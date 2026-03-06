"""
Clinical Metadata Curation and BV-BRC Genome Ingestion Module.

Replaces legacy NCBI Entrez API with direct BV-BRC integration for:
1. Metadata enrichment (collection date, location, host) from BioSample XML
2. Modernity filtering (post-2015 cutoff with 1905 typo rescue)
3. Geographic filtering (India-only clinical isolates)
4. FTP-based genome assembly downloads from BV-BRC
"""

import logging
import shutil
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import urllib.error

import pandas as pd

logger = logging.getLogger(__name__)


class ClinicalMetadataCurator:
    """
    BV-BRC integration engine for clinical strain metadata + genome ingestion.
    
    Workflow:
    1. Read BV-BRC CSV and map to BioSample records
    2. Extract location, collection date, host, isolation source
    3. Apply modernity filter (2015+ with 1905 typo rescue)
    4. Keep only India-anchored strains
    5. Download nucleotide assemblies (.fna) from FTP
    """

    def __init__(
        self,
        email: str = "vihaankulkarni29@gmail.com",
        api_key: str = "c734d921ce3764ef09fcab1c1499e4f41508",
        genomes_dir: Optional[Path] = None,
        results_dir: Optional[Path] = None,
    ):
        """
        Initialize ClinicalMetadataCurator.

        Args:
            email: Email for NCBI API (courteous pacing)
            api_key: NCBI API key
            genomes_dir: Output directory for downloaded .fna files
            results_dir: Output directory for cleaned metadata CSV
        """
        self.email = email
        self.api_key = api_key
        self.genomes_dir = Path(genomes_dir or "data/genomes")
        self.results_dir = Path(results_dir or "data/results")

        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Initialized ClinicalMetadataCurator")
        logger.info(f"  Genomes output: {self.genomes_dir}")
        logger.info(f"  Results output: {self.results_dir}")

    @staticmethod
    def _clean_and_filter_date(date_str: str) -> str:
        """
        Parse date string and determine filter status.

        Rules:
        - 1905 → Rescued (Typo) [known THSTI 2023 data entry error]
        - 2015+ → Pass
        - <2015 → Too Old
        - Invalid/missing → Unknown

        Args:
            date_str: Raw date string from CSV

        Returns:
            Filter status: 'Pass', 'Rescued (Typo)', 'Too Old', 'Unknown'
        """
        if pd.isna(date_str) or date_str == "N/A" or date_str == "missing":
            return "Unknown"

        try:
            year = int(str(date_str)[:4])
            if year == 1905:
                return "Rescued (Typo)"
            if year >= 2015:
                return "Pass"
            return "Too Old"
        except (ValueError, IndexError, TypeError):
            return "Unknown"

    def process_bvbrc_csv(self, csv_path: Path) -> pd.DataFrame:
        """
        Read BV-BRC CSV, apply quality filters, return cleaned metadata.

        Filters applied:
        1. Modernity filter: Keep 2015+, rescue 1905 typo, drop older
        2. Geographic filter: Keep only India-anchored strains
        3. Retain essential columns: Strain, BioSample_ID, Location, Date, Host, Source, Disease

        Args:
            csv_path: Path to BV-BRC metadata CSV (or pre-filtered dataset)

        Returns:
            Cleaned DataFrame with curated clinical strains
        """
        logger.info(f"Processing BV-BRC CSV: {csv_path}")

        df = pd.read_csv(csv_path)
        original_count = len(df)
        logger.info(f"Loaded {original_count} records")

        # Apply modernity filter
        df["Date_Check"] = df["Collection_Date"].apply(
            self._clean_and_filter_date
        )
        keep_status = ["Pass", "Rescued (Typo)", "Unknown"]
        df_modern = df[df["Date_Check"].isin(keep_status)].copy()
        ancient_dropped = len(df[df["Date_Check"] == "Too Old"])

        logger.info(f"Modernity filter: Dropped {ancient_dropped} pre-2015 strains")
        logger.info(f"  After modernity filter: {len(df_modern)} strains")

        # Apply geographic filter (India-only)
        if "Exact_Location" in df_modern.columns:
            df_filtered = df_modern[
                df_modern["Exact_Location"].str.contains("India", na=False, case=False)
            ].copy()
            non_india_dropped = len(df_modern) - len(df_filtered)
            logger.info(f"Geographic filter: Dropped {non_india_dropped} non-India strains")
        else:
            logger.warning("Exact_Location column not found; skipping geographic filter")
            df_filtered = df_modern.copy()

        # Clean up
        df_filtered = df_filtered.drop(columns=["Date_Check"])
        logger.info(f"Final curated dataset: {len(df_filtered)} Indian clinical strains")

        return df_filtered

    def download_bvbrc_genomes(
        self,
        cleaned_df: pd.DataFrame,
        genome_id_column: str = "Genome ID",
        retry_count: int = 3,
        timeout: int = 30,
    ) -> Tuple[int, int]:
        """
        Download nucleotide assemblies (.fna) from BV-BRC FTP for all strains.

        FTP Endpoint: ftp://ftp.bvbrc.org/genomes/{genome_id}/{genome_id}.fna

        Args:
            cleaned_df: Cleaned metadata DataFrame with Genome ID column
            genome_id_column: Name of column containing BV-BRC Genome IDs
            retry_count: Number of retries per failed genome
            timeout: FTP download timeout in seconds

        Returns:
            Tuple of (successful_downloads, failed_downloads)
        """
        logger.info(f"Starting FTP downloads from BV-BRC for {len(cleaned_df)} genomes")

        if genome_id_column not in cleaned_df.columns:
            logger.error(f"Column '{genome_id_column}' not found in DataFrame")
            logger.info("Available columns: " + ", ".join(cleaned_df.columns))
            return 0, 0

        success_count = 0
        fail_count = 0

        for idx, row in cleaned_df.iterrows():
            genome_id = str(row[genome_id_column]).strip()
            strain_name = row.get("Strain", genome_id)

            ftp_url = f"ftp://ftp.bvbrc.org/genomes/{genome_id}/{genome_id}.fna"
            output_path = self.genomes_dir / f"{genome_id}.fna"

            # Skip if already downloaded
            if output_path.exists():
                logger.info(f"[{idx+1}/{len(cleaned_df)}] SKIP (already exists): {strain_name}")
                success_count += 1
                continue

            # Attempt download with retries
            logger.info(f"[{idx+1}/{len(cleaned_df)}] Downloading: {strain_name} ({genome_id})")
            downloaded = False

            for attempt in range(1, retry_count + 1):
                try:
                    with urllib.request.urlopen(ftp_url, timeout=timeout) as response, open(output_path, 'wb') as out_file:
                        shutil.copyfileobj(response, out_file)
                    logger.info(f"  └─ SUCCESS (attempt {attempt}/{retry_count})")
                    success_count += 1
                    downloaded = True
                    break
                except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
                    if attempt < retry_count:
                        logger.warning(f"  └─ RETRY {attempt}/{retry_count}: {type(e).__name__}")
                        time.sleep(1)  # Brief backoff
                    else:
                        logger.error(f"  └─ FAILED after {retry_count} attempts: {e}")
                        fail_count += 1

            # Safe API pacing
            time.sleep(0.15)

        logger.info(f"FTP download complete: {success_count} success, {fail_count} failures")
        return success_count, fail_count

    def run_full_workflow(
        self,
        input_csv: Path,
        genome_id_column: str = "Genome ID",
    ) -> pd.DataFrame:
        """
        Execute complete workflow: read → filter → download → return cleaned dataframe.

        Args:
            input_csv: Path to BV-BRC or pre-filtered metadata CSV
            genome_id_column: Column name for BV-BRC Genome IDs

        Returns:
            Cleaned DataFrame with only successfully-downloaded genomes
        """
        logger.info("="*70)
        logger.info("Clinical Metadata Curation Workflow START")
        logger.info("="*70)

        # Step 1: Process and filter metadata
        cleaned_df = self.process_bvbrc_csv(input_csv)

        # Step 2: Download genomes
        success, fail = self.download_bvbrc_genomes(
            cleaned_df, 
            genome_id_column=genome_id_column
        )

        # Step 3: Filter to only successfully-downloaded genomes
        downloaded_genome_ids = set()
        for fna_file in self.genomes_dir.glob("*.fna"):
            downloaded_genome_ids.add(fna_file.stem)

        final_df = cleaned_df[
            cleaned_df[genome_id_column].astype(str).isin(downloaded_genome_ids)
        ].copy()

        logger.info(f"Final dataset: {len(final_df)} strains with downloaded genomes")

        # Save cleaned metadata
        output_csv = self.results_dir / "final_clinical_dataset.csv"
        final_df.to_csv(output_csv, index=False)
        logger.info(f"Saved cleaned metadata: {output_csv}")

        logger.info("="*70)
        logger.info("Clinical Metadata Curation Workflow COMPLETE")
        logger.info("="*70)

        return final_df
