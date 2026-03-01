"""
NCBI Datasets API v2 handler for genome downloads.

This module provides functionality to query and download bacterial genomes
from NCBI using the modern Datasets API and CLI-based bulk downloading.

Uses Bio.Entrez for text searching and ncbi.datasets.GenomeApi for actual data transfer.
"""

import json
import logging
import re
import shutil
import subprocess
import time
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from io import BytesIO
from pathlib import Path
from typing import List, Dict, Optional, Tuple

from Bio import Entrez
import pandas as pd
import requests

logger = logging.getLogger(__name__)


class BulkGenomeDownloader:
    """
    Enterprise-scale downloader using official ncbi-datasets-cli.

    Workflow:
    - Phase A: Dehydrate genome package from taxon + geo-location
    - Phase B: Unpack dehydrated bundle
    - Phase C: Rehydrate concurrently with datasets CLI
    - Phase D: Organize .fna files into data/genomes as accession-named FASTA files
    """

    def __init__(
        self,
        output_dir: Optional[Path] = None,
        log_file: Optional[Path] = None,
        datasets_binary: str = "datasets",
    ):
        self.output_dir = Path(output_dir or "data/genomes")
        self.log_file = Path(log_file or "data/logs/genome_extractor.log")
        self.datasets_binary = datasets_binary

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        self._setup_logging()

        if not shutil.which(self.datasets_binary):
            raise FileNotFoundError(
                f"'{self.datasets_binary}' was not found in PATH. "
                "Install ncbi-datasets-cli in the runtime environment."
            )

        logger.info("Initialized BulkGenomeDownloader with datasets CLI")
        logger.info(f"Output directory: {self.output_dir}")

    def _setup_logging(self):
        """Configure file-based logging for subprocess progress tracking."""
        existing_file_handlers = [
            handler for handler in logger.handlers if isinstance(handler, logging.FileHandler)
        ]
        if any(getattr(handler, "baseFilename", "") == str(self.log_file.resolve()) for handler in existing_file_handlers):
            return

        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    def _run_command(
        self,
        command: List[str],
        phase_name: str,
        cwd: Optional[Path] = None,
        timeout: int = 0,
    ) -> None:
        """
        Execute a CLI command and route stdout/stderr to logger.

        Raises:
            RuntimeError: If subprocess exits with non-zero status
        """
        command_display = " ".join(command)
        logger.info(f"{phase_name}: Running command: {command_display}")

        completed = subprocess.run(
            command,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            timeout=None if timeout <= 0 else timeout,
            check=False,
        )

        if completed.stdout:
            for line in completed.stdout.splitlines():
                logger.info(f"{phase_name} stdout: {line}")

        if completed.stderr:
            for line in completed.stderr.splitlines():
                logger.warning(f"{phase_name} stderr: {line}")

        if completed.returncode != 0:
            raise RuntimeError(
                f"{phase_name} failed with exit code {completed.returncode}: {command_display}"
            )

    @staticmethod
    def _extract_accession(path_obj: Path) -> Optional[str]:
        """Extract accession ID (GCF_/GCA_) from file path or filename."""
        accession_pattern = re.compile(r"(GC[AF]_\d+\.\d+)")
        path_text = str(path_obj)
        match = accession_pattern.search(path_text)
        if match:
            return match.group(1)

        return None

    def bulk_download_by_accessions(
        self,
        accessions: List[str],
        archive_name: str = "data.zip",
    ) -> Tuple[int, int]:
        """
        Download large genome cohorts via datasets CLI dehydrated + rehydrate flow.

        Args:
            accessions: Explicit list of assembly accessions (GCF_/GCA_)
            archive_name: Filename for dehydrated zip package

        Returns:
            Tuple of (successful_downloads, failed_downloads)
        """
        logger.info(
            f"Starting bulk genome download via datasets CLI using accession input list ({len(accessions)} accessions)"
        )

        if not accessions:
            logger.warning("No accessions provided for bulk download")
            return 0, 0

        successful = 0
        failed = 0

        with tempfile.TemporaryDirectory(prefix="mutation_scan_bulk_") as tmp_dir:
            temp_root = Path(tmp_dir)
            dehydrated_zip = temp_root / archive_name
            unpack_dir = temp_root / "dehydrated"
            accessions_file = temp_root / "accessions.txt"

            accessions_file.write_text("\n".join(accessions), encoding="utf-8")

            try:
                # Phase A: Dehydrate
                dehydrate_command = [
                    self.datasets_binary,
                    "download",
                    "genome",
                    "accession",
                    "--inputfile",
                    str(accessions_file),
                    "--include",
                    "genome",
                    "--dehydrated",
                    "--filename",
                    str(dehydrated_zip),
                ]
                self._run_command(dehydrate_command, "Phase A (Dehydrate)", cwd=temp_root)

                if not dehydrated_zip.exists():
                    raise FileNotFoundError(f"Dehydrated archive not found: {dehydrated_zip}")

                # Phase B: Unpack
                unpack_dir.mkdir(parents=True, exist_ok=True)
                with zipfile.ZipFile(dehydrated_zip, "r") as zip_ref:
                    zip_ref.extractall(unpack_dir)
                logger.info(f"Phase B (Unpack): Extracted dehydrated archive to {unpack_dir}")

                # Phase C: Rehydrate
                rehydrate_command = [
                    self.datasets_binary,
                    "rehydrate",
                    "--directory",
                    str(unpack_dir),
                ]
                self._run_command(rehydrate_command, "Phase C (Rehydrate)", cwd=temp_root)

                # Phase D: Organize
                fna_files = list(unpack_dir.rglob("*.fna"))
                if not fna_files:
                    logger.warning("Phase D (Organize): No .fna files found after rehydration")
                    return 0, 0

                logger.info(f"Phase D (Organize): Found {len(fna_files)} .fna files")

                for fna_file in fna_files:
                    accession = self._extract_accession(fna_file)
                    if not accession:
                        logger.warning(f"Could not infer accession from path: {fna_file}")
                        failed += 1
                        continue

                    destination = self.output_dir / f"{accession}.fasta"
                    if destination.exists():
                        logger.info(f"Overwriting existing genome file: {destination}")

                    shutil.move(str(fna_file), str(destination))
                    successful += 1

                logger.info(
                    f"Bulk download complete via datasets CLI: {successful} successful, {failed} failed"
                )
                return successful, failed

            except Exception as exc:
                logger.error(f"Bulk download failed: {type(exc).__name__}: {exc}")
                raise


class NCBIDatasetsGenomeDownloader:
    """
    Download genomes using NCBI Datasets API v2 (modern, fast, stable).
    
    Supports two modes:
    - Search Mode (A): Query-based (e.g., "Escherichia coli AND Antibiotic Resistant")
    - File Mode (B): Accession list from .txt file
    
    Outputs:
    - Individual FASTA files named {Accession}.fasta
    - metadata_master.csv with clinical metadata
    - genome_extractor.log for error tracking
    """

    DATASETS_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    NCBI_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    
    def __init__(
        self,
        email: str,
        api_key: Optional[str] = None,
        output_dir: Optional[Path] = None,
        log_file: Optional[Path] = None,
    ):
        """
        Initialize the NCBI Datasets genome downloader.

        Args:
            email: Email for NCBI (required by NCBI policy)
            api_key: Optional NCBI API key for faster requests
            output_dir: Directory to save genomes and metadata
            log_file: Path to error log file

        Raises:
            ValueError: If email is not provided
        """
        if not email:
            raise ValueError("Email is required for NCBI API access")
        
        self.email = email
        self.api_key = api_key
        self.output_dir = Path(output_dir or "data/genomes")
        self.log_file = Path(log_file or "data/logs/genome_extractor.log")
        
        # Configure Entrez
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        
        # Create output directories BEFORE logging setup
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self._setup_logging()
        
        logger.info(f"Initialized NCBIDatasetsGenomeDownloader for {email}")
        logger.info(f"Output directory: {self.output_dir}")

    def _setup_logging(self):
        """Configure file-based logging for error tracking."""
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    def search_accessions(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search NCBI Assembly database and resolve to Accession IDs.
        
        Uses Entrez.esearch -> Entrez.esummary for robust text searching.
        Only returns Assembly Accessions (GCF_/GCA_), not Nucleotide UIDs.

        Args:
            query: Search query (e.g., "Escherichia coli AND Antibiotic Resistant")
            max_results: Maximum number of accessions to return

        Returns:
            List of Assembly Accessions (GCF_/GCA_ format)

        Raises:
            Exception: If search fails after retries
        """
        logger.info(f"Searching NCBI Assembly for: {query}")
        
        accessions = []
        retry_count = 0
        max_retries = 3
        
        while retry_count < max_retries:
            try:
                # Step 1: esearch to get Assembly UIDs (returns XML)
                search_params = {
                    "db": "assembly",
                    "term": query,
                    "retmax": max_results,
                    "tool": "mutation_scan",
                    "email": self.email,
                }
                if self.api_key:
                    search_params["api_key"] = self.api_key
                
                response = requests.get(self.NCBI_ESEARCH_URL, params=search_params)
                response.raise_for_status()
                
                # Parse XML response
                root = ET.fromstring(response.text)
                id_list = root.find('.//IdList')
                uids = [id_elem.text for id_elem in id_list.findall('Id')] if id_list is not None else []
                
                if not uids:
                    logger.warning(f"No results found for query: {query}")
                    return []
                
                logger.info(f"Found {len(uids)} Assembly UIDs. Resolving to accessions...")
                
                # Step 2: esummary to get Assembly Accessions (returns XML)
                summary_params = {
                    "db": "assembly",
                    "id": ",".join(uids[:max_results]),
                    "tool": "mutation_scan",
                    "email": self.email,
                }
                if self.api_key:
                    summary_params["api_key"] = self.api_key
                
                response = requests.get(self.NCBI_ESUMMARY_URL, params=summary_params)
                response.raise_for_status()
                
                # Parse XML for assembly accessions
                root = ET.fromstring(response.text)
                for doc_sum in root.findall('.//DocumentSummary'):
                    # Find AssemblyAccession field
                    asm_acc_elem = doc_sum.find(".//AssemblyAccession")
                    if asm_acc_elem is not None:
                        asm_accession = asm_acc_elem.text
                        if asm_accession and (asm_accession.startswith("GCF_") or asm_accession.startswith("GCA_")):
                            accessions.append(asm_accession)
                
                logger.info(f"Resolved {len(accessions)} valid Assembly Accessions")
                return accessions
                
            except requests.exceptions.RequestException as e:
                retry_count += 1
                logger.error(f"Search attempt {retry_count} failed: {e}")
                if retry_count < max_retries:
                    wait_time = 2 ** retry_count
                    logger.info(f"Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"Search failed after {max_retries} retries")
                    raise
            except Exception as e:
                logger.error(f"Unexpected error during search: {e}")
                raise

    def read_accession_file(self, filepath: Path) -> List[str]:
        """
        Read accessions from a text file (one per line).

        Args:
            filepath: Path to accession list file

        Returns:
            List of Accessions

        Raises:
            FileNotFoundError: If file does not exist
        """
        logger.info(f"Reading accessions from: {filepath}")
        
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Accession file not found: {filepath}")
        
        accessions = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):  # Skip empty and comment lines
                    accessions.append(line)
        
        logger.info(f"Loaded {len(accessions)} accessions from file")
        return accessions

    def download_bulk_by_geolocation(self, organism: str, location: str, limit: int = 5000) -> Tuple[int, int]:
        """
        Hybrid geolocation downloader:
        1) Query Entrez for assembly accessions using organism + location
        2) Feed explicit accession list into datasets CLI bulk downloader

        Args:
            organism: Taxon string for Entrez query
            location: Geographic filter token for Entrez text search
            limit: Maximum number of accessions to retrieve from Entrez

        Returns:
            Tuple of (successful_downloads, failed_downloads)
        """
        query = f'"{organism}"[Organism] AND "{location}"[Text Word]'
        logger.info(f"Querying Entrez for Accessions: {query}")

        accessions = self.search_accessions(query, max_results=limit)
        if not accessions:
            logger.warning("No accessions found for hybrid geolocation query")
            return 0, 0

        bulk_downloader = BulkGenomeDownloader(
            output_dir=self.output_dir,
            log_file=self.log_file,
        )
        return bulk_downloader.bulk_download_by_accessions(accessions)

    def download_batch(self, accessions: List[str]) -> Tuple[int, int]:
        """
        Batch download genomes using NCBI Datasets API v2.
        
        Downloads full assemblies (chromosome + plasmids) to avoid missing resistance genes.
        Automatically processes zip and extracts metadata.

        Args:
            accessions: List of Assembly Accessions (GCF_/GCA_)

        Returns:
            Tuple of (successful_downloads, failed_downloads)
        """
        logger.info(f"Starting batch download of {len(accessions)} genomes...")
        
        successful = 0
        failed = 0
        metadata_list = []
        
        for i, accession in enumerate(accessions, 1):
            try:
                logger.info(f"[{i}/{len(accessions)}] Downloading {accession}...")
                
                # Download from Datasets API
                zip_data = self._download_from_datasets(accession)
                if not zip_data:
                    failed += 1
                    continue
                
                # Process zip and extract metadata
                fasta_file, metadata_dict = self._process_zip_and_metadata(
                    zip_data, accession
                )
                
                if fasta_file:
                    successful += 1
                    metadata_list.append(metadata_dict)
                    logger.info(f"Successfully processed {accession}")
                else:
                    failed += 1
                
            except Exception as e:
                logger.error(f"Failed to download {accession}: {e}")
                failed += 1
                continue
        
        # Save master metadata CSV
        if metadata_list:
            self._save_metadata_master(metadata_list)
        
        logger.info(f"Batch download complete: {successful} successful, {failed} failed")
        return successful, failed

    def _download_from_datasets(self, accession: str) -> Optional[BytesIO]:
        """
        Download genome assembly from NCBI Datasets API v2.

        Args:
            accession: Assembly Accession (GCF_/GCA_)

        Returns:
            BytesIO object containing the zip file, or None if failed

        Raises:
            Exception: If download fails after retries
        """
        url = f"{self.DATASETS_API_URL}/genome/accession/{accession}/download"
        
        params = {
            "include_annotation_type": "GENOME_FASTA",  # Include chromosome + plasmids
        }
        
        # Only add api_key if it exists (empty string causes 400 errors)
        if self.api_key:
            params["api_key"] = self.api_key
        
        retry_count = 0
        max_retries = 3
        
        while retry_count < max_retries:
            try:
                response = requests.get(url, params=params, timeout=60, stream=True)
                response.raise_for_status()
                
                # Read the zip file into memory
                zip_buffer = BytesIO()
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        zip_buffer.write(chunk)
                
                zip_buffer.seek(0)
                logger.debug(f"Downloaded {len(zip_buffer.getvalue())} bytes for {accession}")
                return zip_buffer
                
            except requests.exceptions.RequestException as e:
                retry_count += 1
                logger.warning(f"Download attempt {retry_count} failed for {accession}: {e}")
                if retry_count < max_retries:
                    time.sleep(2 ** retry_count)
                else:
                    logger.error(f"Download failed after {max_retries} retries for {accession}")
                    return None

    def _process_zip_and_metadata(
        self, zip_data: BytesIO, accession: str
    ) -> Tuple[Optional[Path], Dict]:
        """
        Extract FASTA from zip and parse metadata_report.jsonl.
        
        Simple QC: Checks if file is populated (> 1KB).

        Args:
            zip_data: BytesIO object of the downloaded zip
            accession: Assembly Accession for naming

        Returns:
            Tuple of (fasta_filepath, metadata_dict)
            Returns (None, metadata_dict) if QC fails
        """
        metadata_dict = {
            "Accession": accession,
            "Organism Name": "N/A",
            "Strain": "N/A",
            "Collection Date": "N/A",
            "Host": "N/A",
            "Isolation Source": "N/A",
            "Geo Location": "N/A",
            "QC Status": "PASS",
            "Downloaded": False,
        }
        
        try:
            with zipfile.ZipFile(zip_data, 'r') as zf:
                file_list = zf.namelist()
                
                # 1. Find and Validate FASTA
                fna_file = next((f for f in file_list if f.endswith(".fna")), None)
                if not fna_file:
                    logger.error(f"No .fna file found for {accession}")
                    return None, metadata_dict
                
                fasta_content = zf.read(fna_file).decode('utf-8')
                
                # Simple QC: Is the file actually populated? (> 1KB)
                if len(fasta_content) < 1000:
                    logger.warning(f"QC FAIL: {accession} is empty or too small.")
                    metadata_dict["QC Status"] = "QC_FAIL"
                    return None, metadata_dict

                # 2. Parse Metadata
                jsonl_file = next((f for f in file_list if "data_report.jsonl" in f), None)
                if jsonl_file:
                    metadata_dict.update(self._parse_jsonl_metadata(zf.read(jsonl_file).decode('utf-8')))
                
                # 3. Save File
                output_file = self.output_dir / f"{accession}.fasta"
                with open(output_file, 'w') as f:
                    f.write(fasta_content)
                
                metadata_dict["Downloaded"] = True
                logger.info(f"Saved {accession}.fasta ({len(fasta_content)} bytes)")
                return output_file, metadata_dict
                
        except Exception as e:
            logger.error(f"Error processing {accession}: {e}")
            return None, metadata_dict

    def _parse_jsonl_metadata(self, jsonl_content: str) -> Dict:
        """
        Parse NCBI data_report.jsonl to extract clinical metadata.
        
        Navigates nested JSON structure: assembly_info -> biosample -> attributes.
        Defaults to "N/A" for missing fields (anti-hallucination).

        Args:
            jsonl_content: Raw JSONL content from data_report.jsonl

        Returns:
            Dictionary with clinical metadata
        """
        metadata = {}
        try:
            # Parse only the first line (first record)
            line = jsonl_content.strip().split('\n')[0]
            record = json.loads(line)
            
            # 1. Extract Assembly Info
            asm_info = record.get("assembly_info", {})
            biosample = asm_info.get("biosample", {})
            
            metadata["Organism Name"] = record.get("organism", {}).get("organism_name", "N/A")
            metadata["Strain"] = biosample.get("strain", "N/A")
            
            # 2. Extract Attributes (The tricky part)
            # attributes is a list: [{'name': 'host', 'value': 'Homo sapiens'}, ...]
            attributes = biosample.get("attributes", [])
            attr_map = {item["name"]: item["value"] for item in attributes}
            
            metadata["Collection Date"] = attr_map.get("collection_date", "N/A")
            metadata["Host"] = attr_map.get("host", "N/A")
            metadata["Isolation Source"] = attr_map.get("isolation_source", "N/A")
            metadata["Geo Location"] = attr_map.get("geo_loc_name", "N/A")
            
        except Exception as e:
            logger.debug(f"Metadata parsing warning: {e}")
            
        return metadata

    def _save_metadata_master(self, metadata_list: List[Dict]) -> Path:
        """
        Save master metadata CSV for all downloaded genomes.

        Args:
            metadata_list: List of metadata dictionaries

        Returns:
            Path to saved CSV file
        """
        df = pd.DataFrame(metadata_list)
        output_file = self.output_dir / "metadata_master.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Saved master metadata to {output_file}")
        return output_file
