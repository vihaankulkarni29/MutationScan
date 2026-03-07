"""
Universal Genome Downloader - Multi-Source Assembly Acquisition Module

Phase 3 of the genomic ingestion pipeline.
Downloads nucleotide assemblies (.fna) from curated metadata using two pathways:
- BV-BRC HTTPS: Native BV-BRC strains via secure FTP
- NCBI Datasets REST API: BioSample-resolved strains via official REST endpoint

Performance: Throttled Handshake with Parallel Payloads pattern
- 8 concurrent workers for maximum throughput
- Thread lock on NCBI API handshakes to enforce 3 requests/second limit
- Unlocked parallel downloads of heavy payloads
"""

import io
import json
import logging
import os
import shutil
import threading
import time
import urllib.parse
import urllib.request
import urllib.error
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple, Optional

import pandas as pd

logger = logging.getLogger(__name__)


class UniversalGenomeDownloader:
    """
    Production-grade genome downloader with dual-source resolution and parallel execution.
    
    Workflow:
    1. Read curated_metadata.csv from Phase 2 (only strains that passed filters)
    2. Spawn 8 concurrent workers for dual-pathway downloading:
       - If Resolved_BioSample == "BVBRC_NATIVE": Download from BV-BRC HTTPS
       - If Resolved_BioSample is standard ID: Use NCBI Datasets REST API
    3. Extract and save .fna files to data/genomes/
    4. Return success/fail statistics for pipeline monitoring
    
    Thread Safety:
    - NCBI handshake (Step A) is protected by ncbi_lock to enforce 0.35s pacing
    - Payload downloads (Step B) run in parallel without lock
    - BV-BRC downloads are lock-free with individual timeout handling
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        genomes_dir: Optional[Path] = None,
        max_retries: int = 3,
    ):
        """
        Initialize UniversalGenomeDownloader.

        Args:
            api_key: Optional NCBI API key for higher rate limits
            genomes_dir: Output directory for genome assemblies (.fna files)
            max_retries: Maximum retry attempts per download (default: 3)
        """
        self.api_key = api_key
        self.genomes_dir = Path(genomes_dir or "data/genomes")
        self.max_retries = max_retries
        self.ncbi_lock = threading.Lock()  # Protects NCBI handshakes (Step A)
        
        # Ensure output directory exists
        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Initialized UniversalGenomeDownloader")
        logger.info(f"  NCBI API Key: {'SET' if self.api_key else 'NOT SET'}")
        logger.info(f"  Genomes Directory: {self.genomes_dir}")
        logger.info(f"  Max Retries: {self.max_retries}")
        logger.info(f"  NCBI Rate Limit: Protected by thread lock (0.35s pacing)")

    def download_curated_genomes(
        self, 
        curated_csv: Optional[Path] = None
    ) -> Tuple[int, int]:
        """
        Main router: Download all genomes from curated metadata using parallel workers.
        
        Implements "Throttled Handshake with Parallel Payloads" pattern:
        - 8 concurrent workers maximize throughput
        - NCBI API handshakes (DNS, auth) protected by lock to enforce 3 req/s limit
        - Large payload downloads (genome zips) run in parallel without contention

        Args:
            curated_csv: Path to curated_metadata.csv (default: data/results/curated_metadata.csv)

        Returns:
            Tuple of (success_count, fail_count)
        """
        logger.info("="*70)
        logger.info("UNIVERSAL GENOME DOWNLOADER - PARALLEL MULTI-SOURCE ACQUISITION")
        logger.info("="*70)
        
        # Load curated metadata from Phase 2
        curated_path = Path(curated_csv or "data/results/curated_metadata.csv")
        
        if not curated_path.exists():
            raise FileNotFoundError(f"Curated metadata not found: {curated_path}")
        
        df = pd.read_csv(curated_path)
        logger.info(f"Loaded {len(df)} curated strains from {curated_path}")
        
        if 'Resolved_BioSample' not in df.columns:
            raise ValueError("Curated metadata must contain 'Resolved_BioSample' column")
        
        logger.info(f"Launching parallel downloads with 8 concurrent workers...")
        
        success_count = 0
        fail_count = 0
        
        # Worker function that processes one row
        def process_row(idx_row_tuple):
            idx, row = idx_row_tuple
            biosample_id = str(row.get('Resolved_BioSample', '')).strip()
            
            if biosample_id == 'BVBRC_NATIVE':
                return self._download_bvbrc_genome(row, idx, len(df))
            else:
                return self._download_ncbi_genome(biosample_id, idx, len(df))
        
        # Execute downloads in parallel
        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = {executor.submit(process_row, (idx, row)): idx 
                      for idx, row in df.iterrows()}
            
            for future in as_completed(futures):
                try:
                    if future.result():
                        success_count += 1
                    else:
                        fail_count += 1
                except Exception as e:
                    logger.error(f"Worker thread crashed: {e}")
                    fail_count += 1
        
        logger.info("="*70)
        logger.info(f"Download complete: {success_count} successful, {fail_count} failed")
        logger.info("="*70)
        
        return success_count, fail_count

    def _download_bvbrc_genome(
        self, 
        row: pd.Series, 
        idx: int, 
        total: int
    ) -> bool:
        """
        Download genome from BV-BRC via HTTPS.

        Args:
            row: DataFrame row containing Genome ID
            idx: Current row index (for logging)
            total: Total row count (for logging)

        Returns:
            True if download succeeded, False otherwise
        """
        try:
            # Extract Genome ID from the row
            if 'Genome ID' not in row or pd.isna(row['Genome ID']):
                logger.error(f"[{idx+1}/{total}] BVBRC_NATIVE row missing Genome ID column")
                return False
            
            genome_id = str(row['Genome ID']).strip()
            safe_id = urllib.parse.quote(genome_id)
            download_url = f"https://ftp.bvbrc.org/genomes/{safe_id}/{safe_id}.fna"
            output_path = self.genomes_dir / f"{genome_id}.fna"
            
            logger.debug(f"[{idx+1}/{total}] BV-BRC URL: {download_url}")
            
            # Retry loop
            for attempt in range(1, self.max_retries + 1):
                try:
                    logger.debug(f"[{idx+1}/{total}] Attempt {attempt}/{self.max_retries}")
                    
                    with urllib.request.urlopen(download_url, timeout=30) as response:
                        with open(output_path, 'wb') as out_file:
                            shutil.copyfileobj(response, out_file)
                    
                    logger.info(f"[{idx+1}/{total}] BV-BRC download successful: {output_path.name}")
                    return True
                    
                except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
                    logger.warning(f"[{idx+1}/{total}] Attempt {attempt} failed: {e}")
                    if attempt < self.max_retries:
                        time.sleep(2 ** attempt)  # Exponential backoff
                    else:
                        logger.error(f"[{idx+1}/{total}] All {self.max_retries} attempts failed for {genome_id}")
                        return False
        
        except Exception as e:
            logger.error(f"[{idx+1}/{total}] Unexpected error in BV-BRC download: {e}")
            return False

    def _download_ncbi_genome(
        self, 
        biosample_id: str, 
        idx: int, 
        total: int
    ) -> bool:
        """
        Download genome from NCBI via Datasets REST API with rate-limit protection.

        Two-step process with differing synchronization strategies:
        - Step A (Handshake): Protected by ncbi_lock to enforce 0.35s pacing
        - Step B (Payload): Unlocked for parallel execution
        
        Args:
            biosample_id: NCBI BioSample accession (e.g., SAMN02604018)
            idx: Current row index (for logging)
            total: Total row count (for logging)

        Returns:
            True if download succeeded, False otherwise
        """
        output_path = self.genomes_dir / f"{biosample_id}.fna"
        
        # Retry loop
        for attempt in range(1, self.max_retries + 1):
            try:
                logger.debug(f"[{idx+1}/{total}] NCBI Attempt {attempt}/{self.max_retries}")
                
                # STEP A: The Handshake (Strictly Locked & Throttled)
                accession = None
                with self.ncbi_lock:
                    logger.debug(f"[{idx+1}/{total}] Acquiring NCBI handshake lock...")
                    api_url = (
                        f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/"
                        f"Escherichia%20coli/dataset_report?filters.biosample={biosample_id}"
                    )
                    if self.api_key:
                        api_url += f"&api_key={self.api_key}"
                    
                    req = urllib.request.Request(api_url, headers={'Accept': 'application/json'})
                    
                    with urllib.request.urlopen(req, timeout=15) as response:
                        data = json.loads(response.read().decode('utf-8'))
                    
                    reports = data.get('reports', [])
                    if not reports:
                        logger.error(f"[{idx+1}/{total}] No assembly found for BioSample {biosample_id}")
                        return False
                    
                    accession = reports[0]['accession']
                    logger.debug(f"[{idx+1}/{total}] Found assembly: {accession}")
                    
                    # Strict NCBI pacing (3 requests per second = 0.35s per request)
                    time.sleep(0.35)
                
                # STEP B: The Payload (Unlocked, Parallel execution)
                logger.debug(f"[{idx+1}/{total}] Released lock, downloading payload in parallel...")
                download_url = (
                    f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
                    f"{accession}/download?include_annotation_type=GENOME_FASTA"
                )
                if self.api_key:
                    download_url += f"&api_key={self.api_key}"
                
                req = urllib.request.Request(download_url, headers={'Accept': 'application/zip'})
                
                with urllib.request.urlopen(req, timeout=60) as response:
                    with zipfile.ZipFile(io.BytesIO(response.read())) as z:
                        # Find the .fna file inside the zip (usually under ncbi_dataset/data/GCF.../*.fna)
                        fna_filename = next((name for name in z.namelist() if name.endswith('.fna')), None)
                        
                        if not fna_filename:
                            logger.error(f"[{idx+1}/{total}] No .fna file found in downloaded zip for {biosample_id}")
                            return False
                        
                        logger.debug(f"[{idx+1}/{total}] Extracting: {fna_filename}")
                        
                        with z.open(fna_filename) as zf:
                            with open(output_path, 'wb') as f:
                                shutil.copyfileobj(zf, f)
                
                logger.info(f"[{idx+1}/{total}] NCBI download successful: {output_path.name}")
                return True
                
            except urllib.error.HTTPError as e:
                logger.warning(f"[{idx+1}/{total}] HTTP error (attempt {attempt}): {e.code} {e.reason}")
                if attempt < self.max_retries:
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    logger.error(f"[{idx+1}/{total}] All {self.max_retries} attempts failed for {biosample_id}")
                    return False
                    
            except urllib.error.URLError as e:
                logger.warning(f"[{idx+1}/{total}] Network error (attempt {attempt}): {e}")
                if attempt < self.max_retries:
                    time.sleep(2 ** attempt)
                else:
                    logger.error(f"[{idx+1}/{total}] All {self.max_retries} attempts failed for {biosample_id}")
                    return False
                    
            except (json.JSONDecodeError, KeyError, zipfile.BadZipFile) as e:
                logger.error(f"[{idx+1}/{total}] Data parsing error: {e}")
                return False
                
            except Exception as e:
                logger.error(f"[{idx+1}/{total}] Unexpected error in NCBI download: {e}")
                return False
        
        return False
