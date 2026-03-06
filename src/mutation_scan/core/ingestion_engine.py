"""
Genomic Ingestion Engine - Flexible Input Resolution Module

Implements "Flexible Input -> Strict Internal Resolution" pattern.
Accepts BioSample IDs, GenBank Accessions, or BV-BRC Genome IDs.
Normalizes all identifiers to a Resolved_BioSample column for downstream processing.
"""

import json
import logging
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)


class GenomicIngestionEngine:
    """
    Production-grade genomic data ingestion with flexible identifier resolution.
    
    Workflow:
    1. Read input CSV with flexible column detection
    2. Normalize column headers (case-insensitive, whitespace-stripped)
    3. Route based on detected identifier type:
       - BioSample ID -> Direct passthrough
       - GenBank Accession -> Resolve to BioSample via NCBI API
       - BV-BRC Genome ID -> Mark as BVBRC_NATIVE for alternate workflow
    4. Return DataFrame with standardized Resolved_BioSample column
    """

    def __init__(
        self,
        email: str = "user@example.com",
        api_key: Optional[str] = None,
    ):
        """
        Initialize GenomicIngestionEngine.

        Args:
            email: Email for NCBI API requests (required for Entrez)
            api_key: Optional NCBI API key for higher rate limits
        """
        self.email = email
        self.api_key = api_key
        
        logger.info("Initialized GenomicIngestionEngine")
        logger.info(f"  NCBI Email: {self.email}")
        logger.info(f"  API Key: {'SET' if self.api_key else 'NOT SET'}")

    def route_and_resolve_input(
        self, 
        input_csv: Union[str, Path]
    ) -> pd.DataFrame:
        """
        Main router: Read CSV, detect identifier type, resolve to BioSample.

        Logic Flow:
        1. If BioSample column exists -> Direct passthrough to Resolved_BioSample
        2. If Accession column exists -> Resolve via NCBI API to BioSample
        3. If Genome ID column exists -> Mark as BVBRC_NATIVE (alternate workflow)
        4. Else -> Raise ValueError (no recognized identifier columns)

        Args:
            input_csv: Path to input CSV file

        Returns:
            DataFrame with Resolved_BioSample column added

        Raises:
            ValueError: If no recognized identifier columns found
        """
        logger.info("="*70)
        logger.info("GENOMIC INGESTION ENGINE - INPUT RESOLUTION")
        logger.info("="*70)
        logger.info(f"Reading input CSV: {input_csv}")
        
        # Read input CSV
        df = pd.read_csv(input_csv)
        original_count = len(df)
        logger.info(f"Loaded {original_count} records")
        
        # Normalize column headers for flexible detection
        col_map = {str(c).strip().lower(): c for c in df.columns}
        logger.info(f"Detected columns: {list(df.columns)}")
        
        # Route based on detected identifier type
        # Priority: BioSample > Accession > Genome ID
        
        # BRANCH 1: BioSample ID (direct passthrough)
        if 'biosample' in col_map or 'biosample_id' in col_map or 'biosample id' in col_map:
            biosample_col = col_map.get('biosample') or col_map.get('biosample_id') or col_map.get('biosample id')
            logger.info(f"Detected BioSample identifier column: '{biosample_col}'")
            logger.info("Route: DIRECT PASSTHROUGH (no resolution needed)")
            
            df['Resolved_BioSample'] = df[biosample_col].astype(str).str.strip()
            logger.info(f"Successfully resolved {len(df)} BioSample identifiers")
            return df
        
        # BRANCH 2: GenBank Accession (resolve via NCBI API)
        elif 'accession' in col_map or 'assembly_accession' in col_map or 'genbank_accession' in col_map:
            accession_col = col_map.get('accession') or col_map.get('assembly_accession') or col_map.get('genbank_accession')
            logger.info(f"Detected GenBank Accession column: '{accession_col}'")
            logger.info("Route: NCBI API RESOLUTION (accession -> BioSample)")
            
            # Resolve each accession to BioSample
            logger.info("Starting NCBI API resolution (this may take several minutes)...")
            df['Resolved_BioSample'] = df[accession_col].apply(self._resolve_accession_to_biosample)
            
            # Drop rows that failed to resolve
            failed_count = df['Resolved_BioSample'].isna().sum()
            if failed_count > 0:
                logger.warning(f"Failed to resolve {failed_count} accessions. Dropping these rows.")
                df = df.dropna(subset=['Resolved_BioSample']).copy()
            
            logger.info(f"Successfully resolved {len(df)}/{original_count} accessions to BioSample IDs")
            return df
        
        # BRANCH 3: BV-BRC Genome ID (alternate workflow marker)
        elif 'genome id' in col_map or 'genome_id' in col_map or 'bvbrc_id' in col_map:
            genome_id_col = col_map.get('genome id') or col_map.get('genome_id') or col_map.get('bvbrc_id')
            logger.info(f"Detected BV-BRC Genome ID column: '{genome_id_col}'")
            logger.info("Route: BVBRC_NATIVE (alternate workflow - no BioSample resolution)")
            
            df['Resolved_BioSample'] = "BVBRC_NATIVE"
            logger.info(f"Marked {len(df)} records as BVBRC_NATIVE workflow")
            return df
        
        # No recognized identifier columns found
        else:
            error_msg = (
                f"No recognized identifier columns found in CSV. "
                f"Detected columns: {list(df.columns)}. "
                f"Expected one of: 'BioSample', 'BioSample_ID', 'Accession', "
                f"'Assembly_Accession', 'GenBank_Accession', 'Genome_ID', 'Genome ID', 'BVBRC_ID'"
            )
            logger.error(error_msg)
            raise ValueError(error_msg)

    def _resolve_accession_to_biosample(self, accession: str) -> Optional[str]:
        """
        Resolve GenBank Assembly Accession to BioSample ID via NCBI Entrez API.

        Uses two-step process:
        1. Search assembly database to get internal ID
        2. Fetch summary to extract BioSample accession

        Args:
            accession: GenBank Assembly Accession (e.g., GCF_000005845.2)

        Returns:
            BioSample ID (e.g., SAMN02604091) or None if resolution fails
        """
        try:
            # Strip version number from accession
            base_acc = str(accession).strip().split('.')[0]
            
            # Step 1: Search assembly database to get internal ID
            search_url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                f"?db=assembly&term={base_acc}&retmode=json&email={self.email}"
            )
            if self.api_key:
                search_url += f"&api_key={self.api_key}"
            
            with urllib.request.urlopen(search_url, timeout=10) as response:
                search_data = json.loads(response.read().decode('utf-8'))
            
            # Extract assembly ID from search results
            id_list = search_data.get('esearchresult', {}).get('idlist', [])
            if not id_list:
                logger.warning(f"No assembly ID found for accession: {base_acc}")
                return None
            
            assembly_id = id_list[0]
            
            # NCBI API pacing (STRICT CONSTRAINT)
            time.sleep(0.15)
            
            # Step 2: Fetch summary to extract BioSample accession
            summary_url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=assembly&id={assembly_id}&retmode=json&email={self.email}"
            )
            if self.api_key:
                summary_url += f"&api_key={self.api_key}"
            
            with urllib.request.urlopen(summary_url, timeout=10) as response:
                summary_data = json.loads(response.read().decode('utf-8'))
            
            # Extract BioSample accession from summary
            result = summary_data.get('result', {})
            biosample_acc = result.get(assembly_id, {}).get('biosampleaccn')
            
            if biosample_acc:
                logger.debug(f"Resolved {base_acc} -> {biosample_acc}")
                return biosample_acc
            else:
                logger.warning(f"No BioSample found in summary for: {base_acc}")
                return None
                
        except urllib.error.URLError as e:
            logger.error(f"Network error resolving {accession}: {e}")
            return None
        except urllib.error.HTTPError as e:
            logger.error(f"HTTP error resolving {accession}: {e}")
            return None
        except (KeyError, json.JSONDecodeError) as e:
            logger.error(f"JSON parsing error resolving {accession}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error resolving {accession}: {e}")
            return None
