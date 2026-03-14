"""
Metadata Interrogator - Epidemiological Data Enrichment and Filtering Module

Phase 2 of the genomic ingestion pipeline.
Fetches metadata from NCBI BioSample, applies geographic/temporal filters,
and partitions data into curated and rejected datasets for scientific integrity.
"""

import logging
import os
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Tuple, Optional

import pandas as pd

logger = logging.getLogger(__name__)


class MetadataInterrogator:
    """
    Production-grade metadata enrichment with scientific filtering.
    
    Workflow:
    1. Accept DataFrame with Resolved_BioSample column from GenomicIngestionEngine
    2. For each BioSample:
       - If BVBRC_NATIVE: Validate existing metadata columns
       - If standard BioSample: Fetch metadata from NCBI XML API
    3. Apply rigorous filters:
       - Geographic: India-only strains
       - Temporal: 2015+ collection dates (with 1905 typo rescue)
    4. Partition into curated (passed) and rejected datasets
    5. Save both datasets with rejection reasons for audit trail
    """

    def __init__(
        self,
        email: str = "user@example.com",
        api_key: Optional[str] = None,
        output_dir: Optional[Path] = None,
        filter_country: Optional[str] = "india",
        filter_min_year: Optional[int] = 2015,
    ):
        """
        Initialize MetadataInterrogator.

        Args:
            email: Email for NCBI API requests (required for Entrez)
            api_key: Optional NCBI API key for higher rate limits
            output_dir: Output directory for curated/rejected datasets
                    filter_country: Geographic filter (case-insensitive). Set to None to disable.
                    filter_min_year: Minimum collection year filter. Set to None to disable.
        """
        self.email = email
        self.api_key = api_key
        self.output_dir = Path(output_dir or "data/results")
        self.filter_country = filter_country.lower() if filter_country else None
        self.filter_min_year = filter_min_year
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Initialized MetadataInterrogator")
        logger.info(f"  NCBI Email: {self.email}")
        logger.info(f"  API Key: {'SET' if self.api_key else 'NOT SET'}")
        logger.info(f"  Output Directory: {self.output_dir}")

    def interrogate_and_filter(
        self, 
        resolved_df: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Main router: Fetch metadata, apply filters, partition datasets.

        Args:
            resolved_df: DataFrame with Resolved_BioSample column from GenomicIngestionEngine

        Returns:
            Tuple of (curated_df, rejected_df) with full metadata and rejection reasons
        """
        logger.info("="*70)
        logger.info("METADATA INTERROGATOR - EPIDEMIOLOGICAL ENRICHMENT & FILTERING")
        logger.info("="*70)
        logger.info(f"Processing {len(resolved_df)} strains")
        
        if 'Resolved_BioSample' not in resolved_df.columns:
            raise ValueError("Input DataFrame must contain 'Resolved_BioSample' column from GenomicIngestionEngine")
        
        # Initialize result lists
        enriched_records = []
        
        for idx, row in resolved_df.iterrows():
            biosample_id = row['Resolved_BioSample']
            
            # BRANCH 1: BVBRC_NATIVE workflow (validate existing columns)
            if str(biosample_id) == "BVBRC_NATIVE":
                logger.debug(f"[{idx+1}/{len(resolved_df)}] BVBRC_NATIVE: Validating existing metadata")
                
                location = row.get('Exact_Location', 'Unknown')
                collection_date = row.get('Collection_Date', 'Unknown')
                host = row.get('Host', 'Unknown')
                isolation_source = row.get('Isolation_Source', 'Unknown')
                
                # Validate with existing metadata
                passed, reason = self._validate_strain(location, collection_date)
                
                enriched_records.append({
                    **row.to_dict(),
                    'BioSample_ID': 'BVBRC_NATIVE',
                    'Geo_Loc_Name': location,
                    'Collection_Date': collection_date,
                    'Host': host,
                    'Isolation_Source': isolation_source,
                    'Filter_Status': reason,
                    'Passed_Filters': passed,
                })
            
            # BRANCH 2: Standard BioSample workflow (fetch from NCBI)
            else:
                logger.debug(f"[{idx+1}/{len(resolved_df)}] Fetching metadata for BioSample: {biosample_id}")
                
                # Fetch metadata from NCBI
                metadata = self._fetch_biosample_metadata(biosample_id)
                
                if metadata is None:
                    # API fetch failed
                    enriched_records.append({
                        **row.to_dict(),
                        'BioSample_ID': biosample_id,
                        'Geo_Loc_Name': 'Unknown',
                        'Collection_Date': 'Unknown',
                        'Host': 'Unknown',
                        'Isolation_Source': 'Unknown',
                        'Filter_Status': 'API Fetch Failed',
                        'Passed_Filters': False,
                    })
                    continue
                
                # Validate fetched metadata
                passed, reason = self._validate_strain(
                    metadata.get('geo_loc_name', 'Unknown'),
                    metadata.get('collection_date', 'Unknown')
                )
                
                enriched_records.append({
                    **row.to_dict(),
                    'BioSample_ID': biosample_id,
                    'Geo_Loc_Name': metadata.get('geo_loc_name', 'Unknown'),
                    'Collection_Date': metadata.get('collection_date', 'Unknown'),
                    'Host': metadata.get('host', 'Unknown'),
                    'Isolation_Source': metadata.get('isolation_source', 'Unknown'),
                    'Filter_Status': reason,
                    'Passed_Filters': passed,
                })
        
        # Create enriched DataFrame
        enriched_df = pd.DataFrame(enriched_records)
        
        # Partition into curated (passed) and rejected datasets
        curated_df = enriched_df[enriched_df['Passed_Filters'] == True].copy()
        rejected_df = enriched_df[enriched_df['Passed_Filters'] == False].copy()
        
        passed_count = len(curated_df)
        rejected_count = len(rejected_df)
        
        logger.info(f"Filtering complete: {passed_count} passed, {rejected_count} rejected")
        
        # Log rejection breakdown
        if rejected_count > 0:
            rejection_summary = rejected_df['Filter_Status'].value_counts()
            logger.info("Rejection breakdown:")
            for reason, count in rejection_summary.items():
                logger.info(f"  {reason}: {count}")
        
        # Save datasets
        curated_path = self.output_dir / "curated_metadata.csv"
        rejected_path = self.output_dir / "rejected_strains.csv"
        
        curated_df.to_csv(curated_path, index=False)
        rejected_df.to_csv(rejected_path, index=False)
        
        logger.info(f"Saved curated dataset: {curated_path} ({passed_count} strains)")
        logger.info(f"Saved rejected dataset: {rejected_path} ({rejected_count} strains)")
        
        return curated_df, rejected_df

    def _fetch_biosample_metadata(self, biosample_id: str) -> Optional[dict]:
        """
        Fetch metadata from NCBI BioSample via efetch XML API.

        Args:
            biosample_id: NCBI BioSample accession (e.g., SAMN02604091)

        Returns:
            Dictionary of metadata attributes or None if fetch fails
        """
        try:
            # Construct efetch URL
            fetch_url = (
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=biosample&id={biosample_id}&retmode=xml&email={self.email}"
            )
            if self.api_key:
                fetch_url += f"&api_key={self.api_key}"
            
            # Fetch and parse XML
            with urllib.request.urlopen(fetch_url, timeout=10) as response:
                tree = ET.parse(response)
                root = tree.getroot()
            
            # Extract attributes from the BioSample XML
            attributes = {}
            for attr in root.findall(".//Attribute"):
                attr_name = attr.get("attribute_name", "").lower()
                attributes[attr_name] = attr.text
            
            # Extract key epidemiological fields
            metadata = {
                'geo_loc_name': attributes.get("geo_loc_name", "Unknown"),
                'collection_date': attributes.get("collection_date", "Unknown"),
                'host': attributes.get("host", "Unknown"),
                'isolation_source': attributes.get("isolation_source", "Unknown"),
            }
            
            # NCBI API pacing (STRICT CONSTRAINT)
            time.sleep(0.15)
            
            logger.debug(f"Successfully fetched metadata for {biosample_id}")
            return metadata
            
        except urllib.error.URLError as e:
            logger.error(f"Network error fetching {biosample_id}: {e}")
            return None
        except urllib.error.HTTPError as e:
            logger.error(f"HTTP error fetching {biosample_id}: {e}")
            return None
        except ET.ParseError as e:
            logger.error(f"XML parsing error for {biosample_id}: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error fetching {biosample_id}: {e}")
            return None

    def _validate_strain(self, location: str, date_str: str) -> Tuple[bool, str]:
        """
        Apply scientific filters: Geographic + Temporal.

        Args:
            location: Geographic location string
            date_str: Collection date string

        Returns:
            Tuple of (Passed: bool, Reason: str)
        """
        # 1. Geographic Filter (if configured)
        if self.filter_country and self.filter_country not in str(location).lower():
            return False, f"Failed Geographic Filter (Not {self.filter_country.title()})"
        
        # 2. Temporal Filter (if configured to a non-zero year)
        if self.filter_min_year is not None and self.filter_min_year > 0:
            if pd.isna(date_str) or str(date_str).lower() in ["n/a", "missing", "unknown"]:
                return False, "Failed Temporal Filter (Missing Date)"

            try:
                year = int(str(date_str)[:4])
                if year == 1905:
                    pass  # Rescue known THSTI typo
                elif year < self.filter_min_year:
                    return False, f"Failed Temporal Filter (Year {year} < {self.filter_min_year})"
            except (ValueError, TypeError):
                return False, "Failed Temporal Filter (Invalid Date Format)"
        
        # All filters passed
        return True, "Passed"
