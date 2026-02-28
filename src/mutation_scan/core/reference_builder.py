"""
Reference Builder Module - Automatic Wild-Type Protein Fetching

This module provides optional automatic fetching of wild-type reference proteins
from NCBI when user-supplied references are missing. User-provided references
always take priority over auto-fetched sequences.
"""

import logging
import time
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class ReferenceBuilder:
    """
    Automatically fetch wild-type reference proteins from NCBI.
    
    Features:
    - RefSeq-first search strategy for high-quality sequences
    - Fallback to non-RefSeq if RefSeq unavailable
    - Length filtering (100-2000 AA) to avoid fragments
    - Retry logic with exponential backoff
    """
    
    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize the ReferenceBuilder.
        
        Args:
            email: Email for NCBI (required by NCBI policy)
            api_key: Optional NCBI API key for higher rate limits
            
        Raises:
            ValueError: If email is not provided
        """
        if not email:
            raise ValueError("Email is required for NCBI API access")
        
        self.email = email
        self.api_key = api_key
        
        # Configure Entrez
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        
        logger.info(f"Initialized ReferenceBuilder for {email}")
    
    def fetch_reference(
        self,
        gene_name: str,
        organism: str,
        output_dir: Path
    ) -> bool:
        """
        Fetch wild-type reference protein from NCBI and save to FASTA.
        
        Search Strategy:
        1. Try RefSeq-filtered search (high quality)
        2. Fallback to non-RefSeq search if RefSeq fails
        3. Apply length filter (100-2000 AA) to avoid fragments
        
        Args:
            gene_name: Gene name (e.g., "oqxB")
            organism: Organism name (e.g., "Klebsiella pneumoniae")
            output_dir: Directory to save reference file
            
        Returns:
            True if reference was successfully fetched and saved, False otherwise
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = output_dir / f"{gene_name}_WT.faa"
        
        logger.info(f"Attempting to fetch reference for {gene_name} from {organism}")
        
        # Strategy 1: Try RefSeq first (highest quality)
        sequence = self._search_and_fetch(
            gene_name=gene_name,
            organism=organism,
            use_refseq=True
        )
        
        # Strategy 2: Fallback to non-RefSeq if RefSeq failed
        if not sequence:
            logger.warning(f"RefSeq search failed for {gene_name}, trying non-RefSeq...")
            sequence = self._search_and_fetch(
                gene_name=gene_name,
                organism=organism,
                use_refseq=False
            )
        
        # Save if successful
        if sequence:
            try:
                with open(output_file, 'w') as f:
                    SeqIO.write(sequence, f, "fasta")
                
                seq_length = len(sequence.seq) if sequence.seq else 0
                logging.info(f"Successfully saved reference: {output_file} ({seq_length} AA)")
                return True
                
            except Exception as e:
                logger.error(f"Failed to save reference for {gene_name}: {e}")
                return False
        else:
            logger.error(f"Could not fetch reference for {gene_name} from NCBI")
            return False
    
    def _search_and_fetch(
        self,
        gene_name: str,
        organism: str,
        use_refseq: bool = True
    ) -> Optional[SeqRecord]:
        """
        Search NCBI Protein database and fetch first valid sequence.
        
        Args:
            gene_name: Gene name to search
            organism: Organism name to filter
            use_refseq: Whether to restrict to RefSeq sequences
            
        Returns:
            SeqRecord if found, None otherwise
        """
        max_retries = 3
        retry_count = 0
        
        # Build search query
        query_parts = [
            f'"{gene_name}"[Gene]',
            f'"{organism}"[Organism]',
            "100:2000[Sequence Length]",  # Filter fragments and large proteins
        ]
        
        if use_refseq:
            query_parts.append('"RefSeq"[Filter]')
        
        query = " AND ".join(query_parts)
        
        while retry_count < max_retries:
            try:
                # Step 1: Search for protein IDs
                logger.debug(f"Searching NCBI Protein with query: {query}")
                
                search_handle = Entrez.esearch(
                    db="protein",
                    term=query,
                    retmax=5,  # Get top 5 results
                    sort="relevance"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                # Extract ID list from results (type: ignore for Entrez.read dict-like return)
                id_list = search_results.get("IdList", []) if isinstance(search_results, dict) else []
                
                if not id_list:
                    logger.debug(f"No results found for query: {query}")
                    return None
                
                logger.debug(f"Found {len(id_list)} protein IDs")
                
                # Step 2: Fetch the first valid protein sequence
                for protein_id in id_list:
                    try:
                        fetch_handle = Entrez.efetch(
                            db="protein",
                            id=protein_id,
                            rettype="fasta",
                            retmode="text"
                        )
                        
                        records = list(SeqIO.parse(fetch_handle, "fasta"))
                        fetch_handle.close()
                        
                        if records:
                            record = records[0]
                            
                            # Validate sequence length
                            if 100 <= len(record.seq) <= 2000:
                                # Clean up header for clarity
                                record.id = f"{gene_name}_WT"
                                record.description = f"Auto-fetched from NCBI | Original: {record.description}"
                                
                                logger.info(f"Fetched {gene_name} reference (ID: {protein_id}, {len(record.seq)} AA)")
                                return record
                            else:
                                logger.debug(f"Skipping protein {protein_id}: length {len(record.seq)} outside 100-2000 AA range")
                        
                    except Exception as e:
                        logger.debug(f"Failed to fetch protein {protein_id}: {e}")
                        continue
                
                # If we got here, no valid sequences were found
                return None
                
            except Exception as e:
                retry_count += 1
                logger.warning(f"Search attempt {retry_count} failed for {gene_name}: {e}")
                
                if retry_count < max_retries:
                    wait_time = 2 ** retry_count
                    logger.info(f"Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"Search failed after {max_retries} retries for {gene_name}")
                    return None
        
        return None
