"""
Sequence Extraction module for MutationScan.

This module extracts DNA sequences from bacterial genomes based on gene coordinates
from GeneFinder, translates them to amino acids, and handles biological edge cases.

CRITICAL COORDINATE SYSTEM NOTES:
- BLAST/ABRicate coordinates are 1-based (inclusive start, inclusive end)
- Python slicing is 0-based (inclusive start, exclusive end)
- Formula: dna_slice = contig_seq[start-1 : end]
- This is the "Off-By-One" (OBO) conversion
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class SequenceExtractor:
    """
    Extract and translate gene sequences from bacterial genomes.
    
    This class takes GeneFinder output (DataFrame with coordinates) and extracts
    the corresponding sequences from genome FASTA files, handling strand orientation
    and bacterial translation properly.
    
    Key Features:
    - Efficient genome loading using Bio.SeqIO.index() (lazy loading)
    - Correct OBO coordinate conversion (1-based to 0-based)
    - Reverse complement for minus strand genes
    - Bacterial translation table (Table 11)
    - Graceful handling of partial genes (no cds=True)
    - Standardized output headers: >GeneName|Accession|Contig|Start-End
    """

    def __init__(self, genomes_dir: Path):
        """
        Initialize SequenceExtractor.

        Args:
            genomes_dir: Path to directory containing genome FASTA files
                        (output from GenomeExtractor module)

        Raises:
            FileNotFoundError: If genomes_dir doesn't exist
        """
        self.genomes_dir = Path(genomes_dir)
        
        if not self.genomes_dir.exists():
            raise FileNotFoundError(f"Genomes directory not found: {self.genomes_dir}")
        
        if not self.genomes_dir.is_dir():
            raise NotADirectoryError(f"Not a directory: {self.genomes_dir}")
        
        logger.info(f"Initialized SequenceExtractor with genomes_dir: {self.genomes_dir}")

    def extract_sequences(
        self, 
        genes_df: pd.DataFrame,
        accession: str,
        output_dir: Path,
        translate: bool = True
    ) -> Tuple[int, int]:
        """
        Extract sequences for all genes in DataFrame from a single genome.

        Args:
            genes_df: DataFrame from GeneFinder with columns:
                     Gene, Contig, Start, End, Strand, Identity, Source
            accession: Genome accession (e.g., 'GCF_000005845')
            output_dir: Output directory for .faa files
            translate: Whether to translate DNA to protein (default: True)

        Returns:
            Tuple of (successful_count, failed_count)

        Raises:
            FileNotFoundError: If genome FASTA file not found
        """
        # Locate genome FASTA file
        fasta_file = self.genomes_dir / f"{accession}.fasta"
        
        if not fasta_file.exists():
            logger.error(f"Genome FASTA not found: {fasta_file}")
            raise FileNotFoundError(f"Genome FASTA not found: {fasta_file}")
        
        # Ensure output directory exists
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load genome EFFICIENTLY using SeqIO.index() (lazy loading)
        logger.info(f"Loading genome: {fasta_file.name}")
        genome_idx = SeqIO.index(str(fasta_file), "fasta")
        
        successful = 0
        failed = 0
        
        logger.info(f"Processing {len(genes_df)} genes from {accession}")
        
        for idx, row in genes_df.iterrows():
            gene_name = row['Gene']
            contig_id = row['Contig']
            start = int(row['Start'])
            end = int(row['End'])
            strand = row['Strand']
            
            try:
                # Extract sequence
                if translate:
                    protein_seq = self._extract_and_translate(
                        genome_idx, contig_id, start, end, strand, gene_name
                    )
                    
                    if protein_seq:
                        # Write protein FASTA
                        output_file = output_dir / f"{accession}_{gene_name}.faa"
                        self._write_fasta(
                            output_file,
                            protein_seq,
                            gene_name,
                            accession,
                            contig_id,
                            start,
                            end
                        )
                        successful += 1
                        logger.debug(f"Extracted {gene_name}: {len(protein_seq)} aa")
                    else:
                        failed += 1
                else:
                    # Extract DNA only (no translation)
                    dna_seq = self._extract_dna(
                        genome_idx, contig_id, start, end, strand
                    )
                    
                    if dna_seq:
                        output_file = output_dir / f"{accession}_{gene_name}.fna"
                        self._write_fasta(
                            output_file,
                            dna_seq,
                            gene_name,
                            accession,
                            contig_id,
                            start,
                            end
                        )
                        successful += 1
                        logger.debug(f"Extracted {gene_name}: {len(dna_seq)} bp")
                    else:
                        failed += 1
                        
            except Exception as e:
                logger.error(f"Failed to extract {gene_name}: {e}")
                failed += 1
        
        genome_idx.close()
        
        logger.info(f"Extraction complete: {successful} succeeded, {failed} failed")
        return (successful, failed)

    def _extract_dna(
        self,
        genome_idx: Dict,
        contig_id: str,
        start: int,
        end: int,
        strand: str
    ) -> Optional[Seq]:
        """
        Extract DNA sequence from genome.

        Args:
            genome_idx: SeqIO.index() genome dictionary
            contig_id: Contig identifier
            start: Start position (1-based, inclusive)
            end: End position (1-based, inclusive)
            strand: '+' or '-'

        Returns:
            Bio.Seq object or None if extraction failed
        """
        # Check if contig exists in genome
        if contig_id not in genome_idx:
            logger.warning(
                f"Contig '{contig_id}' not found in genome. "
                f"Available contigs: {list(genome_idx.keys())[:5]}..."
            )
            return None
        
        try:
            # Get contig sequence
            contig_record = genome_idx[contig_id]
            contig_seq = contig_record.seq
            
            # CRITICAL: OBO conversion (1-based to 0-based)
            # BLAST/ABRicate: 1-based inclusive
            # Python slicing: 0-based start (inclusive), end (exclusive)
            # Formula: contig_seq[start-1 : end]
            dna_seq = contig_seq[start - 1 : end]
            
            # Handle reverse strand
            if strand == '-':
                dna_seq = dna_seq.reverse_complement()
            
            return dna_seq
            
        except Exception as e:
            logger.error(f"Error extracting DNA from {contig_id}:{start}-{end}: {e}")
            return None

    def _extract_and_translate(
        self,
        genome_idx: Dict,
        contig_id: str,
        start: int,
        end: int,
        strand: str,
        gene_name: str
    ) -> Optional[Seq]:
        """
        Extract DNA sequence and translate to protein.

        Args:
            genome_idx: SeqIO.index() genome dictionary
            contig_id: Contig identifier
            start: Start position (1-based, inclusive)
            end: End position (1-based, inclusive)
            strand: '+' or '-'
            gene_name: Gene name (for logging)

        Returns:
            Protein sequence (Bio.Seq) or None if translation failed
        """
        # Extract DNA
        dna_seq = self._extract_dna(genome_idx, contig_id, start, end, strand)
        
        if not dna_seq:
            return None
        
        try:
            # CRITICAL: Use Table 11 (Bacterial/Archaeal/Plant Plastid)
            # Do NOT use cds=True (crashes on partial genes without stop codons)
            # ABRicate finds fragments, so we translate raw sequences
            protein_seq = dna_seq.translate(table=11)
            
            # Trim stop codons (*) from the end if present
            # This is common for partial genes
            protein_str = str(protein_seq).rstrip('*')
            
            if not protein_str:
                logger.warning(f"Empty translation for {gene_name}")
                return None
            
            return Seq(protein_str)
            
        except Exception as e:
            logger.error(f"Translation failed for {gene_name}: {e}")
            return None

    def _write_fasta(
        self,
        output_file: Path,
        sequence: Seq,
        gene_name: str,
        accession: str,
        contig_id: str,
        start: int,
        end: int
    ) -> None:
        """
        Write sequence to FASTA file with standardized header.

        CRITICAL: Header format must be: >GeneName|Accession|Contig|Start-End
        This format ensures downstream tools (VariantCaller) can parse metadata.

        Args:
            output_file: Output FASTA file path
            sequence: Sequence to write (DNA or protein)
            gene_name: Gene name
            accession: Genome accession
            contig_id: Contig identifier
            start: Start position (1-based)
            end: End position (1-based)
        """
        # Standardized header: >GeneName|Accession|Contig|Start-End
        header = f"{gene_name}|{accession}|{contig_id}|{start}-{end}"
        
        record = SeqRecord(
            sequence,
            id=header,
            description=""  # Empty to avoid duplication in header
        )
        
        SeqIO.write(record, str(output_file), "fasta")
        logger.debug(f"Wrote {output_file.name}")

    def extract_all_genomes(
        self,
        genes_df: pd.DataFrame,
        output_dir: Path,
        translate: bool = True
    ) -> Dict[str, Tuple[int, int]]:
        """
        Extract sequences for all genomes in the DataFrame.

        This method groups genes_df by genome accession and processes each
        genome independently. Useful for batch processing multiple genomes.

        Args:
            genes_df: DataFrame from GeneFinder with ALL genes from multiple genomes
            output_dir: Output directory for .faa files
            translate: Whether to translate DNA to protein (default: True)

        Returns:
            Dictionary mapping accession -> (successful_count, failed_count)

        Note:
            Assumes gene coordinates are stored in a 'GeneCoordinates.csv' file
            that includes an 'Accession' column. If not present, you must call
            extract_sequences() for each genome individually.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        # Check if DataFrame has Accession column for batch processing
        if 'Accession' not in genes_df.columns:
            logger.error(
                "DataFrame missing 'Accession' column. "
                "Use extract_sequences() for individual genome processing."
            )
            raise ValueError("DataFrame must contain 'Accession' column for batch processing")
        
        # Group by accession
        grouped = genes_df.groupby('Accession')
        
        logger.info(f"Processing {len(grouped)} genomes")
        
        for accession, genome_genes in grouped:
            logger.info(f"Processing genome: {accession}")
            
            try:
                success, fail = self.extract_sequences(
                    genome_genes,
                    accession,
                    output_dir,
                    translate=translate
                )
                results[accession] = (success, fail)
                
            except Exception as e:
                logger.error(f"Failed to process {accession}: {e}")
                results[accession] = (0, len(genome_genes))
        
        return results

    def get_available_genomes(self) -> List[str]:
        """
        Get list of available genome accessions in genomes_dir.

        Returns:
            List of accession IDs (FASTA filenames without extension)
        """
        fasta_files = list(self.genomes_dir.glob("*.fasta")) + \
                     list(self.genomes_dir.glob("*.fna"))
        
        accessions = [f.stem for f in fasta_files]
        
        logger.info(f"Found {len(accessions)} genomes in {self.genomes_dir}")
        return accessions


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
if __name__ == "__main__":
    """
    Example usage of SequenceExtractor.
    
    Typical workflow:
    1. Run GenomeExtractor → Download genomes to data/genomes/
    2. Run GeneFinder → Generate genes.csv with coordinates
    3. Run SequenceExtractor → Extract and translate sequences
    """
    logging.basicConfig(level=logging.INFO)
    
    # Initialize extractor
    extractor = SequenceExtractor(
        genomes_dir=Path("data/genomes")
    )
    
    # Load GeneFinder output
    genes_df = pd.read_csv("data/genes/GCF_000005845_genes.csv")
    
    # Extract and translate sequences
    success, fail = extractor.extract_sequences(
        genes_df=genes_df,
        accession="GCF_000005845",
        output_dir=Path("data/proteins"),
        translate=True
    )
    
    print(f"Extraction complete: {success} succeeded, {fail} failed")
