"""
Helper module to create housekeeping gene reference databases.

This provides a placeholder function to create a small reference FASTA
containing common bacterial housekeeping genes (gyrA, parC, rpoB, etc.)
"""

import logging
from pathlib import Path
from typing import List, Dict

logger = logging.getLogger(__name__)


def create_housekeeping_reference(output_path: Path) -> Path:
    """
    Create a placeholder housekeeping gene reference database.
    
    This creates a small FASTA file with common bacterial housekeeping genes
    that can be used as a reference for BLASTn searches.

    Args:
        output_path: Path where reference FASTA will be created

    Returns:
        Path to created reference file

    Note:
        This is a placeholder implementation. For production use, download
        actual housekeeping gene sequences from NCBI or use curated databases.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Placeholder sequences (these are NOT real sequences - just examples)
    # In production, these should be downloaded from NCBI or a curated source
    placeholder_genes = {
        "gyrA": "ATGGCAACGTTCGACGATCTGACGCATGACTTCGACCATCGTACGCATGACCTGACGCAT",
        "parC": "ATGAGCGACGTTACGCATGACTTCGATCATCGTACGCATGACCTGACGCATGACTTCGAT",
        "rpoB": "ATGCGTACGCATGACTTCGATCATCGTACGCATGACCTGACGCATGACTTCGATCATCGT",
        "recA": "ATGGATCATCGTACGCATGACCTGACGCATGACTTCGATCATCGTACGCATGACCTGACG",
        "dnaA": "ATGACGCATGACTTCGATCATCGTACGCATGACCTGACGCATGACTTCGATCATCGTACG",
    }
    
    with open(output_path, 'w') as f:
        for gene_name, sequence in placeholder_genes.items():
            f.write(f">{gene_name}\n")
            # Write sequence in 60-character lines (FASTA standard)
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')
    
    logger.info(f"Created placeholder housekeeping reference: {output_path}")
    logger.warning(
        "WARNING: This reference contains placeholder sequences. "
        "For production use, download real sequences from NCBI."
    )
    
    return output_path


def download_housekeeping_genes(gene_names: List[str], output_path: Path) -> Path:
    """
    Download housekeeping gene sequences from NCBI (placeholder).
    
    Args:
        gene_names: List of gene names to download
        output_path: Path where FASTA will be saved

    Returns:
        Path to downloaded FASTA file

    Note:
        This is a placeholder. Actual implementation would use Bio.Entrez
        to fetch sequences from NCBI GenBank.
    """
    logger.warning(
        "download_housekeeping_genes() is not implemented. "
        "Use create_housekeeping_reference() for placeholder data or "
        "implement NCBI Entrez fetching manually."
    )
    
    # For now, just create placeholder
    return create_housekeeping_reference(output_path)


def get_default_housekeeping_genes() -> Dict[str, str]:
    """
    Get default list of housekeeping genes for bacterial analysis.

    Returns:
        Dictionary of gene names and their descriptions
    """
    return {
        "gyrA": "DNA gyrase subunit A (quinolone resistance)",
        "parC": "DNA topoisomerase IV subunit A (quinolone resistance)",
        "rpoB": "RNA polymerase beta subunit (rifampicin resistance)",
        "recA": "Recombinase A (DNA repair)",
        "dnaA": "Chromosomal replication initiator protein",
        "gapA": "Glyceraldehyde-3-phosphate dehydrogenase",
        "16S": "16S ribosomal RNA gene",
    }
