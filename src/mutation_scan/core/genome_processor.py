"""
Genome processing and validation module.

Handles post-download validation, QC checking, and preparation of genome sequences
for downstream analysis.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Tuple

logger = logging.getLogger(__name__)


class GenomeProcessor:
    """Process, validate, and QC-check downloaded genomes."""

    def __init__(self, min_coverage: float = 90.0, min_length: int = 1000000):
        """
        Initialize the genome processor.

        Args:
            min_coverage: Minimum coverage threshold (%) for QC pass
            min_length: Minimum genome length threshold (bp)
        """
        self.min_coverage = min_coverage
        self.min_length = min_length
        logger.info(
            f"Initialized GenomeProcessor: min_coverage={min_coverage}%, "
            f"min_length={min_length}bp"
        )

    def validate_genome(self, filepath: Path) -> Tuple[bool, str]:
        """
        Validate genome FASTA format and length.

        Args:
            filepath: Path to genome FASTA file

        Returns:
            Tuple of (is_valid, status_message)
        """
        try:
            if not filepath.exists():
                return False, f"File not found: {filepath}"
            
            # Parse FASTA
            sequence_length = 0
            is_fasta = False
            
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        is_fasta = True
                    elif is_fasta and line:
                        sequence_length += len(line)
            
            if not is_fasta:
                return False, "Not a valid FASTA file (no header found)"
            
            if sequence_length < self.min_length:
                return False, f"Sequence too short: {sequence_length}bp < {self.min_length}bp"
            
            logger.debug(f"Validated {filepath.name}: {sequence_length}bp")
            return True, f"Valid FASTA: {sequence_length}bp"
            
        except Exception as e:
            logger.error(f"Validation error for {filepath}: {e}")
            return False, f"Validation error: {e}"

    def calculate_coverage(self, filepath: Path, reference_length: Optional[int] = None) -> float:
        """
        Calculate sequence coverage or simply verify presence.
        
        If no reference_length provided, returns 100% for valid files.

        Args:
            filepath: Path to genome FASTA file
            reference_length: Expected reference genome length (bp)

        Returns:
            Coverage percentage (0-100+)
        """
        try:
            sequence_length = 0
            with open(filepath, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        sequence_length += len(line.strip())
            
            if reference_length:
                coverage = (sequence_length / reference_length) * 100
            else:
                coverage = 100.0  # Default to 100% if no reference
            
            logger.debug(f"Coverage for {filepath.name}: {coverage:.1f}%")
            return coverage
            
        except Exception as e:
            logger.error(f"Coverage calculation error: {e}")
            return 0.0

    def extract_metadata(self, filepath: Path) -> Dict:
        """
        Extract basic metadata from genome FASTA header.

        Args:
            filepath: Path to genome FASTA file

        Returns:
            Dictionary with extracted metadata
        """
        metadata = {
            "filename": filepath.name,
            "file_size": filepath.stat().st_size,
            "sequences": 0,
            "first_header": "N/A",
        }
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        metadata["sequences"] += 1
                        if metadata["first_header"] == "N/A":
                            metadata["first_header"] = line.strip()
            
            logger.debug(f"Extracted metadata from {filepath.name}")
            return metadata
            
        except Exception as e:
            logger.error(f"Metadata extraction error: {e}")
            return metadata
