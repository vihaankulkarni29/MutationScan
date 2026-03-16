"""
Core Module - Gene Detection & Sequence Analysis

Contains modules for:
- tblastn_extractor: Translating aligner for frameshift-free protein extraction
- genome_processor: FASTA genome parsing and contig management
- gene_finder: Resistance gene screening logic
- sequence_extractor: Coordinate-based sequence extraction
- coordinate_parser: BLAST HSP coordinate parsing
- translator: DNA codon translation
- reference_builder: Reference sequence index construction
"""

from .tblastn_extractor import TblastnSequenceExtractor
from .genome_processor import GenomeProcessor
from .gene_finder import GeneFinder
from .sequence_extractor import SequenceExtractor
from .coordinate_parser import CoordinateParser
from .translator import SequenceTranslator
from .reference_builder import ReferenceBuilder

__all__ = [
    "TblastnSequenceExtractor",
    "GenomeProcessor",
    "GeneFinder",
    "SequenceExtractor",
    "CoordinateParser",
    "SequenceTranslator",
    "ReferenceBuilder",
]
