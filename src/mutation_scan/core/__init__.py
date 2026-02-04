"""
Core Module - Genome Ingestion & Gene Detection

Contains modules for:
- genome_extractor: NCBI Datasets API v2 genome downloads
- gene_finder: ABRicate + BLAST hybrid gene detection
- sequence_extractor: Coordinate parsing and sequence extraction
"""

from .entrez_handler import NCBIDatasetsGenomeDownloader
from .genome_processor import GenomeProcessor
from .gene_finder import GeneFinder
from .abricate_wrapper import AbricateWrapper
from .blast_wrapper import BLASTWrapper
from .sequence_extractor import SequenceExtractor
from .coordinate_parser import CoordinateParser
from .translator import SequenceTranslator

__all__ = [
    "NCBIDatasetsGenomeDownloader",
    "GenomeProcessor",
    "GeneFinder",
    "AbricateWrapper",
    "BLASTWrapper",
    "SequenceExtractor",
    "CoordinateParser",
    "SequenceTranslator",
]
