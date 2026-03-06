"""
Core Module - Genome Ingestion & Gene Detection

Contains modules for:
- clinical_ingestion: BV-BRC metadata curation and automated FTP downloads
- tblastn_extractor: Translating aligner for frameshift-free protein extraction
- entrez_handler_legacy: Deprecated NCBI Datasets API (Legacy Fallback)
"""

from .entrez_handler_legacy import NCBIDatasetsGenomeDownloader, BulkGenomeDownloader
from .clinical_ingestion import ClinicalMetadataCurator
from .tblastn_extractor import TblastnSequenceExtractor
from .genome_processor import GenomeProcessor
from .gene_finder import GeneFinder
from .abricate_wrapper import AbricateWrapper
from .blast_wrapper import BLASTWrapper
from .sequence_extractor import SequenceExtractor
from .coordinate_parser import CoordinateParser
from .translator import SequenceTranslator
from .reference_builder import ReferenceBuilder

__all__ = [
    "NCBIDatasetsGenomeDownloader",
    "BulkGenomeDownloader",
    "ClinicalMetadataCurator",
    "TblastnSequenceExtractor",
    "GenomeProcessor",
    "GeneFinder",
    "AbricateWrapper",
    "BLASTWrapper",
    "SequenceExtractor",
    "CoordinateParser",
    "SequenceTranslator",
    "ReferenceBuilder",
]
