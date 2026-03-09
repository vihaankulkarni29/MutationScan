"""
Core Module - Genome Ingestion & Gene Detection

Contains modules for:
- ingestion_engine: Flexible identifier resolution engine (BioSample/Accession/BV-BRC)
- metadata_interrogator: Epidemiological metadata enrichment and scientific filtering
- universal_downloader: Multi-source genome assembly downloader (BV-BRC + NCBI Datasets)
- tblastn_extractor: Translating aligner for frameshift-free protein extraction
"""

from .tblastn_extractor import TblastnSequenceExtractor
from .ingestion_engine import GenomicIngestionEngine
from .metadata_interrogator import MetadataInterrogator
from .universal_downloader import UniversalGenomeDownloader
from .genome_processor import GenomeProcessor
from .gene_finder import GeneFinder
from .sequence_extractor import SequenceExtractor
from .coordinate_parser import CoordinateParser
from .translator import SequenceTranslator
from .reference_builder import ReferenceBuilder

__all__ = [
    "TblastnSequenceExtractor",
    "GenomicIngestionEngine",
    "MetadataInterrogator",
    "UniversalGenomeDownloader",
    "GenomeProcessor",
    "GeneFinder",
    "SequenceExtractor",
    "CoordinateParser",
    "SequenceTranslator",
    "ReferenceBuilder",
]
