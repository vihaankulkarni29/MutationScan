"""
Genome Extractor Module

Handles downloading and preprocessing of bacterial genomes from NCBI using Entrez API.
Provides functionality for genome retrieval, validation, and storage.
"""

from .entrez_handler import EntrezGenomeDownloader
from .genome_processor import GenomeProcessor

__all__ = ["EntrezGenomeDownloader", "GenomeProcessor"]
