"""
Genome Extractor Module

Handles downloading and preprocessing of bacterial genomes from NCBI using modern
Datasets API v2. Provides batch download, metadata extraction, and genome validation.
"""

from .entrez_handler import NCBIDatasetsGenomeDownloader
from .genome_processor import GenomeProcessor

__all__ = ["NCBIDatasetsGenomeDownloader", "GenomeProcessor"]
