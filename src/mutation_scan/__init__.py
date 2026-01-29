"""
MutationScan: A Python-based bioinformatics pipeline for AMR gene detection and analysis.

This package provides automated analysis of bacterial genomes, including:
- Genome downloading and processing
- AMR gene screening
- Protein sequence extraction
- Pairwise sequence alignment
- 3D mutation visualization
"""

__version__ = "0.1.0"
__author__ = "Your Organization"
__email__ = "contact@example.com"

from . import genome_extractor
from . import gene_finder
from . import sequence_extractor
from . import variant_caller
from . import visualizer
from . import ml_predictor
from . import utils

__all__ = [
    "genome_extractor",
    "gene_finder",
    "sequence_extractor",
    "variant_caller",
    "visualizer",
    "ml_predictor",
    "utils",
]
