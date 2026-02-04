"""
MutationScan: A Python-based bioinformatics pipeline for AMR gene detection and analysis.

This package provides automated analysis of bacterial genomes, including:
- Genome downloading and processing
- AMR gene screening
- Protein sequence extraction
- Pairwise sequence alignment
- 3D mutation visualization

Refactored v2.0: Modular src-layout architecture
- core: Genome ingestion & gene detection
- analysis: Variant calling & ML predictions
- visualization: 3D structure visualization
- utils: Shared utilities
"""

__version__ = "2.0.0"
__author__ = "MutationScan Team"
__email__ = "contact@example.com"

# Import subpackages
from . import core
from . import analysis
from . import visualization
from . import utils

__all__ = [
    "core",
    "analysis",
    "visualization",
    "utils",
]
