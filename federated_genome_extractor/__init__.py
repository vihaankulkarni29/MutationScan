"""
Federated Genome Extractor

A comprehensive Python package for downloading and extracting genome data
from multiple bioinformatics databases including NCBI, BV-BRC, EnteroBase, and PATRIC.

This package provides unified access to federated genome databases with
standardized metadata extraction and harmonization capabilities.

Key Features:
- Multi-database genome extraction
- Standardized metadata harmonization  
- Parallel download capabilities
- AMR (Antimicrobial Resistance) data integration
- Quality scoring and filtering
- CLI and Python API access

Usage:
    from federated_genome_extractor.clients.ncbi_client import NCBIClient
    
    client = NCBIClient()
    genomes = client.fetch_genomes("Escherichia coli", max_results=10)

Author: Vihaan Kulkarni <vihaankulkarni29@gmail.com>
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Vihaan Kulkarni"
__email__ = "vihaankulkarni29@gmail.com"

# Import main classes for easy access
try:
    from .clients.ncbi_client import NCBIClient
    from .clients.bvbrc_client import BVBRCClient
    from .harmonizer import harmonize_data
    
    __all__ = [
        'NCBIClient',
        'BVBRCClient', 
        'harmonize_data',
        '__version__',
        '__author__',
        '__email__'
    ]
except ImportError:
    # Handle case where optional dependencies might not be available
    __all__ = ['__version__', '__author__', '__email__']