"""
SubScan Analyzer Package

A high-performance mutation detection and co-occurrence analysis engine
for the SubScan bioinformatics pipeline.
"""

__version__ = "1.0.0"
__author__ = "SubScan Development Team"
__email__ = "subscan@genomics.local"

from .engine import MutationAnalyzer

__all__ = ["MutationAnalyzer"]
