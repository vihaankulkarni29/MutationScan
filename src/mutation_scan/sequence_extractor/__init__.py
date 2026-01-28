"""
Sequence Extractor Module

Handles coordinate-based protein sequence extraction from genomic sequences.
Provides translation and sequence processing utilities.
"""

from .coordinate_parser import CoordinateParser
from .translator import SequenceTranslator

__all__ = ["CoordinateParser", "SequenceTranslator"]
