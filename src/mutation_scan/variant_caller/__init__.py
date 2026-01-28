"""
Variant Caller Module

Performs pairwise sequence alignments and identifies variants/mutations
between sequences.
"""

from .alignment import SequenceAligner
from .variant_detector import VariantDetector

__all__ = ["SequenceAligner", "VariantDetector"]
