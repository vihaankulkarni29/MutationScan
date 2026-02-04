"""
Analysis Module - Variant Calling & Machine Learning

Contains modules for:
- variant_caller: Sequence alignment and mutation detection (VariantCaller)
- variant_detector: Mutation detection logic (VariantDetector)
- alignment: Global alignment engine (SequenceAligner)
- ml: Machine learning inference, training, and feature extraction
"""

from .variant_caller import VariantCaller
from .variant_detector import VariantDetector
from .alignment import SequenceAligner

__all__ = [
    "VariantCaller",
    "VariantDetector",
    "SequenceAligner",
]
