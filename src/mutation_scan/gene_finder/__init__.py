"""
Gene Finder Module

Hybrid gene detection using ABRicate (resistance genes) and BLASTn (housekeeping genes).
Returns standardized DataFrames for downstream processing.
"""

from .gene_finder import GeneFinder

__all__ = ["GeneFinder"]

__all__ = ["AbricateWrapper", "BLASTWrapper"]
