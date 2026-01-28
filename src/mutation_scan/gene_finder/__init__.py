"""
Gene Finder Module

Wrapper module for AMR gene detection using ABRicate and local BLASTn.
Provides unified interface for gene screening and annotation.
"""

from .abricate_wrapper import AbricateWrapper
from .blast_wrapper import BLASTWrapper

__all__ = ["AbricateWrapper", "BLASTWrapper"]
