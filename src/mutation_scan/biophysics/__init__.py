"""
Biophysics Module - 3D Structure Analysis and Binding Predictions

This module connects population-level genomics analysis to molecular dynamics
and docking simulations. It bridges the MutationScan framework with external
biophysics engines (AutoScan, AutoDock Vina, OpenMM) to predict the structural
impact of identified mutations.

Contains:
- autoscan_bridge: Integration with AutoScan for ΔΔG calculations
"""

from .autoscan_bridge import AutoScanBridge

__all__ = ["AutoScanBridge"]
