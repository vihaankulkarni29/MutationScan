"""
Module 6: ML Predictor package.

This package provides machine learning-based resistance prediction for novel mutations.
Key components:
- ResistancePredictor: Main inference engine with lazy model loading
- BiophysicalEncoder: Feature engineering for mutation strings

Usage:
    from mutation_scan.ml_predictor import ResistancePredictor
    
    predictor = ResistancePredictor(model_dir="models")
    result = predictor.predict("S83L", antibiotic="Ciprofloxacin")
"""

from .inference import ResistancePredictor
from .features import BiophysicalEncoder

__all__ = ["ResistancePredictor", "BiophysicalEncoder"]
