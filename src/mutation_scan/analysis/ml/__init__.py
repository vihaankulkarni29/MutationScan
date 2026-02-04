"""
ML Submodule - Machine Learning Pipeline

Contains modules for:
- inference: ML model loading and prediction (ResistancePredictor)
- features: Feature extraction from sequences (BiophysicalEncoder)
- trainer: Model training pipeline (ModelZooTrainer - renamed from train_zoo.py)
- dataset: Data pipeline for ML (AMRDataPipeline - renamed from data_pipeline.py)
- benchmark: Model evaluation and benchmarking (ModelBenchmark)
"""

from .inference import ResistancePredictor
from .features import BiophysicalEncoder
from .trainer import ModelZooTrainer
from .dataset import AMRDataPipeline
from .benchmark import ModelBenchmark

__all__ = [
    "ResistancePredictor",
    "BiophysicalEncoder",
    "ModelZooTrainer",
    "AMRDataPipeline",
    "ModelBenchmark",
]
