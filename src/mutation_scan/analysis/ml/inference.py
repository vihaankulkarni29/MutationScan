"""
Resistance Prediction Engine

Unified inference interface for antimicrobial resistance prediction.
Integrates the biophysical feature pipeline with trained model zoo.

Key Feature: Lazy loading of models to minimize memory footprint.
"""

from pathlib import Path
from typing import Dict, Optional
import warnings

import joblib

from .features import BiophysicalEncoder

warnings.filterwarnings('ignore')


class ResistancePredictor:
    """
    Production-ready inference engine for AMR prediction.
    
    Features:
    - Lazy model loading: Models cached only when requested
    - Integrated biophysical encoding pipeline
    - Structured output with risk stratification
    - Comprehensive error handling
    """
    
    def __init__(self, model_dir: str = "models"):
        """
        Initialize the predictor.
        
        Args:
            model_dir: Directory containing trained .pkl models
        """
        self.model_dir = Path(model_dir)
        self.encoder = BiophysicalEncoder()
        self.loaded_models = {}  # Cache for lazy-loaded models
        self._available_models = None
    
    def get_available_models(self) -> list:
        """
        Scan model directory for available antibiotic models.
        
        Returns:
            List of antibiotic names with trained models
        """
        if self._available_models is not None:
            return self._available_models
        
        if not self.model_dir.exists():
            raise FileNotFoundError(f"Models directory not found: {self.model_dir}")
        
        # Find all .pkl files and extract antibiotic names
        model_files = list(self.model_dir.glob("*_predictor.pkl"))
        antibiotics = [f.stem.replace('_predictor', '').replace('_', ' ').title() 
                      for f in model_files]
        
        self._available_models = sorted(antibiotics)
        return self._available_models
    
    def _load_model(self, antibiotic: str):
        """
        Lazy load a model for a specific antibiotic.
        
        Uses caching to avoid reloading the same model multiple times.
        
        Args:
            antibiotic: Name of the antibiotic
        
        Returns:
            Trained RandomForestClassifier model
        
        Raises:
            ValueError: If model doesn't exist
        """
        # Check cache
        if antibiotic in self.loaded_models:
            return self.loaded_models[antibiotic]
        
        # Construct model path
        safe_name = antibiotic.replace(' ', '_').lower()
        model_path = self.model_dir / f"{safe_name}_predictor.pkl"
        
        if not model_path.exists():
            available = self.get_available_models()
            raise ValueError(
                f"Model for '{antibiotic}' not found.\n"
                f"Available models: {', '.join(available)}"
            )
        
        # Load and cache model
        model = joblib.load(model_path)
        self.loaded_models[antibiotic] = model
        
        return model
    
    def _classify_risk(self, probability: float) -> str:
        """
        Stratify risk based on predicted resistance probability.
        
        Args:
            probability: Predicted probability of resistance [0, 1]
        
        Returns:
            Risk level: "High", "Moderate", or "Low"
        """
        if probability > 0.8:
            return "High"
        elif probability > 0.5:
            return "Moderate"
        else:
            return "Low"
    
    def predict(self, mutation: str, antibiotic: str) -> Dict:
        """
        Predict resistance for a given mutation-antibiotic pair.
        
        Args:
            mutation: Mutation string (e.g., "S83L", "gyrA_S83L")
            antibiotic: Antibiotic name (e.g., "Ciprofloxacin")
        
        Returns:
            Dictionary with prediction results:
            {
                'antibiotic': str,
                'mutation': str,
                'wt': str (wild-type amino acid),
                'position': int,
                'mutant': str (mutant amino acid),
                'resistance_prob': float (0-1),
                'risk_level': str ('High'|'Moderate'|'Low'),
                'success': bool,
            }
        
        Raises:
            ValueError: If mutation cannot be parsed or model not found
        """
        # Encode mutation
        features = self.encoder.get_features(mutation)
        
        if features is None:
            return {
                'antibiotic': antibiotic,
                'mutation': mutation,
                'resistance_prob': None,
                'risk_level': None,
                'success': False,
                'error': f"Could not parse mutation: {mutation}"
            }
        
        # Load model (lazy)
        try:
            model = self._load_model(antibiotic)
        except ValueError as e:
            return {
                'antibiotic': antibiotic,
                'mutation': mutation,
                'resistance_prob': None,
                'risk_level': None,
                'success': False,
                'error': str(e)
            }
        
        # Prepare feature vector
        feature_vector = self.encoder.get_features_as_array(mutation)
        
        if feature_vector is None:
            return {
                'antibiotic': antibiotic,
                'mutation': mutation,
                'resistance_prob': None,
                'risk_level': None,
                'success': False,
                'error': f"Could not encode mutation: {mutation}"
            }
        
        # Make prediction
        try:
            resistance_prob = model.predict_proba(feature_vector.reshape(1, -1))[0, 1]
            risk_level = self._classify_risk(resistance_prob)
            
            return {
                'antibiotic': antibiotic,
                'mutation': mutation,
                'wt': features['wt'],
                'position': features['position'],
                'mutant': features['mutant'],
                'resistance_prob': round(float(resistance_prob), 4),
                'risk_level': risk_level,
                'success': True,
            }
        except Exception as e:
            return {
                'antibiotic': antibiotic,
                'mutation': mutation,
                'resistance_prob': None,
                'risk_level': None,
                'success': False,
                'error': f"Prediction failed: {str(e)}"
            }
    
    def batch_predict(self, mutations: list, antibiotic: str) -> list:
        """
        Predict resistance for multiple mutations.
        
        Args:
            mutations: List of mutation strings
            antibiotic: Antibiotic name
        
        Returns:
            List of prediction dictionaries
        """
        results = []
        for mutation in mutations:
            result = self.predict(mutation, antibiotic)
            results.append(result)
        
        return results
