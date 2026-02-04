# ML Models Directory

This directory contains trained machine learning models for Module 6 (ML Predictor).

## Required Files

Place your trained model files here with the following naming convention:

```
{antibiotic_name}_predictor.pkl
{antibiotic_name}_scaler.pkl (optional)
```

### Example Structure

```
models/
├── ciprofloxacin_predictor.pkl
├── ciprofloxacin_scaler.pkl
├── tetracycline_predictor.pkl
├── tetracycline_scaler.pkl
└── README.md (this file)
```

## File Format

- **Model files**: Scikit-learn models saved with `joblib.dump()`
- **Scaler files**: StandardScaler or similar preprocessing objects (if used during training)

## Training Your Own Models

To train models for new antibiotics:

1. Use the training pipeline from Module 6
2. Save the trained RandomForestClassifier with:
   ```python
   import joblib
   joblib.dump(model, "models/ciprofloxacin_predictor.pkl")
   ```
3. If you used feature scaling, save the scaler:
   ```python
   joblib.dump(scaler, "models/ciprofloxacin_scaler.pkl")
   ```

## Supported Antibiotics

The predictor automatically detects available models by scanning this directory for `*_predictor.pkl` files.

To check available models at runtime:
```python
from mutation_scan.analysis.ml import ResistancePredictor

predictor = ResistancePredictor(model_dir="models")
print(predictor.get_available_models())
```

## Model Requirements

Each model file should contain a trained scikit-learn classifier that:
- Accepts 5 numerical features (from BiophysicalEncoder)
- Returns class probabilities with `.predict_proba()`
- Classes: [0=Susceptible, 1=Resistant]

## Notes

- Models are lazy-loaded (only loaded when needed)
- Models are cached after first load
- The ML predictor gracefully handles missing models
- If no models exist, the pipeline falls back to "VUS" (Variant of Unknown Significance)
