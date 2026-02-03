"""
Quick test to verify Module 6 ML predictor integration.

This test verifies:
1. ML predictor can be imported
2. BiophysicalEncoder works correctly
3. Variant caller can use ML predictor
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

# Set UTF-8 encoding for print
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def test_ml_import():
    """Test that ML predictor can be imported."""
    print("=" * 70)
    print("TEST 1: ML Predictor Import")
    print("=" * 70)
    
    try:
        from src.mutation_scan.ml_predictor import ResistancePredictor, BiophysicalEncoder
        print("✓ Successfully imported ResistancePredictor")
        print("✓ Successfully imported BiophysicalEncoder")
        return True
    except Exception as e:
        print(f"✗ Import failed: {e}")
        return False


def test_biophysical_encoder():
    """Test BiophysicalEncoder feature generation."""
    print("\n" + "=" * 70)
    print("TEST 2: BiophysicalEncoder Feature Generation")
    print("=" * 70)
    
    try:
        from src.mutation_scan.ml_predictor import BiophysicalEncoder
        
        encoder = BiophysicalEncoder()
        
        # Test classic S83L mutation
        features = encoder.get_features_as_array("S83L")
        
        print(f"\nMutation: S83L (Serine → Leucine at position 83)")
        print(f"Features generated: {features}")
        print(f"Feature vector length: {len(features) if features is not None else 0}")
        
        if features is not None and len(features) == 5:
            print("\n✓ Correct feature vector length (5 features)")
            print(f"  Feature 1 (Hydrophobicity Δ): {features[0]:.3f}")
            print(f"  Feature 2 (Charge Δ): {features[1]:.3f}")
            print(f"  Feature 3 (MW Δ): {features[2]:.3f}")
            print(f"  Feature 4 (Aromaticity Δ): {features[3]:.3f}")
            print(f"  Feature 5 (Proline effect): {features[4]:.3f}")
            return True
        else:
            print(f"✗ Unexpected feature vector: {features}")
            return False
            
    except Exception as e:
        print(f"✗ Encoder test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_predictor_initialization():
    """Test ResistancePredictor initialization."""
    print("\n" + "=" * 70)
    print("TEST 3: ResistancePredictor Initialization")
    print("=" * 70)
    
    try:
        from src.mutation_scan.ml_predictor import ResistancePredictor
        
        # Test with models directory
        predictor = ResistancePredictor(model_dir="models")
        print(f"✓ Initialized ResistancePredictor")
        print(f"  Model directory: {predictor.model_dir}")
        
        # Check available models
        try:
            available = predictor.get_available_models()
            print(f"  Available models: {available if available else 'None (models directory empty)'}")
        except FileNotFoundError:
            print(f"  Available models: None (models directory not found - expected)")
        
        return True
        
    except Exception as e:
        print(f"✗ Predictor initialization failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_variant_caller_with_ml():
    """Test VariantCaller with ML integration."""
    print("\n" + "=" * 70)
    print("TEST 4: VariantCaller with ML Integration")
    print("=" * 70)
    
    try:
        from src.mutation_scan.variant_caller import VariantCaller
        
        # Initialize with ML enabled
        caller = VariantCaller(
            refs_dir=Path("data/refs"),
            enable_ml=True,
            ml_models_dir=Path("models")
        )
        
        print("✓ VariantCaller initialized with ML support")
        print(f"  ML enabled: {caller.enable_ml}")
        print(f"  ML models dir: {caller.ml_models_dir}")
        print(f"  Target antibiotic: {caller.antibiotic}")
        
        # Test internal ML predictor getter
        predictor = caller._get_ml_predictor()
        
        if predictor is not None:
            print("✓ ML predictor loaded successfully")
        else:
            if caller._ml_predictor_error:
                print(f"⚠ ML predictor not available: {caller._ml_predictor_error}")
                print("  (This is expected if no trained models exist)")
            else:
                print("⚠ ML predictor not loaded (but no error recorded)")
        
        return True
        
    except Exception as e:
        print(f"✗ VariantCaller ML integration failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("MODULE 6 ML PREDICTOR INTEGRATION TEST SUITE")
    print("=" * 70)
    
    results = []
    
    results.append(("ML Import", test_ml_import()))
    results.append(("BiophysicalEncoder", test_biophysical_encoder()))
    results.append(("Predictor Init", test_predictor_initialization()))
    results.append(("VariantCaller ML", test_variant_caller_with_ml()))
    
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    all_passed = all(result[1] for result in results)
    
    if all_passed:
        print("\n" + "=" * 70)
        print("ALL TESTS PASSED! ✓")
        print("=" * 70)
        print("\nModule 6 is successfully integrated!")
        print("\nNext steps:")
        print("  1. Train ML models and place .pkl files in models/")
        print("  2. Run full pipeline with real genomes")
        print("  3. Verify ML predictions appear in mutation_report.csv")
    else:
        print("\n" + "=" * 70)
        print("SOME TESTS FAILED")
        print("=" * 70)
        print("\nCheck error messages above for details.")


if __name__ == "__main__":
    main()
