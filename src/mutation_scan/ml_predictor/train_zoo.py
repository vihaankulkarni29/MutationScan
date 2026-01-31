"""
Model Zoo Trainer for Antimicrobial Resistance Prediction

Trains specialized Random Forest classifiers for each antibiotic using
a "Novel Mutation" evaluation strategy to measure generalization to unseen variants.
"""

from pathlib import Path
from typing import Dict, Tuple, List
import warnings

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import accuracy_score, roc_auc_score
import joblib

warnings.filterwarnings('ignore')


class ModelZooTrainer:
    """
    Trains and manages a zoo of antibiotic-specific resistance predictors.
    
    Key Innovation: Uses GroupKFold with mutation as grouping variable to simulate
    real-world deployment scenario where novel mutations (never seen during training)
    must be predicted at runtime.
    """
    
    # Feature set derived from biophysical encoding
    FEATURE_COLS = [
        'delta_hydrophobicity',
        'delta_charge',
        'delta_mw',
        'is_aromatic_change',
        'is_proline_change'
    ]
    
    MIN_SAMPLES_PER_DRUG = 50  # Filter out rare antibiotics
    
    def __init__(self, data_path: str = "data/processed_features.csv",
                 model_dir: str = "models",
                 report_path: str = "models/zoo_performance.csv"):
        """
        Initialize the Model Zoo Trainer.
        
        Args:
            data_path: Path to processed features CSV
            model_dir: Directory to save models
            report_path: Path to save performance report
        """
        self.data_path = Path(data_path)
        self.model_dir = Path(model_dir)
        self.report_path = Path(report_path)
        
        self.df = None
        self.antibiotics = None
        self.cv_results = []
        self.models = {}
        
        # Create models directory
        self.model_dir.mkdir(parents=True, exist_ok=True)
    
    def load_data(self) -> pd.DataFrame:
        """
        Load processed features from CSV.
        
        Returns:
            Loaded dataframe with features and labels
        
        Raises:
            FileNotFoundError: If data file doesn't exist
        """
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data file not found: {self.data_path}")
        
        self.df = pd.read_csv(self.data_path)
        print(f"✓ Loaded {len(self.df)} samples from {self.data_path}")
        print(f"  Columns: {list(self.df.columns)}")
        
        return self.df
    
    def select_antibiotics(self) -> List[str]:
        """
        Identify unique antibiotics with sufficient data.
        
        Filters out rare antibiotics with <50 samples to ensure robust
        model training and evaluation.
        
        Returns:
            List of antibiotic names passing the filter
        """
        antibiotic_counts = self.df['antibiotic'].value_counts()
        
        # Filter antibiotics with minimum sample count
        self.antibiotics = antibiotic_counts[
            antibiotic_counts >= self.MIN_SAMPLES_PER_DRUG
        ].index.tolist()
        
        print(f"\n📊 Antibiotic Filtering:")
        print(f"  Total unique antibiotics: {len(antibiotic_counts)}")
        print(f"  Antibiotics with ≥{self.MIN_SAMPLES_PER_DRUG} samples: {len(self.antibiotics)}")
        print(f"  Selected: {', '.join(self.antibiotics)}")
        
        return self.antibiotics
    
    def _create_mutation_groups(self, df: pd.DataFrame) -> np.ndarray:
        """
        Create mutation group identifiers for GroupKFold.
        
        Strategy: Combine WT, position, and mutant amino acid into a unique
        mutation identifier. Each unique mutation is treated as a group,
        ensuring the model never sees the same mutation in both training and
        validation sets.
        
        Args:
            df: Dataframe with 'wt', 'position', 'mutant' columns
        
        Returns:
            Array of group labels (one per sample)
        """
        mutations = (df['wt'].astype(str) + 
                    df['position'].astype(str) + 
                    df['mutant'].astype(str))
        
        # Convert to categorical codes for GroupKFold
        mutation_groups = pd.Categorical(mutations).codes
        return mutation_groups
    
    def train_antibiotic_model(self, antibiotic: str) -> Dict[str, float]:
        """
        Train and cross-validate a model for a specific antibiotic.
        
        Uses GroupKFold to simulate real-world scenario:
        - Training set: Models learn from known mutations
        - Validation set: Evaluated on novel (unseen) mutations
        
        This measures the model's ability to predict resistance for new genetic
        variants not encountered during training—critical for clinical utility.
        
        Args:
            antibiotic: Name of the antibiotic
        
        Returns:
            Dictionary with cross-validation metrics:
            - accuracy_mean, accuracy_std
            - roc_auc_mean, roc_auc_std
            - n_mutations (unique mutations in dataset)
            - n_samples (total samples)
        """
        # Filter data for this antibiotic
        df_drug = self.df[self.df['antibiotic'] == antibiotic].copy()
        
        X = df_drug[self.FEATURE_COLS].values
        y = df_drug['phenotype'].values
        
        # Create mutation groups
        groups = self._create_mutation_groups(df_drug)
        n_mutations = len(np.unique(groups))
        
        print(f"\n  Training {antibiotic}...")
        print(f"    Samples: {len(df_drug)}, Unique mutations: {n_mutations}")
        
        # Initialize GroupKFold: 5 folds, shuffle mutations randomly
        gkf = GroupKFold(n_splits=5)
        
        accuracies = []
        roc_aucs = []
        fold = 0
        
        # Cross-validation loop
        for train_idx, val_idx in gkf.split(X, y, groups):
            fold += 1
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            # Train model on training set
            model = RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                class_weight='balanced',
                random_state=42,
                n_jobs=-1
            )
            model.fit(X_train, y_train)
            
            # Evaluate on validation set (novel mutations)
            y_pred = model.predict(X_val)
            y_pred_proba = model.predict_proba(X_val)[:, 1]
            
            acc = accuracy_score(y_val, y_pred)
            roc = roc_auc_score(y_val, y_pred_proba)
            
            accuracies.append(acc)
            roc_aucs.append(roc)
            
            print(f"    Fold {fold}: Acc={acc:.3f}, ROC-AUC={roc:.3f}")
        
        # Compute mean and std
        acc_mean = np.mean(accuracies)
        acc_std = np.std(accuracies)
        roc_mean = np.mean(roc_aucs)
        roc_std = np.std(roc_aucs)
        
        print(f"    ✓ Mean Accuracy: {acc_mean:.3f} ± {acc_std:.3f}")
        print(f"    ✓ Mean ROC-AUC: {roc_mean:.3f} ± {roc_std:.3f}")
        
        return {
            'antibiotic': antibiotic,
            'n_samples': len(df_drug),
            'n_mutations': n_mutations,
            'accuracy_mean': acc_mean,
            'accuracy_std': acc_std,
            'roc_auc_mean': roc_mean,
            'roc_auc_std': roc_std,
        }
    
    def train_final_model(self, antibiotic: str) -> RandomForestClassifier:
        """
        Train final model on all data for this antibiotic.
        
        After cross-validation assessment, train on complete dataset for
        production deployment.
        
        Args:
            antibiotic: Name of the antibiotic
        
        Returns:
            Trained RandomForestClassifier
        """
        df_drug = self.df[self.df['antibiotic'] == antibiotic].copy()
        
        X = df_drug[self.FEATURE_COLS].values
        y = df_drug['phenotype'].values
        
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            class_weight='balanced',
            random_state=42,
            n_jobs=-1
        )
        
        model.fit(X, y)
        self.models[antibiotic] = model
        
        return model
    
    def save_model(self, antibiotic: str, model: RandomForestClassifier) -> Path:
        """
        Save trained model to disk using joblib.
        
        Args:
            antibiotic: Name of the antibiotic
            model: Trained RandomForestClassifier
        
        Returns:
            Path to saved model file
        """
        # Sanitize antibiotic name for filename
        safe_name = antibiotic.replace(' ', '_').replace('/', '_').lower()
        model_path = self.model_dir / f"{safe_name}_predictor.pkl"
        
        joblib.dump(model, model_path)
        print(f"    ✓ Saved model to {model_path}")
        
        return model_path
    
    def save_report(self, results: List[Dict]) -> Path:
        """
        Save cross-validation results to CSV report.
        
        Args:
            results: List of performance dictionaries
        
        Returns:
            Path to saved report
        """
        report_df = pd.DataFrame(results)
        report_df = report_df.sort_values('roc_auc_mean', ascending=False)
        
        self.report_path.parent.mkdir(parents=True, exist_ok=True)
        report_df.to_csv(self.report_path, index=False)
        
        print(f"\n✅ Saved performance report to {self.report_path}")
        print(f"\nModel Zoo Performance Summary:")
        print(report_df.to_string(index=False))
        
        return self.report_path
    
    def train_zoo(self) -> Tuple[List[Dict], Dict[str, RandomForestClassifier]]:
        """
        Train the complete model zoo.
        
        Workflow:
        1. Load data
        2. Identify antibiotics with sufficient samples
        3. For each antibiotic:
           - Cross-validate with GroupKFold (novel mutations strategy)
           - Train final model on all data
           - Save model and results
        
        Returns:
            Tuple of (cv_results_list, models_dict)
        """
        print("="*70)
        print("MODEL ZOO TRAINER - ANTIMICROBIAL RESISTANCE PREDICTION")
        print("="*70)
        
        # Load and prepare data
        self.load_data()
        self.select_antibiotics()
        
        print("\n🚀 Starting Model Zoo Training...")
        print(f"{'='*70}")
        
        # Train model for each antibiotic
        for antibiotic in self.antibiotics:
            # Cross-validation
            cv_result = self.train_antibiotic_model(antibiotic)
            self.cv_results.append(cv_result)
            
            # Train final model
            final_model = self.train_final_model(antibiotic)
            
            # Save model
            self.save_model(antibiotic, final_model)
        
        # Save report
        self.save_report(self.cv_results)
        
        print(f"\n{'='*70}")
        print(f"✅ MODEL ZOO TRAINING COMPLETE")
        print(f"{'='*70}")
        print(f"Total models trained: {len(self.models)}")
        print(f"Models saved to: {self.model_dir}/")
        print(f"Report saved to: {self.report_path}")
        
        return self.cv_results, self.models


def main():
    """Main execution."""
    trainer = ModelZooTrainer(
        data_path="data/processed_features.csv",
        model_dir="models",
        report_path="models/zoo_performance.csv"
    )
    
    cv_results, models = trainer.train_zoo()
    
    print(f"\n📦 Model Zoo Details:")
    print(f"  Antibiotics: {list(models.keys())}")
    print(f"  Total models: {len(models)}")


if __name__ == "__main__":
    main()
