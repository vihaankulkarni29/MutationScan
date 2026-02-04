"""
Benchmark: Biophysical vs Bag-of-Words Models

Validates the hypothesis that biophysical feature encoding outperforms
simple string-based ("one-hot") approaches when predicting novel mutations.

Key Test: Novel mutations (unseen during training) should be predictable
from biophysical properties but unrecognizable to word-based models.
"""

from pathlib import Path
from typing import Dict, Tuple, List
import warnings

import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')


class ModelBenchmark:
    """
    Head-to-head competition between Bag-of-Words and Biophysical models.
    
    Rationale:
    - BOW Model: Learns mutation NAME patterns (e.g., "S83L" as a token)
      → Novel mutations → Unknown tokens → Predicts random (~50%)
    - Biophysical Model: Learns PROPERTY patterns (hydrophobicity change, etc.)
      → Novel mutations → Still have recognizable properties → Predicts well (>65%)
    """
    
    FEATURE_COLS = [
        'delta_hydrophobicity',
        'delta_charge',
        'delta_mw',
        'is_aromatic_change',
        'is_proline_change'
    ]
    
    def __init__(self, data_path: str = "data/processed_features.csv",
                 output_dir: str = "models",
                 antibiotic: str = None):
        """
        Initialize benchmark.
        
        Args:
            data_path: Path to processed features CSV
            output_dir: Directory to save results
            antibiotic: Specific antibiotic to benchmark (auto-selected if None)
        """
        self.data_path = Path(data_path)
        self.output_dir = Path(output_dir)
        self.antibiotic = antibiotic
        
        self.df = None
        self.results = {}
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def load_data(self) -> pd.DataFrame:
        """Load processed features."""
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data file not found: {self.data_path}")
        
        self.df = pd.read_csv(self.data_path)
        print(f"✓ Loaded {len(self.df)} samples")
        
        return self.df
    
    def select_antibiotic(self) -> str:
        """Select top antibiotic by sample count."""
        if self.antibiotic:
            # Verify antibiotic exists
            if self.antibiotic not in self.df['antibiotic'].values:
                raise ValueError(f"Antibiotic '{self.antibiotic}' not found in data")
            return self.antibiotic
        
        # Auto-select top antibiotic
        top_ab = self.df['antibiotic'].value_counts().index[0]
        self.antibiotic = top_ab
        
        return top_ab
    
    def _create_mutation_groups(self, df: pd.DataFrame) -> np.ndarray:
        """Create mutation group identifiers for GroupKFold."""
        mutations = (df['wt'].astype(str) + 
                    df['position'].astype(str) + 
                    df['mutant'].astype(str))
        
        mutation_groups = pd.Categorical(mutations).codes
        return mutation_groups
    
    def _create_mutation_strings(self, df: pd.DataFrame) -> np.ndarray:
        """Create mutation strings for CountVectorizer (BOW model)."""
        mutations = (df['wt'].astype(str) + 
                    df['position'].astype(str) + 
                    df['mutant'].astype(str))
        
        return mutations.values
    
    def benchmark(self) -> Dict[str, Dict]:
        """
        Run benchmark comparing Bag-of-Words vs Biophysical models.
        
        Returns:
            Dictionary with results for each model
        """
        print("="*70)
        print("BENCHMARK: BAG-OF-WORDS vs BIOPHYSICAL MODELS")
        print("="*70)
        
        # Load and filter data
        self.load_data()
        antibiotic = self.select_antibiotic()
        df_ab = self.df[self.df['antibiotic'] == antibiotic].copy()
        
        print(f"\n🎯 Antibiotic: {antibiotic}")
        print(f"   Samples: {len(df_ab)}")
        print(f"   Unique mutations: {df_ab[['wt', 'position', 'mutant']].drop_duplicates().shape[0]}")
        
        # Prepare data
        X_biophysical = df_ab[self.FEATURE_COLS].values
        y = df_ab['phenotype'].values
        mutation_strings = self._create_mutation_strings(df_ab)
        groups = self._create_mutation_groups(df_ab)
        
        # Initialize GroupKFold
        gkf = GroupKFold(n_splits=5)
        
        # Storage for fold results
        bow_accuracies = []
        bow_roc_aucs = []
        bio_accuracies = []
        bio_roc_aucs = []
        
        print(f"\n🔄 Running 5-Fold Cross-Validation (Novel Mutation Strategy)...")
        print(f"{'-'*70}")
        
        fold = 0
        for train_idx, val_idx in gkf.split(X_biophysical, y, groups):
            fold += 1
            
            # Split data
            X_bio_train, X_bio_val = X_biophysical[train_idx], X_biophysical[val_idx]
            mut_train, mut_val = mutation_strings[train_idx], mutation_strings[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            # ========== MODEL A: BAG-OF-WORDS (BASELINE) ==========
            # Vectorize mutations as sparse word features
            vectorizer = CountVectorizer(analyzer='char', ngram_range=(1, 2), max_features=50)
            X_bow_train = vectorizer.fit_transform(mut_train)
            X_bow_val = vectorizer.transform(mut_val)
            
            # Train Logistic Regression (simple baseline)
            bow_model = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
            bow_model.fit(X_bow_train, y_train)
            
            # Evaluate BOW model
            bow_pred = bow_model.predict(X_bow_val)
            bow_pred_proba = bow_model.predict_proba(X_bow_val)[:, 1]
            
            bow_acc = accuracy_score(y_val, bow_pred)
            bow_roc = roc_auc_score(y_val, bow_pred_proba)
            
            bow_accuracies.append(bow_acc)
            bow_roc_aucs.append(bow_roc)
            
            # ========== MODEL B: BIOPHYSICAL (PROPOSED) ==========
            # Use physical feature vectors
            bio_model = RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                class_weight='balanced',
                random_state=42,
                n_jobs=-1
            )
            bio_model.fit(X_bio_train, y_train)
            
            # Evaluate Biophysical model
            bio_pred = bio_model.predict(X_bio_val)
            bio_pred_proba = bio_model.predict_proba(X_bio_val)[:, 1]
            
            bio_acc = accuracy_score(y_val, bio_pred)
            bio_roc = roc_auc_score(y_val, bio_pred_proba)
            
            bio_accuracies.append(bio_acc)
            bio_roc_aucs.append(bio_roc)
            
            # Print fold results
            print(f"Fold {fold}:")
            print(f"  BOW Model  → Acc: {bow_acc:.3f}, ROC-AUC: {bow_roc:.3f}")
            print(f"  Biophysical → Acc: {bio_acc:.3f}, ROC-AUC: {bio_roc:.3f}")
            print(f"  ✓ Improvement: Acc +{(bio_acc - bow_acc)*100:+.1f}%, ROC-AUC +{(bio_roc - bow_roc)*100:+.1f}%")
        
        # Compute aggregate results
        bow_acc_mean = np.mean(bow_accuracies)
        bow_acc_std = np.std(bow_accuracies)
        bow_roc_mean = np.mean(bow_roc_aucs)
        bow_roc_std = np.std(bow_roc_aucs)
        
        bio_acc_mean = np.mean(bio_accuracies)
        bio_acc_std = np.std(bio_accuracies)
        bio_roc_mean = np.mean(bio_roc_aucs)
        bio_roc_std = np.std(bio_roc_aucs)
        
        self.results = {
            'bag_of_words': {
                'accuracy_mean': bow_acc_mean,
                'accuracy_std': bow_acc_std,
                'roc_auc_mean': bow_roc_mean,
                'roc_auc_std': bow_roc_std,
                'accuracies': bow_accuracies,
                'roc_aucs': bow_roc_aucs,
            },
            'biophysical': {
                'accuracy_mean': bio_acc_mean,
                'accuracy_std': bio_acc_std,
                'roc_auc_mean': bio_roc_mean,
                'roc_auc_std': bio_roc_std,
                'accuracies': bio_accuracies,
                'roc_aucs': bio_roc_aucs,
            }
        }
        
        return self.results
    
    def print_report(self):
        """Print formatted comparison report."""
        print(f"\n{'='*70}")
        print(f"BENCHMARK RESULTS")
        print(f"{'='*70}")
        
        bow = self.results['bag_of_words']
        bio = self.results['biophysical']
        
        print(f"\n📊 ACCURACY:")
        print(f"  BOW Model       : {bow['accuracy_mean']:.3f} ± {bow['accuracy_std']:.3f}")
        print(f"  Biophysical     : {bio['accuracy_mean']:.3f} ± {bio['accuracy_std']:.3f}")
        print(f"  ✓ Improvement   : +{(bio['accuracy_mean'] - bow['accuracy_mean'])*100:.1f}%")
        
        print(f"\n📈 ROC-AUC:")
        print(f"  BOW Model       : {bow['roc_auc_mean']:.3f} ± {bow['roc_auc_std']:.3f}")
        print(f"  Biophysical     : {bio['roc_auc_mean']:.3f} ± {bio['roc_auc_std']:.3f}")
        print(f"  ✓ Improvement   : +{(bio['roc_auc_mean'] - bow['roc_auc_mean'])*100:.1f}%")
        
        print(f"\n💡 HYPOTHESIS VALIDATION:")
        if bow['accuracy_mean'] < 0.55:
            print(f"  ✓ BOW model achieves ~random performance on novel mutations (as expected)")
        if bio['accuracy_mean'] > 0.65:
            print(f"  ✓ Biophysical model achieves strong performance on novel mutations")
        if bio['accuracy_mean'] > bow['accuracy_mean']:
            print(f"  ✓ Biophysical model OUTPERFORMS baseline on novel mutation generalization")
        
        print(f"\n{'='*70}")
    
    def save_results(self) -> Path:
        """Save benchmark results to CSV."""
        results_list = [
            {
                'model': 'Bag-of-Words (Baseline)',
                'accuracy_mean': self.results['bag_of_words']['accuracy_mean'],
                'accuracy_std': self.results['bag_of_words']['accuracy_std'],
                'roc_auc_mean': self.results['bag_of_words']['roc_auc_mean'],
                'roc_auc_std': self.results['bag_of_words']['roc_auc_std'],
            },
            {
                'model': 'Biophysical (Proposed)',
                'accuracy_mean': self.results['biophysical']['accuracy_mean'],
                'accuracy_std': self.results['biophysical']['accuracy_std'],
                'roc_auc_mean': self.results['biophysical']['roc_auc_mean'],
                'roc_auc_std': self.results['biophysical']['roc_auc_std'],
            }
        ]
        
        results_df = pd.DataFrame(results_list)
        results_path = self.output_dir / "benchmark_results.csv"
        results_df.to_csv(results_path, index=False)
        
        print(f"\n✅ Saved results to {results_path}")
        
        return results_path
    
    def plot_comparison(self) -> Path:
        """Generate comparison visualization."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle('Benchmark: Bag-of-Words vs Biophysical Model\n(Novel Mutation Generalization)', 
                     fontsize=14, fontweight='bold')
        
        models = ['Bag-of-Words', 'Biophysical']
        bow = self.results['bag_of_words']
        bio = self.results['biophysical']
        
        # Accuracy comparison
        acc_means = [bow['accuracy_mean'], bio['accuracy_mean']]
        acc_stds = [bow['accuracy_std'], bio['accuracy_std']]
        colors = ['#FF6B6B', '#4ECDC4']
        
        ax1 = axes[0]
        bars1 = ax1.bar(models, acc_means, yerr=acc_stds, capsize=10, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
        ax1.set_ylabel('Accuracy', fontsize=12, fontweight='bold')
        ax1.set_ylim([0, 1])
        ax1.set_title('Accuracy on Novel Mutations', fontsize=12, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, (bar, val) in enumerate(zip(bars1, acc_means)):
            ax1.text(bar.get_x() + bar.get_width()/2, val + 0.05, f'{val:.3f}', 
                    ha='center', fontsize=11, fontweight='bold')
        
        # ROC-AUC comparison
        roc_means = [bow['roc_auc_mean'], bio['roc_auc_mean']]
        roc_stds = [bow['roc_auc_std'], bio['roc_auc_std']]
        
        ax2 = axes[1]
        bars2 = ax2.bar(models, roc_means, yerr=roc_stds, capsize=10, color=colors, alpha=0.7, edgecolor='black', linewidth=2)
        ax2.set_ylabel('ROC-AUC', fontsize=12, fontweight='bold')
        ax2.set_ylim([0, 1])
        ax2.set_title('ROC-AUC on Novel Mutations', fontsize=12, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, (bar, val) in enumerate(zip(bars2, roc_means)):
            ax2.text(bar.get_x() + bar.get_width()/2, val + 0.05, f'{val:.3f}', 
                    ha='center', fontsize=11, fontweight='bold')
        
        # Add improvement annotation
        acc_improvement = (bio['accuracy_mean'] - bow['accuracy_mean']) * 100
        roc_improvement = (bio['roc_auc_mean'] - bow['roc_auc_mean']) * 100
        
        fig.text(0.5, 0.02, 
                f'Biophysical Model Improvement: Accuracy +{acc_improvement:.1f}% | ROC-AUC +{roc_improvement:.1f}%',
                ha='center', fontsize=11, fontweight='bold', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout(rect=[0, 0.05, 1, 1])
        
        plot_path = self.output_dir / "benchmark_comparison.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"✅ Saved plot to {plot_path}")
        
        plt.close()
        return plot_path


def main():
    """Main execution."""
    benchmark = ModelBenchmark(
        data_path="data/processed_features.csv",
        output_dir="models"
    )
    
    # Run benchmark
    results = benchmark.benchmark()
    
    # Print and save results
    benchmark.print_report()
    benchmark.save_results()
    benchmark.plot_comparison()
    
    print(f"\n🏆 BENCHMARK COMPLETE")


if __name__ == "__main__":
    main()
