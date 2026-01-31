"""
Data Processing Pipeline for AMR Mutation Analysis

Loads raw mutation data, encodes mutations into biophysical vectors,
validates phenotype labels, and saves processed features.
"""

from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import numpy as np

from .features import BiophysicalEncoder


class AMRDataPipeline:
    """
    ETL pipeline for antimicrobial resistance mutation data.
    
    Flow:
    1. Load raw CSV with mutation strings
    2. Parse and encode mutations using BiophysicalEncoder
    3. Validate phenotype labels (convert to binary 0/1)
    4. Save processed features
    """
    
    def __init__(self, input_path: str = "data/raw/raw_amr.csv", 
                 output_path: str = "data/processed_features.csv"):
        """
        Initialize pipeline with input/output paths.
        
        Args:
            input_path: Path to raw data CSV
            output_path: Path to save processed features
        """
        self.input_path = Path(input_path)
        self.output_path = Path(output_path)
        self.encoder = BiophysicalEncoder()
        self.processed_df = None
    
    def load_raw_data(self) -> pd.DataFrame:
        """
        Load raw mutation data from CSV.
        
        Expected columns:
        - mutation: Mutation string (e.g., "S83L", "gyrA_S83L")
        - antibiotic: Antibiotic name (e.g., "Ciprofloxacin")
        - phenotype: Resistance status (e.g., "R", "S", 1, 0)
        
        Returns:
            Loaded dataframe
        
        Raises:
            FileNotFoundError: If input file doesn't exist
            ValueError: If required columns are missing
        """
        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_path}")
        
        df = pd.read_csv(self.input_path)
        
        required_cols = {'mutation', 'antibiotic', 'phenotype'}
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        print(f"✓ Loaded {len(df)} rows from {self.input_path}")
        return df
    
    def _normalize_phenotype(self, phenotype: any) -> Optional[int]:
        """
        Convert phenotype to binary (0/1).
        
        Args:
            phenotype: Phenotype value (e.g., "R", "S", "Resistant", 1)
        
        Returns:
            0 (Susceptible), 1 (Resistant), or None if invalid
        """
        if isinstance(phenotype, (int, float)):
            if phenotype in [0, 1]:
                return int(phenotype)
            return None
        
        phenotype_str = str(phenotype).strip().upper()
        
        if phenotype_str in ['R', 'RESISTANT', '1']:
            return 1
        elif phenotype_str in ['S', 'SUSCEPTIBLE', '0']:
            return 0
        else:
            return None
    
    def encode_mutations(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Encode mutation strings into biophysical features.
        
        Args:
            df: Dataframe with 'mutation' column
        
        Returns:
            Dataframe with added feature columns:
            - delta_hydrophobicity
            - delta_charge
            - delta_mw
            - is_aromatic_change
            - is_proline_change
        """
        print("\n🧬 Encoding mutations into biophysical vectors...")
        
        features_list = []
        valid_rows = []
        dropped = 0
        
        for idx, row in df.iterrows():
            mutation = row['mutation']
            features = self.encoder.get_features(mutation)
            
            if features is None:
                dropped += 1
                continue
            
            features_list.append(features)
            valid_rows.append(idx)
        
        print(f"   ✓ Successfully encoded {len(features_list)} mutations")
        print(f"   ⚠️  Dropped {dropped} invalid mutations (parsing failed)")
        
        # Create dataframe from encoded features
        features_df = pd.DataFrame(features_list)
        
        # Keep only valid rows from original dataframe
        df_valid = df.iloc[valid_rows].reset_index(drop=True)
        df_encoded = pd.concat([df_valid, features_df], axis=1)
        
        return df_encoded
    
    def validate_phenotypes(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize and validate phenotype column.
        
        Args:
            df: Dataframe with 'phenotype' column
        
        Returns:
            Dataframe with cleaned phenotype column (0/1 binary)
        """
        print("\n🏷️  Validating phenotype labels...")
        
        df = df.copy()
        phenotypes = []
        dropped_indices = []
        
        for idx, phenotype in enumerate(df['phenotype']):
            normalized = self._normalize_phenotype(phenotype)
            if normalized is None:
                dropped_indices.append(idx)
            else:
                phenotypes.append(normalized)
        
        # Remove rows with invalid phenotypes
        df = df.drop(dropped_indices).reset_index(drop=True)
        df['phenotype'] = phenotypes
        
        print(f"   ✓ Valid phenotypes: {len(phenotypes)}")
        print(f"   ⚠️  Dropped {len(dropped_indices)} rows with invalid phenotypes")
        print(f"   Distribution: {np.sum(phenotypes)} Resistant, {len(phenotypes) - np.sum(phenotypes)} Susceptible")
        
        return df
    
    def process(self) -> Tuple[pd.DataFrame, int, int]:
        """
        Execute full ETL pipeline.
        
        Returns:
            Tuple of (processed_dataframe, rows_loaded, rows_processed)
        """
        print("="*70)
        print("AMR MUTATION DATA PIPELINE")
        print("="*70)
        
        # 1. Load
        df = self.load_raw_data()
        initial_count = len(df)
        
        # 2. Encode mutations
        df = self.encode_mutations(df)
        
        # 3. Validate phenotypes
        df = self.validate_phenotypes(df)
        
        # 4. Select output columns
        feature_cols = [
            'antibiotic', 'phenotype',
            'delta_hydrophobicity', 'delta_charge', 'delta_mw',
            'is_aromatic_change', 'is_proline_change',
            'wt', 'position', 'mutant'  # Include for traceability
        ]
        
        df_output = df[[col for col in feature_cols if col in df.columns]]
        
        self.processed_df = df_output
        
        print(f"\n{'='*70}")
        print(f"PIPELINE COMPLETE")
        print(f"{'='*70}")
        print(f"Initial rows: {initial_count}")
        print(f"Final rows: {len(df_output)}")
        print(f"Rows dropped: {initial_count - len(df_output)} ({(1 - len(df_output)/initial_count)*100:.1f}%)")
        
        return df_output, initial_count, len(df_output)
    
    def save(self) -> Path:
        """
        Save processed dataframe to CSV.
        
        Returns:
            Path to saved file
        """
        if self.processed_df is None:
            raise ValueError("No processed data. Run process() first.")
        
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.processed_df.to_csv(self.output_path, index=False)
        
        print(f"\n✅ Saved processed features to {self.output_path}")
        print(f"   Columns: {list(self.processed_df.columns)}")
        print(f"   Shape: {self.processed_df.shape}")
        
        return self.output_path


def main():
    """Main execution."""
    pipeline = AMRDataPipeline(
        input_path="data/raw/raw_amr.csv",
        output_path="data/processed_features.csv"
    )
    
    df_processed, n_initial, n_final = pipeline.process()
    pipeline.save()
    
    # Display sample
    print(f"\nFirst 5 rows of processed data:")
    print(df_processed.head())


if __name__ == "__main__":
    main()
