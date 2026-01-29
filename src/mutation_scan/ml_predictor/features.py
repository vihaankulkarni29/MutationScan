"""
Biophysical Feature Encoder for Antimicrobial Resistance Mutations

Converts mutation strings (e.g., "S83L") into numerical biophysical vectors
capturing structural and chemical changes relevant to protein function.
"""

from typing import Dict, Optional, Tuple
import numpy as np


class BiophysicalEncoder:
    """
    Encodes amino acid mutations into biophysical feature vectors.
    
    Features capture changes in:
    - Hydrophobicity (membrane/protein core interactions)
    - Charge (electrostatic interactions, pH dependence)
    - Molecular weight (steric effects)
    - Aromaticity (π-π interactions, binding specificity)
    - Proline (backbone rigidity, alpha-helix breaking)
    """
    
    # Kyte-Doolittle Hydrophobicity Scale (standard in protein biochemistry)
    # Positive = hydrophobic (prefer lipid environment)
    # Negative = hydrophilic (prefer aqueous environment)
    HYDROPHOBICITY_SCALE: Dict[str, float] = {
        'A': 1.8,   # Alanine (small, hydrophobic)
        'R': -4.5,  # Arginine (long, charged)
        'N': -3.5,  # Asparagine (polar)
        'D': -3.5,  # Aspartate (charged)
        'C': 2.5,   # Cysteine (forms disulfides)
        'Q': -3.5,  # Glutamine (polar)
        'E': -3.5,  # Glutamate (charged)
        'G': -0.4,  # Glycine (no side chain)
        'H': -3.2,  # Histidine (aromatic, can be charged)
        'I': 4.5,   # Isoleucine (highly hydrophobic)
        'L': 3.8,   # Leucine (highly hydrophobic)
        'K': -3.9,  # Lysine (charged, positive)
        'M': 1.9,   # Methionine (nonpolar)
        'F': 2.8,   # Phenylalanine (aromatic, hydrophobic)
        'P': -1.6,  # Proline (unique cyclic structure)
        'S': -0.8,  # Serine (small, polar)
        'T': -0.7,  # Threonine (small, polar)
        'W': -0.9,  # Tryptophan (aromatic, can be polar)
        'Y': -1.3,  # Tyrosine (aromatic, can be polar)
        'V': 4.2,   # Valine (hydrophobic)
    }
    
    # Charge at pH 7.4 (physiological pH)
    # Reflects electrostatic interactions in vivo
    CHARGE_SCALE: Dict[str, float] = {
        'A': 0.0,    # Neutral
        'R': 1.0,    # Positive (pKa ~12.5)
        'N': 0.0,    # Neutral (amide)
        'D': -1.0,   # Negative (pKa ~3.9)
        'C': -0.1,   # Slightly negative (pKa ~8.3)
        'Q': 0.0,    # Neutral (amide)
        'E': -1.0,   # Negative (pKa ~4.3)
        'G': 0.0,    # Neutral
        'H': 0.1,    # Slightly positive (pKa ~6.0, ~50% protonated at pH 7.4)
        'I': 0.0,    # Neutral
        'L': 0.0,    # Neutral
        'K': 1.0,    # Positive (pKa ~10.5)
        'M': 0.0,    # Neutral
        'F': 0.0,    # Neutral
        'P': 0.0,    # Neutral (imino group)
        'S': 0.0,    # Neutral (hydroxyl)
        'T': 0.0,    # Neutral (hydroxyl)
        'W': 0.0,    # Neutral (indole)
        'Y': 0.0,    # Neutral at pH 7.4 (pKa ~10.1)
        'V': 0.0,    # Neutral
    }
    
    # Molecular weight in Daltons (Da)
    # Affects protein folding, turnover, and cellular transport
    MOLECULAR_WEIGHT: Dict[str, float] = {
        'A': 89.09,
        'R': 174.20,
        'N': 132.12,
        'D': 133.10,
        'C': 121.15,
        'Q': 146.15,
        'E': 147.13,
        'G': 75.07,
        'H': 155.16,
        'I': 131.17,
        'L': 131.17,
        'K': 146.19,
        'M': 149.21,
        'F': 165.19,
        'P': 115.13,
        'S': 105.09,
        'T': 119.12,
        'W': 204.23,
        'Y': 181.19,
        'V': 117.15,
    }
    
    # Aromatic amino acids
    # Important for π-π stacking, hydrophobic core, and binding specificity
    AROMATIC_AAS: set = {'F', 'Y', 'W', 'H'}
    
    def __init__(self):
        """Initialize the encoder with biophysical scales."""
        pass
    
    def _parse_mutation(self, mutation_str: str) -> Optional[Tuple[str, int, str]]:
        """
        Parse mutation string into components.
        
        Args:
            mutation_str: Mutation in format "S83L" (WT-Position-Mutant)
        
        Returns:
            Tuple of (wild_type_aa, position, mutant_aa) or None if parsing fails
            
        Examples:
            "S83L" → ('S', 83, 'L')
            "gyrA_S83L" → ('S', 83, 'L')
            "Promoter" → None (invalid)
            "A1X" → None (X is not standard amino acid)
        """
        # Remove common prefixes (e.g., "gyrA_", "rpoB_")
        if '_' in mutation_str:
            mutation_str = mutation_str.split('_')[-1]
        
        # Remove whitespace
        mutation_str = mutation_str.strip().upper()
        
        # Expected format: [WT][Position][Mutant]
        # WT and Mutant are single letters (amino acids)
        # Position is numeric
        
        if len(mutation_str) < 3:
            return None
        
        wt = mutation_str[0]
        mutant = mutation_str[-1]
        
        # Extract position (middle part)
        try:
            position = int(mutation_str[1:-1])
        except ValueError:
            return None
        
        # Validate amino acids (only standard 20 amino acids)
        valid_aas = set(self.HYDROPHOBICITY_SCALE.keys())
        if wt not in valid_aas or mutant not in valid_aas:
            return None
        
        # Reject synonymous mutations (no actual change)
        if wt == mutant:
            return None
        
        return (wt, position, mutant)
    
    def get_features(self, mutation_str: str) -> Optional[Dict[str, float]]:
        """
        Extract biophysical features for a mutation.
        
        Args:
            mutation_str: Mutation string (e.g., "S83L")
        
        Returns:
            Dictionary of features or None if parsing fails:
            - delta_hydrophobicity: Change in hydrophobicity (mutant - WT)
            - delta_charge: Change in charge (mutant - WT)
            - delta_mw: Change in molecular weight (mutant - WT) in Da
            - is_aromatic_change: 1 if aromatic status changes, else 0
            - is_proline_change: 1 if WT or mutant is proline, else 0
            
        Biophysical Significance:
            - Large positive delta_hydro: Increase in core hydrophobicity
              (may alter protein folding, increase antibiotic binding pocket)
            - Large delta_charge: May affect electrostatic interactions
              (critical for antibiotic binding, protein-protein interactions)
            - is_proline_change: Proline breaks alpha helices (structural constraint)
            - is_aromatic_change: Alters binding specificity via π-π interactions
        """
        parsed = self._parse_mutation(mutation_str)
        if parsed is None:
            return None
        
        wt, position, mutant = parsed
        
        # Calculate biophysical changes
        delta_hydrophobicity = self.HYDROPHOBICITY_SCALE[mutant] - self.HYDROPHOBICITY_SCALE[wt]
        delta_charge = self.CHARGE_SCALE[mutant] - self.CHARGE_SCALE[wt]
        delta_mw = self.MOLECULAR_WEIGHT[mutant] - self.MOLECULAR_WEIGHT[wt]
        
        # Aromatic status change (binary: True/False → 1/0)
        wt_is_aromatic = wt in self.AROMATIC_AAS
        mut_is_aromatic = mutant in self.AROMATIC_AAS
        is_aromatic_change = 1 if wt_is_aromatic != mut_is_aromatic else 0
        
        # Proline change (important for secondary structure)
        wt_is_proline = wt == 'P'
        mut_is_proline = mutant == 'P'
        is_proline_change = 1 if (wt_is_proline or mut_is_proline) else 0
        
        return {
            'wt': wt,
            'position': position,
            'mutant': mutant,
            'delta_hydrophobicity': float(delta_hydrophobicity),
            'delta_charge': float(delta_charge),
            'delta_mw': float(delta_mw),
            'is_aromatic_change': int(is_aromatic_change),
            'is_proline_change': int(is_proline_change),
        }
    
    def get_features_as_array(self, mutation_str: str) -> Optional[np.ndarray]:
        """
        Get features as numpy array (convenient for ML pipelines).
        
        Returns:
            Array of [delta_hydro, delta_charge, delta_mw, is_aromatic, is_proline]
            or None if parsing fails
        """
        features = self.get_features(mutation_str)
        if features is None:
            return None
        
        return np.array([
            features['delta_hydrophobicity'],
            features['delta_charge'],
            features['delta_mw'],
            features['is_aromatic_change'],
            features['is_proline_change'],
        ], dtype=np.float32)
