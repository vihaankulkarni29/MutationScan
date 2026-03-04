"""
MutationScorer Module: In-Silico Biochemist for Amino Acid Substitution Severity Assessment

This module evaluates the structural and evolutionary severity of amino acid mutations
using BLOSUM62 (evolutionary distance) and Grantham Distance (chemical properties).

The scoring engine produces:
- Individual mutation severity scores (0-100 scale)
- Categories: BENIGN, UNCERTAIN, PATHOGENIC
- Network-level epistatic severity assessments
"""

import math
import logging
import re
from typing import List, Tuple, Dict, Any

logger = logging.getLogger(__name__)


class MutationScorer:
    """
    In-Silico Biochemist: Scores amino acid substitutions based on evolutionary and chemical properties.
    
    Uses BLOSUM62 matrix for evolutionary distance and Grantham Distance for chemical properties.
    Combines both metrics into a comprehensive severity score (0-100 scale).
    
    Attributes:
        BLOSUM62 (dict): Matrix of evolutionary substitution scores
        AA_PROPERTIES (dict): Amino acid composition, polarity, and volume
    """
    
    def __init__(self):
        """Initialize MutationScorer with BLOSUM62 and amino acid property tables."""
        
        # BLOSUM62 Substitution Matrix (standard)
        self.BLOSUM62 = {
            'A': {'A':  4, 'R': -1, 'N': -2, 'D': -2, 'C':  0, 'Q': -1, 'E': -1, 'G':  0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S':  1, 'T':  0, 'W': -3, 'Y': -2, 'V':  0},
            'R': {'A': -1, 'R':  5, 'N':  0, 'D': -2, 'C': -3, 'Q':  1, 'E':  0, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
            'N': {'A': -2, 'R':  0, 'N':  6, 'D':  1, 'C': -3, 'Q':  0, 'E':  0, 'G':  0, 'H':  1, 'I': -3, 'L': -3, 'K':  0, 'M': -2, 'F': -3, 'P': -2, 'S':  1, 'T':  0, 'W': -4, 'Y': -2, 'V': -3},
            'D': {'A': -2, 'R': -2, 'N':  1, 'D':  6, 'C': -3, 'Q':  0, 'E':  2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
            'C': {'A':  0, 'R': -3, 'N': -3, 'D': -3, 'C':  9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
            'Q': {'A': -1, 'R':  1, 'N':  0, 'D':  0, 'C': -3, 'Q':  5, 'E':  2, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  1, 'M':  0, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
            'E': {'A': -1, 'R':  0, 'N':  0, 'D':  2, 'C': -4, 'Q':  2, 'E':  5, 'G': -2, 'H':  0, 'I': -3, 'L': -3, 'K':  1, 'M': -2, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
            'G': {'A':  0, 'R': -2, 'N':  0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G':  6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S':  0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
            'H': {'A': -2, 'R':  0, 'N':  1, 'D': -1, 'C': -3, 'Q':  0, 'E':  0, 'G': -2, 'H':  8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y':  2, 'V': -3},
            'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I':  4, 'L':  2, 'K': -3, 'M':  1, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V':  3},
            'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I':  2, 'L':  4, 'K': -2, 'M':  2, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V':  1},
            'K': {'A': -1, 'R':  2, 'N':  0, 'D': -1, 'C': -3, 'Q':  1, 'E':  1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K':  5, 'M': -1, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
            'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q':  0, 'E': -2, 'G': -3, 'H': -2, 'I':  1, 'L':  2, 'K': -1, 'M':  5, 'F':  0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V':  1},
            'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I':  0, 'L':  0, 'K': -3, 'M':  0, 'F':  6, 'P': -4, 'S': -3, 'T': -3, 'W':  1, 'Y':  3, 'V': -1},
            'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P':  7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
            'S': {'A':  1, 'R': -1, 'N':  1, 'D':  0, 'C': -1, 'Q':  0, 'E':  0, 'G':  0, 'H': -1, 'I': -2, 'L': -2, 'K':  0, 'M': -1, 'F': -3, 'P': -1, 'S':  4, 'T':  1, 'W': -3, 'Y': -2, 'V': -2},
            'T': {'A':  0, 'R': -1, 'N':  0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -3, 'P': -1, 'S':  1, 'T':  5, 'W': -2, 'Y': -2, 'V':  0},
            'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F':  1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y':  2, 'V': -3},
            'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H':  2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F':  3, 'P': -3, 'S': -2, 'T': -2, 'W':  2, 'Y':  7, 'V': -1},
            'V': {'A':  0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I':  3, 'L':  1, 'K': -2, 'M':  1, 'F': -1, 'P': -2, 'S': -2, 'T':  0, 'W': -3, 'Y': -1, 'V':  4},
        }
        
        # Valid amino acids (including stop codon)
        self.VALID_AA = set("ARNDCQEGHILKMFPSTWYV*")
        
        # Grantham Distance Properties: Composition (c), Polarity (p), Volume (v)
        self.AA_PROPERTIES = {
            'A': {'c': 0.0, 'p': 8.1, 'v': 31.0},
            'R': {'c': 65.6, 'p': 10.5, 'v': 124.0},
            'N': {'c': 133.0, 'p': 11.6, 'v': 56.0},
            'D': {'c': 27.7, 'p': 13.0, 'v': 54.0},
            'C': {'c': 59.9, 'p': 5.5, 'v': 55.0},
            'Q': {'c': 79.7, 'p': 10.5, 'v': 85.0},
            'E': {'c': 43.1, 'p': 12.3, 'v': 83.0},
            'G': {'c': 0.0, 'p': 9.0, 'v': 3.0},
            'H': {'c': 74.0, 'p': 10.4, 'v': 96.0},
            'I': {'c': 0.0, 'p': 5.2, 'v': 111.0},
            'L': {'c': 0.0, 'p': 4.9, 'v': 111.0},
            'K': {'c': 58.0, 'p': 11.3, 'v': 119.0},
            'M': {'c': 0.0, 'p': 5.7, 'v': 105.0},
            'F': {'c': 0.0, 'p': 5.2, 'v': 132.0},
            'P': {'c': 39.8, 'p': 8.0, 'v': 32.5},
            'S': {'c': 142.0, 'p': 9.2, 'v': 32.0},
            'T': {'c': 71.0, 'p': 8.6, 'v': 61.0},
            'W': {'c': 130.0, 'p': 5.4, 'v': 170.0},
            'Y': {'c': 71.0, 'p': 6.2, 'v': 136.0},
            'V': {'c': 0.0, 'p': 5.9, 'v': 84.0},
        }
        
        logger.info("MutationScorer initialized with BLOSUM62 and Grantham properties")
    
    def _get_blosum(self, aa1: str, aa2: str) -> float:
        """
        Retrieve BLOSUM62 substitution score for amino acid pair.
        
        Args:
            aa1 (str): First amino acid (uppercase)
            aa2 (str): Second amino acid (uppercase)
        
        Returns:
            float: BLOSUM62 score, or -4 if not found
        """
        aa1 = aa1.upper()
        aa2 = aa2.upper()
        
        if aa1 == aa2:
            return 0 if aa1 not in self.BLOSUM62 else self.BLOSUM62[aa1][aa1]
        
        try:
            return self.BLOSUM62.get(aa1, {}).get(aa2, -4)
        except (KeyError, TypeError):
            return -4
    
    def _calc_grantham(self, aa1: str, aa2: str) -> float:
        """
        Calculate Grantham Distance between two amino acids.
        
        Formula: distance = sqrt(1.833 * (c1-c2)^2 + 0.1018 * (p1-p2)^2 + 0.000399 * (v1-v2)^2)
        
        Args:
            aa1 (str): First amino acid (uppercase)
            aa2 (str): Second amino acid (uppercase)
        
        Returns:
            float: Grantham distance (0.0 if aa1 == aa2)
        """
        aa1 = aa1.upper()
        aa2 = aa2.upper()
        
        if aa1 == aa2:
            return 0.0
        
        if aa1 not in self.AA_PROPERTIES or aa2 not in self.AA_PROPERTIES:
            return 0.0
        
        props1 = self.AA_PROPERTIES[aa1]
        props2 = self.AA_PROPERTIES[aa2]
        
        c_diff = props1['c'] - props2['c']
        p_diff = props1['p'] - props2['p']
        v_diff = props1['v'] - props2['v']
        
        distance = math.sqrt(
            1.833 * (c_diff ** 2) +
            0.1018 * (p_diff ** 2) +
            0.000399 * (v_diff ** 2)
        )
        
        return distance
    
    def parse_mutation(self, mut_string: str) -> Tuple[str, str, str]:
        """
        Strictly extracts and validates Wild-Type, Position, and Mutant amino acids.
        
        Uses regex to enforce strict format: One letter, one or more digits, one letter/asterisk.
        Validates amino acids against VALID_AA set.
        
        Args:
            mut_string (str): Mutation string like "I174V" or "Q45*"
        
        Returns:
            tuple: (wild_type_aa, position, mutant_aa)
            
        Raises:
            TypeError: If input is not a string
            ValueError: If format or amino acids are invalid
        """
        if not isinstance(mut_string, str):
            raise TypeError(f"Mutation must be a string, got {type(mut_string).__name__}")
        
        mut_string = mut_string.strip().upper()
        
        # Strict regex: One letter, one or more digits, one letter or asterisk (*)
        match = re.match(r'^([A-Z])(\d+)([A-Z\*])$', mut_string)
        if not match:
            raise ValueError(f"Invalid mutation format: '{mut_string}'. Expected format: 'I174V'")
        
        wt, pos, mut = match.groups()
        
        # Validate amino acids
        if wt not in self.AA_PROPERTIES and wt != '*':
            raise ValueError(f"Unknown Wild-Type amino acid: {wt}")
        if mut not in self.AA_PROPERTIES and mut != '*':
            raise ValueError(f"Unknown Mutant amino acid: {mut}")
        
        return wt, pos, mut
    
    def score_single(self, mut_string: str) -> Dict[str, Any]:
        """
        Mathematically bounded and type-safe scoring of a single mutation.
        
        Handles edge cases (synonymous mutations, stop codons) and gracefully
        catches all input validation errors without raising exceptions.
        
        Args:
            mut_string: Mutation string like "I174V", "I174I" (synonymous), or "Q45*" (stop)
        
        Returns:
            dict: Contains Mutation, BLOSUM, Grantham, Severity, Category (and Error if caught)
        """
        # Attempt parsing; return error dict if validation fails
        try:
            wt, pos, mut = self.parse_mutation(mut_string)
        except (ValueError, TypeError) as e:
            return {
                'Mutation': str(mut_string),
                'Error': str(e),
                'Severity': None,
                'Category': 'ERROR'
            }
        
        # Edge Case 1: Synonymous Mutation (Wild-type == Mutant, no change)
        if wt == mut:
            return {
                'Mutation': mut_string,
                'BLOSUM': self._get_blosum(wt, mut),
                'Grantham': 0.0,
                'Severity': 0.0,
                'Category': 'BENIGN'
            }
        
        # Edge Case 2: Nonsense Mutation (Stop Codon = Maximum Severity)
        if mut == '*':
            return {
                'Mutation': mut_string,
                'BLOSUM': -4,
                'Grantham': 215.0,
                'Severity': 100.0,
                'Category': 'PATHOGENIC'
            }
        
        # Standard Scoring Path
        blosum = self._get_blosum(wt, mut)
        grantham = self._calc_grantham(wt, mut)
        
        # Normalize BLOSUM: [-4, 11] -> [0, 1] with strict bounds
        blosum_norm = max(0.0, min((blosum + 4) / 15.0, 1.0))
        
        # Normalize Grantham: [0, 215] -> [0, 1] with strict bounds
        grantham_norm = min(grantham / 215.0, 1.0)
        
        # Calculate severity: 40% evolutionary + 60% chemical penalty
        severity_raw = (0.40 * (1 - blosum_norm) + 0.60 * grantham_norm)
        severity_final = max(0.0, min(severity_raw * 100, 100.0))  # Strict bounds: [0, 100]
        
        # Assign category based on bounded severity
        if severity_final < 30:
            category = "BENIGN"
        elif severity_final < 60:
            category = "UNCERTAIN"
        else:
            category = "PATHOGENIC"
        
        return {
            'Mutation': mut_string,
            'BLOSUM': blosum,
            'Grantham': grantham,
            'Severity': round(severity_final, 2),
            'Category': category
        }
    
    def score_network(self, network_mutations: List[str]) -> Dict[str, Any]:
        """
        Score a network of co-occurring mutations (epistatic interactions).
        
        Filters out error results before calculating network-level statistics.
        
        Args:
            network_mutations (list): List of mutation strings like ["I174V", "K10R"]
        
        Returns:
            dict: Network name, individual scores, Mean_Severity, and Max_Severity
        """
        if not network_mutations:
            raise ValueError("Network must contain at least one mutation")
        
        individual_scores = [self.score_single(mut) for mut in network_mutations]
        
        # Filter out error results (where Severity is None)
        valid_scores = [score for score in individual_scores if score.get('Severity') is not None]
        
        if not valid_scores:
            raise ValueError("No valid mutations in network (all had errors)")
        
        severities = [score['Severity'] for score in valid_scores]
        
        network_name = "_".join(network_mutations)
        
        return {
            'Network': network_name,
            'Individual_Scores': valid_scores,
            'Mean_Severity': round(sum(severities) / len(severities), 2),
            'Max_Severity': round(max(severities), 2)
        }


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    print("MutationScorer module loaded (use via imports or tests/run_scorer_test.py)")
