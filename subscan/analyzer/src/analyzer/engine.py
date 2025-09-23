"""
Core mutation analysis engine for SubScan Analyzer

This module implements sophisticated bioinformatics algorithms for:
1. Reference-driven mutation detection in aligned protein sequences
2. Position-specific amino acid change analysis
3. Statistical co-occurrence pattern analysis across genomes
4. Comprehensive mutation cataloging for antibiotic resistance research
"""

import os
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Set, Optional, Any
from collections import defaultdict, Counter
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MutationAnalyzer:
    """
    High-performance mutation detection and co-occurrence analysis engine.
    
    This class provides core scientific functionality for the SubScan pipeline,
    analyzing aligned protein sequences to identify amino acid mutations and
    their co-occurrence patterns across genome collections.
    
    Key Features:
    - Reference-driven mutation detection
    - Position-specific amino acid change tracking
    - Multi-genome co-occurrence analysis
    - Statistical significance testing
    - Comprehensive mutation cataloging
    """
    
    def __init__(self, reference_pattern: str = "reference", min_coverage: float = 0.8):
        """
        Initialize the MutationAnalyzer.
        
        Args:
            reference_pattern: Pattern to identify reference sequences in alignments
            min_coverage: Minimum sequence coverage required for analysis (0.0-1.0)
        """
        self.reference_pattern = reference_pattern.lower()
        self.min_coverage = min_coverage
        self.mutation_catalog = []
        self.cooccurrence_data = defaultdict(list)
        
        logger.info(f"MutationAnalyzer initialized with reference pattern: '{reference_pattern}'")
        logger.info(f"Minimum coverage threshold: {min_coverage:.1%}")
    
    def identify_mutations(self, alignment_file_path: str) -> List[Dict[str, Any]]:
        """
        Identify amino acid mutations in an aligned FASTA file.
        
        This method performs sophisticated mutation detection by:
        1. Identifying reference sequence(s) in the alignment
        2. Comparing each genome sequence position-by-position
        3. Cataloging substitutions, insertions, and deletions
        4. Recording exact positions and amino acid changes
        
        Args:
            alignment_file_path: Path to aligned FASTA file from WildTypeAligner
            
        Returns:
            List of mutation dictionaries with detailed information
            
        Each mutation dictionary contains:
        - genome_id: Identifier of the genome containing the mutation
        - gene_family: Gene family name (extracted from filename)
        - position: Amino acid position in the alignment (1-based)
        - reference_aa: Reference amino acid at this position
        - mutated_aa: Mutated amino acid at this position
        - mutation_type: Type of mutation (substitution, insertion, deletion)
        - mutation_code: Standard notation (e.g., "A123T" for Ala→Thr at position 123)
        """
        logger.info(f"Analyzing mutations in: {alignment_file_path}")
        
        if not os.path.exists(alignment_file_path):
            raise FileNotFoundError(f"Alignment file not found: {alignment_file_path}")
        
        # Extract gene family from filename
        gene_family = Path(alignment_file_path).stem.replace('_aligned', '')
        logger.info(f"Gene family identified: {gene_family}")
        
        mutations = []
        
        try:
            # Parse aligned sequences
            sequences = list(SeqIO.parse(alignment_file_path, "fasta"))
            if not sequences:
                logger.warning(f"No sequences found in {alignment_file_path}")
                return mutations
            
            logger.info(f"Loaded {len(sequences)} sequences from alignment")
            
            # Identify reference sequence(s)
            reference_seqs = self._identify_reference_sequences(sequences)
            if not reference_seqs:
                logger.warning(f"No reference sequences found in {alignment_file_path}")
                return mutations
            
            logger.info(f"Found {len(reference_seqs)} reference sequence(s)")
            
            # Use first reference sequence as primary reference
            reference_seq = reference_seqs[0]
            reference_str = str(reference_seq.seq).upper()
            
            # Analyze each non-reference sequence
            for seq_record in sequences:
                if seq_record.id in [ref.id for ref in reference_seqs]:
                    continue  # Skip reference sequences
                
                genome_id = seq_record.id
                sequence_str = str(seq_record.seq).upper()
                
                # Check sequence coverage
                coverage = self._calculate_coverage(sequence_str)
                if coverage < self.min_coverage:
                    logger.warning(f"Skipping {genome_id}: coverage {coverage:.1%} below threshold")
                    continue
                
                # Identify mutations in this sequence
                seq_mutations = self._compare_sequences(
                    reference_str, sequence_str, reference_seq.id, genome_id, gene_family
                )
                mutations.extend(seq_mutations)
                
                logger.debug(f"Found {len(seq_mutations)} mutations in {genome_id}")
            
            logger.info(f"Total mutations identified: {len(mutations)}")
            
        except Exception as e:
            logger.error(f"Error analyzing mutations in {alignment_file_path}: {str(e)}")
            raise
        
        return mutations
    
    def analyze_cooccurrence(self, all_mutations: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Analyze co-occurrence patterns of mutations across genomes.
        
        This method performs sophisticated statistical analysis to identify:
        1. Mutations that frequently appear together in the same genome
        2. Gene-level correlation patterns
        3. Potential compensatory mutation relationships
        4. Cross-genome mutation patterns
        
        Args:
            all_mutations: List of all mutations from multiple gene families
            
        Returns:
            Dictionary containing comprehensive co-occurrence analysis:
            - mutation_pairs: Pairwise co-occurrence frequencies
            - gene_correlations: Gene-level correlation matrix
            - genome_patterns: Per-genome mutation patterns
            - statistical_summary: Overall analysis statistics
        """
        logger.info(f"Analyzing co-occurrence patterns for {len(all_mutations)} mutations")
        
        if not all_mutations:
            logger.warning("No mutations provided for co-occurrence analysis")
            return {
                "mutation_pairs": {},
                "gene_correlations": {},
                "genome_patterns": {},
                "statistical_summary": {"total_mutations": 0, "total_genomes": 0}
            }
        
        # Organize mutations by genome
        genome_mutations = defaultdict(list)
        gene_families = set()
        
        for mutation in all_mutations:
            genome_id = mutation["genome_id"]
            gene_family = mutation["gene_family"]
            genome_mutations[genome_id].append(mutation)
            gene_families.add(gene_family)
        
        logger.info(f"Analyzing {len(genome_mutations)} genomes across {len(gene_families)} gene families")
        
        # Calculate pairwise mutation co-occurrence
        mutation_pairs = self._calculate_mutation_cooccurrence(genome_mutations)
        
        # Calculate gene-level correlations
        gene_correlations = self._calculate_gene_correlations(genome_mutations, gene_families)
        
        # Generate genome-specific patterns
        genome_patterns = self._generate_genome_patterns(genome_mutations)
        
        # Calculate statistical summary
        statistical_summary = {
            "total_mutations": len(all_mutations),
            "total_genomes": len(genome_mutations),
            "gene_families": len(gene_families),
            "avg_mutations_per_genome": len(all_mutations) / len(genome_mutations) if genome_mutations else 0,
            "mutation_pair_count": len(mutation_pairs),
            "significant_correlations": sum(1 for corr in gene_correlations.values() if abs(corr.get("correlation", 0)) > 0.5)
        }
        
        logger.info(f"Co-occurrence analysis complete: {statistical_summary['mutation_pair_count']} pairs analyzed")
        
        return {
            "mutation_pairs": mutation_pairs,
            "gene_correlations": gene_correlations,
            "genome_patterns": genome_patterns,
            "statistical_summary": statistical_summary
        }
    
    def _identify_reference_sequences(self, sequences: List[SeqRecord]) -> List[SeqRecord]:
        """Identify reference sequences based on naming pattern."""
        references = []
        for seq in sequences:
            if self.reference_pattern in seq.id.lower() or self.reference_pattern in seq.description.lower():
                references.append(seq)
        return references
    
    def _calculate_coverage(self, sequence: str) -> float:
        """Calculate sequence coverage (non-gap characters)."""
        if not sequence:
            return 0.0
        non_gap_chars = sum(1 for char in sequence if char not in ['-', 'N', 'X'])
        return non_gap_chars / len(sequence)
    
    def _compare_sequences(self, reference: str, sequence: str, ref_id: str, 
                          genome_id: str, gene_family: str) -> List[Dict[str, Any]]:
        """Compare two sequences position-by-position to identify mutations."""
        mutations = []
        
        # Ensure sequences are same length (aligned)
        max_length = max(len(reference), len(sequence))
        reference = reference.ljust(max_length, '-')
        sequence = sequence.ljust(max_length, '-')
        
        for pos in range(max_length):
            ref_aa = reference[pos] if pos < len(reference) else '-'
            seq_aa = sequence[pos] if pos < len(sequence) else '-'
            
            # Skip if both are gaps or identical
            if ref_aa == seq_aa or (ref_aa == '-' and seq_aa == '-'):
                continue
            
            # Skip positions with ambiguous characters
            if ref_aa in ['N', 'X'] or seq_aa in ['N', 'X']:
                continue
            
            # Determine mutation type and create mutation record
            mutation_type = self._classify_mutation(ref_aa, seq_aa)
            mutation_code = self._generate_mutation_code(ref_aa, seq_aa, pos + 1)
            
            mutation = {
                "genome_id": genome_id,
                "gene_family": gene_family,
                "position": pos + 1,  # 1-based position
                "reference_aa": ref_aa,
                "mutated_aa": seq_aa,
                "mutation_type": mutation_type,
                "mutation_code": mutation_code,
                "reference_id": ref_id
            }
            
            mutations.append(mutation)
        
        return mutations
    
    def _classify_mutation(self, ref_aa: str, seq_aa: str) -> str:
        """Classify the type of mutation."""
        if ref_aa == '-' and seq_aa != '-':
            return "insertion"
        elif ref_aa != '-' and seq_aa == '-':
            return "deletion"
        else:
            return "substitution"
    
    def _generate_mutation_code(self, ref_aa: str, seq_aa: str, position: int) -> str:
        """Generate standard mutation notation (e.g., A123T)."""
        # Handle gaps in mutation notation
        ref_display = ref_aa if ref_aa != '-' else 'del'
        seq_display = seq_aa if seq_aa != '-' else 'ins'
        
        return f"{ref_display}{position}{seq_display}"
    
    def _calculate_mutation_cooccurrence(self, genome_mutations: Dict[str, List[Dict]]) -> Dict[str, Dict]:
        """Calculate pairwise mutation co-occurrence frequencies."""
        cooccurrence_counts = defaultdict(int)
        mutation_counts = defaultdict(int)
        
        # Count mutations and co-occurrences
        for genome_id, mutations in genome_mutations.items():
            mutation_codes = [mut["mutation_code"] for mut in mutations]
            
            # Count individual mutations
            for code in mutation_codes:
                mutation_counts[code] += 1
            
            # Count pairwise co-occurrences
            for i, code1 in enumerate(mutation_codes):
                for code2 in mutation_codes[i+1:]:
                    pair = tuple(sorted([code1, code2]))
                    cooccurrence_counts[pair] += 1
        
        # Calculate co-occurrence statistics
        cooccurrence_stats = {}
        total_genomes = len(genome_mutations)
        
        for (mut1, mut2), count in cooccurrence_counts.items():
            if count >= 2:  # Only include pairs that occur together at least twice
                freq1 = mutation_counts[mut1]
                freq2 = mutation_counts[mut2]
                
                # Calculate expected vs observed co-occurrence
                expected = (freq1 * freq2) / total_genomes
                observed = count
                
                cooccurrence_stats[f"{mut1}+{mut2}"] = {
                    "mutation_1": mut1,
                    "mutation_2": mut2,
                    "cooccurrence_count": count,
                    "mutation_1_frequency": freq1,
                    "mutation_2_frequency": freq2,
                    "cooccurrence_frequency": count / total_genomes,
                    "expected_cooccurrence": expected,
                    "enrichment_ratio": observed / expected if expected > 0 else float('inf')
                }
        
        return cooccurrence_stats
    
    def _calculate_gene_correlations(self, genome_mutations: Dict[str, List[Dict]], 
                                   gene_families: Set[str]) -> Dict[str, Dict]:
        """Calculate gene-level mutation correlations."""
        gene_correlations = {}
        
        # Create binary mutation matrix for genes
        gene_mutation_matrix = {}
        for genome_id, mutations in genome_mutations.items():
            genome_genes = set(mut["gene_family"] for mut in mutations)
            gene_mutation_matrix[genome_id] = {gene: 1 if gene in genome_genes else 0 
                                             for gene in gene_families}
        
        # Calculate pairwise gene correlations
        gene_list = list(gene_families)
        for i, gene1 in enumerate(gene_list):
            for gene2 in gene_list[i+1:]:
                # Get binary vectors for both genes
                gene1_vector = [gene_mutation_matrix[gid][gene1] for gid in genome_mutations.keys()]
                gene2_vector = [gene_mutation_matrix[gid][gene2] for gid in genome_mutations.keys()]
                
                # Calculate simple correlation
                correlation = self._calculate_correlation(gene1_vector, gene2_vector)
                
                if abs(correlation) > 0.1:  # Only store meaningful correlations
                    gene_correlations[f"{gene1}+{gene2}"] = {
                        "gene_1": gene1,
                        "gene_2": gene2,
                        "correlation": correlation,
                        "gene_1_mutation_count": sum(gene1_vector),
                        "gene_2_mutation_count": sum(gene2_vector),
                        "total_genomes": len(gene1_vector)
                    }
        
        return gene_correlations
    
    def _calculate_correlation(self, vector1: List[int], vector2: List[int]) -> float:
        """Calculate Pearson correlation coefficient."""
        if len(vector1) != len(vector2) or len(vector1) == 0:
            return 0.0
        
        n = len(vector1)
        sum1 = sum(vector1)
        sum2 = sum(vector2)
        sum1_sq = sum(x * x for x in vector1)
        sum2_sq = sum(x * x for x in vector2)
        sum_products = sum(vector1[i] * vector2[i] for i in range(n))
        
        numerator = n * sum_products - sum1 * sum2
        denominator = ((n * sum1_sq - sum1 * sum1) * (n * sum2_sq - sum2 * sum2)) ** 0.5
        
        if denominator == 0:
            return 0.0
        
        return numerator / denominator
    
    def _generate_genome_patterns(self, genome_mutations: Dict[str, List[Dict]]) -> Dict[str, Dict]:
        """Generate mutation patterns for each genome."""
        genome_patterns = {}
        
        for genome_id, mutations in genome_mutations.items():
            gene_families = set(mut["gene_family"] for mut in mutations)
            mutation_types = Counter(mut["mutation_type"] for mut in mutations)
            
            genome_patterns[genome_id] = {
                "total_mutations": len(mutations),
                "affected_gene_families": list(gene_families),
                "gene_family_count": len(gene_families),
                "mutation_type_distribution": dict(mutation_types),
                "mutation_codes": [mut["mutation_code"] for mut in mutations]
            }
        
        return genome_patterns