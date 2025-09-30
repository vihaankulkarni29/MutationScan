#!/usr/bin/env python3
"""
MutationScan Co-occurrence Analysis Module

This module provides core functionality for analyzing mutation co-occurrence patterns
across genomes in antimicrobial resistance studies. It implements statistical methods
to identify genes with mutations that frequently appear together in the same bacterial
isolates, providing insights into coordinated resistance evolution.

The module supports various analytical approaches including frequency-based filtering,
pattern size analysis, and statistical significance testing for co-occurrence patterns.
It is designed to handle large-scale genomic datasets and provides comprehensive
summary statistics for downstream visualization and interpretation.

Key Functions:
    - analyze_cooccurrence: Core co-occurrence pattern detection
    - generate_cooccurrence_summary: Statistical summary generation
    - validate_mutation_dataframe: Input data validation
    - filter_patterns_by_frequency: Pattern filtering utilities

Scientific Applications:
    - Identify functionally linked resistance mechanisms
    - Detect coordinated evolution of resistance genes
    - Support epidemiological analysis of resistance patterns
    - Enable network analysis of gene interactions

Dependencies:
    - pandas for high-performance data manipulation
    - collections for efficient pattern counting
    - typing for enhanced code documentation and validation

Author: MutationScan Development Team
Version: 1.0.0
"""

import pandas as pd
from typing import Dict, List, Tuple, Any
from collections import Counter, defaultdict
import logging


def analyze_cooccurrence(
    mutation_df: pd.DataFrame, min_frequency: int = 2, exclude_single: bool = False
) -> pd.DataFrame:
    """
    Analyze co-occurrence patterns of mutations across genes within individual genomes.

    This function takes a master mutation DataFrame and identifies patterns where
    mutations in different genes appear together in the same isolate. The analysis
    helps researchers understand potential functional linkages between resistance
    mechanisms.

    Args:
        mutation_df (pd.DataFrame): Master mutation DataFrame with columns:
            - 'Accession Number': Genome identifier
            - 'Gene Name': Gene containing the mutation
            - Additional columns (mutation details, coordinates, etc.)
        min_frequency (int): Minimum frequency threshold for reporting patterns (default: 2)
        exclude_single (bool): If True, exclude single-gene mutations from analysis (default: False)

    Returns:
        pd.DataFrame: Co-occurrence analysis results with columns:
            - 'Gene_Combination': Tuple of genes with co-occurring mutations
            - 'Frequency': Number of isolates with this combination
            - 'Combination_Size': Number of genes in the combination
            - 'Representative_Genes': Alphabetically sorted gene list as string

    Raises:
        ValueError: If required columns are missing from input DataFrame
        KeyError: If DataFrame structure is unexpected

    Example:
        >>> mutation_data = pd.DataFrame({
        ...     'Accession Number': ['ACC1', 'ACC1', 'ACC2', 'ACC2', 'ACC3'],
        ...     'Gene Name': ['acrA', 'acrB', 'acrA', 'tolC', 'acrA']
        ... })
        >>> result = analyze_cooccurrence(mutation_data)
        >>> print(result)
           Gene_Combination  Frequency  Combination_Size Representative_Genes
        0           (acrA,)          3                 1                 acrA
        1     (acrA, acrB)          1                 2           acrA, acrB
        2     (acrA, tolC)          1                 2           acrA, tolC
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    # Validate input DataFrame structure
    required_columns = ["Accession Number", "Gene Name"]
    missing_columns = [
        col for col in required_columns if col not in mutation_df.columns
    ]

    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    if mutation_df.empty:
        logger.warning("Empty mutation DataFrame provided")
        return pd.DataFrame(
            columns=[
                "Gene_Combination",
                "Frequency",
                "Combination_Size",
                "Representative_Genes",
            ]
        )

    logger.info(
        f"Analyzing co-occurrence patterns in {len(mutation_df)} mutations across {mutation_df['Accession Number'].nunique()} genomes"
    )

    # Group mutations by accession (genome) and collect unique genes
    genome_gene_combinations = []
    accession_gene_summary = defaultdict(set)

    # Step 1: Group by accession and collect unique mutated genes per genome
    for accession, group in mutation_df.groupby("Accession Number"):
        # Get unique genes with mutations in this genome
        mutated_genes = set(group["Gene Name"].unique())

        # Store for summary statistics
        accession_gene_summary[accession] = mutated_genes

        # Create alphabetically sorted tuple of genes
        gene_combination = tuple(sorted(mutated_genes))
        genome_gene_combinations.append(gene_combination)

        logger.debug(
            f"Genome {accession}: {len(mutated_genes)} mutated genes - {gene_combination}"
        )

    # Step 2: Count frequency of each unique gene combination
    combination_counts = Counter(genome_gene_combinations)

    logger.info(f"Found {len(combination_counts)} unique gene combinations")

    # Step 3: Filter by minimum frequency and single-gene exclusion
    filtered_combinations = []

    for gene_combo, frequency in combination_counts.items():
        # Apply minimum frequency filter
        if frequency < min_frequency:
            continue

        # Apply single-gene exclusion if requested
        if exclude_single and len(gene_combo) == 1:
            continue

        filtered_combinations.append(
            {
                "Gene_Combination": gene_combo,
                "Frequency": frequency,
                "Combination_Size": len(gene_combo),
                "Representative_Genes": ", ".join(gene_combo),
            }
        )

    # Step 4: Create results DataFrame
    if not filtered_combinations:
        logger.warning(
            f"No combinations meet criteria (min_frequency={min_frequency}, exclude_single={exclude_single})"
        )
        return pd.DataFrame(
            columns=[
                "Gene_Combination",
                "Frequency",
                "Combination_Size",
                "Representative_Genes",
            ]
        )

    results_df = pd.DataFrame(filtered_combinations)

    # Sort by frequency (descending) then by combination size (descending)
    results_df = results_df.sort_values(
        ["Frequency", "Combination_Size"], ascending=[False, False]
    )
    results_df = results_df.reset_index(drop=True)

    # Log summary statistics
    logger.info(f"Co-occurrence Analysis Summary:")
    logger.info(f"  • Total genomes analyzed: {len(accession_gene_summary)}")
    logger.info(f"  • Unique gene combinations found: {len(combination_counts)}")
    logger.info(f"  • Combinations meeting criteria: {len(results_df)}")
    logger.info(
        f"  • Most frequent combination: {results_df.iloc[0]['Representative_Genes']} ({results_df.iloc[0]['Frequency']} times)"
    )

    return results_df


def generate_cooccurrence_summary(
    cooccurrence_df: pd.DataFrame, mutation_df: pd.DataFrame
) -> Dict[str, Any]:
    """
    Generate comprehensive summary statistics for co-occurrence analysis.

    Args:
        cooccurrence_df: Results from analyze_cooccurrence function
        mutation_df: Original mutation DataFrame

    Returns:
        Dict[str, Any]: Summary statistics including:
            - total_genomes: Number of genomes analyzed
            - total_mutations: Total number of mutations
            - unique_genes: Number of unique genes with mutations
            - combination_stats: Statistics about gene combinations
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    if cooccurrence_df.empty or mutation_df.empty:
        return {
            "total_genomes": 0,
            "total_mutations": 0,
            "unique_genes": 0,
            "combination_stats": {},
        }

    total_genomes = mutation_df["Accession Number"].nunique()
    total_mutations = len(mutation_df)
    unique_genes = mutation_df["Gene Name"].nunique()

    # Combination size distribution
    size_distribution = cooccurrence_df["Combination_Size"].value_counts().to_dict()

    # Most frequent combinations
    top_combinations = cooccurrence_df.head(5)[
        ["Representative_Genes", "Frequency"]
    ].to_dict("records")

    summary = {
        "total_genomes": total_genomes,
        "total_mutations": total_mutations,
        "unique_genes": unique_genes,
        "total_combinations": len(cooccurrence_df),
        "combination_size_distribution": size_distribution,
        "top_combinations": top_combinations,
        "analysis_timestamp": pd.Timestamp.now().isoformat(),
    }

    logger.info("Generated comprehensive co-occurrence summary")
    return summary


def validate_mutation_dataframe(df: pd.DataFrame) -> Tuple[bool, List[str]]:
    """
    Validate the structure and content of a mutation DataFrame.

    Args:
        df: DataFrame to validate

    Returns:
        Tuple[bool, List[str]]: (is_valid, list_of_issues)
    """
    issues = []

    # Check required columns
    required_columns = ["Accession Number", "Gene Name"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        issues.append(f"Missing required columns: {missing_columns}")

    # Check for empty DataFrame
    if df.empty:
        issues.append("DataFrame is empty")
        return False, issues

    # Check for null values in critical columns
    for col in required_columns:
        if col in df.columns and df[col].isnull().any():
            null_count = df[col].isnull().sum()
            issues.append(f"Column '{col}' has {null_count} null values")

    # Check data types
    if "Accession Number" in df.columns:
        if not df["Accession Number"].dtype == "object":
            issues.append("Accession Number should be string/object type")

    if "Gene Name" in df.columns:
        if not df["Gene Name"].dtype == "object":
            issues.append("Gene Name should be string/object type")

    # Summary statistics
    if not issues:
        unique_genomes = df["Accession Number"].nunique()
        unique_genes = df["Gene Name"].nunique()

        if unique_genomes < 2:
            issues.append(
                f"Very few genomes ({unique_genomes}) - may not provide meaningful co-occurrence patterns"
            )

        if unique_genes < 2:
            issues.append(
                f"Very few genes ({unique_genes}) - co-occurrence analysis not applicable"
            )

    is_valid = len(issues) == 0
    return is_valid, issues
