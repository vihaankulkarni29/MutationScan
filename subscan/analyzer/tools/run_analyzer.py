#!/usr/bin/env python3
"""
SubScan Analyzer - Command Line Interface

High-performance standalone tool for mutation detection and co-occurrence analysis.
Part of the SubScan bioinformatics pipeline (Domino 5: The Analyzer).

Usage:
    python run_analyzer.py --manifest alignment_manifest.json --output-dir ./analysis_results

Features:
- Processes alignment manifest from WildTypeAligner (Domino 4)
- Performs comprehensive mutation detection across gene families
- Generates detailed CSV reports for downstream analysis
- Creates standardized manifest for pipeline integration
- Supports parallel processing for large genome collections

Scientific Capabilities:
- Reference-driven amino acid mutation detection
- Position-specific change cataloging
- Statistical co-occurrence pattern analysis
- Cross-genome mutation correlation studies
- Antibiotic resistance mutation profiling
"""

import os
import sys
import json
import argparse
import logging
import time
from pathlib import Path
from typing import List, Dict, Any
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# Add the src directory to the Python path for imports
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
)

from analyzer.engine import MutationAnalyzer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def analyze_single_gene_family(
    gene_family_data: Dict[str, Any], analyzer: MutationAnalyzer
) -> List[Dict[str, Any]]:
    """
    Worker function to analyze mutations in a single gene family.

    This function is designed for parallel execution to process multiple
    gene families simultaneously, maximizing computational efficiency.

    Args:
        gene_family_data: Dictionary containing gene family information
        analyzer: Configured MutationAnalyzer instance

    Returns:
        List of mutations detected in this gene family
    """
    gene_family = gene_family_data["gene_family"]
    alignment_file = gene_family_data["alignment_file"]

    logger.info(f"Starting mutation analysis for gene family: {gene_family}")

    try:
        # Analyze mutations in this alignment
        mutations = analyzer.identify_mutations(alignment_file)

        logger.info(
            f"Completed analysis for {gene_family}: {len(mutations)} mutations detected"
        )
        return mutations

    except Exception as e:
        logger.error(f"Error analyzing gene family {gene_family}: {str(e)}")
        return []


def execute_parallel_analysis(
    alignment_manifest: Dict[str, Any], output_dir: str, max_workers: int = None
) -> Dict[str, Any]:
    """
    Execute parallel mutation analysis across all gene families.

    This orchestrator function manages the entire analysis workflow:
    1. Validates input manifest and files
    2. Configures parallel processing workers
    3. Distributes gene families across workers
    4. Aggregates results from all analyses
    5. Generates comprehensive output reports

    Args:
        alignment_manifest: Parsed alignment manifest from Domino 4
        output_dir: Directory for output files
        max_workers: Maximum number of parallel workers

    Returns:
        Dictionary containing analysis results and statistics
    """
    start_time = time.time()

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory created: {output_dir}")

    # Extract gene family information from manifest
    gene_families = []

    # Handle different manifest formats
    if "alignment_results" in alignment_manifest:
        # Expected format from specification
        for result in alignment_manifest["alignment_results"]:
            gene_family = result["gene_family"]
            alignment_file = result["alignment_file"]

            # Validate alignment file exists
            if not os.path.exists(alignment_file):
                logger.warning(f"Alignment file not found: {alignment_file}")
                continue

            gene_families.append(
                {"gene_family": gene_family, "alignment_file": alignment_file}
            )

    elif "alignments" in alignment_manifest:
        # Actual format from WildTypeAligner
        for gene_family, result in alignment_manifest["alignments"].items():
            alignment_file = result["alignment_file_path"]

            # Validate alignment file exists
            if not os.path.exists(alignment_file):
                logger.warning(f"Alignment file not found: {alignment_file}")
                continue

            gene_families.append(
                {"gene_family": gene_family, "alignment_file": alignment_file}
            )

    else:
        logger.error(
            "Invalid manifest format: missing 'alignment_results' or 'alignments'"
        )
        return {"analysis_results": [], "error": "Invalid manifest format"}

    logger.info(f"Processing {len(gene_families)} gene families for mutation analysis")

    if not gene_families:
        logger.error("No valid gene families found for analysis")
        return {"analysis_results": [], "error": "No valid input files"}

    # Initialize analyzer
    analyzer = MutationAnalyzer(reference_pattern="reference", min_coverage=0.8)

    # Execute parallel analysis
    all_mutations = []
    analysis_results = []

    # Determine optimal number of workers
    if max_workers is None:
        max_workers = min(len(gene_families), os.cpu_count() or 4)

    logger.info(f"Starting parallel analysis with {max_workers} workers")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all gene family analysis tasks
        future_to_gene_family = {
            executor.submit(analyze_single_gene_family, gf_data, analyzer): gf_data
            for gf_data in gene_families
        }

        # Process completed tasks with progress bar
        with tqdm(total=len(gene_families), desc="Analyzing gene families") as pbar:
            for future in as_completed(future_to_gene_family):
                gene_family_data = future_to_gene_family[future]
                gene_family = gene_family_data["gene_family"]

                try:
                    mutations = future.result()
                    all_mutations.extend(mutations)

                    # Record analysis result
                    analysis_results.append(
                        {
                            "gene_family": gene_family,
                            "alignment_file": gene_family_data["alignment_file"],
                            "mutation_count": len(mutations),
                            "status": "completed",
                        }
                    )

                    logger.info(
                        f"Gene family {gene_family}: {len(mutations)} mutations"
                    )

                except Exception as e:
                    logger.error(
                        f"Failed to analyze gene family {gene_family}: {str(e)}"
                    )
                    analysis_results.append(
                        {
                            "gene_family": gene_family,
                            "alignment_file": gene_family_data["alignment_file"],
                            "mutation_count": 0,
                            "status": "failed",
                            "error": str(e),
                        }
                    )

                pbar.update(1)

    logger.info(
        f"Parallel analysis completed: {len(all_mutations)} total mutations detected"
    )

    # Perform co-occurrence analysis
    logger.info("Starting co-occurrence analysis...")
    cooccurrence_results = analyzer.analyze_cooccurrence(all_mutations)
    logger.info("Co-occurrence analysis completed")

    # Generate output reports
    mutation_report_path = os.path.join(output_dir, "mutation_report.csv")
    cooccurrence_report_path = os.path.join(output_dir, "cooccurrence_report.csv")

    # Save mutation report
    if all_mutations:
        mutation_df = pd.DataFrame(all_mutations)
        mutation_df.to_csv(mutation_report_path, index=False)
        logger.info(f"Mutation report saved: {mutation_report_path}")
    else:
        logger.warning("No mutations detected - mutation report not generated")

    # Save co-occurrence report
    if cooccurrence_results["mutation_pairs"]:
        cooccurrence_data = []
        for pair_name, pair_data in cooccurrence_results["mutation_pairs"].items():
            cooccurrence_data.append(pair_data)

        cooccurrence_df = pd.DataFrame(cooccurrence_data)
        cooccurrence_df.to_csv(cooccurrence_report_path, index=False)
        logger.info(f"Co-occurrence report saved: {cooccurrence_report_path}")
    else:
        logger.warning(
            "No co-occurrence patterns detected - co-occurrence report not generated"
        )

    # Calculate execution statistics
    execution_time = time.time() - start_time

    # Prepare final results
    final_results = {
        "analysis_results": analysis_results,
        "mutation_report": mutation_report_path if all_mutations else None,
        "cooccurrence_report": (
            cooccurrence_report_path if cooccurrence_results["mutation_pairs"] else None
        ),
        "statistics": {
            "total_gene_families": len(gene_families),
            "successful_analyses": sum(
                1 for r in analysis_results if r["status"] == "completed"
            ),
            "total_mutations": len(all_mutations),
            "execution_time_seconds": execution_time,
            "mutations_per_second": (
                len(all_mutations) / execution_time if execution_time > 0 else 0
            ),
            **cooccurrence_results["statistical_summary"],
        },
        "cooccurrence_data": cooccurrence_results,
    }

    logger.info(f"Analysis completed in {execution_time:.2f} seconds")
    logger.info(
        f"Processing rate: {final_results['statistics']['mutations_per_second']:.1f} mutations/second"
    )

    return final_results


def generate_analysis_manifest(
    input_manifest: Dict[str, Any], analysis_results: Dict[str, Any], output_dir: str
) -> str:
    """
    Generate standardized analysis manifest for pipeline integration.

    This function creates the JSON manifest that will be consumed by
    downstream reporting stages, following the SubScan pipeline standards.

    Args:
        input_manifest: Original alignment manifest from Domino 4
        analysis_results: Results from mutation analysis
        output_dir: Analysis output directory

    Returns:
        Path to the generated analysis manifest file
    """
    logger.info("Generating analysis manifest for pipeline integration")

    # Create analysis manifest based on input manifest
    analysis_manifest = {
        "pipeline_step": "analyzer",
        "version": "1.0.0",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input_manifest": input_manifest,
        "analysis_results": {
            "mutation_report": analysis_results.get("mutation_report"),
            "cooccurrence_report": analysis_results.get("cooccurrence_report"),
            "gene_family_results": analysis_results["analysis_results"],
        },
        "statistics": analysis_results["statistics"],
        "output_directory": output_dir,
        "next_stage": "reporting",
    }

    # Save analysis manifest
    manifest_path = os.path.join(output_dir, "analysis_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(analysis_manifest, f, indent=2)

    logger.info(f"Analysis manifest generated: {manifest_path}")
    return manifest_path


def parse_arguments():
    """Parse command line arguments with comprehensive help."""
    parser = argparse.ArgumentParser(
        description="SubScan Analyzer - Mutation Detection and Co-occurrence Analysis Engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python run_analyzer.py --manifest alignment_manifest.json --output-dir ./results
  
  # High-performance analysis with custom workers
  python run_analyzer.py --manifest alignment_manifest.json --output-dir ./results --max-workers 8
  
  # Analysis with verbose logging
  python run_analyzer.py --manifest alignment_manifest.json --output-dir ./results --verbose

Scientific Output:
  - mutation_report.csv: Comprehensive catalog of all detected mutations
  - cooccurrence_report.csv: Statistical analysis of mutation co-occurrence patterns
  - analysis_manifest.json: Pipeline-ready manifest for downstream processing

Pipeline Integration:
  This tool processes alignment_manifest.json from WildTypeAligner (Domino 4)
  and generates analysis_manifest.json for the reporting stage.
        """,
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to alignment manifest JSON file from WildTypeAligner (Domino 4)",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to store analysis results and reports",
    )

    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Maximum number of parallel workers (default: auto-detect based on CPU cores)",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging for detailed progress tracking",
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.8,
        help="Minimum sequence coverage threshold (0.0-1.0, default: 0.8)",
    )

    return parser.parse_args()


def main():
    """Main entry point for the SubScan Analyzer CLI."""
    print("🧬 SubScan Analyzer - Domino 5: The Analyzer")
    print("Mutation Detection and Co-occurrence Analysis Engine")
    print("=" * 60)

    # Parse command line arguments
    args = parse_arguments()

    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.info("Verbose logging enabled")

    # Validate inputs
    if not os.path.exists(args.manifest):
        logger.error(f"Manifest file not found: {args.manifest}")
        return 1

    if args.min_coverage < 0.0 or args.min_coverage > 1.0:
        logger.error("Minimum coverage must be between 0.0 and 1.0")
        return 1

    try:
        # Load alignment manifest
        logger.info(f"Loading alignment manifest: {args.manifest}")
        with open(args.manifest, "r") as f:
            alignment_manifest = json.load(f)

        # Validate manifest format
        if (
            "alignment_results" not in alignment_manifest
            and "alignments" not in alignment_manifest
        ):
            logger.error(
                "Invalid manifest format: missing 'alignment_results' or 'alignments'"
            )
            return 1

        num_gene_families = len(
            alignment_manifest.get(
                "alignment_results", alignment_manifest.get("alignments", {})
            )
        )
        logger.info(f"Manifest loaded successfully: {num_gene_families} gene families")

        # Execute mutation analysis
        logger.info("Starting mutation analysis pipeline...")
        analysis_results = execute_parallel_analysis(
            alignment_manifest, args.output_dir, max_workers=args.max_workers
        )

        # Generate pipeline manifest
        manifest_path = generate_analysis_manifest(
            alignment_manifest, analysis_results, args.output_dir
        )

        # Print summary
        stats = analysis_results["statistics"]
        print("\n" + "=" * 60)
        print("📊 ANALYSIS SUMMARY")
        print("=" * 60)
        print(f"Gene families processed: {stats['total_gene_families']}")
        print(f"Successful analyses: {stats['successful_analyses']}")
        print(f"Total mutations detected: {stats['total_mutations']}")
        print(f"Execution time: {stats['execution_time_seconds']:.2f} seconds")
        print(f"Processing rate: {stats['mutations_per_second']:.1f} mutations/second")
        print(f"Co-occurrence pairs: {stats.get('mutation_pair_count', 0)}")
        print("\n📁 OUTPUT FILES:")
        if analysis_results.get("mutation_report"):
            print(f"  • Mutation report: {analysis_results['mutation_report']}")
        if analysis_results.get("cooccurrence_report"):
            print(
                f"  • Co-occurrence report: {analysis_results['cooccurrence_report']}"
            )
        print(f"  • Analysis manifest: {manifest_path}")
        print("\n✅ SubScan Analyzer completed successfully!")

        return 0

    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
