#!/usr/bin/env python3
"""
SubScan Pipeline - Domino 6: The Co-occurrence Analyzer

This tool analyzes mutation co-occurrence patterns across genomes to identify
functionally linked mutations in antimicrobial resistance genes.

Purpose:
- Reads analysis_manifest.json from Domino 5 (The Analyzer)
- Loads the subscan_mutation_report.xlsx master mutation list
- Performs statistical analysis of co-occurring mutations in different genes
- Generates cooccurrence_manifest.json for Domino 7 (The Reporter)

Usage:
    python run_cooccurrence_analyzer.py --manifest analysis_manifest.json --output-dir ./cooccurrence_results

Pipeline Flow:
    Domino 5 (Analyzer) → analysis_manifest.json → Domino 6 (Co-occurrence) → cooccurrence_manifest.json → Domino 7 (Reporter)

Scientific Value:
    Identifies patterns where mutations in multiple genes co-occur in the same isolate,
    providing insights into potential functional relationships between resistance mechanisms.

Author: MutationScan Development Team
"""

import argparse
import os
import sys
import json
from typing import Dict, Any, Tuple
import logging
from datetime import datetime

import pandas as pd

# Add the src directory to Python path for imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "src"))

# pylint: disable=wrong-import-position,import-error
from subscan.cooccurrence import (
    analyze_cooccurrence,
    generate_cooccurrence_summary,
    validate_mutation_dataframe,
)


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure the command-line argument parser.

    Returns:
        argparse.ArgumentParser: Configured argument parser for co-occurrence analysis
    """
    parser = argparse.ArgumentParser(
        description="SubScan Domino 6: Co-occurrence Analyzer - Mutation Pattern Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic co-occurrence analysis
  python run_cooccurrence_analyzer.py --manifest analysis_manifest.json --output-dir ./cooccurrence_results
  
  # Analysis with custom output directory
  python run_cooccurrence_analyzer.py --manifest results/analysis_manifest.json --output-dir ./patterns
  
  # Full pipeline integration
  python run_cooccurrence_analyzer.py --manifest analyzer_output/analysis_manifest.json --output-dir ./stage6_cooccurrence

Pipeline Integration:
  This tool processes analysis_manifest.json from the Mutation Analyzer (Domino 5)
  and generates cooccurrence_manifest.json for the HTML Reporter (Domino 7).

Scientific Output:
  - cooccurrence_report.csv: Statistical analysis of mutation co-occurrence patterns
  - cooccurrence_manifest.json: Updated manifest with all pipeline results

Input Requirements:
  - analysis_manifest.json: Contains path to subscan_mutation_report.xlsx
  - subscan_mutation_report.xlsx: Master mutation list from all analyzed genomes
  - All referenced files must be accessible at their manifest paths

Co-occurrence Analysis:
  For each genome isolate, identifies combinations of genes with mutations:
  • Single gene mutations: ('acrA',)
  • Two-gene combinations: ('acrA', 'acrB') 
  • Multi-gene patterns: ('acrA', 'acrB', 'tolC')
  
  Statistical summary provides frequency counts for each unique combination,
  helping researchers identify potential functional linkages between resistance genes.
        """,
    )

    # Required arguments
    required = parser.add_argument_group("Required Arguments")

    required.add_argument(
        "--manifest",
        required=True,
        help="Path to analysis_manifest.json file from the Mutation Analyzer (Domino 5)",
        metavar="FILE",
    )

    required.add_argument(
        "--output-dir",
        required=True,
        help="Directory where co-occurrence analysis results will be saved",
        metavar="DIR",
    )

    # Optional arguments
    optional = parser.add_argument_group("Optional Arguments")

    optional.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging for detailed progress tracking",
    )

    optional.add_argument(
        "--min-frequency",
        type=int,
        default=2,
        help="Minimum frequency threshold for reporting co-occurrence patterns (default: 2)",
        metavar="N",
    )

    optional.add_argument(
        "--exclude-single",
        action="store_true",
        help="Exclude single-gene mutations from co-occurrence analysis (focus on multi-gene patterns)",
    )

    return parser


def validate_inputs(args) -> None:
    """
    Validate input files and directories exist and are accessible.

    Args:
        args: Parsed command-line arguments

    Raises:
        FileNotFoundError: If manifest file is missing
        ValueError: If manifest format is invalid
        OSError: If output directory cannot be created
    """
    # Import shared utilities
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "src"))
    from subscan.utils import validate_manifest_file, validate_output_directory
    
    # Use shared validation functions
    validate_manifest_file(args.manifest)
    validate_output_directory(args.output_dir, create=True)
    
    print(f"📁 Output directory ready: {args.output_dir}")
    print("✅ Input validation completed successfully")


def setup_logging(output_dir: str, verbose: bool = False) -> logging.Logger:
    """
    Configure logging for co-occurrence analysis.

    Args:
        output_dir: Directory for log file output
        verbose: Enable verbose console output

    Returns:
        logging.Logger: Configured logger instance
    """
    log_level = logging.DEBUG if verbose else logging.INFO

    # Create logger
    logger = logging.getLogger("cooccurrence_analyzer")
    logger.setLevel(log_level)

    # Avoid duplicate handlers
    if logger.handlers:
        logger.handlers.clear()

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_format = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

    # File handler
    log_file = os.path.join(output_dir, "cooccurrence_analysis.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s"
    )
    file_handler.setFormatter(file_format)
    logger.addHandler(file_handler)

    return logger


def load_analysis_manifest(manifest_path: str) -> Dict[str, Any]:
    """
    Load and validate the analysis manifest from Domino 5.

    Args:
        manifest_path: Path to analysis_manifest.json

    Returns:
        Dict[str, Any]: Parsed manifest data

    Raises:
        FileNotFoundError: If manifest file doesn't exist
        ValueError: If manifest format is invalid
        KeyError: If required manifest keys are missing
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest_data = json.load(f)

        logger.info(f"Loaded analysis manifest from: {manifest_path}")

        # Validate manifest structure
        required_keys = ["analysis_results", "pipeline_step"]
        missing_keys = [key for key in required_keys if key not in manifest_data]

        if missing_keys:
            raise KeyError(f"Missing required manifest keys: {missing_keys}")

        # Extract mutation report path
        analysis_results = manifest_data.get("analysis_results", {})
        if "mutation_report_path" not in analysis_results:
            raise KeyError("Missing 'mutation_report_path' in analysis_results")

        mutation_report_path = analysis_results["mutation_report_path"]

        # Validate mutation report file exists
        if not os.path.isfile(mutation_report_path):
            raise FileNotFoundError(
                f"Mutation report file not found: {mutation_report_path}"
            )

        logger.info(f"Found mutation report: {mutation_report_path}")
        return manifest_data

    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in manifest file: {e}") from e


def load_mutation_report(report_path: str) -> pd.DataFrame:
    """
    Load the mutation report Excel file from Domino 5.

    Args:
        report_path: Path to subscan_mutation_report.xlsx

    Returns:
        pd.DataFrame: Mutation data with standardized columns

    Raises:
        FileNotFoundError: If report file doesn't exist
        ValueError: If file format is invalid or missing required columns
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    try:
        # Load Excel file
        logger.info(f"Loading mutation report from: {report_path}")
        mutation_df = pd.read_excel(report_path)

        logger.info(f"Loaded {len(mutation_df)} mutation records")

        # Validate DataFrame structure
        is_valid, issues = validate_mutation_dataframe(mutation_df)

        if not is_valid:
            error_msg = f"Invalid mutation report format: {'; '.join(issues)}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Log summary statistics
        unique_genomes = mutation_df["Accession Number"].nunique()
        unique_genes = mutation_df["Gene Name"].nunique()

        logger.info("Mutation report summary:")
        logger.info(f"  • Total mutations: {len(mutation_df)}")
        logger.info(f"  • Unique genomes: {unique_genomes}")
        logger.info(f"  • Unique genes: {unique_genes}")

        return mutation_df

    except Exception as e:
        logger.error(f"Failed to load mutation report: {e}")
        raise


def save_cooccurrence_report(
    cooccurrence_df: pd.DataFrame, summary_stats: Dict[str, Any], output_dir: str
) -> str:
    """
    Save co-occurrence analysis results to CSV file.

    Args:
        cooccurrence_df: Co-occurrence analysis results
        summary_stats: Summary statistics dictionary
        output_dir: Output directory path

    Returns:
        str: Path to saved CSV file
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    # Save main results
    csv_path = os.path.join(output_dir, "cooccurrence_report.csv")
    cooccurrence_df.to_csv(csv_path, index=False)

    logger.info(f"Saved co-occurrence report: {csv_path}")

    # Save summary statistics as JSON
    summary_path = os.path.join(output_dir, "cooccurrence_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary_stats, f, indent=2)

    logger.info(f"Saved summary statistics: {summary_path}")

    return csv_path


def generate_cooccurrence_manifest(
    input_manifest: Dict[str, Any],
    cooccurrence_report_path: str,
    summary_stats: Dict[str, Any],
    output_dir: str,
) -> str:
    """
    Generate the standardized cooccurrence_manifest.json for pipeline integration.

    This function creates the final manifest file that connects Domino 6 to Domino 7,
    preserving all previous analysis results while adding co-occurrence data.

    Args:
        input_manifest: Original analysis_manifest.json from Domino 5
        cooccurrence_report_path: Path to generated cooccurrence_report.csv
        summary_stats: Analysis summary statistics
        output_dir: Output directory for manifest file

    Returns:
        str: Path to generated cooccurrence_manifest.json
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    # Create new manifest based on input manifest
    cooccurrence_manifest = input_manifest.copy()

    # Update pipeline step
    cooccurrence_manifest["pipeline_step"] = "CooccurrenceAnalyzer"
    cooccurrence_manifest["domino_number"] = 6
    cooccurrence_manifest["domino_name"] = "Co-occurrence Analyzer"

    # Update timestamp
    cooccurrence_manifest["timestamp"] = datetime.now().isoformat()

    # Add co-occurrence analysis results to analysis_results section
    if "analysis_results" not in cooccurrence_manifest:
        cooccurrence_manifest["analysis_results"] = {}

    # Add co-occurrence specific results
    cooccurrence_manifest["analysis_results"]["cooccurrence_report_path"] = (
        os.path.abspath(cooccurrence_report_path)
    )
    cooccurrence_manifest["analysis_results"]["cooccurrence_summary"] = summary_stats

    # Add summary path if it exists
    summary_path = os.path.join(output_dir, "cooccurrence_summary.json")
    if os.path.isfile(summary_path):
        cooccurrence_manifest["analysis_results"]["cooccurrence_summary_path"] = (
            os.path.abspath(summary_path)
        )

    # Update processing history
    if "processing_history" not in cooccurrence_manifest:
        cooccurrence_manifest["processing_history"] = []

    cooccurrence_manifest["processing_history"].append(
        {
            "step": "CooccurrenceAnalyzer",
            "domino": 6,
            "timestamp": datetime.now().isoformat(),
            "input_files": [
                input_manifest.get("analysis_results", {}).get(
                    "mutation_report_path", "unknown"
                )
            ],
            "output_files": [os.path.abspath(cooccurrence_report_path)],
            "patterns_found": summary_stats.get("total_combinations", 0),
            "genomes_analyzed": summary_stats.get("total_genomes", 0),
        }
    )

    # Add pipeline status
    cooccurrence_manifest["pipeline_status"] = {
        "current_step": "CooccurrenceAnalyzer",
        "completed_steps": [
            "Harvester",
            "Annotator",
            "Extractor",
            "Aligner",
            "Analyzer",
            "CooccurrenceAnalyzer",
        ],
        "next_step": "Reporter",
        "progress_percentage": 85.7,  # 6/7 steps complete
    }

    # Save manifest to JSON file
    manifest_path = os.path.join(output_dir, "cooccurrence_manifest.json")

    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(cooccurrence_manifest, f, indent=2, ensure_ascii=False)

    logger.info(f"Generated cooccurrence manifest: {manifest_path}")

    # Log manifest summary
    logger.info("Manifest contents summary:")
    logger.info(f"  • Pipeline step: {cooccurrence_manifest['pipeline_step']}")
    logger.info(
        f"  • Total analysis results keys: {len(cooccurrence_manifest.get('analysis_results', {}))}"
    )
    logger.info(
        f"  • Processing history entries: {len(cooccurrence_manifest.get('processing_history', []))}"
    )
    logger.info(
        f"  • Co-occurrence patterns: {summary_stats.get('total_combinations', 0)}"
    )

    return manifest_path


def validate_pipeline_continuity(
    input_manifest: Dict[str, Any], output_manifest_path: str
) -> bool:
    """
    Validate that the generated manifest maintains pipeline continuity.

    Args:
        input_manifest: Original analysis manifest
        output_manifest_path: Path to generated cooccurrence manifest

    Returns:
        bool: True if validation passes
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    try:
        # Load generated manifest
        with open(output_manifest_path, "r", encoding="utf-8") as f:
            output_manifest = json.load(f)

        # Validation checks
        checks = []

        # 1. Pipeline step updated
        checks.append(output_manifest.get("pipeline_step") == "CooccurrenceAnalyzer")

        # 2. Analysis results preserved and extended
        input_results = input_manifest.get("analysis_results", {})
        output_results = output_manifest.get("analysis_results", {})

        # All input keys should be preserved
        for key in input_results.keys():
            checks.append(key in output_results)

        # New co-occurrence key should be added
        checks.append("cooccurrence_report_path" in output_results)

        # 3. Proper file paths
        cooccurrence_path = output_results.get("cooccurrence_report_path", "")
        checks.append(os.path.isfile(cooccurrence_path))

        # 4. Processing history updated
        history = output_manifest.get("processing_history", [])
        checks.append(len(history) > 0)
        checks.append(
            any(entry.get("step") == "CooccurrenceAnalyzer" for entry in history)
        )

        # 5. Pipeline status
        status = output_manifest.get("pipeline_status", {})
        checks.append(status.get("current_step") == "CooccurrenceAnalyzer")
        checks.append("CooccurrenceAnalyzer" in status.get("completed_steps", []))

        validation_passed = all(checks)

        if validation_passed:
            logger.info("✅ Pipeline continuity validation passed")
        else:
            logger.warning(f"⚠️ Pipeline validation issues detected: {checks}")

        return validation_passed

    except (FileNotFoundError, json.JSONDecodeError, KeyError, ValueError) as e:
        logger.error(f"Manifest validation failed: {e}")
        return False


def run_cooccurrence_analysis(args) -> Tuple[str, Dict[str, Any], str]:
    """
    Execute the complete co-occurrence analysis workflow.

    Args:
        args: Parsed command-line arguments

    Returns:
        Tuple[str, Dict[str, Any], str]: (csv_report_path, analysis_summary, manifest_path)
    """
    logger = logging.getLogger("cooccurrence_analyzer")

    # Step 1: Load analysis manifest
    logger.info("🔍 Loading analysis manifest...")
    manifest_data = load_analysis_manifest(args.manifest)

    # Step 2: Load mutation report
    logger.info("📊 Loading mutation report...")
    mutation_report_path = manifest_data["analysis_results"]["mutation_report_path"]
    mutation_df = load_mutation_report(mutation_report_path)

    # Step 3: Perform co-occurrence analysis
    logger.info("🔬 Analyzing co-occurrence patterns...")
    cooccurrence_df = analyze_cooccurrence(
        mutation_df,
        min_frequency=args.min_frequency,
        exclude_single=args.exclude_single,
    )

    if cooccurrence_df.empty:
        logger.warning("No co-occurrence patterns found with current criteria")
        print("⚠️  No patterns found - consider lowering minimum frequency threshold")
    else:
        logger.info(f"Found {len(cooccurrence_df)} co-occurrence patterns")
        print(f"✅ Identified {len(cooccurrence_df)} co-occurrence patterns")

    # Step 4: Generate summary statistics
    logger.info("📈 Generating summary statistics...")
    summary_stats = generate_cooccurrence_summary(cooccurrence_df, mutation_df)

    # Step 5: Save results
    logger.info("💾 Saving analysis results...")
    csv_report_path = save_cooccurrence_report(
        cooccurrence_df, summary_stats, args.output_dir
    )

    # Step 6: Generate pipeline manifest
    logger.info("📄 Generating cooccurrence manifest...")
    manifest_path = generate_cooccurrence_manifest(
        manifest_data, csv_report_path, summary_stats, args.output_dir
    )

    # Step 7: Validate pipeline continuity
    logger.info("🔍 Validating pipeline continuity...")
    validation_passed = validate_pipeline_continuity(manifest_data, manifest_path)

    if not validation_passed:
        logger.warning("Pipeline validation detected potential issues")

    return csv_report_path, summary_stats, manifest_path


if __name__ == "__main__":
    # Parse command-line arguments
    main_parser = create_argument_parser()
    main_args = main_parser.parse_args()

    print("🔬 SubScan Domino 6: Co-occurrence Analyzer")
    print("=" * 55)
    print(f"📅 Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("📊 Configuration:")
    print(f"   • Input manifest: {main_args.manifest}")
    print(f"   • Output directory: {main_args.output_dir}")
    print(f"   • Minimum frequency: {main_args.min_frequency}")
    print(f"   • Exclude singles: {main_args.exclude_single}")
    print(f"   • Verbose mode: {main_args.verbose}")
    print()

    try:
        # Validate inputs
        print("🔍 Validating inputs...")
        validate_inputs(main_args)

        # Setup logging
        main_logger = setup_logging(main_args.output_dir, main_args.verbose)
        main_logger.info("Co-occurrence Analyzer initialized successfully")

        print("✅ CLI initialization completed!")
        print("🚀 Beginning co-occurrence pattern analysis...")
        print()

        # Execute the complete analysis workflow
        main_csv_report_path, main_summary_stats, main_manifest_path = (
            run_cooccurrence_analysis(main_args)
        )

        print("\n" + "=" * 25)
        print("*** DOMINO 6 CO-OCCURRENCE ANALYSIS COMPLETED! ***")
        print("=" * 25)
        print()
        print("📊 ANALYSIS SUMMARY:")
        print("=" * 55)
        print(
            f"⚡ Total genomes analyzed:        {main_summary_stats['total_genomes']}"
        )
        print(
            f"🧬 Total mutations processed:     {main_summary_stats['total_mutations']}"
        )
        print(f"🎯 Unique genes identified:       {main_summary_stats['unique_genes']}")
        print(
            f"🔗 Co-occurrence patterns found:  {main_summary_stats['total_combinations']}"
        )
        print()

        if main_summary_stats["top_combinations"]:
            print("🔬 TOP CO-OCCURRENCE PATTERNS:")
            print("-" * 35)
            for i, combo in enumerate(main_summary_stats["top_combinations"][:5], 1):
                print(
                    f"  {i}. {combo['Representative_Genes']} → {combo['Frequency']} isolates"
                )
        print()

        print("OUTPUT FILES:")
        print("=" * 55)
        print(f"Co-occurrence Report:     {os.path.basename(main_csv_report_path)}")
        print("Summary Statistics:       cooccurrence_summary.json")
        print("Pipeline Manifest:        cooccurrence_manifest.json")
        print()

        print("🚀 PIPELINE INTEGRATION:")
        print("=" * 55)
        print("✅ Input:  analysis_manifest.json (from Domino 5 - Analyzer)")
        print("✅ Output: cooccurrence_manifest.json (ready for Domino 7 - Reporter)")
        print("✅ Status: Pipeline continuity maintained")
        print()

        print("🎯 SCIENTIFIC INSIGHTS:")
        print("-" * 25)
        print("• Identified mutation patterns across resistance genes")
        print("• Quantified co-occurrence frequencies in isolate populations")
        print("• Ready for interactive visualization in HTML dashboard")
        print()

        print(f"📁 All results saved to: {main_args.output_dir}")
        print("🔬 Co-occurrence analysis complete - ready for final reporting!")

        print("\n🧬" + "=" * 54 + "🧬")
        print("  DOMINO 6 SUCCESSFULLY BRIDGES ANALYZER → REPORTER")
        print("🧬" + "=" * 54 + "🧬")

    except (FileNotFoundError, ValueError, OSError, KeyError) as e:
        print(f"❌ Analysis Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"❌ Unexpected Error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(2)
