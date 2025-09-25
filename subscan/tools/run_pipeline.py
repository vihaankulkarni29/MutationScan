#!/usr/bin/env python3
"""
MutationScan Pipeline Orchestrator

This is the master entry point for the complete MutationScan bioinformatics pipeline.
It orchestrates the execution of all seven "Domino" tools in the correct sequence,
managing data handoffs via manifest files to provide a seamless end-to-end experience.

Pipeline Flow:
    Domino 1 (Harvester) → genome_manifest.json →
    Domino 2 (Annotator) → annotation_manifest.json →
    Domino 3 (Extractor) → protein_manifest.json →
    Domino 4 (Aligner) → alignment_manifest.json →
    Domino 5 (Analyzer) → analysis_manifest.json →
    Domino 6 (Co-occurrence) → cooccurrence_manifest.json →
    Domino 7 (Reporter) → Interactive HTML Dashboard

Usage:
    python run_pipeline.py --accessions accessions.txt --gene-list genes.txt --email user@domain.com --output-dir ./results --sepi-species "Escherichia coli"

Author: MutationScan Development Team
"""

import argparse
import os
import sys
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import List


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure the master command-line argument parser.

    This parser collects all parameters needed for the complete pipeline run,
    providing a unified interface for users to configure the entire workflow.

    Returns:
        argparse.ArgumentParser: Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="MutationScan Pipeline Orchestrator - Complete AMR Analysis Workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic pipeline run with SEPI reference mode
  python run_pipeline.py --accessions research_accessions.txt --gene-list resistance_genes.txt --email researcher@university.edu --output-dir ./pipeline_results --sepi-species "Escherichia coli"
  
  # High-performance run with user-provided references
  python run_pipeline.py --accessions genome_list.txt --gene-list target_genes.txt --email user@domain.com --output-dir ./results --user-reference-dir ./my_references --threads 16
  
  # Quick analysis with default settings
  python run_pipeline.py --accessions short_list.txt --gene-list core_genes.txt --email analyst@lab.org --output-dir ./quick_run --sepi-species "Escherichia coli" --threads 8

Pipeline Stages:
  [1/7] Harvester    - Extract genomes from NCBI database
  [2/7] Annotator    - Identify AMR genes using CARD database
  [3/7] Extractor    - Extract protein sequences for target genes
  [4/7] Aligner      - Align sequences to wild-type references
  [5/7] Analyzer     - Detect mutations and variants
  [6/7] Co-occurrence- Analyze mutation co-occurrence patterns
  [7/7] Reporter     - Generate interactive HTML dashboard

Output Structure:
  <output-dir>/
  └── mutationscan_run_YYYYMMDD_HHMMSS/
      ├── 01_harvester_results/
      ├── 02_annotator_results/
      ├── 03_extractor_results/
      ├── 04_aligner_results/
      ├── 05_analyzer_results/
      ├── 06_cooccurrence_results/
      ├── 07_reporter_results/
      └── final_report.html          <- Main result file
        """,
    )

    # Required input arguments
    required = parser.add_argument_group("Required Arguments")

    required.add_argument(
        "--accessions",
        required=True,
        help="Path to text file containing NCBI accession numbers (one per line)",
        metavar="FILE",
    )

    required.add_argument(
        "--gene-list",
        required=True,
        help="Path to text file containing target gene names for analysis (one per line)",
        metavar="FILE",
    )

    required.add_argument(
        "--email",
        required=True,
        help="Valid email address for NCBI API access (required by NCBI guidelines)",
        metavar="EMAIL",
    )

    required.add_argument(
        "--output-dir",
        required=True,
        help="Base directory where all pipeline results will be stored",
        metavar="DIR",
    )

    # Reference mode (mutually exclusive group)
    reference_group = parser.add_mutually_exclusive_group(required=True)

    reference_group.add_argument(
        "--user-reference-dir",
        help="Path to directory containing user-provided wild-type reference sequences",
        metavar="DIR",
    )

    reference_group.add_argument(
        "--sepi-species",
        help="Species name for SEPI wild-type reference lookup (e.g., 'Escherichia coli')",
        metavar="SPECIES",
    )

    # Performance and optional arguments
    optional = parser.add_argument_group("Performance & Optional")

    optional.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of CPU threads to use for parallel processing (default: 4)",
        metavar="N",
    )

    optional.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging for detailed progress information",
    )

    optional.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep all intermediate files for debugging (increases disk usage)",
    )

    optional.add_argument(
        "--open-report",
        action="store_true",
        help="Automatically open the final HTML report in default browser",
    )

    return parser


def validate_input_files(args) -> None:
    """
    Validate that all required input files exist and are readable.

    Args:
        args: Parsed command-line arguments

    Raises:
        FileNotFoundError: If any required file is missing
        ValueError: If file format is invalid
    """
    # Check accessions file
    if not os.path.isfile(args.accessions):
        raise FileNotFoundError(f"Accessions file not found: {args.accessions}")

    # Check gene list file
    if not os.path.isfile(args.gene_list):
        raise FileNotFoundError(f"Gene list file not found: {args.gene_list}")

    # Check user reference directory if specified
    if args.user_reference_dir and not os.path.isdir(args.user_reference_dir):
        raise FileNotFoundError(
            f"User reference directory not found: {args.user_reference_dir}"
        )

    # Validate email format (basic check)
    if "@" not in args.email or "." not in args.email:
        raise ValueError(f"Invalid email format: {args.email}")

    print("✅ Input file validation completed successfully")


def create_run_directory(base_output_dir: str) -> str:
    """
    Create a timestamped subdirectory for this pipeline run.

    Args:
        base_output_dir: Base output directory path

    Returns:
        str: Path to the created run-specific directory
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(base_output_dir, f"mutationscan_run_{timestamp}")

    os.makedirs(run_dir, exist_ok=True)
    print(f"📁 Created pipeline run directory: {run_dir}")

    return run_dir


def find_manifest_file(output_dir: str, manifest_pattern: str) -> str:
    """
    Find the output manifest file from a domino tool execution.

    Args:
        output_dir: Directory where domino tool wrote its results
        manifest_pattern: Expected manifest filename pattern

    Returns:
        str: Full path to the found manifest file

    Raises:
        FileNotFoundError: If manifest file cannot be found
    """
    manifest_path = os.path.join(output_dir, manifest_pattern)

    if not os.path.isfile(manifest_path):
        raise FileNotFoundError(f"Expected manifest file not found: {manifest_path}")

    return manifest_path


def execute_domino_tool(
    tool_script: str,
    tool_args: List[str],
    domino_name: str,
    step_number: int,
    total_steps: int,
    verbose: bool = False,
) -> subprocess.CompletedProcess:
    """
    Execute a single domino tool with proper error handling and enhanced logging.

    Args:
        tool_script: Path to the domino tool script
        tool_args: List of command-line arguments for the tool
        domino_name: Human-readable name of the domino tool
        step_number: Current step number (1-based)
        total_steps: Total number of steps in pipeline
        verbose: Enable verbose subprocess output

    Returns:
        subprocess.CompletedProcess: Result of the subprocess execution

    Raises:
        subprocess.CalledProcessError: If the domino tool fails
    """
    # Enhanced status message with progress bar visualization
    progress_bar = "█" * step_number + "░" * (total_steps - step_number)
    percentage = int((step_number / total_steps) * 100)

    print(f"\n🔧 [{step_number}/{total_steps}] {domino_name}")
    print(f"   Progress: [{progress_bar}] {percentage}%")

    if verbose:
        print(f"   Command: python {tool_script} {' '.join(tool_args)}")

    # Show estimated time for each step
    step_start_time = time.time()

    # Construct full command
    cmd = [sys.executable, tool_script] + tool_args

    # Execute with proper error handling
    env = os.environ.copy()
    env['PYTHONIOENCODING'] = 'utf-8'
    
    result = subprocess.run(
        cmd,
        capture_output=not verbose,  # Show output in verbose mode
        text=True,
        encoding='utf-8',
        errors='replace',
        env=env,
        check=True,  # Raises CalledProcessError on non-zero exit
    )

    # Calculate step execution time
    step_duration = time.time() - step_start_time
    minutes = int(step_duration // 60)
    seconds = int(step_duration % 60)

    print(f"   ✅ Completed in {minutes:02d}:{seconds:02d}")
    return result


def run_pipeline_sequence(args, run_directory: str) -> str:
    """
    Execute the complete pipeline sequence of domino tools.

    This function orchestrates all seven domino tools in the correct order,
    managing manifest file handoffs between each stage.

    Args:
        args: Parsed command-line arguments
        run_directory: Directory for this specific pipeline run

    Returns:
        str: Path to the final HTML report

    Raises:
        subprocess.CalledProcessError: If any domino tool fails
        FileNotFoundError: If expected manifest files are missing
    """
    print("\n" + "🧬" * 25)
    print("🚀 STARTING MUTATIONSCAN PIPELINE EXECUTION")
    print("🧬" * 25)
    print("📋 Pipeline Overview:")
    print("   🔸 Step 1: Harvester    - Extract genomes from NCBI")
    print("   🔸 Step 2: Annotator    - Identify AMR genes with CARD")
    print("   🔸 Step 3: Extractor    - Extract protein sequences")
    print("   🔸 Step 4: Aligner      - Align to wild-type references")
    print("   🔸 Step 5: Analyzer     - Detect mutations and variants")
    print("   🔸 Step 6: Co-occurrence- Analyze mutation patterns")
    print("   🔸 Step 7: Reporter     - Generate interactive dashboard")
    print()

    # Get the tools directory path (relative to this script)
    tools_dir = Path(__file__).parent
    current_manifest = None
    pipeline_start_time = time.time()

    # Step 1: Domino 1 - Harvester (NCBI Genome Extraction)
    print("🧬 PHASE 1: GENOME ACQUISITION")
    print("-" * 35)
    harvester_output = os.path.join(run_directory, "01_harvester_results")
    harvester_args = [
        "--accessions",
        args.accessions,
        "--email",
        args.email,
        "--output-dir",
        harvester_output,
    ]

    execute_domino_tool(
        str(tools_dir / "run_harvester.py"),
        harvester_args,
        "Harvester - NCBI Genome Extraction",
        1,
        7,
        args.verbose,
    )

    # Find genome manifest for next step
    current_manifest = find_manifest_file(harvester_output, "genome_manifest.json")
    print(f"   📄 Manifest Generated: {os.path.basename(current_manifest)}")

    # Step 2: Domino 2 - Annotator (AMR Gene Annotation)
    print("\n🦠 PHASE 2: AMR GENE IDENTIFICATION")
    print("-" * 35)
    annotator_output = os.path.join(run_directory, "02_annotator_results")
    annotator_args = [
        "--manifest",
        current_manifest,
        "--output-dir",
        annotator_output,
        "--threads",
        str(args.threads),
    ]

    execute_domino_tool(
        str(tools_dir / "run_annotator.py"),
        annotator_args,
        "Annotator - AMR Gene Detection (CARD)",
        2,
        7,
        args.verbose,
    )

    # Find annotation manifest for next step
    current_manifest = find_manifest_file(annotator_output, "annotation_manifest.json")
    print(f"   📄 Manifest Generated: {os.path.basename(current_manifest)}")

    # Step 3: Domino 3 - Extractor (Protein Sequence Extraction)
    print("\n🧪 PHASE 3: PROTEIN SEQUENCE EXTRACTION")
    print("-" * 40)
    extractor_output = os.path.join(run_directory, "03_extractor_results")
    extractor_args = [
        "--manifest",
        current_manifest,
        "--gene-list",
        args.gene_list,
        "--output-dir",
        extractor_output,
        "--threads",
        str(args.threads),
    ]

    execute_domino_tool(
        str(tools_dir / "run_extractor.py"),
        extractor_args,
        "Extractor - Protein Sequence Extraction",
        3,
        7,
        args.verbose,
    )

    # Find protein manifest for next step
    current_manifest = find_manifest_file(extractor_output, "protein_manifest.json")
    print(f"   📄 Manifest Generated: {os.path.basename(current_manifest)}")

    # Step 4: Domino 4 - Aligner (Wild-Type Sequence Alignment)
    print("\n🎯 PHASE 4: REFERENCE ALIGNMENT")
    print("-" * 30)
    aligner_output = os.path.join(run_directory, "04_aligner_results")
    aligner_args = [
        "--manifest",
        current_manifest,
        "--output-dir",
        aligner_output,
        "--threads",
        str(args.threads),
    ]

    # Add reference mode arguments
    if args.sepi_species:
        aligner_args.extend(["--sepi-species", args.sepi_species])
        print(f"   🧬 Reference Mode: SEPI Database ({args.sepi_species})")
    else:
        aligner_args.extend(["--user-reference-dir", args.user_reference_dir])
        print(f"   🧬 Reference Mode: User-Provided ({args.user_reference_dir})")

    execute_domino_tool(
        str(tools_dir / "run_aligner.py"),
        aligner_args,
        "Aligner - Wild-Type Reference Alignment",
        4,
        7,
        args.verbose,
    )

    # Find alignment manifest for next step
    current_manifest = find_manifest_file(aligner_output, "alignment_manifest.json")
    print(f"   📄 Manifest Generated: {os.path.basename(current_manifest)}")

    # Step 5: Domino 5 - Analyzer (Mutation Detection & Analysis)
    print("\n🔬 PHASE 5: MUTATION ANALYSIS")
    print("-" * 28)
    analyzer_output = os.path.join(run_directory, "05_analyzer_results")
    analyzer_args = [
        "--manifest",
        current_manifest,
        "--output-dir",
        analyzer_output,
        "--threads",
        str(args.threads),
    ]

    execute_domino_tool(
        str(tools_dir.parent / "analyzer" / "tools" / "run_analyzer.py"),
        analyzer_args,
        "Analyzer - Mutation Detection & Analysis",
        5,
        7,
        args.verbose,
    )

    # Find analysis manifest for next step
    current_manifest = find_manifest_file(analyzer_output, "analysis_manifest.json")
    print(f"   📄 Manifest Generated: {os.path.basename(current_manifest)}")

    # Step 6: Domino 6 - Co-occurrence Analyzer (Mutation Pattern Analysis)
    print("\n🕸️  PHASE 6: PATTERN ANALYSIS")
    print("-" * 28)
    cooccurrence_output = os.path.join(run_directory, "06_cooccurrence_results")
    cooccurrence_args = [
        "--manifest",
        current_manifest,
        "--output-dir",
        cooccurrence_output,
        "--threads",
        str(args.threads),
    ]

    execute_domino_tool(
        str(tools_dir / "run_cooccurrence_analyzer.py"),
        cooccurrence_args,
        "Co-occurrence - Mutation Pattern Analysis",
        6,
        7,
        args.verbose,
    )

    # Find cooccurrence manifest for final step
    current_manifest = find_manifest_file(
        cooccurrence_output, "cooccurrence_manifest.json"
    )
    print(f"   📄 Manifest Generated: {os.path.basename(current_manifest)}")

    # Step 7: Domino 7 - Reporter (Interactive HTML Dashboard)
    print("\n📊 PHASE 7: REPORT GENERATION")
    print("-" * 30)
    reporter_output = os.path.join(run_directory, "07_reporter_results")
    reporter_args = ["--manifest", current_manifest, "--output-dir", reporter_output]

    execute_domino_tool(
        str(tools_dir / "run_reporter.py"),
        reporter_args,
        "Reporter - Interactive HTML Dashboard",
        7,
        7,
        args.verbose,
    )

    # Calculate total pipeline time
    total_pipeline_time = time.time() - pipeline_start_time
    hours = int(total_pipeline_time // 3600)
    minutes = int((total_pipeline_time % 3600) // 60)
    seconds = int(total_pipeline_time % 60)

    print(f"\n⏱️  Pipeline Execution Time: {hours:02d}:{minutes:02d}:{seconds:02d}")

    # Find final HTML report
    final_report = os.path.join(reporter_output, "mutation_analysis_report.html")
    if not os.path.isfile(final_report):
        # Try alternative report name
        final_report = os.path.join(reporter_output, "final_report.html")

    if not os.path.isfile(final_report):
        print("⚠️  Warning: Final HTML report not found in expected location")
        final_report = reporter_output  # Return directory if specific file not found

    return final_report


if __name__ == "__main__":
    # Parse command-line arguments
    parser = create_argument_parser()
    args = parser.parse_args()

    print("🧬 MutationScan Pipeline Orchestrator")
    print("=" * 50)
    print(f"📅 Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📊 Configuration:")
    print(f"   • Accessions file: {args.accessions}")
    print(f"   • Gene list: {args.gene_list}")
    print(f"   • Email: {args.email}")
    print(f"   • Output directory: {args.output_dir}")
    print(f"   • Threads: {args.threads}")

    if args.sepi_species:
        print(f"   • Reference mode: SEPI ({args.sepi_species})")
    else:
        print(f"   • Reference mode: User-provided ({args.user_reference_dir})")

    print()

    try:
        # Validate input files
        print("🔍 Validating input files...")
        validate_input_files(args)

        # Create run directory
        print("📁 Setting up output directory...")
        run_directory = create_run_directory(args.output_dir)

        # Execute the complete pipeline sequence
        print("🚀 Beginning pipeline execution...")
        start_time = time.time()

        final_report_path = run_pipeline_sequence(args, run_directory)

        # Calculate execution time
        end_time = time.time()
        execution_time = end_time - start_time
        hours = int(execution_time // 3600)
        minutes = int((execution_time % 3600) // 60)
        seconds = int(execution_time % 60)

        # GRAND FINALE - Professional Success Message
        print("\n" + "🎉" * 25)
        print("� MUTATIONSCAN PIPELINE COMPLETED SUCCESSFULLY! �")
        print("�🎉" * 25)
        print()
        print("📊 EXECUTION SUMMARY:")
        print("=" * 50)
        print(f"⏱️  Total Runtime:     {hours:02d}:{minutes:02d}:{seconds:02d}")
        print(f"📁 Results Directory:  {run_directory}")
        print(f"� Configuration:")
        print(f"   • Input Genomes:    {os.path.basename(args.accessions)}")
        print(f"   • Target Genes:     {os.path.basename(args.gene_list)}")
        print(f"   • CPU Threads:      {args.threads}")
        if args.sepi_species:
            print(f"   • Reference Mode:   SEPI ({args.sepi_species})")
        else:
            print(f"   • Reference Mode:   User-provided")
        print()

        print("🔬 ANALYSIS PIPELINE RESULTS:")
        print("=" * 50)
        print("✅ Genome Extraction     - NCBI genomes successfully retrieved")
        print("✅ AMR Gene Annotation   - Resistance genes identified with CARD")
        print("✅ Protein Extraction    - Target protein sequences extracted")
        print("✅ Reference Alignment   - Wild-type alignments completed")
        print("✅ Mutation Analysis     - Variants and mutations detected")
        print("✅ Pattern Analysis      - Co-occurrence relationships mapped")
        print("✅ Report Generation     - Interactive dashboard created")
        print()

        print("🎯 YOUR RESULTS ARE READY!")
        print("=" * 50)

        if os.path.isfile(final_report_path):
            print(f"📊 Interactive Report:  {final_report_path}")
            print("🌐 Open this HTML file in your web browser to explore:")
            print("   • Mutation frequency distributions")
            print("   • Gene-by-gene resistance profiles")
            print("   • Co-occurrence network visualizations")
            print("   • Statistical summaries and trends")
        else:
            print(f"📁 Results Directory:   {final_report_path}")
            print("📊 Check the reporter results folder for output files")

        print()
        print("💡 NEXT STEPS:")
        print("   1. Open the HTML report in your favorite browser")
        print("   2. Explore the interactive visualizations")
        print("   3. Export specific datasets for further analysis")
        print("   4. Share findings with your research team")
        print()
        print("🙏 Thank you for using MutationScan!")
        print("   For support: https://github.com/vihaankulkarni29/MutationScan")
        print()

        # Optionally open report in browser
        if args.open_report and os.path.isfile(final_report_path):
            try:
                import webbrowser

                webbrowser.open(f"file://{os.path.abspath(final_report_path)}")
                print("🚀 Opening report in your default browser...")
                print("   (If it doesn't open automatically, copy the file path above)")
            except Exception as e:
                print(f"⚠️  Could not auto-open browser: {e}")
                print("   Please manually open the HTML file in your browser.")

        print("🧬" * 25)

    except subprocess.CalledProcessError as e:
        print(f"\n❌ Pipeline Error: Domino tool failed", file=sys.stderr)
        print(f"   Command: {' '.join(e.cmd)}", file=sys.stderr)
        print(f"   Exit code: {e.returncode}", file=sys.stderr)
        if e.stdout:
            print(f"   Output: {e.stdout}", file=sys.stderr)
        if e.stderr:
            print(f"   Error: {e.stderr}", file=sys.stderr)
        print(
            f"\n💡 Check the individual domino tool logs for detailed error information.",
            file=sys.stderr,
        )
        sys.exit(3)

    except FileNotFoundError as e:
        print(f"\n❌ File System Error: {e}", file=sys.stderr)
        print(
            "💡 Ensure all domino tools are present and manifest files are generated correctly.",
            file=sys.stderr,
        )
        sys.exit(4)

    except (ValueError, OSError) as e:
        print(f"\n❌ Configuration Error: {e}", file=sys.stderr)
        sys.exit(1)

    except KeyboardInterrupt:
        print(f"\n⏹️  Pipeline interrupted by user", file=sys.stderr)
        print(
            "🧹 Clean up any partial results in the output directory if needed.",
            file=sys.stderr,
        )
        sys.exit(130)  # Standard exit code for Ctrl+C

    except Exception as e:
        print(f"\n❌ Unexpected Error: {e}", file=sys.stderr)
        print(f"📍 Error type: {type(e).__name__}", file=sys.stderr)
        import traceback

        print("📋 Full traceback:", file=sys.stderr)
        traceback.print_exc()
        sys.exit(2)
