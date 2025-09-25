#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SubScan Pipeline - Domino 2: The Annotator (High-Performance Parallel Version)

This script serves as an integration wrapper for the ABRicate-Automator tool,
acting as the second domino in the SubScan bioinformatics pipeline with
parallel execution for maximum performance.

Purpose:
- Reads genome_manifest.json from Domino 1 (The Harvester)
- Executes ABRicate-Automator for AMR gene annotation in parallel
- Generates annotation_manifest.json for Domino 3 (FastaAAExtractor)

Author: SubScan Pipeline Development Team
"""

import argparse
import json
import os
import subprocess
import sys

# Set UTF-8 encoding for Windows console
if sys.platform.startswith('win'):
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')
import sys
import multiprocessing
from datetime import datetime
from multiprocessing import Pool
from typing import Dict, List, Tuple, Any

from tqdm import tqdm

# Add src directory to Python path for imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src'))
from subscan.utils import load_json_manifest, save_json_manifest, extract_genomes_from_manifest  # pylint: disable=import-error


def create_argument_parser():
    """Create and configure the command-line argument parser"""
    parser = argparse.ArgumentParser(
        description="SubScan Domino 2: High-Performance AMR Gene Annotation Wrapper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_annotator.py --manifest genome_manifest.json --output-dir ./annotations
  python run_annotator.py --manifest results/genome_manifest.json --output-dir ./amr_results --threads 8

Pipeline Flow:
  Domino 1 (Harvester) → genome_manifest.json → Domino 2 (Annotator) → annotation_manifest.json → Domino 3 (FastaAAExtractor)

Performance Features:
  - Parallel execution using multiprocessing.Pool
  - Individual genome processing for maximum efficiency
  - Real-time progress tracking with tqdm
  - Configurable thread count for optimal performance
        """,
    )

    # Required arguments
    parser.add_argument(
        "--manifest",
        type=str,
        required=True,
        help="Path to the JSON manifest file from the Harvester",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Path to the directory where Abricate results will be saved",
    )

    # Optional arguments
    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help=f"Number of genomes to process in parallel (default: {multiprocessing.cpu_count()} CPU cores)",
    )

    return parser


def validate_arguments(args):
    """Validate command-line arguments"""
    # Check if manifest file exists
    if not os.path.isfile(args.manifest):
        print(f"Error: Manifest file not found: {args.manifest}", file=sys.stderr)
        sys.exit(1)

    # Check if manifest file is JSON
    if not args.manifest.lower().endswith(".json"):
        print(
            f"Warning: Manifest file should be a JSON file: {args.manifest}",
            file=sys.stderr,
        )

    # Validate thread count
    if args.threads < 1:
        print(
            f"Error: Thread count must be at least 1, got: {args.threads}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"[OK] Manifest file: {args.manifest}")
    print(f"[OK] Output directory: {args.output_dir}")
    print(f"[OK] Parallel threads: {args.threads}")


# REFACTORED: Replaced with shared utility function from subscan.utils
# def load_genome_manifest(manifest_path): - REMOVED
# This function is now handled by load_json_manifest() and extract_genomes_from_manifest()
# from the shared utilities module, eliminating code duplication.


def annotate_single_genome(work_item: Tuple[Dict[str, Any], str, List[str], bool]) -> Dict[str, Any]:
    """
    Worker function to annotate a single genome using ABRicate-Automator.
    
    This function is designed to be used with multiprocessing.Pool to enable
    parallel annotation of multiple genomes simultaneously.
    
    Args:
        work_item: Tuple containing:
            - genome_data (Dict[str, Any]): Dictionary with genome information including 'fasta_path'
            - output_dir (str): Directory where annotation results will be saved
            - abricate_command_base (List[str]): Base command for ABRicate-Automator execution
            - use_mock (bool): Whether to use mock tool for testing
    
    Returns:
        Dict[str, Any]: Dictionary containing annotation results:
            - genome_id (str): Genome identifier (accession or filename)
            - success (bool): Whether annotation completed successfully
            - error (str): Error message if annotation failed (None if successful)
            - card_results_path (str): Path to generated TSV results file (None if failed)
    
    Raises:
        No exceptions are raised - all errors are captured and returned in the result dictionary
    """
    genome_data, output_dir, abricate_command_base, use_mock = work_item

    accession = genome_data.get("accession", "")
    if not accession:
        # Extract accession from FASTA filename
        fasta_path = genome_data.get("fasta_path", "")
        if fasta_path:
            filename = os.path.basename(fasta_path)
            accession = filename.replace(".fasta", "").replace(".fa", "")

    genome_id = accession or "unknown"
    fasta_path = genome_data.get("fasta_path", "")

    if not os.path.isfile(fasta_path):
        return {
            "genome_id": genome_id,
            "success": False,
            "error": f"FASTA file not found: {fasta_path}",
            "card_results_path": None,
        }

    try:
        # Define output TSV file path with proper accession naming
        output_tsv_file = f"{accession}.tsv"
        tsv_path = os.path.join(output_dir, output_tsv_file)

        if use_mock:
            # Use mock tool with individual file processing
            command = abricate_command_base + [
                "--input",
                fasta_path,
                "--output",
                tsv_path,
                "--db",
                "card",
                "--mincov",
                "80.0",
                "--minid",
                "90.0",
            ]
        else:
            # Use ABRicate directly
            command = [
                "abricate",
                "--db", "card",
                "--mincov", "80.0",
                "--minid", "90.0",
                fasta_path
            ]

        # Execute ABRicate command
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace',
            timeout=300,  # 5 minute timeout per genome
        )

        # Save the output to TSV file
        with open(tsv_path, 'w', encoding='utf-8') as f:
            f.write(result.stdout)

        # Check if output was created and has content
        if os.path.isfile(tsv_path) and os.path.getsize(tsv_path) > 0:
            return {
                "genome_id": genome_id,
                "success": True,
                "error": None,
                "card_results_path": os.path.abspath(tsv_path),
            }

        return {
                "genome_id": genome_id,
                "success": False,
                "error": f"Expected TSV output not created or empty: {output_tsv_file}",
                "card_results_path": None,
            }

    except subprocess.TimeoutExpired:
        return {
            "genome_id": genome_id,
            "success": False,
            "error": "ABRicate processing timed out (5 minutes)",
            "card_results_path": None,
        }
    except subprocess.CalledProcessError as e:
        return {
            "genome_id": genome_id,
            "success": False,
            "error": f"ABRicate failed: {e}",
            "card_results_path": None,
        }
    except (OSError, PermissionError) as e:
        return {
            "genome_id": genome_id,
            "success": False,
            "error": f"Process error: {str(e)}",
            "card_results_path": None,
        }


def execute_parallel_annotation(genomes: List[Dict[str, Any]], output_dir: str, threads: int) -> List[Dict[str, Any]]:
    """
    Execute ABRicate-Automator annotation in parallel across multiple genomes.
    
    This function orchestrates high-performance parallel processing to annotate
    multiple genomes simultaneously using ABRicate-Automator tool for AMR gene
    identification against the CARD database.
    
    Args:
        genomes (List[Dict[str, Any]]): List of genome dictionaries, each containing:
            - 'fasta_path': Path to the genome FASTA file
            - 'accession': Genome accession number (optional)
        output_dir (str): Directory where annotation results will be saved
        threads (int): Number of parallel worker processes to use
    
    Returns:
        List[Dict[str, Any]]: List of annotation results from worker functions:
            - Each result contains 'genome_id', 'success', 'error', 'card_results_path'
    
    Raises:
        Exception: If parallel processing setup or execution fails
    """
    print("\n🚀 Step 3: Running parallel ABRicate-Automator annotation...")
    print(f"    Output directory: {output_dir}")
    print(f"    Parallel workers: {threads}")
    print(f"    Genomes to process: {len(genomes)}")

    # Prepare base command for ABRicate-Automator (check for mock tool first)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    mock_tool_path = os.path.join(script_dir, "mock_abricate_automator.py")

    if os.path.isfile(mock_tool_path):
        print("🧪 Using mock ABRicate Automator for testing")
        abricate_command_base = ["python", mock_tool_path]
        use_mock = True
    else:
        # Test if abricate is available
        try:
            subprocess.run(
                ["abricate", "--version"],
                capture_output=True,
                check=True,
                encoding='utf-8',
                errors='replace',
                timeout=10,
            )
            print("✅ Using ABRicate directly")
            use_mock = False
        except (
            subprocess.CalledProcessError,
            FileNotFoundError,
            subprocess.TimeoutExpired,
        ):
            print("❌ ABRicate command not found!")
            print("🔧 Falling back to MOCK mode for demonstration purposes...")
            use_mock = True

    # Create work items for parallel processing
    work_items = [
        (genome, output_dir, abricate_command_base, use_mock) for genome in genomes
    ]

    # Execute parallel annotation
    print("\n" + "=" * 50)
    print("Parallel ABRicate-Automator Execution:")
    print("=" * 50)

    results = []
    try:
        with Pool(processes=threads) as pool:
            # Use tqdm for progress tracking
            results = list(
                tqdm(
                    pool.imap(annotate_single_genome, work_items),
                    total=len(work_items),
                    desc="Annotating Genomes",
                    unit="genome",
                )
            )
    except (OSError, MemoryError, KeyboardInterrupt) as e:
        print(f"❌ Parallel execution failed: {e}")
        return []

    # Analyze results
    successful = [r for r in results if r["success"]]
    failed = [r for r in results if not r["success"]]

    print("=" * 50)
    print("✅ Parallel annotation completed!")
    print(f"📊 Results: {len(successful)} successful, {len(failed)} failed")

    if failed:
        print("\n⚠️  Failed annotations:")
        for failure in failed[:5]:  # Show first 5 failures
            print(f"    - {failure['genome_id']}: {failure['error']}")
        if len(failed) > 5:
            print(f"    ... and {len(failed) - 5} more failures")

    return results


def _calculate_annotation_stats(parallel_results):
    """Calculate annotation statistics from parallel results."""
    total_genomes = len(parallel_results)
    successful_annotations = sum(1 for r in parallel_results if r["success"])

    return {
        "total_genomes": total_genomes,
        "successful_annotations": successful_annotations,
        "failed_annotations": total_genomes - successful_annotations,
        "success_rate": (
            f"{(successful_annotations/total_genomes*100):.1f}%"
            if total_genomes > 0
            else "0%"
        ),
        "parallel_execution": True,
        "processing_method": "multiprocessing.Pool with individual genome workers",
    }


def _update_genomes_with_results(genomes_list, parallel_results):
    """Update genome entries with annotation results."""
    # Create mapping from parallel results for efficient lookup
    results_mapping = {r["genome_id"]: r for r in parallel_results}
    updated_count = 0

    for genome in genomes_list:
        accession = genome.get("accession", "")

        # If no accession in manifest, try to extract from FASTA filename
        if not accession:
            fasta_path = genome.get("fasta_path", "")
            if fasta_path:
                filename = os.path.basename(fasta_path)
                accession = filename.replace(".fasta", "").replace(".fa", "")
                # Add the extracted accession to the genome entry
                genome["accession"] = accession

        # Add AMR CARD results path from parallel results
        if accession and accession in results_mapping:
            result = results_mapping[accession]
            if result["success"] and result["card_results_path"]:
                genome["amr_card_results"] = result["card_results_path"]
                updated_count += 1
                print(f"    ✅ Added AMR results for: {accession}")
            else:
                genome["amr_card_results"] = ""
                genome["annotation_error"] = result.get("error", "Unknown error")
                print(f"    ❌ Failed annotation for: {accession} - {result.get('error', 'Unknown error')}")
        else:
            genome["amr_card_results"] = ""
            print(f"    ⚠️  No parallel result found for: {accession}")

    return updated_count


def generate_annotation_manifest(original_manifest, parallel_results, output_dir):
    """Generate the annotation manifest for Domino 3 using parallel results"""
    print("\n📄 Step 4: Generating annotation manifest from parallel results...")

    # Create a deep copy of the original manifest
    annotation_manifest = json.loads(json.dumps(original_manifest))

    # Update pipeline metadata
    annotation_manifest["pipeline_step"] = "Annotator"
    annotation_manifest["generated_by"] = (
        "SubScan Domino 2: The Annotator (High-Performance Parallel)"
    )
    annotation_manifest["generation_timestamp"] = datetime.now().isoformat()

    # Calculate statistics from parallel results
    annotation_manifest["annotation_stats"] = _calculate_annotation_stats(parallel_results)

    # Extract and update genomes with AMR results
    genomes_list = extract_genomes_from_manifest(annotation_manifest)
    updated_count = _update_genomes_with_results(genomes_list, parallel_results)

    print(f"✅ Updated {updated_count} genome entries with AMR results")

    # Save the annotation manifest
    manifest_path = os.path.join(output_dir, "annotation_manifest.json")

    try:
        save_json_manifest(annotation_manifest, manifest_path)

        print(f"✅ Annotation manifest saved: {manifest_path}")
        print("📊 High-Performance Statistics:")
        print(f"    - Total genomes: {annotation_manifest['annotation_stats']['total_genomes']}")
        print(f"    - Successful annotations: {annotation_manifest['annotation_stats']['successful_annotations']}")
        print(f"    - Success rate: {annotation_manifest['annotation_stats']['success_rate']}")
        print("    - Processing method: Parallel (individual genome workers)")

        return manifest_path

    except (OSError, PermissionError, json.JSONDecodeError) as e:
        print(f"❌ Failed to save annotation manifest: {e}")
        return None


def main() -> int:
    """
    Main entry point for the high-performance AMR gene annotation tool.
    
    This function orchestrates the complete annotation workflow:
    1. Loads and validates the input genome manifest from Domino 1 (Harvester)
    2. Validates that all FASTA files exist and are accessible
    3. Executes parallel ABRicate-Automator annotation using multiprocessing
    4. Aggregates results and generates annotation manifest for Domino 3 (Extractor)
    
    The tool uses high-performance parallel processing to annotate multiple genomes
    simultaneously against the CARD database for antimicrobial resistance gene identification.
    
    Returns:
        int: Exit code (0 for success, 1 for error, 130 for user interruption)
    
    Command Line Usage:
        python run_annotator.py --manifest genome_manifest.json --output-dir ./annotations --threads 8
    
    Pipeline Integration:
        Domino 1 (Harvester) → genome_manifest.json → Domino 2 (Annotator) → annotation_manifest.json → Domino 3 (Extractor)
    """
    parser = create_argument_parser()
    args = parser.parse_args()

    # Validate arguments
    validate_arguments(args)

    print("\n" + "=" * 60)
    print("SubScan Pipeline - Domino 2: The Annotator (High-Performance)")
    print("=" * 60)
    print(f"Input manifest: {args.manifest}")
    print(f"Output directory: {args.output_dir}")
    print(f"Parallel workers: {args.threads} threads")
    print("🚀 Performance Mode: ENABLED")
    print("=" * 60)

    # High-performance parallel processing workflow
    try:
        # Step 1: Load and validate the genome manifest
        print("\n📋 Step 1: Loading genome manifest...")
        manifest = load_json_manifest(args.manifest, required_keys=["output_files"])
        genomes = extract_genomes_from_manifest(manifest)

        print(f"[OK] Loaded manifest: {len(genomes)} genomes found")

        if not genomes:
            print("❌ No genomes found in manifest")
            return 1

        # Step 2: Validate FASTA files exist
        print("\n🔗 Step 2: Validating input files...")
        valid_genomes = []
        for genome in genomes:
            fasta_path = genome.get("fasta_path", "")
            if os.path.isfile(fasta_path):
                valid_genomes.append(genome)
                print(f"    ✅ Valid: {os.path.basename(fasta_path)}")
            else:
                print(f"    ❌ Missing: {fasta_path}")

        if not valid_genomes:
            print("❌ No valid FASTA files found to process")
            return 1

        print(
            f"✅ Found {len(valid_genomes)} valid genomes out of {len(genomes)} total"
        )

        # Step 3: Execute parallel ABRicate-Automator annotation
        parallel_results = execute_parallel_annotation(
            valid_genomes, args.output_dir, args.threads
        )

        if not parallel_results:
            print("❌ Parallel annotation execution failed")
            return 1

        # Step 4: Generate annotation manifest from parallel results
        manifest_path = generate_annotation_manifest(
            manifest, parallel_results, args.output_dir
        )

        if not manifest_path:
            print("❌ Failed to generate annotation manifest")
            return 1

        # Final success message
        print("\n" + "=" * 60)
        print("🎉 DOMINO 2 COMPLETE: The Annotator (High-Performance)")
        print("=" * 60)
        print("✅ All phases completed successfully with parallel processing!")
        print(f"📄 Input manifest: {args.manifest}")
        print(f"📁 ABRicate results: {args.output_dir}")
        print(f"📋 Output manifest: {manifest_path}")
        print(f"⚡ Performance: {args.threads} parallel workers")
        print("\n🔄 Pipeline Status:")
        print("    Domino 1 (Harvester) → ✅ Complete")
        print("    Domino 2 (Annotator) → ✅ Complete (High-Performance)")
        print("    Domino 3 (FastaAAExtractor) → 🚀 Ready to proceed")
        print("=" * 60)

        return 0

    except KeyboardInterrupt:
        print("\n⚠️  Process interrupted by user")
        return 1
    except (OSError, PermissionError, ImportError, ValueError) as e:
        print(f"\n❌ Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
