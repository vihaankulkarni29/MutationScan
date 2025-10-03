#!/usr/bin/env python3
"""
MutationScan Domino 3: Protein Sequence Extractor

This module implements the third domino in the MutationScan pipeline, responsible
for extracting protein sequences from annotated genomes based on AMR gene coordinates.
It processes annotation results from Domino 2 and generates FASTA files containing
protein sequences for target genes specified in the user-provided gene list.

The extractor interfaces with the FastaAAExtractor tool for high-performance protein
sequence extraction and creates standardized JSON manifests for downstream alignment
and mutation analysis workflows.

Key Features:
- High-performance parallel processing of multiple genomes
- Precise protein extraction based on genome coordinates
- Gene filtering using user-defined target gene lists
- Comprehensive error handling and progress tracking
- Mock mode support for testing environments
- Standardized manifest generation for pipeline integration

Usage:
    python run_extractor.py --manifest annotation_manifest.json --gene-list targets.txt --output-dir ./proteins

Integration:
    - Input: annotation_manifest.json from Annotator (Domino 2) + gene list file
    - Output: protein_manifest.json + extracted protein FASTA files
    - Next Domino: run_aligner.py (sequence alignment to references)

Dependencies:
    - FastaAAExtractor tool for coordinate-based protein extraction
    - Valid annotation TSV files with gene coordinates
    - Target gene list file for filtering relevant proteins

Author: MutationScan Development Team
Version: 1.0.0
"""

import argparse
import os
import sys
import multiprocessing
import subprocess
import json
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
from tqdm import tqdm


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure command-line argument parser for protein extraction.
    
    Sets up argument parsing for the extractor with validation for annotation
    manifest input, gene list specification, and output directory configuration
    for high-throughput protein sequence extraction workflows.

    Returns:
        argparse.ArgumentParser: Configured parser for extractor arguments with:
            - manifest: Path to annotation_manifest.json from Domino 2 (Annotator)
            - gene_list: Path to text file containing target genes to extract
            - output_dir: Directory for extracted protein FASTA files and manifest
            - threads: Number of parallel processes for extraction

    Example:
        >>> parser = create_argument_parser()
        >>> args = parser.parse_args(['--manifest', 'annotation_manifest.json', '--gene-list', 'targets.txt'])
    """
    parser = argparse.ArgumentParser(
        description="SubScan Domino 3: High-Performance Protein Extraction Wrapper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_extractor.py --manifest annotation_manifest.json --gene-list targets.txt --output-dir ./proteins
  python run_extractor.py --manifest results/annotation_manifest.json --gene-list resistance_genes.txt --output-dir ./extracted --threads 16

Pipeline Flow:
  Domino 2 (Annotator) → annotation_manifest.json → Domino 3 (Extractor) → protein_manifest.json → Domino 4 (Aligner)

High-Performance Features:
  - Parallel protein extraction using multiprocessing.Pool
  - Individual genome processing for maximum efficiency
  - Real-time progress tracking with tqdm
  - Configurable thread count for optimal performance
  - Precise gene filtering based on user-defined target lists

Input Requirements:
  - annotation_manifest.json: Contains genome FASTA paths and TSV coordinate files
  - gene_list.txt: User-defined list of target genes/proteins to extract
  - Functional fasta_aa_extractor tool available in PATH
        """,
    )

    # Required arguments
    parser.add_argument(
        "--manifest",
        type=str,
        required=True,
        help="Path to the annotation_manifest.json file from the Annotator",
    )

    parser.add_argument(
        "--gene-list",
        type=str,
        required=True,
        help="Path to the user's target gene list file",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Path to save the extracted protein FASTA files",
    )

    # Optional arguments
    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help=f"Number of parallel processes (default: {multiprocessing.cpu_count()} CPU cores)",
    )

    return parser


def validate_arguments(args):
    """Validate command-line arguments and inputs"""
    # Check if annotation manifest file exists
    if not os.path.isfile(args.manifest):
        print(
            f"Error: Annotation manifest file not found: {args.manifest}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check if manifest file is JSON
    if not args.manifest.lower().endswith(".json"):
        print(
            f"Warning: Manifest file should be a JSON file: {args.manifest}",
            file=sys.stderr,
        )

    # Check if gene list file exists
    if not os.path.isfile(args.gene_list):
        print(f"Error: Gene list file not found: {args.gene_list}", file=sys.stderr)
        sys.exit(1)

    # Validate thread count
    if args.threads < 1:
        print(
            f"Error: Thread count must be at least 1, got: {args.threads}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"✓ Annotation manifest: {args.manifest}")
    print(f"✓ Gene list file: {args.gene_list}")
    print(f"✓ Output directory: {args.output_dir}")
    print(f"✓ Parallel processes: {args.threads}")


def extract_proteins_for_single_genome(work_item):
    """
    Worker function to extract proteins for a single genome using FastaAAExtractor

    This function is designed to be the core 'work unit' for parallel processing.
    It handles all aspects of protein extraction for one genome, from input validation
    to command execution and result reporting.

    Args:
        work_item (tuple): Contains (genome_info_dict, gene_list_path, output_dir)
            - genome_info_dict: Dictionary containing genome information including:
                - 'accession': Genome accession number
                - 'fasta_path': Path to genome FASTA file
                - 'amr_card_results': Path to TSV coordinate file
            - gene_list_path: Path to user's target gene list file
            - output_dir: Directory where extracted proteins should be saved

    Returns:
        tuple or None: On success, returns (accession, protein_fasta_path)
                      On failure, returns None
    """
    try:
        # Unpack the work item
        genome_info_dict, gene_list_path, output_dir = work_item

        # Extract required information from genome dictionary
        accession = genome_info_dict.get("accession", "")
        genome_fasta_path = genome_info_dict.get("fasta_path", "")
        tsv_coordinate_path = genome_info_dict.get("amr_card_results", "")

        # Validate required inputs
        if not accession:
            print(f"⚠️  Warning: No accession found in genome info")
            return None

        if not os.path.isfile(genome_fasta_path):
            print(
                f"⚠️  Warning: Genome FASTA not found for {accession}: {genome_fasta_path}"
            )
            return None

        if not os.path.isfile(tsv_coordinate_path):
            print(
                f"⚠️  Warning: TSV coordinate file not found for {accession}: {tsv_coordinate_path}"
            )
            return None

        if not os.path.isfile(gene_list_path):
            print(f"⚠️  Warning: Gene list file not found: {gene_list_path}")
            return None

        # Construct output filename for extracted proteins
        output_filename = f"{accession}_proteins.faa"
        output_protein_path = os.path.join(output_dir, output_filename)

        # Construct the full command-line call to FastaAAExtractor
        # Expected command format: fasta_aa_extractor --genome FASTA --coordinates TSV --genes GENELIST --output OUTPUT

        # For testing, check if mock tool exists in the same directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        mock_tool_path = os.path.join(script_dir, "mock_fasta_aa_extractor.py")

        if os.path.isfile(mock_tool_path):
            # Use mock tool for testing
            command = [
                "python",
                mock_tool_path,
                "--genome",
                genome_fasta_path,
                "--coordinates",
                tsv_coordinate_path,
                "--genes",
                gene_list_path,
                "--output",
                output_protein_path,
            ]
        else:
            # Use real tool
            command = [
                "fasta_aa_extractor",
                "--genome",
                genome_fasta_path,
                "--coordinates",
                tsv_coordinate_path,
                "--genes",
                gene_list_path,
                "--output",
                output_protein_path,
            ]

        # Execute the FastaAAExtractor command
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout per genome
        )

        # Verify that the output file was created
        if os.path.isfile(output_protein_path):
            return (accession, os.path.abspath(output_protein_path))
        else:
            print(
                f"⚠️  Warning: Expected output file not created for {accession}: {output_filename}"
            )
            return None

    except subprocess.TimeoutExpired:
        print(f"❌ Timeout: FastaAAExtractor timed out for {accession} (5 minutes)")
        return None

    except subprocess.CalledProcessError as e:
        print(f"❌ Error: FastaAAExtractor failed for {accession}: {e.stderr}")
        return None

    except FileNotFoundError:
        print(f"❌ Error: fasta_aa_extractor command not found in PATH")
        return None

    except Exception as e:
        print(f"❌ Unexpected error processing {accession}: {str(e)}")
        return None


def execute_parallel_extraction(
    annotation_manifest_path, gene_list_path, output_dir, num_processes
):
    """
    Orchestrates parallel protein extraction across multiple genomes using multiprocessing

    This function manages the entire parallel processing pipeline, from loading the
    annotation manifest to distributing work across worker processes and collecting results.

    Args:
        annotation_manifest_path (str): Path to annotation_manifest.json from Domino 2
        gene_list_path (str): Path to user's target gene list file
        output_dir (str): Directory where extracted proteins should be saved
        num_processes (int): Number of parallel worker processes to use

    Returns:
        list: List of successful extractions as (accession, protein_fasta_path) tuples
    """
    try:
        # Load the annotation manifest produced by Domino 2
        with open(annotation_manifest_path, "r") as f:
            manifest_data = json.load(f)

        # Extract genome list from manifest
        genome_list = manifest_data.get("genomes", [])
        if not genome_list and "output_files" in manifest_data:
            genome_list = manifest_data["output_files"].get("genomes", [])

        if not genome_list:
            print("❌ Error: No genomes found in annotation manifest")
            return []

        print(f"📊 Processing {len(genome_list)} genomes for protein extraction...")

        # Prepare work items for parallel processing
        # Each work item is a tuple: (genome_info_dict, gene_list_path, output_dir)
        work_items = []
        for genome_info in genome_list:
            work_item = (genome_info, gene_list_path, output_dir)
            work_items.append(work_item)

        # Execute parallel processing with progress tracking
        successful_extractions = []

        with multiprocessing.Pool(processes=num_processes) as pool:
            # Initialize progress tracking
            with tqdm(
                total=len(work_items), desc="🧬 Extracting proteins", unit="genome"
            ) as pbar:

                # Submit all work items for parallel processing
                results = []
                for work_item in work_items:
                    result = pool.apply_async(
                        extract_proteins_for_single_genome, (work_item,)
                    )
                    results.append(result)

                # Collect results as they complete
                for result in results:
                    try:
                        extraction_result = result.get(
                            timeout=360
                        )  # 6 minute timeout per genome
                        if extraction_result is not None:
                            successful_extractions.append(extraction_result)
                        pbar.update(1)
                    except multiprocessing.TimeoutError:
                        print("⏱️  Warning: Worker process timed out")
                        pbar.update(1)
                    except Exception as e:
                        print(f"⚠️  Warning: Worker process error: {str(e)}")
                        pbar.update(1)

        # Report processing summary
        success_count = len(successful_extractions)
        total_count = len(genome_list)
        success_rate = (success_count / total_count * 100) if total_count > 0 else 0

        print(f"\n🎯 Protein extraction summary:")
        print(
            f"   ✅ Successful: {success_count}/{total_count} genomes ({success_rate:.1f}%)"
        )

        if success_count < total_count:
            failed_count = total_count - success_count
            print(f"   ❌ Failed: {failed_count} genomes")

        return successful_extractions

    except FileNotFoundError:
        print(f"❌ Error: Annotation manifest not found: {annotation_manifest_path}")
        return []

    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON in annotation manifest: {str(e)}")
        return []

    except Exception as e:
        print(f"❌ Unexpected error in parallel extraction: {str(e)}")
        return []


def generate_protein_manifest(
    successful_extractions, output_dir, gene_list_path, annotation_manifest_path
):
    """
    Generates protein_manifest.json for handoff to Domino 4

    This function creates a standardized JSON manifest that documents all successfully
    extracted protein files, their locations, and metadata for downstream processing.

    Args:
        successful_extractions (list): List of (accession, protein_fasta_path) tuples
        output_dir (str): Directory where protein files were saved
        gene_list_path (str): Path to the user's gene list file
        annotation_manifest_path (str): Path to input annotation manifest

    Returns:
        str: Path to the generated protein_manifest.json file
    """
    try:
        # Prepare manifest data structure
        manifest_data = {
            "metadata": {
                "tool": "SubScan-Extractor",
                "domino": 3,
                "description": "Protein extraction results from FastaAAExtractor",
                "input_manifest": os.path.abspath(annotation_manifest_path),
                "gene_list": os.path.abspath(gene_list_path),
                "output_directory": os.path.abspath(output_dir),
                "total_genomes_processed": len(successful_extractions),
                "timestamp": None,  # Will be set by subprocess if available
            },
            "extracted_proteins": [],
            "summary": {
                "total_protein_files": len(successful_extractions),
                "extraction_success_rate": "100%" if successful_extractions else "0%",
            },
        }

        # Add individual protein file entries
        for accession, protein_fasta_path in successful_extractions:
            protein_entry = {
                "accession": accession,
                "protein_fasta_path": os.path.abspath(protein_fasta_path),
                "file_size_bytes": None,
                "protein_count": None,
            }

            # Try to get file size if file exists
            if os.path.isfile(protein_fasta_path):
                try:
                    protein_entry["file_size_bytes"] = os.path.getsize(
                        protein_fasta_path
                    )
                except OSError:
                    protein_entry["file_size_bytes"] = "unknown"

                # Try to count proteins in FASTA file
                try:
                    with open(protein_fasta_path, "r") as f:
                        protein_count = sum(1 for line in f if line.startswith(">"))
                    protein_entry["protein_count"] = protein_count
                except (IOError, UnicodeDecodeError):
                    protein_entry["protein_count"] = "unknown"

            manifest_data["extracted_proteins"].append(protein_entry)

        # Try to get current timestamp
        try:
            import datetime

            manifest_data["metadata"]["timestamp"] = datetime.datetime.now().isoformat()
        except ImportError:
            manifest_data["metadata"]["timestamp"] = "unavailable"

        # Write manifest to output directory
        manifest_path = os.path.join(output_dir, "protein_manifest.json")

        with open(manifest_path, "w") as f:
            json.dump(manifest_data, f, indent=2)

        print(f"📄 Generated protein manifest: {os.path.basename(manifest_path)}")
        print(f"   📁 Location: {os.path.dirname(manifest_path)}")
        print(f"   🔢 Protein files documented: {len(successful_extractions)}")

        return os.path.abspath(manifest_path)

    except Exception as e:
        print(f"❌ Error generating protein manifest: {str(e)}")
        return None


def main():
    """Main entry point for the high-performance Extractor wrapper"""
    parser = create_argument_parser()
    args = parser.parse_args()

    # Validate the Input Manifest: Check if manifest file exists before processing
    if not os.path.isfile(args.manifest):
        print(f"ERROR: Input manifest not found at '{args.manifest}'. This is likely because the previous pipeline step (Annotator) failed.", file=sys.stderr)
        sys.exit(1)

    # Validate arguments
    validate_arguments(args)

    print("\n" + "=" * 60)
    print("SubScan Pipeline - Domino 3: The Extractor")
    print("=" * 60)
    print(f"✓ Annotation manifest: {args.manifest}")
    print(f"✓ Gene list: {args.gene_list}")
    print(f"✓ Output directory: {args.output_dir}")
    print(f"✓ Parallel processes: {args.threads}")
    print("🚀 High-Performance Mode: ENABLED")
    print("=" * 60)

    try:
        # Execute the parallel protein extraction
        successful_extractions = execute_parallel_extraction(
            args.manifest, args.gene_list, args.output_dir, args.threads
        )

        # Generate the protein manifest for Domino 4
        manifest_path = generate_protein_manifest(
            successful_extractions, args.output_dir, args.gene_list, args.manifest
        )

        # Final summary
        if successful_extractions:
            print(f"\n🎉 Protein extraction completed successfully!")
            print(f"   ✅ {len(successful_extractions)} genomes processed")
            print(
                f"   � Manifest ready for Domino 4: {os.path.basename(manifest_path) if manifest_path else 'Failed'}"
            )
            return 0
        else:
            print(f"\n⚠️  No proteins were successfully extracted")
            print(f"   🔍 Check gene list and annotation results")
            return 1

    except KeyboardInterrupt:
        print(f"\n🛑 Process interrupted by user")
        return 1
    except Exception as e:
        print(f"\n❌ Fatal error: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
