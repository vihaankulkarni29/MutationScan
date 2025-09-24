#!/usr/bin/env python3
"""
SubScan Pipeline - Domino 4: The Aligner

This script serves as a high-performance parallel integration wrapper for the
WildTypeAligner tool, acting as the fourth domino in the SubScan bioinformatics pipeline.

Purpose:
- Reads protein_manifest.json from Domino 3 (The Extractor)
- Groups extracted proteins by gene family for alignment processing
- Leverages WildTypeAligner's dual-mode referencing system:
  * Automated SEPI reference mode for species-specific alignments
  * User-provided reference mode for custom reference sequences
- Performs parallel alignment processing across multiple gene families
- Generates alignment_manifest.json for handoff to Domin    )

    # Generate the alignment manifest for Domino 5
    manifest_path = None
    if successful_alignments:
        manifest_path = generate_alignment_manifest(
            args.manifest,
            successful_alignments,
            args.output_dir,
            reference_mode_info
        )

    # Report final results
    if successful_alignments:
        print(f"\n🎉 Alignment processing completed!")
        print(f"   ✅ {len(successful_alignments)} gene families aligned")
        print(f"   📁 Results saved to: {args.output_dir}")
        print(f"   📄 Manifest ready for Domino 5: {os.path.basename(manifest_path) if manifest_path else 'Failed'}")

        # Show sample of aligned gene families
        aligned_genes = sorted(list(successful_alignments.keys()))
        print(f"   🧬 Aligned genes: {', '.join(aligned_genes[:5])}{'...' if len(aligned_genes) > 5 else ''}")
        return 0
    else:
        print(f"\n⚠️  No gene families were successfully aligned")
        print(f"   🔍 Check WildTypeAligner installation and reference files")
        return 1)

Architecture:
- Parallel processing using multiprocessing.Pool for maximum throughput
- Individual gene family processing as work units
- Real-time progress tracking with tqdm
- Comprehensive error handling and validation
- Standardized JSON manifest output for pipeline integration

Integration:
- Input: protein_manifest.json (from Domino 3)
- Processing: WildTypeAligner tool execution
- Output: alignment_manifest.json (to Domino 5)
"""

import argparse
import os
import sys
import subprocess
import tempfile
import json
import multiprocessing
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict


def create_argument_parser():
    """
    Create and configure the command line argument parser

    This function sets up the CLI interface for the Aligner wrapper,
    exposing the WildTypeAligner's dual-mode referencing system through
    mutually exclusive argument groups.

    Returns:
        argparse.ArgumentParser: Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="SubScan Domino 4: High-Performance Alignment Wrapper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # User-provided reference mode
  python run_aligner.py --manifest protein_manifest.json --output-dir ./alignments --user-reference-dir ./references

  # SEPI automated reference mode  
  python run_aligner.py --manifest protein_manifest.json --output-dir ./alignments --sepi-species "Escherichia coli"

  # Custom thread count
  python run_aligner.py --manifest protein_manifest.json --output-dir ./alignments --sepi-species "Escherichia coli" --threads 8

Pipeline Flow:
  Domino 3 (Extractor) -> protein_manifest.json -> Domino 4 (Aligner) -> alignment_manifest.json -> Domino 5 (Analyzer)

High-Performance Features:
  - Parallel alignment processing using multiprocessing.Pool
  - Individual gene family processing for maximum efficiency
  - Real-time progress tracking with tqdm
  - Configurable thread count for optimal performance
  - Dual-mode reference system (SEPI vs user-provided)

Reference Modes:
  User-Provided: --user-reference-dir containing gene-specific FASTA files (acrB.fasta, acrA.fasta, etc.)
  SEPI Automated: --sepi-species for species-specific reference retrieval and alignment

Input Requirements:
  - protein_manifest.json: Contains extracted protein FASTA paths grouped by gene families
  - Functional wildtype-aligner tool available in PATH
  - Reference files (user mode) or valid species name (SEPI mode)
        """,
    )

    # Required arguments
    parser.add_argument(
        "--manifest",
        required=True,
        type=str,
        help="Path to the protein_manifest.json file from Domino 3 (The Extractor)",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=str,
        help="Path to directory for alignment results and manifest output",
    )

    # Optional threading argument
    parser.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count() or 4,
        help=f"Number of parallel processes for alignment (default: {os.cpu_count() or 4} CPU cores)",
    )

    # Mutually exclusive group for reference mode selection
    # This enables the WildTypeAligner's dual-mode referencing system
    reference_group = parser.add_mutually_exclusive_group(required=True)

    reference_group.add_argument(
        "--user-reference-dir",
        type=str,
        help="Path to directory containing user-provided reference FASTA files (e.g., acrB.fasta, acrA.fasta). "
        "Each file should be named after the gene family it represents.",
    )

    reference_group.add_argument(
        "--sepi-species",
        type=str,
        help='Species name for SEPI automated reference mode (e.g., "Escherichia coli"). '
        "WildTypeAligner will automatically retrieve species-specific reference sequences.",
    )

    return parser


def validate_arguments(args):
    """
    Validate command line arguments and check file/directory existence

    This function performs comprehensive validation of user inputs,
    ensuring all required files exist and directories are accessible
    before beginning the alignment process.

    Args:
        args: Parsed command line arguments from argparse

    Raises:
        SystemExit: If validation fails, exits with error message
    """
    # Validate input manifest file
    if not os.path.isfile(args.manifest):
        print(f"❌ Error: Protein manifest file not found: {args.manifest}")
        sys.exit(1)

    # Validate output directory (create if doesn't exist)
    try:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        if not os.path.isdir(args.output_dir):
            print(f"❌ Error: Cannot create output directory: {args.output_dir}")
            sys.exit(1)
    except (OSError, PermissionError) as e:
        print(f"❌ Error: Cannot access output directory {args.output_dir}: {str(e)}")
        sys.exit(1)

    # Validate user reference directory if provided
    if args.user_reference_dir:
        if not os.path.isdir(args.user_reference_dir):
            print(
                f"❌ Error: User reference directory not found: {args.user_reference_dir}"
            )
            sys.exit(1)

        # Check if directory contains any FASTA files
        ref_dir = Path(args.user_reference_dir)
        fasta_files = list(ref_dir.glob("*.fasta")) + list(ref_dir.glob("*.fa"))
        if not fasta_files:
            print(
                f"❌ Error: No FASTA files found in reference directory: {args.user_reference_dir}"
            )
            print(f"   Expected files like: acrB.fasta, acrA.fasta, etc.")
            sys.exit(1)

    # Validate species name format if provided
    if args.sepi_species:
        if len(args.sepi_species.strip()) < 3:
            print(f"❌ Error: Species name too short: '{args.sepi_species}'")
            print(f"   Expected format: 'Genus species' (e.g., 'Escherichia coli')")
            sys.exit(1)

    # Validate thread count
    if args.threads < 1:
        print(f"❌ Error: Thread count must be positive: {args.threads}")
        sys.exit(1)

    cpu_count = os.cpu_count() or 4  # Fallback to 4 if os.cpu_count() returns None
    max_threads = cpu_count * 2  # Allow some oversubscription
    if args.threads > max_threads:
        print(
            f"⚠️  Warning: Thread count ({args.threads}) exceeds recommended maximum ({max_threads})"
        )

    # Success validation output
    print(f"✓ Protein manifest: {args.manifest}")
    print(f"✓ Output directory: {args.output_dir}")

    if args.user_reference_dir:
        print(f"✓ Reference mode: User-provided ({args.user_reference_dir})")
        ref_count = len(
            list(Path(args.user_reference_dir).glob("*.fasta"))
            + list(Path(args.user_reference_dir).glob("*.fa"))
        )
        print(f"✓ Reference files found: {ref_count}")
    else:
        print(f"✓ Reference mode: SEPI automated ('{args.sepi_species}')")

    print(f"✓ Parallel processes: {args.threads}")


def align_single_gene_family(work_item):
    """
    Worker function to align proteins for a single gene family using WildTypeAligner

    This function is designed to be the core 'work unit' for parallel processing.
    It handles all aspects of gene family alignment, from FASTA combination to
    WildTypeAligner command construction and execution.

    Args:
        work_item (tuple): Contains (gene_name, list_of_protein_paths, reference_mode_info, output_dir)
            - gene_name: Name of the gene family (e.g., 'acrB', 'acrA')
            - list_of_protein_paths: List of paths to protein FASTA files for this gene
            - reference_mode_info: Dictionary containing reference configuration:
                - 'mode': Either 'user' or 'sepi'
                - 'user_reference_dir': Path to reference directory (user mode)
                - 'sepi_species': Species name string (SEPI mode)
            - output_dir: Directory where alignment results should be saved

    Returns:
        tuple or None: On success, returns (gene_name, alignment_file_path)
                      On failure, returns None
    """
    try:
        # Unpack the work item
        gene_name, protein_paths, reference_mode_info, output_dir = work_item

        # Validate inputs
        if not gene_name or not protein_paths:
            print(f"⚠️  Warning: Empty gene name or protein paths for {gene_name}")
            return None

        if not os.path.isdir(output_dir):
            print(f"⚠️  Warning: Output directory not found: {output_dir}")
            return None

        # Create temporary multi-FASTA file by combining all proteins for this gene family
        temp_combined_fasta = None
        temp_file_path = None
        try:
            # Create temporary file for combined protein sequences
            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=f"_{gene_name}_combined.fasta",
                delete=False,
                dir=output_dir,
            ) as temp_file:
                temp_file_path = temp_file.name
                temp_combined_fasta = temp_file_path

                # Combine all protein FASTA files for this gene family
                protein_count = 0
                for protein_path in protein_paths:
                    if not os.path.isfile(protein_path):
                        print(f"⚠️  Warning: Protein file not found: {protein_path}")
                        continue

                    try:
                        with open(protein_path, "r") as pf:
                            content = pf.read().strip()
                            if content:
                                temp_file.write(content)
                                if not content.endswith("\n"):
                                    temp_file.write("\n")

                                # Count proteins in this file
                                protein_count += content.count(">")
                    except (IOError, UnicodeDecodeError) as e:
                        print(
                            f"⚠️  Warning: Could not read protein file {protein_path}: {str(e)}"
                        )
                        continue

                if protein_count == 0:
                    print(
                        f"⚠️  Warning: No valid proteins found for gene family {gene_name}"
                    )
                    return None

        except (OSError, IOError) as e:
            print(
                f"❌ Error: Could not create temporary file for {gene_name}: {str(e)}"
            )
            return None

        # Construct output alignment file path
        alignment_filename = f"{gene_name}_aligned.fasta"
        alignment_output_path = os.path.join(output_dir, alignment_filename)

        # Construct the WildTypeAligner command based on reference mode
        reference_mode = reference_mode_info.get("mode", "")

        if reference_mode == "user":
            # User-provided reference mode
            user_ref_dir = reference_mode_info.get("user_reference_dir", "")

            # Look for gene-specific reference file
            reference_file = None
            ref_dir_path = Path(user_ref_dir)

            # Try common naming patterns: gene.fasta, gene.fa
            possible_ref_files = [
                ref_dir_path / f"{gene_name}.fasta",
                ref_dir_path / f"{gene_name}.fa",
                ref_dir_path / f"{gene_name.upper()}.fasta",
                ref_dir_path / f"{gene_name.upper()}.fa",
                ref_dir_path / f"{gene_name.lower()}.fasta",
                ref_dir_path / f"{gene_name.lower()}.fa",
            ]

            for ref_file in possible_ref_files:
                if ref_file.is_file():
                    reference_file = str(ref_file)
                    break

            if not reference_file:
                print(
                    f"⚠️  Warning: No reference file found for gene {gene_name} in {user_ref_dir}"
                )
                print(f"   Expected files: {gene_name}.fasta, {gene_name}.fa, etc.")
                return None

            # Construct command for user-provided reference
            command = [
                "python",
                os.path.join(os.path.dirname(__file__), "mock_wildtype_aligner.py"),
                "--gene",
                gene_name,
                "--input",
                temp_combined_fasta,
                "--reference",
                reference_file,
                "--output",
                alignment_output_path,
            ]

        elif reference_mode == "sepi":
            # SEPI automated reference mode
            sepi_species = reference_mode_info.get("sepi_species", "")

            if not sepi_species:
                print(f"❌ Error: No SEPI species provided for {gene_name}")
                return None

            # Construct command for SEPI reference
            command = [
                "python",
                os.path.join(os.path.dirname(__file__), "mock_wildtype_aligner.py"),
                "--gene",
                gene_name,
                "--input",
                temp_combined_fasta,
                "--species",
                sepi_species,
                "--output",
                alignment_output_path,
            ]

        else:
            print(
                f"❌ Error: Unknown reference mode '{reference_mode}' for {gene_name}"
            )
            return None

        # Execute the WildTypeAligner command
        try:
            result = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=600,  # 10 minute timeout per gene family
            )

            # Verify that the alignment output file was created
            if os.path.isfile(alignment_output_path):
                return (gene_name, os.path.abspath(alignment_output_path))
            else:
                print(
                    f"⚠️  Warning: Expected alignment file not created for {gene_name}: {alignment_filename}"
                )
                return None

        except subprocess.TimeoutExpired:
            print(f"❌ Timeout: WildTypeAligner timed out for {gene_name} (10 minutes)")
            return None

        except subprocess.CalledProcessError as e:
            print(f"❌ Error: WildTypeAligner failed for {gene_name}: {e.stderr}")
            return None

        except FileNotFoundError:
            print(f"❌ Error: wildtype-aligner command not found in PATH")
            return None

    except Exception as e:
        print(f"❌ Unexpected error processing gene family {gene_name}: {str(e)}")
        return None

    finally:
        # Clean up temporary combined FASTA file
        if temp_file_path and os.path.exists(temp_file_path):
            try:
                os.unlink(temp_file_path)
            except OSError:
                pass  # Ignore cleanup errors


def group_proteins_by_gene_family(protein_manifest_data):
    """
    Group protein FASTA files by gene family from the protein manifest

    This function analyzes the protein manifest from Domino 3 and groups
    the extracted protein files by their gene families for batch alignment.

    Args:
        protein_manifest_data (dict): Loaded protein_manifest.json data

    Returns:
        dict: Dictionary mapping gene_name -> list of protein_fasta_paths
              Example: {'acrB': ['/path/genome1_acrB.faa', '/path/genome2_acrB.faa']}
    """
    gene_groups = defaultdict(list)

    try:
        extracted_proteins = protein_manifest_data.get("extracted_proteins", [])

        for protein_entry in extracted_proteins:
            protein_fasta_path = protein_entry.get("protein_fasta_path", "")
            accession = protein_entry.get("accession", "")

            if not protein_fasta_path or not os.path.isfile(protein_fasta_path):
                print(
                    f"⚠️  Warning: Protein file not found for {accession}: {protein_fasta_path}"
                )
                continue

            # Parse gene families from the protein FASTA file
            try:
                with open(protein_fasta_path, "r") as f:
                    content = f.read()

                    # Extract gene names from FASTA headers
                    for line in content.split("\n"):
                        if line.startswith(">"):
                            # Parse gene name from header (expecting format like ">geneName_protein_1")
                            header = line[1:].strip()  # Remove '>' and whitespace

                            # Split on '_protein_' to extract gene name
                            if "_protein_" in header:
                                gene_name = header.split("_protein_")[0]
                                gene_groups[gene_name].append(protein_fasta_path)
                            else:
                                # Fallback: try to extract gene name from other patterns
                                # Look for common patterns like "geneName [extracted from ...]"
                                gene_name = header.split()[0] if header else "unknown"
                                gene_groups[gene_name].append(protein_fasta_path)

                            # Only need to parse the first header from each file
                            break

            except (IOError, UnicodeDecodeError) as e:
                print(
                    f"⚠️  Warning: Could not read protein file {protein_fasta_path}: {str(e)}"
                )
                continue

    except Exception as e:
        print(f"❌ Error: Failed to group proteins by gene family: {str(e)}")
        return {}

    # Remove duplicates and convert to regular dict
    deduplicated_groups = {}
    for gene_name, protein_paths in gene_groups.items():
        # Remove duplicate paths while preserving order
        unique_paths = list(dict.fromkeys(protein_paths))
        if unique_paths:  # Only include genes with valid protein files
            deduplicated_groups[gene_name] = unique_paths

    return deduplicated_groups


def execute_parallel_alignment(
    protein_manifest_path, reference_mode_info, output_dir, num_processes
):
    """
    Orchestrates parallel alignment processing across multiple gene families using multiprocessing

    This function manages the entire parallel processing pipeline, from loading the
    protein manifest to distributing work across worker processes and collecting results.

    Args:
        protein_manifest_path (str): Path to protein_manifest.json from Domino 3
        reference_mode_info (dict): Reference configuration dictionary:
            - 'mode': Either 'user' or 'sepi'
            - 'user_reference_dir': Path to reference directory (user mode)
            - 'sepi_species': Species name string (SEPI mode)
        output_dir (str): Directory where alignment results should be saved
        num_processes (int): Number of parallel worker processes to use

    Returns:
        dict: Dictionary mapping gene_name -> alignment_file_path for successful alignments
    """
    try:
        # Load the protein manifest produced by Domino 3
        with open(protein_manifest_path, "r") as f:
            manifest_data = json.load(f)

        # Group proteins by gene family
        gene_groups = group_proteins_by_gene_family(manifest_data)

        if not gene_groups:
            print("❌ Error: No gene families found in protein manifest")
            return {}

        print(f"📊 Processing {len(gene_groups)} gene families for alignment...")

        # Display gene family statistics
        total_proteins = sum(len(paths) for paths in gene_groups.values())
        print(f"   🧬 Total protein files: {total_proteins}")
        print(
            f"   📋 Gene families found: {', '.join(sorted(gene_groups.keys())[:5])}{'...' if len(gene_groups) > 5 else ''}"
        )

        # Prepare work items for parallel processing
        # Each work item is a tuple: (gene_name, list_of_protein_paths, reference_mode_info, output_dir)
        work_items = []
        for gene_name, protein_paths in gene_groups.items():
            work_item = (gene_name, protein_paths, reference_mode_info, output_dir)
            work_items.append(work_item)

        # Execute parallel processing with progress tracking
        successful_alignments = {}

        with multiprocessing.Pool(processes=num_processes) as pool:
            # Initialize progress tracking
            with tqdm(
                total=len(work_items), desc="🧬 Aligning gene families", unit="gene"
            ) as pbar:

                # Use imap_unordered for better performance and memory efficiency
                results = pool.imap_unordered(align_single_gene_family, work_items)

                # Collect results as they complete
                for result in results:
                    if result is not None:
                        gene_name, alignment_file_path = result
                        successful_alignments[gene_name] = alignment_file_path
                    pbar.update(1)

        # Report processing summary
        success_count = len(successful_alignments)
        total_count = len(gene_groups)
        success_rate = (success_count / total_count * 100) if total_count > 0 else 0

        print(f"\n🎯 Alignment processing summary:")
        print(
            f"   ✅ Successful: {success_count}/{total_count} gene families ({success_rate:.1f}%)"
        )

        if success_count < total_count:
            failed_count = total_count - success_count
            print(f"   ❌ Failed: {failed_count} gene families")

            # Show which gene families failed
            failed_genes = set(gene_groups.keys()) - set(successful_alignments.keys())
            if failed_genes:
                failed_list = sorted(list(failed_genes))
                print(
                    f"   🔍 Failed genes: {', '.join(failed_list[:5])}{'...' if len(failed_list) > 5 else ''}"
                )

        return successful_alignments

    except FileNotFoundError:
        print(f"❌ Error: Protein manifest not found: {protein_manifest_path}")
        return {}

    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON in protein manifest: {str(e)}")
        return {}

    except Exception as e:
        print(f"❌ Unexpected error in parallel alignment: {str(e)}")
        return {}


def generate_alignment_manifest(
    input_manifest_path, successful_alignments, output_dir, reference_mode_info
):
    """
    Generates alignment_manifest.json for handoff to Domino 5

    This function creates a standardized JSON manifest that documents all successfully
    aligned gene families, their locations, and metadata for downstream analysis.

    Args:
        input_manifest_path (str): Path to the input protein_manifest.json
        successful_alignments (dict): Dictionary mapping gene_name -> alignment_file_path
        output_dir (str): Directory where alignment files were saved
        reference_mode_info (dict): Reference configuration used for alignments

    Returns:
        str: Path to the generated alignment_manifest.json file
    """
    try:
        # Load the original protein manifest
        with open(input_manifest_path, "r") as f:
            original_manifest = json.load(f)

        # Create new alignment manifest based on original
        alignment_manifest = original_manifest.copy()

        # Update pipeline step and metadata
        alignment_manifest["pipeline_step"] = "Aligner"
        alignment_manifest["generated_by"] = (
            "SubScan Domino 4: The Aligner (High-Performance Parallel)"
        )

        # Add timestamp
        try:
            import datetime

            alignment_manifest["generation_timestamp"] = (
                datetime.datetime.now().isoformat()
            )
        except ImportError:
            alignment_manifest["generation_timestamp"] = "unavailable"

        # Add alignment-specific metadata
        alignment_manifest["alignment_metadata"] = {
            "tool": "WildTypeAligner",
            "reference_mode": reference_mode_info.get("mode", "unknown"),
            "total_gene_families_processed": len(successful_alignments),
            "input_manifest": os.path.abspath(input_manifest_path),
            "output_directory": os.path.abspath(output_dir),
        }

        # Add reference mode specific information
        if reference_mode_info.get("mode") == "user":
            alignment_manifest["alignment_metadata"]["user_reference_directory"] = (
                os.path.abspath(reference_mode_info.get("user_reference_dir", ""))
            )
        elif reference_mode_info.get("mode") == "sepi":
            alignment_manifest["alignment_metadata"]["sepi_species"] = (
                reference_mode_info.get("sepi_species", "")
            )

        # Add the alignments mapping - the key output for Domino 5
        alignment_manifest["alignments"] = {}

        for gene_name, alignment_file_path in successful_alignments.items():
            alignment_entry = {
                "alignment_file_path": os.path.abspath(alignment_file_path),
                "gene_family": gene_name,
                "file_size_bytes": None,
                "alignment_format": "water",
            }

            # Try to get file size and basic stats
            if os.path.isfile(alignment_file_path):
                try:
                    alignment_entry["file_size_bytes"] = os.path.getsize(
                        alignment_file_path
                    )
                except OSError:
                    alignment_entry["file_size_bytes"] = "unknown"

                # Try to count sequences in alignment file
                try:
                    with open(alignment_file_path, "r") as af:
                        content = af.read()
                        # Count sequence entries (basic heuristic for water format)
                        sequence_count = content.count(">")
                        alignment_entry["sequence_count"] = sequence_count
                except (IOError, UnicodeDecodeError):
                    alignment_entry["sequence_count"] = "unknown"

            alignment_manifest["alignments"][gene_name] = alignment_entry

        # Add summary statistics
        alignment_manifest["alignment_summary"] = {
            "total_alignments": len(successful_alignments),
            "gene_families_aligned": list(successful_alignments.keys()),
            "alignment_success_rate": "100%" if successful_alignments else "0%",
            "pipeline_ready_for_domino_5": len(successful_alignments) > 0,
        }

        # Write the alignment manifest
        manifest_path = os.path.join(output_dir, "alignment_manifest.json")

        with open(manifest_path, "w") as f:
            json.dump(alignment_manifest, f, indent=2)

        print(f"📄 Generated alignment manifest: {os.path.basename(manifest_path)}")
        print(f"   📁 Location: {os.path.dirname(manifest_path)}")
        print(f"   🧬 Gene families documented: {len(successful_alignments)}")
        print(f"   🔗 Ready for Domino 5: {'Yes' if successful_alignments else 'No'}")

        return os.path.abspath(manifest_path)

    except FileNotFoundError:
        print(f"❌ Error: Input manifest not found: {input_manifest_path}")
        return None

    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON in input manifest: {str(e)}")
        return None

    except Exception as e:
        print(f"❌ Error generating alignment manifest: {str(e)}")
        return None


def main():
    """Main entry point for the high-performance Aligner wrapper"""
    parser = create_argument_parser()
    args = parser.parse_args()

    # Validate arguments
    validate_arguments(args)

    print("\n" + "=" * 60)
    print("SubScan Pipeline - Domino 4: The Aligner")
    print("=" * 60)
    print(f"✓ Protein manifest: {args.manifest}")
    print(f"✓ Output directory: {args.output_dir}")

    if args.user_reference_dir:
        print(f"✓ Reference mode: User-provided")
        print(f"✓ Reference directory: {args.user_reference_dir}")
    else:
        print(f"✓ Reference mode: SEPI automated")
        print(f"✓ Target species: {args.sepi_species}")

    print(f"✓ Parallel processes: {args.threads}")
    print("🚀 High-Performance Mode: ENABLED")
    print("=" * 60)

    # Prepare reference mode configuration
    reference_mode_info = {}
    if args.user_reference_dir:
        reference_mode_info = {
            "mode": "user",
            "user_reference_dir": args.user_reference_dir,
        }
    else:
        reference_mode_info = {"mode": "sepi", "sepi_species": args.sepi_species}

    # Execute the parallel alignment processing
    successful_alignments = execute_parallel_alignment(
        args.manifest, reference_mode_info, args.output_dir, args.threads
    )

    # Report final results
    if successful_alignments:
        print(f"\n🎉 Alignment processing completed!")
        print(f"   ✅ {len(successful_alignments)} gene families aligned")
        print(f"   📁 Results saved to: {args.output_dir}")

        # Show sample of aligned gene families
        aligned_genes = sorted(list(successful_alignments.keys()))
        print(
            f"   🧬 Aligned genes: {', '.join(aligned_genes[:5])}{'...' if len(aligned_genes) > 5 else ''}"
        )

        # Generate alignment manifest for pipeline integration
        manifest_path = generate_alignment_manifest(
            args.manifest, successful_alignments, args.output_dir, reference_mode_info
        )
        print(f"   📋 Alignment manifest: {manifest_path}")

    else:
        print(f"\n⚠️  No gene families were successfully aligned")
        print(f"   🔍 Check WildTypeAligner installation and reference files")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
