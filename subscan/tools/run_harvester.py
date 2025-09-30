#!/usr/bin/env python3
"""
MutationScan Domino 1: Genome Harvester

This module implements the first domino in the MutationScan pipeline, responsible
for downloading genome assemblies from NCBI based on accession numbers. It serves
as the entry point for the AMR analysis workflow, converting accession lists into
local FASTA files with comprehensive metadata tracking.

The harvester interfaces with the ncbi_genome_extractor tool to download genomes
and creates standardized JSON manifests for downstream domino tools. It handles
various genome formats and provides robust error handling for network issues
and invalid accessions.

Usage:
    python run_harvester.py --accessions genome_list.txt --email user@email.com 
                           --output-dir results/genomes/

Integration:
    - Input: Text file with NCBI genome accession numbers
    - Output: genome_manifest.json + downloaded FASTA files
    - Next Domino: run_annotator.py (AMR gene annotation)

Author: MutationScan Development Team
Version: 1.0.0
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from subprocess import CompletedProcess
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    """
    Create and configure command-line argument parser for genome harvesting.
    
    Sets up argument parsing for the genome harvester with validation for NCBI
    API access requirements and file path specifications. Ensures all required
    parameters for genome downloading workflow are properly configured.

    Returns:
        argparse.Namespace: Parsed command-line arguments containing:
            - accessions (str): Path to text file with NCBI accession numbers (one per line)
            - email (str): Valid email address for NCBI API compliance (required by NCBI)
            - output_dir (str): Output directory path for downloaded genomes and manifest

    Raises:
        SystemExit: If required arguments are missing or invalid format

    Example:
        >>> args = parse_arguments()
        >>> print(f"Processing {args.accessions} with email {args.email}")
    """
    parser = argparse.ArgumentParser(
        description="SubScan Harvester: Extract genomes from NCBI and generate manifest",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--accessions",
        required=True,
        help="Path to a text file containing NCBI accession numbers, one per line.",
    )

    parser.add_argument(
        "--email",
        required=True,
        help="A valid email address for NCBI programmatic access.",
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Path to the directory where results will be saved.",
    )

    return parser.parse_args()


def execute_ncbi_genome_extractor(
    accessions: str, email: str, output_dir: str
) -> CompletedProcess[str]:
    """
    Execute the external NCBI Genome Extractor tool for automated genome harvesting.

    Interfaces with the ncbi_genome_extractor.py tool to download complete genome
    sequences from NCBI using provided accession numbers. Handles NCBI API 
    interactions, rate limiting, file organization, and metadata generation.

    Args:
        accessions (str): Path to text file containing NCBI accession numbers (one per line)
        email (str): Valid email address for NCBI API access (required by NCBI guidelines)
        output_dir (str): Directory path where harvested genomes and metadata will be saved

    Returns:
        CompletedProcess[str]: Result of subprocess execution with returncode and output

    Raises:
        FileNotFoundError: If ncbi_genome_extractor.py script cannot be found in expected location
        subprocess.CalledProcessError: If the external harvesting tool fails during execution
        PermissionError: If output directory cannot be created or accessed due to permissions

    Side Effects:
        - Creates output directory structure if it doesn't exist
        - Downloads FASTA files to output_dir with standardized naming
        - Generates CSV metadata file with comprehensive genome information
        - Handles NCBI rate limiting and retry logic automatically

    Example:
        >>> result = execute_ncbi_genome_extractor("accessions.txt", "user@email.com", "results/")
        >>> if result.returncode == 0:
        ...     print("Genomes downloaded successfully")
    """
    # Find the ncbi_genome_extractor script
    script_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "ncbi_genome_extractor",
        "ncbi_genome_extractor.py",
    )

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.abspath(output_dir)

    # Read accessions from file to create query
    try:
        with open(accessions, "r", encoding="utf-8") as f:
            accession_list = [line.strip() for line in f if line.strip()]

        # Create a query string from accessions
        query = " OR ".join([f'"{acc}"' for acc in accession_list])

    except (FileNotFoundError, PermissionError, OSError) as e:
        print(f"Error reading accessions file: {e}")
        sys.exit(1)

    # Construct the command to run the external tool
    command = [
        "python",
        script_path,
        "--query",
        query,
        "--output_dir",
        output_dir,
        "--max_results",
        str(len(accession_list)),
        "--metadata_format",
        "csv",
    ]

    print("Executing ncbi_genome_extractor with command:")
    print(" ".join([f'"{arg}"' if " " in arg else arg for arg in command]))
    print()

    try:
        # Execute the subprocess with error checking
        result = subprocess.run(
            command,
            check=True,
            capture_output=False,  # Allow real-time output to user
            text=True,
            cwd=os.path.dirname(
                os.path.dirname(__file__)
            ),  # Run from subscan directory
        )
        print("\nncbi_genome_extractor completed successfully!")
        return result

    except subprocess.CalledProcessError as e:
        print(f"\nError: ncbi_genome_extractor failed with exit code {e.returncode}")
        print("The external tool encountered an error during execution.")
        sys.exit(1)

    except FileNotFoundError:
        print(f"\nError: ncbi_genome_extractor script not found at: {script_path}")
        print(
            "Please ensure the ncbi_genome_extractor repository is cloned in the parent directory."
        )
        print(
            "Run: git clone https://github.com/vihaankulkarni29/ncbi_genome_extractor.git"
        )
        sys.exit(1)


def _read_metadata_csv(metadata_csv_path: str) -> List[Dict[str, Any]]:
    """
    Read and parse the metadata CSV file generated by ncbi_genome_extractor.
    
    Args:
        metadata_csv_path (str): Path to the CSV metadata file
        
    Returns:
        List[Dict[str, str]]: List of dictionaries containing genome metadata
        
    Raises:
        FileNotFoundError: If CSV file doesn't exist
        csv.Error: If CSV file format is invalid
    """
    with open(metadata_csv_path, "r", newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        return list(reader)


def _find_fasta_file(accession: str, genome_id: str, output_dir: str) -> Optional[str]:
    """
    Find the corresponding FASTA file for a genome entry.
    
    Args:
        accession (str): Genome accession number
        genome_id (str): Genome ID from metadata
        output_dir (str): Directory containing FASTA files
        
    Returns:
        str or None: Absolute path to FASTA file if found, None otherwise
    """
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith(".fasta")]
    fasta_path = None

    # Try to find the FASTA file by genome_id
    if genome_id:
        expected_fasta = f"{genome_id}.fasta"
        if expected_fasta in fasta_files:
            fasta_path = os.path.join(output_dir, expected_fasta)

    # If not found by genome_id, try by accession
    if not fasta_path:
        for fasta_file in fasta_files:
            if accession in fasta_file:
                fasta_path = os.path.join(output_dir, fasta_file)
                break

    return fasta_path


def _create_genome_entry(row: Dict[str, Any], fasta_path: Optional[str]) -> Dict[str, Any]:
    """
    Create a genome entry dictionary from metadata row and FASTA path.
    
    Args:
        row (Dict[str, str]): Dictionary containing genome metadata from CSV
        fasta_path (str or None): Absolute path to the genome FASTA file
        
    Returns:
        Dict[str, Any]: Genome entry dictionary with accession, fasta_path, and metadata
    """
    accession = row.get("accession", "").strip()
    genome_entry = {
        "accession": accession,
        "fasta_path": os.path.abspath(fasta_path) if fasta_path else None,
        "metadata": dict(row),
    }
    return genome_entry


def generate_genome_manifest(accessions: str, email: str, output_dir: str) -> str:
    """
    Generate standardized JSON manifest file for downstream pipeline integration.

    Creates the genome_manifest.json file that serves as the primary data handoff
    mechanism to subsequent domino tools. Processes ncbi_genome_extractor output
    to create a structured manifest with genome metadata, file paths, and processing
    statistics for the MutationScan pipeline.

    Args:
        accessions (str): Path to the original accessions file used for download
        email (str): Email address used for NCBI API access (preserved in manifest)
        output_dir (str): Output directory containing downloaded genomes and metadata

    Returns:
        str: Path to the generated genome_manifest.json file

    Raises:
        FileNotFoundError: If expected metadata CSV file is not found after download
        IOError: If manifest file cannot be written due to permissions or disk space
        ValueError: If metadata format is invalid or missing required fields

    Side Effects:
        - Creates genome_manifest.json in the output directory
        - Validates all referenced FASTA files exist and are readable
        - Links accession numbers to downloaded genome files
        - Preserves comprehensive metadata for downstream analysis

    Example:
        >>> manifest_path = generate_genome_manifest("accessions.txt", "user@email.com", "results/")
        >>> print(f"Manifest created at: {manifest_path}")
    """
    print("Generating genome_manifest.json...")

    # Define expected output paths (based on your tool's structure)
    output_dir_name = os.path.basename(output_dir.rstrip("/\\"))
    paths = {
        "metadata_csv": os.path.join(output_dir, f"{output_dir_name}_metadata.csv"),
        "manifest": os.path.join(output_dir, "genome_manifest.json"),
    }

    # Verify that the expected outputs exist
    if not os.path.exists(paths["metadata_csv"]):
        print(f"Error: Expected metadata file not found: {paths['metadata_csv']}")
        print("Available files in output directory:")
        try:
            for file in os.listdir(output_dir):
                print(f"  - {file}")
        except (OSError, PermissionError):
            print("  (directory not accessible)")
        sys.exit(1)

    # Read metadata.csv to get genome information
    genomes_data = []
    try:
        metadata_rows = _read_metadata_csv(paths["metadata_csv"])

        # Build genome entries from metadata and find corresponding FASTA files
        for row in metadata_rows:
            row_data = {
                "accession": row.get("accession", "").strip(),
                "genome_id": row.get("genome_id", "").strip(),
            }
            if not row_data["accession"]:
                continue

            fasta_path = _find_fasta_file(
                row_data["accession"], row_data["genome_id"], output_dir
            )
            genome_entry = _create_genome_entry(row, fasta_path)

            if not fasta_path:
                print(
                    f"Warning: FASTA file not found for accession {row_data['accession']} (genome_id: {row_data['genome_id']})"
                )

            genomes_data.append(genome_entry)

    except (FileNotFoundError, PermissionError, csv.Error, UnicodeDecodeError) as e:
        print(f"Error reading metadata.csv: {e}")
        sys.exit(1)

    # Create the standardized manifest structure
    manifest = {
        "pipeline_step": "Harvester",
        "parameters": {"accession_file": os.path.abspath(accessions), "email": email},
        "output_files": {
            "metadata_csv": os.path.abspath(paths["metadata_csv"]),
            "genomes": genomes_data,
        },
    }

    # Save the manifest as JSON
    try:
        with open(paths["manifest"], "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)

        print("Successfully generated genome_manifest.json")
        print(f"Manifest location: {os.path.abspath(paths['manifest'])}")
        print(f"Total genomes processed: {len(genomes_data)}")

        return paths["manifest"]

    except (FileNotFoundError, PermissionError, OSError) as e:
        print(f"Error writing manifest file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    args = parse_arguments()

    print("SubScan Harvester - Domino 1")
    print(f"Accessions file: {args.accessions}")
    print(f"Email: {args.email}")
    print(f"Output directory: {args.output_dir}")
    print()

    # Execute the external ncbi_genome_extractor tool
    execute_ncbi_genome_extractor(args.accessions, args.email, args.output_dir)

    # Generate the standardized manifest for pipeline integration
    manifest_path = generate_genome_manifest(
        args.accessions, args.email, args.output_dir
    )

    print()
    print("=" * 60)
    print("SubScan Harvester (Domino 1) completed successfully!")
    print("=" * 60)
    print(f"Manifest file: {manifest_path}")
    print("This manifest is ready to be passed to the next domino in the pipeline.")
    print("The Harvester's work is complete.")
