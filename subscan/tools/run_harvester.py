#!/usr/bin/env python3
"""
MutationScan Domino 1: Genome Harvester

This module implements the first domino in the MutationScan pipeline, responsible
for downloading genome assemblies from NCBI based on accession numbers. It serves
as the entry point for the AMR analysis workflow, converting accession lists into
local FASTA files with comprehensive metadata tracking.

The harvester interfaces with the federated_genome_extractor tool to download genomes
from multiple databases (NCBI, BV-BRC, EnteroBase, PATRIC) and creates standardized 
JSON manifests for downstream domino tools. It handles various genome formats and 
provides robust error handling for network issues and invalid accessions.

Usage:
    python run_harvester.py --accessions genome_list.txt --email user@email.com 
                           --output-dir results/genomes/ [--database ncbi|all]

Integration:
    - Input: Text file with genome accession numbers (one per line)
    - Output: genome_manifest.json + downloaded FASTA files with source database metadata
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

    parser.add_argument(
        "--database",
        choices=["ncbi", "bvbrc", "enterobase", "patric", "all"],
        default="all",
        help="Database source(s) to query for genomes (default: all). "
             "Use 'ncbi' for NCBI-only, or 'all' for federated multi-database search.",
    )

    return parser.parse_args()


def execute_federated_genome_extractor(
    accessions: str, email: str, output_dir: str, database: str = "all"
) -> CompletedProcess[str]:
    """
    Execute the federated genome extractor for multi-database genome harvesting.

    Interfaces with the federated_genome_extractor/harvester.py tool to download complete genome
    sequences from NCBI, BV-BRC, EnteroBase, or PATRIC databases using provided accession numbers.
    Handles multi-database API interactions, rate limiting, file organization, and metadata generation
    with source database tracking.

    Args:
        accessions (str): Path to text file containing genome accession numbers (one per line)
        email (str): Valid email address for API access (required by NCBI guidelines)
        output_dir (str): Directory path where harvested genomes and metadata will be saved
        database (str): Database source - 'ncbi', 'bvbrc', 'enterobase', 'patric', or 'all' (default: 'all')

    Returns:
        CompletedProcess[str]: Result of subprocess execution with returncode and output

    Raises:
        FileNotFoundError: If federated_genome_extractor/harvester.py script cannot be found in expected location
        subprocess.CalledProcessError: If the external harvesting tool fails during execution
        PermissionError: If output directory cannot be created or accessed due to permissions

    Side Effects:
        - Creates output directory structure if it doesn't exist
        - Downloads FASTA files to output_dir with standardized naming
        - Generates JSON/CSV metadata file with comprehensive genome information including source database
        - Handles multi-database rate limiting and retry logic automatically

    Example:
        >>> result = execute_federated_genome_extractor("accessions.txt", "user@email.com", "results/", "all")
        >>> if result.returncode == 0:
        ...     print("Genomes downloaded successfully from federated databases")
    """
    # Find the federated_genome_extractor script
    # Look for it in a few common locations:
    # 1. Adjacent to subscan (for development)
    # 2. In a tools/external directory
    # 3. As a system command (if installed globally)
    
    possible_paths = [
        # Option 1: Adjacent to MutationScan directory (common development setup)
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
            "federated_genome_extractor", 
            "harvester.py"
        ),
        # Option 2: Within subscan external tools directory
        os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "external", 
            "federated_genome_extractor", 
            "harvester.py"
        ),
        # Option 3: Within subscan as submodule
        os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "federated_genome_extractor", 
            "harvester.py"
        ),
    ]
    
    script_path = None
    for path in possible_paths:
        if os.path.exists(path):
            script_path = path
            break
    
    if not script_path:
        print(f"Error: federated_genome_extractor not found in expected locations:")
        for path in possible_paths:
            print(f"  - {path}")
        print()
        print("Please install federated_genome_extractor by:")
        print("1. Clone: git clone https://github.com/vihaankulkarni29/federated_genome_extractor.git")
        print("2. Place it adjacent to MutationScan directory, or")
        print("3. Place it in subscan/external/ directory")
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.abspath(output_dir)

    # Read accessions from file to create query
    try:
        with open(accessions, "r", encoding="utf-8") as f:
            accession_list = [line.strip() for line in f if line.strip()]

        # Create comma-separated list for download_only mode
        query = ",".join(accession_list)

    except (FileNotFoundError, PermissionError, OSError) as e:
        print(f"Error reading accessions file: {e}")
        sys.exit(1)

    # Construct the command to run the federated genome extractor
    command = [
        "python",
        script_path,
        "--source",
        database,
        "--query",
        query,
        "--output_dir",
        output_dir,
        "--download_only",  # Download FASTAs (metadata will be created manually)
        "--metadata_format",
        "csv",
    ]

    print(f"Executing federated genome extractor (database: {database})...")
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
        print("\nFederated genome extractor completed successfully!")
        
        # Create metadata CSV if it doesn't exist (for download_only mode)
        output_dir_name = os.path.basename(output_dir.rstrip("/\\"))
        expected_metadata_path = os.path.join(output_dir, f"{output_dir_name}_metadata.csv")
        
        if not os.path.exists(expected_metadata_path):
            _create_minimal_metadata_csv(output_dir, expected_metadata_path, database)
            
        return result

    except subprocess.CalledProcessError as e:
        print(f"\nError: Federated genome extractor failed with exit code {e.returncode}")
        print("The external tool encountered an error during execution.")
        sys.exit(1)

    except FileNotFoundError:
        print(f"\nError: Federated genome extractor script not found at: {script_path}")
        print(
            "Please ensure the federated_genome_extractor repository is installed via pip."
        )
        print(
            "Run: pip install git+https://github.com/vihaankulkarni29/federated_genome_extractor.git"
        )
        sys.exit(1)


def _create_minimal_metadata_csv(output_dir: str, metadata_csv_path: str, database: str) -> None:
    """
    Create minimal metadata CSV from downloaded FASTA files when federated extractor 
    uses --download_only mode and doesn't generate metadata.
    
    Args:
        output_dir (str): Directory containing downloaded FASTA files
        metadata_csv_path (str): Path where metadata CSV should be created
        database (str): Source database name
    """
    import glob
    
    print(f"Creating minimal metadata CSV from downloaded FASTA files...")
    
    # Find all FASTA files in output directory
    fasta_files = glob.glob(os.path.join(output_dir, "*.fasta")) + \
                  glob.glob(os.path.join(output_dir, "*.fna")) + \
                  glob.glob(os.path.join(output_dir, "*.fa"))
    
    if not fasta_files:
        print("Warning: No FASTA files found to create metadata")
        return
    
    # Create minimal metadata CSV
    with open(metadata_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        writer.writerow(['accession', 'genome_id', 'database', 'organism', 'title'])
        
        # Write minimal data for each FASTA file
        for fasta_file in fasta_files:
            filename = os.path.basename(fasta_file)
            # Extract accession from filename (remove extension)
            accession = os.path.splitext(filename)[0]
            
            # Create minimal metadata row
            writer.writerow([
                accession,           # accession
                accession,           # genome_id (same as accession)
                database.upper(),    # database (NCBI, BV-BRC, etc.)
                f"Unknown organism", # organism (placeholder)
                f"{accession} genome sequence"  # title (placeholder)
            ])
    
    print(f"Created minimal metadata CSV with {len(fasta_files)} entries: {metadata_csv_path}")


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
    print(f"Database source: {args.database}")
    print()

    # Execute the external federated genome extractor tool
    execute_federated_genome_extractor(args.accessions, args.email, args.output_dir, args.database)

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
