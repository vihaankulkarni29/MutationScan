#!/usr/bin/env python3
"""
SubScan Harvester - Integration Wrapper for ncbi_genome_extractor

This module serves as "Domino 1" in the SubScan bioinformatics pipeline.
It wraps the ncbi_genome_extractor tool and produces a standardized
genome_manifest.json file for downstream processing.
"""

import argparse
import csv
import json
import os
import subprocess
import sys


def parse_arguments():
    """Parse command-line arguments for the Harvester wrapper."""
    parser = argparse.ArgumentParser(
        description="SubScan Harvester: Extract genomes from NCBI and generate manifest",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--accessions",
        required=True,
        help="Path to a text file containing NCBI accession numbers, one per line."
    )
    
    parser.add_argument(
        "--email",
        required=True,
        help="A valid email address for NCBI programmatic access."
    )
    
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Path to the directory where results will be saved."
    )
    
    return parser.parse_args()


def execute_ncbi_genome_extractor(accessions, email, output_dir):
    """
    Execute the external ncbi_genome_extractor tool with the provided arguments.
    
    Args:
        accessions (str): Path to the accessions file
        email (str): Email address for NCBI access
        output_dir (str): Output directory path
    
    Raises:
        subprocess.CalledProcessError: If the external tool fails
    """
    # Find the ncbi_genome_extractor script
    script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ncbi_genome_extractor", "ncbi_genome_extractor.py")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.abspath(output_dir)
    
    # Read accessions from file to create query
    try:
        with open(accessions, 'r') as f:
            accession_list = [line.strip() for line in f if line.strip()]
        
        # Create a query string from accessions
        query = " OR ".join([f'"{acc}"' for acc in accession_list])
        
    except Exception as e:
        print(f"Error reading accessions file: {e}")
        sys.exit(1)
    
    # Construct the command to run the external tool
    command = [
        "python", script_path,
        "--query", query,
        "--output_dir", output_dir,
        "--max_results", str(len(accession_list)),
        "--metadata_format", "csv"
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
            cwd=os.path.dirname(os.path.dirname(__file__))  # Run from subscan directory
        )
        print("\nncbi_genome_extractor completed successfully!")
        return result
        
    except subprocess.CalledProcessError as e:
        print(f"\nError: ncbi_genome_extractor failed with exit code {e.returncode}")
        print("The external tool encountered an error during execution.")
        sys.exit(1)
        
    except FileNotFoundError:
        print(f"\nError: ncbi_genome_extractor script not found at: {script_path}")
        print("Please ensure the ncbi_genome_extractor repository is cloned in the parent directory.")
        print("Run: git clone https://github.com/vihaankulkarni29/ncbi_genome_extractor.git")
        sys.exit(1)


def generate_genome_manifest(accessions, email, output_dir):
    """
    Generate the standardized genome_manifest.json file for pipeline integration.
    
    This function reads the output from ncbi_genome_extractor and creates a structured
    manifest that serves as input for the next domino in the SubScan pipeline.
    
    Args:
        accessions (str): Path to the original accessions file
        email (str): Email used for NCBI access
        output_dir (str): Output directory where results were saved
    """
    print("Generating genome_manifest.json...")
    
    # Define expected output paths (based on your tool's structure)
    # The tool creates files with the pattern: {output_dir}_metadata.csv
    output_dir_name = os.path.basename(output_dir.rstrip('/\\'))
    metadata_csv_path = os.path.join(output_dir, f"{output_dir_name}_metadata.csv")
    manifest_path = os.path.join(output_dir, "genome_manifest.json")
    
    # Verify that the expected outputs exist
    if not os.path.exists(metadata_csv_path):
        print(f"Error: Expected metadata file not found: {metadata_csv_path}")
        print("Available files in output directory:")
        try:
            for file in os.listdir(output_dir):
                print(f"  - {file}")
        except:
            print("  (directory not accessible)")
        sys.exit(1)
    
    # Read metadata.csv to get genome information
    genomes_data = []
    try:
        with open(metadata_csv_path, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            metadata_rows = list(reader)
            
        # Build genome entries from metadata and find corresponding FASTA files
        for row in metadata_rows:
            accession = row.get('accession', '').strip()
            genome_id = row.get('genome_id', '').strip()
            if not accession:
                continue
            
            # Look for corresponding FASTA file by genome_id (how your tool names files)
            fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
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
            
            if fasta_path:
                genome_entry = {
                    "accession": accession,
                    "fasta_path": os.path.abspath(fasta_path),
                    "metadata": dict(row)  # Include all metadata columns
                }
                genomes_data.append(genome_entry)
            else:
                print(f"Warning: FASTA file not found for accession {accession} (genome_id: {genome_id})")
                # Still include the entry but without fasta_path
                genome_entry = {
                    "accession": accession,
                    "fasta_path": None,
                    "metadata": dict(row)
                }
                genomes_data.append(genome_entry)
                
    except Exception as e:
        print(f"Error reading metadata.csv: {e}")
        sys.exit(1)
    
    # Create the standardized manifest structure
    manifest = {
        "pipeline_step": "Harvester",
        "parameters": {
            "accession_file": os.path.abspath(accessions),
            "email": email
        },
        "output_files": {
            "metadata_csv": os.path.abspath(metadata_csv_path),
            "genomes": genomes_data
        }
    }
    
    # Save the manifest as JSON
    try:
        with open(manifest_path, 'w', encoding='utf-8') as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)
            
        print(f"Successfully generated genome_manifest.json")
        print(f"Manifest location: {os.path.abspath(manifest_path)}")
        print(f"Total genomes processed: {len(genomes_data)}")
        
        return manifest_path
        
    except Exception as e:
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
    manifest_path = generate_genome_manifest(args.accessions, args.email, args.output_dir)
    
    print()
    print("="*60)
    print("SubScan Harvester (Domino 1) completed successfully!")
    print("="*60)
    print(f"Manifest file: {manifest_path}")
    print("This manifest is ready to be passed to the next domino in the pipeline.")
    print("The Harvester's work is complete.")