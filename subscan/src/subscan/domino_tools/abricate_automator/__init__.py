#!/usr/bin/env python3
"""
ABRicate Automator - Integrated version for MutationScan
"""

import os
import sys
import argparse
import subprocess
import glob
import re
import logging


def get_accession(fasta_file):
    """Extract accession from FASTA header."""
    try:
        with open(fasta_file, 'r') as f:
            line = f.readline().strip()
            if line.startswith('>'):
                header_content = line[1:].strip()
                if header_content:
                    potential_accession = header_content.split()[0]
                    if potential_accession and is_valid_accession(potential_accession):
                        return potential_accession
                return os.path.splitext(os.path.basename(fasta_file))[0]
            else:
                return os.path.splitext(os.path.basename(fasta_file))[0]
    except Exception as e:
        logging.warning(f"Error parsing {fasta_file}: {e}")
        return os.path.splitext(os.path.basename(fasta_file))[0]


def is_valid_accession(accession):
    """Check if the string looks like a valid NCBI accession."""
    pattern = r'^[A-Z]{2,4}_[A-Z0-9]+(\.[0-9]+)?$'
    return bool(re.match(pattern, accession))


def find_abricate_path():
    """Find the full path to abricate executable."""
    possible_paths = [
        'abricate',  # In PATH
        'wsl abricate',  # WSL installation
        '/usr/local/bin/abricate',  # System install
        os.path.expanduser('~/miniconda3/bin/abricate'),  # Miniconda
        os.path.expanduser('~/miniconda3/envs/abricate-env/bin/abricate'),  # Conda env
        os.path.expanduser('~/anaconda3/bin/abricate'),  # Anaconda
    ]

    for path in possible_paths:
        try:
            result = subprocess.run(path.split() + ['--version'], capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                return path
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            continue

    return None


def run_abricate(fasta_file, output_file, db='card'):
    """Run ABRicate on a FASTA file."""
    abricate_path = find_abricate_path()
    if not abricate_path:
        logging.error("ABRicate not found. Please install ABRicate and ensure it's in PATH.")
        return False

    # Handle WSL paths
    if abricate_path.startswith('wsl '):
        wsl_fasta = fasta_file.replace('C:', '/mnt/c').replace('\\', '/')
        wsl_output = output_file.replace('C:', '/mnt/c').replace('\\', '/')
        cmd = abricate_path.split() + ['--db', db, wsl_fasta]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                with open(output_file, 'w') as out:
                    out.write(result.stdout)
                return True
            else:
                logging.error(f"ABRicate failed for {fasta_file}: {result.stderr.strip()}")
                return False
        except Exception as e:
            logging.error(f"Error running ABRicate on {fasta_file}: {e}")
            return False
    else:
        # Regular execution
        cmd = [abricate_path, '--db', db, fasta_file]

        try:
            with open(output_file, 'w') as out:
                result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                logging.error(f"ABRicate failed for {fasta_file}: {result.stderr.strip()}")
                return False
            return True
        except Exception as e:
            logging.error(f"Error running ABRicate on {fasta_file}: {e}")
            return False


def process_directory(input_dir, output_dir, db='card'):
    """Process all FASTA files in a directory with ABRicate."""
    # Find FASTA files
    fasta_files = glob.glob(os.path.join(input_dir, '*.fasta')) + \
                  glob.glob(os.path.join(input_dir, '*.fa')) + \
                  glob.glob(os.path.join(input_dir, '*.fna'))

    if not fasta_files:
        logging.warning(f"No FASTA files found in {input_dir}")
        return []

    logging.info(f"Found {len(fasta_files)} FASTA files.")
    os.makedirs(output_dir, exist_ok=True)

    successful = []
    for fasta_file in fasta_files:
        accession = get_accession(fasta_file)
        output_file = os.path.join(output_dir, f"{accession}.tsv")
        logging.info(f"Processing {os.path.basename(fasta_file)} -> {accession}.tsv")
        
        if run_abricate(fasta_file, output_file, db):
            successful.append(output_file)
        else:
            logging.error(f"Failed to process {fasta_file}")

    logging.info(f"Successfully processed {len(successful)}/{len(fasta_files)} files.")
    return successful


# Main function for backward compatibility
def run_abricate_analysis(input_dir: str, output_dir: str, database: str = 'card') -> list:
    """Simple function interface for ABRicate analysis"""
    return process_directory(input_dir, output_dir, database)


if __name__ == '__main__':
    # Simple CLI interface
    parser = argparse.ArgumentParser(description='Run ABRicate on multiple genome FASTA files')
    parser.add_argument('--input-dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('--output-dir', required=True, help='Directory to save .tsv output files')
    parser.add_argument('--db', default='card', help='ABRicate database (default: card)')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    run_abricate_analysis(args.input_dir, args.output_dir, args.db)