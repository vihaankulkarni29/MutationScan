#!/usr/bin/env python3
"""
FASTA AA Extractor - Integrated version for MutationScan
"""

import os
import sys
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def extract_protein_sequences(input_file, output_file, gene_list=None):
    """Extract protein sequences from a FASTA file."""
    try:
        sequences = []
        with open(input_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # If gene_list is provided, filter sequences
                if gene_list:
                    if any(gene.lower() in record.description.lower() for gene in gene_list):
                        sequences.append(record)
                else:
                    sequences.append(record)
        
        # Write filtered sequences
        with open(output_file, 'w') as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")
        
        logging.info(f"Extracted {len(sequences)} sequences to {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Error extracting sequences from {input_file}: {e}")
        return False


def process_directory(input_dir, output_dir, gene_list=None):
    """Process all FASTA files in a directory."""
    import glob
    
    fasta_files = glob.glob(os.path.join(input_dir, '*.fasta')) + \
                  glob.glob(os.path.join(input_dir, '*.fa')) + \
                  glob.glob(os.path.join(input_dir, '*.fna'))
    
    if not fasta_files:
        logging.warning(f"No FASTA files found in {input_dir}")
        return []
    
    os.makedirs(output_dir, exist_ok=True)
    successful = []
    
    for fasta_file in fasta_files:
        basename = os.path.splitext(os.path.basename(fasta_file))[0]
        output_file = os.path.join(output_dir, f"{basename}_proteins.fasta")
        
        if extract_protein_sequences(fasta_file, output_file, gene_list):
            successful.append(output_file)
    
    return successful


# Main function for backward compatibility
def extract_aa_sequences(input_dir: str, output_dir: str, genes: list = None) -> list:
    """Simple function interface for protein extraction"""
    return process_directory(input_dir, output_dir, genes)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract protein sequences from FASTA files')
    parser.add_argument('--input-dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('--output-dir', required=True, help='Directory to save protein FASTA files')
    parser.add_argument('--genes', nargs='*', help='List of genes to filter for')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    extract_aa_sequences(args.input_dir, args.output_dir, args.genes)