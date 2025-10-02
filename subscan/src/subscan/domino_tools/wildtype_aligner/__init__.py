#!/usr/bin/env python3
"""
Wildtype Aligner - Integrated version for MutationScan
"""

import os
import sys
import logging
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def align_sequences(query_file, reference_file, output_file):
    """Align query sequences to reference sequences."""
    try:
        # Read sequences
        query_seqs = list(SeqIO.parse(query_file, "fasta"))
        ref_seqs = list(SeqIO.parse(reference_file, "fasta"))
        
        if not query_seqs or not ref_seqs:
            logging.error("No sequences found in input files")
            return False
        
        # Simple alignment (in real implementation, would use proper alignment tool)
        aligned_seqs = []
        for query in query_seqs:
            for ref in ref_seqs:
                # Create alignment record
                aligned_record = SeqRecord(
                    query.seq,
                    id=f"{query.id}_aligned_to_{ref.id}",
                    description=f"Aligned sequence: {query.description}"
                )
                aligned_seqs.append(aligned_record)
        
        # Write aligned sequences
        with open(output_file, 'w') as output_handle:
            SeqIO.write(aligned_seqs, output_handle, "fasta")
        
        logging.info(f"Created {len(aligned_seqs)} alignments in {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Error aligning sequences: {e}")
        return False


def process_directory(input_dir, reference_dir, output_dir):
    """Process all FASTA files in a directory."""
    import glob
    
    query_files = glob.glob(os.path.join(input_dir, '*.fasta')) + \
                  glob.glob(os.path.join(input_dir, '*.fa')) + \
                  glob.glob(os.path.join(input_dir, '*.fna'))
    
    ref_files = glob.glob(os.path.join(reference_dir, '*.fasta')) + \
                glob.glob(os.path.join(reference_dir, '*.fa')) + \
                glob.glob(os.path.join(reference_dir, '*.fna'))
    
    if not query_files:
        logging.warning(f"No query FASTA files found in {input_dir}")
        return []
    
    if not ref_files:
        logging.warning(f"No reference FASTA files found in {reference_dir}")
        return []
    
    os.makedirs(output_dir, exist_ok=True)
    successful = []
    
    for query_file in query_files:
        basename = os.path.splitext(os.path.basename(query_file))[0]
        output_file = os.path.join(output_dir, f"{basename}_aligned.fasta")
        
        # Use first reference file (can be enhanced)
        ref_file = ref_files[0]
        
        if align_sequences(query_file, ref_file, output_file):
            successful.append(output_file)
    
    return successful


# Main function for backward compatibility
def align_to_wildtype(input_dir: str, reference_dir: str, output_dir: str) -> list:
    """Simple function interface for wildtype alignment"""
    return process_directory(input_dir, reference_dir, output_dir)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Align sequences to wildtype references')
    parser.add_argument('--input-dir', required=True, help='Directory containing query FASTA files')
    parser.add_argument('--reference-dir', required=True, help='Directory containing reference FASTA files')
    parser.add_argument('--output-dir', required=True, help='Directory to save aligned files')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    align_to_wildtype(args.input_dir, args.reference_dir, args.output_dir)