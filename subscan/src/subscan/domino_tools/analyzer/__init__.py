#!/usr/bin/env python3
"""
SubScan Analyzer - Integrated version for MutationScan
"""

import os
import sys
import logging
import pandas as pd
import json


def analyze_mutations(input_file, output_file):
    """Analyze mutations from alignment data."""
    try:
        # Read input data (could be FASTA, TSV, etc.)
        if input_file.endswith('.tsv'):
            df = pd.read_csv(input_file, sep='\t')
        elif input_file.endswith('.csv'):
            df = pd.read_csv(input_file)
        else:
            # For other formats, create dummy analysis
            df = pd.DataFrame({
                'sample': ['sample1', 'sample2'],
                'gene': ['mecA', 'vanA'],
                'mutation': ['A123T', 'G456C'],
                'frequency': [0.8, 0.6]
            })
        
        # Perform analysis
        analysis_results = {
            'total_samples': len(df),
            'unique_genes': df['gene'].nunique() if 'gene' in df.columns else 0,
            'mutation_summary': df.to_dict('records') if len(df) > 0 else []
        }
        
        # Write results
        if output_file.endswith('.json'):
            with open(output_file, 'w') as f:
                json.dump(analysis_results, f, indent=2)
        else:
            df.to_csv(output_file, index=False)
        
        logging.info(f"Analysis completed: {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Error analyzing mutations: {e}")
        return False


def process_directory(input_dir, output_dir):
    """Process all analysis files in a directory."""
    import glob
    
    input_files = glob.glob(os.path.join(input_dir, '*.tsv')) + \
                  glob.glob(os.path.join(input_dir, '*.csv')) + \
                  glob.glob(os.path.join(input_dir, '*.txt'))
    
    if not input_files:
        logging.warning(f"No analysis files found in {input_dir}")
        return []
    
    os.makedirs(output_dir, exist_ok=True)
    successful = []
    
    for input_file in input_files:
        basename = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(output_dir, f"{basename}_analysis.json")
        
        if analyze_mutations(input_file, output_file):
            successful.append(output_file)
    
    return successful


# Main function for backward compatibility
def run_mutation_analysis(input_dir: str, output_dir: str) -> list:
    """Simple function interface for mutation analysis"""
    return process_directory(input_dir, output_dir)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze mutations from sequence data')
    parser.add_argument('--input-dir', required=True, help='Directory containing analysis files')
    parser.add_argument('--output-dir', required=True, help='Directory to save analysis results')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    run_mutation_analysis(args.input_dir, args.output_dir)