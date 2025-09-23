#!/usr/bin/env python3
"""
SubScan Pipeline - Domino 2: The Annotator

This script serves as an integration wrapper for the ABRicate-Automator tool,
acting as the second domino in the SubScan bioinformatics pipeline.

Purpose:
- Reads genome_manifest.json from Domino 1 (The Harvester)
- Executes ABRicate-Automator for AMR gene annotation in parallel        # St            # Step 3: Execute ABRicate-Automator
            print("\n🚀 Step 3: Running ABRicate-Automator...") 2: Create temporary directory with symlinks
        pri        print("\n🔗 Step 2: Validating input files...")t("\n🔗 Step 2: Preparing input files...")- Generates annotation_manifest.json for Domino 3 (FastaAAExtractor)

Author: SubScan Pipeline Development Team
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import multiprocessing
from pathlib import Path
import glob
from datetime import datetime
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm


def create_argument_parser():
    """Create and configure the command-line argument parser"""
    parser = argparse.ArgumentParser(
        description="SubScan Domino 2: AMR Gene Annotation Wrapper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_annotator.py --manifest genome_manifest.json --output-dir ./annotations
  python run_annotator.py --manifest results/genome_manifest.json --output-dir ./amr_results --threads 8

Pipeline Flow:
  Domino 1 (Harvester) → genome_manifest.json → Domino 2 (Annotator) → annotation_manifest.json → Domino 3 (FastaAAExtractor)
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--manifest',
        type=str,
        required=True,
        help='Path to the JSON manifest file from the Harvester'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Path to the directory where Abricate results will be saved'
    )
    
    # Optional arguments
    parser.add_argument(
        '--threads',
        type=int,
        default=multiprocessing.cpu_count(),
        help=f'Number of genomes to process in parallel (default: {multiprocessing.cpu_count()} CPU cores)'
    )
    
    return parser


def validate_arguments(args):
    """Validate command-line arguments"""
    # Check if manifest file exists
    if not os.path.isfile(args.manifest):
        print(f"Error: Manifest file not found: {args.manifest}", file=sys.stderr)
        sys.exit(1)
    
    # Check if manifest file is JSON
    if not args.manifest.lower().endswith('.json'):
        print(f"Warning: Manifest file should be a JSON file: {args.manifest}", file=sys.stderr)
    
    # Validate thread count
    if args.threads < 1:
        print(f"Error: Thread count must be at least 1, got: {args.threads}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"✓ Manifest file: {args.manifest}")
    print(f"✓ Output directory: {args.output_dir}")
    print(f"✓ Parallel threads: {args.threads}")


def load_genome_manifest(manifest_path):
    """Load and parse the genome manifest JSON file"""
    try:
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
        
        print(f"✓ Loaded manifest: {len(manifest.get('genomes', []))} genomes found")
        
        # Validate manifest structure
        if 'genomes' not in manifest:
            raise ValueError("Manifest missing 'genomes' key")
        
        if not isinstance(manifest['genomes'], list):
            raise ValueError("'genomes' must be a list")
        
        # Validate genome entries
        for i, genome in enumerate(manifest['genomes']):
            if not isinstance(genome, dict):
                raise ValueError(f"Genome entry {i} must be a dictionary")
            
            if 'fasta_path' not in genome:
                raise ValueError(f"Genome entry {i} missing 'fasta_path'")
            
            if not os.path.isfile(genome['fasta_path']):
                print(f"Warning: FASTA file not found: {genome['fasta_path']}")
        
        return manifest
        
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in manifest file: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading manifest: {e}", file=sys.stderr)
        sys.exit(1)


def create_symlink_directory(genomes, temp_dir):
    """Create symbolic links to genome FASTA files in a temporary directory"""
    print(f"\n🔗 Creating symbolic links in temporary directory: {temp_dir}")
    
    symlink_count = 0
    for i, genome in enumerate(genomes):
        fasta_path = genome['fasta_path']
        
        if not os.path.isfile(fasta_path):
            print(f"⚠️  Skipping missing file: {fasta_path}")
            continue
        
        # Create symlink with original filename
        original_filename = os.path.basename(fasta_path)
        symlink_path = os.path.join(temp_dir, original_filename)
        
        try:
            # Convert to absolute paths for reliable symlinks
            abs_fasta_path = os.path.abspath(fasta_path)
            
            # Create symlink (works on Windows with proper permissions)
            if os.name == 'nt':  # Windows
                # On Windows, try creating symlink, fallback to copy if needed
                try:
                    os.symlink(abs_fasta_path, symlink_path)
                except OSError:
                    # If symlink fails (permissions), use hard link or copy
                    import shutil
                    shutil.copy2(abs_fasta_path, symlink_path)
                    print(f"    📄 Copied (fallback): {original_filename}")
            else:  # Unix-like systems
                os.symlink(abs_fasta_path, symlink_path)
            
            symlink_count += 1
            print(f"    🔗 Linked: {original_filename}")
            
        except Exception as e:
            print(f"⚠️  Failed to link {original_filename}: {e}")
    
    print(f"✓ Created {symlink_count} file links")
    return symlink_count


def execute_abricate_automator(input_dir, output_dir, threads):
    """Execute the ABRicate-Automator tool via subprocess"""
    print(f"\n🚀 Executing ABRicate-Automator...")
    print(f"    Input directory: {input_dir}")
    print(f"    Output directory: {output_dir}")
    print(f"    Threads: {threads}")
    
    # Construct the command
    command = [
        'abricate-automator',
        '--input-dir', input_dir,
        '--output-dir', output_dir,
        '--threads', str(threads)
    ]
    
    print(f"    Command: {' '.join(command)}")
    
    try:
        # Execute the command with real-time output
        print("\n" + "="*50)
        print("ABRicate-Automator Output:")
        print("="*50)
        
        result = subprocess.run(
            command,
            check=True,
            capture_output=False,  # Allow real-time output
            text=True
        )
        
        print("="*50)
        print("✅ ABRicate-Automator completed successfully!")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ ABRicate-Automator failed with exit code {e.returncode}")
        print(f"Error: {e}")
        return False
        
    except FileNotFoundError:
        print("\n❌ ABRicate-Automator command not found!")
        print("Please ensure 'abricate-automator' is installed and available in PATH")
        return False
        
    except Exception as e:
        print(f"\n❌ Unexpected error running ABRicate-Automator: {e}")
        return False


def annotate_single_genome(genome_data, output_dir, abricate_command_base):
    """Worker function to annotate a single genome using ABRicate-Automator"""
    genome_id = genome_data.get('accession', 'unknown')
    fasta_path = genome_data.get('fasta_path', '')
    
    if not os.path.isfile(fasta_path):
        return {
            'genome_id': genome_id,
            'success': False,
            'error': f"FASTA file not found: {fasta_path}",
            'card_results_path': None
        }
    
    try:
        # Create individual temporary directory for this genome
        with tempfile.TemporaryDirectory() as temp_genome_dir:
            # Create symlink to FASTA file in temp directory
            original_filename = os.path.basename(fasta_path)
            symlink_path = os.path.join(temp_genome_dir, original_filename)
            
            # Create symlink (with Windows fallback)
            try:
                abs_fasta_path = os.path.abspath(fasta_path)
                if os.name == 'nt':  # Windows
                    try:
                        os.symlink(abs_fasta_path, symlink_path)
                    except OSError:
                        import shutil
                        shutil.copy2(abs_fasta_path, symlink_path)
                else:  # Unix-like systems
                    os.symlink(abs_fasta_path, symlink_path)
            except Exception as e:
                return {
                    'genome_id': genome_id,
                    'success': False,
                    'error': f"Failed to create symlink: {e}",
                    'card_results_path': None
                }
            
            # Create individual output directory for this genome
            genome_output_dir = os.path.join(output_dir, f"temp_{genome_id}")
            os.makedirs(genome_output_dir, exist_ok=True)
            
            # Construct command for this specific genome
            command = abricate_command_base + [
                '--input-dir', temp_genome_dir,
                '--output-dir', genome_output_dir,
                '--threads', '1'  # Single thread per worker to avoid conflicts
            ]
            
            # Execute ABRicate-Automator for this genome
            result = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout per genome
            )
            
            # Find the generated CARD results file
            expected_card_file = f"{genome_id}_card.tsv"
            card_path = os.path.join(genome_output_dir, expected_card_file)
            
            if os.path.isfile(card_path):
                # Move result file to main output directory
                final_card_path = os.path.join(output_dir, expected_card_file)
                import shutil
                shutil.move(card_path, final_card_path)
                
                # Clean up temporary genome output directory
                shutil.rmtree(genome_output_dir, ignore_errors=True)
                
                return {
                    'genome_id': genome_id,
                    'success': True,
                    'error': None,
                    'card_results_path': os.path.abspath(final_card_path)
                }
            else:
                return {
                    'genome_id': genome_id,
                    'success': False,
                    'error': f"Expected output file not found: {expected_card_file}",
                    'card_results_path': None
                }
                
    except subprocess.TimeoutExpired:
        return {
            'genome_id': genome_id,
            'success': False,
            'error': "ABRicate-Automator timed out (5 minutes)",
            'card_results_path': None
        }
    except subprocess.CalledProcessError as e:
        return {
            'genome_id': genome_id,
            'success': False,
            'error': f"ABRicate-Automator failed: {e.stderr}",
            'card_results_path': None
        }
    except Exception as e:
        return {
            'genome_id': genome_id,
            'success': False,
            'error': f"Unexpected error: {str(e)}",
            'card_results_path': None
        }


def execute_parallel_annotation(genomes, output_dir, threads):
    """Execute ABRicate-Automator in parallel across multiple genomes"""
    print(f"\n🚀 Step 3: Running parallel ABRicate-Automator annotation...")
    print(f"    Output directory: {output_dir}")
    print(f"    Parallel workers: {threads}")
    print(f"    Genomes to process: {len(genomes)}")
    
    # Prepare base command for ABRicate-Automator (without input/output specific args)
    abricate_command_base = ['abricate-automator']
    
    # Test if abricate-automator is available
    try:
        subprocess.run(['abricate-automator', '--help'], 
                      capture_output=True, check=True, timeout=10)
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        print("❌ ABRicate-Automator command not found or not working!")
        print("Please ensure 'abricate-automator' is installed and available in PATH")
        return []
    
    # Create worker function with fixed parameters
    worker_func = partial(annotate_single_genome, 
                         output_dir=output_dir, 
                         abricate_command_base=abricate_command_base)
    
    # Execute parallel annotation
    print("\n" + "="*50)
    print("Parallel ABRicate-Automator Execution:")
    print("="*50)
    
    results = []
    try:
        with Pool(processes=threads) as pool:
            # Use tqdm for progress tracking
            results = list(tqdm(
                pool.imap(worker_func, genomes),
                total=len(genomes),
                desc="Annotating Genomes",
                unit="genome"
            ))
    except Exception as e:
        print(f"❌ Parallel execution failed: {e}")
        return []
    
    # Analyze results
    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]
    
    print("="*50)
    print(f"✅ Parallel annotation completed!")
    print(f"📊 Results: {len(successful)} successful, {len(failed)} failed")
    
    if failed:
        print("\n⚠️  Failed annotations:")
        for failure in failed[:5]:  # Show first 5 failures
            print(f"    - {failure['genome_id']}: {failure['error']}")
        if len(failed) > 5:
            print(f"    ... and {len(failed) - 5} more failures")
    
    return results


def find_abricate_results(output_dir, genomes):
    """Find and map ABRicate result files to genome entries"""
    print(f"\n🔍 Step 4: Locating ABRicate result files...")
    
    # Look for TSV files in the output directory
    tsv_pattern = os.path.join(output_dir, "*.tsv")
    found_files = glob.glob(tsv_pattern)
    
    print(f"Found {len(found_files)} TSV files in output directory")
    
    # Create mapping of accession to result files
    result_mapping = {}
    
    for genome in genomes:
        accession = genome.get('accession', '')
        if not accession:
            # Try to extract accession from FASTA filename
            fasta_path = genome.get('fasta_path', '')
            if fasta_path:
                filename = os.path.basename(fasta_path)
                # Remove .fasta extension to get accession
                accession = filename.replace('.fasta', '').replace('.fa', '')
        
        if accession:
            # Look for corresponding CARD results file
            card_pattern = f"{accession}_card.tsv"
            matching_files = [f for f in found_files if os.path.basename(f) == card_pattern]
            
            if matching_files:
                # Use absolute path for maximum portability
                abs_path = os.path.abspath(matching_files[0])
                result_mapping[accession] = abs_path
                print(f"    ✓ Found: {card_pattern}")
            else:
                print(f"    ⚠️  Missing: {card_pattern}")
                result_mapping[accession] = None
    
    return result_mapping


def generate_annotation_manifest(original_manifest, parallel_results, output_dir):
    """Generate the annotation manifest for Domino 3 using parallel results"""
    print(f"\n📄 Step 4: Generating annotation manifest from parallel results...")
    
    # Create a deep copy of the original manifest
    annotation_manifest = json.loads(json.dumps(original_manifest))
    
    # Update pipeline metadata
    annotation_manifest['pipeline_step'] = 'Annotator'
    annotation_manifest['generated_by'] = 'SubScan Domino 2: The Annotator (Parallel)'
    annotation_manifest['generation_timestamp'] = datetime.now().isoformat()
    
    # Calculate statistics from parallel results
    total_genomes = len(parallel_results)
    successful_annotations = sum(1 for r in parallel_results if r['success'])
    
    annotation_manifest['annotation_stats'] = {
        'total_genomes': total_genomes,
        'successful_annotations': successful_annotations,
        'failed_annotations': total_genomes - successful_annotations,
        'success_rate': f"{(successful_annotations/total_genomes*100):.1f}%" if total_genomes > 0 else "0%",
        'parallel_execution': True,
        'processing_method': 'multiprocessing.Pool'
    }
    
    # Create mapping from genome_id to results
    results_mapping = {r['genome_id']: r for r in parallel_results}
    
    # Update each genome entry with AMR CARD results
    updated_count = 0
    for genome in annotation_manifest['genomes']:
        accession = genome.get('accession', '')
        
        # If no accession in manifest, try to extract from FASTA filename
        if not accession:
            fasta_path = genome.get('fasta_path', '')
            if fasta_path:
                filename = os.path.basename(fasta_path)
                accession = filename.replace('.fasta', '').replace('.fa', '')
                # Add the extracted accession to the genome entry
                genome['accession'] = accession
        
        # Add AMR CARD results path from parallel results
        if accession and accession in results_mapping:
            result = results_mapping[accession]
            if result['success'] and result['card_results_path']:
                genome['amr_card_results'] = result['card_results_path']
                updated_count += 1
                print(f"    ✅ Added AMR results for: {accession}")
            else:
                genome['amr_card_results'] = ""
                genome['annotation_error'] = result.get('error', 'Unknown error')
                print(f"    ❌ Failed annotation for: {accession} - {result.get('error', 'Unknown error')}")
        else:
            genome['amr_card_results'] = ""
            print(f"    ⚠️  No parallel result found for: {accession}")
    
    print(f"✅ Updated {updated_count} genome entries with AMR results")
    
    # Save the annotation manifest
    manifest_path = os.path.join(output_dir, "annotation_manifest.json")
    
    try:
        with open(manifest_path, 'w') as f:
            json.dump(annotation_manifest, f, indent=2)
        
        print(f"✅ Annotation manifest saved: {manifest_path}")
        print(f"📊 Parallel Processing Statistics:")
        print(f"    - Total genomes: {total_genomes}")
        print(f"    - Successful annotations: {successful_annotations}")
        print(f"    - Success rate: {annotation_manifest['annotation_stats']['success_rate']}")
        print(f"    - Processing method: Parallel (multiprocessing.Pool)")
        
        return manifest_path
        
    except Exception as e:
        print(f"❌ Failed to save annotation manifest: {e}")
        return None


def main():
    """Main entry point for the Annotator wrapper"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Validate arguments
    validate_arguments(args)
    
    print("\n" + "="*60)
    print("SubScan Pipeline - Domino 2: The Annotator")
    print("="*60)
    print(f"Input manifest: {args.manifest}")
    print(f"Output directory: {args.output_dir}")
    print(f"Parallel processing: {args.threads} threads")
    print("="*60)
    
    # Complete parallel processing workflow
    try:
        # Step 1: Load and validate the genome manifest
        print("\n📋 Step 1: Loading genome manifest...")
        manifest = load_genome_manifest(args.manifest)
        
        genomes = manifest['genomes']
        if not genomes:
            print("❌ No genomes found in manifest")
            return 1
        
        # Step 2: Create temporary directory with symlinks
        print("\n� Step 2: Preparing input files...")
        with tempfile.TemporaryDirectory() as temp_dir:
            print(f"Temporary directory: {temp_dir}")
            
            # Create symlinks to all FASTA files
            symlink_count = create_symlink_directory(genomes, temp_dir)
            
            if symlink_count == 0:
                print("❌ No valid FASTA files found to process")
                return 1
            
            # Step 3: Execute ABRicate-Automator
            print("\n� Step 3: Running ABRicate-Automator...")
            success = execute_abricate_automator(temp_dir, args.output_dir, args.threads)
            
            if not success:
                print("❌ ABRicate-Automator execution failed")
                return 1
        
        # Step 4: Find ABRicate result files
        result_mapping = find_abricate_results(args.output_dir, genomes)
        
        # Step 5: Generate annotation manifest
        manifest_path = generate_annotation_manifest(manifest, result_mapping, args.output_dir)
        
        if not manifest_path:
            print("❌ Failed to generate annotation manifest")
            return 1
        
        # Final success message
        print("\n" + "="*60)
        print("🎉 DOMINO 2 COMPLETE: The Annotator (High-Performance)")
        print("="*60)
        print("✅ All phases completed successfully with parallel processing!")
        print(f"📄 Input manifest: {args.manifest}")
        print(f"📁 ABRicate results: {args.output_dir}")
        print(f"📋 Output manifest: {manifest_path}")
        print(f"⚡ Performance: {args.threads} parallel workers")
        print("\n🔄 Pipeline Status:")
        print("    Domino 1 (Harvester) → ✅ Complete")
        print("    Domino 2 (Annotator) → ✅ Complete (Parallel)")
        print("    Domino 3 (FastaAAExtractor) → 🚀 Ready to proceed")
        print("="*60)
        
        return 0
        
    except KeyboardInterrupt:
        print("\n⚠️  Process interrupted by user")
        return 1
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())