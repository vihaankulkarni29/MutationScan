#!/usr/bin/env python3
"""
MutationScan Final Comprehensive Test Suite
Validates pipeline with parameter variations, output quality, and performance
"""

import os
import sys
import json
import time
import shutil
from pathlib import Path
from typing import List, Dict, Tuple, Any

# ANSI color codes for better output
class Colors:
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    END = '\033[0m'

def color_text(text: str, color: str) -> str:
    """Apply color to text for terminal output"""
    return f"{color}{text}{Colors.END}"

def print_header(title: str) -> None:
    """Print a formatted section header"""
    print("\n" + "="*70)
    print(color_text(title, Colors.BLUE + Colors.BOLD))
    print("="*70)

def create_test_accession_files() -> Dict[str, Path]:
    """Create different test accession files for parameter variation"""
    test_dir = Path("comprehensive_test_data")
    test_dir.mkdir(exist_ok=True)
    
    accession_files = {}
    
    # Small set - 2 genomes (existing sample)
    small_file = test_dir / "small_accessions.txt"
    with open(small_file, 'w') as f:
        f.write("GCF_000005825.2\n")  # E. coli
        f.write("GCF_000009605.1\n")  # E. coli
    accession_files['small'] = small_file
    
    # Medium set - 4 genomes  
    medium_file = test_dir / "medium_accessions.txt"
    with open(medium_file, 'w') as f:
        f.write("GCF_000005825.2\n")  # E. coli K-12
        f.write("GCF_000009605.1\n")  # E. coli
        f.write("NC_000913.3\n")      # E. coli K-12 MG1655
        f.write("GCF_000008865.2\n")  # E. coli
    accession_files['medium'] = medium_file
    
    # Large set - 6 genomes (for performance testing)
    large_file = test_dir / "large_accessions.txt"
    with open(large_file, 'w') as f:
        f.write("GCF_000005825.2\n")  # E. coli K-12
        f.write("GCF_000009605.1\n")  # E. coli
        f.write("NC_000913.3\n")      # E. coli K-12 MG1655
        f.write("GCF_000008865.2\n")  # E. coli
        f.write("GCF_000006945.2\n")  # Salmonella
        f.write("GCF_000009685.1\n")  # Klebsiella
    accession_files['large'] = large_file
    
    return accession_files

def create_test_gene_files() -> Dict[str, Path]:
    """Create different gene lists for testing"""
    test_dir = Path("comprehensive_test_data")
    test_dir.mkdir(exist_ok=True)
    
    gene_files = {}
    
    # Basic AMR genes
    basic_file = test_dir / "basic_genes.txt"
    with open(basic_file, 'w') as f:
        f.write("mecA\n")     # MRSA resistance
        f.write("vanA\n")     # Vancomycin resistance
    gene_files['basic'] = basic_file
    
    # Extended AMR panel
    extended_file = test_dir / "extended_genes.txt"
    with open(extended_file, 'w') as f:
        f.write("mecA\n")     # MRSA resistance
        f.write("vanA\n")     # Vancomycin resistance
        f.write("blaTEM\n")   # Beta-lactamase
        f.write("qnrS\n")     # Quinolone resistance
        f.write("tetA\n")     # Tetracycline resistance
    gene_files['extended'] = extended_file
    
    # Comprehensive panel
    comprehensive_file = test_dir / "comprehensive_genes.txt"
    with open(comprehensive_file, 'w') as f:
        f.write("mecA\n")     # MRSA resistance
        f.write("vanA\n")     # Vancomycin resistance  
        f.write("blaTEM\n")   # Beta-lactamase
        f.write("qnrS\n")     # Quinolone resistance
        f.write("tetA\n")     # Tetracycline resistance
        f.write("aph\n")      # Aminoglycoside resistance
        f.write("cat\n")      # Chloramphenicol resistance
        f.write("sul\n")      # Sulfonamide resistance
    gene_files['comprehensive'] = comprehensive_file
    
    return gene_files

def test_parameter_variations() -> bool:
    """Test pipeline with different parameter combinations"""
    print_header("Parameter Variation Testing")
    
    # Create test files
    accession_files = create_test_accession_files()
    gene_files = create_test_gene_files()
    
    # Test combinations
    test_cases = [
        {
            'name': 'Small Dataset + Basic Genes',
            'accessions': accession_files['small'],
            'genes': gene_files['basic'],
            'species': 'Escherichia coli'
        },
        {
            'name': 'Medium Dataset + Extended Genes', 
            'accessions': accession_files['medium'],
            'genes': gene_files['extended'],
            'species': 'Escherichia coli'
        },
        {
            'name': 'Large Dataset + Comprehensive Genes',
            'accessions': accession_files['large'], 
            'genes': gene_files['comprehensive'],
            'species': 'Mixed bacterial species'
        }
    ]
    
    success_count = 0
    
    for i, test_case in enumerate(test_cases, 1):
        test_name = f"Test Case {i}: {test_case['name']}"
        print(f"\n{color_text(test_name, Colors.CYAN)}")
        print(f"  Accessions: {test_case['accessions'].name}")
        print(f"  Genes: {test_case['genes'].name}")
        print(f"  Species: {test_case['species']}")
        
        # Run dry-run test (faster for validation)
        output_dir = f"param_test_{i}"
        
        try:
            import subprocess
            cmd = [
                sys.executable, "run_pipeline.py",
                "--accessions", str(test_case['accessions']),
                "--gene-list", str(test_case['genes']),
                "--email", "test@example.com",
                "--output-dir", output_dir,
                "--sepi-species", test_case['species'],
                "--dry-run"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())
            
            if result.returncode == 0:
                print(f"  {color_text('✓ PASSED', Colors.GREEN)} - Dry-run successful")
                success_count += 1
            else:
                print(f"  {color_text('✗ FAILED', Colors.RED)} - Exit code: {result.returncode}")
                print(f"  Error: {result.stderr[:200]}...")
                
        except Exception as e:
            print(f"  {color_text('✗ ERROR', Colors.RED)} - {str(e)}")
    
    print(f"\n{color_text('Parameter Variation Results:', Colors.BOLD)}")
    print(f"  Passed: {success_count}/{len(test_cases)}")
    
    return success_count == len(test_cases)

def test_database_variations() -> bool:
    """Test different database sources with federated extractor"""
    print_header("Database Source Testing")
    
    databases = ['ncbi', 'all']  # Test NCBI-only vs federated
    
    test_dir = Path("comprehensive_test_data")
    small_accessions = test_dir / "small_accessions.txt"
    basic_genes = test_dir / "basic_genes.txt"
    
    success_count = 0
    
    for db in databases:
        print(f"\n{color_text(f'Testing Database: {db.upper()}', Colors.CYAN)}")
        
        try:
            import subprocess
            cmd = [
                sys.executable, "run_harvester.py", 
                "--accessions", str(small_accessions),
                "--email", "test@example.com",
                "--output-dir", f"db_test_{db}",
                "--database", db
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())
            
            if result.returncode == 0:
                print(f"  {color_text('✓ PASSED', Colors.GREEN)} - Database {db} working")
                success_count += 1
                
                # Clean up
                shutil.rmtree(f"db_test_{db}", ignore_errors=True)
            else:
                print(f"  {color_text('✗ FAILED', Colors.RED)} - Database {db} failed")
                
        except Exception as e:
            print(f"  {color_text('✗ ERROR', Colors.RED)} - {str(e)}")
    
    print(f"\n{color_text('Database Testing Results:', Colors.BOLD)}")
    print(f"  Passed: {success_count}/{len(databases)}")
    
    return success_count == len(databases)

def test_performance_characteristics() -> bool:
    """Test pipeline performance with different dataset sizes"""
    print_header("Performance Characteristics Testing")
    
    accession_files = create_test_accession_files()
    basic_genes = Path("comprehensive_test_data") / "basic_genes.txt"
    
    performance_results = {}
    
    for size, acc_file in accession_files.items():
        print(f"\n{color_text(f'Performance Test: {size.upper()} dataset', Colors.CYAN)}")
        
        start_time = time.time()
        
        try:
            import subprocess
            cmd = [
                sys.executable, "run_pipeline.py",
                "--accessions", str(acc_file),
                "--gene-list", str(basic_genes),
                "--email", "test@example.com", 
                "--output-dir", f"perf_test_{size}",
                "--sepi-species", "Escherichia coli",
                "--dry-run"  # Use dry-run for consistent timing
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())
            end_time = time.time()
            
            duration = end_time - start_time
            performance_results[size] = {
                'success': result.returncode == 0,
                'duration': duration,
                'genome_count': len(open(acc_file).readlines())
            }
            
            if result.returncode == 0:
                print(f"  {color_text('✓ PASSED', Colors.GREEN)} - {duration:.2f}s")
            else:
                print(f"  {color_text('✗ FAILED', Colors.RED)} - {duration:.2f}s")
                
        except Exception as e:
            print(f"  {color_text('✗ ERROR', Colors.RED)} - {str(e)}")
            performance_results[size] = {'success': False, 'duration': 0, 'genome_count': 0}
    
    # Print performance analysis
    print(f"\n{color_text('Performance Analysis:', Colors.BOLD)}")
    for size, result in performance_results.items():
        if result['success']:
            per_genome = result['duration'] / result['genome_count'] if result['genome_count'] > 0 else 0
            print(f"  {size.capitalize():12} | {result['genome_count']:2d} genomes | {result['duration']:5.2f}s total | {per_genome:.2f}s/genome")
    
    # All tests should pass for performance validation
    success_count = sum(1 for r in performance_results.values() if r['success'])
    total_tests = len(performance_results)
    
    print(f"\n{color_text('Performance Results:', Colors.BOLD)}")
    print(f"  Passed: {success_count}/{total_tests}")
    
    return success_count == total_tests

def validate_output_structure() -> bool:
    """Validate the structure of pipeline outputs"""
    print_header("Output Structure Validation")
    
    # Run a single test to generate outputs for validation
    test_dir = Path("comprehensive_test_data")
    small_accessions = test_dir / "small_accessions.txt"
    basic_genes = test_dir / "basic_genes.txt"
    
    print(f"{color_text('Running pipeline to generate outputs...', Colors.CYAN)}")
    
    try:
        import subprocess
        cmd = [
            sys.executable, "run_pipeline.py",
            "--accessions", str(small_accessions),
            "--gene-list", str(basic_genes), 
            "--email", "test@example.com",
            "--output-dir", "output_validation_test",
            "--sepi-species", "Escherichia coli",
            "--dry-run"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())
        
        if result.returncode != 0:
            print(f"  {color_text('✗ Pipeline failed', Colors.RED)}")
            return False
            
        # Find the run directory
        run_dirs = list(Path("output_validation_test").glob("mutationscan_run_*"))
        if not run_dirs:
            print(f"  {color_text('✗ No run directory found', Colors.RED)}")
            return False
            
        run_dir = run_dirs[0]
        
        # Check expected outputs
        expected_outputs = [
            "01_harvester_results/genome_manifest.json",
            "02_annotator_results/annotation_manifest.json", 
            "03_extractor_results/protein_manifest.json",
            "04_aligner_results/alignment_manifest.json",
            "05_analyzer_results/analysis_manifest.json",
            "06_cooccurrence_results/cooccurrence_manifest.json",
            "07_reporter_results/final_report.html"
        ]
        
        validation_results = []
        
        for expected in expected_outputs:
            output_path = run_dir / expected
            exists = output_path.exists()
            
            status = color_text('✓ EXISTS', Colors.GREEN) if exists else color_text('✗ MISSING', Colors.RED)
            print(f"  {expected:45} | {status}")
            
            validation_results.append(exists)
        
        # Clean up
        shutil.rmtree("output_validation_test", ignore_errors=True)
        
        success_count = sum(validation_results)
        total_expected = len(expected_outputs)
        
        print(f"\n{color_text('Output Validation Results:', Colors.BOLD)}")
        print(f"  Found: {success_count}/{total_expected} expected outputs")
        
        return success_count == total_expected
        
    except Exception as e:
        print(f"  {color_text('✗ ERROR', Colors.RED)} - {str(e)}")
        return False

def main():
    """Run comprehensive final testing suite"""
    print(color_text("MutationScan Final Comprehensive Test Suite", Colors.BOLD + Colors.BLUE))
    print("="*70)
    
    # Change to tools directory
    tools_dir = Path(__file__).parent
    os.chdir(tools_dir)
    
    test_results = {}
    
    # Run all test suites
    test_suites = [
        ("Parameter Variations", test_parameter_variations),
        ("Database Sources", test_database_variations), 
        ("Performance Characteristics", test_performance_characteristics),
        ("Output Structure", validate_output_structure)
    ]
    
    for suite_name, test_func in test_suites:
        try:
            result = test_func()
            test_results[suite_name] = result
        except Exception as e:
            print(f"\n{color_text(f'ERROR in {suite_name}: {e}', Colors.RED)}")
            test_results[suite_name] = False
    
    # Final summary
    print_header("Final Comprehensive Test Results")
    
    passed = sum(test_results.values())
    total = len(test_results)
    
    for suite_name, result in test_results.items():
        status = color_text('PASS', Colors.GREEN) if result else color_text('FAIL', Colors.RED)
        print(f"  {suite_name:30} | {status}")
    
    print(f"\n{color_text('OVERALL RESULTS:', Colors.BOLD)}")
    print(f"  Test Suites Passed: {passed}/{total}")
    
    if passed == total:
        print(f"\n{color_text('🎉 ALL COMPREHENSIVE TESTS PASSED!', Colors.GREEN + Colors.BOLD)}")
        print(f"{color_text('The MutationScan pipeline is ready for production deployment!', Colors.GREEN)}")
        return True
    else:
        print(f"\n{color_text('⚠️  Some tests failed. Review results above.', Colors.YELLOW + Colors.BOLD)}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)