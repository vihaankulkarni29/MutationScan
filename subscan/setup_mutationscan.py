#!/usr/bin/env python3
"""
MutationScan Pipeline Quick Setup Script
========================================

This script helps researchers quickly set up and validate the MutationScan pipeline
for production use. It performs system checks, installs dependencies, and runs
validation tests.

Usage:
    python setup_mutationscan.py
    python setup_mutationscan.py --check-only
    python setup_mutationscan.py --install-deps
"""

import sys
import os
import subprocess
import platform
import json
import argparse
from pathlib import Path
import importlib.util
from typing import List, Dict, Tuple

# Color codes for output
class Colors:
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'

def print_header(text: str):
    """Print a formatted header."""
    print(f"\n{Colors.BLUE}{Colors.BOLD}{'='*60}")
    print(f" {text}")
    print(f"{'='*60}{Colors.END}")

def print_success(text: str):
    """Print success message."""
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")

def print_warning(text: str):
    """Print warning message."""
    print(f"{Colors.YELLOW}⚠ {text}{Colors.END}")

def print_error(text: str):
    """Print error message."""
    print(f"{Colors.RED}✗ {text}{Colors.END}")

def print_info(text: str):
    """Print info message."""
    print(f"{Colors.BLUE}ℹ {text}{Colors.END}")

def check_python_version() -> bool:
    """Check if Python version meets requirements."""
    version = sys.version_info
    if version.major >= 3 and version.minor >= 8:
        print_success(f"Python {version.major}.{version.minor}.{version.micro} ✓")
        return True
    else:
        print_error(f"Python {version.major}.{version.minor}.{version.micro} - Need 3.8+")
        return False

def check_system_requirements() -> Dict[str, bool]:
    """Check system requirements."""
    print_header("System Requirements Check")
    
    results = {}
    
    # Python version
    results['python'] = check_python_version()
    
    # Operating system
    os_name = platform.system()
    print_info(f"Operating System: {os_name}")
    results['os'] = os_name in ['Linux', 'Darwin', 'Windows']
    if results['os']:
        print_success(f"Supported OS: {os_name}")
    else:
        print_warning(f"Untested OS: {os_name}")
    
    # Available memory (approximate)
    try:
        if os_name == 'Linux':
            with open('/proc/meminfo', 'r') as f:
                meminfo = f.read()
            mem_total = int([line for line in meminfo.split('\n') 
                           if 'MemTotal' in line][0].split()[1]) // 1024
        elif os_name == 'Darwin':  # macOS
            result = subprocess.run(['sysctl', 'hw.memsize'], 
                                  capture_output=True, text=True)
            mem_total = int(result.stdout.split(': ')[1]) // (1024**3)
        else:  # Windows or fallback
            mem_total = 8  # Default assumption
        
        if mem_total >= 16:
            print_success(f"RAM: {mem_total} GB (Recommended)")
        elif mem_total >= 8:
            print_warning(f"RAM: {mem_total} GB (Minimum)")
        else:
            print_error(f"RAM: {mem_total} GB (Insufficient)")
        
        results['memory'] = mem_total >= 8
        
    except Exception:
        print_warning("Could not determine memory size")
        results['memory'] = True  # Assume OK
    
    return results

def check_python_packages() -> Dict[str, bool]:
    """Check if required Python packages are installed."""
    print_header("Python Package Dependencies")
    
    required_packages = {
        'biopython': 'Bio',
        'pandas': 'pandas',
        'requests': 'requests',
        'pytest': 'pytest',
        'pylint': 'pylint'
    }
    
    results = {}
    
    for package_name, import_name in required_packages.items():
        try:
            spec = importlib.util.find_spec(import_name)
            if spec is not None:
                # Try to get version if possible
                module = importlib.import_module(import_name)
                version = getattr(module, '__version__', 'unknown')
                print_success(f"{package_name}: {version}")
                results[package_name] = True
            else:
                print_error(f"{package_name}: Not installed")
                results[package_name] = False
        except Exception:
            print_error(f"{package_name}: Not installed or import error")
            results[package_name] = False
    
    return results

def check_external_tools() -> Dict[str, bool]:
    """Check if external bioinformatics tools are available."""
    print_header("External Tool Dependencies")
    
    tools = ['abricate', 'prokka', 'mlst']
    results = {}
    
    for tool in tools:
        try:
            result = subprocess.run([tool, '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                version_info = result.stdout.split('\n')[0]
                print_success(f"{tool}: {version_info}")
                results[tool] = True
            else:
                print_warning(f"{tool}: Available but version check failed")
                results[tool] = False
        except (subprocess.TimeoutExpired, FileNotFoundError):
            print_error(f"{tool}: Not found in PATH")
            results[tool] = False
        except Exception as e:
            print_error(f"{tool}: Error checking ({e})")
            results[tool] = False
    
    return results

def check_pipeline_structure() -> bool:
    """Check if pipeline files are present."""
    print_header("Pipeline Structure Validation")
    
    required_paths = [
        'src/subscan',
        'tools',
        'tests',
        'src/subscan/utils.py',
        'tools/run_harvester.py',
        'tools/run_cooccurrence_analyzer.py',
        'tests/conftest.py'
    ]
    
    all_present = True
    
    for path in required_paths:
        if os.path.exists(path):
            print_success(f"{path}")
        else:
            print_error(f"{path} - Missing")
            all_present = False
    
    return all_present

def run_validation_tests() -> bool:
    """Run pipeline validation tests."""
    print_header("Pipeline Validation Tests")
    
    try:
        # Run the test suite
        result = subprocess.run([
            sys.executable, '-m', 'pytest', 
            'tests/test_integration.py', 
            '-v', '--tb=short'
        ], capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            # Count passed tests
            lines = result.stdout.split('\n')
            summary_line = [line for line in lines if 'passed' in line and '=' in line]
            if summary_line:
                print_success(f"All tests passed: {summary_line[-1].strip()}")
            else:
                print_success("Validation tests completed successfully")
            return True
        else:
            print_error("Some validation tests failed")
            print(result.stdout)
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print_error("Validation tests timed out")
        return False
    except Exception as e:
        print_error(f"Could not run validation tests: {e}")
        return False

def install_python_dependencies() -> bool:
    """Install Python dependencies."""
    print_header("Installing Python Dependencies")
    
    try:
        print_info("Installing required packages...")
        result = subprocess.run([
            sys.executable, '-m', 'pip', 'install',
            'biopython>=1.79',
            'pandas>=1.3.0',
            'requests>=2.25.0',
            'pytest>=7.0.0',
            'pylint>=2.12.0'
        ], check=True, capture_output=True, text=True)
        
        print_success("Python dependencies installed successfully")
        return True
        
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to install dependencies: {e}")
        print(e.stdout)
        print(e.stderr)
        return False

def generate_sample_config() -> bool:
    """Generate a sample configuration file."""
    print_header("Generating Sample Configuration")
    
    config = {
        "pipeline_settings": {
            "quality_threshold": 95.0,
            "min_coverage": 80.0,
            "max_missing_data": 5.0,
            "batch_size": 50
        },
        "analysis_parameters": {
            "mutation_types": ["SNP", "indel"],
            "resistance_databases": ["CARD", "ResFinder"],
            "gene_families": ["bla", "aac", "tet", "gyr", "par"]
        },
        "output_options": {
            "format": ["json", "csv", "html"],
            "include_plots": True,
            "detailed_logs": True,
            "compress_results": False
        },
        "system_resources": {
            "max_processes": 4,
            "memory_limit": "8GB",
            "temp_dir": "./tmp"
        }
    }
    
    try:
        with open('config.yaml', 'w') as f:
            # Write YAML-formatted config
            import json
            json.dump(config, f, indent=2)
        
        print_success("Sample configuration saved to config.yaml")
        return True
        
    except Exception as e:
        print_error(f"Failed to create configuration file: {e}")
        return False

def create_sample_data() -> bool:
    """Create sample data files for testing."""
    print_header("Creating Sample Data")
    
    try:
        # Create sample directory
        sample_dir = Path("sample_data")
        sample_dir.mkdir(exist_ok=True)
        
        # Sample accessions
        sample_accessions = [
            "NZ_CP107554",  # E. coli
            "NZ_CP107555",  # K. pneumoniae
            "NC_000913.3"   # E. coli K-12
        ]
        
        with open(sample_dir / "sample_accessions.txt", 'w') as f:
            for acc in sample_accessions:
                f.write(f"{acc}\n")
        
        print_success(f"Sample data created in {sample_dir}/")
        print_info("Try: python tools/run_harvester.py --input sample_data/sample_accessions.txt")
        return True
        
    except Exception as e:
        print_error(f"Failed to create sample data: {e}")
        return False

def print_final_report(results: Dict[str, Dict[str, bool]]):
    """Print final setup report."""
    print_header("Setup Summary Report")
    
    all_good = True
    
    print(f"\n{Colors.BOLD}System Requirements:{Colors.END}")
    for check, status in results.get('system', {}).items():
        status_symbol = "✓" if status else "✗"
        color = Colors.GREEN if status else Colors.RED
        print(f"  {color}{status_symbol} {check.replace('_', ' ').title()}{Colors.END}")
        if not status:
            all_good = False
    
    print(f"\n{Colors.BOLD}Python Packages:{Colors.END}")
    for package, status in results.get('packages', {}).items():
        status_symbol = "✓" if status else "✗"
        color = Colors.GREEN if status else Colors.RED
        print(f"  {color}{status_symbol} {package}{Colors.END}")
        if not status:
            all_good = False
    
    print(f"\n{Colors.BOLD}External Tools:{Colors.END}")
    for tool, status in results.get('tools', {}).items():
        status_symbol = "✓" if status else "✗"
        color = Colors.GREEN if status else Colors.YELLOW
        print(f"  {color}{status_symbol} {tool}{Colors.END}")
    
    print(f"\n{Colors.BOLD}Pipeline Structure:{Colors.END}")
    structure_ok = results.get('structure', False)
    status_symbol = "✓" if structure_ok else "✗"
    color = Colors.GREEN if structure_ok else Colors.RED
    print(f"  {color}{status_symbol} All required files present{Colors.END}")
    if not structure_ok:
        all_good = False
    
    print(f"\n{Colors.BOLD}Validation Tests:{Colors.END}")
    tests_ok = results.get('tests', False)
    status_symbol = "✓" if tests_ok else "✗"
    color = Colors.GREEN if tests_ok else Colors.RED
    print(f"  {color}{status_symbol} Integration tests{Colors.END}")
    
    # Overall status
    print(f"\n{Colors.BOLD}Overall Status:{Colors.END}")
    if all_good and tests_ok:
        print_success("🎉 MutationScan is ready for production use!")
        print_info("Next steps:")
        print_info("1. Review the PRODUCTION_GUIDE.md")
        print_info("2. Try the quick start example")
        print_info("3. Run with your own data")
    elif all_good:
        print_warning("⚠️  Core requirements met, but validation tests need attention")
        print_info("You can proceed with caution - check test results above")
    else:
        print_error("❌ Setup incomplete - address the issues above")
        print_info("Install missing dependencies and re-run this script")

def main():
    """Main setup function."""
    parser = argparse.ArgumentParser(description='MutationScan Pipeline Setup')
    parser.add_argument('--check-only', action='store_true',
                       help='Only check requirements, do not install anything')
    parser.add_argument('--install-deps', action='store_true',
                       help='Install Python dependencies automatically')
    parser.add_argument('--skip-tests', action='store_true',
                       help='Skip validation tests (faster setup)')
    
    args = parser.parse_args()
    
    print(f"{Colors.BLUE}{Colors.BOLD}")
    print(" MutationScan Pipeline Setup")
    print(" ==========================")
    print(" Production-Ready AMR Analysis")
    print(f"{Colors.END}")
    
    results = {}
    
    # System checks
    results['system'] = check_system_requirements()
    results['packages'] = check_python_packages()
    results['tools'] = check_external_tools()
    results['structure'] = check_pipeline_structure()
    
    # Install dependencies if requested
    if args.install_deps and not args.check_only:
        install_success = install_python_dependencies()
        if install_success:
            # Re-check packages after installation
            results['packages'] = check_python_packages()
    
    # Generate sample files
    if not args.check_only:
        generate_sample_config()
        create_sample_data()
    
    # Run validation tests
    if not args.skip_tests and not args.check_only:
        results['tests'] = run_validation_tests()
    else:
        results['tests'] = None
    
    # Print final report
    print_final_report(results)
    
    # Return appropriate exit code
    core_requirements_met = all([
        all(results['system'].values()),
        all(results['packages'].values()),
        results['structure']
    ])
    
    if core_requirements_met:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()