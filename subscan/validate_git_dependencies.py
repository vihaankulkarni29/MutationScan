#!/usr/bin/env python3
"""
Dependency Validation Script for MutationScan

This script validates that all Git-based dependencies are properly pip-installable
before attempting to install the main MutationScan package.

Usage:
    python validate_git_dependencies.py
    python validate_git_dependencies.py --fix  # Attempt to fix issues
    python validate_git_dependencies.py --verbose  # Detailed output

Author: MutationScan Development Team
"""

import argparse
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path
from typing import List, Dict, Tuple, Optional

# Git-based dependencies from setup.cfg
GIT_DEPENDENCIES = {
    "federated-genome-extractor": {
        "url": "https://github.com/vihaankulkarni29/federated_genome_extractor.git",
        "required_files": ["pyproject.toml", "setup.py", "__init__.py"],
        "import_name": "federated_genome_extractor",
    },
    "abricate-automator": {
        "url": "https://github.com/vihaankulkarni29/ABRicate-Automator.git",
        "required_files": ["pyproject.toml", "setup.py", "__init__.py"],
        "import_name": "abricate_automator",
    },
    "fasta-aa-extractor": {
        "url": "https://github.com/vihaankulkarni29/FastaAAExtractor.git",
        "required_files": ["pyproject.toml", "setup.py", "__init__.py"],
        "import_name": "fasta_aa_extractor",
    },
    "wildtype-aligner": {
        "url": "https://github.com/vihaankulkarni29/wildtype-aligner.git",
        "required_files": ["pyproject.toml", "setup.py", "__init__.py"],
        "import_name": "wildtype_aligner",
    },
    "subscan-analyzer": {
        "url": "https://github.com/vihaankulkarni29/SubScan.git",
        "required_files": ["setup.py", "__init__.py"],
        "import_name": "subscan_analyzer",
    },
}


class Colors:
    """ANSI color codes for terminal output"""
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BLUE = "\033[94m"
    BOLD = "\033[1m"
    END = "\033[0m"


def print_header(text: str) -> None:
    """Print formatted header"""
    print(f"\n{Colors.BLUE}{Colors.BOLD}{'='*70}")
    print(f" {text}")
    print(f"{'='*70}{Colors.END}")


def print_success(text: str) -> None:
    """Print success message"""
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")


def print_warning(text: str) -> None:
    """Print warning message"""
    print(f"{Colors.YELLOW}⚠ {text}{Colors.END}")


def print_error(text: str) -> None:
    """Print error message"""
    print(f"{Colors.RED}✗ {text}{Colors.END}")


def print_info(text: str) -> None:
    """Print info message"""
    print(f"{Colors.BLUE}ℹ {text}{Colors.END}")


def check_git_installed() -> bool:
    """Check if git is installed"""
    try:
        result = subprocess.run(
            ["git", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def check_pip_installed() -> bool:
    """Check if pip is installed"""
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def clone_repository(url: str, temp_dir: Path, verbose: bool = False) -> Tuple[bool, Optional[Path]]:
    """Clone a Git repository to temporary directory"""
    try:
        clone_path = temp_dir / Path(url).stem
        
        if verbose:
            print_info(f"Cloning {url}...")
        
        result = subprocess.run(
            ["git", "clone", "--depth=1", url, str(clone_path)],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0:
            return True, clone_path
        else:
            if verbose:
                print_error(f"Clone failed: {result.stderr}")
            return False, None
            
    except subprocess.TimeoutExpired:
        if verbose:
            print_error("Clone timed out")
        return False, None
    except Exception as e:
        if verbose:
            print_error(f"Clone error: {e}")
        return False, None


def check_required_files(repo_path: Path, required_files: List[str], verbose: bool = False) -> Tuple[bool, List[str]]:
    """Check if required packaging files exist"""
    missing_files = []
    
    for file in required_files:
        file_path = repo_path / file
        if not file_path.exists():
            missing_files.append(file)
            if verbose:
                print_warning(f"Missing file: {file}")
    
    return len(missing_files) == 0, missing_files


def test_pip_install(url: str, temp_dir: Path, verbose: bool = False) -> Tuple[bool, str]:
    """Test if package can be pip installed"""
    try:
        venv_path = temp_dir / "test_venv"
        
        if verbose:
            print_info(f"Creating test environment...")
        
        # Create virtual environment
        result = subprocess.run(
            [sys.executable, "-m", "venv", str(venv_path)],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            return False, f"Failed to create venv: {result.stderr}"
        
        # Determine pip executable path
        if sys.platform == "win32":
            pip_exe = venv_path / "Scripts" / "pip.exe"
        else:
            pip_exe = venv_path / "bin" / "pip"
        
        if not pip_exe.exists():
            return False, f"Pip not found in venv: {pip_exe}"
        
        if verbose:
            print_info(f"Testing pip install from {url}...")
        
        # Try pip install
        result = subprocess.run(
            [str(pip_exe), "install", f"git+{url}"],
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode == 0:
            return True, "Successfully installed"
        else:
            error_msg = result.stderr if result.stderr else result.stdout
            return False, f"Pip install failed: {error_msg[:200]}"
            
    except subprocess.TimeoutExpired:
        return False, "Pip install timed out"
    except Exception as e:
        return False, f"Pip install error: {str(e)}"


def validate_dependency(
    name: str,
    config: Dict,
    temp_dir: Path,
    verbose: bool = False,
    test_install: bool = True
) -> Dict:
    """Validate a single Git-based dependency"""
    
    result = {
        "name": name,
        "url": config["url"],
        "status": "unknown",
        "issues": [],
        "warnings": [],
    }
    
    print(f"\n{Colors.BOLD}Validating: {name}{Colors.END}")
    print(f"Repository: {config['url']}")
    
    # Step 1: Clone repository
    clone_success, repo_path = clone_repository(config["url"], temp_dir, verbose)
    
    if not clone_success:
        result["status"] = "failed"
        result["issues"].append("Failed to clone repository")
        print_error("✗ Clone failed")
        return result
    
    print_success("✓ Repository cloned")
    
    # Step 2: Check required files
    files_ok, missing_files = check_required_files(
        repo_path,
        config["required_files"],
        verbose
    )
    
    if not files_ok:
        result["status"] = "failed"
        result["issues"].append(f"Missing required files: {', '.join(missing_files)}")
        print_error(f"✗ Missing files: {', '.join(missing_files)}")
        return result
    
    print_success(f"✓ All required files present: {', '.join(config['required_files'])}")
    
    # Step 3: Test pip install (optional, slower)
    if test_install:
        install_ok, install_msg = test_pip_install(config["url"], temp_dir, verbose)
        
        if not install_ok:
            result["status"] = "failed"
            result["issues"].append(f"Pip install failed: {install_msg}")
            print_error(f"✗ Pip install test failed")
            if verbose:
                print_error(f"   {install_msg}")
            return result
        
        print_success("✓ Pip install test passed")
    
    # All checks passed
    result["status"] = "passed"
    print_success(f"✓ {name} validation PASSED")
    
    return result


def main():
    """Main validation function"""
    parser = argparse.ArgumentParser(
        description="Validate MutationScan Git-based dependencies",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output with detailed error messages"
    )
    
    parser.add_argument(
        "--no-install-test",
        action="store_true",
        help="Skip pip install testing (faster, less thorough)"
    )
    
    parser.add_argument(
        "--dependency", "-d",
        choices=list(GIT_DEPENDENCIES.keys()),
        help="Test only a specific dependency"
    )
    
    args = parser.parse_args()
    
    print_header("MutationScan Git Dependency Validator")
    
    # Check prerequisites
    print("\nChecking prerequisites...")
    
    if not check_git_installed():
        print_error("Git is not installed or not in PATH")
        sys.exit(1)
    print_success("Git is installed")
    
    if not check_pip_installed():
        print_error("Pip is not installed or not available")
        sys.exit(1)
    print_success("Pip is installed")
    
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="mutationscan_dep_check_"))
    
    if args.verbose:
        print_info(f"Using temporary directory: {temp_dir}")
    
    try:
        # Select dependencies to test
        if args.dependency:
            deps_to_test = {args.dependency: GIT_DEPENDENCIES[args.dependency]}
        else:
            deps_to_test = GIT_DEPENDENCIES
        
        print_header(f"Validating {len(deps_to_test)} Dependencies")
        
        results = []
        for name, config in deps_to_test.items():
            result = validate_dependency(
                name,
                config,
                temp_dir,
                verbose=args.verbose,
                test_install=not args.no_install_test
            )
            results.append(result)
        
        # Print summary
        print_header("Validation Summary")
        
        passed = sum(1 for r in results if r["status"] == "passed")
        failed = sum(1 for r in results if r["status"] == "failed")
        
        print(f"\nTotal dependencies tested: {len(results)}")
        print_success(f"Passed: {passed}")
        
        if failed > 0:
            print_error(f"Failed: {failed}")
        
        # Print failed dependencies
        if failed > 0:
            print("\n" + Colors.BOLD + "Failed Dependencies:" + Colors.END)
            for result in results:
                if result["status"] == "failed":
                    print(f"\n{Colors.RED}✗ {result['name']}{Colors.END}")
                    for issue in result["issues"]:
                        print(f"  - {issue}")
        
        # Exit with appropriate code
        if failed > 0:
            print(f"\n{Colors.RED}Validation FAILED - {failed} dependencies have issues{Colors.END}")
            sys.exit(1)
        else:
            print(f"\n{Colors.GREEN}All dependencies validated successfully!{Colors.END}")
            sys.exit(0)
    
    finally:
        # Cleanup
        if temp_dir.exists():
            if args.verbose:
                print_info(f"Cleaning up temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
