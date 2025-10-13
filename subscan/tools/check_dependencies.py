#!/usr/bin/env python3
"""
Dependency Validator

Validates that all required dependencies for the MutationScan pipeline are available,
including Python version, core libraries, GitHub-based tools, and platform-specific
external tools like ABRicate.

Run this from the repository root:
    python subscan/tools/check_dependencies.py

Exit codes:
    0 - All required dependencies met (warnings allowed for optional/platform-specific)
    1 - One or more critical dependencies missing
"""
import os
import sys
import platform
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple, Optional

# ANSI color codes
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
BLUE = "\033[94m"
RESET = "\033[0m"

def color_text(text: str, color: str) -> str:
    """Wrap text in color codes if terminal supports it."""
    if sys.stdout.isatty():
        return f"{color}{text}{RESET}"
    return text

def get_platform_info() -> dict:
    """Detect platform details."""
    return {
        'system': platform.system(),  # Windows, Linux, Darwin
        'is_windows': platform.system() == 'Windows',
        'is_linux': platform.system() == 'Linux',
        'is_macos': platform.system() == 'Darwin',
        'is_wsl': 'microsoft' in platform.uname().release.lower() if platform.system() == 'Linux' else False,
    }

def check_python_version() -> Tuple[bool, str]:
    """Check if Python version meets requirements (>=3.9)."""
    version = sys.version_info
    version_str = f"{version.major}.{version.minor}.{version.micro}"
    
    if version.major >= 3 and version.minor >= 9:
        return True, f"[OK] Python {version_str} (>= 3.9 required)"
    else:
        return False, f"[FAIL] Python {version_str} (>= 3.9 required, found {version_str})"

def check_import(module_name: str, package_name: Optional[str] = None) -> Tuple[bool, str]:
    """Check if a Python module can be imported."""
    try:
        __import__(module_name)
        return True, f"[OK] {module_name}"
    except ImportError:
        install_name = package_name or module_name
        return False, f"[FAIL] {module_name} (install: pip install {install_name})"

def check_command(cmd: str, args: List[str] = None) -> Tuple[bool, str, str]:
    """Check if a command-line tool is available.
    
    Returns: (success, message, version_output)
    """
    if args is None:
        args = ['--version']
    
    exe_path = shutil.which(cmd)
    if not exe_path:
        return False, f"[FAIL] '{cmd}' not found in PATH", ""
    
    try:
        result = subprocess.run(
            [cmd] + args,
            capture_output=True,
            text=True,
            timeout=5
        )
        
        version_out = (result.stdout + result.stderr).strip()
        if result.returncode == 0:
            # Extract first line for concise display
            first_line = version_out.split('\n')[0] if version_out else "available"
            return True, f"[OK] {cmd}: {first_line}", version_out
        else:
            return False, f"[WARN] '{cmd}' found but returned error", version_out
    except subprocess.TimeoutExpired:
        return False, f"[WARN] '{cmd}' timed out", ""
    except Exception as e:
        return False, f"[FAIL] Error checking '{cmd}': {e}", ""

def check_wsl_command(cmd: str, args: List[str] = None) -> Tuple[bool, str, str]:
    """Check if a command is available via WSL on Windows."""
    if args is None:
        args = ['--version']
    
    try:
        result = subprocess.run(
            ['wsl', cmd] + args,
            capture_output=True,
            text=True,
            timeout=10
        )
        
        version_out = (result.stdout + result.stderr).strip()
        if result.returncode == 0:
            first_line = version_out.split('\n')[0] if version_out else "available"
            return True, f"[OK] {cmd} (via WSL): {first_line}", version_out
        else:
            return False, f"[WARN] WSL '{cmd}' returned error", version_out
    except FileNotFoundError:
        return False, f"[FAIL] WSL not available on this system", ""
    except subprocess.TimeoutExpired:
        return False, f"[WARN] WSL '{cmd}' timed out", ""
    except Exception as e:
        return False, f"[FAIL] Error checking WSL '{cmd}': {e}", ""

def main() -> int:
    """Run all dependency validation checks."""
    print(color_text("\n" + "=" * 70, BLUE))
    print(color_text("MutationScan Dependency Validator", BLUE))
    print(color_text("=" * 70 + "\n", BLUE))
    
    # Detect platform
    plat = get_platform_info()
    print(f"Platform: {plat['system']}")
    if plat['is_wsl']:
        print(color_text("  (WSL detected - Ubuntu on Windows)", BLUE))
    print()
    
    checks: List[Tuple[bool, str]] = []
    warnings: List[str] = []
    critical_failures: List[str] = []
    
    # ========== Python Version ==========
    print(color_text("Checking Python version...", BLUE))
    passed, msg = check_python_version()
    checks.append((passed, msg))
    print(color_text(msg, GREEN if passed else RED))
    if not passed:
        critical_failures.append(msg)
    print()
    
    # ========== Core Python Libraries ==========
    print(color_text("Checking core Python libraries...", BLUE))
    
    core_libs = [
        ('pandas', 'pandas'),
        ('Bio', 'biopython'),
        ('numpy', 'numpy'),
        ('requests', 'requests'),
        ('tqdm', 'tqdm'),
        ('jinja2', 'jinja2'),
        ('plotly', 'plotly'),
        ('openpyxl', 'openpyxl'),
        ('colorama', 'colorama'),
    ]
    
    for module, package in core_libs:
        passed, msg = check_import(module, package)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else RED))
        if not passed:
            critical_failures.append(msg)
    
    print()
    
    # ========== GitHub-based Domino Tools ==========
    print(color_text("Checking GitHub-based pipeline tools...", BLUE))
    
    github_tools = [
        ('federated_genome_extractor', 'Federated Genome Harvester'),
        ('fasta_aa_extractor', 'Protein Extractor'),
        ('wildtype_aligner', 'Sequence Aligner'),
        ('subscan', 'SubScan Core'),
    ]
    
    for module, desc in github_tools:
        passed, msg = check_import(module)
        msg_with_desc = msg.replace(module, f"{module} ({desc})")
        checks.append((passed, msg_with_desc))
        print(color_text(msg_with_desc, GREEN if passed else YELLOW))
        if not passed:
            warnings.append(f"{desc}: {msg}")
    
    # Special check for subscan-analyzer (might be under different import name)
    try:
        # Try common variations
        import subscan_analyzer
        print(color_text("[OK] subscan_analyzer (Mutation Analyzer)", GREEN))
        checks.append((True, "subscan_analyzer imported"))
    except ImportError:
        try:
            import analyzer
            print(color_text("[OK] analyzer (Mutation Analyzer - alternate import)", GREEN))
            checks.append((True, "analyzer imported"))
        except ImportError:
            msg = "[WARN] subscan_analyzer (Mutation Analyzer) - not found"
            print(color_text(msg, YELLOW))
            checks.append((False, msg))
            warnings.append("Mutation Analyzer: install from GitHub")
    
    print()
    
    # ========== FastaAAExtractor CLI ==========
    print(color_text("Checking FastaAAExtractor CLI availability...", BLUE))
    
    extractor_cli = shutil.which('fasta_aa_extractor')
    if extractor_cli:
        print(color_text(f"[OK] fasta_aa_extractor console script: {extractor_cli}", GREEN))
        checks.append((True, "fasta_aa_extractor CLI available"))
    else:
        # Try python -m fallback
        try:
            result = subprocess.run(
                [sys.executable, '-m', 'fasta_aa_extractor', '--help'],
                capture_output=True,
                timeout=5
            )
            if result.returncode == 0:
                print(color_text("[OK] fasta_aa_extractor via 'python -m' (CLI fallback)", GREEN))
                checks.append((True, "fasta_aa_extractor module available"))
            else:
                msg = "[WARN] fasta_aa_extractor CLI not functional"
                print(color_text(msg, YELLOW))
                checks.append((False, msg))
                warnings.append("FastaAAExtractor: reinstall or use mock for testing")
        except Exception:
            msg = "[WARN] fasta_aa_extractor not available (console or module)"
            print(color_text(msg, YELLOW))
            checks.append((False, msg))
            warnings.append("FastaAAExtractor: install from GitHub or use mock")
    
    print()
    
    # ========== ABRicate (Platform-specific) ==========
    print(color_text("Checking ABRicate (Annotator dependency)...", BLUE))
    
    if plat['is_linux'] or plat['is_wsl']:
        # Native Linux or WSL - should have ABRicate natively
        passed, msg, ver = check_command('abricate', ['--version'])
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else YELLOW))
        
        if passed:
            # Check if databases are set up
            passed_db, msg_db, _ = check_command('abricate', ['--list'])
            if passed_db:
                print(color_text("[OK] ABRicate databases accessible", GREEN))
                checks.append((True, "ABRicate databases OK"))
            else:
                warn_msg = "[WARN] ABRicate databases may not be set up (run: abricate --setupdb)"
                print(color_text(warn_msg, YELLOW))
                checks.append((False, warn_msg))
                warnings.append("Run 'abricate --setupdb' to initialize databases")
        else:
            warnings.append("ABRicate not installed. Install: sudo apt-get install abricate")
            print(color_text("  Fix: sudo apt-get install abricate && abricate --setupdb", YELLOW))
    
    elif plat['is_windows']:
        # Windows native - check if WSL is available
        print(color_text("  (Windows detected - checking WSL for ABRicate...)", BLUE))
        passed, msg, ver = check_wsl_command('abricate', ['--version'])
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else YELLOW))
        
        if passed:
            print(color_text("[OK] ABRicate available via WSL - full runs supported", GREEN))
            checks.append((True, "ABRicate via WSL OK"))
        else:
            warn_msg = "[WARN] ABRicate not available (WSL or native)"
            print(color_text(warn_msg, YELLOW))
            checks.append((False, warn_msg))
            warnings.append("For full runs, install WSL + Ubuntu and run: sudo apt-get install abricate")
            print(color_text("  Note: Windows can run dry-runs and use mock annotator", BLUE))
            print(color_text("  For production: Use WSL Ubuntu or native Linux", BLUE))
    
    else:
        # macOS or other
        passed, msg, ver = check_command('abricate', ['--version'])
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else YELLOW))
        if not passed:
            warnings.append("ABRicate: Install via Homebrew (brew install brewsci/bio/abricate)")
    
    print()
    
    # ========== Optional External Tools ==========
    print(color_text("Checking optional external tools...", BLUE))
    
    optional_tools = [
        ('git', ['--version'], "Version control"),
        ('wget', ['--version'], "File download (harvester fallback)"),
        ('curl', ['--version'], "File download (alternative)"),
    ]
    
    for tool, args, desc in optional_tools:
        passed, msg, ver = check_command(tool, args)
        # Don't add to main checks list (optional)
        print(color_text(f"  {msg} - {desc}", GREEN if passed else BLUE))
    
    print()
    
    # ========== Build Tools ==========
    print(color_text("Checking Python build tools...", BLUE))
    
    build_modules = [
        ('pip', None),
        ('setuptools', None),
        ('wheel', None),
    ]
    
    for module, _ in build_modules:
        passed, msg = check_import(module)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else RED))
        if not passed:
            critical_failures.append(msg)
    
    print()
    
    # ========== Summary ==========
    print(color_text("=" * 70, BLUE))
    print(color_text("SUMMARY", BLUE))
    print(color_text("=" * 70, BLUE))
    
    total_checks = len(checks)
    passed_checks = sum(1 for passed, _ in checks if passed)
    failed_checks = total_checks - passed_checks
    
    print(f"\nTotal checks: {total_checks}")
    print(color_text(f"Passed: {passed_checks}", GREEN))
    
    if failed_checks > 0:
        print(color_text(f"Failed/Warnings: {failed_checks}", YELLOW))
    
    # Critical failures
    if critical_failures:
        print()
        print(color_text("CRITICAL FAILURES (must fix):", RED))
        for msg in critical_failures:
            print(f"  - {msg}")
        print()
        print(color_text("Remediation:", YELLOW))
        print("  Install missing core dependencies:")
        print("    pip install -U pip setuptools wheel")
        print("    pip install -e subscan/")
        print()
        return 1
    
    # Warnings
    if warnings:
        print()
        print(color_text("WARNINGS (recommended to fix for full functionality):", YELLOW))
        for msg in warnings:
            print(f"  - {msg}")
        print()
        print(color_text("Remediation:", YELLOW))
        print("  Install GitHub-based tools:")
        print("    pip install -U \\")
        print("      'git+https://github.com/vihaankulkarni29/federated_genome_extractor.git' \\")
        print("      'git+https://github.com/vihaankulkarni29/FastaAAExtractor.git' \\")
        print("      'git+https://github.com/vihaankulkarni29/wildtype-aligner.git' \\")
        print("      'git+https://github.com/vihaankulkarni29/SubScan.git'")
        print()
        if plat['is_windows'] and not any('WSL' in w for w in warnings):
            print("  For ABRicate (required for Annotator):")
            print("    Windows: Install WSL Ubuntu, then: sudo apt-get install abricate && abricate --setupdb")
        print()
    
    if not critical_failures:
        print(color_text("✓ CORE DEPENDENCIES MET - Pipeline ready for development!", GREEN))
        if warnings:
            print(color_text("  (Address warnings above for full production capability)", YELLOW))
        print()
        return 0
    
    return 1

if __name__ == "__main__":
    sys.exit(main())
