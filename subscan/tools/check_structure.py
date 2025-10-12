#!/usr/bin/env python3
"""
Pipeline Structure Validator

Validates that the MutationScan repository has all required files, directories,
and entry points for the 7-domino pipeline to function correctly.

Run this from the repository root:
    python subscan/tools/check_structure.py

Exit codes:
    0 - All checks passed
    1 - One or more checks failed
"""
import os
import sys
from pathlib import Path
from typing import List, Tuple

# ANSI color codes for terminal output (safe fallback for non-color terminals)
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

def check_exists(path: Path, description: str) -> Tuple[bool, str]:
    """Check if a path exists and return (passed, message)."""
    if path.exists():
        return True, f"[OK] {description}: {path}"
    else:
        return False, f"[FAIL] {description} missing: {path}"

def main() -> int:
    """Run all structure validation checks."""
    print(color_text("\n" + "=" * 70, BLUE))
    print(color_text("MutationScan Pipeline Structure Validator", BLUE))
    print(color_text("=" * 70 + "\n", BLUE))
    
    # Detect repository root
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]  # tools -> subscan -> MutationScan
    
    print(f"Repository root: {repo_root}\n")
    
    if not (repo_root / "subscan").exists():
        print(color_text("[ERROR] Cannot find 'subscan' directory. Please run from repo root.", RED))
        return 1
    
    checks: List[Tuple[bool, str]] = []
    
    # ========== Core Repository Files ==========
    print(color_text("Checking core repository files...", BLUE))
    
    core_files = [
        (repo_root / "README.md", "Main README"),
        (repo_root / "LICENSE", "Root LICENSE"),
        (repo_root / "requirements.txt", "Root requirements"),
        (repo_root / "CONTRIBUTING.md", "Contributing guide"),
        (repo_root / "SECURITY.md", "Security policy"),
    ]
    
    for path, desc in core_files:
        passed, msg = check_exists(path, desc)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else RED))
    
    print()
    
    # ========== SubScan Package Structure ==========
    print(color_text("Checking subscan package structure...", BLUE))
    
    subscan_root = repo_root / "subscan"
    subscan_files = [
        (subscan_root / "pyproject.toml", "SubScan pyproject.toml"),
        (subscan_root / "LICENSE", "SubScan LICENSE"),
        (subscan_root / "README.md", "SubScan README"),
        (subscan_root / "setup.cfg", "SubScan setup.cfg"),
    ]
    
    for path, desc in subscan_files:
        passed, msg = check_exists(path, desc)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else RED))
    
    print()
    
    # ========== Domino Tool Scripts ==========
    print(color_text("Checking domino tool scripts...", BLUE))
    
    tools_dir = subscan_root / "tools"
    domino_scripts = [
        ("run_harvester.py", "Domino 1: Harvester"),
        ("run_annotator.py", "Domino 2: Annotator"),
        ("run_extractor.py", "Domino 3: Extractor"),
        ("run_aligner.py", "Domino 4: Aligner"),
        ("run_analyzer.py", "Domino 5: Analyzer"),
        ("run_cooccurrence_analyzer.py", "Domino 6: CoOccurrence"),
        ("run_reporter.py", "Domino 7: Reporter"),
        ("run_pipeline.py", "Orchestrator"),
    ]
    
    for script, desc in domino_scripts:
        path = tools_dir / script
        passed, msg = check_exists(path, desc)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else RED))
    
    print()
    
    # ========== Support Scripts ==========
    print(color_text("Checking support scripts...", BLUE))
    
    support_scripts = [
        (tools_dir / "smoke_test_extractor.py", "Extractor smoke test"),
        (tools_dir / "mock_fasta_aa_extractor.py", "Mock extractor (for testing)"),
    ]
    
    for path, desc in support_scripts:
        passed, msg = check_exists(path, desc)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else YELLOW))
    
    print()
    
    # ========== Data Directories ==========
    print(color_text("Checking data directories...", BLUE))
    
    data_dirs = [
        (repo_root / "data_input", "Input data directory"),
        (repo_root / "data_output", "Output data directory"),
        (subscan_root / "sample_data", "Sample data directory"),
        (repo_root / "examples", "Examples directory"),
    ]
    
    for path, desc in data_dirs:
        passed, msg = check_exists(path, desc)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else YELLOW))
    
    print()
    
    # ========== Documentation ==========
    print(color_text("Checking documentation...", BLUE))
    
    docs_files = [
        (repo_root / "docs" / "README.md", "Technical documentation"),
        (repo_root / "docs" / "TEST_PLAN.md", "Test plan"),
    ]
    
    for path, desc in docs_files:
        passed, msg = check_exists(path, desc)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else YELLOW))
    
    print()
    
    # ========== Sample Data Files ==========
    print(color_text("Checking sample data files...", BLUE))
    
    sample_files = [
        (subscan_root / "sample_data" / "sample_accessions.txt", "Sample accessions"),
        (subscan_root / "sample_data" / "gene_list.txt", "Sample gene list"),
    ]
    
    for path, desc in sample_files:
        passed, msg = check_exists(path, desc)
        checks.append((passed, msg))
        print(color_text(msg, GREEN if passed else YELLOW))
    
    print()
    
    # ========== Entry Point Validation ==========
    print(color_text("Checking entry points in pyproject.toml...", BLUE))
    
    pyproject = subscan_root / "pyproject.toml"
    if pyproject.exists():
        try:
            with open(pyproject, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Check for mutationscan entry point
            if 'mutationscan =' in content:
                print(color_text("[OK] Entry point 'mutationscan' defined", GREEN))
                checks.append((True, "Entry point 'mutationscan' defined"))
            else:
                print(color_text("[FAIL] Entry point 'mutationscan' not found in pyproject.toml", RED))
                checks.append((False, "Entry point 'mutationscan' not found"))
            
            # Check for project.scripts section
            if '[project.scripts]' in content:
                print(color_text("[OK] [project.scripts] section present", GREEN))
                checks.append((True, "[project.scripts] section present"))
            else:
                print(color_text("[FAIL] [project.scripts] section missing", RED))
                checks.append((False, "[project.scripts] section missing"))
                
        except Exception as e:
            print(color_text(f"[WARN] Could not parse pyproject.toml: {e}", YELLOW))
    
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
        print(color_text(f"Failed: {failed_checks}", RED))
        print()
        print(color_text("FAILED CHECKS:", RED))
        for passed, msg in checks:
            if not passed:
                print(f"  - {msg}")
        
        print()
        print(color_text("Remediation:", YELLOW))
        print("  1. Ensure you're on the latest 'main' branch:")
        print("     git fetch origin && git switch main && git pull")
        print()
        print("  2. If files are missing after pulling, check for:")
        print("     - Nested clone directories (e.g., MutationScan/MutationScan)")
        print("     - Incorrect working directory")
        print()
        print("  3. For missing sample data, create placeholder files:")
        print("     mkdir -p subscan/sample_data")
        print("     echo 'GCF_000005825.2' > subscan/sample_data/sample_accessions.txt")
        print("     echo 'mecA' > subscan/sample_data/gene_list.txt")
        print()
        return 1
    else:
        print(color_text("\n✓ ALL CHECKS PASSED - Repository structure is valid!\n", GREEN))
        return 0

if __name__ == "__main__":
    sys.exit(main())
