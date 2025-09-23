# MutationScan Repository Cleanup

This directory contains cleanup scripts to prepare the MutationScan repository for its first official public release.

## Files Created

### Phase 1: .gitignore
- **`.gitignore`** - Comprehensive gitignore file that excludes:
  - Python cache files (`__pycache__/`, `*.pyc`)
  - Output directories (`*_output/`, `*_results/`)
  - Temporary files (`*.tmp`, `debug_*.py`)
  - IDE configurations (`.vscode/`, `.idea/`)
  - System files (`.DS_Store`, `Thumbs.db`)

### Phase 2: Cleanup Scripts
- **`cleanup_repository.sh`** - Bash script for Linux/macOS
- **`cleanup_repository.ps1`** - PowerShell script for Windows

## How to Run Cleanup

### On Windows (PowerShell):
```powershell
# Navigate to repository root
cd C:\Users\Vihaan\Desktop\GenomeAnalyzer

# Run the PowerShell cleanup script
.\cleanup_repository.ps1
```

### On Linux/macOS (Bash):
```bash
# Navigate to repository root
cd /path/to/GenomeAnalyzer

# Make script executable
chmod +x cleanup_repository.sh

# Run the cleanup script
./cleanup_repository.sh
```

## What Gets Removed

✅ **Safe to Remove:**
- Python bytecode files (`__pycache__/`, `*.pyc`)
- Test output directories (`harvester_test_results/`, `annotator_test_results/`)
- Temporary files (`test_html_reporter.py`, `debug_paths.py`)
- Pipeline output directories (`harvester_output/`, `test_parallel_output/`)
- Log files and system files

❌ **Preserved:**
- Source code in `subscan/src/`
- Essential tools in `subscan/tools/`
- Test code in `subscan/tests/`, `subscan/analyzer/tests/`, etc.
- Documentation and configuration files
- Core pipeline components

## After Cleanup

The repository will be ready for:
1. **Public Release** on GitHub
2. **Clean Development** environment
3. **Professional Presentation** to users

The cleanup ensures only essential source code and documentation are tracked by Git, making the repository fast, clean, and professional.