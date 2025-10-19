# Git Dependency Troubleshooting Guide

## Overview

MutationScan uses several Git-based Python packages as dependencies. This guide helps you diagnose and fix common installation issues.

---

## 🔍 Quick Diagnosis

Run the automated validator:
```bash
cd subscan
python validate_git_dependencies.py --verbose
```

This will check all Git dependencies for:
- ✅ Repository accessibility
- ✅ Required packaging files (`pyproject.toml`, `setup.py`, `__init__.py`)
- ✅ Pip installability

---

## 📦 Git Dependencies

MutationScan depends on these Git-based packages:

| Package | Repository | Purpose |
|---------|-----------|---------|
| `federated-genome-extractor` | [federated_genome_extractor](https://github.com/vihaankulkarni29/federated_genome_extractor) | Multi-database genome download |
| `abricate-automator` | [ABRicate-Automator](https://github.com/vihaankulkarni29/ABRicate-Automator) | AMR gene annotation |
| `fasta-aa-extractor` | [FastaAAExtractor](https://github.com/vihaankulkarni29/FastaAAExtractor) | Protein sequence extraction |
| `wildtype-aligner` | [wildtype-aligner](https://github.com/vihaankulkarni29/wildtype-aligner) | Sequence alignment |
| `subscan-analyzer` | [SubScan](https://github.com/vihaankulkarni29/SubScan) | Mutation analysis engine |

---

## ❌ Common Errors

### Error: "does not appear to be a Python project"

**Full Error:**
```
ERROR: <package> @ git+https://github.com/... does not appear to be a Python project: 
neither 'setup.py' nor 'pyproject.toml' found.
```

**Cause:** The repository is missing packaging files.

**Fix:**
1. Verify the repository has these files in the root:
   - `pyproject.toml` (modern, recommended)
   - `setup.py` (legacy fallback)
   - `__init__.py` (makes it a package)

2. Test the specific dependency:
   ```bash
   pip install git+https://github.com/vihaankulkarni29/<repo>.git
   ```

3. If it fails, the repository needs packaging files added:
   ```bash
   # Clone the repo
   git clone https://github.com/vihaankulkarni29/<repo>.git
   cd <repo>
   
   # Check for packaging files
   ls -la | grep -E "pyproject|setup.py|__init__.py"
   ```

---

### Error: Circular dependency

**Cause:** A dependency is trying to install itself (e.g., SubScan installing `subscan-analyzer @ git+SubScan.git`).

**Fix:**
1. Check the dependency's `pyproject.toml` or `setup.cfg`
2. Remove self-referencing dependencies
3. Ensure package names are distinct (e.g., `subscan` vs `subscan-analyzer`)

---

### Error: Git clone failed

**Full Error:**
```
ERROR: Error [Errno 2] No such file or directory: 'git'
```

**Cause:** Git is not installed or not in PATH.

**Fix:**

**Ubuntu/Debian:**
```bash
sudo apt-get install git
```

**Windows:**
- Download from https://git-scm.com/download/win
- Ensure "Git from the command line" is selected during installation
- Restart terminal after installation

**Verify:**
```bash
git --version
```

---

### Error: SSL certificate verification failed

**Full Error:**
```
fatal: unable to access 'https://github.com/...': SSL certificate problem
```

**Temporary Fix (not recommended for production):**
```bash
git config --global http.sslVerify false
pip install ...
git config --global http.sslVerify true
```

**Proper Fix:**
```bash
# Update certificates (Ubuntu/Debian)
sudo apt-get install ca-certificates
sudo update-ca-certificates

# Windows: Reinstall Git with proper certificate bundle
```

---

### Error: Permission denied (publickey)

**Cause:** Using SSH URL instead of HTTPS.

**Fix:** Git dependencies must use HTTPS URLs:
```python
# ❌ Wrong (SSH)
"package @ git+git@github.com:user/repo.git"

# ✅ Correct (HTTPS)
"package @ git+https://github.com/user/repo.git"
```

---

## 🔧 Manual Dependency Testing

Test each dependency individually:

```bash
# Test federated-genome-extractor
pip install git+https://github.com/vihaankulkarni29/federated_genome_extractor.git
python -c "import federated_genome_extractor; print('OK')"

# Test abricate-automator
pip install git+https://github.com/vihaankulkarni29/ABRicate-Automator.git
python -c "import abricate_automator; print('OK')"

# Test fasta-aa-extractor
pip install git+https://github.com/vihaankulkarni29/FastaAAExtractor.git
python -c "import fasta_aa_extractor; print('OK')"

# Test wildtype-aligner
pip install git+https://github.com/vihaankulkarni29/wildtype-aligner.git
python -c "import wildtype_aligner; print('OK')"

# Test subscan-analyzer
pip install git+https://github.com/vihaankulkarni29/SubScan.git
python -c "import subscan_analyzer; print('OK')"
```

---

## 🛠️ Fixing a Broken Dependency

If a dependency repository is missing packaging files:

### 1. Clone the repository
```bash
git clone https://github.com/vihaankulkarni29/<repo>.git
cd <repo>
```

### 2. Create `pyproject.toml`
```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "package-name"
version = "1.0.0"
description = "Package description"
requires-python = ">=3.8"
license = "MIT"
authors = [{name = "Vihaan Kulkarni", email = "vihaankulkarni29@gmail.com"}]
dependencies = [
    "biopython>=1.79,<2.0",
]

[project.urls]
Homepage = "https://github.com/vihaankulkarni29/<repo>"
Repository = "https://github.com/vihaankulkarni29/<repo>.git"

[tool.setuptools.packages.find]
where = ["."]
include = ["*"]
exclude = ["tests*", "docs*"]
```

### 3. Create `setup.py` (fallback)
```python
from setuptools import setup, find_packages

setup(
    name="package-name",
    version="1.0.0",
    description="Package description",
    author="Vihaan Kulkarni",
    author_email="vihaankulkarni29@gmail.com",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=["biopython>=1.79,<2.0"],
)
```

### 4. Create `__init__.py`
```python
"""Package description"""
__version__ = "1.0.0"
__author__ = "Vihaan Kulkarni"
```

### 5. Test locally
```bash
pip install -e .
```

### 6. Commit and push
```bash
git add pyproject.toml setup.py __init__.py
git commit -m "Add pip packaging support"
git push origin main
```

### 7. Verify from GitHub
```bash
pip install git+https://github.com/vihaankulkarni29/<repo>.git
```

---

## 🤖 Automated CI Validation

MutationScan includes GitHub Actions CI that automatically validates all Git dependencies:

**Workflow:** `.github/workflows/validate_dependencies.yml`

**What it checks:**
- Repository accessibility
- Presence of packaging files
- Pip install success
- Import test
- Cross-platform compatibility (Ubuntu + Windows)
- Python version compatibility (3.9-3.12)

**Triggers:**
- Every push to `main` or `develop`
- Every pull request
- Daily at 2 AM UTC (catches upstream issues)
- Manual workflow dispatch

**On failure:**
- Creates GitHub issue with diagnostic information
- Labels: `automated`, `dependencies`, `ci-failure`, `high-priority`

---

## 📝 Dependency Checklist

Before releasing a dependency update:

- [ ] `pyproject.toml` exists in repository root
- [ ] `setup.py` exists (fallback compatibility)
- [ ] `__init__.py` exists (makes it a package)
- [ ] No circular dependencies in `dependencies` list
- [ ] Test local install: `pip install -e .`
- [ ] Test GitHub install: `pip install git+<url>`
- [ ] Test import: `python -c "import <package>"`
- [ ] Run validator: `python validate_git_dependencies.py -d <package>`
- [ ] Version bumped if making changes
- [ ] Changes committed and pushed to GitHub

---

## 🚀 Full Installation Verification

After fixing dependencies, verify the complete installation:

```bash
# 1. Validate all Git dependencies
cd subscan
python validate_git_dependencies.py --verbose

# 2. Install MutationScan
pip install -e .

# 3. Run structure validation
python tools/check_structure.py

# 4. Run dependency check
python tools/check_dependencies.py

# 5. Run smoke tests
python tools/smoke_test_dominos.py

# 6. Quick validation
python tools/quick_validation.py
```

---

## 📞 Getting Help

If you encounter persistent dependency issues:

1. **Check GitHub Issues:** https://github.com/vihaankulkarni29/MutationScan/issues
2. **Run diagnostic:** `python validate_git_dependencies.py --verbose`
3. **Check CI status:** View latest workflow run
4. **Contact:** vihaankulkarni29@gmail.com

---

## 🔄 Dependency Update Workflow

When updating a dependency:

1. **Update dependency repository:**
   ```bash
   cd /path/to/dependency-repo
   # Make changes
   git add .
   git commit -m "Update: <description>"
   git push
   ```

2. **Verify packaging:**
   ```bash
   pip install git+https://github.com/vihaankulkarni29/<repo>.git
   ```

3. **Update MutationScan (if needed):**
   ```bash
   cd /path/to/MutationScan/subscan
   # Version constraints in setup.cfg or pyproject.toml
   ```

4. **Test integration:**
   ```bash
   python validate_git_dependencies.py
   pip install -e .
   python tools/quick_validation.py
   ```

5. **Commit and push MutationScan:**
   ```bash
   git add .
   git commit -m "Update dependency: <package>"
   git push
   ```

---

*Last updated: 2025-10-19*