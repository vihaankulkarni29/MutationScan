# Contributing to MutationScan 🧬

Thank you for your interest in contributing to MutationScan! This document provides guidelines for contributing to our antimicrobial resistance analysis pipeline.

## 🎯 How to Contribute

### 🐛 Reporting Bugs

Before creating bug reports, please check the [issue tracker](https://github.com/vihaankulkarni29/MutationScan/issues) to see if the problem has already been reported.

When creating a bug report, please include:

- **Clear description** of the issue
- **Steps to reproduce** the behavior
- **Expected behavior** vs. actual behavior
- **Environment details** (OS, Python version, dependencies)
- **Log files** or error messages if applicable
- **Sample data** if the issue is data-specific

### 💡 Suggesting Features

We welcome feature requests! Please:

1. Check existing [feature requests](https://github.com/vihaankulkarni29/MutationScan/issues?q=is%3Aissue+is%3Aopen+label%3Aenhancement)
2. Create a detailed issue describing:
   - The problem you're trying to solve
   - Your proposed solution
   - Alternative solutions considered
   - Additional context or examples

### 🔧 Development Workflow

#### 1. Fork and Clone

```bash
# Fork the repository on GitHub, then clone your fork
git clone https://github.com/YOUR-USERNAME/MutationScan.git
cd MutationScan

# Add upstream remote
git remote add upstream https://github.com/vihaankulkarni29/MutationScan.git
```

#### 2. Set Up Development Environment

```bash
# Create a virtual environment
python -m venv mutationscan-dev
source mutationscan-dev/bin/activate  # On Windows: mutationscan-dev\Scripts\activate

# Install in development mode
pip install -e .

# Install development dependencies
pip install -r requirements-dev.txt

# Install pre-commit hooks
pre-commit install
```

#### 3. Create a Feature Branch

```bash
# Update your fork
git checkout main
git pull upstream main

# Create feature branch
git checkout -b feature/your-feature-name
```

#### 4. Make Your Changes

- Follow the [coding standards](#coding-standards)
- Write tests for new functionality
- Update documentation as needed
- Ensure all tests pass

#### 5. Test Your Changes

```bash
# Run the full test suite
python -m pytest subscan/tests/

# Run specific test modules
python -m pytest subscan/tests/test_pipeline.py

# Check code coverage
python -m pytest --cov=subscan subscan/tests/

# Run linting
flake8 subscan/
black --check subscan/
mypy subscan/
```

#### 6. Submit a Pull Request

```bash
# Commit your changes
git add .
git commit -m "feat: add your feature description"

# Push to your fork
git push origin feature/your-feature-name
```

Then create a pull request on GitHub with:
- Clear title and description
- Reference to related issues
- Description of changes made
- Test coverage information

## 📝 Coding Standards

### Python Style Guide

We follow [PEP 8](https://pep8.org/) with these specific guidelines:

- **Line length**: 88 characters (Black formatter default)
- **Imports**: Use absolute imports, group by standard library, third-party, local
- **Docstrings**: Google-style docstrings for all public functions and classes
- **Type hints**: Required for all function signatures
- **Variable naming**: `snake_case` for variables and functions, `PascalCase` for classes

### Code Quality Tools

- **Formatter**: [Black](https://black.readthedocs.io/)
- **Linter**: [Flake8](https://flake8.pycqa.org/)
- **Type checker**: [MyPy](https://mypy.readthedocs.io/)
- **Import sorting**: [isort](https://pycqa.github.io/isort/)

### Example Code Style

```python
from typing import List, Dict, Optional
import os
import sys

from subscan.utils import load_json_manifest


def process_genomes(
    manifest_path: str, 
    output_dir: str, 
    threads: int = 4
) -> Dict[str, str]:
    """
    Process genomes from manifest file.
    
    Args:
        manifest_path: Path to the genome manifest JSON file
        output_dir: Directory for output files
        threads: Number of parallel threads to use
        
    Returns:
        Dictionary mapping genome IDs to output file paths
        
    Raises:
        FileNotFoundError: If manifest file doesn't exist
        ValueError: If manifest format is invalid
    """
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")
    
    manifest_data = load_json_manifest(manifest_path)
    # ... implementation
    
    return results
```

## 🧪 Testing Guidelines

### Test Structure

- **Unit tests**: Test individual functions in isolation
- **Integration tests**: Test pipeline components working together
- **End-to-end tests**: Test complete workflows with sample data

### Writing Tests

```python
import pytest
from pathlib import Path
from subscan.tools.run_harvester import parse_arguments


class TestHarvester:
    """Test suite for the genome harvester module."""
    
    def test_parse_arguments_valid_input(self):
        """Test argument parsing with valid inputs."""
        args = parse_arguments([
            "--accessions", "NZ_CP107554,NZ_CP107555",
            "--output-dir", "/tmp/test"
        ])
        assert args.accessions == "NZ_CP107554,NZ_CP107555"
        assert args.output_dir == "/tmp/test"
    
    @pytest.mark.parametrize("invalid_accession", [
        "",
        "invalid-format",
        "123456789"
    ])
    def test_parse_arguments_invalid_accessions(self, invalid_accession):
        """Test argument parsing with invalid accession formats."""
        with pytest.raises(SystemExit):
            parse_arguments([
                "--accessions", invalid_accession,
                "--output-dir", "/tmp/test"
            ])
```

### Test Data

- Place test data in `subscan/tests/data/`
- Use small, representative datasets
- Include both positive and negative test cases
- Mock external API calls when possible

## 📚 Documentation

### Docstring Format

Use Google-style docstrings:

```python
def analyze_mutations(alignment_file: str, reference_file: str) -> Dict[str, Any]:
    """
    Analyze mutations in aligned sequences.
    
    This function compares aligned sequences to reference sequences and identifies
    mutations including SNPs, insertions, and deletions.
    
    Args:
        alignment_file: Path to the alignment file in FASTA format
        reference_file: Path to the reference sequence file
        
    Returns:
        Dictionary containing mutation analysis results with keys:
            - 'mutations': List of detected mutations
            - 'statistics': Summary statistics
            - 'quality_score': Overall alignment quality
            
    Raises:
        FileNotFoundError: If input files don't exist
        ValueError: If file formats are invalid
        
    Example:
        >>> results = analyze_mutations("aligned.fasta", "reference.fasta")
        >>> print(f"Found {len(results['mutations'])} mutations")
    """
```

### README Updates

When adding new features, update relevant documentation:
- Main README.md
- Module-specific README files
- API documentation
- Usage examples

## 🏷️ Commit Message Guidelines

We use [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>[optional scope]: <description>

[optional body]

[optional footer(s)]
```

### Types

- **feat**: New feature
- **fix**: Bug fix
- **docs**: Documentation changes
- **style**: Code style changes (formatting, etc.)
- **refactor**: Code refactoring
- **test**: Adding or updating tests
- **chore**: Maintenance tasks

### Examples

```
feat(analyzer): add support for complex variant detection

- Implement detection of complex structural variants
- Add support for tandem repeats and inversions
- Update test suite with new variant types

Closes #123
```

```
fix(harvester): handle NCBI API rate limiting

- Add exponential backoff for API requests
- Implement retry logic with maximum attempts
- Log rate limiting events for debugging

Fixes #456
```

## 🌟 Recognition

Contributors will be recognized in:
- `CONTRIBUTORS.md` file
- Release notes for their contributions
- GitHub contributors graph

## 📞 Getting Help

- **Discussions**: [GitHub Discussions](https://github.com/vihaankulkarni29/MutationScan/discussions)
- **Issues**: [GitHub Issues](https://github.com/vihaankulkarni29/MutationScan/issues)
- **Email**: vihaankulkarni29@gmail.com

## 📋 Checklist for Pull Requests

Before submitting a pull request, ensure:

- [ ] Code follows the style guide
- [ ] All tests pass
- [ ] New functionality has tests
- [ ] Documentation is updated
- [ ] Commit messages follow conventions
- [ ] No merge conflicts with main branch
- [ ] PR description is clear and complete

Thank you for contributing to MutationScan! Your efforts help advance antimicrobial resistance research. 🧬