# Contributing to MutationScan

Thank you for your interest in contributing to MutationScan! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

We are committed to providing a welcoming and inclusive environment for all contributors. Please be respectful and constructive in all interactions.

## How to Contribute

### 1. Reporting Bugs

Found a bug? Please open an issue with:
- **Clear title** describing the problem
- **Detailed description** of the bug
- **Steps to reproduce** the issue
- **Expected vs. actual behavior**
- **Environment details** (OS, Python version, etc.)
- **Relevant code snippets** or error messages

### 2. Suggesting Features

Have an idea? Please open an issue with:
- **Clear title** for the feature
- **Detailed description** of what and why
- **Use cases** and benefits
- **Potential implementation approach** (optional)

### 3. Code Contributions

#### Setup Development Environment

```bash
# Clone repository
git clone https://github.com/yourusername/MutationScan.git
cd MutationScan

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install development dependencies
pip install -e ".[dev]"

# Pre-commit hooks
pip install pre-commit
pre-commit install
```

#### Development Workflow

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make changes** following code standards (see below)

3. **Write/update tests**
   - Unit tests in `tests/unit/`
   - Integration tests in `tests/integration/`
   - Run tests: `pytest -v`

4. **Run quality checks**
   ```bash
   # Format code
   black src/ tests/

   # Sort imports
   isort src/ tests/

   # Lint code
   flake8 src/ tests/

   # Type checking
   mypy src/

   # Coverage report
   pytest --cov=src --cov-report=term-missing tests/
   ```

5. **Commit changes** with clear messages
   ```bash
   git commit -m "feat: add genome parallel processing"
   ```

6. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

7. **Open a Pull Request** with description of changes

## Code Standards

### Style Guide

- **Language**: Python 3.8+
- **Formatter**: Black (line length: 100)
- **Linter**: Flake8
- **Sorting**: isort
- **Type hints**: Required for public APIs

### Naming Conventions

```python
# Classes: PascalCase
class GenomeDownloader:
    pass

# Functions/methods: snake_case
def download_genome():
    pass

# Constants: UPPER_SNAKE_CASE
MAX_RETRIES = 3

# Private members: _leading_underscore
_internal_function = None
```

### Docstring Format

Use Google-style docstrings:

```python
def example_function(param1: str, param2: int = 10) -> bool:
    """
    Brief description of what the function does.

    Longer description explaining the purpose, behavior, and any important
    considerations when using this function.

    Args:
        param1: Description of param1
        param2: Description of param2 (default: 10)

    Returns:
        Description of return value

    Raises:
        ValueError: Description of when this is raised
        TypeError: Description of when this is raised

    Example:
        >>> result = example_function("test", param2=20)
        >>> print(result)
        True
    """
    return True
```

### Type Hints

Always include type hints:

```python
from typing import Optional, List, Dict, Tuple

def process_sequences(
    sequences: List[str],
    config: Dict[str, any],
    output_file: Optional[str] = None
) -> Tuple[List[str], Dict[str, int]]:
    """Process sequences with optional output."""
    pass
```

### Logging

Use structured logging:

```python
import logging

logger = logging.getLogger(__name__)

logger.debug("Detailed information for debugging")
logger.info("General informational message")
logger.warning("Warning about potential issues")
logger.error("Error occurred but continuing")
logger.critical("Critical error, may need to stop")
```

### Error Handling

Handle errors gracefully:

```python
try:
    result = download_genome(accession)
except FileNotFoundError:
    logger.error(f"Genome {accession} not found")
    raise
except Exception as e:
    logger.error(f"Unexpected error downloading genome: {e}")
    raise
```

## Testing Guidelines

### Writing Tests

```python
import pytest
from mutation_scan.sequence_extractor import SequenceTranslator

class TestSequenceTranslator:
    """Test suite for SequenceTranslator class."""

    def setup_method(self):
        """Setup test fixtures."""
        self.translator = SequenceTranslator()

    def test_translate_standard_sequence(self):
        """Test translation of standard DNA sequence."""
        dna = "ATGAAAGCG"
        protein = self.translator.translate(dna)
        assert protein == "MAA"

    def test_translate_with_frame(self):
        """Test translation with different reading frame."""
        dna = "ATGAAAGCG"
        protein = self.translator.translate(dna, frame=1)
        assert len(protein) > 0

    @pytest.mark.parametrize("frame", [0, 1, 2])
    def test_translate_all_frames(self, frame):
        """Test translation in all reading frames."""
        dna = "ATGAAAGCGTAA"
        protein = self.translator.translate(dna, frame=frame)
        assert isinstance(protein, str)
```

### Test Organization

```
tests/
├── unit/                          # Fast, isolated tests
│   ├── test_genome_extractor.py
│   ├── test_gene_finder.py
│   └── test_sequence_extractor.py
├── integration/                   # Slower, end-to-end tests
│   ├── test_pipeline_flow.py
│   └── test_api_integration.py
└── fixtures/                      # Test data
    ├── sample_genomes/
    └── sample_results/
```

### Coverage Requirements

- Aim for >80% code coverage
- All public APIs must be tested
- Test both success and error cases
- Run coverage report: `pytest --cov=src tests/`

## Documentation

### Docstring Requirements

- All modules must have module-level docstrings
- All classes must have class-level docstrings
- All public methods must have docstrings
- Private methods should have docstrings

### Documentation Updates

When adding features:
1. Update relevant docstrings
2. Update README.md if needed
3. Add examples in docstrings
4. Update CHANGELOG.md

## Commit Message Format

Follow conventional commits:

```
type(scope): subject

body

footer
```

**Types**: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

**Examples**:
```
feat(genome_extractor): add retry logic for NCBI API

fix(visualizer): correct PyMOL image generation resolution

docs(README): add advanced usage examples

refactor(variant_caller): optimize alignment algorithm
```

## Pull Request Process

1. **Update documentation** if needed
2. **Ensure tests pass**: `pytest -v`
3. **Ensure code quality**: `black`, `flake8`, `mypy`
4. **Add descriptive PR description** with:
   - What changes were made
   - Why these changes were needed
   - Related issues/PRs
   - Testing performed
5. **Request reviewers** from core team
6. **Address feedback** in subsequent commits
7. **Squash commits** before merge if requested

## Release Process

1. Update version in `setup.py` and `__init__.py`
2. Update CHANGELOG.md
3. Tag release: `git tag -a v0.1.0 -m "Release version 0.1.0"`
4. Push tag: `git push origin --tags`

## Additional Resources

- [Python Style Guide (PEP 8)](https://pep8.org/)
- [Google Docstring Style](https://google.github.io/styleguide/pyguide.html)
- [Conventional Commits](https://www.conventionalcommits.org/)
- [Pytest Documentation](https://docs.pytest.org/)

## Questions?

- Open an issue for questions or discussions
- Check existing issues for similar questions
- Join project discussions on GitHub

Thank you for contributing to MutationScan! 🙏
