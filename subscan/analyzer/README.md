# SubScan Analyzer

**SubScan Pipeline - Domino 5: The Analyzer**

A high-performance mutation detection and co-occurrence analysis engine for the SubScan bioinformatics pipeline.

## Overview

The SubScan Analyzer is the fifth and final analytical domino in the SubScan pipeline, designed to perform core scientific discovery by identifying and cataloging amino acid mutations relative to reference sequences. It also performs sophisticated co-occurrence analysis to find patterns of mutations across different genes within the same genome.

## Features

### 🧬 **Mutation Detection Engine**
- **Reference-driven analysis**: Compares aligned sequences to reference sequences
- **Comprehensive mutation types**: Detects substitutions, insertions, and deletions
- **Position-specific tracking**: Records exact amino acid positions and changes
- **Multi-genome analysis**: Processes mutations across entire genome collections

### 🔍 **Co-occurrence Analysis**
- **Pattern recognition**: Identifies mutations that appear together in genomes
- **Gene-level correlation**: Analyzes relationships between mutated genes
- **Statistical significance**: Provides quantitative co-occurrence metrics
- **Cross-genome patterns**: Discovers mutation patterns across sample collections

### 📊 **Scientific Output**
- **Detailed CSV reports**: Machine-readable mutation catalogs
- **Co-occurrence matrices**: Quantitative relationship analysis
- **Standardized manifest**: Pipeline-ready JSON output for reporting
- **Quality metrics**: Comprehensive analysis statistics

## Pipeline Integration

```
Domino 4 (Aligner) → alignment_manifest.json → Domino 5 (Analyzer) → analysis_manifest.json → Reporting Stage
```

### Input
- `alignment_manifest.json`: Contains paths to reference-aligned protein sequences
- Aligned FASTA files from WildTypeAligner processing

### Output
- `mutation_report.csv`: Comprehensive catalog of all detected mutations
- `cooccurrence_report.csv`: Analysis of mutation co-occurrence patterns
- `analysis_manifest.json`: Standardized manifest for downstream reporting

## Installation

```bash
# Development installation
pip install -e .

# Production installation
pip install subscan-analyzer

# With development dependencies
pip install -e ".[dev]"
```

## Usage

### Command Line Interface

```bash
# Basic usage
subscan-analyzer --manifest alignment_manifest.json --output-dir ./analysis_results

# Using the tools script directly
python tools/run_analyzer.py --manifest alignment_manifest.json --output-dir ./results
```

### Python API

```python
from analyzer.engine import MutationAnalyzer

# Initialize analyzer
analyzer = MutationAnalyzer()

# Analyze mutations in an alignment
mutations = analyzer.identify_mutations("path/to/alignment.fasta")

# Analyze co-occurrence patterns
cooccurrence = analyzer.analyze_cooccurrence(all_mutations)
```

## Architecture

### Core Components

- **`src/analyzer/engine.py`**: Core mutation analysis engine
- **`tools/run_analyzer.py`**: Command-line interface
- **`tests/`**: Comprehensive test suite

### Scientific Algorithms

1. **Reference Identification**: Automatically identifies reference sequences in alignments
2. **Position-wise Comparison**: Iterates through every position for mutation detection
3. **Mutation Classification**: Categorizes substitutions, insertions, and deletions
4. **Co-occurrence Matrix**: Builds quantitative relationship matrices
5. **Statistical Analysis**: Provides significance testing for patterns

## Development

### Setup Development Environment

```bash
# Clone repository
git clone https://github.com/subscan/analyzer.git
cd analyzer

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=analyzer --cov-report=html

# Run specific test file
pytest tests/test_engine.py
```

### Code Quality

```bash
# Format code
black src/ tools/ tests/

# Check code style
flake8 src/ tools/ tests/

# Type checking
mypy src/ tools/
```

## Scientific Background

### Mutation Analysis
The analyzer uses sophisticated bioinformatics algorithms to identify amino acid changes that may affect protein function, particularly in antibiotic resistance contexts.

### Co-occurrence Analysis
Statistical analysis of mutation patterns helps identify:
- Compensatory mutations that occur together
- Resistance mechanisms involving multiple genes
- Evolutionary relationships between mutations

## License

MIT License - see LICENSE file for details.

## Citation

If you use SubScan Analyzer in your research, please cite:

```
SubScan: A High-Performance Bioinformatics Pipeline for Antibiotic Resistance Analysis
[Citation details to be added]
```

## Support

- **Documentation**: [https://subscan-analyzer.readthedocs.io/](https://subscan-analyzer.readthedocs.io/)
- **Issues**: [https://github.com/subscan/analyzer/issues](https://github.com/subscan/analyzer/issues)
- **Discussions**: [https://github.com/subscan/analyzer/discussions](https://github.com/subscan/analyzer/discussions)