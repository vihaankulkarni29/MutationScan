# MutationScan 🧬

**A Comprehensive Bioinformatics Pipeline for Antimicrobial Resistance Analysis**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Pipeline Status](https://img.shields.io/badge/pipeline-production%20ready-green.svg)](https://github.com/vihaankulkarni29/MutationScan)

MutationScan is a state-of-the-art bioinformatics pipeline designed for comprehensive antimicrobial resistance (AMR) analysis. It provides end-to-end functionality from genome acquisition to interactive visualization, making it an essential tool for researchers studying bacterial resistance mechanisms.

## 🎯 Overview

MutationScan implements a sophisticated "Domino Effect" pipeline architecture, where each stage seamlessly feeds into the next, ensuring comprehensive and reproducible AMR analysis:

```
🧬 Genome Acquisition → 🔍 Gene Annotation → 🧪 Sequence Extraction → 
📏 Wild-type Alignment → 🔬 Mutation Analysis → 📊 Co-occurrence Analysis → 📄 Interactive Reports
```

## ✨ Key Features

### 🚀 **Complete Pipeline Integration**
- **Domino 1**: NCBI Genome Harvester with intelligent quality scoring
- **Domino 2**: ABRicate Automator for CARD database integration
- **Domino 3**: FastaAA Extractor for targeted sequence extraction
- **Domino 4**: Wild-type Aligner using EMBOSS alignment tools
- **Domino 5**: Mutation Analyzer with sophisticated detection algorithms
- **Domino 6**: Co-occurrence Analysis for resistance pattern discovery
- **Domino 7**: HTML Report Generator with interactive visualizations

### 🔬 **Advanced Analytics**
- **Smart Mutation Detection**: Identifies SNPs, insertions, deletions, and complex variants
- **Co-occurrence Pattern Analysis**: Discovers resistance gene relationships
- **Quality Scoring**: Automated genome quality assessment and ranking
- **Statistical Analysis**: Comprehensive mutation frequency and distribution analysis

### 📊 **Interactive Visualizations**
- **Plotly.js Charts**: Interactive mutation distributions, gene family analysis, and network visualizations
- **MSAViewer Integration**: In-browser multiple sequence alignment viewer
- **Responsive Dashboard**: Professional HTML reports with real-time data exploration
- **Publication-Ready Outputs**: High-quality figures and comprehensive data tables

## 🛠️ Installation

### Prerequisites
- Python 3.8+
- Git
- EMBOSS Suite (for sequence alignment)
- Internet connection (for NCBI API access)

### Quick Start

```bash
# Clone the repository
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan

# Install dependencies
pip install -r requirements.txt

# Install EMBOSS (Ubuntu/Debian)
sudo apt-get install emboss

# Install EMBOSS (macOS with Homebrew)
brew install brewsci/bio/emboss

# Install EMBOSS (Windows)
# Download from: http://emboss.open-bio.org/
```

## 🚀 Quick Start Guide

### 1. Basic Pipeline Run

```bash
# Step 1: Harvest genomes from NCBI
python subscan/tools/run_harvester.py \
  --accessions "NZ_CP107554,NZ_CP107555,NZ_CP107556" \
  --output-dir harvester_results

# Step 2: Annotate resistance genes
python subscan/tools/run_annotator.py \
  --manifest harvester_results/genome_manifest.json \
  --output-dir annotator_results

# Step 3: Extract target sequences
python subscan/tools/run_extractor.py \
  --manifest annotator_results/annotation_manifest.json \
  --output-dir extractor_results

# Step 4: Align to wild-type sequences
python subscan/tools/run_aligner.py \
  --manifest extractor_results/extraction_manifest.json \
  --output-dir aligner_results

# Step 5: Analyze mutations
python subscan/tools/run_analyzer.py \
  --manifest aligner_results/alignment_manifest.json \
  --output-dir analyzer_results

# Step 6: Generate interactive report
python subscan/tools/run_reporter.py \
  --manifest analyzer_results/analysis_manifest.json \
  --output-dir reports \
  --open-browser
```

### 2. Automated End-to-End Analysis

```bash
# Run complete pipeline with a single command
python subscan/pipeline.py \
  --accessions-file research_accessions.txt \
  --output-dir complete_analysis \
  --threads 8 \
  --generate-report
```

## 📁 Project Structure

```
MutationScan/
├── subscan/
│   ├── src/subscan/           # Core pipeline modules
│   ├── tools/                 # Command-line interfaces for each domino
│   ├── analyzer/              # Advanced mutation analysis engine
│   ├── ncbi_genome_extractor/ # NCBI integration tools
│   ├── tests/                 # Comprehensive test suite
│   └── pyproject.toml         # Package configuration
├── docs/                      # Documentation and tutorials
├── examples/                  # Example datasets and workflows
├── .gitignore                 # Git ignore rules
└── README.md                  # This file
```

## 🔧 Advanced Usage

### Custom Gene Lists

Create custom gene lists for targeted analysis:

```python
# custom_genes.txt
gyrA
parC
rpoB
katG
inhA
```

```bash
python subscan/tools/run_extractor.py \
  --gene-list custom_genes.txt \
  --manifest annotation_manifest.json
```

### Parallel Processing

Leverage multi-core processing for large datasets:

```bash
python subscan/tools/run_analyzer.py \
  --threads 16 \
  --batch-size 100 \
  --manifest alignment_manifest.json
```

### Report Customization

Generate customized reports with specific visualizations:

```bash
python subscan/tools/run_reporter.py \
  --template custom_template.html \
  --include-msa \
  --export-data \
  --manifest analysis_manifest.json
```

## 📊 Output Examples

### Interactive Dashboard
![Dashboard Preview](docs/images/dashboard_preview.png)

### Mutation Analysis Report
- **Total Genomes Analyzed**: 150
- **Mutations Detected**: 5,284
- **Co-occurrence Pairs**: 171,379
- **Gene Families**: 12

### Visualization Components
- 📈 Interactive mutation distribution charts
- 🔗 Gene co-occurrence network diagrams
- 📋 Comprehensive data tables with filtering
- 🧬 Multiple sequence alignment viewer

## 🧪 Testing

Run the comprehensive test suite:

```bash
# Run all tests
python -m pytest subscan/tests/

# Run specific test modules
python -m pytest subscan/analyzer/tests/test_engine.py

# Run with coverage report
python -m pytest --cov=subscan subscan/tests/
```

## 📚 Documentation

- **[User Guide](docs/user_guide.md)**: Comprehensive usage instructions
- **[API Reference](docs/api_reference.md)**: Detailed module documentation
- **[Tutorial](docs/tutorial.md)**: Step-by-step walkthrough with example data
- **[Contributing](docs/contributing.md)**: Guidelines for contributors

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone for development
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan

# Install in development mode
pip install -e .

# Install development dependencies
pip install -r requirements-dev.txt

# Run pre-commit hooks
pre-commit install
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **NCBI**: For providing comprehensive genomic databases
- **CARD Database**: For curated antimicrobial resistance gene data
- **EMBOSS Suite**: For robust sequence alignment tools
- **Plotly.js**: For interactive visualization capabilities
- **MSAViewer**: For multiple sequence alignment visualization

## 📧 Contact

- **Author**: Vihaan Kulkarni
- **Email**: [vihaankulkarni29@gmail.com](mailto:vihaankulkarni29@gmail.com)
- **GitHub**: [@vihaankulkarni29](https://github.com/vihaankulkarni29)

## 🚀 Citation

If you use MutationScan in your research, please cite:

```bibtex
@software{kulkarni2025mutationscan,
  title={MutationScan: A Comprehensive Bioinformatics Pipeline for Antimicrobial Resistance Analysis},
  author={Kulkarni, Vihaan},
  year={2025},
  url={https://github.com/vihaankulkarni29/MutationScan},
  version={1.0.0}
}
```

---

<div align="center">

**🧬 Advancing AMR Research Through Computational Excellence 🧬**

[Report Issues](https://github.com/vihaankulkarni29/MutationScan/issues) • [Request Features](https://github.com/vihaankulkarni29/MutationScan/issues) • [Join Discussions](https://github.com/vihaankulkarni29/MutationScan/discussions)

</div>