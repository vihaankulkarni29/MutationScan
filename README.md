# MutationScan

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python: 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Status: Alpha](https://img.shields.io/badge/Status-Alpha-orange.svg)

A professional, modular bioinformatics pipeline for automated detection, analysis, and visualization of antimicrobial resistance (AMR) genes in bacterial genomes.

## 📚 Documentation

All documentation is organized in the [docs/](docs/) folder:

- **[Contributing Guide](docs/CONTRIBUTING.md)** - How to contribute to the project
- **[Project Structure](docs/PROJECT_STRUCTURE.md)** - Detailed directory architecture
- **[Quick Reference](docs/QUICK_REFERENCE.md)** - Quick lookup guide
- **[Documentation Index](docs/INDEX.md)** - Complete documentation index
- **[Setup Details](docs/SETUP_COMPLETE.md)** - Installation and setup verification

## Overview

MutationScan streamlines the complex process of genomic analysis by automating:

- **Genome Acquisition**: Automated downloading of bacterial genomes from NCBI
- **AMR Screening**: Detection of antimicrobial resistance genes using ABRicate and BLASTn
- **Sequence Processing**: Coordinate-based protein extraction and translation
- **Variant Analysis**: Pairwise alignment and mutation calling
- **3D Visualization**: PyMOL-based mutation mapping and structure visualization

## Research Runs (User-Provided Inputs)

For scientific runs, MutationScan expects the user to provide the following inputs up front:

- Genome data: either a local genome FASTA file or an organism query you want analyzed.
- Target genes: a focused list of genes to analyze (for example, acrR).
- NCBI credentials: email required by NCBI policy and an optional API key for higher throughput.

The preferred execution path is a single Docker (or WSL) command. See the setup guide at [docs/SETUP_COMPLETE.md](docs/SETUP_COMPLETE.md) for the exact one-line run instructions.

## Key Features

✨ **Modular Architecture**: Six independent, well-documented modules for maximum flexibility

🔧 **Configuration-Driven**: YAML-based configuration for easy customization without code changes

📊 **Comprehensive Logging**: Structured logging throughout the pipeline for debugging and tracking

🐳 **Docker Support**: Full containerization for reproducible, portable deployments

🧪 **Production-Ready**: Extensive test coverage with CI/CD integration

📚 **Well-Documented**: Complete docstrings, type hints, and comprehensive documentation

## Project Structure

```
MutationScan/
├── src/
│   └── mutation_scan/
│       ├── __init__.py
│       ├── genome_extractor/      # NCBI genome downloading (Entrez API)
│       ├── gene_finder/           # AMR gene detection (ABRicate, BLASTn)
│       ├── sequence_extractor/    # Protein extraction & translation
│       ├── variant_caller/        # Pairwise alignment & mutation calling
│       ├── visualizer/            # PyMOL 3D visualization
│       └── utils/                 # Logging, config, file I/O
├── tests/
│   ├── unit/                      # Unit tests for individual modules
│   ├── integration/               # Integration tests for workflows
│   └── fixtures/                  # Test data and fixtures
├── config/
│   └── config.yaml                # Main configuration file
├── data/
│   ├── genomes/                   # Downloaded genome sequences
│   ├── databases/                 # Custom BLAST databases
│   └── output/                    # Pipeline results
├── docker/
│   └── Dockerfile                 # Container definition
├── .github/
│   └── workflows/
│       └── tests.yml              # CI/CD workflow
├── docs/                          # Documentation
├── setup.py                       # Package setup configuration
├── README.md                      # This file
├── CONTRIBUTING.md                # Contribution guidelines
├── LICENSE                        # MIT License
└── requirements.txt               # Python dependencies
```

## Module Overview

### 🧬 genome_extractor
Handles NCBI Entrez API interactions for automated genome downloading and basic validation.

**Key Classes:**
- `EntrezGenomeDownloader`: Query and download genomes from NCBI
- `GenomeProcessor`: Validate and process downloaded sequences

**Usage:**
```python
from mutation_scan.core import NCBIDatasetsGenomeDownloader

downloader = NCBIDatasetsGenomeDownloader(email="your.email@example.com")
genomes = downloader.search_genomes("Escherichia coli", max_results=10)
```

### 🔍 gene_finder
Wrapper modules for AMR gene detection using industry-standard tools.

**Key Classes:**
- `AbricateWrapper`: ABRicate database screening
- `BLASTWrapper`: Custom BLASTn searches

**Usage:**
```python
from mutation_scan.core import AbricateWrapper

abricate = AbricateWrapper(database="card")
results = abricate.screen_genome("genome.fasta", "results.tsv")
```

### 🧪 sequence_extractor
Coordinate-based sequence extraction and DNA-to-protein translation.

**Key Classes:**
- `CoordinateParser`: Parse GFF/GenBank annotations
- `SequenceTranslator`: DNA to protein translation

**Usage:**
```python
from mutation_scan.core import SequenceTranslator

translator = SequenceTranslator()
protein = translator.translate("ATGAAAGCG...", frame=0)
```

### 🔀 variant_caller
Pairwise sequence alignment and mutation detection.

**Key Classes:**
- `SequenceAligner`: Needleman-Wunsch and Smith-Waterman alignment
- `VariantDetector`: Mutation identification and classification

**Usage:**
```python
from mutation_scan.analysis import SequenceAligner, VariantDetector

aligner = SequenceAligner()
result = aligner.global_alignment(seq1, seq2)

detector = VariantDetector()
variants = detector.detect_variants(result["alignment1"], result["alignment2"])
```

### 🎨 visualizer
PyMOL automation for 3D protein structure visualization.

**Key Classes:**
- `PyMOLAutomator`: PyMOL control and image/movie generation
- `StructureMapper`: Map mutations to 3D coordinates

**Usage:**
```python
from mutation_scan.visualization import PyMOLAutomator

pymol = PyMOLAutomator(headless=True)
pymol.load_structure("protein.pdb")
pymol.highlight_mutations(mutations, color="red")
pymol.generate_image("output.png", dpi=300)
```

### 🛠️ utils
Shared utilities for logging, configuration, and file handling.

**Key Classes:**
- `setup_logger()`: Configure structured logging
- `ConfigParser`: Parse YAML configuration files
- `FileHandler`: Safe file I/O operations

## Installation

### Option 1: Development Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/MutationScan.git
cd MutationScan

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with all dependencies
pip install -e ".[dev]"
```

### Option 2: Production Installation

```bash
pip install mutation-scan
```

### Option 3: Docker Installation

```bash
# Build the container
docker build -t mutation-scan:latest .

# Run the container
docker run -v $(pwd)/data:/app/data mutation-scan:latest
```

## Configuration

MutationScan uses YAML-based configuration for all parameters. Edit `config/config.yaml` to customize:

```yaml
# NCBI Settings
ncbi:
  email: "your.email@example.com"
  api_key: null  # Optional, for faster requests
  max_retries: 3

# ABRicate Settings
abricate:
  database: "card"
  min_coverage: 80
  min_identity: 90

# BLAST Settings
blast:
  database_path: "/path/to/blast/db"
  evalue: 1e-5
  num_threads: 4

# Visualization Settings
visualization:
  pymol_headless: true
  dpi: 300
  color_scheme: "default"

# Output Settings
output:
  results_dir: "data/output"
  format: "json"  # json or csv
```

## Quick Start

### 1. Basic Pipeline Execution

```python
from mutation_scan.utils import setup_logging
from mutation_scan.core import NCBIDatasetsGenomeDownloader
from mutation_scan.core import AbricateWrapper

# Setup logging
from mutation_scan.utils import get_logger
logger = get_logger("mutation_scan")

# Step 1: Download genome
downloader = EntrezGenomeDownloader(config.get("ncbi.email"))
# ... download logic

# Step 2: Screen for AMR genes
abricate = AbricateWrapper(database=config.get("abricate.database"))
# ... screening logic

logger.info("Pipeline completed successfully")
```

### 2. Running Tests

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=src tests/

# Run specific test module
pytest tests/unit/test_sequence_extractor.py
```

### 3. Building Documentation

```bash
cd docs
make html
# Documentation will be in docs/_build/html/
```

## External Dependencies

### Required Tools
- **ABRicate**: `conda install -c bioconda abricate` or `apt-get install abricate`
- **BLAST+**: `conda install -c bioconda blast` or `apt-get install ncbi-blast+`
- **PyMOL** (optional): `conda install -c schrodinger pymol-open-source`

### Python Packages
See `requirements.txt` for the complete list. Main dependencies:
- `biopython`: Biological sequence analysis
- `pyyaml`: YAML configuration parsing
- `pandas`: Data manipulation
- `numpy/scipy`: Scientific computing
- `matplotlib`: Data visualization

## Usage Examples

### Example 1: Download and Screen a Single Genome

```python
from mutation_scan.core import NCBIDatasetsGenomeDownloader
from mutation_scan.core import AbricateWrapper
from mutation_scan.utils import get_logger

logger = get_logger("example1")

# Initialize downloader
downloader = NCBIDatasetsGenomeDownloader(email="user@example.com")

# Search for E. coli genomes
genomes = downloader.search_genomes("Escherichia coli[organism]", max_results=5)
logger.info(f"Found {len(genomes)} genomes")

# Download first genome
if genomes:
    accession = genomes[0]["accession"]
    success = downloader.download_genome(accession, f"data/genomes/{accession}.fasta")
    
    if success:
        # Screen for AMR genes
        abricate = AbricateWrapper(database="card")
        results = abricate.screen_genome(
            f"data/genomes/{accession}.fasta",
            "data/output/amr_results.tsv"
        )
        logger.info("AMR screening completed")
```

### Example 2: Sequence Alignment and Mutation Detection

```python
from mutation_scan.core import SequenceTranslator
from mutation_scan.analysis import SequenceAligner, VariantDetector

# Translate sequences
translator = SequenceTranslator()
seq1_protein = translator.translate(dna_seq1)
seq2_protein = translator.translate(dna_seq2)

# Perform alignment
aligner = SequenceAligner(match_score=2, mismatch_penalty=-1, gap_penalty=-1)
alignment = aligner.global_alignment(seq1_protein, seq2_protein)

# Detect variants
detector = VariantDetector()
variants = detector.detect_variants(
    alignment["alignment1"],
    alignment["alignment2"],
    ref_seq="seq1"
)

# Analyze results
stats = detector.calculate_statistics(variants)
print(f"Found {stats['total_variants']} variants")
```

### Example 3: 3D Visualization

```python
from mutation_scan.visualization import PyMOLAutomator, StructureMapper

# Initialize visualizer
pymol = PyMOLAutomator(headless=True)
pymol.load_structure("protein.pdb")

# Map mutations to structure
mapper = StructureMapper()
mapped_mutations = mapper.map_mutations_to_structure(variants, "protein.pdb")

# Highlight and visualize
pymol.highlight_mutations(mapped_mutations, color="red")
pymol.generate_image("mutation_visualization.png", dpi=300)
```

## Testing

The project includes comprehensive tests organized by type:

```bash
# Run all tests
pytest -v

# Run only unit tests
pytest tests/unit/ -v

# Run only integration tests
pytest tests/integration/ -v

# Generate coverage report
pytest --cov=src --cov-report=html tests/

# Run specific test
pytest tests/unit/test_sequence_extractor.py::TestSequenceTranslator -v
```

## CI/CD Pipeline

GitHub Actions automatically run:
- ✅ Unit and integration tests
- ✅ Code coverage analysis
- ✅ Linting (flake8, black)
- ✅ Type checking (mypy)

See `.github/workflows/tests.yml` for workflow configuration.

## Contributing

We welcome contributions! Please see [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) for:
- Code standards and style guide
- Pull request process
- Bug reporting guidelines
- Feature request templates

## Troubleshooting

### Common Issues

**Issue**: NCBI Entrez connection timeout
```
Solution: Add API key to config.yaml or increase max_retries
```

**Issue**: ABRicate database not found
```
Solution: Run 'abricate --setupdb' to initialize databases
```

**Issue**: PyMOL module not found
```
Solution: Install with: pip install "mutation-scan[pymol]"
```

## Performance Considerations

- **Parallel Processing**: Configure BLAST thread count in config.yaml
- **Memory Usage**: Large genomes may require 8GB+ RAM
- **Network**: NCBI downloads depend on internet connection speed
- **Storage**: Allocate sufficient space for genomes and databases (~50GB+ recommended)

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use MutationScan in your research, please cite:

```bibtex
@software{mutationscan2024,
  author = {Your Organization},
  title = {MutationScan: A Modular Bioinformatics Pipeline for AMR Gene Analysis},
  year = {2024},
  url = {https://github.com/yourusername/MutationScan}
}
```

## Support & Contact

- 📧 **Email**: contact@example.com
- 🐛 **Issues**: [GitHub Issues](https://github.com/yourusername/MutationScan/issues)
- 📖 **Documentation**: [GitHub Wiki](https://github.com/yourusername/MutationScan/wiki)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/yourusername/MutationScan/discussions)

## Acknowledgments

MutationScan builds upon excellent open-source bioinformatics tools:
- [Biopython](https://biopython.org/)
- [ABRicate](https://github.com/tseemann/abricate)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/)
- [PyMOL](https://pymol.org/)

## Changelog

See [CHANGELOG.md](docs/CHANGELOG.md) for a history of changes and releases.

---

**Last Updated**: January 2024 | **Current Version**: 0.1.0 (Alpha)
