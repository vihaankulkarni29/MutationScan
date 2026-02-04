# MutationScan Architecture Documentation

## Directory Structure

```
.
├── .github/workflows/          # CI/CD pipelines
├── config/                     # Configuration files
│   ├── drug_mapping.json       # Gene-to-drug mapping for phenotype interpretation
│   ├── settings.yaml           # Global application settings
│   ├── target_genes.txt        # Default target genes list
│   └── README.md               # Configuration documentation
├── data/                       # Data directory (runtime)
│   ├── genomes/               # Downloaded genome FASTA files
│   ├── proteins/              # Extracted protein sequences
│   ├── refs/                  # Reference sequences and databases
│   ├── results/               # Analysis output (CSV, visualizations)
│   └── logs/                  # Pipeline execution logs
├── docs/                       # Documentation
│   ├── SETUP_COMPLETE.md      # Setup and deployment guide
│   └── API.md                 # Module API documentation
├── models/                     # Trained ML models (.pkl files)
├── src/
│   └── mutation_scan/          # Main Package
│       ├── __init__.py
│       ├── __main__.py         # Module entry point (python -m mutation_scan.main)
│       ├── main.py             # Pipeline orchestrator
│       ├── core/               # Genome Ingestion & Gene Detection
│       │   ├── __init__.py
│       │   ├── entrez_handler.py      # NCBI Datasets API downloader
│       │   ├── genome_processor.py     # Genome validation and preprocessing
│       │   ├── gene_finder.py          # Hybrid ABRicate + BLAST gene detection
│       │   ├── abricate_wrapper.py     # ABRicate resistance gene wrapper
│       │   ├── blast_wrapper.py        # BLAST housekeeping gene wrapper
│       │   ├── reference_builder.py    # Reference FASTA builder
│       │   ├── sequence_extractor.py   # Coordinate-based sequence extraction
│       │   ├── coordinate_parser.py    # GFF/GenBank parser
│       │   └── translator.py           # DNA → Protein translation
│       ├── analysis/           # Variant Calling & ML
│       │   ├── __init__.py
│       │   ├── variant_caller.py       # Main variant calling orchestrator
│       │   ├── variant_detector.py     # Mutation detection logic
│       │   ├── alignment.py            # Global alignment engine
│       │   └── ml/                     # ML Sub-package
│       │       ├── __init__.py
│       │       ├── inference.py        # ML model loading and prediction
│       │       ├── features.py         # Feature extraction from sequences
│       │       ├── trainer.py          # Model training pipeline
│       │       ├── dataset.py          # Data pipeline for ML
│       │       └── benchmark.py        # Model evaluation
│       ├── visualization/      # 3D Structure Visualization
│       │   ├── __init__.py
│       │   ├── pymol_viz.py            # PyMOL automation
│       │   ├── structure_mapper.py     # PDB structure mapping
│       │   └── pymol_automation.py     # PyMOL scripting helpers
│       └── utils/              # Shared Utilities
│           ├── __init__.py
│           ├── logger.py               # Centralized logging
│           ├── config_parser.py        # Configuration file parser
│           └── file_handler.py         # File I/O utilities
├── tests/                      # Unit and Integration Tests
│   ├── unit/
│   └── integration/
├── Dockerfile                  # Production container
├── docker-compose.yml          # Multi-container orchestration
├── requirements.txt            # Python dependencies
├── setup.py                    # Package installation
├── README.md                   # Project overview
└── SIMPLIFIED_RUN_GUIDE.md    # Quick start guide
```

## Package Organization

### Core Module (`mutation_scan.core`)
**Purpose**: Genome ingestion and gene detection

**Key Components**:
- **entrez_handler**: Downloads genomes from NCBI using Datasets API v2
- **gene_finder**: Hybrid detection using ABRicate (resistance genes) + BLAST (housekeeping genes)
- **sequence_extractor**: Extracts genomic regions based on coordinates
- **translator**: Performs DNA → Protein translation with codon table support

**Exports**:
```python
from mutation_scan.core import (
    NCBIDatasetsGenomeDownloader,
    GenomeProcessor,
    GeneFinder,
    SequenceExtractor
)
```

### Analysis Module (`mutation_scan.analysis`)
**Purpose**: Variant calling and machine learning predictions

**Key Components**:
- **variant_caller**: Main orchestrator for mutation detection
- **variant_detector**: Implements residue counter algorithm for gap-aware mutation calling
- **alignment**: Global alignment using BioPython PairwiseAligner + BLOSUM62
- **ml/**: Machine learning subpackage for resistance prediction

**ML Subpackage** (`mutation_scan.analysis.ml`):
- **inference**: Loads trained models and makes predictions
- **features**: Extracts sequence features (AAC, CTD, physicochemical properties)
- **trainer**: Model training pipeline (Random Forest, SVM, XGBoost)
- **dataset**: Data pipeline for model training
- **benchmark**: Model evaluation metrics

**Exports**:
```python
from mutation_scan.analysis import VariantCaller
from mutation_scan.analysis.ml import MLPredictor, FeatureExtractor
```

### Visualization Module (`mutation_scan.visualization`)
**Purpose**: 3D structure visualization with PyMOL

**Key Components**:
- **pymol_viz**: Main visualization orchestrator
- **structure_mapper**: Maps mutations to PDB structures
- **pymol_automation**: PyMOL scripting utilities

**Exports**:
```python
from mutation_scan.visualization import PyMOLVisualizer
```

### Utils Module (`mutation_scan.utils`)
**Purpose**: Shared utilities across all modules

**Key Components**:
- **logger**: Centralized logging configuration
- **config_parser**: YAML/JSON configuration file parsing
- **file_handler**: File I/O utilities with error handling

**Exports**:
```python
from mutation_scan.utils import setup_logging, get_logger
```

## Entry Points

### Command-Line Execution
```bash
# Run as module (recommended)
python -m mutation_scan.main --email user@example.com --genome data/genomes/genome.fasta

# Run in Docker
docker run -v $(pwd)/data:/app/data mutationscan:v1 --email user@example.com --organism "E. coli"
```

### Docker Entry Point
The Dockerfile uses `python -m mutation_scan.main` as the entry point, ensuring proper module resolution.

## Configuration

### Config Directory (`config/`)
- **drug_mapping.json**: Maps gene names to antibiotics
- **settings.yaml**: Global application settings
- **target_genes.txt**: Default genes to analyze

### Environment Variables
- `NCBI_EMAIL`: Email for NCBI API (required)
- `NCBI_API_KEY`: Optional API key for higher rate limits

## Import Conventions

### Old Structure (Deprecated)
```python
from mutation_scan.genome_extractor import NCBIDatasetsGenomeDownloader
from mutation_scan.gene_finder import GeneFinder
from mutation_scan.visualizer import PyMOLVisualizer
```

### New Structure (Current)
```python
from mutation_scan.core.entrez_handler import NCBIDatasetsGenomeDownloader
from mutation_scan.core.gene_finder import GeneFinder
from mutation_scan.visualization.pymol_viz import PyMOLVisualizer
```

## Design Principles

1. **Modular Architecture**: Each package has a single responsibility
2. **Clear Separation of Concerns**: Core (ingestion) → Analysis (detection) → Visualization (output)
3. **Industry-Standard Layout**: Follows src-layout packaging conventions
4. **Configuration-Driven**: No hardcoded paths, all configs in `config/` directory
5. **Pythonic Imports**: Use relative imports within packages, absolute imports across packages
6. **Docker-First**: Optimized for containerized deployment

## Testing Strategy

### Unit Tests (`tests/unit/`)
- Test individual functions in isolation
- Mock external dependencies (NCBI API, ABRicate, BLAST)

### Integration Tests (`tests/integration/`)
- Test full pipeline workflows
- Use small reference datasets
- Validate end-to-end functionality

## CI/CD Pipeline

### GitHub Actions (`.github/workflows/`)
- **Linting**: flake8, black, mypy
- **Testing**: pytest with coverage reports
- **Docker Build**: Automated container builds
- **Release**: Automated versioning and tagging

## Deployment

### Docker Deployment (Recommended)
```bash
docker build -t mutationscan:v2 .
docker run -v $(pwd)/data:/app/data mutationscan:v2 --help
```

### Native Deployment
```bash
pip install -e .
python -m mutation_scan.main --help
```

## Migration Guide

### Updating Imports
If upgrading from the old structure, update imports as follows:

| Old Import | New Import |
|------------|------------|
| `from mutation_scan.genome_extractor import NCBIDatasetsGenomeDownloader` | `from mutation_scan.core.entrez_handler import NCBIDatasetsGenomeDownloader` |
| `from mutation_scan.gene_finder import GeneFinder` | `from mutation_scan.core.gene_finder import GeneFinder` |
| `from mutation_scan.sequence_extractor import SequenceExtractor` | `from mutation_scan.core.sequence_extractor import SequenceExtractor` |
| `from mutation_scan.variant_caller import VariantCaller` | `from mutation_scan.analysis.variant_caller import VariantCaller` |
| `from mutation_scan.visualizer import PyMOLVisualizer` | `from mutation_scan.visualization.pymol_viz import PyMOLVisualizer` |
| `from mutation_scan.ml_predictor import MLPredictor` | `from mutation_scan.analysis.ml.inference import MLPredictor` |

### Updating Configuration Paths
- Old: `data/config/drug_mapping.json`
- New: `config/drug_mapping.json`

Update all references in code and scripts accordingly.

## Version History

### v2.0 (Current) - February 2026
- Refactored to src-layout structure
- Moved configurations to root `config/` directory
- Renamed ML files: `train_zoo.py` → `trainer.py`, `data_pipeline.py` → `dataset.py`
- Renamed visualization: `visualizer.py` → `pymol_viz.py`
- Updated Docker entry point to use module execution
- Centralized logging in `utils/logger.py`

### v1.0 - January 2026
- Initial release with flat structure
- Basic pipeline functionality
- Docker deployment support
