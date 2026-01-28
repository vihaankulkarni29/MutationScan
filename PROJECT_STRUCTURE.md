# MutationScan - Project Structure

```
MutationScan/
│
├── 📁 src/
│   └── 📁 mutation_scan/                    # Main package directory
│       ├── __init__.py                      # Package initialization
│       │
│       ├── 📁 genome_extractor/             # Module 1: NCBI genome acquisition
│       │   ├── __init__.py
│       │   ├── entrez_handler.py            # NCBI Entrez API wrapper
│       │   └── genome_processor.py          # Genome validation & processing
│       │
│       ├── 📁 gene_finder/                  # Module 2: AMR gene detection
│       │   ├── __init__.py
│       │   ├── abricate_wrapper.py          # ABRicate database screening
│       │   └── blast_wrapper.py             # Local BLASTn wrapper
│       │
│       ├── 📁 sequence_extractor/           # Module 3: Sequence processing
│       │   ├── __init__.py
│       │   ├── coordinate_parser.py         # GFF/GenBank parsing
│       │   └── translator.py                # DNA → Protein translation
│       │
│       ├── 📁 variant_caller/               # Module 4: Mutation analysis
│       │   ├── __init__.py
│       │   ├── alignment.py                 # Pairwise alignment (global/local)
│       │   └── variant_detector.py          # Mutation detection & classification
│       │
│       ├── 📁 visualizer/                   # Module 5: 3D visualization
│       │   ├── __init__.py
│       │   ├── pymol_automation.py          # PyMOL control & rendering
│       │   └── structure_mapper.py          # Map mutations to 3D coordinates
│       │
│       └── 📁 utils/                        # Module 6: Shared utilities
│           ├── __init__.py
│           ├── logger.py                    # Structured logging setup
│           ├── config_parser.py             # YAML config parsing
│           └── file_handler.py              # Safe file I/O operations
│
├── 📁 tests/
│   ├── __init__.py
│   ├── conftest.py                          # Pytest fixtures & configuration
│   │
│   ├── 📁 unit/                             # Unit tests
│   │   ├── __init__.py
│   │   ├── test_sequence_extractor.py       # Tests for sequence translation
│   │   └── test_utils.py                    # Tests for utilities
│   │
│   ├── 📁 integration/                      # Integration tests
│   │   ├── __init__.py
│   │   └── test_pipeline.py                 # End-to-end pipeline tests
│   │
│   └── 📁 fixtures/                         # Test data
│       └── README.md                        # Test data guidance
│
├── 📁 config/
│   └── config.yaml                          # Main configuration file (YAML)
│
├── 📁 data/
│   ├── README.md                            # Data directory guide
│   ├── 📁 genomes/                          # Downloaded genome sequences
│   │   └── .gitkeep
│   ├── 📁 databases/                        # BLAST & custom databases
│   │   └── .gitkeep
│   └── 📁 output/                           # Pipeline results
│       └── .gitkeep
│
├── 📁 docker/
│   └── README.md                            # Docker utilities guide
│
├── 📁 docs/
│   └── README.md                            # Documentation guide
│
├── 📁 .github/
│   └── 📁 workflows/
│       └── tests.yml                        # CI/CD pipeline (GitHub Actions)
│
├── 📄 setup.py                              # Package setup configuration
├── 📄 README.md                             # Project README (main documentation)
├── 📄 CONTRIBUTING.md                       # Contribution guidelines
├── 📄 LICENSE                               # MIT License
├── 📄 requirements.txt                      # Python dependencies
├── 📄 Dockerfile                            # Container definition (multi-stage)
└── 📄 .gitignore                            # Git ignore patterns

```

## Key Statistics

| Metric | Value |
|--------|-------|
| **Total Modules** | 6 specialized modules |
| **Test Files** | 4 (unit + integration) |
| **Python Packages** | 8 (main + 6 modules + utils) |
| **Configuration Files** | 3 (YAML, requirements.txt, setup.py) |
| **CI/CD Stages** | 6 (test, code-quality, security, docker, integration, docs) |
| **Documentation Files** | 4 (README, CONTRIBUTING, docs, API refs) |
| **Lines of Code (Est.)** | 2000+ (including docstrings) |

## Directory Purposes

### Core Source (`src/mutation_scan/`)
- **genome_extractor**: NCBI Entrez API interaction and genome processing
- **gene_finder**: ABRicate and BLASTn wrappers for AMR gene screening  
- **sequence_extractor**: Coordinate parsing and DNA→Protein translation
- **variant_caller**: Pairwise alignment and mutation detection
- **visualizer**: PyMOL automation and 3D structure visualization
- **utils**: Logging, configuration, and file handling

### Testing (`tests/`)
- **unit/**: Focused tests for individual components
- **integration/**: End-to-end workflow testing
- **fixtures/**: Test data and sample files
- **conftest.py**: Shared pytest fixtures and configuration

### Configuration (`config/`)
- **config.yaml**: Central configuration with sections for NCBI, ABRicate, BLAST, alignment, visualization, output, logging

### Data (`data/`)
- **genomes/**: Downloaded bacterial genome sequences (git-ignored, local only)
- **databases/**: BLAST and custom search databases (git-ignored, local only)
- **output/**: Analysis results and pipeline outputs (git-ignored, local only)

### DevOps
- **Dockerfile**: Multi-stage Docker build (development, testing, production)
- **.github/workflows/tests.yml**: GitHub Actions CI/CD with matrix testing (multi-OS, multi-Python)

## Module Dependencies

```
genome_extractor
    ├─ Biopython (NCBI Entrez)
    ├─ requests
    └─ Pathlib (file handling)

gene_finder
    ├─ subprocess (tool execution)
    └─ External tools (ABRicate, BLAST+)

sequence_extractor
    ├─ Biopython (sequence handling)
    └─ Standard genetic code tables

variant_caller
    ├─ numpy (matrix operations for alignment)
    └─ pandas (results handling)

visualizer
    ├─ pymol-open-source (optional)
    └─ numpy (coordinate manipulation)

utils
    ├─ PyYAML (config parsing)
    ├─ logging (structured logging)
    └─ pathlib (file operations)
```

## Quick Navigation

| Task | File | Location |
|------|------|----------|
| Install package | [setup.py](setup.py) | Root |
| Configure pipeline | [config.yaml](config/config.yaml) | config/ |
| View main docs | [README.md](README.md) | Root |
| Contribute code | [CONTRIBUTING.md](CONTRIBUTING.md) | Root |
| Run tests | [tests.yml](.github/workflows/tests.yml) | .github/workflows/ |
| Deploy container | [Dockerfile](Dockerfile) | Root |
| Main code | [__init__.py](src/mutation_scan/__init__.py) | src/mutation_scan/ |

