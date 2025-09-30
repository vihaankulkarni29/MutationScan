# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Repository Overview

**MutationScan** is a production-ready bioinformatics pipeline for comprehensive antimicrobial resistance (AMR) analysis. It implements a sophisticated "7-Domino Effect" architecture where each pipeline stage seamlessly feeds data to the next stage via JSON manifest files.

The pipeline processes bacterial genome data from NCBI through complete mutation analysis and interactive visualization, designed specifically for antimicrobial resistance research.

## Common Development Commands

### Environment Setup
```powershell
# Windows PowerShell setup
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
pip install -r subscan/requirements.txt

# Install in development mode
cd subscan
pip install -e .
```

```bash
# Linux/macOS setup
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -r subscan/requirements.txt

# Install in development mode
cd subscan
pip install -e .
```

### Testing Commands

```bash
# Run all tests with coverage
python -m pytest subscan/tests/ --cov=src/subscan --cov-report=html --cov-report=term-missing

# Run specific test categories
python -m pytest subscan/tests/ -m "unit"              # Unit tests only
python -m pytest subscan/tests/ -m "integration"      # Integration tests only
python -m pytest subscan/tests/ -m "not slow"         # Skip slow tests

# Run analyzer-specific tests
python -m pytest subscan/analyzer/tests/

# Run tests for a specific domino tool
python -m pytest subscan/tests/test_harvester.py
```

### Pipeline Execution

```bash
# Complete pipeline execution (automated)
python subscan/tools/run_pipeline.py \
  --accessions demo_accessions.txt \
  --gene-list subscan/sample_data/gene_list.txt \
  --email your.email@domain.com \
  --output-dir ./results \
  --sepi-species "Escherichia coli" \
  --threads 8

# Individual domino execution
python subscan/tools/run_harvester.py --accessions demo_accessions.txt --email user@domain.com --output-dir harvester_results
python subscan/tools/run_annotator.py --manifest harvester_results/genome_manifest.json --output-dir annotator_results
python subscan/tools/run_extractor.py --manifest annotator_results/annotation_manifest.json --gene-list genes.txt --output-dir extractor_results
python subscan/tools/run_aligner.py --manifest extractor_results/protein_manifest.json --sepi-species "Escherichia coli" --output-dir aligner_results
python subscan/analyzer/tools/run_analyzer.py --manifest aligner_results/alignment_manifest.json --output-dir analyzer_results
python subscan/tools/run_cooccurrence_analyzer.py --manifest analyzer_results/analysis_manifest.json --output-dir cooccurrence_results
python subscan/tools/run_reporter.py --manifest cooccurrence_results/cooccurrence_manifest.json --output-dir reports
```

### Development Tools

```bash
# Code quality checks
pylint subscan/src/subscan/
black subscan/src/ subscan/tools/ subscan/tests/
isort subscan/src/ subscan/tools/ subscan/tests/
mypy subscan/src/subscan/

# Setup validation and dependency check
python subscan/setup_mutationscan.py
python subscan/setup_mutationscan.py --check-only
python subscan/setup_mutationscan.py --install-deps

# Repository cleanup (remove temp files)
# Windows:
.\cleanup_repository.ps1
# Linux/macOS:
./cleanup_repository.sh
```

## Architecture Overview

### 7-Domino Pipeline Pattern

MutationScan implements a unique "Domino Effect" architecture where each tool produces a JSON manifest that becomes the input for the next tool:

```
📊 Input: Accession List
    ↓
🧬 Domino 1: Harvester (NCBI Genome Extraction)
    ↓ genome_manifest.json
🔍 Domino 2: Annotator (AMR Gene Detection via CARD/ABRicate)
    ↓ annotation_manifest.json  
🧪 Domino 3: Extractor (Protein Sequence Extraction)
    ↓ protein_manifest.json
📏 Domino 4: Aligner (Wild-type Reference Alignment via EMBOSS)
    ↓ alignment_manifest.json
🔬 Domino 5: Analyzer (Mutation Detection & Analysis)
    ↓ analysis_manifest.json
📊 Domino 6: Co-occurrence (Mutation Pattern Analysis)  
    ↓ cooccurrence_manifest.json
📄 Domino 7: Reporter (Interactive HTML Dashboard)
    ↓ Final Report
```

### Key Architectural Components

**Manifest System**: Every domino produces a standardized JSON manifest containing:
- Execution metadata (timestamps, versions, parameters)
- Input/output file paths
- Processing statistics
- Error handling information
- Data handoff structure for next domino

**Utility Layer (`subscan/src/subscan/utils.py`)**:
- Centralized manifest loading/saving functions
- Standard error handling patterns
- Logging configuration
- File validation utilities
- Used by all domino tools to ensure consistency

**Modular Tool Structure**: Each domino tool follows a consistent pattern:
- `subscan/tools/run_*.py` - CLI interfaces for each domino
- Standardized argument parsing and validation
- Progress reporting and status messages
- Error handling with descriptive messages

**Analyzer Engine (`subscan/analyzer/src/analyzer/engine.py`)**:
- Sophisticated mutation detection algorithms
- Reference-driven sequence comparison
- Statistical co-occurrence analysis
- Parallel processing support for large datasets

## Project Structure

```
MutationScan/
├── subscan/                        # Main pipeline package
│   ├── src/subscan/               # Core pipeline modules
│   │   ├── utils.py               # Shared utilities (manifest handling, logging)
│   │   ├── cooccurrence.py        # Co-occurrence analysis logic
│   │   └── reporting.py           # HTML report generation
│   ├── tools/                     # CLI interfaces for each domino
│   │   ├── run_harvester.py       # Domino 1: NCBI genome extraction
│   │   ├── run_annotator.py       # Domino 2: AMR gene annotation
│   │   ├── run_extractor.py       # Domino 3: Protein sequence extraction
│   │   ├── run_aligner.py         # Domino 4: Wild-type alignment
│   │   ├── run_cooccurrence_analyzer.py # Domino 6: Pattern analysis
│   │   ├── run_reporter.py        # Domino 7: Report generation
│   │   └── run_pipeline.py        # Master orchestrator
│   ├── analyzer/                  # Advanced mutation analysis engine
│   │   ├── src/analyzer/engine.py # Core analysis algorithms
│   │   └── tools/run_analyzer.py  # Domino 5: Mutation detection
│   ├── ncbi_genome_extractor/     # NCBI integration tools
│   ├── tests/                     # Comprehensive test suite
│   ├── sample_data/              # Example datasets
│   ├── config.yaml               # Pipeline configuration
│   ├── pyproject.toml            # Package metadata
│   └── requirements.txt          # Python dependencies
├── debug_test/                   # Development test data
├── demo_accessions.txt          # Sample genome accessions
└── requirements.txt             # Root-level dependencies
```

## Configuration and Data Flow

### Configuration Management
- Primary config: `subscan/config.yaml` (JSON format despite .yaml extension)
- Contains pipeline settings, analysis parameters, output options
- System resource limits and performance tuning parameters

### Manifest File Patterns
Each domino produces standardized manifest files:
- **genome_manifest.json**: Harvester → Annotator
- **annotation_manifest.json**: Annotator → Extractor  
- **protein_manifest.json**: Extractor → Aligner
- **alignment_manifest.json**: Aligner → Analyzer
- **analysis_manifest.json**: Analyzer → Co-occurrence
- **cooccurrence_manifest.json**: Co-occurrence → Reporter

### External Dependencies
- **EMBOSS Suite**: Required for sequence alignment (Domino 4)
- **ABRicate/CARD Database**: For resistance gene annotation (Domino 2)
- **NCBI Entrez**: For genome downloading (Domino 1)
- **Python 3.9+**: Core runtime requirement
- **BioPython**: For sequence manipulation
- **Plotly.js**: For interactive visualizations

## Development Guidelines

### Adding New Domino Tools
1. Follow the established CLI pattern in `subscan/tools/`
2. Use `subscan.utils` for manifest handling and logging
3. Implement proper error handling and progress reporting
4. Add comprehensive tests in `subscan/tests/`
5. Update the master orchestrator (`run_pipeline.py`)

### Manifest File Standards
- Always include execution metadata (timestamps, versions)
- Provide absolute file paths for cross-platform compatibility
- Include error handling information and processing statistics
- Follow the established JSON schema patterns

### Testing Strategy
- Unit tests for individual functions and classes
- Integration tests for domino tool interactions
- End-to-end pipeline tests with sample data
- Performance tests for large-scale genomic analyses
- Use pytest markers: `@pytest.mark.unit`, `@pytest.mark.integration`, `@pytest.mark.slow`

### Error Handling Patterns
- Use custom exceptions from `subscan.utils.ManifestError`
- Provide descriptive error messages with context
- Log errors with appropriate levels (ERROR, WARNING, INFO)
- Ensure graceful degradation when possible

## Key Development Notes

**Bioinformatics Context**: This pipeline specifically targets antimicrobial resistance research, focusing on bacterial genome analysis and mutation detection in resistance genes.

**Performance Considerations**: The pipeline is designed for large-scale analyses (hundreds to thousands of genomes) with parallel processing support and memory-efficient algorithms.

**External Tool Integration**: Several dominos integrate with external bioinformatics tools (EMBOSS, ABRicate) - ensure these are properly installed and accessible in PATH.

**Platform Compatibility**: Designed to work on Windows, Linux, and macOS with appropriate dependency management and path handling.

**Production Ready**: This is a mature pipeline with comprehensive error handling, logging, and validation suitable for research publication and clinical applications.
