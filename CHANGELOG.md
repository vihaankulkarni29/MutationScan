# Changelog

All notable changes to MutationScan will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive documentation overhaul with Google-style docstrings
- Complete type hints across all modules and functions
- Professional repository structure and GitHub integration

### Changed
- Enhanced CLI argument parsing in reporter tool
- Improved error handling and validation across pipeline components
- Updated project structure for better maintainability

### Fixed
- Reporter CLI output path handling
- Aligner documentation formatting issues
- Import path resolution for better IDE support

## [1.0.0] - 2025-10-01

### Added
- Complete 7-domino pipeline architecture for AMR analysis
- **Domino 1**: NCBI Genome Harvester with quality scoring
- **Domino 2**: ABRicate Automator for CARD database integration  
- **Domino 3**: FastaAA Extractor for targeted sequence extraction
- **Domino 4**: Wild-type Aligner using EMBOSS alignment tools
- **Domino 5**: Mutation Analyzer with sophisticated detection algorithms
- **Domino 6**: Co-occurrence Analysis for resistance pattern discovery
- **Domino 7**: HTML Report Generator with interactive visualizations
- Comprehensive test suite with pytest framework
- Professional documentation and user guides
- MIT License for open-source distribution
- GitHub integration with issue templates and workflows

### Features
- **Smart Mutation Detection**: SNPs, insertions, deletions, and complex variants
- **Co-occurrence Pattern Analysis**: Resistance gene relationship discovery
- **Quality Scoring**: Automated genome quality assessment and ranking
- **Interactive Visualizations**: Plotly.js charts and MSAViewer integration
- **Parallel Processing**: Multi-core support for large-scale analysis
- **Flexible Input**: Support for NCBI accessions, local files, and custom gene lists
- **Publication-Ready Outputs**: High-quality figures and comprehensive reports

### Pipeline Components
- Robust error handling and progress tracking
- Manifest-based data flow between pipeline stages
- Configurable parameters for different analysis requirements
- Integration with external tools (EMBOSS, ABRicate, CARD database)
- Mock mode support for testing and development

### Technical Implementation
- Python 3.8+ compatibility
- Modern packaging with pyproject.toml
- Comprehensive logging and debugging features
- Memory-efficient processing for large datasets
- Cross-platform support (Windows, macOS, Linux)

### Documentation
- Comprehensive README with installation and usage guides
- API reference documentation
- Step-by-step tutorials with example data
- Contributing guidelines for developers
- Security policy and vulnerability reporting

## [0.9.0] - 2025-09-15

### Added
- Initial pipeline architecture design
- Core domino tool implementations
- Basic testing framework setup
- Preliminary documentation structure

### Changed
- Refactored codebase for modular design
- Improved error handling mechanisms
- Enhanced logging capabilities

## [0.5.0] - 2025-08-01

### Added
- Proof-of-concept implementations
- Basic genome processing capabilities
- Initial alignment and analysis tools
- Foundation testing infrastructure

### Features
- Basic mutation detection algorithms
- Simple report generation
- Command-line interface prototypes
- Integration with NCBI databases

---

## Release Notes

### Version 1.0.0 Highlights

This major release represents the first production-ready version of MutationScan, featuring:

🧬 **Complete Pipeline**: Full 7-domino architecture from genome acquisition to interactive reports
📊 **Advanced Analytics**: Sophisticated mutation detection and co-occurrence analysis
🎨 **Interactive Visualizations**: Professional HTML reports with Plotly.js and MSAViewer
🚀 **Performance**: Multi-core processing and memory-efficient algorithms
📚 **Documentation**: Comprehensive guides, tutorials, and API references
🧪 **Testing**: Robust test suite ensuring reliability and reproducibility

### Breaking Changes

None for initial release.

### Migration Guide

This is the initial production release. No migration is required.

### Dependencies

- Python 3.8+
- EMBOSS Suite for sequence alignment
- Required Python packages listed in requirements.txt
- Optional: ABRicate for enhanced resistance gene annotation

### Known Issues

- Large genome datasets (>1000 genomes) may require significant computational resources
- NCBI API rate limiting may affect genome harvesting for very large batches
- Some visualization features require modern web browsers with JavaScript enabled

### Future Roadmap

- Support for additional resistance databases
- Enhanced machine learning-based mutation prediction
- Cloud deployment options (AWS, Google Cloud)
- Real-time analysis capabilities
- Integration with laboratory information systems (LIMS)

---

For complete details on any release, see the [GitHub releases page](https://github.com/vihaankulkarni29/MutationScan/releases).