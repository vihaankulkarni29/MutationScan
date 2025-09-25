# MutationScan Pipeline: Production Guide for Researchers

## Overview

The **MutationScan Pipeline** is a production-ready bioinformatics tool designed for antimicrobial resistance (AMR) analysis in bacterial genomes. This guide provides comprehensive instructions for researchers and project investigators to deploy and use the pipeline effectively.

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Pipeline Architecture](#pipeline-architecture)
5. [Usage Instructions](#usage-instructions)
6. [Input Data Requirements](#input-data-requirements)
7. [Output Interpretation](#output-interpretation)
8. [Troubleshooting](#troubleshooting)
9. [Advanced Configuration](#advanced-configuration)
10. [Performance Optimization](#performance-optimization)
11. [Citation and Support](#citation-and-support)

## System Requirements

### Minimum Requirements
- **Operating System**: Linux, macOS, or Windows 10+
- **Python**: Version 3.8 or higher
- **RAM**: 8 GB minimum, 16 GB recommended
- **Storage**: 10 GB free space minimum
- **Internet**: Required for initial setup and genome downloads

### Recommended Requirements
- **Python**: Version 3.9-3.12
- **RAM**: 32 GB for large-scale analyses
- **Storage**: 100 GB+ for extensive datasets
- **CPU**: Multi-core processor (4+ cores recommended)

### Required Software Dependencies
```bash
# Core Python packages (auto-installed)
biopython>=1.79
pandas>=1.3.0
requests>=2.25.0
pytest>=7.0.0
pylint>=2.12.0

# External tools (manual installation required)
abricate>=1.0.0
prokka>=1.14.0
mlst>=2.19.0
```

## Installation

### Method 1: Direct Download (Recommended)
```bash
# Clone the repository
git clone https://github.com/yourusername/MutationScan.git
cd MutationScan/subscan

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate     # Windows

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -m pytest tests/ -v
```

### Method 2: Development Setup
```bash
# For contributors and advanced users
git clone https://github.com/yourusername/MutationScan.git
cd MutationScan/subscan

# Create development environment
python -m venv dev_env
source dev_env/bin/activate

# Install in development mode
pip install -e .
pip install -r requirements-dev.txt

# Run full test suite
python -m pytest tests/ --cov=src --cov-report=html
```

## Quick Start

### 30-Second Demo
```bash
# Navigate to pipeline directory
cd MutationScan/subscan

# Run with sample data
python tools/run_harvester.py --input data/sample_accessions.txt --output results/

# View results
ls results/
cat results/genome_manifest.json
```

### Basic Analysis Workflow
```bash
# Step 1: Prepare your data
echo "NZ_CP107554" > my_genomes.txt
echo "NZ_CP107555" >> my_genomes.txt

# Step 2: Run the pipeline
python tools/run_complete_pipeline.py \
  --input my_genomes.txt \
  --output_dir results/ \
  --email your.email@institution.edu

# Step 3: Check results
ls results/
python tools/run_reporter.py --manifest results/final_manifest.json
```

## Pipeline Architecture

The MutationScan pipeline follows a **7-domino architecture**:

```
[Input] → [D1: Harvester] → [D2: Annotator] → [D3: Mutator] → 
[D4: Aggregator] → [D5: Filter] → [D6: Co-occurrence] → [D7: Reporter] → [Output]
```

### Domino Functions
1. **Harvester**: Downloads and validates genome sequences
2. **Annotator**: Identifies genes and resistance markers
3. **Mutator**: Detects mutations and variants
4. **Aggregator**: Combines mutation data across genomes
5. **Filter**: Applies quality and relevance filters
6. **Co-occurrence**: Analyzes mutation pattern associations
7. **Reporter**: Generates publication-ready reports

### Data Flow
- **Input**: Genome accession lists, FASTA files, or custom datasets
- **Intermediate**: JSON manifests tracking pipeline progress
- **Output**: Comprehensive reports, tables, and visualizations

## Usage Instructions

### Command-Line Interface

#### Basic Usage
```bash
# Single domino execution
python tools/run_harvester.py --input accessions.txt --output results/

# Full pipeline execution
python tools/run_complete_pipeline.py --config pipeline_config.yaml

# Generate reports
python tools/run_reporter.py --manifest results/final_manifest.json --format html
```

#### Advanced Options
```bash
# Custom configuration
python tools/run_complete_pipeline.py \
  --input my_genomes.txt \
  --output_dir /path/to/results/ \
  --config custom_config.yaml \
  --threads 8 \
  --memory 16GB \
  --quality_threshold 95.0 \
  --email your.email@domain.com \
  --verbose

# Restart from specific domino
python tools/run_complete_pipeline.py \
  --resume_from domino_3 \
  --manifest results/domino_2_manifest.json
```

### Python API

```python
from src.subscan.pipeline import MutationScanPipeline
from src.subscan.utils import load_json_manifest

# Initialize pipeline
pipeline = MutationScanPipeline(
    output_dir="results/",
    config_file="config.yaml"
)

# Load input data
accessions = ["NZ_CP107554", "NZ_CP107555"]
pipeline.set_input(accessions)

# Run analysis
results = pipeline.run_complete_analysis()

# Access results
print(f"Found {results.total_mutations} mutations")
print(f"Analyzed {results.total_genomes} genomes")

# Generate custom reports
pipeline.generate_report(format="pdf", output="analysis_report.pdf")
```

## Input Data Requirements

### Genome Accession Lists
```
# Format: One accession per line
# Supported formats:
NZ_CP107554.1
GCF_123456789.1
NC_000913.3

# Comments supported
NZ_CP107554  # E. coli strain
NZ_CP107555  # K. pneumoniae strain
```

### FASTA Files
- **Format**: Standard FASTA format
- **Naming**: Consistent sequence identifiers
- **Quality**: Complete or high-quality draft genomes preferred
- **Size**: No strict limits, but >1MB per genome recommended

### Configuration Files
```yaml
# config.yaml
pipeline_settings:
  quality_threshold: 95.0
  min_coverage: 80.0
  max_missing_data: 5.0
  
analysis_parameters:
  mutation_types: ["SNP", "indel"]
  resistance_databases: ["CARD", "ResFinder"]
  gene_families: ["bla", "aac", "tet"]
  
output_options:
  format: ["json", "csv", "html"]
  include_plots: true
  detailed_logs: true
```

## Output Interpretation

### Primary Outputs

#### 1. Final Manifest (final_manifest.json)
```json
{
  "pipeline_info": {
    "version": "1.0.0",
    "execution_date": "2025-01-14",
    "total_genomes": 100,
    "successful_analyses": 98
  },
  "mutation_summary": {
    "total_mutations": 1547,
    "unique_genes": 23,
    "resistance_mechanisms": ["efflux", "enzymatic", "target_modification"]
  },
  "quality_metrics": {
    "average_coverage": 96.2,
    "average_identity": 98.7
  }
}
```

#### 2. Mutation Matrix (mutations_matrix.csv)
```csv
Genome,Gene,Position,Reference,Alternative,Mutation_Type,Frequency,Resistance_Impact
NZ_CP107554,blaOXA-1,123,A,G,SNP,1,High
NZ_CP107555,blaSHV-1,456,C,T,SNP,1,Medium
```

#### 3. Co-occurrence Report (cooccurrence_analysis.html)
- Interactive visualization of mutation patterns
- Statistical significance testing
- Network analysis of gene interactions

### Secondary Outputs

#### Quality Control Reports
- **Coverage Statistics**: Per-genome and per-gene coverage metrics
- **Identity Scores**: Sequence similarity to reference databases
- **Processing Logs**: Detailed execution logs for troubleshooting

#### Visualization Files
- **Mutation Heatmaps**: Visual representation of mutation patterns
- **Phylogenetic Trees**: Genomic relationships based on mutations
- **Resistance Profiles**: AMR gene distribution across genomes

## Troubleshooting

### Common Issues and Solutions

#### 1. Installation Problems
```bash
# Issue: Package conflicts
# Solution: Use clean virtual environment
python -m venv fresh_env
source fresh_env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

# Issue: Missing external tools
# Solution: Install using conda
conda install -c bioconda abricate prokka mlst
```

#### 2. Memory Issues
```bash
# Issue: Out of memory errors
# Solution: Reduce batch size
python tools/run_complete_pipeline.py --batch_size 10 --max_memory 8GB

# Issue: Large genome files
# Solution: Enable compression
python tools/run_complete_pipeline.py --compress_intermediate
```

#### 3. Network Issues
```bash
# Issue: Download failures
# Solution: Enable retry mechanism
python tools/run_harvester.py --retry_count 3 --timeout 300

# Issue: NCBI rate limiting
# Solution: Provide email for higher limits
python tools/run_harvester.py --email your.email@domain.com --delay 1
```

#### 4. Data Quality Issues
```bash
# Issue: Low-quality genomes
# Solution: Adjust quality thresholds
python tools/run_complete_pipeline.py --min_quality 80.0 --exclude_drafts false

# Issue: Missing annotations
# Solution: Force re-annotation
python tools/run_annotator.py --force_reannotation --database_update
```

### Error Code Reference

| Error Code | Description | Solution |
|------------|-------------|----------|
| E001 | Invalid accession format | Check accession list format |
| E002 | Network connection failed | Check internet connection |
| E003 | Insufficient disk space | Free up storage space |
| E004 | Missing dependencies | Install required software |
| E005 | Corrupted input file | Re-download or repair file |

### Getting Help

1. **Check Documentation**: Review this guide and inline help
2. **Run Diagnostics**: `python tools/diagnose_installation.py`
3. **Enable Verbose Logging**: Add `--verbose` flag to commands
4. **Contact Support**: Create issue on GitHub repository

## Advanced Configuration

### High-Performance Computing (HPC)

#### SLURM Configuration
```bash
#!/bin/bash
#SBATCH --job-name=mutationscan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --memory=64GB
#SBATCH --time=24:00:00
#SBATCH --output=mutationscan_%j.log

module load python/3.9
source venv/bin/activate

python tools/run_complete_pipeline.py \
  --input large_dataset.txt \
  --output_dir $SCRATCH/mutationscan_results/ \
  --threads 16 \
  --memory 64GB \
  --config hpc_config.yaml
```

#### Cluster Configuration (cluster_config.yaml)
```yaml
cluster_settings:
  scheduler: "slurm"
  max_jobs: 50
  queue: "normal"
  walltime: "24:00:00"
  
resource_allocation:
  harvester: {memory: "8GB", cpu: 2, time: "2:00:00"}
  annotator: {memory: "16GB", cpu: 4, time: "6:00:00"}
  mutator: {memory: "32GB", cpu: 8, time: "12:00:00"}
```

### Database Configuration

#### Custom Resistance Database
```python
# custom_db_config.py
CUSTOM_DATABASES = {
    "my_resistance_db": {
        "path": "/path/to/my_database.fasta",
        "type": "nucleotide",
        "tool": "blast",
        "parameters": {
            "evalue": 1e-5,
            "identity": 80.0,
            "coverage": 70.0
        }
    }
}
```

### Workflow Customization

#### Custom Analysis Pipeline
```python
from src.subscan.pipeline import BasePipeline
from src.subscan.dominoes import CustomDomino

class CustomMutationScanPipeline(BasePipeline):
    def __init__(self, config):
        super().__init__(config)
        self.add_domino(CustomDomino("my_analysis"))
    
    def run_custom_analysis(self):
        # Your custom analysis logic
        results = self.process_genomes()
        return self.generate_custom_report(results)
```

## Performance Optimization

### System Optimization

#### Memory Management
```python
# config.yaml
memory_settings:
  max_genome_cache: 10
  use_memory_mapping: true
  compress_intermediate: true
  cleanup_frequency: 5
```

#### Parallel Processing
```python
# Enable parallel processing
parallelization:
  enable_multiprocessing: true
  max_processes: 8
  chunk_size: 100
  use_thread_pool: false
```

### Performance Benchmarks

| Dataset Size | Processing Time | Memory Usage | Recommended Hardware |
|--------------|-----------------|--------------|----------------------|
| 10 genomes   | 30 minutes      | 4 GB        | Standard laptop      |
| 100 genomes  | 4 hours         | 16 GB       | Workstation         |
| 1000 genomes | 2 days          | 64 GB       | HPC cluster         |
| 10000 genomes| 2 weeks         | 256 GB      | Large HPC cluster   |

### Optimization Tips

1. **Use SSD Storage**: Significantly improves I/O performance
2. **Optimize Batch Sizes**: Balance memory usage and processing efficiency
3. **Enable Compression**: Reduces storage requirements for large datasets
4. **Use Local Databases**: Avoid repeated downloads of reference data
5. **Monitor Resource Usage**: Use system monitoring tools to identify bottlenecks

## Citation and Support

### How to Cite

If you use MutationScan in your research, please cite:

```
MutationScan: A comprehensive pipeline for antimicrobial resistance mutation analysis
[Authors], [Year]
[Journal], [Volume(Issue)], [Pages]
DOI: [DOI]
```

### Bibtex Entry
```bibtex
@article{mutationscan2025,
  title={MutationScan: A comprehensive pipeline for antimicrobial resistance mutation analysis},
  author={[Authors]},
  journal={[Journal]},
  year={2025},
  volume={X},
  number={Y},
  pages={Z--Z},
  doi={[DOI]}
}
```

### License

This software is licensed under the MIT License. See LICENSE file for details.

### Support and Community

- **GitHub Repository**: https://github.com/yourusername/MutationScan
- **Documentation**: https://mutationscan.readthedocs.io
- **Issue Tracker**: https://github.com/yourusername/MutationScan/issues
- **Discussions**: https://github.com/yourusername/MutationScan/discussions

### Acknowledgments

We thank the bioinformatics community for tools and databases that make this pipeline possible:
- NCBI for genome databases
- ABRicate developers for resistance gene detection
- BioPython community for sequence analysis tools
- All beta testers and contributors

---

**Version**: 1.0.0  
**Last Updated**: January 14, 2025  
**Maintainers**: [Your Research Team]  
**Institution**: [Your Institution]