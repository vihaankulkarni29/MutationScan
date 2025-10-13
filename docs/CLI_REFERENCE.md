# MutationScan CLI Reference

## 🎯 Overview

MutationScan provides both **pipeline orchestration** and **individual domino tools** with comprehensive CLI interfaces. The federated multi-database architecture allows flexible genome source selection for maximum research coverage.

---

## 🚀 Pipeline Orchestrator

### `run_pipeline.py` - Complete AMR Analysis

**Purpose**: Execute the full 7-domino MutationScan pipeline with automated manifest handoff.

```bash
python subscan/tools/run_pipeline.py [REQUIRED_ARGS] [OPTIONAL_ARGS]
```

#### Required Arguments

| Flag | Description | Example |
|------|-------------|---------|
| `--accessions` | Path to genome accession list file | `data_input/accession_list.txt` |
| `--gene-list` | Path to target gene list file | `data_input/gene_list.txt` |
| `--email` | Valid email for NCBI API compliance | `researcher@university.edu` |
| `--output-dir` | Directory for all pipeline results | `data_output/run_20241215` |
| `--sepi-species` OR `--user-reference-dir` | Species name or custom reference directory | `"Escherichia coli"` |

#### Optional Arguments

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--threads` | integer | `4` | CPU threads for parallel processing |
| `--verbose` | flag | `False` | Detailed logging output |
| `--open-report` | flag | `False` | Auto-open HTML report after completion |
| `--dry-run` | flag | `False` | Structure validation with placeholder data |

#### Examples

```bash
# Basic E. coli AMR analysis
python subscan/tools/run_pipeline.py \
  --accessions examples/demo_accessions.txt \
  --gene-list examples/gene_list.txt \
  --email researcher@lab.org \
  --output-dir results/ecoli_study \
  --sepi-species "Escherichia coli" \
  --threads 8

# Dry run for structure validation
python subscan/tools/run_pipeline.py \
  --accessions examples/demo_accessions.txt \
  --gene-list examples/gene_list.txt \
  --email test@example.com \
  --output-dir validation_test \
  --sepi-species "Escherichia coli" \
  --dry-run

# Custom reference sequences
python subscan/tools/run_pipeline.py \
  --accessions strain_collection.txt \
  --gene-list resistance_genes.txt \
  --email researcher@hospital.org \
  --output-dir clinical_isolates \
  --user-reference-dir custom_references/ \
  --verbose --open-report
```

---

## 🧬 Individual Domino Tools

### 1. `run_harvester.py` - Multi-Database Genome Extraction

**Purpose**: Download genome assemblies from federated microbial databases.

```bash
python subscan/tools/run_harvester.py --accessions FILE --email EMAIL --output-dir DIR [--database DB]
```

#### Arguments

| Flag | Required | Choices | Default | Description |
|------|----------|---------|---------|-------------|
| `--accessions` | ✅ | file path | - | Genome accession list (one per line) |
| `--email` | ✅ | valid email | - | Email for NCBI API compliance |
| `--output-dir` | ✅ | directory | - | Output directory for genomes and manifest |
| `--database` | ❌ | `ncbi`, `bvbrc`, `enterobase`, `patric`, `all` | `all` | Database source selection |

#### Database Selection Guide

| Option | Coverage | Use Case |
|--------|----------|----------|
| `all` | **Maximum coverage** (recommended) | Comprehensive research, novel isolates |
| `ncbi` | Global reference genomes | Legacy workflows, curated assemblies |
| `bvbrc` | NIAID priority pathogens | Biodefense research, pathogen studies |
| `enterobase` | Enterobacteriaceae specialists | E. coli, Salmonella, Shigella research |
| `patric` | Legacy bacterial genomes | Historical data access |

#### Examples

```bash
# Multi-database extraction (recommended)
python subscan/tools/run_harvester.py \
  --accessions strain_list.txt \
  --email researcher@university.edu \
  --output-dir genomes_2024 \
  --database all

# NCBI-only extraction (legacy mode)
python subscan/tools/run_harvester.py \
  --accessions ncbi_accessions.txt \
  --email lab@institution.org \
  --output-dir ncbi_genomes \
  --database ncbi

# Enterobacteriaceae-focused research
python subscan/tools/run_harvester.py \
  --accessions ecoli_collection.txt \
  --email microbiologist@research.edu \
  --output-dir entero_study \
  --database enterobase
```

#### Output Structure

```
output-dir/
├── genome_manifest.json          # Standardized metadata with database sources
└── genomes/
    ├── GENOME001.fasta           # Downloaded assemblies
    ├── GENOME002.fasta
    └── ...
```

### 2. `run_annotator.py` - AMR Gene Annotation

**Purpose**: Identify antimicrobial resistance genes using ABRicate.

```bash
python subscan/tools/run_annotator.py --input-dir DIR --gene-list FILE --output-dir DIR [--threads N]
```

#### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--input-dir` | ✅ | Directory containing genome_manifest.json |
| `--gene-list` | ✅ | Target gene list file |
| `--output-dir` | ✅ | Output directory for annotation results |
| `--threads` | ❌ | CPU threads (default: 4) |

#### Example

```bash
python subscan/tools/run_annotator.py \
  --input-dir genomes_2024 \
  --gene-list resistance_genes.txt \
  --output-dir annotations \
  --threads 8
```

### 3. `run_extractor.py` - Gene Sequence Extraction

**Purpose**: Extract target gene sequences from annotated genomes.

```bash
python subscan/tools/run_extractor.py --input-dir DIR --output-dir DIR [--threads N]
```

#### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--input-dir` | ✅ | Directory containing annotation_manifest.json |
| `--output-dir` | ✅ | Output directory for protein sequences |
| `--threads` | ❌ | CPU threads (default: 4) |

### 4. `run_aligner.py` - Wild-type Sequence Alignment

**Purpose**: Align extracted sequences to reference wild-type genes.

```bash
python subscan/tools/run_aligner.py --input-dir DIR --reference-dir DIR --output-dir DIR [--threads N]
```

#### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--input-dir` | ✅ | Directory containing protein_manifest.json |
| `--reference-dir` | ✅ | Directory with reference sequences |
| `--output-dir` | ✅ | Output directory for alignments |
| `--threads` | ❌ | CPU threads (default: 4) |

### 5. `run_analyzer.py` - Mutation Detection

**Purpose**: Identify SNPs, insertions, and deletions in aligned sequences.

```bash
python subscan/tools/run_analyzer.py --input-dir DIR --output-dir DIR [--threads N]
```

### 6. `run_cooccurrence_analyzer.py` - Resistance Pattern Analysis

**Purpose**: Discover co-occurrence patterns between resistance genes.

```bash
python subscan/tools/run_cooccurrence_analyzer.py --input-dir DIR --output-dir DIR [--min-support FLOAT]
```

#### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--input-dir` | ✅ | - | Directory containing analysis_manifest.json |
| `--output-dir` | ✅ | - | Output directory for co-occurrence results |
| `--min-support` | ❌ | `0.1` | Minimum support threshold (0.0-1.0) |

### 7. `run_reporter.py` - Interactive Report Generation

**Purpose**: Create comprehensive HTML reports with visualizations.

```bash
python subscan/tools/run_reporter.py --input-dir DIR --output-dir DIR [--report-title STR]
```

#### Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--input-dir` | ✅ | - | Directory containing cooccurrence_manifest.json |
| `--output-dir` | ✅ | - | Output directory for final report |
| `--report-title` | ❌ | `"MutationScan Analysis Report"` | Custom report title |

---

## 🔧 Diagnostic and Validation Tools

### Structure Validation
```bash
# Check pipeline component structure
python subscan/check_structure.py

# Expected output: "✅ All 29 domino script files found"
```

### Dependency Verification
```bash
# Verify Python packages and system tools
python subscan/check_dependencies.py

# Platform-specific behavior:
# Windows: Validates Python packages, mocks bioinformatics tools
# Linux: Validates full toolchain including ABRicate, BLAST+, MAFFT
```

### Smoke Testing
```bash
# Test all domino tools with --help flags
python subscan/smoke_test_dominos.py

# Validates CLI interfaces for 8 tools + orchestrator
```

### Comprehensive Testing
```bash
# Full end-to-end validation with multiple test scenarios
python subscan/comprehensive_final_tests.py

# Tests:
# - Parameter variations (dataset sizes, gene panels)  
# - Database source variations (NCBI vs federated)
# - Performance characteristics (scaling analysis)
# - Output structure validation (manifest schemas)
```

### Quick Validation
```bash
# Run all validation tests in sequence
python subscan/quick_validation.py

# Orchestrates: structure → dependencies → smoke → full → comprehensive
```

---

## 🌐 Database-Specific CLI Patterns

### NCBI-Focused Workflows
```bash
# Traditional NCBI-centric research
python subscan/tools/run_harvester.py --database ncbi --accessions ncbi_list.txt --email lab@edu --output-dir ncbi_study

# Continue with standard pipeline
python subscan/tools/run_pipeline.py --accessions ncbi_list.txt --gene-list amr_genes.txt --email lab@edu --output-dir ncbi_study --sepi-species "Klebsiella pneumoniae"
```

### Multi-Database Research
```bash
# Maximum genome coverage approach
python subscan/tools/run_harvester.py --database all --accessions diverse_collection.txt --email researcher@org --output-dir federated_study

# Specialized pathogen research
python subscan/tools/run_harvester.py --database bvbrc --accessions biodefense_isolates.txt --email defense@lab.gov --output-dir pathogen_analysis
```

### Enterobacteriaceae Studies
```bash
# E. coli/Salmonella specialist workflows
python subscan/tools/run_harvester.py --database enterobase --accessions entero_strains.txt --email micro@university.edu --output-dir enterobacteria

# Species-specific pipeline continuation
python subscan/tools/run_pipeline.py --accessions entero_strains.txt --gene-list enteric_resistance.txt --email micro@university.edu --output-dir enterobacteria --sepi-species "Salmonella enterica"
```

---

## 📊 Output Manifest Schemas

### Genome Manifest (`genome_manifest.json`)
```json
{
  "version": "1.0.0",
  "timestamp": "2024-12-15T10:30:00Z",
  "database_sources": ["ncbi", "bvbrc", "enterobase"],
  "genomes": [
    {
      "accession": "NZ_CP123456",
      "source_database": "ncbi",
      "species": "Escherichia coli",
      "file_path": "genomes/NZ_CP123456.fasta",
      "download_timestamp": "2024-12-15T10:31:15Z",
      "metadata": {
        "assembly_level": "Complete Genome",
        "organism": "Escherichia coli strain ABC123"
      }
    }
  ]
}
```

### Annotation Manifest (`annotation_manifest.json`)
```json
{
  "version": "1.0.0",
  "timestamp": "2024-12-15T10:45:00Z",
  "input_genomes": 15,
  "target_genes": ["blaTEM", "blaCTX-M", "qnrS"],
  "annotation_results": [
    {
      "genome": "NZ_CP123456",
      "genes_found": ["blaTEM-1", "qnrS1"],
      "abricate_output": "abricate_results.tsv"
    }
  ]
}
```

---

## 🚨 Error Handling & Exit Codes

### Common Exit Codes

| Code | Meaning | Troubleshooting |
|------|---------|-----------------|
| `0` | Success | Pipeline completed successfully |
| `1` | General error | Check input files and arguments |
| `2` | File not found | Verify file paths and permissions |
| `3` | Network error | Check internet connectivity for genome downloads |
| `4` | Database error | Try alternative database with `--database` flag |
| `5` | Tool not found | Install missing bioinformatics tools (Linux) |

### Error Message Patterns

```bash
# Input validation errors
Error: Accessions file missing: /path/to/file.txt
Solution: Verify file path and ensure file exists

# Database connectivity issues  
Error: Failed to download from NCBI, trying BV-BRC...
Solution: This is normal federated fallback behavior

# Permission errors
Error: Cannot write to output directory: /restricted/path
Solution: Ensure write permissions or choose different output directory

# Email validation
Error: Invalid email (must contain '@')
Solution: Provide valid email address for NCBI API compliance
```

---

## 🔍 Advanced Usage Patterns

### Batch Processing
```bash
# Process multiple studies in parallel
for study in study1 study2 study3; do
  python subscan/tools/run_pipeline.py \
    --accessions "${study}_accessions.txt" \
    --gene-list "${study}_genes.txt" \
    --email batch@research.org \
    --output-dir "batch_results/${study}" \
    --sepi-species "Escherichia coli" &
done
wait  # Wait for all background jobs to complete
```

### Performance Optimization
```bash
# High-performance configuration
python subscan/tools/run_pipeline.py \
  --accessions large_dataset.txt \
  --gene-list comprehensive_genes.txt \
  --email hpc@supercomputer.edu \
  --output-dir hpc_analysis \
  --sepi-species "Klebsiella pneumoniae" \
  --threads 32 \
  --verbose
```

### Custom Reference Integration
```bash
# Use laboratory-specific reference sequences
python subscan/tools/run_pipeline.py \
  --accessions clinical_isolates.txt \
  --gene-list hospital_resistance_panel.txt \
  --email clinicallab@hospital.org \
  --output-dir clinical_analysis \
  --user-reference-dir lab_references/ \
  --threads 16
```

---

## 📞 CLI Support Resources

### Quick Help
```bash
# General pipeline help
python subscan/tools/run_pipeline.py --help

# Individual domino help
python subscan/tools/run_harvester.py --help
python subscan/tools/run_annotator.py --help
# ... etc for all domino tools
```

### Validation Commands
```bash
# Validate installation
python subscan/check_structure.py && echo "✅ Structure OK"
python subscan/check_dependencies.py && echo "✅ Dependencies OK"

# Test pipeline components
python subscan/smoke_test_dominos.py && echo "✅ CLI interfaces OK"
```

### Performance Testing
```bash
# Benchmark your system
time python subscan/comprehensive_final_tests.py
# Expected output includes performance metrics and pass/fail status
```

---

*MutationScan CLI Reference - Complete command-line interface documentation for federated AMR analysis*