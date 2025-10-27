# 🧬 MutationScan

**Comprehensive Antimicrobial Resistance Mutation Analysis Pipeline for researchers**

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](VERSION)
[![Status](https://img.shields.io/badge/status-Production%20Ready-green.svg)](RELEASE_NOTES.md)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](requirements.txt)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 🚀 Quick Start

### 🎯 Choose Your Interface

MutationScan offers **three ways** to run analyses - pick what works best for you:

#### 1. 🖥️ **Web GUI** (Recommended for Beginners)
```bash
python subscan/gui/launch.py
```
- 🖱️ Drag-and-drop file uploads
- 📊 Real-time visual progress tracking
- 🎨 No command-line knowledge needed
- 📋 Pre-configured templates for common AMR genes
- 🔍 Built-in system health checks
- 📥 One-click report viewing and downloads

**Perfect for:** First-time users, visual learners, researchers unfamiliar with command-line tools

👉 **[See Full Web GUI Guide](docs/GUI_USER_GUIDE.md)**

#### 2. 🧙 **CLI Wizard** (Best for Remote Servers)
```bash
python subscan/tools/run_wizard.py
```
- ⌨️ Interactive terminal prompts (no complex commands)
- 📝 Common gene templates (Fluoroquinolone, Beta-Lactamase, etc.)
- 💾 Generates reproducible command files
- 🖥️ Works over SSH (no GUI needed)
- ✅ Input validation and helpful error messages

**Perfect for:** HPC clusters, remote servers, automated workflows, terminal power users

👉 **[See Full CLI Wizard Guide](docs/CLI_WIZARD_GUIDE.md)**

#### 3. ⚙️ **Manual Command-Line** (Advanced Users)
```bash
python subscan/tools/run_pipeline.py \
  --accessions data_input/accession_list.txt \
  --gene-list data_input/gene_list.txt \
  --email you@example.com \
  --output-dir data_output/run1 \
  --sepi-species "Escherichia coli" \
  --threads 4
```
**Perfect for:** Scripting, automation, fine-grained control, integration with other tools

---

### 🏃 Quickest Path (5 Minutes)

1. **Clone the repository:**
   ```bash
   git clone https://github.com/vihaankulkarni29/MutationScan.git
   cd MutationScan
   ```
2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Launch the Web GUI:**
   ```bash
   python subscan/gui/launch.py
   ```
   Browser opens automatically at `http://localhost:5000`

4. **Follow the visual setup:**
   - Check system health (automatic)
   - Select "Quick Start" template
   - Drag-drop your accession file (or use example)
   - Enter your email
   - Click "Start Pipeline"
   - Watch live progress!

**No files ready?** Use the included demo:
```bash
# Copy example files
cp examples/demo_accessions.txt data_input/accession_list.txt

# Launch GUI and select "Quick Start" template
python subscan/gui/launch.py
```

---

### 🛠️ System Requirements

**Minimum:**
- Python 3.8+
- 4 GB RAM
- 10 GB disk space
- Internet connection (for genome downloads)

**Required Software:**
- **EMBOSS** (sequence alignment)
  - Ubuntu/WSL: `sudo apt install emboss`
  - macOS: `brew install emboss`
- **ABRicate** (AMR annotation)
  - `conda install -c bioconda abricate`

**Verify Installation:**
```bash
python subscan/tools/health_check.py
```
All checks should show ✓ (green). The GUI also displays health status automatically.

---

### 📚 Example Workflows

#### Workflow 1: Fluoroquinolone Resistance Analysis (GUI)

1. Launch GUI: `python subscan/gui/launch.py`
2. Select **"Fluoroquinolone Resistance"** template
3. Upload your accession list
4. Enter email, click Start
5. View results when complete

**Genes analyzed:** gyrA, parC, gyrB, parE (quinolone targets)

#### Workflow 2: Beta-Lactamase Study (CLI Wizard)

```bash
python subscan/tools/run_wizard.py

# Follow prompts:
# → Upload accession file
# → Select "Beta-Lactamase" template
# → Enter email
# → Choose database: NCBI
# → Run immediately or save command
```

**Genes analyzed:** blaKPC, blaNDM, blaOXA-48, blaCTX-M, blaTEM, blaSHV

#### Workflow 3: Custom Gene List (Manual CLI)

```bash
# Create custom gene list
cat > my_genes.txt << EOF
mecA
vanA
ermA
tetM
EOF

# Run pipeline
   ```bash
   python subscan/tools/run_pipeline.py --accessions data_input/accession_list.txt --gene-list data_input/gene_list.txt --email you@example.com --output-dir data_output/run1 --sepi-species "Escherichia coli" --threads 4
   ```
```

---

### 🌐 Multi-Database Support

MutationScan supports federated genome extraction from multiple databases:
   - **NCBI** (National Center for Biotechnology Information)  
   - **BV-BRC** (Bacterial and Viral Bioinformatics Resource Center)
   - **EnteroBase** (Salmonella, Escherichia/Shigella, and other Enterobacteriaceae)
   - **PATRIC** (Pathosystems Resource Integration Center)

   The pipeline automatically uses all databases by default for maximum genome coverage. Individual domino tools can specify database sources with the `--database` flag:
   ```bash
   # Use all databases (default, recommended)
   python subscan/tools/run_harvester.py --accessions list.txt --email you@example.com --output-dir results/ --database all
   
   # Use NCBI only (legacy behavior)
   python subscan/tools/run_harvester.py --accessions list.txt --email you@example.com --output-dir results/ --database ncbi
   
   # Use specific database
   python subscan/tools/run_harvester.py --accessions list.txt --email you@example.com --output-dir results/ --database bvbrc
   ```

---

## 🛠️ Installation & Setup

### Quick Setup
```bash
cd subscan/
pip install -e .
```

### 🌍 Multi-Database Genome Extraction
MutationScan uses a **federated genome extractor** supporting multiple microbial databases:

| Database | Coverage | Specialization |
|----------|----------|----------------|
| **NCBI** | Global reference genomes | Comprehensive, curated assemblies |
| **BV-BRC** | Bacterial & viral pathogens | NIAID priority pathogens |
| **EnteroBase** | Enterobacteriaceae | Salmonella, E. coli, Shigella specialization |
| **PATRIC** | Bacterial pathogens | Legacy bacterial genomics platform |

The federated extractor is automatically integrated and requires no additional setup. It gracefully handles:
- **Cross-database accession lookup** (finds genomes across all sources)
- **Metadata standardization** (unified JSON manifest format)
- **Network resilience** (automatic fallback between databases)
- **Platform compatibility** (Windows development, Linux production)

### Domino Tool Architecture
All domino tools are included in the MutationScan package:
- `subscan/tools/run_harvester.py` - Multi-database genome extraction
- `subscan/tools/run_annotator.py` - AMR gene annotation  
- `subscan/tools/run_extractor.py` - Gene sequence extraction
- `subscan/tools/run_aligner.py` - Wild-type alignment
- `subscan/tools/run_analyzer.py` - Mutation detection
- `subscan/tools/run_cooccurrence_analyzer.py` - Resistance pattern analysis
- `subscan/tools/run_reporter.py` - Interactive report generation

---

## 🧩 Pipeline Overview

MutationScan is a 7-stage federated AMR analysis pipeline:
1. **🌐 Genome Harvester**: Downloads genomes from NCBI, BV-BRC, EnteroBase, PATRIC
2. **🧬 Gene Annotator**: Identifies AMR genes using ABRicate
3. **✂️ Sequence Extractor**: Extracts target gene sequences  
4. **🎯 Wild-type Aligner**: Aligns sequences to reference genes
5. **🔍 Mutation Analyzer**: Detects SNPs, insertions, deletions
6. **🔗 Co-occurrence Analyzer**: Discovers resistance gene patterns
7. **📊 Report Generator**: Creates interactive HTML reports

**Output:**
- Interactive HTML dashboard
- CSV mutation summary
- Co-occurrence matrix
- Sequence alignments

---

## 🧪 Example Usage & Results

### Basic Pipeline Run
```bash
python subscan/tools/run_pipeline.py --accessions data_input/accession_list.txt --gene-list data_input/gene_list.txt --email you@example.com --output-dir data_output/run1 --sepi-species "Escherichia coli" --threads 4
```

### Individual Domino with Database Selection
```bash
# Multi-database genome harvesting (recommended)
python subscan/tools/run_harvester.py --accessions examples/demo_accessions.txt --email you@example.com --output-dir test_results/ --database all

# NCBI-only harvesting (legacy mode)  
python subscan/tools/run_harvester.py --accessions examples/demo_accessions.txt --email you@example.com --output-dir test_results/ --database ncbi
```

### Expected Output Structure
```
data_output/run1/
├── 01_harvester_results/
│   ├── genome_manifest.json     # Multi-database genome metadata
│   └── genomes/                 # Downloaded FASTA files
├── 02_annotator_results/
│   ├── annotation_manifest.json # AMR gene annotations
│   └── abricate_results.tsv     # ABRicate output
├── ...                          # Intermediate domino outputs
└── 07_reporter_results/
    ├── final_report.html        # 📊 Interactive dashboard
    └── mutation_summary.csv     # 📋 Structured results
```

---

## 🆘 Troubleshooting

### Common Issues

**🌐 Database Connection Problems**
```bash
# Test federated extractor connectivity
python -c "from subscan.tools.run_harvester import execute_federated_genome_extractor; print('Federated extractor ready')"
```

**📁 Missing Output Files**
- Check input file formats (one accession/gene per line)
- Verify email format (required for NCBI API compliance)
- Ensure output directories are writable

**🔧 Platform-Specific Issues**
- **Windows**: ABRicate is mocked for development (see `mock_abricate.py`)
- **Linux**: Full tool support including native ABRicate
- **macOS**: Use Linux-compatible tooling via Docker/containers

**🚫 Permission Errors**
```powershell
# Windows PowerShell
New-Item -Path "data_output" -ItemType Directory -Force
```

**🌍 Network Issues**
- Multi-database extraction provides automatic fallback
- Check internet connectivity for genome downloads
- Verify email is valid for NCBI API access

**⚠️ Advanced Debugging**
```bash
# Validate pipeline structure
python subscan/check_structure.py

# Test dependencies
python subscan/check_dependencies.py

# Run smoke tests
python subscan/smoke_test_dominos.py
```

---

## ❓ FAQ

**Q: Do I need bioinformatics experience?**
A: No. MutationScan provides a user-friendly pipeline for all researchers.

**Q: What databases are supported?**
A: NCBI, BV-BRC, EnteroBase, and PATRIC via our federated genome extractor. The pipeline automatically searches all databases for maximum genome coverage.

**Q: Can I use my own genome files?**
A: Currently, the pipeline uses accession-based downloading from public databases. Local FASTA file support is planned for future releases.

**Q: Which operating systems are supported?**
A: **Windows** (development with mocked tools), **Linux** (full production support), **macOS** (via containerization). The federated extractor works on all platforms.

**Q: How do I install dependencies?**
A: Simple: `cd subscan && pip install -e .` - All domino tools are included.

**Q: What if a database is down?**
A: The federated extractor automatically tries alternative databases. Use `--database all` for maximum resilience.

**Q: Where do I get help?**
A: Check this README, run diagnostic tools (`check_*.py`), or open a GitHub Issue.

---

## 👨‍💻 Developer & Advanced Info

- See `docs/README.md` for advanced configuration, developer setup, API reference, and contribution guidelines.
- All technical details, performance tuning, and security info are consolidated there.

---

## 📧 Contact & Support

---

## 📖 Documentation

**User Guides:**
- 🖥️ **[Web GUI User Guide](docs/GUI_USER_GUIDE.md)** - Complete guide to the browser interface
   - Drag-and-drop workflows
   - Real-time progress tracking
   - Template system
   - Troubleshooting
   - FAQ

- 🧙 **[CLI Wizard User Guide](docs/CLI_WIZARD_GUIDE.md)** - Interactive terminal interface
   - Step-by-step walkthrough
   - Gene templates
   - Command file generation
   - HPC cluster integration
   - Examples

- ⚙️ **[Main README](README.md)** - You are here!
   - Quick start for all interfaces
   - Installation & setup
   - Pipeline overview
   - Troubleshooting

**Developer Resources:**
- 📘 **[Developer Documentation](docs/README.md)** - Advanced configuration, API reference
- 🤝 **[Contributing Guide](CONTRIBUTING.md)** - How to contribute code
- 🔒 **[Security Policy](SECURITY.md)** - Reporting vulnerabilities
- 📝 **[Changelog](CHANGELOG.md)** - Version history

**Tools & Utilities:**
- `subscan/tools/health_check.py` - System diagnostics
   ```bash
   python subscan/tools/health_check.py --verbose
   ```
- `subscan/tools/run_wizard.py` - Interactive CLI wizard
- `subscan/gui/launch.py` - Web GUI launcher

---

## 🆘 Getting Help

### 1. Check Documentation
- **GUI Users:** See [Web GUI Guide](docs/GUI_USER_GUIDE.md)
- **CLI Users:** See [CLI Wizard Guide](docs/CLI_WIZARD_GUIDE.md)
- **Errors:** Check troubleshooting sections in respective guides

### 2. Run Diagnostics
```bash
# Comprehensive health check
python subscan/tools/health_check.py --verbose

# Verify all dependencies
# → Python version (3.8+)
# → Required packages (BioPython, pandas, Flask)
# → EMBOSS installation
# → ABRicate installation
# → Disk space (10+ GB)
# → Internet connectivity
```

### 3. Review Logs
```bash
# GUI logs: Check browser console (F12)
# Pipeline logs: Located in results/run_*/pipeline.log
tail -n 100 results/run_*/pipeline.log
```

### 4. Report Issues
If problems persist:
1. Open a [GitHub Issue](https://github.com/vihaankulkarni29/MutationScan/issues)
2. Include:
    - Error messages
    - Health check output
    - OS version and Python version
    - Steps to reproduce
3. Tag appropriately: `bug`, `gui`, `wizard`, `pipeline`

### 5. Contact Maintainers
- **Author:** Vihaan Kulkarni
- **Email:** vihaankulkarni29@gmail.com
- **GitHub:** [@vihaankulkarni29](https://github.com/vihaankulkarni29)
- **Issues:** [GitHub Issues](https://github.com/vihaankulkarni29/MutationScan/issues)

---


*MutationScan: Professional, user-centric AMR mutation analysis for everyone.*
