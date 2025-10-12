# 🧬 MutationScan

**Comprehensive Antimicrobial Resistance Mutation Analysis Pipeline for researchers**

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](VERSION)
[![Status](https://img.shields.io/badge/status-Production%20Ready-green.svg)](RELEASE_NOTES.md)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](requirements.txt)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 🚀 Quick Start

1. **Clone the repository:**
   ```bash
   git clone https://github.com/vihaankulkarni29/MutationScan.git
   cd MutationScan
   ```
2. **Install dependencies:**
   ```bash
   cd subscan/
   pip install -e .
   ```
   
   > 📹 **New User?** - For step-by-step video installation guide, see our [Installation Video Guide](docs/INSTALLATION_GUIDE.md)

3. **Prepare input files:**
   - `data_input/accession_list.txt`: NCBI accession numbers (one per line)
   - `data_input/gene_list.txt`: Target gene names (one per line)
4. **Run the pipeline:**
   ```bash
   python subscan/tools/run_pipeline.py --accessions data_input/accession_list.txt --gene-list data_input/gene_list.txt --email you@example.com --output-dir data_output/run1 --sepi-species "Escherichia coli" --threads 4
   ```

   Note on Domino 3 (Extractor): The pipeline now uses the CLI-capable FastaAAExtractor. Before running full analyses, you can verify it is available with either of the following in your active environment:

   ```bash
   fasta_aa_extractor --help
   # or
   python -m fasta_aa_extractor --help
   ```

   If neither works, install the extractor in your environment and ensure the console script is on PATH.

---

## 🛠️ Installation & Domino Tool Setup

MutationScan requires several domino tools. If you see an error about missing scripts, follow these steps:

1. **Clone required domino tools:**
   ```bash
   cd subscan
   # Example for NCBI Genome Extractor
   mkdir -p ncbi_genome_extractor
   cd ncbi_genome_extractor
   git clone https://github.com/vihaankulkarni29/ncbi_genome_extractor.git .
   # Repeat for other tools as needed (see pyproject.toml)
   ```
2. **Verify file locations:**
   - `subscan/ncbi_genome_extractor/ncbi_genome_extractor.py`
   - `subscan/tools/run_harvester.py`, `run_annotator.py`, `run_extractor.py`, `run_aligner.py`, `run_analyzer.py`, `run_cooccurrence_analyzer.py`, `run_reporter.py`
3. **Re-run the pipeline.**

---

## 🧩 Pipeline Overview

MutationScan is a 7-stage pipeline:
1. **Genome Harvester**: Downloads genomes from NCBI
2. **Gene Annotator**: Identifies AMR genes
3. **Sequence Extractor**: Extracts target gene sequences
4. **Wild-type Aligner**: Aligns sequences to reference genes
5. **Mutation Analyzer**: Detects SNPs, insertions, deletions
6. **Co-occurrence Analyzer**: Discovers resistance gene patterns
7. **Report Generator**: Creates interactive HTML reports

**Output:**
- Interactive HTML dashboard
- CSV mutation summary
- Co-occurrence matrix
- Sequence alignments

---

## 🧪 Example Usage

```bash
python subscan/tools/run_pipeline.py --accessions data_input/accession_list.txt --gene-list data_input/gene_list.txt --email you@example.com --output-dir data_output/run1 --sepi-species "Escherichia coli" --threads 4
```

**Results:**
- `data_output/run1/07_reporter_results/final_report.html`
- All intermediate manifests and CSVs in stage folders

---

## 🆘 Troubleshooting

- **Missing domino tool error:** See "Domino Tool Setup" above.
- **No results:** Check input file paths and formats.
- **Permission errors:** Ensure folders are writable.
- **Network errors:** Verify internet connection for NCBI downloads.
- **Other issues:** See FAQ or report on GitHub Issues.

---

## ❓ FAQ

**Q: Do I need bioinformatics experience?**
A: No. MutationScan is designed for all users.

**Q: Can I use my own genome files?**
A: Currently, only NCBI accessions are supported. Local file support is planned.

**Q: How do I install all dependencies?**
A: Use `pip install -e .[dev]` and clone all domino tools as described above.

**Q: Where do I get help?**
A: Read this README, check the Troubleshooting section, or open a GitHub Issue.

---

## 👨‍💻 Developer & Advanced Info

- See `docs/README.md` for advanced configuration, developer setup, API reference, and contribution guidelines.
- All technical details, performance tuning, and security info are consolidated there.

---

## 📧 Contact & Support

- **Author:** Vihaan Kulkarni
- **Email:** vihaankulkarni29@gmail.com
- **GitHub:** [@vihaankulkarni29](https://github.com/vihaankulkarni29)
- **Issues:** [GitHub Issues](https://github.com/vihaankulkarni29/MutationScan/issues)

---


*MutationScan: Professional, user-centric AMR mutation analysis for everyone.*
