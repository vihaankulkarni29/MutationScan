# MutationScan Web GUI - User Guide

**Version 1.0** | **Last Updated: October 2025**

---

## Table of Contents

1. [Overview](#overview)
2. [System Requirements](#system-requirements)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Detailed Walkthrough](#detailed-walkthrough)
6. [Features & Interface](#features--interface)
7. [Templates & Presets](#templates--presets)
8. [Troubleshooting](#troubleshooting)
9. [FAQ](#faq)

---

## Overview

The MutationScan Web GUI provides a user-friendly browser interface for running the complete AMR (Antimicrobial Resistance) mutation analysis pipeline. No command-line expertise required!

**What it does:**
- Downloads bacterial genomes from NCBI
- Extracts and aligns target AMR genes
- Identifies mutations compared to reference sequences
- Analyzes co-occurrence patterns
- Generates interactive HTML reports

**Key Benefits:**
- 🖱️ **Drag-and-drop** file uploads
- 📊 **Real-time progress** tracking across 7 pipeline stages
- 🎨 **Visual interface** - no terminal commands needed
- 📝 **Pre-configured templates** for common AMR genes
- 🔍 **System health checks** before pipeline execution
- 📥 **One-click downloads** of results and reports

---

## System Requirements

### Minimum Requirements

- **Operating System:** Linux (Ubuntu 20.04+), macOS (10.15+), or Windows 10/11 with WSL2
- **Python:** 3.8 or higher
- **RAM:** 4 GB minimum, 8 GB recommended
- **Disk Space:** 10 GB free space minimum (more for large datasets)
- **Internet:** Required for downloading genomes from NCBI

### Required Software

The following tools must be installed:

1. **Python 3.8+** with pip
2. **EMBOSS** (for sequence alignment)
   - Ubuntu/WSL: `sudo apt install emboss`
   - macOS: `brew install emboss`
3. **ABRicate** (for AMR gene annotation)
   - Installation: `conda install -c bioconda abricate` or follow [ABRicate installation guide](https://github.com/tseemann/abricate)

### Required Python Packages

Install via pip:

```bash
pip install biopython pandas flask
```

**Optional (recommended):**
```bash
pip install psutil  # For memory monitoring
```

---

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan
```

### Step 2: Install Python Dependencies

```bash
pip install -r requirements.txt
```

### Step 3: Verify Installation

Run the health check to ensure all dependencies are installed:

```bash
python subscan/tools/health_check.py
```

You should see all checks marked with ✓ (green checkmarks). Warnings are acceptable, but errors must be resolved.

### Step 4: Launch the GUI

```bash
python subscan/gui/launch.py
```

The browser should automatically open to `http://localhost:5000`. If not, manually navigate to this URL.

---

## Quick Start

### 5-Minute Tutorial

1. **Launch the GUI:**
   ```bash
   python subscan/gui/launch.py
   ```

2. **Check System Health:**
   - The setup page automatically displays system status
   - Ensure critical components (Python, EMBOSS, ABRicate) show ✓

3. **Choose a Template:**
   - Click **"Quick Start"** for a demo analysis
   - Pre-configured with common AMR genes: `gyrA`, `parC`, `blaKPC`

4. **Upload Accessions:**
   - Create a text file with NCBI accession numbers (one per line)
   - Example: `NZ_CP012345.1`, `NC_000913.3`
   - Drag and drop into the **Accession IDs** zone

5. **Upload Genes (Optional):**
   - If using a template, this is pre-filled
   - Or upload a custom gene list (one gene name per line)

6. **Configure Settings:**
   - Email: Required for NCBI downloads (use your real email)
   - Database: Select NCBI (recommended) or alternatives
   - Advanced options: Usually defaults are fine

7. **Start Pipeline:**
   - Click **"Start Pipeline"**
   - You'll be redirected to the live progress page

8. **Monitor Progress:**
   - Watch real-time updates across 7 stages
   - View live logs in the terminal viewer
   - See estimated time remaining

9. **View Results:**
   - Automatically redirected when complete
   - Click a result to view the interactive HTML report
   - Download report or CSV data files

**Total Time:** 5-30 minutes depending on dataset size

---

## Detailed Walkthrough

### Page 1: Setup & Configuration

#### System Health Dashboard

Located at the top of the setup page, this displays:

- **Python:** Version and status
- **Python Packages:** BioPython, pandas, Flask
- **EMBOSS:** Sequence alignment toolkit
- **ABRicate:** AMR database tool
- **Disk Space:** Available storage
- **Memory:** Available RAM (if psutil installed)
- **Internet:** NCBI connectivity

**Status Indicators:**
- ✓ Green = Ready
- ⚠ Yellow = Warning (may limit some features)
- ✗ Red = Error (must fix before proceeding)

#### File Upload Zones

**Accession IDs File:**
- Format: Plain text, one accession per line
- Supported: RefSeq (NZ_*), GenBank (NC_*, CP_*, etc.)
- Example content:
  ```
  NZ_CP012345.1
  NZ_CP067890.1
  NC_000913.3
  ```
- **Upload Methods:**
  - Drag and drop file onto the blue zone
  - Click "Browse" to select file manually

**Target Genes File:**
- Format: Plain text, one gene name per line
- Example content:
  ```
  gyrA
  parC
  blaKPC
  blaNDM
  ```
- **Note:** Template selection auto-fills this

#### Template Selection

Pre-configured gene sets for common AMR analysis:

1. **Quick Start** (Beginner-friendly)
   - Genes: `gyrA`, `parC`, `blaKPC`
   - Database: NCBI
   - Use case: Demo or general fluoroquinolone + beta-lactamase analysis

2. **Fluoroquinolone Resistance**
   - Genes: `gyrA`, `parC`, `gyrB`, `parE`
   - Target: DNA gyrase and topoisomerase mutations
   - Use case: Studying quinolone resistance mechanisms

3. **Beta-Lactamase**
   - Genes: `blaKPC`, `blaNDM`, `blaOXA-48`, `blaCTX-M`
   - Target: Major carbapenemase and ESBL genes
   - Use case: Carbapenem and cephalosporin resistance

**Selecting a Template:**
- Click the template card
- Gene list is automatically populated
- Database is pre-selected
- You can still modify settings as needed

#### Advanced Configuration

Click **"Advanced Options"** to expand:

- **Database Selection:**
  - `NCBI`: National Center for Biotechnology Information (recommended)
  - `BV-BRC`: Bacterial and Viral Bioinformatics Resource Center
  - `EnteroBase`: Specialized for Enterobacteriaceae
  - `PATRIC`: Legacy bacterial database

- **Number of Threads:**
  - Default: 4 (suitable for most systems)
  - Increase for faster processing if you have more CPU cores
  - Decrease if system becomes unresponsive

- **Quality Threshold:**
  - Default: 30 (high quality)
  - Range: 0-50 (higher = stricter quality filtering)
  - Lower for older/lower-quality assemblies

- **Output Directory:**
  - Default: `results/`
  - Results saved in timestamped subdirectories
  - Example: `results/run_20251027_143022/`

- **Email Address:**
  - **Required** for NCBI API compliance
  - Used to track download requests
  - No spam - only for NCBI's logging

#### Starting the Pipeline

1. Ensure all required files are uploaded (green checkmarks)
2. Verify system health shows no critical errors
3. Click **"Start Pipeline"** button
4. Browser automatically navigates to progress page

---

### Page 2: Pipeline Execution & Progress

#### Visual Stage Tracker

The pipeline executes in 7 sequential stages:

```
┌─────────────┐   ┌─────────────┐   ┌─────────────┐   ┌─────────────┐
│  Harvester  │ → │  Annotator  │ → │  Extractor  │ → │   Aligner   │
└─────────────┘   └─────────────┘   └─────────────┘   └─────────────┘
       ↓                                                      ↓
┌─────────────┐   ┌─────────────┐   ┌─────────────┐
│   Reporter  │ ← │Co-Occurrence│ ← │   Analyzer  │
└─────────────┘   └─────────────┘   └─────────────┘
```

**Stage Details:**

1. **🔽 Genome Harvester** (Stage 1)
   - Downloads genome assemblies from NCBI/database
   - Validates accession numbers
   - Stores in `genomes/` directory
   - Progress: Per-genome download status

2. **🏷️ Gene Annotator** (Stage 2)
   - Uses ABRicate to annotate AMR genes
   - Searches against multiple resistance databases
   - Generates annotation reports
   - Progress: Annotation completion percentage

3. **✂️ Sequence Extractor** (Stage 3)
   - Extracts target gene sequences from genomes
   - Creates FASTA files for each gene
   - Handles multiple copies/variants
   - Progress: Per-gene extraction status

4. **📏 Sequence Aligner** (Stage 4)
   - Aligns extracted sequences to references using EMBOSS Needle
   - Pairwise global alignment for each gene
   - Generates alignment files
   - Progress: Alignment completion per genome

5. **🧬 Mutation Analyzer** (Stage 5)
   - Identifies SNPs and amino acid changes
   - Calculates mutation frequencies
   - Maps to protein domains
   - Progress: Analysis completion percentage

6. **🔗 Co-occurrence Analyzer** (Stage 6)
   - Finds mutation patterns that appear together
   - Statistical correlation analysis
   - Network visualization data
   - Progress: Pattern analysis status

7. **📄 Report Generator** (Stage 7)
   - Compiles all results into interactive HTML
   - Creates downloadable CSV tables
   - Generates summary statistics
   - Progress: Report compilation status

#### Stage Status Indicators

Each stage card displays:

- **Icon:** Stage-specific visual indicator
- **Status Badge:**
  - `Pending` (gray): Not yet started
  - `Running` (blue, animated): Currently executing
  - `Completed` (green): Successfully finished
  - `Error` (red): Failed (see logs for details)
- **Progress Bar:** 0-100% completion
- **Current Item:** Which genome/gene is being processed

#### Live Log Viewer

Real-time terminal output appears in the log viewer:

- **Color-coded messages:**
  - Blue = Info
  - Green = Success
  - Yellow = Warning
  - Red = Error

- **Auto-scroll:** Automatically scrolls to latest message
- **Searchable:** Use browser's Find (Ctrl+F) to search logs
- **Copy-friendly:** Select and copy for troubleshooting

#### Progress Statistics

Top of the page shows:

- **Elapsed Time:** How long the pipeline has been running
- **Current Stage:** Which of 7 stages is active
- **Overall Progress:** Combined progress across all stages
- **Estimated Time Remaining:** Calculated from current pace

#### Stopping the Pipeline

- Click **"Stop Pipeline"** button to halt execution
- Partial results may be saved depending on stage
- Safer to let pipeline complete when possible

---

### Page 3: Results & Reports

#### Results List

All completed pipeline runs are listed with:

- **Run ID:** Timestamped identifier (e.g., `run_20251027_143022`)
- **Date & Time:** When the run completed
- **Status:** Success, Error, or Stopped
- **Actions:**
  - `View Report`: Opens interactive HTML report
  - `Download`: Get all result files as ZIP

#### Interactive Report Viewer

Clicking "View Report" opens the HTML report in an embedded iframe:

**Report Sections:**

1. **Summary Dashboard**
   - Total genomes analyzed
   - Genes examined
   - Mutations identified
   - Key statistics

2. **Mutation Tables**
   - Sortable, filterable tables
   - Columns: Position, Reference, Alternate, Frequency, Effect
   - Click column headers to sort
   - Search box to filter rows

3. **Visualizations**
   - Mutation frequency heatmaps
   - Co-occurrence network graphs
   - Geographic/temporal distribution (if metadata provided)
   - Interactive charts (hover for details)

4. **Genome Details**
   - Per-genome mutation profiles
   - Resistance gene presence/absence
   - Sequence quality metrics

5. **Gene-Specific Analysis**
   - Separate section for each target gene
   - Alignment visualizations
   - Protein domain mapping
   - Known resistance mutations highlighted

**Viewer Controls:**

- **Fullscreen Toggle:** Expand report to full browser window
- **Download Report:** Save HTML file locally
- **Download CSV:** Export tables as spreadsheet
- **Download Logs:** Get pipeline execution logs
- **Print:** Browser print function (Ctrl+P)

#### File Downloads

**Available Downloads:**

1. **HTML Report:** Complete interactive report
   - Filename: `mutation_report.html`
   - Opens in any web browser
   - Self-contained (includes all CSS/JS)

2. **CSV Data Files:**
   - `mutations_summary.csv`: All mutations across genomes
   - `cooccurrence_matrix.csv`: Mutation pair frequencies
   - `gene_presence.csv`: Which genes found in each genome

3. **Alignment Files:**
   - Located in `alignments/` subdirectory
   - FASTA format alignments
   - One file per gene per genome

4. **Log Files:**
   - `pipeline.log`: Complete execution log
   - Useful for troubleshooting errors
   - Includes timestamps and detailed messages

#### Managing Results

- Results are saved in timestamped directories under `results/`
- Each run is independent (won't overwrite previous results)
- Manually delete old results to free disk space
- Export important results before cleaning

---

## Features & Interface

### Drag-and-Drop Upload

**How it works:**
1. Prepare your file (accessions.txt or genes.txt)
2. Drag file from file explorer
3. Hover over upload zone (turns darker blue)
4. Release mouse button to drop
5. File automatically uploads and validates

**Visual Feedback:**
- Gray → Blue: Hover effect
- Green checkmark: Upload successful
- Red X: Upload failed (see error message)
- Preview: First few lines displayed for verification

### Real-Time Progress Tracking

**Server-Sent Events (SSE) Technology:**
- Live updates without page refresh
- Sub-second latency
- Automatic reconnection if connection drops
- No polling overhead

**What Updates in Real-Time:**
- Stage status changes (pending → running → completed)
- Progress bars (0% → 100%)
- Log messages (appended as they occur)
- Current genome/gene being processed
- Estimated time remaining

### Responsive Design

**Works on Different Screen Sizes:**
- **Desktop (1920×1080+):** Full multi-column layout
- **Laptop (1366×768):** Optimized spacing
- **Tablet (768×1024):** Stacked columns, touch-friendly
- **Mobile (375×667+):** Single column, enlarged buttons

**Accessibility Features:**
- High contrast colors for readability
- Keyboard navigation support (Tab key)
- Screen reader compatible (ARIA labels)
- Font size adjustable via browser zoom

### Template System

**Purpose:** Simplify setup for common use cases

**How Templates Work:**
1. Pre-defined gene lists curated by AMR experts
2. Optimal database selection for gene type
3. One-click configuration
4. Still fully customizable after selection

**Available Templates:**

| Template | Genes | Use Case |
|----------|-------|----------|
| Quick Start | gyrA, parC, blaKPC | General demo/testing |
| Fluoroquinolone | gyrA, parC, gyrB, parE | Quinolone resistance |
| Beta-Lactamase | blaKPC, blaNDM, blaOXA-48, blaCTX-M | Carbapenem resistance |

**Custom Templates:**
- Users can request new templates via GitHub issues
- Copy template structure for in-house gene sets

### Health Monitoring

**Checks Performed:**

1. **Python Version:** Requires 3.8+
2. **Required Packages:** BioPython, pandas, Flask
3. **Optional Packages:** psutil (for memory stats)
4. **EMBOSS:** `needle` command availability
5. **ABRicate:** Installation and database setup
6. **Disk Space:** 10+ GB recommended
7. **Memory:** 2+ GB available
8. **Internet:** NCBI API reachability

**Health Status Levels:**

- **Good (Green):** All systems operational
- **Warning (Yellow):** Minor issues, pipeline may run with limitations
- **Error (Red):** Critical failure, must fix before proceeding

**Command-Line Health Check:**

```bash
# Human-readable output
python subscan/tools/health_check.py

# Detailed diagnostics
python subscan/tools/health_check.py --verbose

# JSON format (for programmatic use)
python subscan/tools/health_check.py --json
```

---

## Templates & Presets

### Creating Custom Templates

To add your own template to the GUI:

1. **Edit `subscan/gui/app.py`**

   Locate the `get_templates()` function and add your template:

   ```python
   {
       'id': 'my-custom-template',
       'name': 'My Custom Analysis',
       'description': 'Analysis of custom AMR genes',
       'genes': ['geneA', 'geneB', 'geneC'],
       'database': 'ncbi'
   }
   ```

2. **Restart the GUI**

   ```bash
   # Stop current server (Ctrl+C)
   # Restart
   python subscan/gui/launch.py
   ```

3. **Verify**

   Template should appear on the setup page

### Recommended Gene Sets

**For Different Organisms:**

- **E. coli:** gyrA, parC, blaKPC, blaCTX-M, mcr-1
- **K. pneumoniae:** blaKPC, blaNDM, blaOXA-48, rmtB
- **S. aureus:** mecA, vanA, ermA, fusA
- **P. aeruginosa:** oprD, mexB, blaVIM, blaPDC
- **M. tuberculosis:** rpoB, katG, inhA, embB

**For Different Resistance Classes:**

- **Fluoroquinolones:** gyrA, parC, gyrB, parE
- **Beta-lactams:** blaKPC, blaNDM, blaOXA, blaCTX-M, blaTEM, blaSHV
- **Aminoglycosides:** aac(6')-Ib, aph(3')-Ia, ant(2'')-Ia
- **Tetracyclines:** tetA, tetB, tetC, tetM
- **Macrolides:** ermA, ermB, mphA, msrA

---

## Troubleshooting

### Common Issues & Solutions

#### 1. GUI Won't Start

**Symptom:** `python subscan/gui/launch.py` shows error

**Possible Causes:**

- Flask not installed
- Port 5000 already in use
- Python version too old

**Solutions:**

```bash
# Install Flask
pip install flask

# Use different port
python subscan/gui/launch.py --port 8080

# Check Python version (need 3.8+)
python --version
```

#### 2. Browser Doesn't Auto-Open

**Symptom:** GUI starts but browser doesn't open

**Solutions:**

```bash
# Disable auto-open and manually navigate
python subscan/gui/launch.py --no-browser
# Then open http://localhost:5000 manually

# Check if running on remote server
# If yes, use SSH port forwarding:
ssh -L 5000:localhost:5000 user@server
# Then access http://localhost:5000 on your local machine
```

#### 3. File Upload Fails

**Symptom:** Red X appears after dropping file

**Possible Causes:**

- File too large (>10 MB)
- Wrong file format (must be .txt or .csv)
- File permissions issue

**Solutions:**

```bash
# Check file size
ls -lh accessions.txt

# Verify format (should be plain text)
file accessions.txt

# Fix permissions
chmod 644 accessions.txt

# Try smaller file first (split large files)
head -n 10 large_accessions.txt > test_accessions.txt
```

#### 4. Health Check Shows Errors

**Symptom:** Red X on EMBOSS or ABRicate

**Solutions:**

```bash
# For EMBOSS
# Ubuntu/WSL:
sudo apt update
sudo apt install emboss

# macOS:
brew install emboss

# Verify installation
needle -version

# For ABRicate
# Using conda (recommended):
conda install -c bioconda abricate

# Or using homebrew (macOS):
brew install abricate

# Verify installation
abricate --version

# Update ABRicate databases
abricate --setupdb
```

#### 5. Pipeline Fails at Stage 1 (Harvester)

**Symptom:** Error during genome download

**Possible Causes:**

- Invalid accession numbers
- NCBI connection issue
- Email not provided
- Rate limiting by NCBI

**Solutions:**

```bash
# Verify accession format
# Correct: NZ_CP012345.1, NC_000913.3
# Incorrect: CP012345 (missing NZ_ prefix)

# Test NCBI connectivity
curl -I https://www.ncbi.nlm.nih.gov

# Check if email is set in config
# Always provide a valid email address

# Reduce batch size if rate limited
# Split large accession lists into smaller batches
```

#### 6. Pipeline Fails at Stage 2 (Annotator)

**Symptom:** ABRicate annotation errors

**Solutions:**

```bash
# Update ABRicate databases
abricate --setupdb

# Check available databases
abricate --list

# Test ABRicate manually
abricate --db ncbi test_genome.fasta

# Verify genome files exist
ls genomes/
```

#### 7. Pipeline Fails at Stage 4 (Aligner)

**Symptom:** EMBOSS needle errors

**Solutions:**

```bash
# Verify EMBOSS installation
needle -version

# Test needle manually
needle -asequence seq1.fasta -bsequence seq2.fasta -outfile test.align

# Check if sequences are valid FASTA
grep ">" extracted_sequences/*.fasta
```

#### 8. Out of Memory

**Symptom:** Pipeline crashes, system freezes

**Solutions:**

```bash
# Reduce number of threads
# In GUI: Advanced Options → Threads: 2

# Process fewer genomes at once
# Split accession list into batches

# Check available memory
free -h  # Linux
vm_stat  # macOS

# Close other applications
# Increase swap space (Linux):
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

#### 9. Disk Space Error

**Symptom:** "No space left on device"

**Solutions:**

```bash
# Check disk usage
df -h

# Clean old results
rm -rf results/run_20251001_*

# Clean temporary files
rm -rf /tmp/mutationscan_*

# Move large result files to external storage
mv results/run_20251027_* /external_drive/
```

#### 10. Results Not Appearing

**Symptom:** Pipeline completes but no results shown

**Solutions:**

```bash
# Check results directory
ls -la results/

# Verify report file exists
find results/ -name "mutation_report.html"

# Check pipeline logs
tail -n 50 results/run_*/pipeline.log

# Manually open report
firefox results/run_20251027_*/mutation_report.html

# Restart GUI
# Sometimes cache needs clearing
```

### Getting Help

If issues persist:

1. **Check Logs:**
   - Browser Console (F12 → Console tab)
   - GUI terminal output
   - Pipeline log files in `results/run_*/`

2. **Run Health Check:**
   ```bash
   python subscan/tools/health_check.py --verbose
   ```

3. **Test Components Individually:**
   ```bash
   # Test genome download
   python subscan/tools/run_harvester.py --help
   
   # Test annotation
   abricate --version
   
   # Test alignment
   needle -version
   ```

4. **Report Issue:**
   - Open issue on GitHub: https://github.com/vihaankulkarni29/MutationScan/issues
   - Include: Error messages, health check output, pipeline logs
   - Specify: OS version, Python version, installation method

---

## FAQ

### General Questions

**Q: Do I need programming experience to use the GUI?**  
A: No! The GUI is designed for researchers with no coding background. Just drag, drop, and click.

**Q: Can I run multiple analyses simultaneously?**  
A: Currently, only one pipeline can run at a time per GUI instance. Start multiple instances on different ports if needed.

**Q: How long does a typical analysis take?**  
A: Depends on dataset size:
- 10 genomes, 5 genes: 5-10 minutes
- 100 genomes, 10 genes: 30-60 minutes
- 1000 genomes, 20 genes: 2-4 hours

**Q: Are my data uploaded to any server?**  
A: No! Everything runs locally on your computer. Only NCBI accession downloads require internet.

**Q: Can I use genomes I already downloaded?**  
A: Yes, but currently you need to use the CLI. GUI support for local genomes is planned for future versions.

### Technical Questions

**Q: What browsers are supported?**  
A: Chrome, Firefox, Edge, Safari (latest versions). Chrome recommended for best performance.

**Q: Can I run this on a high-performance computing (HPC) cluster?**  
A: Yes, but you'll need SSH port forwarding to access the GUI from your local machine.

**Q: Is there an API for programmatic access?**  
A: The GUI has REST endpoints (`/api/*`) documented in `subscan/gui/app.py`. For full CLI automation, use `run_wizard.py`.

**Q: Can I customize the HTML report template?**  
A: Yes! Edit `subscan/src/subscan/report_template.html`.

**Q: How do I cite MutationScan in publications?**  
A: Citation information is in the main README.md. Include the GitHub repository URL.

### Data & Results

**Q: What file formats are accepted?**  
A: Plain text (.txt) or CSV (.csv) for accessions and gene lists.

**Q: Can I analyze viruses or fungi?**  
A: The pipeline is optimized for bacteria. Viral/fungal genomes may work but haven't been extensively tested.

**Q: How are reference sequences selected?**  
A: First occurrence of each gene in your dataset is used as reference, or you can provide custom references.

**Q: What databases does ABRicate search?**  
A: NCBI, CARD, ResFinder, ARG-ANNOT, and others. See `abricate --list` for full list.

**Q: Can I export results to Excel?**  
A: CSV files can be opened directly in Excel. For formatted reports, copy tables from HTML report.

### Troubleshooting

**Q: Why is the progress bar stuck at 0%?**  
A: Some stages (like downloading) may appear stuck initially. Check the log viewer for activity. If truly stuck (no logs for >5 min), stop and restart.

**Q: What if I get "Connection lost" message?**  
A: The GUI server may have crashed. Check the terminal for errors, then restart with `python subscan/gui/launch.py`.

**Q: Can I pause and resume a pipeline?**  
A: No direct pause feature yet. Stopping requires restarting from beginning. Plan analyses to complete in one session.

**Q: Why are some mutations marked "Unknown effect"?**  
A: The mutation wasn't in reference databases. Doesn't mean it's not important - may require literature review.

---

## Best Practices

### Before Starting

1. ✅ Run health check first
2. ✅ Start with small test dataset (5-10 genomes)
3. ✅ Verify accession numbers are correct
4. ✅ Ensure sufficient disk space (estimate: 100 MB per genome)
5. ✅ Use a template if analyzing common resistance genes

### During Analysis

1. ✅ Monitor the log viewer for errors
2. ✅ Don't close browser window (can minimize)
3. ✅ Avoid putting computer to sleep mid-pipeline
4. ✅ Keep internet connection stable (for downloads)

### After Completion

1. ✅ Download and save HTML report externally
2. ✅ Review mutation tables for quality control
3. ✅ Export CSV files for further statistical analysis
4. ✅ Clean old results to free disk space
5. ✅ Document parameters used for reproducibility

---

## Advanced Usage

### Running on Remote Server

```bash
# On server:
python subscan/gui/launch.py --host 0.0.0.0 --port 5000 --no-browser

# On local machine:
ssh -L 5000:localhost:5000 user@server.address
# Then navigate to http://localhost:5000
```

### Custom Port

```bash
# Use port 8080 instead of default 5000
python subscan/gui/launch.py --port 8080
```

### Debug Mode

```bash
# Enable Flask debug mode (shows detailed errors)
python subscan/gui/launch.py --debug
```

### Integration with Other Tools

The GUI generates standard file formats compatible with:

- **Excel/R/Python:** CSV outputs
- **Jalview/MEGA:** FASTA alignments
- **Cytoscape:** Co-occurrence networks (CSV matrix)
- **Tableau/PowerBI:** Mutation tables (CSV)

---

## Keyboard Shortcuts

- **Ctrl+F:** Search in page (find in logs)
- **F5:** Refresh page (if needed)
- **Ctrl+P:** Print report (when viewing results)
- **Ctrl+S:** Save report HTML (when viewing)
- **Esc:** Exit fullscreen mode
- **Tab:** Navigate between form fields

---

## Version History

**v1.0 (October 2025)**
- Initial GUI release
- 7-stage pipeline visualization
- Drag-and-drop file upload
- Real-time progress tracking
- Template system
- Interactive HTML reports

---

## Additional Resources

- **Main Documentation:** `README.md`
- **CLI Wizard Guide:** `CLI_WIZARD_GUIDE.md`
- **Contributing:** `CONTRIBUTING.md`
- **GitHub Issues:** https://github.com/vihaankulkarni29/MutationScan/issues
- **Example Files:** `examples/` directory

---

**Need more help?** Open an issue on GitHub or consult the CLI Wizard Guide for terminal-based usage.
