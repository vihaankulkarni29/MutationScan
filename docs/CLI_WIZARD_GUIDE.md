# MutationScan CLI Wizard - User Guide

**Version 1.0** | **Last Updated: October 2025**

---

## Table of Contents

1. [Overview](#overview)
2. [When to Use the CLI Wizard](#when-to-use-the-cli-wizard)
3. [Installation & Requirements](#installation--requirements)
4. [Quick Start](#quick-start)
5. [Step-by-Step Walkthrough](#step-by-step-walkthrough)
6. [Gene Templates](#gene-templates)
7. [Command Files](#command-files)
8. [Advanced Usage](#advanced-usage)
9. [Troubleshooting](#troubleshooting)
10. [Examples](#examples)

---

## Overview

The **MutationScan CLI Wizard** is an interactive command-line tool that guides you through configuring the AMR mutation analysis pipeline without requiring extensive command-line expertise.

**What it does:**
- Asks simple questions to collect pipeline parameters
- Provides pre-configured gene templates for common analyses
- Validates your inputs in real-time
- Generates a complete pipeline command
- Saves the command for reproducibility
- Optionally runs the pipeline immediately

**Why use the CLI Wizard over the GUI?**
- Remote server access (no GUI display needed)
- Scriptable/automatable workflows
- Lower resource usage (no web browser)
- Better for HPC/cluster environments
- Preferred by command-line power users

---

## When to Use the CLI Wizard

### ✅ Use CLI Wizard When:

- Working on remote server via SSH
- Running on HPC cluster without X11 forwarding
- Automating analyses via shell scripts
- Resources constrained (no GUI overhead)
- Prefer keyboard-driven workflows
- Need to document exact commands used
- Running headless/server systems

### ❌ Use Web GUI Instead When:

- New to command-line tools
- Want visual progress tracking
- Prefer drag-and-drop file upload
- Need to monitor multiple analyses
- Want interactive report viewing
- Working on local desktop/laptop with display

---

## Installation & Requirements

### Prerequisites

Same as main pipeline:

- **Python 3.8+**
- **EMBOSS** (for alignment)
- **ABRicate** (for annotation)
- **Internet connection** (for NCBI downloads)

### Verification

```bash
# Check Python version
python --version  # Should be 3.8 or higher

# Verify EMBOSS
needle -version

# Verify ABRicate
abricate --version

# Run health check
python subscan/tools/health_check.py
```

### No Additional Installation

The wizard uses only Python standard library, so if your MutationScan installation works, the wizard is ready to use!

---

## Quick Start

### 30-Second Demo

```bash
# Navigate to MutationScan directory
cd MutationScan

# Launch wizard
python subscan/tools/run_wizard.py

# Follow the prompts:
# 1. Upload or create accession file
# 2. Select gene template (e.g., "Quick Start")
# 3. Enter your email
# 4. Choose species/reference option
# 5. Select database (default: NCBI)
# 6. Configure output directory
# 7. Review and run!
```

**What happens:**
1. Wizard asks 6 simple questions
2. Shows color-coded configuration summary
3. Generates complete pipeline command
4. Saves command to timestamped `.sh` file
5. Asks if you want to run now or later

---

## Step-by-Step Walkthrough

### Launching the Wizard

```bash
python subscan/tools/run_wizard.py
```

You'll see a colorful banner:

```
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║         🧬  MUTATION SCAN - INTERACTIVE CLI WIZARD  🧬          ║
║                                                                  ║
║        AMR Mutation Analysis Pipeline Configuration             ║
║                      Version 1.0                                 ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝
```

---

### Step 1: Accession IDs

**Question:** "How would you like to provide accession IDs?"

**Options:**

1. **Upload existing file** - You already have a `.txt` file with accessions
2. **Create new file from input** - Type accessions directly in terminal

#### Option 1: Upload Existing File

```
Enter path to accessions file: examples/demo_accessions.txt
```

**Tips:**
- Use Tab completion for file paths
- Relative or absolute paths both work
- File must exist and be readable

#### Option 2: Create from Input

```
Enter accession IDs (one per line, Ctrl+D or Ctrl+Z when done):
NZ_CP012345.1
NZ_CP067890.1
NC_000913.3
^D  # Press Ctrl+D (Linux/Mac) or Ctrl+Z (Windows)
```

**Tips:**
- Paste multiple lines at once
- Blank lines are ignored
- Press Ctrl+D (Unix) or Ctrl+Z (Windows) to finish

**Validation:**
- Wizard checks file exists and is readable
- Shows preview of first few accessions
- Asks for confirmation

---

### Step 2: Target Genes

**Question:** "How would you like to specify target genes?"

**Options:**

1. **Use a common gene template** - Pre-configured gene lists
2. **Upload existing genes file** - Your own `.txt` file
3. **Create new file from input** - Type genes in terminal

#### Option 1: Gene Templates (Recommended)

```
Available gene templates:
  1. Fluoroquinolone Resistance
     Genes: gyrA, parC, gyrB, parE
     
  2. Beta-Lactamase
     Genes: blaKPC, blaNDM, blaOXA-48, blaCTX-M, blaTEM, blaSHV
     
  3. Aminoglycoside Resistance
     Genes: aac(6')-Ib, aph(3')-Ia, ant(2'')-Ia, rmtB
     
  4. Tetracycline Resistance
     Genes: tetA, tetB, tetC, tetM
     
  5. Common AMR Mix
     Genes: gyrA, parC, blaKPC, blaNDM, aac(6')-Ib, tetA

Select a template (1-5): 2
```

**See [Gene Templates](#gene-templates) section for full details.**

#### Option 2: Upload File

```
Enter path to genes file: data_input/gene_list.txt
```

#### Option 3: Create from Input

```
Enter gene names (one per line, Ctrl+D or Ctrl+Z when done):
gyrA
parC
blaKPC
^D
```

---

### Step 3: Email Address

**Question:** "Enter your email address (required for NCBI downloads):"

```
Email: researcher@university.edu
```

**Why needed:**
- NCBI requires email for API rate tracking
- Used only for NCBI's internal logging
- No spam or marketing emails

**Validation:**
- Must contain `@` symbol
- Must look like valid email format
- Not sent anywhere except NCBI API

---

### Step 4: Species/Reference Selection

**Question:** "Do you want to auto-download species references or provide your own?"

**Options:**

1. **Auto-download using sepi-species** - Automatic reference genome download
2. **Provide custom reference directory** - Use your own reference sequences

#### Option 1: Auto-Download (Recommended for Beginners)

```
Select option (1-2): 1
Enter species name: Escherichia coli
```

**Supported species:**
- Escherichia coli
- Klebsiella pneumoniae
- Staphylococcus aureus
- Pseudomonas aeruginosa
- Salmonella enterica
- Mycobacterium tuberculosis
- And many more...

**What happens:**
- Wizard will download reference genome for that species
- Creates reference directory automatically
- Uses latest complete genome from NCBI

#### Option 2: Custom Reference

```
Select option (1-2): 2
Enter path to reference directory: /path/to/my/references/
```

**Requirements:**
- Directory must exist
- Should contain reference FASTA files for your target genes
- One file per gene (e.g., `gyrA.fasta`, `parC.fasta`)

---

### Step 5: Database Selection

**Question:** "Select database source:"

**Options:**

1. **NCBI** (National Center for Biotechnology Information) - Default, recommended
2. **BV-BRC** (Bacterial and Viral Bioinformatics Resource Center)
3. **EnteroBase** (Specialized for Enterobacteriaceae)
4. **PATRIC** (Legacy database)

```
Select database (1-4) [default: 1]: 1
```

**Recommendations:**
- Use NCBI for most analyses (most comprehensive)
- Use EnteroBase for E. coli, Salmonella, Yersinia
- Use BV-BRC for specialized bacterial datasets
- Press Enter to accept default (NCBI)

---

### Step 6: Output & Advanced Options

#### Output Directory

```
Enter output directory [default: results]: my_analysis_output
```

**Tips:**
- Creates directory if it doesn't exist
- Timestamped subdirectory created automatically
- Press Enter for default `results/`

#### Number of Threads

```
Enter number of CPU threads [default: 4]: 8
```

**Guidelines:**
- More threads = faster processing
- Don't exceed your CPU core count
- Check with: `nproc` (Linux) or `sysctl -n hw.ncpu` (Mac)
- Reduce if system becomes unresponsive

#### Verbose Output

```
Enable verbose output? (Y/n) [default: yes]: y
```

**Verbose mode:**
- Shows detailed progress messages
- Useful for troubleshooting
- Slightly slower due to extra logging
- Recommended for first runs

#### Auto-Open Report

```
Open HTML report automatically when complete? (Y/n) [default: yes]: y
```

**If yes:**
- Report opens in default browser when pipeline finishes
- Requires GUI display (disable on servers)
- Only works on local machines

---

### Step 7: Configuration Summary

After answering all questions, wizard shows a colored summary:

```
╔══════════════════════════════════════════════════════════════════╗
║                  CONFIGURATION SUMMARY                           ║
╚══════════════════════════════════════════════════════════════════╝

📋 Accession File:    examples/demo_accessions.txt (10 accessions)
🧬 Genes File:        /tmp/mutationscan_genes_12345.txt (5 genes)
📧 Email:             researcher@university.edu
🔬 Species:           Escherichia coli (auto-download)
💾 Database:          NCBI
📁 Output Directory:  results/
🔧 Threads:           8
📊 Verbose:           Yes
🌐 Open Report:       Yes

✓ Configuration ready!

Do you want to proceed? (Y/n):
```

**Review carefully before proceeding!**

---

### Step 8: Command Generation & Execution

#### Generated Command

Wizard creates a complete pipeline command:

```bash
python subscan/tools/run_pipeline.py \
  --accessions examples/demo_accessions.txt \
  --genes /tmp/mutationscan_genes_12345.txt \
  --email researcher@university.edu \
  --species "Escherichia coli" \
  --database ncbi \
  --output results \
  --threads 8 \
  --verbose \
  --open-report
```

#### Saved Command File

Command saved to timestamped shell script:

```
✓ Command saved to: run_command_20251027_143522.sh
```

**File contents:**

```bash
#!/bin/bash
# MutationScan Pipeline Command
# Generated: 2025-10-27 14:35:22
# Configuration: Quick Start template

python subscan/tools/run_pipeline.py \
  --accessions examples/demo_accessions.txt \
  --genes /tmp/mutationscan_genes_12345.txt \
  --email researcher@university.edu \
  --species "Escherichia coli" \
  --database ncbi \
  --output results \
  --threads 8 \
  --verbose \
  --open-report
```

**Benefits:**
- Reproducibility: Re-run exact same analysis
- Documentation: Record of parameters used
- Sharing: Send to collaborators
- Version control: Commit to Git

#### Run Now or Later?

```
Do you want to run the pipeline now? (Y/n):
```

**Option 1: Run Now (Y)**
- Pipeline starts immediately
- Shows real-time output
- Blocks until completion
- Ctrl+C to interrupt

**Option 2: Run Later (n)**
- Command saved to file
- Run when convenient:
  ```bash
  bash run_command_20251027_143522.sh
  ```
- Can modify command file if needed

---

## Gene Templates

### 1. Fluoroquinolone Resistance

**Genes:** `gyrA`, `parC`, `gyrB`, `parE`

**Targets:**
- DNA gyrase subunits (gyrA, gyrB)
- Topoisomerase IV subunits (parC, parE)

**Mutations to watch:**
- gyrA: S83L, D87N (E. coli numbering)
- parC: S80I, E84K
- QRDR (Quinolone Resistance-Determining Region) variants

**Use cases:**
- Fluoroquinolone/quinolone resistance profiling
- Ciprofloxacin, levofloxacin resistance studies
- High-frequency mutation analysis

---

### 2. Beta-Lactamase

**Genes:** `blaKPC`, `blaNDM`, `blaOXA-48`, `blaCTX-M`, `blaTEM`, `blaSHV`

**Targets:**
- Carbapenemases: KPC, NDM, OXA-48
- ESBLs: CTX-M, TEM, SHV families

**Resistance profile:**
- KPC: Carbapenems (meropenem, imipenem)
- NDM: Broad-spectrum beta-lactams + carbapenems
- OXA-48: Carbapenems, weak ESBL activity
- CTX-M: Extended-spectrum cephalosporins

**Use cases:**
- Carbapenem resistance surveillance
- ESBL detection
- Hospital outbreak investigations

---

### 3. Aminoglycoside Resistance

**Genes:** `aac(6')-Ib`, `aph(3')-Ia`, `ant(2'')-Ia`, `rmtB`

**Targets:**
- Acetyltransferases: aac enzymes
- Phosphotransferases: aph enzymes
- Nucleotidyltransferases: ant enzymes
- 16S rRNA methyltransferases: rmtB

**Resistance profile:**
- aac(6')-Ib: Amikacin, tobramycin
- aph(3')-Ia: Kanamycin, neomycin
- rmtB: High-level pan-aminoglycoside resistance

**Use cases:**
- Aminoglycoside resistance mechanisms
- Co-resistance with beta-lactams (common)

---

### 4. Tetracycline Resistance

**Genes:** `tetA`, `tetB`, `tetC`, `tetM`

**Targets:**
- Efflux pumps: tetA, tetB, tetC
- Ribosomal protection: tetM

**Resistance profile:**
- tetA/B/C: Tetracycline, doxycycline (efflux-mediated)
- tetM: Tetracycline (ribosomal protection)

**Use cases:**
- Tetracycline resistance surveillance
- Mobile genetic element tracking (often on plasmids)

---

### 5. Common AMR Mix

**Genes:** `gyrA`, `parC`, `blaKPC`, `blaNDM`, `aac(6')-Ib`, `tetA`

**Rationale:**
- Representative genes from multiple resistance classes
- Good for initial/broad screening
- Covers most common resistance mechanisms

**Use cases:**
- First-pass AMR profiling
- Unknown resistance patterns
- Multi-drug resistance (MDR) analysis

---

## Command Files

### Structure

Generated `.sh` files are standard bash scripts:

```bash
#!/bin/bash
# MutationScan Pipeline Command
# Generated: 2025-10-27 14:35:22
# Configuration: Beta-Lactamase template

python subscan/tools/run_pipeline.py \
  --accessions data_input/accession_list.txt \
  --genes /tmp/mutationscan_genes_98765.txt \
  --email user@example.com \
  --species "Klebsiella pneumoniae" \
  --database ncbi \
  --output results \
  --threads 4 \
  --verbose \
  --open-report
```

### Running Saved Commands

```bash
# Method 1: Direct execution
bash run_command_20251027_143522.sh

# Method 2: Make executable
chmod +x run_command_20251027_143522.sh
./run_command_20251027_143522.sh

# Method 3: Source it (runs in current shell)
source run_command_20251027_143522.sh
```

### Modifying Commands

Edit the file before running:

```bash
# Open in text editor
nano run_command_20251027_143522.sh

# Common modifications:
# - Change --threads value
# - Update --output directory
# - Add/remove --verbose flag
# - Change database
```

### Batch Processing

Create multiple command files and run sequentially:

```bash
#!/bin/bash
# batch_analysis.sh

# Run analysis 1
bash run_command_001.sh

# Run analysis 2
bash run_command_002.sh

# Run analysis 3
bash run_command_003.sh
```

### Version Control

Track command files in Git:

```bash
# Add to repository
git add run_command_*.sh

# Commit with descriptive message
git commit -m "Add pipeline commands for E. coli quinolone study"

# Push to remote
git push origin main
```

---

## Advanced Usage

### Non-Interactive Mode (Future Feature)

For fully automated workflows, you can modify the wizard to accept environment variables:

```bash
export MUTATIONSCAN_ACCESSIONS="data_input/accessions.txt"
export MUTATIONSCAN_GENES="data_input/genes.txt"
export MUTATIONSCAN_EMAIL="auto@example.com"
# ... etc

python subscan/tools/run_wizard.py --auto
```

*Note: This feature is planned but not yet implemented.*

### Custom Gene Lists

Create your own gene list file:

```bash
# my_custom_genes.txt
rpoB
katG
inhA
embB
pncA
gyrA
gyrB
```

Then use in wizard:

```
Step 2 → Option 2 (Upload existing file) → my_custom_genes.txt
```

### Scripting with Wizard Output

Parse the generated command file:

```bash
# Extract parameters from command file
ACCESSIONS=$(grep -oP '(?<=--accessions )\S+' run_command_*.sh)
GENES=$(grep -oP '(?<=--genes )\S+' run_command_*.sh)

echo "Accessions file: $ACCESSIONS"
echo "Genes file: $GENES"
```

### Integration with HPC Job Schedulers

#### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=mutationscan
#SBATCH --output=mutationscan_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Load modules
module load python/3.9
module load emboss
module load abricate

# Run wizard-generated command
bash run_command_20251027_143522.sh
```

#### PBS Example

```bash
#!/bin/bash
#PBS -N mutationscan
#PBS -l nodes=1:ppn=16
#PBS -l mem=32gb
#PBS -l walltime=04:00:00
#PBS -j oe

cd $PBS_O_WORKDIR

# Load environment
source ~/miniconda3/bin/activate mutationscan_env

# Execute command
bash run_command_20251027_143522.sh
```

### Logging

Capture all output to log file:

```bash
# Run with full logging
bash run_command_20251027_143522.sh 2>&1 | tee analysis.log

# Run in background with logging
nohup bash run_command_20251027_143522.sh > analysis.log 2>&1 &

# Check progress
tail -f analysis.log
```

---

## Troubleshooting

### Wizard Won't Start

**Error:** `python: can't open file 'subscan/tools/run_wizard.py'`

**Solution:**

```bash
# Ensure you're in MutationScan root directory
cd /path/to/MutationScan

# Or use absolute path
python /full/path/to/MutationScan/subscan/tools/run_wizard.py
```

---

### File Not Found During Input

**Error:** `File not found: data_input/accessions.txt`

**Solutions:**

```bash
# Check if file exists
ls -la data_input/accessions.txt

# Use absolute path
/home/user/MutationScan/data_input/accessions.txt

# Or create the file first
cat > data_input/accessions.txt << EOF
NZ_CP012345.1
NC_000913.3
EOF
```

---

### Email Validation Fails

**Error:** `Invalid email format`

**Solution:**

- Must contain `@` symbol
- Examples of valid emails:
  - `user@university.edu`
  - `researcher@institution.org`
  - `name@company.com`

---

### Cannot Create Output Directory

**Error:** `Permission denied: /protected/directory`

**Solutions:**

```bash
# Use directory you have write access to
Enter output directory: ~/mutationscan_results

# Or create directory first with proper permissions
mkdir -p ~/mutationscan_results
chmod 755 ~/mutationscan_results
```

---

### Pipeline Fails After Generation

**Error:** Command runs but pipeline fails

**Diagnostic steps:**

```bash
# 1. Check health status
python subscan/tools/health_check.py --verbose

# 2. Verify files exist
cat data_input/accession_list.txt
cat /tmp/mutationscan_genes_*.txt

# 3. Test pipeline components individually
python subscan/tools/run_harvester.py --help

# 4. Check logs
tail -n 100 results/run_*/pipeline.log
```

---

### Keyboard Interrupt (Ctrl+C)

**What happens:**

- Wizard catches Ctrl+C gracefully
- Displays: "Setup cancelled by user"
- Temporary files cleaned up
- Safe to restart wizard

**To fully exit:**

- Press Ctrl+C once
- Wizard exits cleanly

---

### Colors Not Displaying

**Symptom:** See escape codes like `\033[92m` instead of colors

**Cause:** Terminal doesn't support ANSI colors

**Solutions:**

```bash
# Windows: Use Windows Terminal, PowerShell, or Git Bash
# (Not cmd.exe)

# Linux/Mac: Should work by default

# SSH: Ensure TERM is set
echo $TERM  # Should be xterm-256color or similar
export TERM=xterm-256color
```

---

### Gene Template Not Listed

**Problem:** Want a template that doesn't exist

**Solutions:**

1. **Use existing template as starting point**
   - Select closest match
   - Modify genes file afterward

2. **Request new template**
   - Open GitHub issue with gene list
   - We'll add popular templates

3. **Create custom template**
   - Edit `subscan/tools/run_wizard.py`
   - Add to `suggest_common_genes()` function
   - See [Contributing Guide](../CONTRIBUTING.md)

---

## Examples

### Example 1: Quick Fluoroquinolone Analysis

```bash
$ python subscan/tools/run_wizard.py

# Step 1: Accessions
How would you like to provide accession IDs?
  1. Upload existing file
  2. Create new file from input
Select (1-2): 1
Enter path: examples/demo_accessions.txt

# Step 2: Genes
How would you like to specify target genes?
  1. Use a common gene template
  2. Upload existing file
  3. Create new file from input
Select (1-3): 1

Available templates:
  1. Fluoroquinolone Resistance
  ...
Select template (1-5): 1

# Step 3: Email
Enter email: researcher@university.edu

# Step 4: Reference
Auto-download or custom?
  1. Auto-download using sepi-species
  2. Provide custom reference
Select (1-2): 1
Species name: Escherichia coli

# Step 5: Database
Select database:
  1. NCBI
  ...
Select (1-4) [default: 1]: [Enter]

# Step 6: Output
Output directory [default: results]: fluoroquinolone_study
Threads [default: 4]: 8
Verbose (Y/n): y
Open report (Y/n): y

# Review and run
Do you want to proceed? Y
Run now? Y
```

---

### Example 2: Custom Gene List

```bash
# Create custom gene file first
cat > my_genes.txt << EOF
mecA
vanA
ermA
tetM
EOF

# Run wizard
python subscan/tools/run_wizard.py

# Select "Upload existing file" for genes
# Point to my_genes.txt
# Continue with other steps...
```

---

### Example 3: Batch Analysis with Templates

Create a script:

```bash
#!/bin/bash
# analyze_all.sh - Run multiple template analyses

for template in "Fluoroquinolone" "Beta-Lactamase" "Tetracycline"
do
    echo "Analyzing with $template template..."
    
    # Run wizard non-interactively (manual for now)
    # Or modify command files and run sequentially
    
    # Wait for completion
    wait
    
    echo "$template analysis complete!"
done
```

---

### Example 4: Remote Server Workflow

```bash
# SSH to server
ssh user@research-server.edu

# Navigate to project
cd /data/projects/mutationscan

# Load environment
module load python/3.9 emboss abricate

# Run wizard
python subscan/tools/run_wizard.py

# Important: Disable auto-open report (no GUI on server)
Open report (Y/n): n

# Run in background with logging
Run now? n
nohup bash run_command_20251027_143522.sh > analysis.log 2>&1 &

# Monitor progress
tail -f analysis.log

# Download results to local machine (from local terminal)
scp -r user@server:/data/projects/mutationscan/results/run_* .
```

---

## Best Practices

### Before Running Wizard

1. ✅ Prepare input files beforehand
2. ✅ Run health check to verify dependencies
3. ✅ Estimate disk space needed (100 MB/genome)
4. ✅ Choose appropriate template or create gene list
5. ✅ Have NCBI email ready

### During Wizard Session

1. ✅ Read questions carefully
2. ✅ Use Tab completion for file paths
3. ✅ Review configuration summary before proceeding
4. ✅ Save command file even if running immediately
5. ✅ Use verbose mode for first runs

### After Command Generation

1. ✅ Backup command file (version control)
2. ✅ Document parameters in lab notebook
3. ✅ Run small test dataset first
4. ✅ Monitor first few stages for errors
5. ✅ Save results with descriptive names

---

## Keyboard Shortcuts

While in wizard:

- **Ctrl+C:** Cancel and exit
- **Ctrl+D:** Finish multi-line input (Unix/Mac)
- **Ctrl+Z:** Finish multi-line input (Windows)
- **Tab:** Autocomplete file paths (if shell supports)
- **↑ / ↓ arrows:** Navigate command history

---

## Comparison: Wizard vs GUI vs Manual CLI

| Feature | Wizard | Web GUI | Manual CLI |
|---------|--------|---------|------------|
| Ease of use | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ |
| Remote server | ✅ | ❌ | ✅ |
| Visual progress | ❌ | ✅ | ❌ |
| Resource usage | Low | Medium | Low |
| Reproducibility | ✅✅ | ⭐ | ✅✅ |
| Automation | ⭐ | ❌ | ✅✅ |
| Learning curve | Low | Very Low | High |
| HPC/cluster | ✅ | ❌ | ✅ |

**Legend:**
- ⭐ = Stars (rating out of 5)
- ✅ = Fully supported
- ✅✅ = Excellent support
- ⭐ = Partial support
- ❌ = Not supported

---

## Additional Resources

- **Web GUI Guide:** `docs/GUI_USER_GUIDE.md`
- **Main README:** `README.md`
- **Health Check Tool:** `subscan/tools/health_check.py --help`
- **Pipeline Help:** `python subscan/tools/run_pipeline.py --help`
- **GitHub Issues:** https://github.com/vihaankulkarni29/MutationScan/issues

---

## Version History

**v1.0 (October 2025)**
- Initial CLI wizard release
- 6-step interactive configuration
- 5 pre-configured gene templates
- Input validation and error handling
- Command file generation
- Immediate execution option

---

## Contributing

Help improve the wizard:

- **Report bugs:** GitHub issues
- **Request templates:** Open an issue with your gene list
- **Submit code:** See `CONTRIBUTING.md`
- **Improve docs:** Pull requests welcome!

---

**Need help?** Check the troubleshooting section or open a GitHub issue with the `wizard` tag.
