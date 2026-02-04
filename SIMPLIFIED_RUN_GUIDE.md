# MutationScan - Simplified Run Guide

## Quick Start (One-Command Execution)

### Docker (Recommended)

**Windows PowerShell:**
```powershell
# Analyze local genome for specific genes
.\run_docker.ps1 -Email your@email.com -Genome data/genomes/GCF_000005845.2.fasta -Targets data/config/acrR_targets.txt

# Download organism and analyze
.\run_docker.ps1 -Email your@email.com -Organism "Escherichia coli" -Targets data/config/acrR_targets.txt -ApiKey YOUR_KEY
```

**Linux/Mac:**
```bash
# Analyze local genome for specific genes
./run_docker.sh --email your@email.com --genome data/genomes/GCF_000005845.2.fasta --targets data/config/acrR_targets.txt

# Download organism and analyze
./run_docker.sh --email your@email.com --organism "Escherichia coli" --targets data/config/acrR_targets.txt --api-key YOUR_KEY
```

### WSL (Windows Subsystem for Linux)

```bash
# Install ABRicate in WSL
wsl -e sudo apt update
wsl -e sudo apt install abricate ncbi-blast+ python3-pip

# Run analysis
wsl -e python3 src/main.py --email your@email.com --genome data/genomes/GCF_000005845.2.fasta --targets data/config/acrR_targets.txt
```

---

## Required Inputs

### 1. Email (Always Required)
NCBI requires an email for all API requests.

### 2. Genome Source (Choose One)
- **Local genome file** (`--genome`): Use if you already have a FASTA file
- **Organism query** (`--organism`): Download from NCBI

### 3. Target Genes (Optional but Recommended)
Create a text file with one gene per line:

**Example: `data/config/acrR_targets.txt`**
```
acrR
```

**Example: `data/config/mtb_targets.txt`**
```
katG
rpoB
embB
pncA
gyrA
```

If not provided, all resistance genes detected by ABRicate will be analyzed.

---

## Optional Parameters

- `--api-key YOUR_KEY`: NCBI API key for higher download limits (recommended)
- `--limit N`: Number of genomes to download when using `--organism` (default: 1)
- `--visualize`: Generate 3D PyMOL visualizations
- `--no-ml`: Disable machine learning predictions

---

## Output

Results are saved to `data/results/`:
- `mutation_report.csv`: All detected mutations with resistance phenotypes
- `visualizations/`: 3D structure images (if `--visualize` used)

---

## Complete Examples

### Research Run: E. coli acrR Analysis

**With local genome:**
```bash
# Docker (PowerShell)
.\run_docker.ps1 -Email user@example.com -Genome data/genomes/GCF_000005845.2.fasta -Targets data/config/acrR_targets.txt

# WSL
wsl -e python3 src/main.py --email user@example.com --genome data/genomes/GCF_000005845.2.fasta --targets data/config/acrR_targets.txt
```

**Download from NCBI:**
```bash
# Docker (PowerShell)
.\run_docker.ps1 -Email user@example.com -Organism "Escherichia coli" -Targets data/config/acrR_targets.txt -ApiKey YOUR_KEY -Limit 5

# WSL
wsl -e python3 src/main.py --email user@example.com --organism "Escherichia coli" --targets data/config/acrR_targets.txt --api-key YOUR_KEY --limit 5
```

### Research Run: M. tuberculosis

Create `data/config/mtb_targets.txt`:
```
katG
rpoB
embB
pncA
gyrA
```

Run analysis:
```bash
# Docker
.\run_docker.ps1 -Email user@example.com -Organism "Mycobacterium tuberculosis" -Targets data/config/mtb_targets.txt -ApiKey YOUR_KEY

# WSL
wsl -e python3 src/main.py --email user@example.com --organism "Mycobacterium tuberculosis" --targets data/config/mtb_targets.txt --api-key YOUR_KEY
```

---

## Troubleshooting

### "Email is required"
Add `--email` or `-Email` with your email address.

### "Either --genome or --organism is required"
Provide one of:
- `--genome path/to/file.fasta` (local genome)
- `--organism "Organism name"` (download from NCBI)

### "NCBI download failed"
- Add `--api-key` for better reliability
- Check your internet connection
- Verify organism name is correct

### "No genes found"
- Check that your genome file is valid FASTA format
- Try without `--targets` to see all detected genes first
- Verify target gene names match CARD database nomenclature

---

## Getting an NCBI API Key (Recommended)

1. Go to https://www.ncbi.nlm.nih.gov/account/
2. Sign in or create account
3. Go to Settings → API Key Management
4. Create a new API key
5. Use it with `--api-key YOUR_KEY`

Benefits: 10x higher request limits, better reliability
