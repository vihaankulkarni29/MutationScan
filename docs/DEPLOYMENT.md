# MutationScan - Production Deployment Guide

**Status:** ✅ PRODUCTION READY (February 3, 2026)

## Overview

MutationScan is a clinical-grade end-to-end bacterial genome analysis pipeline:

```
NCBI Genomes → ABRicate Gene Finding → Sequence Extraction → 
Variant Calling (with ML) → 3D Visualization
```

**Key Features:**
- 6 integrated Python modules (GenomeExtractor, GeneFinder, SequenceExtractor, VariantCaller, PyMOLVisualizer, ML Predictor)
- Trained ML models for novel mutation prediction (82% accuracy on holdout set)
- Windows-compatible (pathlib.Path, Docker)
- One-click deployment via Docker Compose
- CI/CD pipeline ready (GitHub Actions)

---

## Quick Start (Docker)

### Option 1: Docker Run (Linux/Mac)

```bash
# Build image
docker build -t mutationscan:v1 .

# Run pipeline
docker run -v $(pwd)/data:/app/data mutationscan:v1 \
  --email your.email@example.com \
  --query "Escherichia coli" \
  --limit 5 \
  --visualize
```

### Option 2: Docker Compose (Windows PowerShell)

```powershell
# Set email environment variable
$env:NCBI_EMAIL = "your.email@example.com"

# Run via docker-compose
docker-compose up

# OR with custom parameters
docker-compose run --rm mutationscan --email your.email@example.com --query "Klebsiella pneumoniae" --limit 2
```

### Option 3: GitHub Release (Auto-Build)

```bash
# Tag a release
git tag v1.0
git push origin v1.0

# GitHub Actions automatically:
# 1. Builds Docker image
# 2. Creates GitHub Release
# 3. Tagged as mutationscan:v1.0 and mutationscan:latest
```

---

## Command-Line Options

```
python src/main.py --email YOUR_EMAIL [OPTIONS]

Required:
  --email EMAIL              NCBI email (required by NCBI policy)

Optional:
  --query QUERY              Search query (default: "Escherichia coli")
  --limit LIMIT              Max genomes to download (default: 10)
  --api-key API_KEY          NCBI API key (for faster downloads)
  --visualize                Generate 3D structure visualizations
  --no-ml                    Disable ML predictions for novel mutations
  --antibiotic ANTIBIOTIC    Target antibiotic (default: Ciprofloxacin)
```

---

## Outputs

### Successful Run Generates:

```
data/
├── genomes/                    # Downloaded FASTA files
├── proteins/                   # Extracted protein sequences
├── results/
│   ├── mutation_report.csv     # Key output: Mutations with predictions
│   └── visualizations/         # PNG images of structures (if --visualize)
└── logs/
    └── pipeline.log            # Execution log (DEBUG level)
```

### mutation_report.csv Columns:

| Column | Meaning |
|--------|---------|
| Accession | Genome accession ID |
| Gene | Gene name (e.g., gyrA) |
| Mutation | Amino acid change (e.g., S83L) |
| Status | Resistant / Predicted High Risk / VUS |
| Phenotype | Clinical interpretation |
| prediction_score | 0.0-1.0 (confidence) |
| prediction_source | Clinical DB or AI Model |

---

## Docker Image Details

**Name:** `mutationscan:v1`  
**Size:** ~2.94 GB  
**Base:** staphb/abricate:latest (ABRicate + BLAST+ + Python)

**Includes:**
- ✅ ABRicate + CARD/NCBI databases (pre-downloaded, offline-ready)
- ✅ BLAST+ (BLASTn for sequence alignment)
- ✅ PyMOL (3D structure visualization)
- ✅ Python 3 + ML stack (scikit-learn, joblib, pandas)
- ✅ Pre-trained ML models (Ciprofloxacin resistance predictor)
- ✅ All MutationScan modules (src/)

**Pre-installed Python Packages:**
- biopython, pandas, numpy
- scikit-learn, joblib (ML engine)
- matplotlib, seaborn (visualization)
- ncbi-datasets-pylib (NCBI access)

---

## Architecture

### 5-Step Pipeline (Orchestrated by src/main.py)

```
Step 1: Download Genomes
├── Input: NCBI email, search query
├── Tool: NCBIDatasetsGenomeDownloader
└── Output: data/genomes/*.fasta

Step 2: Find Resistance Genes
├── Input: data/genomes/
├── Tool: ABRicate (CARD database)
└── Output: Gene coordinates (DataFrame)

Step 3: Extract & Translate Sequences
├── Input: Gene coordinates + genomes
├── Tool: SequenceExtractor
└── Output: data/proteins/*.faa (translated)

Step 4: Call Variants & Predict Resistance
├── Input: data/proteins/ + data/refs/
├── Tools: VariantCaller + ML Predictor
├── ML: Lazy loads trained models from models/
└── Output: data/results/mutation_report.csv

Step 5: Visualize (Optional)
├── Input: mutation_report.csv
├── Tool: PyMOLVisualizer
├── Condition: Only if --visualize flag
└── Output: data/results/visualizations/*.png
```

### Error Handling

Each step wrapped in try/except:
- **Step 1 fails** → Pipeline stops (no genomes, no point continuing)
- **Step 2 fails** → Pipeline stops (no genes found)
- **Step 3 fails** → Pipeline stops (no sequences extracted)
- **Step 4 fails** → Pipeline stops (variant calling critical)
- **Step 5 fails** → Warning only (visualization optional)

---

## Dependency Checks

On startup, the pipeline verifies:

| Tool | Status | Failure Mode |
|------|--------|--------------|
| ABRicate | CRITICAL | Fails with Docker/WSL instructions |
| BLAST+ | WARNING | Uses ABRicate alone if missing |
| PyMOL | WARNING | Skips visualization if missing |

Example error on Windows (without Docker/WSL):
```
[CRITICAL] ABRicate not found.

On Windows, ABRicate is typically not available natively.

SOLUTION: Run MutationScan via Docker:
  docker build -t mutationscan .
  docker run -v $(pwd)/data:/app/data mutationscan --email YOUR_EMAIL --query "E. coli"
```

---

## CI/CD Pipeline (GitHub Actions)

### Automatic Build on Release Tag

```bash
# Create release
git tag v1.0 -m "Production Release v1.0"
git push origin v1.0

# GitHub Actions automatically:
# 1. Checks out code
# 2. Builds Dockerfile with Docker Buildx
# 3. Creates GitHub Release with build info
# 4. Tags image as mutationscan:v1.0
```

**Workflow Configuration:** `.github/workflows/publish.yml`

To enable Docker Hub push:
1. Add Docker Hub credentials as GitHub Secrets
2. Set `push: true` in workflow YAML
3. Configure image registry tags

---

## ML Integration (Module 6)

### Prediction Accuracy

| Metric | Value |
|--------|-------|
| Novel Mutation Accuracy | 82.0% ± 11.7% |
| ROC-AUC Score | 89.4% ± 9.1% |
| Training Strategy | Leave-One-Mutation-Out CV |
| Model Type | RandomForest (5 features) |
| Features | Biophysical (hydrophobicity, charge, MW, etc.) |

### How It Works

1. **Unknown mutations** (not in resistance DB) → Route to ML
2. **ML loads models** from `models/ciprofloxacin_predictor.pkl`
3. **Encodes mutation** into 5 biophysical features
4. **Predicts resistance** probability (0.0-1.0)
5. **Outputs results** with `prediction_source = "AI Model"`

### Disable ML (Database-Only Mode)

```bash
docker run mutationscan:v1 --email user@email.com --no-ml
```

---

## Troubleshooting

### "ABRicate not found on Windows"

**Solution:** Use Docker (no native Windows support for ABRicate)

```powershell
docker-compose up
```

### "No genomes found for query"

**Check:**
1. Verify NCBI email is correct: `--email YOUR_EMAIL@example.com`
2. Try simpler query: `--query "E. coli"`
3. Try increasing limit: `--limit 20`
4. Check internet connection

### "PyMOL not found (visualization skipped)"

**Solution (Optional):** Install PyMOL locally

```bash
conda install -c schrodinger pymol-bundle
```

Or just skip: pipeline continues without visualization

### "No mutations detected"

**Possible causes:**
1. Wild-type reference sequences missing in `data/refs/`
2. Sequence query doesn't align well to references
3. All sequences are perfect matches to reference

---

## Production Deployment

### Recommended Setup

1. **Use Docker image** (standardized environment)
2. **Mount data volume** (preserve results)
3. **Set up logging** (pipeline.log for debugging)
4. **Use NCBI API key** (faster downloads)
5. **Run with timeout** (set max execution time)

### Example Kubernetes Deployment

```yaml
apiVersion: batch/v1
kind: Job
metadata:
  name: mutationscan-analysis
spec:
  template:
    spec:
      containers:
      - name: mutationscan
        image: mutationscan:v1.0
        command:
        - python3
        - src/main.py
        args:
        - --email
        - $(NCBI_EMAIL)
        - --query
        - "Klebsiella pneumoniae"
        - --limit
        - "10"
        - --visualize
        env:
        - name: NCBI_EMAIL
          valueFrom:
            secretKeyRef:
              name: ncbi-credentials
              key: email
        volumeMounts:
        - name: data
          mountPath: /app/data
      volumes:
      - name: data
        persistentVolumeClaim:
          claimName: mutationscan-data
      restartPolicy: Never
```

---

## Performance

| Operation | Time |
|-----------|------|
| Docker build | ~10-15 min (first time) |
| Download 10 genomes | ~2-5 min |
| Find genes (10 genomes) | ~1-2 min |
| Extract sequences | ~30 sec |
| Call variants | ~1-2 min |
| Generate visualizations | ~2-3 min |
| **Total (10 genomes)** | **~10-15 minutes** |

---

## Support & Citation

**Repository:** https://github.com/vihaankulkarni29/MutationScan  
**Documentation:** See DEVELOPMENT_LOG.md for detailed module docs  
**Issues:** Report via GitHub Issues

---

**MutationScan v2.0 - Production Ready**  
**Last Updated:** February 3, 2026
