# MutationScan Deployment Guide

## 🚀 Quick Deploy

### Windows Development Setup
```powershell
# Clone repository
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan

# Install Python dependencies
cd subscan
pip install -e .

# Verify installation
python check_structure.py
python check_dependencies.py

# Test with sample data
python tools/run_pipeline.py --accessions examples/demo_accessions.txt --gene-list examples/gene_list.txt --email you@example.com --output-dir test_output --sepi-species "Escherichia coli" --dry-run
```

### Linux Production Setup
```bash
# System dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y python3 python3-pip git abricate blast+ mafft

# Clone and install
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan/subscan
pip3 install -e .

# Production test
python3 tools/run_pipeline.py --accessions examples/demo_accessions.txt --gene-list examples/gene_list.txt --email your.email@lab.org --output-dir production_test --sepi-species "Escherichia coli"
```

---

## 🌍 Multi-Platform Architecture

### Development vs Production Environments

| Component | Windows (Dev) | Linux (Production) |
|-----------|---------------|-------------------|
| **Python** | 3.9+ | 3.9+ |
| **ABRicate** | Mocked (`mock_abricate.py`) | Native installation |
| **BLAST+** | Optional | Required |
| **MAFFT** | Optional | Required |
| **Genome Extraction** | Federated (all platforms) | Federated (all platforms) |
| **Performance** | Limited by mocks | Full bioinformatics suite |

### Cross-Platform Compatibility Matrix

| Feature | Windows | Linux | macOS | Docker |
|---------|---------|-------|--------|--------|
| Pipeline Orchestrator | ✅ | ✅ | ✅ | ✅ |
| Federated Extractor | ✅ | ✅ | ✅ | ✅ |
| AMR Gene Annotation | 🟡 Mock | ✅ Native | 🔶 Container | ✅ |
| Sequence Alignment | 🟡 Mock | ✅ Native | 🔶 Container | ✅ |
| Report Generation | ✅ | ✅ | ✅ | ✅ |

**Legend:** ✅ Full Support | 🟡 Development Mode | 🔶 Requires Container

---

## 🗃️ Database Configuration

### Federated Multi-Database Support

MutationScan's federated genome extractor supports 4 major microbial databases:

```python
# Database priority and fallback order
DATABASES = {
    "ncbi": {
        "name": "NCBI GenBank",
        "api": "https://eutils.ncbi.nlm.nih.gov/",
        "coverage": "Global reference genomes",
        "specialization": "Curated assemblies"
    },
    "bvbrc": {
        "name": "BV-BRC", 
        "api": "https://www.bv-brc.org/api/",
        "coverage": "NIAID priority pathogens",
        "specialization": "Bacterial and viral"
    },
    "enterobase": {
        "name": "EnteroBase",
        "api": "https://enterobase.warwick.ac.uk/api/",
        "coverage": "Enterobacteriaceae",
        "specialization": "Salmonella, E. coli, Shigella"
    },
    "patric": {
        "name": "PATRIC",
        "api": "https://patricbrc.org/api/",
        "coverage": "Bacterial pathogens",
        "specialization": "Legacy genomics platform"
    }
}
```

### Database Selection Strategies

#### 1. **Multi-Database (Recommended)**
```bash
# Maximum genome coverage with automatic fallback
python tools/run_harvester.py --database all --accessions genomes.txt --email you@lab.org --output-dir results/
```

#### 2. **NCBI-Only (Legacy)**
```bash
# Traditional NCBI-focused approach
python tools/run_harvester.py --database ncbi --accessions genomes.txt --email you@lab.org --output-dir results/
```

#### 3. **Specialized Databases**
```bash
# E. coli/Salmonella research
python tools/run_harvester.py --database enterobase --accessions genomes.txt --email you@lab.org --output-dir results/

# NIAID priority pathogen research  
python tools/run_harvester.py --database bvbrc --accessions genomes.txt --email you@lab.org --output-dir results/
```

---

## 🔧 Environment Setup

### Python Environment Management

#### Conda Setup (Recommended)
```bash
# Create isolated environment
conda create -n mutationscan python=3.11
conda activate mutationscan

# Install MutationScan
cd MutationScan/subscan
pip install -e .

# Linux: Add bioinformatics tools
conda install -c bioconda abricate blast mafft
```

#### Virtual Environment Setup
```bash
# Create virtual environment
python -m venv mutationscan_env

# Activate (Linux/macOS)
source mutationscan_env/bin/activate

# Activate (Windows)
mutationscan_env\Scripts\activate

# Install
cd MutationScan/subscan
pip install -e .
```

### System Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get install -y \
    python3-dev \
    python3-pip \
    git \
    abricate \
    ncbi-blast+ \
    mafft \
    curl \
    wget
```

#### CentOS/RHEL
```bash
sudo yum install -y \
    python3-devel \
    python3-pip \
    git \
    curl \
    wget

# Install bioconda for bioinformatics tools
curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
conda install -c bioconda abricate blast mafft
```

#### macOS
```bash
# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install python git
brew install --cask docker

# Use Docker for bioinformatics tools
docker pull quay.io/biocontainers/abricate:1.0.1--ha8f3691_1
```

---

## 🐳 Docker Deployment

### Production Container
```dockerfile
FROM ubuntu:22.04

# System dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    abricate \
    ncbi-blast+ \
    mafft \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install MutationScan
WORKDIR /app
COPY . .
RUN cd subscan && pip3 install -e .

# Create data directories
RUN mkdir -p /data/input /data/output

# Entry point
ENTRYPOINT ["python3", "subscan/tools/run_pipeline.py"]
```

### Build and Run
```bash
# Build image
docker build -t mutationscan:latest .

# Run with mounted data
docker run -v $(pwd)/data:/data \
    mutationscan:latest \
    --accessions /data/input/accessions.txt \
    --gene-list /data/input/genes.txt \
    --email you@lab.org \
    --output-dir /data/output \
    --sepi-species "Escherichia coli"
```

---

## 📊 Performance Tuning

### Resource Configuration

#### CPU Optimization
```bash
# Utilize all CPU cores
python tools/run_pipeline.py --threads $(nproc) [other args]

# Conservative threading for shared systems
python tools/run_pipeline.py --threads 4 [other args]
```

#### Memory Management
```python
# Large dataset optimization (config.py)
BATCH_SIZE = 100          # Process genomes in batches
MAX_CONCURRENT = 8        # Concurrent downloads
CHUNK_SIZE = 1024 * 1024  # 1MB file chunks
```

#### Network Optimization
```python
# Federated extractor resilience settings
TIMEOUT_SECONDS = 300     # 5-minute timeout per genome
RETRY_ATTEMPTS = 3        # Database fallback tries
CONNECTION_POOL = 10      # HTTP connection pooling
```

### Scaling Guidelines

| Dataset Size | Recommended Resources | Expected Runtime |
|--------------|----------------------|------------------|
| **Small** (1-10 genomes) | 2 CPU, 4GB RAM | 5-15 minutes |
| **Medium** (10-100 genomes) | 4 CPU, 8GB RAM | 30-90 minutes |
| **Large** (100-1000 genomes) | 8+ CPU, 16GB+ RAM | 2-8 hours |
| **Enterprise** (1000+ genomes) | HPC cluster recommended | Days |

---

## 🔒 Security & Compliance

### Data Privacy
- **Local Processing**: All genome analysis occurs on your infrastructure
- **No Data Upload**: Genomes are downloaded from public databases only
- **Temporary Files**: Automatic cleanup of intermediate files
- **Access Control**: Standard Unix permissions for output directories

### Network Security
```bash
# Firewall configuration for outbound API calls
# NCBI
iptables -A OUTPUT -p tcp --dport 443 -d eutils.ncbi.nlm.nih.gov -j ACCEPT

# BV-BRC
iptables -A OUTPUT -p tcp --dport 443 -d www.bv-brc.org -j ACCEPT

# EnteroBase  
iptables -A OUTPUT -p tcp --dport 443 -d enterobase.warwick.ac.uk -j ACCEPT

# PATRIC
iptables -A OUTPUT -p tcp --dport 443 -d patricbrc.org -j ACCEPT
```

### Audit Logging
```python
# Enable comprehensive logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('mutationscan.log'),
        logging.StreamHandler()
    ]
)
```

---

## 🚀 Production Checklist

### Pre-Deployment Validation
```bash
# 1. Structural validation
python subscan/check_structure.py

# 2. Dependency verification  
python subscan/check_dependencies.py

# 3. Smoke tests
python subscan/smoke_test_dominos.py

# 4. Full pipeline test
python subscan/full_test_dominos.py

# 5. Comprehensive validation
python subscan/comprehensive_final_tests.py
```

### Monitoring & Maintenance
- **Log Rotation**: Configure `logrotate` for `mutationscan.log`
- **Disk Space**: Monitor output directories for growth
- **Database Updates**: ABRicate database updates via `abricate --setupdb`
- **Security Patches**: Regular OS and Python package updates

### Performance Benchmarks
Run the comprehensive test suite to establish baseline performance:
```bash
python subscan/comprehensive_final_tests.py
# Expected: 0.04-0.08 seconds per genome processing time
# Expected: All 4 test categories pass (Parameter/Database/Performance/Output)
```

---

## 📞 Support & Troubleshooting

### Quick Diagnostics
```bash
# Platform compatibility check
python -c "import platform; print(f'Platform: {platform.system()} {platform.release()}')"

# Python environment
python -c "import sys; print(f'Python: {sys.version}')"

# MutationScan installation
python -c "import subscan; print('MutationScan installed successfully')"

# Database connectivity
python subscan/tools/run_harvester.py --help
```

### Community Support
- **Issues**: [GitHub Issues](https://github.com/vihaankulkarni29/MutationScan/issues)
- **Documentation**: [User Guide](README.md) 
- **Email**: vihaankulkarni29@gmail.com

---

*MutationScan Deployment Guide - Professional multi-platform AMR analysis pipeline*