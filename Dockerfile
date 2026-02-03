# =============================================================================
# MutationScan - Production Dockerfile (Final)
# =============================================================================
# Optimized for "Democratization": Low bandwidth, fast deployment, offline-ready
# 
# Architecture:
# - Base: StaPH-B ABRicate image (proven bioinformatics base)
# - Layer 1: System tools + Python 3 + PyMOL
# - Layer 2: ABRicate database pre-download (offline strategy)
# - Layer 3: ML libraries (scikit-learn, joblib, matplotlib)
# - Layer 4: Python dependencies (pandas, biopython, NCBI Datasets)
# - Layer 5: Source code + ML models + permissions
#
# Result: Container is 100% ready on first run with ML predictions enabled
# Entry Point: src/main.py (Orchestrator)
# =============================================================================

FROM staphb/abricate:latest

# METADATA
LABEL maintainer="MutationScan Team"
LABEL version="2.0"
LABEL description="Clinical-Grade AMR Pipeline: NCBI -> ABRicate -> ML Prediction -> 3D Visualization"

# 1. SWITCH TO ROOT
# The base image defaults to a non-root user, but we need root to install dependencies
USER root

# 2. INSTALL SYSTEM TOOLS & PYMOL
# Added 'pymol' and 'libglew-dev' for the visualizer module
# We clean apt lists immediately to keep the image small (democratization = low bandwidth)
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    curl \
    unzip \
    zip \
    pymol \
    libglew-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. SETUP ABRICATE DATABASES (Offline Strategy)
# We force the download of CARD and NCBI databases during the build.
# This ensures the container works 100% offline once pulled (democratization principle).
RUN abricate --setupdb

# 4. INSTALL PYTHON DEPENDENCIES
# Copy requirements first to leverage Docker caching.
# Includes: pandas, biopython, scikit-learn, joblib, matplotlib (ML stack)
WORKDIR /app
COPY requirements.txt .
RUN pip3 install --no-cache-dir --break-system-packages -r requirements.txt

# 5. COPY APPLICATION CODE & ML MODELS
# Copy source code and pre-trained models into container
COPY src /app/src
COPY models /app/models

# 6. CREATE DATA DIRECTORIES
# Pre-create standard folders to prevent runtime "DirectoryNotFound" errors
RUN mkdir -p /app/data/genomes \
    /app/data/raw \
    /app/data/proteins \
    /app/data/refs \
    /app/data/results \
    /app/data/results/visualizations \
    /app/data/logs

# 7. PERMISSIONS (Security Best Practice)
# Create a dedicated non-root user for execution
RUN useradd -m bioinfo && \
    chown -R bioinfo:bioinfo /app
USER bioinfo

# 8. ENTRYPOINT (The Orchestrator)
# Default: Show help menu
# User can override: docker run mutationscan:v1 --email user@example.com --query "E. coli"
ENTRYPOINT ["python3", "src/main.py"]
CMD ["--help"]
