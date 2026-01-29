# =============================================================================
# MutationScan - Production Dockerfile
# =============================================================================
# Optimized for "Democratization": Low bandwidth, fast deployment, offline-ready
# 
# Architecture:
# - Base: StaPH-B ABRicate image (proven bioinformatics base)
# - Layer 1: Python 3 installation
# - Layer 2: ABRicate database pre-download (20-30 min at build time)
# - Layer 3: Python dependencies (pandas, biopython, NCBI Datasets)
# - Layer 4: Directory structure & permissions
#
# Result: Container is 100% ready on first run (no runtime database downloads)
# =============================================================================

FROM staphb/abricate:latest

# METADATA
LABEL maintainer="MutationScan Team"
LABEL version="1.0"
LABEL description="Democratized AMR Pipeline: NCBI Ingestion -> ABRicate/BLAST -> Variant Calling"

# 1. SWITCH TO ROOT
# The base image defaults to a non-root user, but we need root to install Python.
USER root

# 2. INSTALL PYTHON & SYSTEM TOOLS
# We clean apt lists immediately to keep the image small (Democratization = Low Bandwidth)
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    curl \
    unzip \
    zip \
    && rm -rf /var/lib/apt/lists/*

# 3. SETUP ABRICATE DATABASES (The "Offline" Fix)
# We force the download of CARD and NCBI databases during the build.
# This ensures the container works 100% offline once pulled.
RUN abricate --setupdb

# 4. INSTALL PYTHON DEPENDENCIES
# We copy requirements first to leverage Docker caching.
WORKDIR /app
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# 5. SETUP DIRECTORY STRUCTURE
# Create standard folders so the script never fails on "DirectoryNotFound"
RUN mkdir -p /app/data/genomes \
    /app/data/raw \
    /app/data/results \
    /app/logs \
    /app/refs

# 6. COPY SOURCE CODE
# Copy the entire src folder
COPY src /app/src

# 7. PERMISSIONS (Security Best Practice)
# Create a dedicated user so we don't run as root
RUN useradd -m bioinfo && \
    chown -R bioinfo:bioinfo /app
USER bioinfo

# 8. ENTRYPOINT
# Default command prints the python version to verify install, 
# but user can override this to run specific scripts.
CMD ["python3", "--version"]
