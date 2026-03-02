# BASE IMAGE: State Public Health Bioinformatics (StaPH-B)
FROM staphb/abricate:latest

# METADATA
LABEL maintainer="MutationScan Team"
LABEL version="2.1"
LABEL description="Democratized AMR Pipeline (Genomics + Biophysics)"

# 1. SWITCH TO ROOT & INSTALL SYSTEM DEPENDENCIES
USER root
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    bzip2 \
    unzip \
    zip \
    libgl1 \
    libegl1 \
    libxrandr2 \
    libxss1 \
    libxcursor1 \
    libxcomposite1 \
    libasound2t64 \
    libxi6 \
    libxtst6 \
    && rm -rf /var/lib/apt/lists/*

# 2. INSTALL MINICONDA
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# 3. CONDA PRE-LOAD (Physics + Data Science Stack)
# We load ncbi-datasets-pylib and ML tools here to bypass pip compiler crashes
RUN /opt/conda/bin/conda create -n biophysics --override-channels -c conda-forge -c bioconda python=3.10 \
    openbabel vina pymol-open-source openmm pdbfixer \
    pandas numpy biopython scikit-learn joblib matplotlib seaborn ncbi-datasets-pylib -y

# 4. ACTIVATE ENVIRONMENT
ENV PATH=/opt/conda/envs/biophysics/bin:$CONDA_DIR/bin:$PATH

# 5. INSTALL NCBI DATASETS CLI
RUN curl -L -o /usr/local/bin/datasets \
    "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets" && \
    chmod +x /usr/local/bin/datasets

# 6. SETUP ABRICATE OFFLINE DATABASES
RUN abricate --setupdb

# 7. INSTALL PYTHON PIPELINE DEPENDENCIES
WORKDIR /app
COPY requirements.txt .
# Pip will now gracefully skip the heavy packages and only install AutoScan
RUN pip install --no-cache-dir setuptools==69.5.1 wheel pbr && \
    pip install --no-cache-dir -r requirements.txt

# 8. COPY APPLICATION ARCHITECTURE
COPY src /app/src
COPY config /app/config
COPY models /app/models
COPY pyproject.toml /app/

# 9. CONFIGURE ENVIRONMENT
ENV PYTHONPATH="${PYTHONPATH}:/app/src"

# 10. PERMISSIONS
RUN useradd -m bioinfo && \
    chown -R bioinfo:bioinfo /app
USER bioinfo

# 11. ENTRYPOINT
ENTRYPOINT ["python", "-m", "mutation_scan.main"]
CMD ["--help"]
