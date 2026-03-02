# BASE IMAGE: StaPH-B ABRicate (includes NCBI tools)
FROM staphb/abricate:latest

LABEL maintainer="MutationScan Team"
LABEL version="2.1"
LABEL description="Grand Unified MutationScan: Genomic + Biophysics Pipeline with Native OpenMM/Vina"

# 1. SWITCH TO ROOT & INSTALL SYSTEM DEPENDENCIES
USER root

RUN apt-get update && apt-get install -y \
    curl wget bzip2 unzip zip \
    libgl1 libegl1 libxrandr2 libxss1 libxcursor1 libxcomposite1 \
    libasound2t64 libxi6 libxtst6 \
    && rm -rf /var/lib/apt/lists/*

# 2. INSTALL MINICONDA
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# 3. INSTALL BIOPHYSICS STACK (Isolated Environment)
RUN /opt/conda/bin/conda create -n biophysics --override-channels -c conda-forge python=3.10 openbabel vina pymol-open-source openmm pdbfixer pandas -y
ENV PATH=/opt/conda/envs/biophysics/bin:$CONDA_DIR/bin:$PATH

# 4. INSTALL NCBI DATASETS CLI
RUN curl -L -o /usr/local/bin/datasets \
    "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets" && \
    chmod +x /usr/local/bin/datasets

# 5. SETUP ABRicate DATABASES
RUN abricate --setupdb

# 6. APPLICATION SETUP
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir "setuptools<70.0.0" wheel && \
    pip install --no-cache-dir --no-build-isolation -r requirements.txt

COPY src /app/src
COPY config /app/config
COPY models /app/models
COPY pyproject.toml /app/

# 7. CONFIGURE PYTHON PATH
ENV PYTHONPATH="${PYTHONPATH}:/app/src"

# 8. PERMISSIONS & USER
RUN useradd -m bioinfo && \
    chown -R bioinfo:bioinfo /app

USER bioinfo

# 9. ENTRYPOINT
ENTRYPOINT ["python", "-m", "mutation_scan.main"]
CMD ["--help"]

# 7. ENTRYPOINT (Updated for Module Execution)
# We now run it as a module (-m) instead of a script file
ENTRYPOINT ["python3", "-m", "mutation_scan.main"]
CMD ["--help"]
