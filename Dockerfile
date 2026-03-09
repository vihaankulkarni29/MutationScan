# Use Miniforge for fast, robust Conda environment solving
FROM condaforge/miniforge3:latest

# Set working directory
WORKDIR /app

# Install system-level dependencies required by SMINA and OpenMM
RUN apt-get update && apt-get install -y \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the environment definition
COPY environment.yml .

# Create the conda environment and clean up caches to keep the image small
RUN mamba env create -f environment.yml && mamba clean -afy

# Activate the conda environment globally for all subsequent commands
ENV PATH="/opt/conda/envs/mutationscan/bin:$PATH"

# Copy the entire codebase into the container
COPY . .

# Set Snakemake as the default entrypoint
# Set Snakemake as the default command (entrypoint can be overridden if needed)
CMD ["snakemake", "--help"]
