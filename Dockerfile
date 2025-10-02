# Use a stable Conda base image
FROM continuumio/miniconda3:latest

# Set maintainer information
LABEL maintainer="vihaankulkarni29"
LABEL description="MutationScan: Complete AMR mutation analysis pipeline"
LABEL version="1.0.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_ENV_NAME=mutationscan_env
ENV PATH="/opt/conda/envs/${CONDA_ENV_NAME}/bin:$PATH"

# Set working directory
WORKDIR /app

# Update system and install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Configure Conda channels with strict priority
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels defaults && \
    conda config --set channel_priority strict

# Create a new Conda environment and install system bioinformatics tools
RUN conda create -n ${CONDA_ENV_NAME} python=3.9 -y && \
    conda install -n ${CONDA_ENV_NAME} -c conda-forge -c bioconda \
    git \
    abricate \
    -y && \
    conda clean -all -y

# Copy project files (excluding items in .dockerignore)
COPY . /app/

# Change to the subscan directory where setup.cfg is located
WORKDIR /app/subscan

# Install the MutationScan pipeline and all dependencies
# This will:
# 1. Install the subscan package itself
# 2. Install all Python dependencies from install_requires
# 3. Use git to clone and install all seven Domino tools from GitHub
RUN /opt/conda/envs/${CONDA_ENV_NAME}/bin/pip install .

# Verify installation by checking if key tools are available
RUN /opt/conda/envs/${CONDA_ENV_NAME}/bin/python -c "import subscan; print('SubScan installed successfully')" && \
    /opt/conda/envs/${CONDA_ENV_NAME}/bin/python -c "import ncbi_genome_extractor; print('NCBI Genome Extractor installed')" && \
    /opt/conda/envs/${CONDA_ENV_NAME}/bin/abricate --version

# Set the working directory back to the main app directory
WORKDIR /app

# Create directories for input and output data
RUN mkdir -p /app/data_input /app/data_output /app/results

# Set proper permissions
RUN chmod +x /app/subscan/tools/run_pipeline.py

# Set the entrypoint to the master orchestrator script
# This makes the container itself an executable pipeline
ENTRYPOINT ["/opt/conda/envs/mutationscan_env/bin/python", "/app/subscan/tools/run_pipeline.py"]

# Default command (can be overridden)
CMD ["--help"]

# Add health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD /opt/conda/envs/${CONDA_ENV_NAME}/bin/python -c "import subscan" || exit 1