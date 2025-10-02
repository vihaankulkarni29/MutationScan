# Use a specific, stable version of miniconda3 as the base image for reproducibility
FROM continuumio/miniconda3:4.12.0

# Add metadata labels for professionalism and traceability
LABEL maintainer="Vihaan Kulkarni <vihaankulkarni29@gmail.com>"
LABEL version="1.0.0"
LABEL description="MutationScan: A modular pipeline for targeted, reference-driven mutation analysis of AMR genes."

# Set the working directory inside the container
WORKDIR /app

# Copy the .dockerignore file first to leverage caching
COPY .dockerignore .dockerignore

# Copy all project files into the container, respecting the .dockerignore
COPY . .

# Configure Conda channels with strict priority, create the environment, and install Abricate.
# This is a critical step for ensuring dependencies are resolved correctly.
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda create -n mutationscan_env python=3.9 abricate -y && \
    conda clean -afy

# Install all Python package dependencies from the setup.cfg into the new Conda environment.
# This includes all of our standalone "Domino" tools from GitHub.
RUN /opt/conda/envs/mutationscan_env/bin/pip install --no-cache-dir .

# Set the default command to run when the container starts.
# This makes the container act like a dedicated executable for our master orchestrator.
ENTRYPOINT ["conda", "run", "-n", "mutationscan_env", "python", "tools/run_pipeline.py"]