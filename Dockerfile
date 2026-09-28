# MutationScan.
#
#   docker build -t mutationscan .
#   docker run --rm -v "$PWD/data:/app/data" mutationscan \
#       run --genomes data/genomes --references data/references --out data/output
#
# Miniforge because BLAST+ comes from bioconda; everything else is pure Python.
FROM condaforge/miniforge3:latest

WORKDIR /app

# Environment first, so a code change does not invalidate the conda layer.
COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean -afy
ENV PATH="/opt/conda/envs/mutationscan/bin:$PATH"

COPY . .
RUN pip install --no-deps -e .

# `docker run mutationscan <args>` is `mutationscan <args>`.
ENTRYPOINT ["mutationscan"]
CMD ["--help"]
