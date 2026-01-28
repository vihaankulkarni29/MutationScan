# Multi-stage Docker build for MutationScan

# Stage 1: Base image with dependencies
FROM python:3.10-slim as base

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    wget \
    ncbi-blast+ \
    abricate \
    && rm -rf /var/lib/apt/lists/*

# Create application directory
WORKDIR /app

# Stage 2: Development image
FROM base as development

# Install development dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    vim \
    less \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements
COPY requirements.txt .

# Install Python dependencies with dev extras
RUN pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir \
    pytest>=6.2.5 \
    pytest-cov>=2.12.1 \
    black>=21.7b0 \
    flake8>=3.9.2 \
    mypy>=0.910 \
    isort>=5.9.3

# Copy source code
COPY . .

# Install package in development mode
RUN pip install -e .

# Stage 3: Production image
FROM base as production

# Copy requirements
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY . .

# Install package
RUN pip install .

# Create non-root user
RUN useradd -m -u 1000 mutationscan && \
    chown -R mutationscan:mutationscan /app

USER mutationscan

# Set up volumes for data
VOLUME ["/app/data", "/app/logs"]

# Default command
ENTRYPOINT ["python", "-m", "mutation_scan"]

# Stage 4: Testing image
FROM development as testing

# Set PYTHONPATH for testing
ENV PYTHONPATH=/app/src

# Run tests as default command
CMD ["pytest", "-v", "--cov=src", "tests/"]

# Final: Production ready image
FROM production as final

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import mutation_scan; print('healthy')" || exit 1

LABEL maintainer="Your Organization <contact@example.com>" \
      version="0.1.0" \
      description="MutationScan: Bioinformatics pipeline for AMR gene detection"
