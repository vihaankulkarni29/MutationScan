# BASE IMAGE: State Public Health Bioinformatics (StaPH-B)
FROM staphb/abricate:latest

# METADATA
LABEL maintainer="MutationScan Team"
LABEL version="1.0"
LABEL description="Democratized AMR Pipeline"

# 1. SWITCH TO ROOT & INSTALL SYSTEM TOOLS
USER root
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

# 2. SETUP OFFLINE DATABASES
RUN abricate --setupdb

# 3. INSTALL PYTHON DEPENDENCIES
WORKDIR /app
COPY requirements.txt .
# We explicitly install the package in editable mode so 'src' layout works
RUN pip3 install --no-cache-dir --break-system-packages -r requirements.txt

# 4. COPY APPLICATION ARCHITECTURE
# We must copy the specific folders to match the new structure
COPY src /app/src
COPY config /app/config
COPY models /app/models
COPY pyproject.toml /app/

# 5. CONFIGURE ENVIRONMENT
# This line is CRITICAL: It tells Python to look inside /app/src for modules
ENV PYTHONPATH="${PYTHONPATH}:/app/src"

# 6. PERMISSIONS
RUN useradd -m bioinfo && \
    chown -R bioinfo:bioinfo /app
USER bioinfo

# 7. ENTRYPOINT (Updated for Module Execution)
# We now run it as a module (-m) instead of a script file
ENTRYPOINT ["python3", "-m", "mutation_scan.main"]
CMD ["--help"]
