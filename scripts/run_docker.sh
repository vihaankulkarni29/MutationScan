#!/bin/bash
# MutationScan - Simplified Docker Run Script
# Usage: ./run_docker.sh --email YOUR_EMAIL [--genome FILE | --organism "NAME"] [--targets FILE]

set -e

# Parse arguments
ARGS="$@"

# Validate minimum requirements
if [[ ! "$ARGS" =~ --email ]]; then
    echo "ERROR: --email is required"
    echo ""
    echo "Usage examples:"
    echo "  ./run_docker.sh --email user@example.com --genome data/genomes/my_genome.fasta --targets data/config/my_genes.txt"
    echo "  ./run_docker.sh --email user@example.com --organism \"Escherichia coli\" --targets data/config/acrR_targets.txt --api-key YOUR_KEY"
    exit 1
fi

if [[ ! "$ARGS" =~ --genome ]] && [[ ! "$ARGS" =~ --organism ]]; then
    echo "ERROR: Either --genome or --organism is required"
    echo ""
    echo "Usage examples:"
    echo "  ./run_docker.sh --email user@example.com --genome data/genomes/my_genome.fasta"
    echo "  ./run_docker.sh --email user@example.com --organism \"Escherichia coli\""
    exit 1
fi

# Build image if it doesn't exist
if [[ "$(docker images -q mutationscan:latest 2> /dev/null)" == "" ]]; then
    echo "Building MutationScan Docker image..."
    docker build -t mutationscan:latest .
fi

# Run container with mounted data directory
echo "Running MutationScan analysis..."
docker run --rm \
    -v "$(pwd)/data:/app/data" \
    -v "$(pwd)/models:/app/models:ro" \
    mutationscan:latest \
    $ARGS

echo ""
echo "Analysis complete. Check data/results/ for output files."
