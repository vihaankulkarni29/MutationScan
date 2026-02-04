# MutationScan - Simplified PowerShell Run Script
# Usage: .\run_docker.ps1 -Email YOUR_EMAIL [-Genome FILE | -Organism "NAME"] [-Targets FILE]

param(
    [Parameter(Mandatory=$true)]
    [string]$Email,
    
    [Parameter(ParameterSetName='LocalGenome')]
    [string]$Genome,
    
    [Parameter(ParameterSetName='DownloadOrganism')]
    [string]$Organism,
    
    [string]$Targets,
    [string]$ApiKey,
    [int]$Limit = 1,
    [switch]$Visualize,
    [switch]$NoML
)

# Build argument list
$args = @("--email", $Email)

if ($Genome) {
    $args += @("--genome", $Genome)
} elseif ($Organism) {
    $args += @("--organism", $Organism)
    $args += @("--limit", $Limit)
} else {
    Write-Host "ERROR: Either -Genome or -Organism is required" -ForegroundColor Red
    Write-Host ""
    Write-Host "Usage examples:"
    Write-Host "  .\run_docker.ps1 -Email user@example.com -Genome data/genomes/my_genome.fasta"
    Write-Host "  .\run_docker.ps1 -Email user@example.com -Organism 'Escherichia coli' -Targets data/config/acrR_targets.txt"
    exit 1
}

if ($Targets) { $args += @("--targets", $Targets) }
if ($ApiKey) { $args += @("--api-key", $ApiKey) }
if ($Visualize) { $args += @("--visualize") }
if ($NoML) { $args += @("--no-ml") }

# Check if image exists
$imageExists = docker images -q mutationscan:latest 2>$null
if (-not $imageExists) {
    Write-Host "Building MutationScan Docker image..." -ForegroundColor Cyan
    docker build -t mutationscan:latest .
}

# Run container
Write-Host "Running MutationScan analysis..." -ForegroundColor Cyan
docker run --rm `
    -v "${PWD}/data:/app/data" `
    -v "${PWD}/models:/app/models:ro" `
    mutationscan:latest `
    $args

Write-Host ""
Write-Host "Analysis complete. Check data/results/ for output files." -ForegroundColor Green
