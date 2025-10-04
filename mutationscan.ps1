# MutationScan Pipeline Launcher
# Professional AMR Analysis Pipeline for Windows
# Compatible with Windows PowerShell 5.1+ and PowerShell Core 7+

param(
    [Parameter(ValueFromRemainingArguments=$true)]
    [string[]]$Arguments
)

Write-Host ""
Write-Host "MutationScan Pipeline" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host "Professional AMR Analysis Pipeline" -ForegroundColor Gray
Write-Host ""

# Get script directory (root of MutationScan repository)
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir

# Check if this is the correct directory structure
if (-not (Test-Path "subscan\tools\run_pipeline.py")) {
    Write-Host "[ERROR] Invalid directory structure" -ForegroundColor Red
    Write-Host "   Please ensure you're running this from the MutationScan root directory" -ForegroundColor Yellow
    Write-Host "   Expected files:" -ForegroundColor Yellow
    Write-Host "     - mutationscan.ps1 (this file)" -ForegroundColor Yellow
    Write-Host "     - subscan\tools\run_pipeline.py" -ForegroundColor Yellow
    Write-Host "     - requirements.txt" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "   Current directory: $(Get-Location)" -ForegroundColor Gray
    Read-Host "Press Enter to exit"
    exit 1
}

Write-Host "[OK] Found MutationScan directory structure" -ForegroundColor Green

# Check for virtual environment
$VenvPath = "mutationscan_env\Scripts\Activate.ps1"
if (-not (Test-Path $VenvPath)) {
    Write-Host "[ERROR] Virtual environment not found" -ForegroundColor Red
    Write-Host "   Expected: $VenvPath" -ForegroundColor Yellow
    Write-Host "   Please run the installer first:" -ForegroundColor Yellow
    Write-Host "     .\install.bat" -ForegroundColor White
    Write-Host ""
    Read-Host "Press Enter to exit"
    exit 1
}

Write-Host "[OK] Found virtual environment" -ForegroundColor Green

# Activate virtual environment
Write-Host "Activating Python environment..." -ForegroundColor Yellow
try {
    & $VenvPath
    Write-Host "[OK] Virtual environment activated" -ForegroundColor Green
} catch {
    Write-Host "[ERROR] Could not activate virtual environment: $_" -ForegroundColor Red
    Write-Host "   Try running: .\install.bat" -ForegroundColor Yellow
    Read-Host "Press Enter to exit"
    exit 1
}

# Verify Python
try {
    $pythonVersion = python --version 2>&1
    Write-Host "[OK] Python: $pythonVersion" -ForegroundColor Green
} catch {
    Write-Host "[ERROR] Python not found or not working" -ForegroundColor Red
    Write-Host "   Please ensure Python is installed and in PATH" -ForegroundColor Yellow
    Read-Host "Press Enter to exit"
    exit 1
}

# Show help if no arguments provided
if ($Arguments.Count -eq 0) {
    Write-Host ""
    Write-Host "MutationScan Pipeline Help" -ForegroundColor Cyan
    Write-Host "=============================" -ForegroundColor Cyan
    Write-Host ""
    python subscan\tools\run_pipeline.py --help
    Write-Host ""
    Write-Host "Quick Start Examples:" -ForegroundColor Cyan
    Write-Host "   # Basic analysis with sample data:" -ForegroundColor White
    Write-Host "   .\mutationscan.ps1 --accessions examples\demo_accessions.txt --gene-list examples\gene_list.txt --email your@email.com --output-dir results --sepi-species ""Escherichia coli""" -ForegroundColor Gray
    Write-Host ""
    Write-Host "   # Analysis with your own data:" -ForegroundColor White
    Write-Host "   .\mutationscan.ps1 --accessions data_input\accession_list.txt --gene-list data_input\gene_list.txt --email your@email.com --output-dir results" -ForegroundColor Gray
    Write-Host ""
    Write-Host "For more information, see README.md" -ForegroundColor Yellow
    Write-Host ""
    Read-Host "Press Enter to close"
    exit 0
}

# Run the pipeline with arguments
Write-Host ""
Write-Host "Starting MutationScan Pipeline..." -ForegroundColor Green
Write-Host "Command: python subscan\tools\run_pipeline.py $($Arguments -join ' ')" -ForegroundColor Gray
Write-Host ""

$startTime = Get-Date

try {
    python subscan\tools\run_pipeline.py @Arguments
    $exitCode = $LASTEXITCODE
    $endTime = Get-Date
    $duration = $endTime - $startTime
    
    Write-Host ""
    Write-Host "Pipeline completed in: $($duration.Minutes):$($duration.Seconds.ToString('00'))" -ForegroundColor Cyan
    
    if ($exitCode -eq 0) {
        Write-Host ""
        Write-Host "Pipeline completed successfully!" -ForegroundColor Green
        Write-Host "   Check the data_output directory for results" -ForegroundColor Gray
    } else {
        Write-Host ""
        Write-Host "Pipeline finished with errors (exit code: $exitCode)" -ForegroundColor Red
        Write-Host "   Check the error messages above for details" -ForegroundColor Yellow
        Write-Host "   For troubleshooting, see README.md or report issues" -ForegroundColor Yellow
    }
} catch {
    Write-Host ""
    Write-Host "Error running pipeline: $_" -ForegroundColor Red
    Write-Host "   Please check your input files and try again" -ForegroundColor Yellow
}

Write-Host ""
Read-Host "Press Enter to close"