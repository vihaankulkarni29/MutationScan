# MutationScan PowerShell Launcher
# Compatible with Windows PowerShell and PowerShell Core

param(
    [Parameter(ValueFromRemainingArguments=$true)]
    [string[]]$Arguments
)

Write-Host "MutationScan Pipeline Launcher" -ForegroundColor Cyan
Write-Host "================================" -ForegroundColor Cyan
Write-Host ""

# Get script directory
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

# Check if MutationScan directory exists
$MutationScanDir = Join-Path $ScriptDir "MutationScan"
if (-not (Test-Path $MutationScanDir)) {
    Write-Host "[ERROR] MutationScan directory not found at $MutationScanDir" -ForegroundColor Red
    Write-Host "        Make sure you're running this from the correct directory" -ForegroundColor Yellow
    Write-Host "        Expected structure:" -ForegroundColor Yellow
    Write-Host "          Current Directory" -ForegroundColor Yellow
    Write-Host "          +-- mutationscan.ps1 (this file)" -ForegroundColor Yellow
    Write-Host "          +-- MutationScan/" -ForegroundColor Yellow
    Read-Host "Press Enter to exit"
    exit 1
}

# Change to MutationScan directory
Set-Location $MutationScanDir
Write-Host "[OK] Found MutationScan directory: $(Get-Location)" -ForegroundColor Green

# Check if virtual environment exists
$VenvPath = Join-Path $MutationScanDir "mutationscan_env\Scripts\Activate.ps1"
if (-not (Test-Path $VenvPath)) {
    Write-Host "[ERROR] Virtual environment not found" -ForegroundColor Red
    Write-Host "        Expected: mutationscan_env\Scripts\Activate.ps1" -ForegroundColor Yellow
    Write-Host "        Please run the installer again: install.bat" -ForegroundColor Yellow
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
    Read-Host "Press Enter to exit"
    exit 1
}

# Check if pipeline script exists
$PipelineScript = Join-Path $MutationScanDir "subscan\tools\run_pipeline.py"
if (-not (Test-Path $PipelineScript)) {
    Write-Host "[ERROR] Pipeline script not found" -ForegroundColor Red
    Write-Host "        Expected: subscan\tools\run_pipeline.py" -ForegroundColor Yellow
    Write-Host "        Current directory: $(Get-Location)" -ForegroundColor Yellow
    if (Test-Path "subscan\tools\") {
        Get-ChildItem "subscan\tools\" | Format-Table Name
    }
    Read-Host "Press Enter to exit"
    exit 1
}

Write-Host "[OK] Found pipeline script" -ForegroundColor Green

# Show help if no arguments provided
if ($Arguments.Count -eq 0) {
    Write-Host "Showing help - add your arguments to run analysis" -ForegroundColor Yellow
    Write-Host ""
    python subscan\tools\run_pipeline.py --help
    Write-Host ""
    Write-Host "Example usage:" -ForegroundColor Cyan
    Write-Host "   .\mutationscan.ps1 --accessions data_input\accessions.txt --genes data_input\genes.txt --output results" -ForegroundColor White
    Write-Host ""
    Read-Host "Press Enter to close"
    exit 0
}

# Run the pipeline with arguments
Write-Host "Running MutationScan pipeline..." -ForegroundColor Yellow
Write-Host "Command: python subscan\tools\run_pipeline.py $($Arguments -join ' ')" -ForegroundColor Gray
Write-Host ""

try {
    python subscan\tools\run_pipeline.py @Arguments
    if ($LASTEXITCODE -eq 0) {
        Write-Host ""
        Write-Host "Pipeline completed successfully!" -ForegroundColor Green
    } else {
        Write-Host ""
        Write-Host "Pipeline finished with errors (exit code: $LASTEXITCODE)" -ForegroundColor Red
        Write-Host "Check the error messages above for details" -ForegroundColor Yellow
    }
} catch {
    Write-Host ""
    Write-Host "Error running pipeline: $_" -ForegroundColor Red
}

Write-Host ""
Read-Host "Press Enter to close"