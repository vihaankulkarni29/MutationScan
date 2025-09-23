# ============================================================================
# MutationScan Repository Cleanup Script (PowerShell Version)
# ============================================================================
# This script removes all files and directories that should be ignored
# according to the new comprehensive .gitignore file.
# 
# IMPORTANT: This script prepares the repository for its first public release.
# It will remove temporary files, test outputs, and cache directories while
# preserving essential source code and legitimate test directories.
#
# Usage: .\cleanup_repository.ps1
# ============================================================================

Write-Host "🧹 MutationScan Repository Cleanup - Phase 2" -ForegroundColor Cyan
Write-Host "==============================================" -ForegroundColor Cyan
Write-Host ""

# Function to safely remove files/directories if they exist
function Safe-Remove {
    param([string]$Target)
    
    if (Test-Path $Target) {
        Write-Host "🗑️  Removing: $Target" -ForegroundColor Yellow
        Remove-Item -Path $Target -Recurse -Force -ErrorAction SilentlyContinue
    } else {
        Write-Host "✅ Already clean: $Target" -ForegroundColor Green
    }
}

Write-Host "📋 Step 1: Un-tracking all files and re-adding with new .gitignore rules" -ForegroundColor Blue
Write-Host "----------------------------------------------------------------------" -ForegroundColor Blue

# Un-track all files while keeping them locally
Write-Host "🔄 Un-tracking all currently tracked files..." -ForegroundColor Yellow
try {
    git rm -r --cached . 2>$null
} catch {
    # Ignore errors, continue
}

# Re-add all files (this will now respect the new .gitignore)
Write-Host "➕ Re-adding files (respecting new .gitignore rules)..." -ForegroundColor Yellow
git add .

Write-Host "✅ Git tracking updated successfully!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 2: Removing Python cache and compiled files" -ForegroundColor Blue
Write-Host "---------------------------------------------------" -ForegroundColor Blue

# Remove __pycache__ directories and .pyc files
Get-ChildItem -Path . -Recurse -Directory -Name "__pycache__" -ErrorAction SilentlyContinue | ForEach-Object {
    Safe-Remove $_
}

Get-ChildItem -Path . -Recurse -File -Include "*.pyc", "*.pyo", "*.pyd" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

Write-Host "✅ Python cache files cleaned!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 3: Removing output directories and test results" -ForegroundColor Blue
Write-Host "------------------------------------------------------" -ForegroundColor Blue

# Root-level output directories
Safe-Remove "data_output"
Safe-Remove "demo_output"
Safe-Remove "test_output"
Safe-Remove "results_output"
Safe-Remove "analysis_output"

# Test result directories at root level
Safe-Remove "harvester_test_results"
Safe-Remove "annotator_test_results"

# Subscan-level output directories
Safe-Remove "subscan\harvester_output"
Safe-Remove "subscan\test_parallel_output"
Safe-Remove "subscan\harvester_test"
Safe-Remove "subscan\harvester_test_results"

# Pipeline-specific output directories
Safe-Remove "subscan\annotator_output"
Safe-Remove "subscan\extractor_output"
Safe-Remove "subscan\aligner_output"
Safe-Remove "subscan\analyzer_output"
Safe-Remove "subscan\reporter_output"

Write-Host "✅ Output directories cleaned!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 4: Removing temporary and debug files" -ForegroundColor Blue
Write-Host "--------------------------------------------" -ForegroundColor Blue

# Root-level temporary test files
Safe-Remove "test_html_reporter.py"
Safe-Remove "test_accessions.txt"

# Subscan-level debug and temporary files
Safe-Remove "subscan\debug_paths.py"
Safe-Remove "subscan\test_manifest_fix.py"

# Tools-level temporary files
Safe-Remove "subscan\tools\test_genome_manifest.json"
Safe-Remove "subscan\tools\test_alignment.water"

# Find and remove other temporary files
Get-ChildItem -Path . -Recurse -File -Include "*.tmp", "*.temp" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

Get-ChildItem -Path . -Recurse -File -Include "debug_*.py", "temp_*.py", "mock_*.py" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

Write-Host "✅ Temporary files cleaned!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 5: Removing log files and validation guides" -ForegroundColor Blue
Write-Host "--------------------------------------------------" -ForegroundColor Blue

# Remove log files
Get-ChildItem -Path . -Recurse -File -Include "*.log", "log_*.txt", "debug_*.txt" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

# Remove validation guides
Safe-Remove "VALIDATION_RUN_GUIDE.md"
Get-ChildItem -Path . -Recurse -File -Include "validation_*.md" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

# Remove disabled test files
Get-ChildItem -Path . -Recurse -File -Include "disabled_test_*.py", "*_disabled_test.py" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

Write-Host "✅ Log files and validation guides cleaned!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 6: Removing IDE and editor configuration files" -ForegroundColor Blue
Write-Host "-----------------------------------------------------" -ForegroundColor Blue

# IDE directories
Safe-Remove ".vscode"
Safe-Remove ".idea"

# Editor temporary files
Get-ChildItem -Path . -Recurse -File -Include "*.swp", "*.swo", "*~" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

Write-Host "✅ IDE configuration cleaned!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 7: Removing bioinformatics temporary files" -ForegroundColor Blue
Write-Host "-------------------------------------------------" -ForegroundColor Blue

# Temporary sequence files
Get-ChildItem -Path . -Recurse -File -Include "*.fasta.tmp", "*.fastq.tmp", "*.fa.tmp", "*.fq.tmp" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

# Temporary alignment files
Get-ChildItem -Path . -Recurse -File -Include "*.aln.tmp", "*.alignment.tmp", "*.water", "*.needle" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

# Temporary database files
Get-ChildItem -Path . -Recurse -File -Include "*.db.tmp", "*.sqlite.tmp" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

Write-Host "✅ Bioinformatics temporary files cleaned!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 8: Removing system files" -ForegroundColor Blue
Write-Host "-------------------------------" -ForegroundColor Blue

# Windows system files
Get-ChildItem -Path . -Recurse -File -Include "Thumbs.db", "Desktop.ini", ".DS_Store" -ErrorAction SilentlyContinue | ForEach-Object {
    Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue
}

Write-Host "✅ System files cleaned!" -ForegroundColor Green
Write-Host ""

Write-Host "📋 Step 9: Verifying essential directories are preserved" -ForegroundColor Blue
Write-Host "------------------------------------------------------" -ForegroundColor Blue

# Check that essential directories still exist
$essentialDirs = @(
    "subscan\src",
    "subscan\tools",
    "subscan\tests",
    "subscan\analyzer\tests",
    "subscan\ncbi_genome_extractor\tests"
)

foreach ($dir in $essentialDirs) {
    if (Test-Path $dir) {
        Write-Host "✅ Preserved: $dir" -ForegroundColor Green
    } else {
        Write-Host "⚠️  Warning: $dir not found (this may be expected)" -ForegroundColor Yellow
    }
}

Write-Host ""
Write-Host "📋 Step 10: Final git status check" -ForegroundColor Blue
Write-Host "---------------------------------" -ForegroundColor Blue

Write-Host "📊 Current repository status:" -ForegroundColor Yellow
git status --porcelain | Select-Object -First 20

Write-Host ""
Write-Host "🎉 Repository Cleanup Complete!" -ForegroundColor Green
Write-Host "===============================" -ForegroundColor Green
Write-Host ""
Write-Host "📁 The following actions were completed:" -ForegroundColor Cyan
Write-Host "   ✅ Updated Git tracking to respect new .gitignore" -ForegroundColor Green
Write-Host "   ✅ Removed Python cache files and compiled bytecode" -ForegroundColor Green
Write-Host "   ✅ Removed all output directories and test results" -ForegroundColor Green
Write-Host "   ✅ Removed temporary and debug files" -ForegroundColor Green
Write-Host "   ✅ Removed log files and validation guides" -ForegroundColor Green
Write-Host "   ✅ Removed IDE configuration files" -ForegroundColor Green
Write-Host "   ✅ Removed bioinformatics temporary files" -ForegroundColor Green
Write-Host "   ✅ Removed system files" -ForegroundColor Green
Write-Host "   ✅ Preserved essential source code and test directories" -ForegroundColor Green
Write-Host ""
Write-Host "🚀 Your MutationScan repository is now ready for public release!" -ForegroundColor Cyan
Write-Host ""
Write-Host "📝 Next steps:" -ForegroundColor Yellow
Write-Host "   1. Review the changes: git status" -ForegroundColor White
Write-Host "   2. Commit the cleanup: git add ." -ForegroundColor White
Write-Host "      Then: git commit -m 'Repository cleanup for public release'" -ForegroundColor White
Write-Host "   3. Push to GitHub: git push origin main" -ForegroundColor White
Write-Host ""
Write-Host "🔗 Repository ready for: https://github.com/vihaankulkarni29/MutationScan" -ForegroundColor Cyan