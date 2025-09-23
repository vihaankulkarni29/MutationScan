#!/bin/bash
# ============================================================================
# MutationScan Repository Cleanup Script
# ============================================================================
# This script removes all files and directories that should be ignored
# according to the new comprehensive .gitignore file.
# 
# IMPORTANT: This script prepares the repository for its first public release.
# It will remove temporary files, test outputs, and cache directories while
# preserving essential source code and legitimate test directories.
#
# Usage: ./cleanup_repository.sh
# ============================================================================

echo "🧹 MutationScan Repository Cleanup - Phase 2"
echo "=============================================="
echo ""

# Set error handling
set -e

# Function to safely remove files/directories if they exist
safe_remove() {
    local target="$1"
    if [ -e "$target" ]; then
        echo "🗑️  Removing: $target"
        rm -rf "$target"
    else
        echo "✅ Already clean: $target"
    fi
}

# Function to remove files using git (for tracked files)
git_remove() {
    local target="$1"
    if [ -e "$target" ]; then
        echo "📝 Git removing: $target"
        git rm -rf --cached "$target" 2>/dev/null || true
        rm -rf "$target"
    fi
}

echo "📋 Step 1: Un-tracking all files and re-adding with new .gitignore rules"
echo "----------------------------------------------------------------------"

# Un-track all files while keeping them locally
echo "🔄 Un-tracking all currently tracked files..."
git rm -r --cached . 2>/dev/null || true

# Re-add all files (this will now respect the new .gitignore)
echo "➕ Re-adding files (respecting new .gitignore rules)..."
git add .

echo "✅ Git tracking updated successfully!"
echo ""

echo "📋 Step 2: Removing Python cache and compiled files"
echo "---------------------------------------------------"

# Remove __pycache__ directories
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -name "*.pyc" -delete 2>/dev/null || true
find . -name "*.pyo" -delete 2>/dev/null || true
find . -name "*.pyd" -delete 2>/dev/null || true

echo "✅ Python cache files cleaned!"
echo ""

echo "📋 Step 3: Removing output directories and test results"
echo "------------------------------------------------------"

# Root-level output directories
safe_remove "data_output"
safe_remove "demo_output" 
safe_remove "test_output"
safe_remove "results_output"
safe_remove "analysis_output"

# Test result directories at root level
safe_remove "harvester_test_results"
safe_remove "annotator_test_results"

# Subscan-level output directories
safe_remove "subscan/harvester_output"
safe_remove "subscan/test_parallel_output"
safe_remove "subscan/harvester_test"
safe_remove "subscan/harvester_test_results"

# Pipeline-specific output directories
safe_remove "subscan/harvester_output"
safe_remove "subscan/annotator_output"
safe_remove "subscan/extractor_output"
safe_remove "subscan/aligner_output"
safe_remove "subscan/analyzer_output"
safe_remove "subscan/reporter_output"

echo "✅ Output directories cleaned!"
echo ""

echo "📋 Step 4: Removing temporary and debug files"
echo "--------------------------------------------"

# Root-level temporary test files
safe_remove "test_html_reporter.py"
safe_remove "test_accessions.txt"

# Subscan-level debug and temporary files
safe_remove "subscan/debug_paths.py"
safe_remove "subscan/test_manifest_fix.py"

# Tools-level temporary files
safe_remove "subscan/tools/test_genome_manifest.json"
safe_remove "subscan/tools/test_alignment.water"

# Find and remove other temporary files
find . -name "*.tmp" -delete 2>/dev/null || true
find . -name "*.temp" -delete 2>/dev/null || true
find . -name "debug_*.py" -delete 2>/dev/null || true
find . -name "temp_*.py" -delete 2>/dev/null || true
find . -name "mock_*.py" -delete 2>/dev/null || true

echo "✅ Temporary files cleaned!"
echo ""

echo "📋 Step 5: Removing log files and validation guides"
echo "--------------------------------------------------"

# Remove log files
find . -name "*.log" -delete 2>/dev/null || true
find . -name "log_*.txt" -delete 2>/dev/null || true
find . -name "debug_*.txt" -delete 2>/dev/null || true

# Remove validation guides
safe_remove "VALIDATION_RUN_GUIDE.md"
find . -name "validation_*.md" -delete 2>/dev/null || true

# Remove disabled test files
find . -name "disabled_test_*.py" -delete 2>/dev/null || true
find . -name "*_disabled_test.py" -delete 2>/dev/null || true

echo "✅ Log files and validation guides cleaned!"
echo ""

echo "📋 Step 6: Removing IDE and editor configuration files"
echo "-----------------------------------------------------"

# IDE directories
safe_remove ".vscode"
safe_remove ".idea"

# Editor temporary files
find . -name "*.swp" -delete 2>/dev/null || true
find . -name "*.swo" -delete 2>/dev/null || true
find . -name "*~" -delete 2>/dev/null || true

echo "✅ IDE configuration cleaned!"
echo ""

echo "📋 Step 7: Removing bioinformatics temporary files"
echo "-------------------------------------------------"

# Temporary sequence files
find . -name "*.fasta.tmp" -delete 2>/dev/null || true
find . -name "*.fastq.tmp" -delete 2>/dev/null || true
find . -name "*.fa.tmp" -delete 2>/dev/null || true
find . -name "*.fq.tmp" -delete 2>/dev/null || true

# Temporary alignment files
find . -name "*.aln.tmp" -delete 2>/dev/null || true
find . -name "*.alignment.tmp" -delete 2>/dev/null || true
find . -name "*.water" -delete 2>/dev/null || true
find . -name "*.needle" -delete 2>/dev/null || true

# Temporary database files
find . -name "*.db.tmp" -delete 2>/dev/null || true
find . -name "*.sqlite.tmp" -delete 2>/dev/null || true

echo "✅ Bioinformatics temporary files cleaned!"
echo ""

echo "📋 Step 8: Removing system files"
echo "-------------------------------"

# macOS system files
find . -name ".DS_Store" -delete 2>/dev/null || true
find . -name ".AppleDouble" -delete 2>/dev/null || true
find . -name "._*" -delete 2>/dev/null || true

# Windows system files
find . -name "Thumbs.db" -delete 2>/dev/null || true
find . -name "Desktop.ini" -delete 2>/dev/null || true

echo "✅ System files cleaned!"
echo ""

echo "📋 Step 9: Verifying essential directories are preserved"
echo "------------------------------------------------------"

# Check that essential directories still exist
essential_dirs=(
    "subscan/src"
    "subscan/tools"
    "subscan/tests"
    "subscan/analyzer/tests"
    "subscan/ncbi_genome_extractor/tests"
)

for dir in "${essential_dirs[@]}"; do
    if [ -d "$dir" ]; then
        echo "✅ Preserved: $dir"
    else
        echo "⚠️  Warning: $dir not found (this may be expected)"
    fi
done

echo ""
echo "📋 Step 10: Final git status check"
echo "---------------------------------"

echo "📊 Current repository status:"
git status --porcelain | head -20

echo ""
echo "🎉 Repository Cleanup Complete!"
echo "==============================="
echo ""
echo "📁 The following actions were completed:"
echo "   ✅ Updated Git tracking to respect new .gitignore"
echo "   ✅ Removed Python cache files and compiled bytecode"
echo "   ✅ Removed all output directories and test results"
echo "   ✅ Removed temporary and debug files"
echo "   ✅ Removed log files and validation guides"
echo "   ✅ Removed IDE configuration files"
echo "   ✅ Removed bioinformatics temporary files"
echo "   ✅ Removed system files"
echo "   ✅ Preserved essential source code and test directories"
echo ""
echo "🚀 Your MutationScan repository is now ready for public release!"
echo ""
echo "📝 Next steps:"
echo "   1. Review the changes: git status"
echo "   2. Commit the cleanup: git add . && git commit -m 'Repository cleanup for public release'"
echo "   3. Push to GitHub: git push origin main"
echo ""
echo "🔗 Repository ready for: https://github.com/vihaankulkarni29/MutationScan"