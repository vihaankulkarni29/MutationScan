@echo off
setlocal enabledelayedexpansion

REM MutationScan Pipeline Launcher
REM Professional AMR Analysis Pipeline for Windows
REM Compatible with Windows 10/11

echo.
echo 🧬 MutationScan Pipeline
echo ================================
echo Professional AMR Analysis Pipeline
echo.

REM Get script directory (root of MutationScan repository)
cd /d "%~dp0"

REM Check if this is the correct directory structure
if not exist "subscan\tools\run_pipeline.py" (
    echo [❌ ERROR] Invalid directory structure
    echo    Please ensure you're running this from the MutationScan root directory
    echo    Expected files:
    echo      ✓ mutationscan.bat (this file)
    echo      ✓ subscan\tools\run_pipeline.py
    echo      ✓ requirements.txt
    echo.
    echo    Current directory: %CD%
    pause
    exit /b 1
)

echo [✅ OK] Found MutationScan directory structure

REM Check for virtual environment
if not exist "mutationscan_env\Scripts\activate.bat" (
    echo [❌ ERROR] Virtual environment not found
    echo    Expected: mutationscan_env\Scripts\activate.bat
    echo    Please run the installer first:
    echo      .\install.bat
    echo.
    pause
    exit /b 1
)

echo [✅ OK] Found virtual environment

REM Activate virtual environment
echo 🔧 Activating Python environment...
call mutationscan_env\Scripts\activate.bat
if errorlevel 1 (
    echo [❌ ERROR] Could not activate virtual environment
    echo    Try running: .\install.bat
    pause
    exit /b 1
)

echo [✅ OK] Virtual environment activated

REM Verify Python and key packages
python --version >nul 2>&1
if errorlevel 1 (
    echo [❌ ERROR] Python not found or not working
    echo    Please ensure Python is installed and in PATH
    pause
    exit /b 1
)

python -c "import tqdm; print('tqdm OK')" >nul 2>&1
if errorlevel 1 (
    echo [⚠️  WARNING] Missing dependencies detected
    echo    Installing missing packages...
    pip install -r requirements.txt
    if errorlevel 1 (
        echo [❌ ERROR] Failed to install dependencies
        pause
        exit /b 1
    )
    echo [✅ OK] Dependencies installed
) else (
    echo [✅ OK] Python dependencies verified
)

REM Show help if no arguments provided
if "%~1"=="" (
    echo.
    echo 📖 MutationScan Pipeline Help
    echo =============================
    echo.
    python subscan\tools\run_pipeline.py --help
    echo.
    echo 💡 Quick Start Examples:
    echo    # Basic analysis with sample data:
    echo    .\mutationscan.bat --accessions examples\demo_accessions.txt --gene-list examples\gene_list.txt --email your@email.com --output-dir results --sepi-species "Escherichia coli"
    echo.
    echo    # Analysis with your own data:
    echo    .\mutationscan.bat --accessions data_input\accession_list.txt --gene-list data_input\gene_list.txt --email your@email.com --output-dir results
    echo.
    echo 📚 For more information, see README.md
    echo.
    pause
    exit /b 0
)

REM Run the pipeline with arguments
echo.
echo 🚀 Starting MutationScan Pipeline...
echo Command: python subscan\tools\run_pipeline.py %*
echo.

set start_time=%time%

python subscan\tools\run_pipeline.py %*
set exit_code=%errorlevel%

set end_time=%time%

echo.
echo ⏱️  Pipeline completed

if %exit_code% equ 0 (
    echo.
    echo 🎉 Pipeline completed successfully!
    echo    Check the data_output directory for results
) else (
    echo.
    echo ❌ Pipeline finished with errors (exit code: %exit_code%)
    echo    Check the error messages above for details
    echo    For troubleshooting, see README.md or report issues
)

echo.
pause