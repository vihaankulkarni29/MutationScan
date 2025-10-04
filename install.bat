@echo off
setlocal enabledelayedexpansion

REM MutationScan One-Click Installer for Windows
REM Professional AMR Analysis Pipeline Installation
REM Compatible with Windows 10/11

echo.
echo 🧬 MutationScan Pipeline Installer
echo =====================================
echo Professional AMR Analysis Pipeline
echo Installing dependencies and setting up environment...
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo [❌ ERROR] Python not found
    echo    Please install Python 3.9+ from https://python.org
    echo    Make sure to check "Add Python to PATH" during installation
    echo.
    echo    After installing Python, run this installer again
    pause
    exit /b 1
)

for /f "tokens=2" %%i in ('python --version 2^>^&1') do set python_version=%%i
echo [✅ OK] Python found: %python_version%

REM Check if Git is installed (optional for cloning, but good to have)
git --version >nul 2>&1
if errorlevel 1 (
    echo [⚠️  WARNING] Git not found (optional)
    echo    Git is recommended for updates but not required for basic usage
) else (
    for /f "tokens=3" %%i in ('git --version 2^>^&1') do set git_version=%%i
    echo [✅ OK] Git found: %git_version%
)

REM Create virtual environment
echo.
echo 🔧 Setting up Python environment...
if exist "mutationscan_env" (
    echo [⚠️  WARNING] Virtual environment already exists
    echo    Removing old environment...
    rmdir /s /q mutationscan_env
)

python -m venv mutationscan_env
if errorlevel 1 (
    echo [❌ ERROR] Failed to create virtual environment
    echo    Please check your Python installation
    pause
    exit /b 1
)

echo [✅ OK] Virtual environment created

REM Activate virtual environment
call mutationscan_env\Scripts\activate.bat
if errorlevel 1 (
    echo [❌ ERROR] Could not activate virtual environment
    pause
    exit /b 1
)

REM Upgrade pip
echo.
echo 📦 Upgrading pip...
python -m pip install --upgrade pip --quiet
if errorlevel 1 (
    echo [❌ ERROR] Failed to upgrade pip
    pause
    exit /b 1
)

REM Install dependencies
echo.
echo 📚 Installing Python dependencies...
echo    This may take a few minutes...
pip install -r requirements.txt
if errorlevel 1 (
    echo [❌ ERROR] Failed to install dependencies
    echo    Please check your internet connection and try again
    pause
    exit /b 1
)

echo [✅ OK] Dependencies installed successfully

REM Install the subscan package in development mode
echo.
echo 🔧 Installing MutationScan package...
cd subscan
pip install -e .
if errorlevel 1 (
    echo [❌ ERROR] Failed to install MutationScan package
    pause
    exit /b 1
)

cd ..

REM Verify installation
echo.
echo 🧪 Verifying installation...
python -c "import subscan; import tqdm; import pandas; import biopython; print('✅ All packages imported successfully')" 2>nul
if errorlevel 1 (
    echo [⚠️  WARNING] Some packages may not be working correctly
    echo    The pipeline should still function, but you may encounter issues
) else (
    echo [✅ OK] Installation verified successfully
)

REM Test the pipeline help
echo.
echo 🧪 Testing pipeline...
python subscan\tools\run_pipeline.py --help >nul 2>&1
if errorlevel 1 (
    echo [❌ ERROR] Pipeline test failed
    echo    Please check the installation and try again
    pause
    exit /b 1
)

echo [✅ OK] Pipeline test passed

REM Create example data if it doesn't exist
if not exist "examples\demo_accessions.txt" (
    echo.
    echo 📁 Creating example data...
    echo NC_013791 > examples\demo_accessions.txt
    echo acrB > examples\gene_list.txt
    echo acrA >> examples\gene_list.txt
    echo mdtH >> examples\gene_list.txt
    echo [✅ OK] Example data created
)

REM Create data directories
if not exist "data_input" mkdir data_input
if not exist "data_output" mkdir data_output

echo.
echo 🎉 Installation Complete!
echo ==========================
echo.
echo MutationScan is ready to use!
echo.
echo 💡 Quick Start:
echo    PowerShell users: .\mutationscan.ps1 --help
echo    Command Prompt users: .\mutationscan.bat --help
echo.
echo 📚 Example analysis:
echo    .\mutationscan.ps1 --accessions examples\demo_accessions.txt --gene-list examples\gene_list.txt --email your@email.com --output-dir results --sepi-species "Escherichia coli"
echo.
echo 📖 Documentation: README.md
echo 🐛 Issues: https://github.com/vihaankulkarni29/MutationScan/issues
echo.
echo Thank you for using MutationScan!
echo.
pause