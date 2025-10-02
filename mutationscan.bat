@echo off
REM MutationScan Launcher - Enhanced with error handling
echo 🧬 MutationScan Pipeline Launcher
echo ================================

REM Get the directory where this batch file is located
set SCRIPT_DIR=%~dp0

REM Check if MutationScan directory exists
if not exist "%SCRIPT_DIR%MutationScan" (
    echo ❌ ERROR: MutationScan directory not found at %SCRIPT_DIR%MutationScan
    echo    Make sure you're running this from the correct directory
    echo    Expected structure:
    echo      📁 Current Directory
    echo      ├── mutationscan.bat (this file)
    echo      └── 📁 MutationScan/
    pause
    exit /b 1
)

REM Change to MutationScan directory
cd /d "%SCRIPT_DIR%MutationScan"
if errorlevel 1 (
    echo ❌ ERROR: Could not change to MutationScan directory
    pause
    exit /b 1
)

echo ✅ Found MutationScan directory: %CD%

REM Check if virtual environment exists
if not exist "mutationscan_env\Scripts\activate.bat" (
    echo ❌ ERROR: Virtual environment not found
    echo    Expected: mutationscan_env\Scripts\activate.bat
    echo    Please run the installer again: install.bat
    pause
    exit /b 1
)

echo ✅ Found virtual environment

REM Activate virtual environment
echo 🐍 Activating Python environment...
call mutationscan_env\Scripts\activate.bat
if errorlevel 1 (
    echo ❌ ERROR: Could not activate virtual environment
    pause
    exit /b 1
)

echo ✅ Virtual environment activated

REM Check if subscan directory exists
if not exist "subscan\tools\run_pipeline.py" (
    echo ❌ ERROR: Pipeline script not found
    echo    Expected: subscan\tools\run_pipeline.py
    echo    Current directory: %CD%
    dir subscan\tools\ 2>nul
    pause
    exit /b 1
)

echo ✅ Found pipeline script

REM Show help if no arguments provided
if "%~1"=="" (
    echo 📚 Showing help - add your arguments to run analysis
    echo.
    python subscan\tools\run_pipeline.py --help
    echo.
    echo 💡 Example usage:
    echo    mutationscan.bat --accessions data_input\accessions.txt --genes data_input\genes.txt --output results
    echo.
    echo Press any key to close...
    pause >nul
    exit /b 0
)

REM Run the pipeline with arguments
echo 🚀 Running MutationScan pipeline...
echo Command: python subscan\tools\run_pipeline.py %*
echo.

python subscan\tools\run_pipeline.py %*

REM Check exit code
if errorlevel 1 (
    echo.
    echo ❌ Pipeline finished with errors (exit code: %errorlevel%)
    echo    Check the error messages above for details
) else (
    echo.
    echo ✅ Pipeline completed successfully!
)

echo.
echo Press any key to close...
pause >nul
