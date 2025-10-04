@echo off
REM MutationScan Debug Launcher - Shows all details
echo MutationScan Debug Launcher
echo =============================

echo System Information:
echo Current Directory: %CD%
echo Script Directory: %~dp0
echo Python Version:
python --version 2>nul || echo ERROR: Python not found in PATH

echo.
echo Directory Structure Check:
echo Looking for MutationScan folder...
if exist "MutationScan" (
    echo [OK] MutationScan folder found
    dir MutationScan /w
) else (
    echo [ERROR] MutationScan folder NOT found
    echo Current directory contents:
    dir /w
)

echo.
echo Virtual Environment Check:
if exist "MutationScan\mutationscan_env" (
    echo [OK] Virtual environment folder found
    dir "MutationScan\mutationscan_env" /w
) else (
    echo [ERROR] Virtual environment NOT found
)

echo.
echo Pipeline Script Check:
if exist "MutationScan\subscan\tools\run_pipeline.py" (
    echo [OK] Pipeline script found
) else (
    echo [ERROR] Pipeline script NOT found
    if exist "MutationScan\subscan" (
        echo Subscan directory contents:
        dir "MutationScan\subscan" /w
    )
)

echo.
echo Testing Python Import:
cd /d "MutationScan" 2>nul
if exist "mutationscan_env\Scripts\activate.bat" (
    call mutationscan_env\Scripts\activate.bat
    echo Testing subscan import...
    python -c "import subscan; print('[OK] SubScan import successful')" 2>nul || echo "[ERROR] SubScan import failed"
    python -c "from subscan.domino_tools.ncbi_extractor import NCBIExtractor; print('[OK] NCBI Extractor import successful')" 2>nul || echo "[ERROR] NCBI Extractor import failed"
) else (
    echo [ERROR] Cannot test - virtual environment not found
)

echo.
echo Troubleshooting Tips:
echo 1. Make sure you ran install.bat successfully
echo 2. Make sure you're in the correct directory (where install.bat created files)
echo 3. Try running install.bat again if something is missing
echo 4. Check if antivirus blocked any files

echo.
echo Press any key to close...
pause >nul