@echo off
REM MutationScan One-Click Installer for Windows
REM Usage: Save as install.bat and double-click

echo MutationScan One-Click Installer for Windows
echo ================================================

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found. Please install Python 3.9+ from https://python.org
    echo    Make sure to check "Add Python to PATH" during installation
    pause
    exit /b 1
)

echo [OK] Python found

REM Check if Git is installed
git --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Git not found. Please install Git from https://git-scm.com
    pause
    exit /b 1
)

echo [OK] Git found

REM Create virtual environment
echo Setting up Python environment...
python -m venv mutationscan_env
call mutationscan_env\Scripts\activate.bat

REM Upgrade pip
python -m pip install --upgrade pip

REM Clone and install MutationScan
echo Installing MutationScan...
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan\subscan
pip install -e .

REM Verify installation
python -c "import subscan; print('[OK] MutationScan installed successfully!')"

REM Create batch launcher
echo Creating launcher...
cd ..

REM Create batch file launcher
echo @echo off > mutationscan.bat
echo REM MutationScan Launcher >> mutationscan.bat
echo echo MutationScan Pipeline Launcher >> mutationscan.bat
echo echo ================================ >> mutationscan.bat
echo echo. >> mutationscan.bat
echo cd /d "%%~dp0MutationScan" >> mutationscan.bat
echo if errorlevel 1 ( >> mutationscan.bat
echo     echo ERROR: Could not find MutationScan directory >> mutationscan.bat
echo     pause >> mutationscan.bat
echo     exit /b 1 >> mutationscan.bat
echo ^) >> mutationscan.bat
echo call mutationscan_env\Scripts\activate.bat >> mutationscan.bat
echo if errorlevel 1 ( >> mutationscan.bat
echo     echo ERROR: Could not activate virtual environment >> mutationscan.bat
echo     pause >> mutationscan.bat
echo     exit /b 1 >> mutationscan.bat
echo ^) >> mutationscan.bat
echo if "%%~1"=="" ( >> mutationscan.bat
echo     python subscan\tools\run_pipeline.py --help >> mutationscan.bat
echo     echo. >> mutationscan.bat
echo     echo Example: .\mutationscan.bat --accessions file.txt --genes genes.txt --output results\ >> mutationscan.bat
echo     pause >> mutationscan.bat
echo ^) else ( >> mutationscan.bat
echo     python subscan\tools\run_pipeline.py %%* >> mutationscan.bat
echo     pause >> mutationscan.bat
echo ^) >> mutationscan.bat

REM Create PowerShell launcher (already exists as mutationscan.ps1)

echo.
echo Installation Complete!
echo =========================
echo.
echo MutationScan is ready to use!
echo.
echo Quick Start:
echo    PowerShell users: .\mutationscan.ps1 --help
echo    Command Prompt users: .\mutationscan.bat --help
echo    
echo    Example analysis:
echo    .\mutationscan.ps1 --accessions file.txt --genes genes.txt --output results\
echo.
echo Documentation: https://github.com/vihaankulkarni29/MutationScan
echo Issues: https://github.com/vihaankulkarni29/MutationScan/issues
echo.
echo Thank you for using MutationScan!
pause