# 📹 MutationScan Installation Video Guide

## Quick Installation Video Tutorial

> **Coming Soon:** Video installation guide for non-technical users

### For Now: Step-by-Step Installation with Screenshots

## 🖥️ Windows Installation

### Step 1: Install Python
1. Download Python 3.9+ from [python.org](https://python.org)
2. **IMPORTANT:** Check "Add Python to PATH" during installation
3. Verify installation:
   ```cmd
   python --version
   ```

### Step 2: Install Git (Optional but Recommended)
1. Download Git from [git-scm.com](https://git-scm.com)
2. Use default installation options
3. Verify installation:
   ```cmd
   git --version
   ```

### Step 3: Download MutationScan
**Option A: Using Git (Recommended)**
```cmd
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan
```

**Option B: Download ZIP**
1. Go to the GitHub repository
2. Click "Code" → "Download ZIP"
3. Extract to desired location
4. Open Command Prompt and navigate to folder

### Step 4: Install MutationScan
```cmd
cd subscan
pip install -e .
```

### Step 5: Verify Installation
```cmd
mutationscan --help
```

## 🐧 Linux/macOS Installation

### Prerequisites
```bash
# Ubuntu/Debian
sudo apt update
sudo apt install python3 python3-pip git

# macOS (with Homebrew)
brew install python git
```

### Installation Steps
```bash
git clone https://github.com/vihaankulkarni29/MutationScan.git
cd MutationScan/subscan
pip3 install -e .
```

### Verify Installation
```bash
mutationscan --help
```

## 🆘 Troubleshooting Installation

### Common Issues

#### "Python is not recognized"
- **Solution:** Add Python to PATH or reinstall Python with "Add to PATH" checked

#### "pip is not recognized"
- **Solution:** 
  ```cmd
  python -m pip install -e .
  ```

#### Permission denied errors
- **Windows:** Run Command Prompt as Administrator
- **Linux/macOS:** Use `sudo` or virtual environment

#### Package installation failures
- **Solution:** Update pip first:
  ```bash
  python -m pip install --upgrade pip
  ```

### Virtual Environment (Recommended for Advanced Users)
```bash
# Create virtual environment
python -m venv mutationscan_env

# Activate (Windows)
mutationscan_env\Scripts\activate

# Activate (Linux/macOS)
source mutationscan_env/bin/activate

# Install MutationScan
cd subscan
pip install -e .
```

## 📧 Need Help?

If you encounter installation issues:
1. Check the [main troubleshooting guide](../README.md#troubleshooting)
2. Ensure you have Python 3.9+
3. Try using a virtual environment
4. Create an issue on GitHub with:
   - Your operating system
   - Python version (`python --version`)
   - Error message (full text)

## 🎯 What's Next?

After successful installation:
1. Try the [Quick Start guide](../README.md#quick-start)
2. Review [example workflows](../examples/README.md)
3. Read the [complete documentation](README.md)

---

**Installation Video Coming Soon!** 📹  
We're working on a comprehensive video tutorial for non-technical users.