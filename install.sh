#!/bin/bash
# MutationScan One-Click Installer
# Usage: curl -sSL https://raw.githubusercontent.com/vihaankulkarni29/MutationScan/main/install.sh | bash

set -e

echo "🧬 MutationScan One-Click Installer"
echo "==================================="

# Detect OS
OS=$(uname -s)
ARCH=$(uname -m)

echo "🖥️  Detected: $OS $ARCH"

# Install dependencies based on OS
install_system_deps() {
    echo "📦 Installing system dependencies..."
    
    if [[ "$OS" == "Linux" ]]; then
        if command -v apt-get &> /dev/null; then
            sudo apt-get update -qq
            sudo apt-get install -y python3 python3-pip python3-venv git abricate
        elif command -v yum &> /dev/null; then
            sudo yum install -y python3 python3-pip git
        elif command -v conda &> /dev/null; then
            conda install -c bioconda abricate -y
        fi
    elif [[ "$OS" == "Darwin" ]]; then
        if command -v brew &> /dev/null; then
            brew install python3 git
            brew install brewsci/bio/abricate
        fi
    fi
}

# Install Python and create environment
setup_python_env() {
    echo "🐍 Setting up Python environment..."
    
    # Create virtual environment
    python3 -m venv mutationscan_env
    source mutationscan_env/bin/activate
    
    # Upgrade pip
    pip install --upgrade pip
}

# Install MutationScan
install_mutationscan() {
    echo "🧬 Installing MutationScan..."
    
    # Clone repository
    git clone https://github.com/vihaankulkarni29/MutationScan.git
    cd MutationScan/subscan
    
    # Install the package
    pip install -e .
    
    # Verify installation
    python -c "import subscan; print('✅ MutationScan installed successfully!')"
}

# Create launcher script
create_launcher() {
    echo "🚀 Creating launcher script..."
    
    cat > ~/mutationscan << 'EOF'
#!/bin/bash
# MutationScan Launcher
cd ~/MutationScan
source mutationscan_env/bin/activate
python subscan/tools/run_pipeline.py "$@"
EOF
    
    chmod +x ~/mutationscan
    
    # Add to PATH if not already there
    if [[ ":$PATH:" != *":$HOME:"* ]]; then
        echo 'export PATH="$HOME:$PATH"' >> ~/.bashrc
        echo 'export PATH="$HOME:$PATH"' >> ~/.zshrc 2>/dev/null || true
    fi
}

# Main installation
main() {
    echo "🚀 Starting installation..."
    
    # Check if already installed
    if [ -d "MutationScan" ] && [ -f "mutationscan_env/bin/activate" ]; then
        echo "✅ MutationScan appears to be already installed!"
        echo "   To reinstall, delete the MutationScan directory and run again."
        exit 0
    fi
    
    install_system_deps
    setup_python_env
    install_mutationscan
    create_launcher
    
    echo ""
    echo "🎉 Installation Complete!"
    echo "========================="
    echo ""
    echo "✅ MutationScan is ready to use!"
    echo ""
    echo "🚀 Quick Start:"
    echo "   1. Restart your terminal (or run: source ~/.bashrc)"
    echo "   2. Run: mutationscan --help"
    echo "   3. Analyze genomes: mutationscan --accessions file.txt --genes genes.txt --output results/"
    echo ""
    echo "📚 Documentation: https://github.com/vihaankulkarni29/MutationScan"
    echo "🐛 Issues: https://github.com/vihaankulkarni29/MutationScan/issues"
    echo ""
    echo "Thank you for using MutationScan! 🧬"
}

# Run main function
main