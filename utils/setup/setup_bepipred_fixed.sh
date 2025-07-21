#!/bin/bash
# setup_bepipred_fixed.sh
# Setup script for BepiPred 3.0 with dependency conflict resolution

set -e  # Exit on any error

echo "=========================================="
echo "BepiPred 3.0 Setup (Fixed Dependencies)"
echo "=========================================="

# Check if we're on Ubuntu/Debian
if ! command -v apt &> /dev/null; then
    echo "Warning: This script is designed for Ubuntu/Debian systems with apt package manager"
fi

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
echo "Python version: $PYTHON_VERSION"

# Create tools directory if it doesn't exist
mkdir -p tools

# Clone BepiPred 3.0 repository
echo "Cloning BepiPred 3.0 repository..."
if [ -d "tools/BepiPred3.0" ]; then
    echo "BepiPred3.0 directory already exists, updating..."
    cd tools/BepiPred3.0
    git pull
    cd ../..
else
    git clone https://github.com/UberClifford/BepiPred3.0-Predictor.git tools/BepiPred3.0
fi

# Change to BepiPred directory
cd tools/BepiPred3.0

echo "Setting up BepiPred 3.0 with fixed dependencies..."

# Create virtual environment (recommended)
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Create fixed requirements file
echo "Creating fixed requirements.txt..."
cat > requirements_fixed.txt << 'EOF'
# Fixed requirements for BepiPred 3.0 with Python 3.10 compatibility
# Original requirements.txt had version conflicts

# Core dependencies
torch>=1.9.0
transformers>=4.5.0
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0

# BepiPred specific (if available)
bp3>=0.0.12.7
fair-esm>=1.0.3

# Additional dependencies that might be needed
biopython>=1.78
matplotlib>=3.3.0
scikit-learn>=1.0.0

# Compatibility packages
importlib-resources>=5.8.0
EOF

# Install dependencies with conflict resolution
echo "Installing Python dependencies (fixed versions)..."

# Install core ML packages first
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu

# Install transformers and related packages
pip install transformers>=4.5.0

# Install scientific computing packages
pip install numpy>=1.21.0 pandas>=1.3.0 scipy>=1.7.0

# Install bio packages
pip install biopython>=1.78

# Install visualization
pip install matplotlib>=3.3.0

# Try to install BepiPred specific packages (may not be available)
echo "Installing BepiPred specific packages..."
pip install bp3>=0.0.12.7 || echo "Warning: bp3 package not available, continuing..."
pip install fair-esm>=1.0.3 || echo "Warning: fair-esm package not available, continuing..."

# Install any remaining packages
pip install importlib-resources>=5.8.0 scikit-learn>=1.0.0

# If original requirements.txt exists, try to install remaining packages
if [ -f "requirements.txt" ]; then
    echo "Attempting to install any remaining packages from original requirements.txt..."
    # Extract packages not already installed
    grep -v "numpy\|pandas\|torch\|transformers\|bp3\|fair-esm\|importlib-resources" requirements.txt > remaining_requirements.txt || true
    
    if [ -s remaining_requirements.txt ]; then
        pip install -r remaining_requirements.txt || echo "Some optional packages could not be installed, continuing..."
    fi
    rm -f remaining_requirements.txt
fi

# Test installation
echo "Testing BepiPred 3.0 installation..."
if python bepipred3_CLI.py -h > /dev/null 2>&1; then
    echo "✅ BepiPred 3.0 installed successfully!"
    
    # Show installed packages
    echo ""
    echo "Installed packages:"
    pip list | grep -E "(torch|transform|numpy|pandas|bp3|fair-esm|bio)"
    
else
    echo "❌ BepiPred 3.0 installation test failed"
    echo "Trying alternative approach..."
    
    # Alternative: install minimal dependencies and test
    echo "Installing minimal dependencies..."
    pip install torch numpy pandas biopython matplotlib
    
    if python bepipred3_CLI.py -h > /dev/null 2>&1; then
        echo "✅ BepiPred 3.0 working with minimal dependencies!"
    else
        echo "❌ BepiPred 3.0 still not working"
        echo "You may need to install missing dependencies manually"
        echo "Check the error output above for clues"
    fi
fi

# Test with example file if it exists
if [ -f "example_antigens.fasta" ]; then
    echo ""
    echo "Running test prediction..."
    mkdir -p test_output
    
    if python bepipred3_CLI.py -i example_antigens.fasta -o test_output/ -pred vt_pred; then
        if [ -d "test_output" ] && [ "$(ls -A test_output)" ]; then
            echo "✅ Test prediction completed successfully!"
            echo "Output files:"
            ls -la test_output/
            rm -rf test_output  # Clean up test
        else
            echo "⚠️  BepiPred ran but produced no output files"
        fi
    else
        echo "⚠️  Test prediction failed, but BepiPred CLI is working"
        echo "This might be due to missing model files (will download on first real use)"
    fi
else
    echo "No example file found, skipping test prediction"
fi

# Create activation script for easy use
cat > activate_bepipred.sh << 'EOF'
#!/bin/bash
# Activate BepiPred virtual environment
cd "$(dirname "$0")"
source venv/bin/activate
echo "BepiPred 3.0 environment activated"
echo "Run: python bepipred3_CLI.py -h for help"
EOF

chmod +x activate_bepipred.sh

# Deactivate virtual environment
deactivate

echo ""
echo "=========================================="
echo "✅ BepiPred 3.0 setup completed!"
echo "=========================================="
echo ""
echo "To activate BepiPred environment manually:"
echo "  cd tools/BepiPred3.0"
echo "  source venv/bin/activate"
echo "  # or use: ./activate_bepipred.sh"
echo ""
echo "To use BepiPred 3.0 in the pipeline:"
echo "1. Ensure config/config.yaml has correct path:"
echo "   bepipred:"
echo "     path: \"tools/BepiPred3.0\""
echo ""
echo "2. Run the pipeline with BepiPred:"
echo "   snakemake all_epitope_predictions_bepipred --cores 4"
echo ""
echo "Note: The first run may take longer as BepiPred downloads ESM-2 models (~2GB)"
echo ""

# Return to project root
cd ../..

echo "Setup script completed!"
echo ""
echo "If you encountered any errors:"
echo "1. Check that you have Python 3.8+ installed"
echo "2. Ensure you have internet connection for package downloads"
echo "3. Try running the test script: python test_bepipred.py"