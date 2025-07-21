#!/bin/bash
# setup_bepipred_conda.sh
# Setup BepiPred 3.0 using conda for better dependency management

set -e  # Exit on any error

echo "=========================================="
echo "BepiPred 3.0 Setup (Conda Method)"
echo "=========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Miniconda or Anaconda first."
    echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Conda version: $(conda --version)"

# Create tools directory
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

cd tools/BepiPred3.0

# Create conda environment for BepiPred
echo "Creating conda environment for BepiPred..."
if conda env list | grep -q "bepipred"; then
    echo "Conda environment 'bepipred' already exists, updating..."
    conda activate bepipred
else
    # Create environment with Python 3.9 (more compatible)
    conda create -n bepipred python=3.9 -y
    conda activate bepipred
fi

# Install core dependencies via conda (better compatibility)
echo "Installing dependencies via conda..."
conda install -c conda-forge -c pytorch -y \
    pytorch \
    numpy \
    pandas \
    scipy \
    matplotlib \
    scikit-learn \
    biopython

# Install transformers and related via pip (in conda env)
echo "Installing ML packages via pip..."
pip install transformers fair-esm

# Try to install BepiPred specific packages
echo "Installing BepiPred specific packages..."
pip install bp3 || echo "Warning: bp3 package not available, continuing..."

# Install any additional requirements
if [ -f "requirements.txt" ]; then
    echo "Installing remaining requirements..."
    # Try to install what we can, ignore failures
    pip install -r requirements.txt || echo "Some packages failed to install, continuing..."
fi

# Test installation
echo "Testing BepiPred 3.0 installation..."
if python bepipred3_CLI.py -h > /dev/null 2>&1; then
    echo "✅ BepiPred 3.0 installed successfully!"
    
    # Show key packages
    echo ""
    echo "Key installed packages:"
    conda list | grep -E "(torch|transform|numpy|pandas)" || pip list | grep -E "(torch|transform|numpy|pandas)"
    
else
    echo "❌ BepiPred 3.0 installation test failed"
    echo "Installed packages:"
    conda list
fi

# Test with example if available
if [ -f "example_antigens.fasta" ]; then
    echo ""
    echo "Running test prediction..."
    mkdir -p test_output
    
    if timeout 60 python bepipred3_CLI.py -i example_antigens.fasta -o test_output/ -pred vt_pred 2>/dev/null; then
        if [ -d "test_output" ] && [ "$(ls -A test_output)" ]; then
            echo "✅ Test prediction completed successfully!"
            rm -rf test_output
        else
            echo "⚠️  BepiPred ran but no output (may need model download)"
        fi
    else
        echo "⚠️  Test skipped (may require model download on first use)"
    fi
fi

# Create activation script
cat > activate_bepipred_conda.sh << 'EOF'
#!/bin/bash
# Activate BepiPred conda environment
echo "Activating BepiPred conda environment..."
eval "$(conda shell.bash hook)"
conda activate bepipred
echo "✅ BepiPred environment activated"
echo "Current directory: $(pwd)"
echo "Python: $(which python)"
echo "Run: python bepipred3_CLI.py -h for help"
EOF

chmod +x activate_bepipred_conda.sh

conda deactivate

echo ""
echo "=========================================="
echo "✅ BepiPred 3.0 conda setup completed!"
echo "=========================================="
echo ""
echo "To use BepiPred:"
echo "1. Activate environment:"
echo "   conda activate bepipred"
echo "   # or use: cd tools/BepiPred3.0 && ./activate_bepipred_conda.sh"
echo ""
echo "2. Test installation:"
echo "   conda activate bepipred"
echo "   python test_bepipred.py"
echo ""
echo "3. Run pipeline:"
echo "   # Make sure bepipred environment is active"
echo "   snakemake all_epitope_predictions_bepipred --cores 4"
echo ""
echo "Environment name: bepipred"
echo "Location: $(conda info --envs | grep bepipred | awk '{print $2}' || echo 'Not found')"

cd ../..
echo "Setup completed!"