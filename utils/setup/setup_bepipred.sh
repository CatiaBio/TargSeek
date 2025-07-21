#!/bin/bash
# setup_bepipred.sh
# Setup script for BepiPred 3.0 installation on Ubuntu

set -e  # Exit on any error

echo "=========================================="
echo "BepiPred 3.0 Setup for PureMilk Pipeline"
echo "=========================================="

# Check if we're on Ubuntu/Debian
if ! command -v apt &> /dev/null; then
    echo "Warning: This script is designed for Ubuntu/Debian systems with apt package manager"
fi

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
echo "Python version: $PYTHON_VERSION"

if [ "$(printf '%s\n' "3.8" "$PYTHON_VERSION" | sort -V | head -n1)" != "3.8" ]; then
    echo "Warning: BepiPred 3.0 recommends Python 3.8.8, you have $PYTHON_VERSION"
fi

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

echo "Installing BepiPred 3.0 dependencies..."

# Create virtual environment (recommended)
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies from requirements.txt
echo "Installing Python dependencies..."
pip install -r requirements.txt

# Test installation
echo "Testing BepiPred 3.0 installation..."
if python bepipred3_CLI.py -h > /dev/null 2>&1; then
    echo "✅ BepiPred 3.0 installed successfully!"
else
    echo "❌ BepiPred 3.0 installation test failed"
    exit 1
fi

# Test with example file if it exists
if [ -f "example_antigens.fasta" ]; then
    echo "Running test prediction..."
    mkdir -p test_output
    python bepipred3_CLI.py -i example_antigens.fasta -o test_output/ -pred vt_pred
    
    if [ -d "test_output" ] && [ "$(ls -A test_output)" ]; then
        echo "✅ Test prediction completed successfully!"
        echo "Output files:"
        ls -la test_output/
        rm -rf test_output  # Clean up test
    else
        echo "❌ Test prediction failed"
        exit 1
    fi
fi

# Deactivate virtual environment
deactivate

echo ""
echo "=========================================="
echo "✅ BepiPred 3.0 setup completed!"
echo "=========================================="
echo ""
echo "To use BepiPred 3.0 in the pipeline:"
echo "1. Update config/config.yaml if needed:"
echo "   bepipred:"
echo "     path: \"tools/BepiPred3.0\""
echo ""
echo "2. Run the pipeline with BepiPred:"
echo "   snakemake all_epitope_predictions_bepipred --cores 4"
echo ""
echo "3. Or run for specific group:"
echo "   snakemake results/epitope_predictions_bepipred/analysis_1_params_1_gram_positive --cores 4"
echo ""
echo "Note: The first run may take longer as BepiPred downloads ESM-2 models"
echo ""

# Return to project root
cd ../..

echo "Setup script completed!"