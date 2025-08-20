#!/bin/bash
# TargSeek Pipeline Activation Script
# This script sets up the proper environment and library paths for running the TargSeek pipeline

# Activate targseek conda environment
source /home/catia/tools/miniconda3/bin/activate targseek

# Set library path for BLAST tools
export LD_LIBRARY_PATH=/home/catia/tools/miniconda3/envs/targseek/lib:$LD_LIBRARY_PATH

# Set PATH for ConSurf tools
export PATH=/home/catia/tools/rostlab_consurf:/home/catia/tools/stand_alone_consurf/stand_alone_consurf-1.00:$PATH

echo "✅ TargSeek pipeline environment activated!"
echo "📊 Available tools:"
echo "   - Snakemake: $(which snakemake)"
echo "   - MAFFT: $(which mafft)"
echo "   - ClipKIT: $(which clipkit)"
echo "   - PyMOL: $(which pymol)"
echo "   - BLAST: /home/catia/tools/ncbi-blast-2.17.0+/bin/blastp"
echo "   - ConSurf: /home/catia/tools/stand_alone_consurf/stand_alone_consurf-1.00/stand_alone_consurf.py"
echo "   - BepiPred: /home/catia/tools/BepiPred3.0-Predictor/bepipred3_CLI.py (use 'conda activate bepipred' first)"
echo ""
echo "🚀 To run the analysis pipeline:"
echo "   snakemake -s Snakefile_analysis all_analysis --cores 4"
echo ""