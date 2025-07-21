#!/usr/bin/env python3
"""
Test script for BepiPred 3.0 epitope prediction
Runs BepiPred on a sample gene from your existing MSA sequences
"""

import sys
from pathlib import Path
sys.path.append('scripts')

from predict_epitopes_bepipred import predict_epitopes_for_gene
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_bepipred():
    """Test BepiPred on existing data"""
    
    # Test parameters
    gene = "dnaK"  # Use dnaK as test gene
    msa_sequences_dir = "results/msa_sequences/analysis_1_params_1_gram_positive"
    conservation_dir = "results/conservation/analysis_1_params_1_gram_positive"  
    structures_dir = "results/3d_structures/analysis_1_params_1_gram_positive"
    output_dir = "test_epitope_predictions_bepipred"
    bepipred_path = "tools/BepiPred3.0"
    
    print("Testing BepiPred 3.0 epitope prediction...")
    print(f"Gene: {gene}")
    print(f"BepiPred path: {bepipred_path}")
    print()
    
    # Check if required directories exist
    required_dirs = [msa_sequences_dir, conservation_dir, structures_dir]
    for dir_path in required_dirs:
        if not Path(dir_path).exists():
            print(f"❌ Required directory not found: {dir_path}")
            print("Please ensure you have run the pipeline up to conservation analysis")
            return False
    
    # Check if BepiPred is installed
    bepipred_script = Path(bepipred_path) / "bepipred3_CLI.py"
    if not bepipred_script.exists():
        print(f"❌ BepiPred script not found: {bepipred_script}")
        print("Please run ./setup_bepipred.sh first")
        return False
    
    # Check if test gene exists
    gene_msa_file = Path(msa_sequences_dir) / f"{gene}.fasta"
    if not gene_msa_file.exists():
        print(f"❌ MSA file not found for {gene}: {gene_msa_file}")
        print("Available genes:")
        msa_files = Path(msa_sequences_dir).glob("*.fasta")
        for msa_file in sorted(msa_files):
            print(f"  - {msa_file.stem}")
        return False
    
    try:
        # Run epitope prediction
        summary = predict_epitopes_for_gene(
            gene, msa_sequences_dir, conservation_dir, structures_dir, output_dir, bepipred_path
        )
        
        if summary:
            print("✅ BepiPred prediction successful!")
            print(f"Results summary for {gene}:")
            print(f"  - Total epitopes: {summary['total_epitopes']}")
            print(f"  - High-score epitopes: {summary['high_score_epitopes']}")
            print(f"  - Average score: {summary['avg_score']:.3f}")
            print(f"  - Sequence length: {summary['sequence_length']}")
            print(f"  - Results saved to: {output_dir}")
            return True
        else:
            print("❌ BepiPred prediction failed - no results returned")
            return False
            
    except Exception as e:
        print(f"❌ BepiPred prediction failed with error: {e}")
        return False

if __name__ == "__main__":
    success = test_bepipred()
    sys.exit(0 if success else 1)