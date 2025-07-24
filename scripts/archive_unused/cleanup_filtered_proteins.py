#!/usr/bin/env python3
"""
Clean up data for proteins that were filtered out by surface accessibility
"""

import os
import shutil
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define filtered proteins for each group
filtered_proteins = {
    'positive': ['cls', 'dnaK', 'era', 'glyS', 'pstS'],
    'negative': ['bamD', 'bamE', 'cls', 'comL', 'copB', 'dnaK', 'era', 
                 'flgF', 'flgG', 'fliI', 'fliL', 'glyS', 'lptA', 'omlA', 
                 'pstS', 'yfiO']
}

def cleanup_protein_directories(base_dir, analysis, paramset, group):
    """Remove directories for filtered proteins"""
    
    proteins_to_remove = filtered_proteins[group]
    removed_count = 0
    
    for protein in proteins_to_remove:
        protein_dir = base_dir / protein
        if protein_dir.exists() and protein_dir.is_dir():
            logging.info(f"Removing {protein_dir}")
            shutil.rmtree(protein_dir)
            removed_count += 1
    
    return removed_count

def main():
    # Define analysis and paramset
    analysis = "analysis_1"
    paramset = "params_1"
    
    # Define directories to clean
    directories_to_clean = [
        "results/protein_fasta",
        "results/proteins_to_download",
        "results/msa_sequences", 
        "results/msa_alignments",
        "results/msa_trimmed",
        "results/msa_quality",
        "results/conservation",
        "results/epitope_predictions",
        "results/epitope_predictions_bepipred",
        "results/3d_structures"
    ]
    
    # Clean each directory
    for dir_path in directories_to_clean:
        logging.info(f"\nCleaning {dir_path}...")
        
        for group in ['positive', 'negative']:
            group_dir = Path(dir_path) / f"{analysis}_{paramset}_gram_{group}"
            
            if group_dir.exists():
                removed = cleanup_protein_directories(group_dir, analysis, paramset, group)
                logging.info(f"  Removed {removed} protein directories from gram-{group}")
            else:
                logging.info(f"  Directory {group_dir} does not exist, skipping")
    
    logging.info("\nCleanup complete!")

if __name__ == "__main__":
    main()