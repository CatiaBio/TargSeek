#!/usr/bin/env python3
"""
Create Sequence References for MSA
==================================

This script creates reference files (filelists) for both regular sequences and 3D structures.
It generates:
- {gene}_filelist.txt: Regular sequences from data/proteins_fasta/
- {gene}_3d_filelist.txt: 3D structure sequences from data/proteins_3d_structure/
"""

import sys
import logging
from pathlib import Path
from select_msa_proteins import (
    load_gene_species_lists,
    find_available_sequences,
    create_3d_structure_filelist
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    try:
        gene_lists_dir = Path(snakemake.input.gene_lists)
        output_dir = Path(snakemake.output[0])
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
    except NameError:
        gene_lists_dir = Path(sys.argv[1])
        output_dir = Path(sys.argv[2])
        analysis = "analysis_1"
        paramset = "params_1"
        group = "gram_negative"

    output_dir.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"=== Creating Sequence References ===")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    logging.info(f"Gene lists directory: {gene_lists_dir}")
    logging.info(f"Output directory: {output_dir}")

    gene_species_mapping = load_gene_species_lists(gene_lists_dir)
    
    if not gene_species_mapping:
        logging.error("No gene species lists found")
        return
    
    logging.info(f"Processing {len(gene_species_mapping)} genes")

    for gene_idx, (gene_name, species_list) in enumerate(gene_species_mapping.items(), 1):
        logging.info(f"\n[{gene_idx}/{len(gene_species_mapping)}] Processing gene: {gene_name}")
        
        # Find available sequences (use correct data directory)
        protein_sequences_dir = "data/protein_sequences"  # From config
        sequences = find_available_sequences(gene_name, species_list, protein_sequences_dir)
        logging.info(f"  Found {len(sequences)} sequences for {gene_name}")

        # Create gene directory
        gene_dir = output_dir / gene_name
        gene_dir.mkdir(parents=True, exist_ok=True)

        # Write regular sequences filelist
        filelist_path = gene_dir / f"{gene_name}_filelist.txt"
        with open(filelist_path, 'w') as f:
            for species, fasta_path in sequences.items():
                f.write(f"{fasta_path}\n")
        logging.info(f"  Created filelist: {len(sequences)} sequences")

        # Create 3D structure filelist
        # Note: We pass empty sequence_references since we're only creating the 3D filelist
        protein_structures_dir = "data/protein_structures"  # From config
        create_3d_structure_filelist(gene_name, gene_dir, [], protein_structures_dir)
    
    logging.info(f"\n=== Sequence reference creation completed successfully! ===")

if __name__ == "__main__":
    main()
