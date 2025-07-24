#!/usr/bin/env python3
"""
Create Gene-Specific Species Lists
=================================

This script creates individual species lists for each gene based on the 
proteins_to_study TSV files. This helps with organization and allows
for easier debugging and monitoring of download progress.

Input: results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv
Output: results/proteins_to_download/{group}/{gene}.txt

Each output file contains the species names (one per line) that should
be searched for that specific gene.
"""

import pandas as pd
import os
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def safe_filename(name):
    """Create safe filename from gene name"""
    return name.replace("/", "_").replace("\\", "_").replace(" ", "_")

def create_gene_species_lists(proteins_file, species_file, output_dir, group):
    """
    Create individual species list files for each gene
    
    Args:
        proteins_file: TSV file with gene information
        species_file: Text file with species names (one per line)
        output_dir: Output directory for gene-specific species lists
        group: Group name (positive/negative) for logging
    """
    
    logging.info(f"Creating gene-specific species lists for gram {group}")
    logging.info(f"Proteins file: {proteins_file}")
    logging.info(f"Species file: {species_file}")
    logging.info(f"Output directory: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read proteins data to get list of genes
    try:
        proteins_df = pd.read_csv(proteins_file, sep='\t')
        genes = proteins_df['gene'].unique().tolist()
        logging.info(f"Found {len(genes)} unique genes in proteins file")
    except Exception as e:
        logging.error(f"Error reading proteins file: {e}")
        raise
    
    # Read species list
    try:
        with open(species_file, 'r') as f:
            species_list = [line.strip() for line in f if line.strip()]
        logging.info(f"Found {len(species_list)} species in species file")
    except Exception as e:
        logging.error(f"Error reading species file: {e}")
        raise
    
    # Create a file for each gene containing all species
    created_files = 0
    for gene in genes:
        # Create safe filename for gene
        safe_gene_name = safe_filename(gene)
        gene_file = os.path.join(output_dir, f"{safe_gene_name}.txt")
        
        # Write species list for this gene
        with open(gene_file, 'w') as f:
            for species in species_list:
                f.write(f"{species}\n")
        
        created_files += 1
        logging.info(f"Created species list for gene '{gene}': {gene_file} ({len(species_list)} species)")
    
    # Create summary file
    summary_file = os.path.join(output_dir, "_gene_summary.txt")
    with open(summary_file, 'w') as f:
        f.write(f"Gene-specific species lists for gram {group}\n")
        f.write(f"Generated on: {pd.Timestamp.now()}\n")
        f.write(f"Total genes: {len(genes)}\n")
        f.write(f"Total species per gene: {len(species_list)}\n")
        f.write(f"Files created: {created_files}\n\n")
        f.write("Genes:\n")
        for gene in sorted(genes):
            f.write(f"  {gene}.txt\n")
    
    logging.info(f"Summary written to: {summary_file}")
    logging.info(f"Successfully created {created_files} gene-specific species list files")
    
    return created_files

def main():
    """Main function for Snakemake integration"""
    logging.info("Gene-Specific Species Lists Creator")
    logging.info("=" * 50)
    
    try:
        # Get inputs from Snakemake
        proteins_file = snakemake.input.proteins
        species_file = snakemake.input.species
        output_dir = snakemake.output[0]
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        
        # Create gene-specific species lists
        created_files = create_gene_species_lists(proteins_file, species_file, output_dir, group)
        
        logging.info(f"Process completed successfully: {created_files} files created")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        
        # Example test parameters
        proteins_file = "results/proteins_to_study/analysis_1_params_1_gram_positive.tsv"
        species_file = "data/bacdive/analysis_1/gram_positive.txt"
        output_dir = "results/proteins_to_download/test_positive"
        group = "positive"
        
        if os.path.exists(proteins_file) and os.path.exists(species_file):
            created_files = create_gene_species_lists(proteins_file, species_file, output_dir, group)
            logging.info(f"Test completed: {created_files} files created")
        else:
            logging.error("Test files not found")
            logging.info(f"Expected proteins file: {proteins_file}")
            logging.info(f"Expected species file: {species_file}")

if __name__ == "__main__":
    main()