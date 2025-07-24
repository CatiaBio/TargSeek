#!/usr/bin/env python3
"""
Create Gene-Specific Species Lists from Coverage Data
====================================================

This script creates individual species lists for each gene based on the 
coverage analysis data. It uses the actual species that have each gene
(from the coverage TSV files) rather than the full BacDive species list.

This ensures we only try to download proteins for species that actually
have the gene according to the coverage analysis.

Input: 
- results/coverage/{analysis}_{paramset}_gram_{group}_coverage_count.tsv
- results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv

Output: results/proteins_to_download/{analysis}_{paramset}_gram_{group}/{gene}.txt
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

def create_gene_species_lists_from_coverage(coverage_file, proteins_file, output_dir, group):
    """
    Create individual species list files for each gene based on coverage data
    
    Args:
        coverage_file: TSV file with gene coverage data (includes species column)
        proteins_file: TSV file with surface-accessible proteins 
        output_dir: Output directory for gene-specific species lists
        group: Group name (positive/negative) for logging
    """
    
    logging.info(f"Creating gene-specific species lists from coverage data for gram {group}")
    logging.info(f"Coverage file: {coverage_file}")
    logging.info(f"Proteins file: {proteins_file}")
    logging.info(f"Output directory: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read unified coverage data to get gene-species mappings
    try:
        coverage_df = pd.read_csv(coverage_file, sep='\t')
        logging.info(f"Found {len(coverage_df)} gene-gram combinations in unified coverage file")
        
        # Filter by Gram type for this specific analysis
        if 'gram' in coverage_df.columns:
            # Convert group name to gram type (positive/negative)
            gram_type = group  # positive or negative from params
            filtered_df = coverage_df[coverage_df['gram'] == gram_type].copy()
            logging.info(f"Filtered to {len(filtered_df)} entries for Gram-{gram_type}")
            coverage_df = filtered_df
        else:
            logging.warning("No 'gram' column found in coverage file, using all data")
            
    except Exception as e:
        logging.error(f"Error reading coverage file: {e}")
        raise
    
    # Read proteins data to get list of target genes (only surface-accessible ones)
    try:
        # Check if it's a TSV or TXT file
        if proteins_file.endswith('.txt'):
            # Read as simple text file with gene names
            with open(proteins_file, 'r') as f:
                target_genes = set(line.strip() for line in f if line.strip())
            logging.info(f"Read {len(target_genes)} target genes from TXT file")
        else:
            # Read as TSV file
            proteins_df = pd.read_csv(proteins_file, sep='\t')
            # Handle both 'gene' and 'protein' column names
            if 'protein' in proteins_df.columns:
                target_genes = set(proteins_df['protein'].unique())
            elif 'gene' in proteins_df.columns:
                target_genes = set(proteins_df['gene'].unique())
            else:
                raise ValueError("Proteins file must have 'protein' or 'gene' column")
            logging.info(f"Found {len(target_genes)} target genes in TSV file")
    except Exception as e:
        logging.error(f"Error reading proteins file: {e}")
        raise
    
    # Filter coverage data to only target genes
    target_coverage = coverage_df[coverage_df['gene'].isin(target_genes)]
    logging.info(f"Filtered to {len(target_coverage)} target genes with coverage data")
    
    # Create a file for each gene containing only species that have that gene
    created_files = 0
    total_species_mapped = 0
    
    for _, row in target_coverage.iterrows():
        gene = row['gene']
        species_str = row['species_names_with_gene']  # Correct column name from unified coverage
        count = row['species_with_gene']  # Correct column name from unified coverage
        coverage_pct = row['coverage_percentage']
        
        # Parse comma-separated species list
        if pd.isna(species_str) or species_str.strip() == '':
            logging.warning(f"No species data for gene '{gene}', skipping")
            continue
            
        species_list = [s.strip() for s in species_str.split(',') if s.strip()]
        
        if not species_list:
            logging.warning(f"Empty species list for gene '{gene}', skipping")
            continue
        
        # Create safe filename for gene
        safe_gene_name = safe_filename(gene)
        gene_file = os.path.join(output_dir, f"{safe_gene_name}.txt")
        
        # Write species list for this gene
        with open(gene_file, 'w') as f:
            for species in species_list:
                f.write(f"{species}\n")
        
        created_files += 1
        total_species_mapped += len(species_list)
        logging.info(f"Created species list for gene '{gene}': {gene_file} ({len(species_list)} species, {coverage_pct:.1f}% coverage)")
    
    # Create summary file with detailed statistics
    summary_file = os.path.join(output_dir, "_gene_summary.txt")
    with open(summary_file, 'w') as f:
        f.write(f"Gene-specific species lists from coverage data for gram {group}\n")
        f.write(f"Generated on: {pd.Timestamp.now()}\n")
        f.write(f"Source coverage file: {coverage_file}\n")
        f.write(f"Source proteins file: {proteins_file}\n")
        f.write(f"Total target genes: {len(target_genes)}\n")
        f.write(f"Genes with coverage data: {len(target_coverage)}\n")
        f.write(f"Files created: {created_files}\n")
        f.write(f"Total species-gene mappings: {total_species_mapped}\n")
        f.write(f"Average species per gene: {total_species_mapped/created_files:.1f}\n\n")
        
        f.write("Gene details:\n")
        f.write("Gene\tSpecies_Count\tCoverage_%\tFile\n")
        for _, row in target_coverage.iterrows():
            gene = row['gene']
            count = row['species_with_gene']  # Use correct column name
            coverage_pct = row['coverage_percentage']
            safe_gene_name = safe_filename(gene)
            f.write(f"{gene}\t{count}\t{coverage_pct:.1f}\t{safe_gene_name}.txt\n")
    
    logging.info(f"Summary written to: {summary_file}")
    logging.info(f"Successfully created {created_files} gene-specific species list files")
    logging.info(f"Total species-gene mappings: {total_species_mapped}")
    logging.info(f"Average species per gene: {total_species_mapped/created_files:.1f}")
    
    return created_files

def main():
    """Main function for Snakemake integration"""
    logging.info("Gene-Specific Species Lists Creator (from Coverage Data)")
    logging.info("=" * 60)
    
    try:
        # Get inputs from Snakemake
        coverage_file = snakemake.input.coverage
        proteins_file = snakemake.input.genes
        output_dir = snakemake.output[0]
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        
        # Create gene-specific species lists from coverage data
        created_files = create_gene_species_lists_from_coverage(coverage_file, proteins_file, output_dir, group)
        
        logging.info(f"Process completed successfully: {created_files} files created")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        
        # Example test parameters
        coverage_file = "results/coverage/analysis_1_params_1_gram_positive_coverage_count.tsv"
        proteins_file = "results/analysis_1_params_1/proteins_to_study/gram_positive.tsv"
        output_dir = "results/proteins_to_download/test_coverage_positive"
        group = "positive"
        
        if os.path.exists(coverage_file) and os.path.exists(proteins_file):
            created_files = create_gene_species_lists_from_coverage(coverage_file, proteins_file, output_dir, group)
            logging.info(f"Test completed: {created_files} files created")
        else:
            logging.error("Test files not found")
            logging.info(f"Expected coverage file: {coverage_file}")
            logging.info(f"Expected proteins file: {proteins_file}")

if __name__ == "__main__":
    main()