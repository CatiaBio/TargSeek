#!/usr/bin/env python3
"""
Process Custom Gene List for PureMilk Pipeline

This script processes a user-provided gene list and creates the necessary
files for the custom pipeline.

Input:
- Custom gene list file (one gene symbol per line)

Output:
- gene_symbols.txt: Cleaned gene symbols
- gene_aliases.txt: Gene aliases (same as symbols for custom genes)
- surface_accessible_proteins.txt: All genes (assumes surface accessibility)
- proteins_to_be_tested.txt: Final gene list for testing

Author: Generated for PureMilk Custom Pipeline
"""

import os
import sys
from pathlib import Path

def main():
    # Get parameters from snakemake
    gene_file = snakemake.input.gene_list
    paramset = snakemake.params.paramset
    
    # Output files
    genes_file = snakemake.output.genes
    aliases_file = snakemake.output.aliases
    surface_genes_file = snakemake.output.surface_genes
    proteins_to_test_file = snakemake.output.proteins_to_test
    
    print(f"Processing custom gene list for paramset: {paramset}")
    print(f"Input file: {gene_file}")
    
    # Ensure output directories exist
    for output_file in [genes_file, aliases_file, surface_genes_file, proteins_to_test_file]:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Read gene list
    try:
        with open(gene_file, 'r', encoding='utf-8') as f:
            gene_list = [line.strip() for line in f if line.strip()]
        
        print(f"Found {len(gene_list)} genes in input file")
        
        if not gene_list:
            raise ValueError("No genes found in input file")
            
    except Exception as e:
        print(f"Error reading gene file: {e}")
        sys.exit(1)
    
    # Clean and deduplicate gene names
    cleaned_genes = []
    seen_genes = set()
    
    for gene in gene_list:
        gene = gene.strip()
        if not gene or gene in seen_genes:
            continue
        
        # Basic cleaning - remove common prefixes/suffixes that might cause issues
        cleaned_gene = gene.replace(' ', '_').replace('-', '_')
        
        cleaned_genes.append(cleaned_gene)
        seen_genes.add(gene)
    
    print(f"Cleaned gene list: {len(cleaned_genes)} unique genes")
    
    # Write gene symbols file
    with open(genes_file, 'w', encoding='utf-8') as f:
        for gene in cleaned_genes:
            f.write(f"{gene}\n")
    
    # Write gene aliases file (for custom genes, aliases are the same as symbols)
    # Format: gene_symbol\taliases (tab-separated)
    with open(aliases_file, 'w', encoding='utf-8') as f:
        f.write("gene_symbol\taliases\n")  # Header
        for gene in cleaned_genes:
            # For custom genes, use the gene name as its own alias
            f.write(f"{gene}\t{gene}\n")
    
    # Write surface accessible proteins (assume all custom genes are surface accessible)
    with open(surface_genes_file, 'w', encoding='utf-8') as f:
        for gene in cleaned_genes:
            f.write(f"{gene}\n")
    
    # Write proteins to be tested (final list)
    with open(proteins_to_test_file, 'w', encoding='utf-8') as f:
        for gene in cleaned_genes:
            f.write(f"{gene}\n")
    
    # Print summary
    print(f"\nProcessing complete!")
    print(f"Total genes processed: {len(cleaned_genes)}")
    print(f"Files created:")
    print(f"  - {genes_file}")
    print(f"  - {aliases_file}")
    print(f"  - {surface_genes_file}")
    print(f"  - {proteins_to_test_file}")

if __name__ == "__main__":
    main()