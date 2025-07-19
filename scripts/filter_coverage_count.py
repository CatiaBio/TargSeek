#!/usr/bin/env python3
"""
Coverage Count Summary Script
============================

This script reads a coverage file and creates a summary showing:
- gene: gene name
- count: number of species that have this gene (count > 0)
- species: comma-separated list of species that have the gene
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict

def main():
    """Main function"""
    print("Starting Coverage Count Summary...")
    
    # Get inputs from Snakemake
    try:
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        threshold = snakemake.params.threshold
        print(f"Input: {input_file}")
        print(f"Output: {output_file}")
        print(f"Threshold: {threshold}")
    except NameError:
        print("Snakemake object not available, using test values")
        input_file = "results/coverage/analysis_1_params_1_gram_positive_coverage.tsv"
        output_file = "results/coverage/analysis_1_params_1_gram_positive_coverage_count.tsv"
        threshold = 25  # Test with gram_positive threshold
    
    # Load coverage data
    print("Loading coverage data...")
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Loaded {len(df)} records from {input_file}")
        
        # Ensure the species_with_gene column is numeric
        df['species_with_gene'] = pd.to_numeric(df['species_with_gene'], errors='coerce')
        
        # Remove rows with NaN counts
        df = df.dropna(subset=['species_with_gene'])
        print(f"After cleaning: {len(df)} records")
        
    except Exception as e:
        print(f"Error loading coverage data: {e}")
        return
    
    # Group by gene and summarize
    print("Grouping by gene and counting species...")
    gene_summary = {}
    
    # Get total number of species for coverage calculation
    total_species = df['total_species'].iloc[0] if 'total_species' in df.columns else len(df['gene'].unique())
    print(f"Total species count: {total_species}")
    
    for gene in df['gene'].unique():
        gene_data = df[df['gene'] == gene]
        
        # Filter species with species_with_gene > 0 (species that have the gene)
        species_with_gene = gene_data[gene_data['species_with_gene'] > 0]
        
        if len(species_with_gene) > 0:
            species_count = species_with_gene['species_with_gene'].iloc[0]
            coverage_percentage = species_with_gene['coverage_percentage'].iloc[0]
            species_string = species_with_gene['species_names_with_gene'].iloc[0]
            
            # Filter by coverage percentage (0.25 = 25%)
            if coverage_percentage >= 25.0:
                gene_summary[gene] = {
                    'count': species_count,
                    'coverage_percentage': coverage_percentage,
                    'species': species_string
                }
    
    print(f"Found {len(gene_summary)} genes with species above threshold")
    
    # Create output directory
    output_path = Path(output_file)
    
    # Try multiple methods to create directory and save file
    saved_successfully = False
    
    # Method 1: pathlib with mkdir
    if not saved_successfully:
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open('w') as f:
                f.write("gene\tcount\tcoverage_percentage\tspecies\n")
                
                # Sort genes by coverage percentage (descending) then by gene name
                sorted_genes = sorted(gene_summary.items(), 
                                    key=lambda x: (-x[1]['coverage_percentage'], x[0]))
                
                for gene, data in sorted_genes:
                    f.write(f"{gene}\t{data['count']}\t{data['coverage_percentage']:.2f}\t{data['species']}\n")
            
            print(f"Results saved to {output_file}")
            saved_successfully = True
        except Exception as e:
            print(f"Method 1 failed: {e}")
    
    # Method 2: os.makedirs with string paths
    if not saved_successfully:
        try:
            import os
            os.makedirs(str(output_path.parent), exist_ok=True)
            with open(str(output_path), 'w') as f:
                f.write("gene\tcount\tcoverage_percentage\tspecies\n")
                
                # Sort genes by coverage percentage (descending) then by gene name
                sorted_genes = sorted(gene_summary.items(), 
                                    key=lambda x: (-x[1]['coverage_percentage'], x[0]))
                
                for gene, data in sorted_genes:
                    f.write(f"{gene}\t{data['count']}\t{data['coverage_percentage']:.2f}\t{data['species']}\n")
            
            print(f"Results saved with os.makedirs: {output_file}")
            saved_successfully = True
        except Exception as e:
            print(f"Method 2 failed: {e}")
    
    # Method 3: Create directory step by step
    if not saved_successfully:
        try:
            # Create results directory first
            results_dir = Path("results")
            results_dir.mkdir(exist_ok=True)
            print(f"Created results directory: {results_dir.absolute()}")
            
            # Create coverage directory
            coverage_dir = results_dir / "coverage"
            coverage_dir.mkdir(exist_ok=True)
            print(f"Created coverage directory: {coverage_dir.absolute()}")
            
            # Now save file
            final_path = coverage_dir / output_path.name
            with final_path.open('w') as f:
                f.write("gene\tcount\tcoverage_percentage\tspecies\n")
                
                # Sort genes by coverage percentage (descending) then by gene name
                sorted_genes = sorted(gene_summary.items(), 
                                    key=lambda x: (-x[1]['coverage_percentage'], x[0]))
                
                for gene, data in sorted_genes:
                    f.write(f"{gene}\t{data['count']}\t{data['coverage_percentage']:.2f}\t{data['species']}\n")
            
            print(f"Results saved with step-by-step creation: {final_path.absolute()}")
            saved_successfully = True
        except Exception as e:
            print(f"Method 3 failed: {e}")
    
    # Method 4: Save to current directory with warning
    if not saved_successfully:
        try:
            fallback_file = Path(output_path.name)
            with fallback_file.open('w') as f:
                f.write("gene\tcount\tcoverage_percentage\tspecies\n")
                
                # Sort genes by coverage percentage (descending) then by gene name
                sorted_genes = sorted(gene_summary.items(), 
                                    key=lambda x: (-x[1]['coverage_percentage'], x[0]))
                
                for gene, data in sorted_genes:
                    f.write(f"{gene}\t{data['count']}\t{data['coverage_percentage']:.2f}\t{data['species']}\n")
            
            print(f"Results saved to current directory: {fallback_file.absolute()}")
            print(f"WARNING: File saved to current directory instead of {output_file}")
            saved_successfully = True
        except Exception as e:
            print(f"Method 4 failed: {e}")
    
    if not saved_successfully:
        print("ERROR: Could not save results file with any method!")
    else:
        print(f"\nSummary complete!")
        print(f"Total genes processed: {len(gene_summary)}")
        if gene_summary:
            max_count = max(data['count'] for data in gene_summary.values())
            min_count = min(data['count'] for data in gene_summary.values())
            print(f"Species count range: {min_count} - {max_count}")

if __name__ == "__main__":
    main()