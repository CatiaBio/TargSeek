#!/usr/bin/env python3
"""
Create proteins_to_study file from coverage results.

This script takes the coverage TSV file (pivot format) and converts it to the 
proteins_to_study format by:
1. Calculating coverage statistics per gene
2. Adding placeholder GO cellular component information
3. Filtering by minimum coverage threshold
4. Creating the final TSV with the required format
"""

import pandas as pd
import sys
from pathlib import Path

# Default GO cellular component mappings for common genes
GO_CELLULAR_COMPONENTS = {
    "eno": "cell surface; extracellular region; phosphopyruvate hydratase complex",
    "bamA": "cell outer membrane; pore complex",
    "htrA": "periplasmic space",
    "dnaK": "cytoplasm",
    "era": "cytoplasm; ribosome",
    "cls": "membrane",
    "pal": "cell outer membrane",
    "lamB": "cell outer membrane; pore complex",
    "pstS": "periplasmic space",
    "ushA": "extracellular region; membrane; outer membrane-bounded periplasmic space",
    "pepN": "cytoplasm; extracellular space; membrane",
    "cheZ": "bacterial-type flagellum; cytoplasm",
    "ompA": "cell outer membrane; pore complex",
    "traA": "extracellular region; plasma membrane",
    "ail": "cell outer membrane",
    # Add more as needed
}

def process_coverage_to_proteins_to_study(coverage_file, output_file, min_coverage=50.0):
    """
    Convert coverage TSV to proteins_to_study format
    
    Args:
        coverage_file: Path to coverage TSV file (long format: gene, species_with_gene, total_species, coverage_percentage, species_names_with_gene)
        output_file: Path to output proteins_to_study TSV
        min_coverage: Minimum coverage percentage to include gene
    """
    print(f"Reading coverage file: {coverage_file}")
    
    # Read the coverage TSV (long format)
    df = pd.read_csv(coverage_file, sep='\t')
    
    print(f"Loaded coverage data: {len(df)} genes")
    
    results = []
    
    for _, row in df.iterrows():
        gene = row['gene']
        coverage_percentage = row['coverage_percentage']
        count = row['species_with_gene']
        species_list = row['species_names_with_gene']
        
        print(f"  {gene}: {count} species ({coverage_percentage:.1f}%)")
        
        # Only include genes above minimum coverage
        if coverage_percentage >= min_coverage:
            # Get GO cellular component (use default or placeholder)
            go_component = GO_CELLULAR_COMPONENTS.get(gene, "unknown cellular component")
            
            # Add to results
            results.append({
                "gene": gene,
                "go_cellular_component": go_component,
                "count": count,
                "coverage_percentage": round(coverage_percentage, 2),
                "species": species_list
            })
    
    # Create DataFrame and sort by coverage percentage (descending)
    if results:
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('coverage_percentage', ascending=False)
        
        # Save to TSV
        results_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"\nSaved {len(results_df)} genes to {output_file}")
        print(f"Coverage range: {results_df['coverage_percentage'].min():.1f}% - {results_df['coverage_percentage'].max():.1f}%")
    else:
        # Create empty DataFrame with correct columns
        results_df = pd.DataFrame(columns=["gene", "go_cellular_component", "count", "coverage_percentage", "species"])
        results_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"\nNo genes met the minimum coverage threshold of {min_coverage}%")
        print(f"Saved empty results to {output_file}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python create_proteins_to_study_from_coverage.py <coverage_tsv> <output_tsv> [min_coverage]")
        print("  coverage_tsv: Input coverage TSV file (pivot format)")
        print("  output_tsv: Output proteins_to_study TSV file")
        print("  min_coverage: Minimum coverage percentage (default: 50.0)")
        sys.exit(1)
    
    coverage_file = sys.argv[1]
    output_file = sys.argv[2]
    min_coverage = float(sys.argv[3]) if len(sys.argv) > 3 else 50.0
    
    # Verify input file exists
    if not Path(coverage_file).exists():
        print(f"Error: Coverage file not found: {coverage_file}")
        sys.exit(1)
    
    # Create output directory if needed
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Process the data
    process_coverage_to_proteins_to_study(coverage_file, output_file, min_coverage)

if __name__ == "__main__":
    # Support both command line and Snakemake usage
    if 'snakemake' in globals():
        # Running from Snakemake
        coverage_file = snakemake.input.coverage
        output_file = snakemake.output[0]
        min_coverage = getattr(snakemake.params, 'min_coverage', 50.0)
        
        process_coverage_to_proteins_to_study(coverage_file, output_file, min_coverage)
    else:
        # Running from command line
        main()