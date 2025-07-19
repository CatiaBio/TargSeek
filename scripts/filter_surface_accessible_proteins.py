#!/usr/bin/env python3
"""
Filter Surface-Accessible Proteins
==================================

This script filters proteins based on their GO cellular component location.
Only proteins with GO cellular components matching the surface-accessible list
are retained for further study.
"""

import pandas as pd
from pathlib import Path
import sys


def load_surface_accessible_terms(file_path):
    """
    Load the list of surface-accessible GO cellular component terms
    
    Args:
        file_path: Path to surface_accessible.txt file
        
    Returns:
        set: Set of surface-accessible terms (lowercase for case-insensitive matching)
    """
    try:
        with open(file_path, 'r') as f:
            terms = [line.strip().lower() for line in f if line.strip()]
        print(f"Loaded {len(terms)} surface-accessible terms")
        return set(terms)
    except Exception as e:
        print(f"Error loading surface-accessible terms: {e}")
        sys.exit(1)


def is_surface_accessible(go_cellular_component, surface_terms):
    """
    Check if any GO cellular component matches surface-accessible terms
    
    Args:
        go_cellular_component: String containing GO cellular components (semicolon-separated)
        surface_terms: Set of surface-accessible terms
        
    Returns:
        bool: True if any component matches surface terms
    """
    if pd.isna(go_cellular_component) or not go_cellular_component:
        return False
    
    if go_cellular_component in ['Requires manual check', 'Not found in UniProt bacteria database',
                                  'Network error - requires manual check', 'Error - requires manual check']:
        return False
    
    # Split by semicolon and check each component
    components = [comp.strip().lower() for comp in go_cellular_component.split(';')]
    
    for component in components:
        if component in surface_terms:
            return True
    
    return False


def filter_proteins_by_location(input_file, output_file, surface_terms_file):
    """
    Filter proteins based on GO cellular component location
    
    Args:
        input_file: Path to enriched coverage file with GO cellular component
        output_file: Path to output filtered proteins file
        surface_terms_file: Path to surface_accessible.txt
    """
    
    print(f"\nProcessing {input_file}...")
    
    # Load surface-accessible terms
    surface_terms = load_surface_accessible_terms(surface_terms_file)
    
    # Read enriched coverage data
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Loaded {len(df)} genes from coverage file")
    except Exception as e:
        print(f"Error reading coverage file: {e}")
        return
    
    # Check if go_cellular_component column exists
    if 'go_cellular_component' not in df.columns:
        print("Error: 'go_cellular_component' column not found in input file")
        print("Available columns:", df.columns.tolist())
        return
    
    # Filter proteins with surface-accessible GO cellular components
    print("\nFiltering proteins by surface accessibility...")
    df['is_surface_accessible'] = df['go_cellular_component'].apply(
        lambda x: is_surface_accessible(x, surface_terms)
    )
    
    # Keep only surface-accessible proteins
    filtered_df = df[df['is_surface_accessible']].copy()
    filtered_df.drop('is_surface_accessible', axis=1, inplace=True)
    
    # Create output directory if needed
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save filtered results
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nFiltered results saved to: {output_file}")
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"Total genes analyzed: {len(df)}")
    print(f"Surface-accessible proteins: {len(filtered_df)} ({len(filtered_df)/len(df)*100:.1f}%)")
    
    # Show examples of filtered proteins
    if len(filtered_df) > 0:
        print(f"\nExamples of surface-accessible proteins:")
        for idx, row in filtered_df.head(5).iterrows():
            print(f"  - {row['gene']}: {row['go_cellular_component']}")
    
    # Show what was filtered out
    excluded_df = df[~df['is_surface_accessible']]
    if len(excluded_df) > 0:
        print(f"\nExamples of excluded proteins:")
        for idx, row in excluded_df.head(5).iterrows():
            go_comp = row['go_cellular_component']
            if pd.isna(go_comp) or not go_comp:
                go_comp = "No GO cellular component"
            print(f"  - {row['gene']}: {go_comp}")


def main():
    """Main function for Snakemake integration"""
    print("Surface-Accessible Protein Filter")
    print("="*50)
    
    try:
        # Get inputs from Snakemake
        input_file = snakemake.input.coverage
        surface_terms_file = snakemake.input.surface_terms
        output_file = snakemake.output[0]
        
        print(f"Input coverage file: {input_file}")
        print(f"Surface terms file: {surface_terms_file}")
        print(f"Output file: {output_file}")
        
        # Filter proteins
        filter_proteins_by_location(input_file, output_file, surface_terms_file)
        
    except NameError:
        print("Running in test mode (no Snakemake)")
        
        # Test with gram positive file
        test_files = [
            {
                'input': 'results/uniprot_info/analysis_1_params_1_gram_positive_uniprot_info/analysis_1_params_1_gram_positive_coverage_count_location.tsv',
                'output': 'results/proteins_to_study/gram_positive.tsv'
            },
            {
                'input': 'results/uniprot_info/analysis_1_params_1_gram_negative_uniprot_info/analysis_1_params_1_gram_negative_coverage_count_location.tsv',
                'output': 'results/proteins_to_study/gram_negative.tsv'
            }
        ]
        
        surface_terms_file = 'config/quickgo/surface_accessible.txt'
        
        for test_case in test_files:
            if Path(test_case['input']).exists():
                filter_proteins_by_location(
                    test_case['input'],
                    test_case['output'],
                    surface_terms_file
                )
                break
        else:
            print("No test files found")
            print("\nLooking for enriched coverage files:")
            for dir_path in Path("results/uniprot_info").glob("*"):
                if dir_path.is_dir():
                    for file in dir_path.glob("*_location.tsv"):
                        print(f"  - {file}")


if __name__ == "__main__":
    main()