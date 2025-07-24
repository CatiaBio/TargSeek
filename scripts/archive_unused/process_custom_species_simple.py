#!/usr/bin/env python3
"""
Simple Species List Processing for Custom Pipeline

This script processes a user-provided species list without Gram classification.

Input:
- Custom species list file (one species per line)

Output:
- all_species.txt: All species from the input

Author: Generated for PureMilk Custom Pipeline
"""

import os
import sys
from pathlib import Path

def main():
    # Get parameters from snakemake
    species_file = snakemake.input.species_list
    analysis = snakemake.params.analysis
    
    # Output files
    all_species_file = snakemake.output.all_species
    
    print(f"Processing custom species list for analysis: {analysis}")
    print(f"Input file: {species_file}")
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(all_species_file), exist_ok=True)
    
    # Read species list
    try:
        with open(species_file, 'r', encoding='utf-8') as f:
            species_list = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
        
        print(f"Found {len(species_list)} species in input file")
        
        if not species_list:
            raise ValueError("No species found in input file")
            
    except Exception as e:
        print(f"Error reading species file: {e}")
        sys.exit(1)
    
    # Write all species file
    with open(all_species_file, 'w', encoding='utf-8') as f:
        for species in species_list:
            f.write(f"{species.strip()}\n")
    
    # Print summary
    print(f"\nProcessing complete!")
    print(f"Total species: {len(species_list)}")
    print(f"File created: {all_species_file}")

if __name__ == "__main__":
    main()