#!/usr/bin/env python3
"""
Process Custom Species List for PureMilk Pipeline

This script processes a user-provided species list and creates the necessary
files for the custom pipeline, including Gram classification.

Input:
- Custom species list file (one species per line)

Output:
- all_identified.txt: All species from the input
- updated_gram.tsv: Gram classification table
- gram_positive.txt: List of Gram-positive species
- gram_negative.txt: List of Gram-negative species

Author: Generated for PureMilk Custom Pipeline
"""

import os
import sys
import pandas as pd
import random
from pathlib import Path

def classify_gram_by_genus(species_name, default_classification="positive"):
    """
    Simple Gram classification based on common genus patterns.
    
    Args:
        species_name (str): Scientific name of the species
        default_classification (str): Default classification if genus not recognized
    
    Returns:
        str: "positive" or "negative"
    """
    
    # Extract genus (first word)
    genus = species_name.split()[0].lower()
    
    # Known Gram-positive genera (common ones)
    gram_positive_genera = {
        'lactobacillus', 'streptococcus', 'staphylococcus', 'bacillus',
        'clostridium', 'enterococcus', 'bifidobacterium', 'propionibacterium',
        'corynebacterium', 'mycobacterium', 'actinomyces', 'listeria',
        'leuconostoc', 'pediococcus', 'lactococcus', 'weissella'
    }
    
    # Known Gram-negative genera (common ones)
    gram_negative_genera = {
        'escherichia', 'salmonella', 'shigella', 'klebsiella', 'enterobacter',
        'citrobacter', 'serratia', 'proteus', 'yersinia', 'vibrio',
        'pseudomonas', 'acinetobacter', 'haemophilus', 'neisseria',
        'helicobacter', 'campylobacter', 'bacteroides', 'prevotella',
        'fusobacterium', 'porphyromonas', 'veillonella', 'alcaligenes',
        'bordetella', 'burkholderia', 'stenotrophomonas', 'aeromonas'
    }
    
    if genus in gram_positive_genera:
        return "positive"
    elif genus in gram_negative_genera:
        return "negative"
    else:
        # Use default classification for unknown genera
        if default_classification == "mixed":
            return random.choice(["positive", "negative"])
        else:
            return default_classification

def main():
    # Get parameters from snakemake
    species_file = snakemake.input.species_list
    analysis = snakemake.params.analysis
    default_gram = snakemake.params.default_gram_classification
    
    # Output files
    all_species_file = snakemake.output.all_species
    gram_classification_file = snakemake.output.gram_classification
    gram_positive_file = snakemake.output.gram_positive
    gram_negative_file = snakemake.output.gram_negative
    
    print(f"Processing custom species list for analysis: {analysis}")
    print(f"Input file: {species_file}")
    print(f"Default Gram classification: {default_gram}")
    
    # Ensure output directories exist
    for output_file in [all_species_file, gram_classification_file, 
                       gram_positive_file, gram_negative_file]:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
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
    
    # Process each species
    gram_data = []
    gram_positive = []
    gram_negative = []
    
    for species in species_list:
        # Clean species name
        species = species.strip()
        if not species:
            continue
            
        # Classify by Gram stain
        gram_classification = classify_gram_by_genus(species, default_gram)
        
        # Add to data
        gram_data.append({
            'species': species,
            'gram_stain': gram_classification
        })
        
        # Add to appropriate lists
        if gram_classification == "positive":
            gram_positive.append(species)
        else:
            gram_negative.append(species)
    
    # Create DataFrame for gram classification
    gram_df = pd.DataFrame(gram_data)
    
    # Write all species file
    with open(all_species_file, 'w', encoding='utf-8') as f:
        for species in species_list:
            f.write(f"{species}\n")
    
    # Write gram classification TSV
    gram_df.to_csv(gram_classification_file, sep='\t', index=False)
    
    # Write gram-positive species
    with open(gram_positive_file, 'w', encoding='utf-8') as f:
        for species in gram_positive:
            f.write(f"{species}\n")
    
    # Write gram-negative species  
    with open(gram_negative_file, 'w', encoding='utf-8') as f:
        for species in gram_negative:
            f.write(f"{species}\n")
    
    # Print summary
    print(f"\nProcessing complete!")
    print(f"Total species: {len(species_list)}")
    print(f"Gram-positive: {len(gram_positive)}")
    print(f"Gram-negative: {len(gram_negative)}")
    print(f"Files created:")
    print(f"  - {all_species_file}")
    print(f"  - {gram_classification_file}")
    print(f"  - {gram_positive_file}")
    print(f"  - {gram_negative_file}")

if __name__ == "__main__":
    main()