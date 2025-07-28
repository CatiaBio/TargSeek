#!/usr/bin/env python3

"""
Script to create epitope reports for only the used structures
Creates TSV files in gene-specific BepiPred folders with the exact format requested
"""

import pandas as pd
import os
import sys
from pathlib import Path

def load_used_structures_mapping(mapping_file):
    """Load the used structures mapping file"""
    try:
        df = pd.read_csv(mapping_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading used structures mapping: {e}")
        sys.exit(1)

def load_structure_ranges(ranges_file):
    """Load structure ranges TSV file"""
    try:
        df = pd.read_csv(ranges_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading structure ranges: {e}")
        sys.exit(1)

def load_bepipred_output(csv_file):
    """Load BepiPred raw output CSV file"""
    try:
        df = pd.read_csv(csv_file)
        # Add position column (1-indexed)
        df['Position'] = range(1, len(df) + 1)
        return df
    except Exception as e:
        print(f"Error loading BepiPred output: {e}")
        return None

def filter_by_structure_range(df, start_pos, end_pos):
    """Filter dataframe by structure position range"""
    return df[(df['Position'] >= start_pos) & (df['Position'] <= end_pos)].copy()

def identify_linear_epitopes(df, min_length=5, min_score=0.1):
    """
    Identify linear epitopes: 5+ consecutive residues with BepiPred-3.0 linear epitope score ≥ 0.1
    Returns epitopes with additional metrics
    """
    epitopes = []
    current_epitope = []
    
    for idx, row in df.iterrows():
        score = row['BepiPred-3.0 linear epitope score']
        
        if score >= min_score:
            current_epitope.append({
                'position': row['Position'],
                'residue': row['Residue'],
                'linear_score': score,
                'raw_score': row['BepiPred-3.0 score']
            })
        else:
            # End of potential epitope
            if len(current_epitope) >= min_length:
                epitopes.append(current_epitope)
            current_epitope = []
    
    # Check last epitope
    if len(current_epitope) >= min_length:
        epitopes.append(current_epitope)
    
    return epitopes

def calculate_epitope_score(epitope):
    """
    Calculate overall epitope score (average of linear scores)
    """
    linear_scores = [aa['linear_score'] for aa in epitope]
    return sum(linear_scores) / len(linear_scores)

def create_epitope_tsv(epitopes, gene, structure_id, chain, output_file):
    """
    Create TSV file with the requested format
    """
    
    # Calculate scores and sort epitopes by score (descending)
    epitope_data = []
    for epitope in epitopes:
        start_pos = epitope[0]['position']
        end_pos = epitope[-1]['position']
        sequence = ''.join([aa['residue'] for aa in epitope])
        length = len(epitope)
        score = calculate_epitope_score(epitope)
        
        epitope_data.append({
            'start': start_pos,
            'end': end_pos,
            'sequence': sequence,
            'length': length,
            'score': score
        })
    
    # Sort by score (descending)
    epitope_data.sort(key=lambda x: x['score'], reverse=True)
    
    # Create TSV content
    with open(output_file, 'w') as f:
        f.write(f"Predicted Linear Epitope(s) {structure_id}:\n")
        f.write("No.\tChain\tStart\tEnd\tPeptide\tNumber of residues\tScore\n")
        
        for i, epitope in enumerate(epitope_data, 1):
            f.write(f"{i}\t{chain}\t{epitope['start']}\t{epitope['end']}\t{epitope['sequence']}\t{epitope['length']}\t{epitope['score']:.3f}\n")
    
    return len(epitope_data)

def process_used_structures(used_structures_file, structure_ranges_file, bepipred_dir):
    """Process only the used structures and create epitope reports"""
    
    # Load data files
    used_structures_df = load_used_structures_mapping(used_structures_file)
    ranges_df = load_structure_ranges(structure_ranges_file)
    
    print(f"Processing {len(used_structures_df)} used structures")
    
    processed_count = 0
    
    for _, row in used_structures_df.iterrows():
        gene = row['Gene']
        structure_id = row['Structure_ID']
        chain = row['Chain']
        
        print(f"Processing: {gene} - {structure_id} - Chain {chain}")
        
        # Find BepiPred raw output file
        bepipred_file = os.path.join(bepipred_dir, gene, f"{gene}_{structure_id}_raw_output.csv")
        
        if not os.path.exists(bepipred_file):
            print(f"  Warning: BepiPred file not found: {bepipred_file}")
            continue
        
        # Load BepiPred output
        bepipred_df = load_bepipred_output(bepipred_file)
        if bepipred_df is None:
            continue
        
        # Find structure range
        range_match = ranges_df[
            (ranges_df['Gene'] == gene) & 
            (ranges_df['Structure_ID'] == structure_id) & 
            (ranges_df['Chain'] == chain)
        ]
        
        if range_match.empty:
            print(f"  Warning: No structure range found for {gene}/{structure_id}/{chain}")
            continue
        
        start_pos = int(range_match.iloc[0]['Start'])
        end_pos = int(range_match.iloc[0]['End'])
        
        print(f"  Structure range: {start_pos}-{end_pos}")
        
        # Filter by structure range
        filtered_df = filter_by_structure_range(bepipred_df, start_pos, end_pos)
        
        if filtered_df.empty:
            print(f"  Warning: No data in range {start_pos}-{end_pos}")
            continue
        
        # Identify linear epitopes
        epitopes = identify_linear_epitopes(filtered_df)
        
        print(f"  Found {len(epitopes)} predicted epitopes")
        
        # Create output file in gene-specific folder
        gene_dir = os.path.join(bepipred_dir, gene)
        os.makedirs(gene_dir, exist_ok=True)
        
        output_file = os.path.join(gene_dir, f"{gene}_{structure_id}_epitopes.tsv")
        
        # Create TSV report
        num_epitopes = create_epitope_tsv(epitopes, gene, structure_id, chain, output_file)
        
        print(f"  Report saved: {output_file} ({num_epitopes} epitopes)")
        processed_count += 1
    
    print(f"\nProcessing complete! Generated {processed_count} epitope reports")

def main():
    # Define file paths
    used_structures_file = "results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/used_structures_mapping.tsv"
    structure_ranges_file = "data/protein_structures/structure_ranges.tsv"
    bepipred_dir = "results/analysis1_params1/protein_analysis/epitope_predictions_bepipred"
    
    # Validate inputs
    if not os.path.exists(used_structures_file):
        print(f"Error: Used structures mapping not found: {used_structures_file}")
        print("Run create_used_structures_mapping.py first")
        sys.exit(1)
        
    if not os.path.exists(structure_ranges_file):
        print(f"Error: Structure ranges file not found: {structure_ranges_file}")
        sys.exit(1)
        
    if not os.path.exists(bepipred_dir):
        print(f"Error: BepiPred directory not found: {bepipred_dir}")
        sys.exit(1)
    
    print(f"Used structures mapping: {used_structures_file}")
    print(f"Structure ranges: {structure_ranges_file}")
    print(f"BepiPred directory: {bepipred_dir}")
    print()
    
    process_used_structures(used_structures_file, structure_ranges_file, bepipred_dir)

if __name__ == "__main__":
    main()