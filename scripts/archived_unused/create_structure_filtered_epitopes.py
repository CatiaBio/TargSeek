#!/usr/bin/env python3

"""
Script to filter BepiPred epitope predictions based on 3D structure ranges
and identify predicted linear epitopes (5+ consecutive AA with score ≥0.1)
"""

import pandas as pd
import argparse
import os
from pathlib import Path
import sys

def load_structure_ranges(tsv_file):
    """Load structure ranges from TSV file"""
    try:
        ranges_df = pd.read_csv(tsv_file, sep='\t')
        return ranges_df
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
        sys.exit(1)

def filter_by_structure_range(df, start_pos, end_pos):
    """Filter dataframe by structure position range"""
    return df[(df['Position'] >= start_pos) & (df['Position'] <= end_pos)].copy()

def identify_linear_epitopes(df, min_length=5, min_score=0.1):
    """
    Identify linear epitopes: 5+ consecutive residues with BepiPred-3.0 linear epitope score ≥ 0.1
    """
    epitopes = []
    current_epitope = []
    
    for idx, row in df.iterrows():
        score = row['BepiPred-3.0 linear epitope score']
        
        if score >= min_score:
            current_epitope.append({
                'position': row['Position'],
                'residue': row['Residue'],
                'score': score,
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

def format_epitope_report(epitopes, gene, structure_id, chain, structure_range):
    """Format epitope results into a report"""
    report_lines = []
    report_lines.append(f"# BepiPred-3.0 Epitope Prediction Report")
    report_lines.append(f"# Gene: {gene}")
    report_lines.append(f"# Structure: {structure_id}")
    report_lines.append(f"# Chain: {chain}")
    report_lines.append(f"# Structure Range: {structure_range}")
    report_lines.append(f"# Generated: {pd.Timestamp.now()}")
    report_lines.append(f"#")
    report_lines.append(f"# Criteria: Linear epitope score ≥ 0.1, minimum 5 consecutive residues")
    report_lines.append(f"# Total predicted epitopes: {len(epitopes)}")
    report_lines.append("")
    
    if not epitopes:
        report_lines.append("No linear epitopes found meeting the criteria.")
        return "\n".join(report_lines)
    
    report_lines.append("Position\tSequence\tLength\tAvg_Linear_Score\tMax_Linear_Score\tAvg_Raw_Score\tMax_Raw_Score")
    
    for i, epitope in enumerate(epitopes, 1):
        start_pos = epitope[0]['position']
        end_pos = epitope[-1]['position']
        sequence = ''.join([aa['residue'] for aa in epitope])
        length = len(epitope)
        
        linear_scores = [aa['score'] for aa in epitope]
        raw_scores = [aa['raw_score'] for aa in epitope]
        
        avg_linear = sum(linear_scores) / len(linear_scores)
        max_linear = max(linear_scores)
        avg_raw = sum(raw_scores) / len(raw_scores)
        max_raw = max(raw_scores)
        
        position_range = f"{start_pos}-{end_pos}"
        
        report_lines.append(f"{position_range}\t{sequence}\t{length}\t{avg_linear:.4f}\t{max_linear:.4f}\t{avg_raw:.4f}\t{max_raw:.4f}")
        
        # Add detailed residue information
        report_lines.append(f"# Epitope {i} details:")
        for aa in epitope:
            report_lines.append(f"#   {aa['position']}: {aa['residue']} (Linear: {aa['score']:.4f}, Raw: {aa['raw_score']:.4f})")
        report_lines.append("")
    
    return "\n".join(report_lines)

def process_structure_epitopes(bepipred_dir, structure_ranges_file, output_dir):
    """Process all BepiPred outputs with structure filtering"""
    
    # Load structure ranges
    ranges_df = load_structure_ranges(structure_ranges_file)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all raw_output.csv files
    bepipred_path = Path(bepipred_dir)
    raw_output_files = list(bepipred_path.rglob("*_raw_output.csv"))
    
    if not raw_output_files:
        print(f"No raw_output.csv files found in {bepipred_dir}")
        return
    
    print(f"Found {len(raw_output_files)} BepiPred output files")
    
    processed_count = 0
    
    for csv_file in raw_output_files:
        try:
            # Extract gene and structure info from filename
            filename = csv_file.stem  # Remove .csv extension
            # Expected format: gene_structureID_raw_output
            parts = filename.replace('_raw_output', '').split('_')
            
            if len(parts) < 2:
                print(f"Warning: Cannot parse filename {filename}, skipping")
                continue
                
            gene = parts[0]
            structure_id = '_'.join(parts[1:])  # Handle structure IDs with underscores
            
            print(f"Processing: {gene} - {structure_id}")
            
            # Load BepiPred output
            bepipred_df = load_bepipred_output(csv_file)
            
            # Find matching structure ranges
            structure_matches = ranges_df[
                (ranges_df['Gene'] == gene) & 
                (ranges_df['Structure_ID'] == structure_id)
            ]
            
            if structure_matches.empty:
                print(f"  Warning: No structure range found for {gene}/{structure_id}")
                continue
            
            # Process each chain
            for _, row in structure_matches.iterrows():
                chain = row['Chain']
                start_pos = int(row['Start'])
                end_pos = int(row['End'])
                structure_range = row['Range']
                
                print(f"  Processing chain {chain}: {structure_range}")
                
                # Filter by structure range
                filtered_df = filter_by_structure_range(bepipred_df, start_pos, end_pos)
                
                if filtered_df.empty:
                    print(f"    Warning: No data in range {structure_range}")
                    continue
                
                # Identify linear epitopes
                epitopes = identify_linear_epitopes(filtered_df)
                
                print(f"    Found {len(epitopes)} predicted epitopes")
                
                # Generate report
                report = format_epitope_report(epitopes, gene, structure_id, chain, structure_range)
                
                # Save report
                output_filename = f"{gene}_{structure_id}_chain{chain}_epitopes.txt"
                output_path = os.path.join(output_dir, output_filename)
                
                with open(output_path, 'w') as f:
                    f.write(report)
                
                print(f"    Report saved: {output_filename}")
                processed_count += 1
                
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            continue
    
    print(f"\nProcessing complete! Generated {processed_count} epitope reports in {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Filter BepiPred epitopes by 3D structure ranges')
    parser.add_argument('--bepipred-dir', 
                       default='results/analysis1_params1/protein_analysis/epitope_predictions_bepipred',
                       help='Directory containing BepiPred raw_output.csv files')
    parser.add_argument('--structure-ranges', 
                       default='data/protein_structures/structure_ranges.tsv',
                       help='TSV file with structure ranges')
    parser.add_argument('--output-dir', 
                       default='results/analysis1_params1/protein_analysis/structure_filtered_epitopes',
                       help='Output directory for epitope reports')
    parser.add_argument('--min-score', type=float, default=0.1,
                       help='Minimum BepiPred-3.0 linear epitope score (default: 0.1)')
    parser.add_argument('--min-length', type=int, default=5,
                       help='Minimum epitope length (default: 5)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.bepipred_dir):
        print(f"Error: BepiPred directory not found: {args.bepipred_dir}")
        sys.exit(1)
        
    if not os.path.exists(args.structure_ranges):
        print(f"Error: Structure ranges file not found: {args.structure_ranges}")
        sys.exit(1)
    
    print(f"BepiPred directory: {args.bepipred_dir}")
    print(f"Structure ranges: {args.structure_ranges}")
    print(f"Output directory: {args.output_dir}")
    print(f"Minimum score: {args.min_score}")
    print(f"Minimum length: {args.min_length}")
    print()
    
    process_structure_epitopes(args.bepipred_dir, args.structure_ranges, args.output_dir)

if __name__ == "__main__":
    main()