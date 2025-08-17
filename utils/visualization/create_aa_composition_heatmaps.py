#!/usr/bin/env python3
"""
Create amino acid composition heatmaps for all bamA epitopes.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def parse_consurf_aa_data(consurf_file):
    """Parse ConSurf amino acid percentage data"""
    
    with open(consurf_file, 'r') as f:
        lines = f.readlines()
    
    # Find the header line
    header_line = None
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith('pos,'):
            header_line = line.strip()
            data_start = i + 1
            break
    
    # Parse header to get column names
    columns = [col.strip('"') for col in header_line.split(',')]
    
    # Read data lines
    data_rows = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        values = [val.strip('"') for val in line.split(',')]
        if len(values) == len(columns):
            data_rows.append(values)
    
    # Create DataFrame
    df = pd.DataFrame(data_rows, columns=columns)
    
    # Convert position to numeric
    df['pos'] = df['pos'].astype(int)
    
    # Convert amino acid percentages to numeric
    aa_columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    for col in aa_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    return df, aa_columns

def get_all_epitope_regions():
    """Get all bamA epitope regions from integrated analysis file"""
    integrated_file = "results/analysis1_params1/protein_analysis/integrated_epitope_analysis/gram_negative/bamA/bamA_5OR1_integrated_epitope_analysis.tsv"
    
    epitopes = []
    with open(integrated_file, 'r') as f:
        lines = f.readlines()
        
    for line in lines[1:]:  # Skip header
        if line.strip():
            parts = line.strip().split('\t')
            epitope_no = int(parts[3])
            start = int(parts[5])
            end = int(parts[6])
            peptide = parts[7]
            epitopes.append({
                "no": epitope_no, 
                "start": start, 
                "end": end, 
                "peptide": peptide
            })
    
    # Sort by epitope number
    epitopes.sort(key=lambda x: x['no'])
    return epitopes

def create_epitope_aa_composition_heatmap(epitope, df, aa_columns, output_dir):
    """Create amino acid composition heatmap for a single epitope"""
    
    # Get epitope data
    epitope_data = df[(df['pos'] >= epitope['start']) & (df['pos'] <= epitope['end'])].copy()
    
    if epitope_data.empty:
        print(f"No amino acid data found for epitope {epitope['no']}")
        return
    
    # Prepare amino acid percentage matrix
    aa_matrix = []
    positions = epitope_data['pos'].values
    
    for _, row in epitope_data.iterrows():
        aa_percentages = [float(row[aa]) for aa in aa_columns]
        aa_matrix.append(aa_percentages)
    
    aa_matrix = np.array(aa_matrix).T  # Transpose so amino acids are rows, positions are columns
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(10, len(positions) * 0.8), 8))
    
    # Create heatmap
    im = ax.imshow(aa_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=100)
    
    # Set x-axis (positions)
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(positions, fontsize=10)
    ax.set_xlabel('Amino Acid Position', fontsize=12, fontweight='bold')
    
    # Set y-axis (amino acids)
    ax.set_yticks(range(len(aa_columns)))
    ax.set_yticklabels(aa_columns, fontsize=10)
    ax.set_ylabel('Amino Acid', fontsize=12, fontweight='bold')
    
    # Add title
    ax.set_title(f'Epitope {epitope["no"]} Amino Acid Composition\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add percentage values on heatmap
    for i in range(len(aa_columns)):
        for j in range(len(positions)):
            value = aa_matrix[i, j]
            if value > 0:  # Only show non-zero values
                text_color = 'white' if value > 50 else 'black'
                ax.text(j, i, f'{value:.0f}%', ha='center', va='center', 
                       fontsize=8, color=text_color, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.05)
    cbar.set_label('Amino Acid Percentage (%)', fontsize=12, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    output_file = f"{output_dir}/bamA_epitope_{epitope['no']}_aa_composition.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    print(f"Epitope {epitope['no']} amino acid composition heatmap saved to: {output_file}")

def create_all_aa_composition_heatmaps():
    """Create amino acid composition heatmaps for all epitopes"""
    
    # Load ConSurf amino acid data
    consurf_file = "results/analysis1_params1/protein_analysis/consurf_analysis/gram_negative/bamA/msa_aa_variety_percentage.csv"
    df, aa_columns = parse_consurf_aa_data(consurf_file)
    epitopes = get_all_epitope_regions()
    
    # Create output directory
    output_dir = "results/analysis1_params1/protein_analysis"
    
    print(f"Creating amino acid composition heatmaps for {len(epitopes)} epitopes...")
    
    # Create heatmap for each epitope
    for epitope in epitopes:
        create_epitope_aa_composition_heatmap(epitope, df, aa_columns, output_dir)
    
    print("\n" + "="*60)
    print("AMINO ACID COMPOSITION HEATMAPS CREATED:")
    print("="*60)
    for epitope in epitopes:
        print(f"✅ Epitope {epitope['no']}: {epitope['peptide']} (positions {epitope['start']}-{epitope['end']})")
    print("="*60)
    print("\nFeatures of each heatmap:")
    print("• Y-axis: 20 standard amino acids (A-Y)")
    print("• X-axis: Amino acid positions within epitope")
    print("• Color intensity: Percentage of sequences with that amino acid")
    print("• Text overlay: Percentage values for non-zero entries")
    print("• Higher percentages = more conserved amino acid at that position")

if __name__ == "__main__":
    create_all_aa_composition_heatmaps()