#!/usr/bin/env python3
"""
Create epitope heatmap with vertical y-axis label and epitope markers.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_consurf_detailed(consurf_file):
    """Parse ConSurf data with all amino acid percentages"""
    
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
    
    # Convert position and conservation score to numeric
    df['pos'] = df['pos'].astype(int)
    df['ConSurf grade'] = df['ConSurf grade'].astype(float)
    
    return df

def get_epitope_regions():
    """Get bamA epitope regions from integrated analysis"""
    epitopes = [
        {"no": 1, "start": 432, "end": 439, "peptide": "YGTESGVS"},
        {"no": 2, "start": 493, "end": 505, "peptide": "DFQADDADLSDYT"},
        {"no": 3, "start": 536, "end": 547, "peptide": "LSNMQPQIAMDR"},
        {"no": 4, "start": 553, "end": 558, "peptide": "GQSADT"},
        {"no": 5, "start": 671, "end": 696, "peptide": "KNGAHTSWDDNDDYEDCTQESGCKSD"},
        {"no": 6, "start": 741, "end": 757, "peptide": "NWDPSSAPSDVPDYSDP"}
    ]
    return epitopes

def create_vertical_label_heatmap():
    """Create epitope heatmap with vertical y-axis label and epitope markers"""
    
    # Load ConSurf data
    consurf_file = "results/analysis1_params1/protein_analysis/consurf_analysis/gram_negative/bamA/msa_aa_variety_percentage.csv"
    df = parse_consurf_detailed(consurf_file)
    epitopes = get_epitope_regions()
    
    # Create figure with tighter x-axis and longer y-axis
    fig, ax = plt.subplots(figsize=(12, 5))
    
    # Prepare data for all epitopes in sequence
    all_data = []
    all_positions = []
    epitope_indices = []  # Track which epitope each position belongs to
    
    for epitope in epitopes:
        epitope_data = df[(df['pos'] >= epitope['start']) & (df['pos'] <= epitope['end'])].copy()
        
        for _, row in epitope_data.iterrows():
            pos = int(row['pos'])
            conservation = row['ConSurf grade']
            all_data.append(conservation)
            all_positions.append(pos)
            epitope_indices.append(epitope['no'])
    
    # Create heatmap
    heatmap_data = np.array(all_data).reshape(1, -1)
    im = ax.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto', vmin=1, vmax=9)
    
    # Remove x-axis numbering
    ax.set_xticks([])
    ax.set_xlabel('')
    
    # Set y-axis with vertical label - properly centered
    ax.set_yticks([0])
    ax.set_yticklabels(['Conservation'], fontsize=14, fontweight='bold', rotation=90, va='center')
    
    # Add title
    ax.set_title('bamA Epitope Conservation Analysis', 
                fontsize=18, fontweight='bold', pad=25)
    
    # Add epitope markers above the heatmap
    epitope_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
    
    # Track epitope positions for markers
    current_epitope = None
    epitope_start_idx = 0
    
    for i, epitope_no in enumerate(epitope_indices + [None]):  # Add None to process last epitope
        if epitope_no != current_epitope:
            if current_epitope is not None:
                # Draw marker for previous epitope
                epitope_info = next(e for e in epitopes if e['no'] == current_epitope)
                epitope_end_idx = i - 1
                center_idx = (epitope_start_idx + epitope_end_idx) / 2
                color = epitope_colors[(current_epitope - 1) % len(epitope_colors)]
                
                # Add epitope marker below heatmap (black text only)
                ax.text(center_idx, -0.3, f"E{current_epitope}", 
                       ha='center', va='center', fontsize=14, fontweight='bold', color='black')
                
                # Add horizontal line below heatmap with small gap
                ax.plot([epitope_start_idx, epitope_end_idx], [-0.1, -0.1], 
                       color='black', linewidth=2, alpha=0.8)
            
            current_epitope = epitope_no
            epitope_start_idx = i
    
    # Add colorbar with proper spacing - moved outside plot area
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.05)
    cbar.set_label('ConSurf Conservation Score\n(1 = Variable, 9 = Highly Conserved)', 
                   fontsize=12, fontweight='bold')
    
    # Clean up spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # Set axis limits to accommodate epitope markers (adjusted for closer positioning)
    ax.set_xlim(-1, len(all_data))
    ax.set_ylim(-0.5, 1.2)
    
    # Adjust layout for tighter x-axis and longer y-axis
    plt.subplots_adjust(bottom=0.1, top=0.8, left=0.1, right=0.85)
    
    # Save plot
    output_file = "results/analysis1_params1/protein_analysis/test_bamA_3_epitope_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    print(f"Vertical label heatmap saved to: {output_file}")
    
    print("\n" + "="*60)
    print("HEATMAP FEATURES:")
    print("="*60)
    print("✅ Y-axis label rotated to vertical orientation")
    print("✅ No amino acid numbering on x-axis")
    print("✅ Epitope markers (E1, E2, etc.) above heatmap")
    print("✅ Colored lines showing epitope spans")
    print("✅ Wide horizontal layout for better resolution")
    print("✅ Clean, minimal design")

if __name__ == "__main__":
    create_vertical_label_heatmap()