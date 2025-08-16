#!/usr/bin/env python3
"""
Create plots showing amino acid presence across species from ConSurf r4s.res data.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re

def parse_r4s_file(r4s_file):
    """Parse ConSurf r4s.res file to extract amino acid and species count data"""
    
    data = []
    with open(r4s_file, 'r') as f:
        lines = f.readlines()
    
    # Find data lines (skip comments and headers)
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#') and not line.startswith('The alpha'):
            # Parse lines like: "    1     M  -1.337   [-1.361,-1.323] 0.02233  153/157"
            parts = line.split()
            if len(parts) >= 6:
                try:
                    pos = int(parts[0])
                    aa = parts[1]
                    score = float(parts[2])
                    # Extract species count from last part (e.g., "153/157")
                    species_info = parts[-1]
                    if '/' in species_info:
                        present, total = species_info.split('/')
                        present = int(present)
                        total = int(total)
                        percentage = (present / total) * 100
                        
                        data.append({
                            'position': pos,
                            'amino_acid': aa,
                            'conservation_score': score,
                            'species_present': present,
                            'total_species': total,
                            'percentage_present': percentage
                        })
                except (ValueError, IndexError):
                    continue
    
    return pd.DataFrame(data)

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

def create_separate_epitope_plots(epitope, df, output_dir):
    """Create separate plots for a single epitope in its own subfolder"""
    
    # Get epitope data
    epitope_data = df[(df['position'] >= epitope['start']) & (df['position'] <= epitope['end'])].copy()
    
    if epitope_data.empty:
        print(f"No species count data found for epitope {epitope['no']}")
        return
    
    # Create epitope-specific subfolder
    epitope_folder = f"{output_dir}/epitope{epitope['no']}_{epitope['peptide']}"
    os.makedirs(epitope_folder, exist_ok=True)
    
    positions = epitope_data['position'].values
    amino_acids = epitope_data['amino_acid'].values
    species_present = epitope_data['species_present'].values
    total_species = epitope_data['total_species'].values
    percentages = epitope_data['percentage_present'].values
    conservation_scores = epitope_data['conservation_score'].values
    
    # Plot 1: Species count (absolute numbers)
    fig1, ax1 = plt.subplots(figsize=(max(10, len(epitope_data) * 0.8), 6))
    bars1 = ax1.bar(range(len(positions)), species_present, 
                    color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axhline(y=total_species[0], color='red', linestyle='--', 
                label=f'Total species ({total_species[0]})')
    ax1.set_ylabel('Number of Species', fontsize=12, fontweight='bold')
    ax1.set_title(f'Epitope {epitope["no"]} - Species Count per Position\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                  fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add amino acid labels and values on bars
    for i, (bar, aa, count) in enumerate(zip(bars1, amino_acids, species_present)):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{aa}\n{count}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Set x-axis
    ax1.set_xticks(range(len(positions)))
    ax1.set_xticklabels(positions, fontsize=10)
    ax1.set_xlabel('Position', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    output1 = f"{epitope_folder}/species_count_absolute.png"
    plt.savefig(output1, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    # Plot 2: Percentage present
    fig2, ax2 = plt.subplots(figsize=(max(10, len(epitope_data) * 0.8), 6))
    bars2 = ax2.bar(range(len(positions)), percentages, 
                    color='forestgreen', alpha=0.7, edgecolor='black')
    ax2.set_ylabel('Percentage of Species (%)', fontsize=12, fontweight='bold')
    ax2.set_title(f'Epitope {epitope["no"]} - Percentage of Species with Each Amino Acid\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                  fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 100)
    
    # Add percentage values on bars
    for i, (bar, aa, pct) in enumerate(zip(bars2, amino_acids, percentages)):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{aa}\n{pct:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Set x-axis
    ax2.set_xticks(range(len(positions)))
    ax2.set_xticklabels(positions, fontsize=10)
    ax2.set_xlabel('Position', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    output2 = f"{epitope_folder}/species_count_percentage.png"
    plt.savefig(output2, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    # Plot 3: Conservation score vs species presence
    fig3, ax3 = plt.subplots(figsize=(10, 8))
    scatter = ax3.scatter(percentages, conservation_scores, 
                         c=range(len(positions)), cmap='viridis', 
                         s=100, alpha=0.7, edgecolors='black')
    ax3.set_xlabel('Percentage of Species with Amino Acid (%)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Normalized Rate4Site Conservation Score', fontsize=12, fontweight='bold')
    ax3.set_title(f'Epitope {epitope["no"]} - Normalized Rate4Site Conservation Score vs Species Presence\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                  fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Add amino acid labels
    for i, (pct, score, aa, pos) in enumerate(zip(percentages, conservation_scores, amino_acids, positions)):
        ax3.annotate(f'{aa}\n{pos}', (pct, score),
                    xytext=(5, 5), textcoords='offset points', 
                    fontsize=9, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Position Index', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    output3 = f"{epitope_folder}/conservation_vs_presence.png"
    plt.savefig(output3, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    print(f"Epitope {epitope['no']} plots saved to folder: {epitope_folder}/")
    print(f"  - {output1}")
    print(f"  - {output2}")
    print(f"  - {output3}")

def create_all_species_count_plots():
    """Create separate species count plots for all epitopes in individual subfolders"""
    
    # Load ConSurf r4s data
    r4s_file = "results/analysis1_params1/protein_analysis/consurf_analysis/gram_negative/bamA/r4s.res"
    df = parse_r4s_file(r4s_file)
    epitopes = get_all_epitope_regions()
    
    # Create output directory
    output_dir = "results/analysis1_params1/protein_analysis"
    
    print(f"Creating separate species count plots for {len(epitopes)} epitopes...")
    print(f"Total species in analysis: {df['total_species'].iloc[0] if not df.empty else 'Unknown'}")
    
    # Create separate plots for each epitope
    for epitope in epitopes:
        create_separate_epitope_plots(epitope, df, output_dir)
    
    print("\n" + "="*60)
    print("SEPARATE SPECIES COUNT PLOTS CREATED:")
    print("="*60)
    for epitope in epitopes:
        print(f"✅ Epitope {epitope['no']}: {epitope['peptide']} (positions {epitope['start']}-{epitope['end']})")
        print(f"   📁 Folder: epitope{epitope['no']}_{epitope['peptide']}/")
    print("="*60)
    print("\nEach epitope folder contains 3 separate plots:")
    print("• species_count_absolute.png: Absolute number of species with each amino acid")
    print("• species_count_percentage.png: Percentage of species with each amino acid")
    print("• conservation_vs_presence.png: Normalized Rate4Site conservation score vs species presence correlation")
    print("• Higher species counts = more conserved amino acid")
    print("• Lower Rate4Site scores = more conserved (negative values indicate high conservation)")

if __name__ == "__main__":
    create_all_species_count_plots()