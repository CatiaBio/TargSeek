#!/usr/bin/env python3
"""
Integrated Epitope Analysis Script

Combines data from:
1. BepiPred 3.0 epitope predictions
2. ConSurf conservation analysis  
3. DeepTMHMM topology predictions

Creates comprehensive epitope reports with conservation scores, topology annotations,
and visualizations.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
import sys
from pathlib import Path
import argparse
from typing import Dict, List, Tuple, Optional

def parse_topology_data(topology_file: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    Parse DeepTMHMM topology results from markdown file.
    
    Returns:
        Dictionary mapping sequence regions to their topology annotations
    """
    topology_regions = []
    
    with open(topology_file, 'r') as f:
        content = f.read()
    
    # Extract GFF3 section
    gff_section = re.search(r'```\n##gff-version 3\n(.*?)\n```', content, re.DOTALL)
    if not gff_section:
        print(f"Warning: Could not parse topology data from {topology_file}")
        return {}
    
    gff_lines = gff_section.group(1).strip().split('\n')
    
    for line in gff_lines:
        if line.startswith('#') or not line.strip():
            continue
            
        parts = line.split('\t')
        if len(parts) >= 4:
            seq_id = parts[0]
            region_type = parts[1]
            try:
                start = int(parts[2])
                end = int(parts[3])
                topology_regions.append((start, end, region_type))
            except ValueError:
                print(f"Warning: Could not parse topology line: {line.strip()}")
                continue
    
    return {'topology_regions': topology_regions}

def parse_consurf_data(consurf_file: str) -> pd.DataFrame:
    """
    Parse ConSurf conservation analysis results.
    
    Returns:
        DataFrame with position, conservation scores, and amino acid variety
    """
    # Read the file line by line to handle the complex header structure
    with open(consurf_file, 'r') as f:
        lines = f.readlines()
    
    # Find the header line with column names
    header_line = None
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith('pos,'):
            header_line = line.strip()
            data_start = i + 1
            break
    
    if header_line is None:
        raise ValueError(f"Could not find header line in {consurf_file}")
    
    # Parse header to get column names
    columns = [col.strip('"') for col in header_line.split(',')]
    
    # Read data lines
    data_rows = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        # Split by comma and clean up values
        values = [val.strip('"') for val in line.split(',')]
        if len(values) == len(columns):
            data_rows.append(values)
    
    # Create DataFrame
    df = pd.DataFrame(data_rows, columns=columns)
    
    # Extract useful columns
    conservation_df = df[['pos', 'ConSurf grade', 'MAX ACID']].copy()
    conservation_df.columns = ['position', 'conservation_score', 'max_residue']
    
    # Ensure proper data types
    conservation_df['position'] = conservation_df['position'].astype(int)
    conservation_df['conservation_score'] = conservation_df['conservation_score'].astype(float)
    
    return conservation_df

def parse_bepipred_data(bepipred_file: str) -> pd.DataFrame:
    """
    Parse BepiPred 3.0 epitope predictions.
    
    Returns:
        DataFrame with epitope predictions
    """
    epitopes = []
    
    with open(bepipred_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header and comment lines
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('No.'):
            continue
            
        parts = line.split('\t')
        if len(parts) >= 7:
            try:
                epitopes.append({
                    'epitope_no': int(parts[0]),
                    'chain': parts[1],
                    'start': int(parts[2]),
                    'end': int(parts[3]),
                    'peptide': parts[4],
                    'length': int(parts[5]),
                    'score': float(parts[6])
                })
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line: {line}")
                continue
    
    return pd.DataFrame(epitopes)

def annotate_epitope_topology(epitope_start: int, epitope_end: int, 
                             topology_regions: List[Tuple[int, int, str]]) -> List[str]:
    """
    Annotate epitope with topology information.
    
    Returns:
        List of topology regions that overlap with the epitope
    """
    overlapping_regions = []
    
    for region_start, region_end, region_type in topology_regions:
        # Check for overlap
        if not (epitope_end < region_start or epitope_start > region_end):
            overlapping_regions.append(region_type)
    
    return list(set(overlapping_regions))  # Remove duplicates

def calculate_epitope_conservation(epitope_start: int, epitope_end: int, 
                                 conservation_df: pd.DataFrame) -> Dict[str, float]:
    """
    Calculate conservation statistics for an epitope region.
    
    Returns:
        Dictionary with conservation statistics
    """
    epitope_conservation = conservation_df[
        (conservation_df['position'] >= epitope_start) & 
        (conservation_df['position'] <= epitope_end)
    ]
    
    if epitope_conservation.empty:
        return {
            'mean_conservation': 0.0,
            'min_conservation': 0.0,
            'max_conservation': 0.0,
            'conservation_std': 0.0
        }
    
    return {
        'mean_conservation': epitope_conservation['conservation_score'].mean(),
        'min_conservation': epitope_conservation['conservation_score'].min(),
        'max_conservation': epitope_conservation['conservation_score'].max(),
        'conservation_std': epitope_conservation['conservation_score'].std()
    }

def create_integrated_analysis(gene: str, structure_id: str, gram_type: str,
                             bepipred_file: str, consurf_file: str, 
                             topology_file: str, output_dir: str):
    """
    Create integrated epitope analysis combining all three tools.
    """
    print(f"Processing {gene} ({structure_id}) - {gram_type}")
    
    # Parse input data
    epitopes_df = parse_bepipred_data(bepipred_file)
    conservation_df = parse_consurf_data(consurf_file)
    topology_data = parse_topology_data(topology_file)
    
    if epitopes_df.empty:
        print(f"Warning: No epitopes found for {gene}")
        return
    
    # Create integrated analysis
    integrated_results = []
    
    for _, epitope in epitopes_df.iterrows():
        # Calculate conservation statistics
        conservation_stats = calculate_epitope_conservation(
            epitope['start'], epitope['end'], conservation_df
        )
        
        # Annotate topology
        topology_annotations = annotate_epitope_topology(
            epitope['start'], epitope['end'], 
            topology_data.get('topology_regions', [])
        )
        
        # Combine all data
        result = {
            'gene': gene,
            'structure_id': structure_id,
            'gram_type': gram_type,
            'epitope_no': epitope['epitope_no'],
            'chain': epitope['chain'],
            'start': epitope['start'],
            'end': epitope['end'],
            'peptide': epitope['peptide'],
            'length': epitope['length'],
            'bepipred_score': epitope['score'],
            'mean_conservation': conservation_stats['mean_conservation'],
            'min_conservation': conservation_stats['min_conservation'],
            'max_conservation': conservation_stats['max_conservation'],
            'conservation_std': conservation_stats['conservation_std'],
            'topology_regions': ';'.join(topology_annotations) if topology_annotations else 'Unknown',
            'accessibility': 'Surface' if any('outside' in region.lower() for region in topology_annotations) else 'Buried'
        }
        
        integrated_results.append(result)
    
    # Create output DataFrame
    results_df = pd.DataFrame(integrated_results)
    
    # Add composite scores
    results_df['conservation_rank'] = results_df['mean_conservation'].rank(ascending=False)
    results_df['bepipred_rank'] = results_df['bepipred_score'].rank(ascending=False)
    
    # Composite score: combine conservation and BepiPred scores
    # Normalize scores to 0-1 range
    norm_conservation = (results_df['mean_conservation'] - results_df['mean_conservation'].min()) / \
                       (results_df['mean_conservation'].max() - results_df['mean_conservation'].min())
    norm_bepipred = (results_df['bepipred_score'] - results_df['bepipred_score'].min()) / \
                   (results_df['bepipred_score'].max() - results_df['bepipred_score'].min())
    
    # Weight conservation more heavily (70% conservation, 30% BepiPred)
    results_df['composite_score'] = 0.7 * norm_conservation + 0.3 * norm_bepipred
    results_df['final_rank'] = results_df['composite_score'].rank(ascending=False)
    
    # Sort by composite score
    results_df = results_df.sort_values('composite_score', ascending=False)
    
    # Save integrated results
    os.makedirs(output_dir, exist_ok=True)
    
    # Save main results table
    output_file = os.path.join(output_dir, f"{gene}_{structure_id}_integrated_epitope_analysis.tsv")
    results_df.to_csv(output_file, sep='\t', index=False)
    
    # Create visualizations
    create_epitope_visualizations(results_df, gene, structure_id, output_dir)
    
    print(f"Integrated analysis saved to {output_file}")
    return results_df

def create_epitope_visualizations(results_df: pd.DataFrame, gene: str, 
                                structure_id: str, output_dir: str):
    """
    Create visualizations for epitope analysis.
    """
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'Epitope Analysis for {gene} ({structure_id})', fontsize=16, fontweight='bold')
    
    # 1. Epitope ranking by composite score
    ax1 = axes[0, 0]
    bars = ax1.bar(range(len(results_df)), results_df['composite_score'], 
                   color='steelblue', alpha=0.7)
    ax1.set_xlabel('Epitope Rank')
    ax1.set_ylabel('Composite Score')
    ax1.set_title('Epitope Ranking by Composite Score')
    ax1.set_xticks(range(len(results_df)))
    ax1.set_xticklabels([f"E{i}" for i in results_df['epitope_no']], rotation=45)
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9)
    
    # 2. Conservation vs BepiPred score scatter
    ax2 = axes[0, 1]
    scatter = ax2.scatter(results_df['bepipred_score'], results_df['mean_conservation'], 
                         c=results_df['composite_score'], cmap='viridis', 
                         s=100, alpha=0.7)
    ax2.set_xlabel('BepiPred Score')
    ax2.set_ylabel('Mean Conservation Score')
    ax2.set_title('Conservation vs BepiPred Score')
    
    # Add epitope labels
    for _, row in results_df.iterrows():
        ax2.annotate(f"E{row['epitope_no']}", 
                    (row['bepipred_score'], row['mean_conservation']),
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.colorbar(scatter, ax=ax2, label='Composite Score')
    
    # 3. Epitope positions and conservation
    ax3 = axes[1, 0]
    positions = [(row['start'] + row['end']) / 2 for _, row in results_df.iterrows()]
    colors = ['red' if 'outside' in row['topology_regions'].lower() else 'blue' 
              for _, row in results_df.iterrows()]
    
    scatter3 = ax3.scatter(positions, results_df['mean_conservation'], 
                          c=colors, s=results_df['length']*10, alpha=0.7)
    ax3.set_xlabel('Epitope Position (Center)')
    ax3.set_ylabel('Mean Conservation Score')
    ax3.set_title('Epitope Position vs Conservation\n(Red=Surface, Blue=Buried, Size=Length)')
    
    # Add epitope labels
    for _, row in results_df.iterrows():
        pos = (row['start'] + row['end']) / 2
        ax3.annotate(f"E{row['epitope_no']}", (pos, row['mean_conservation']),
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    # 4. Epitope properties summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    # Create summary table
    summary_data = []
    for _, row in results_df.iterrows():
        summary_data.append([
            f"E{row['epitope_no']}",
            f"{row['start']}-{row['end']}",
            f"{row['length']} aa",
            f"{row['bepipred_score']:.3f}",
            f"{row['mean_conservation']:.1f}",
            f"{row['composite_score']:.3f}",
            row['accessibility']
        ])
    
    headers = ['Epitope', 'Position', 'Length', 'BepiPred', 'Conservation', 'Composite', 'Location']
    
    table = ax4.table(cellText=summary_data, colLabels=headers,
                     cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    # Color code by composite score rank
    for i in range(len(summary_data)):
        for j in range(len(headers)):
            if i == 0:  # Best epitope
                table[(i+1, j)].set_facecolor('#90EE90')
            elif i == 1:  # Second best
                table[(i+1, j)].set_facecolor('#FFE4B5')
    
    ax4.set_title('Epitope Summary Table\n(Green=Best, Orange=Second Best)', pad=20)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = os.path.join(output_dir, f"{gene}_{structure_id}_epitope_analysis.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Visualization saved to {plot_file}")

def main():
    parser = argparse.ArgumentParser(description='Integrate epitope, conservation, and topology analysis')
    parser.add_argument('--gene', required=True, help='Gene name')
    parser.add_argument('--structure-id', required=True, help='Structure ID (e.g., 5OR1)')
    parser.add_argument('--gram-type', required=True, choices=['positive', 'negative'], 
                       help='Gram type (positive/negative)')
    parser.add_argument('--analysis', default='analysis1', help='Analysis name')
    parser.add_argument('--paramset', default='params1', help='Parameter set')
    parser.add_argument('--output-dir', help='Output directory (optional)')
    
    args = parser.parse_args()
    
    # Construct file paths
    base_dir = f"results/{args.analysis}_{args.paramset}/protein_analysis"
    
    bepipred_file = f"{base_dir}/bepipred_epitope_predictions/{args.gene}/{args.gene}_{args.structure_id}_linear_epitopes.tsv"
    consurf_file = f"{base_dir}/consurf_analysis/gram_{args.gram_type}/{args.gene}/msa_aa_variety_percentage.csv"
    topology_file = f"{base_dir}/topology_analysis/gram_{args.gram_type}/{args.gene}/biolib_results/deeptmhmm_results.md"
    
    # Set output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = f"{base_dir}/integrated_epitope_analysis/gram_{args.gram_type}/{args.gene}"
    
    # Check if input files exist
    for file_path, file_type in [(bepipred_file, "BepiPred"), 
                                (consurf_file, "ConSurf"), 
                                (topology_file, "Topology")]:
        if not os.path.exists(file_path):
            print(f"Error: {file_type} file not found: {file_path}")
            sys.exit(1)
    
    # Run integrated analysis
    try:
        results_df = create_integrated_analysis(
            args.gene, args.structure_id, args.gram_type,
            bepipred_file, consurf_file, topology_file, output_dir
        )
        
        print(f"\nIntegrated analysis completed successfully!")
        print(f"Best epitope: E{results_df.iloc[0]['epitope_no']} "
              f"({results_df.iloc[0]['peptide']}) - "
              f"Composite Score: {results_df.iloc[0]['composite_score']:.3f}")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()