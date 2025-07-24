#!/usr/bin/env python3
"""
Amino Acid Conservation Analysis
===============================

This script analyzes amino acid conservation in trimmed alignments to identify:
- Highly conserved positions (potential functional sites)
- Conservation scores per position
- Consensus sequences
- Logomaker-style conservation plots
"""

from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import logging
from collections import Counter
import json

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_shannon_entropy(column):
    """Calculate Shannon entropy for an alignment column"""
    if len(column) == 0:
        return 0
    
    # Count amino acids (excluding gaps)
    aa_counts = Counter([aa for aa in column if aa != '-'])
    
    if len(aa_counts) == 0:
        return 0  # All gaps
    
    total = sum(aa_counts.values())
    entropy = 0
    
    for count in aa_counts.values():
        if count > 0:
            freq = count / total
            entropy -= freq * np.log2(freq)
    
    return entropy

def calculate_conservation_score(column):
    """Calculate conservation score (0 = low conservation, 1 = high conservation)"""
    max_entropy = np.log2(20)  # Maximum possible entropy for 20 amino acids
    entropy = calculate_shannon_entropy(column)
    
    # Convert to conservation score (inverse of normalized entropy)
    conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
    return conservation

def get_consensus_aa(column, min_frequency=0.5):
    """Get consensus amino acid for a position"""
    aa_counts = Counter([aa for aa in column if aa != '-'])
    
    if not aa_counts:
        return '-'
    
    most_common_aa, count = aa_counts.most_common(1)[0]
    frequency = count / len([aa for aa in column if aa != '-'])
    
    # Return consensus only if it meets minimum frequency threshold
    return most_common_aa if frequency >= min_frequency else 'X'

def analyze_alignment_conservation(alignment_file):
    """Analyze conservation for a single alignment"""
    try:
        # Read alignment
        alignment = AlignIO.read(alignment_file, "fasta")
        
        if len(alignment) < 2:
            logging.warning(f"Alignment {alignment_file} has fewer than 2 sequences")
            return None
        
        alignment_length = alignment.get_alignment_length()
        num_sequences = len(alignment)
        
        # Calculate conservation metrics for each position
        conservation_data = []
        consensus_sequence = []
        
        for i in range(alignment_length):
            column = alignment[:, i]  # Get column (all AAs at position i)
            
            # Calculate metrics
            conservation_score = calculate_conservation_score(column)
            entropy = calculate_shannon_entropy(column)
            consensus_aa = get_consensus_aa(column)
            gap_frequency = column.count('-') / len(column)
            
            # Count amino acids
            aa_counts = Counter([aa for aa in column if aa != '-'])
            most_common = aa_counts.most_common(1)[0] if aa_counts else ('-', 0)
            
            conservation_data.append({
                'position': i + 1,  # 1-indexed
                'conservation_score': conservation_score,
                'entropy': entropy,
                'consensus_aa': consensus_aa,
                'gap_frequency': gap_frequency,
                'most_common_aa': most_common[0],
                'most_common_count': most_common[1],
                'most_common_freq': most_common[1] / num_sequences if num_sequences > 0 else 0
            })
            
            consensus_sequence.append(consensus_aa)
        
        # Create summary statistics
        conservation_scores = [pos['conservation_score'] for pos in conservation_data]
        
        summary = {
            'gene': alignment_file.stem.replace('_trimmed', '').replace('_aligned', ''),
            'num_sequences': num_sequences,
            'alignment_length': alignment_length,
            'mean_conservation': np.mean(conservation_scores),
            'median_conservation': np.median(conservation_scores),
            'highly_conserved_positions': len([s for s in conservation_scores if s > 0.8]),
            'moderately_conserved_positions': len([s for s in conservation_scores if 0.5 < s <= 0.8]),
            'variable_positions': len([s for s in conservation_scores if s <= 0.5]),
            'consensus_sequence': ''.join(consensus_sequence)
        }
        
        return {
            'summary': summary,
            'position_data': conservation_data
        }
        
    except Exception as e:
        logging.error(f"Error analyzing {alignment_file}: {e}")
        return None

def create_conservation_plot(conservation_data, output_file, gene_name):
    """Create conservation plot for a gene"""
    try:
        positions = [pos['position'] for pos in conservation_data]
        scores = [pos['conservation_score'] for pos in conservation_data]
        
        plt.figure(figsize=(12, 6))
        
        # Create bar plot colored by conservation level
        colors = ['red' if s > 0.8 else 'orange' if s > 0.5 else 'lightblue' for s in scores]
        plt.bar(positions, scores, color=colors, alpha=0.7)
        
        # Add horizontal lines for thresholds
        plt.axhline(y=0.8, color='red', linestyle='--', alpha=0.5, label='High conservation (>0.8)')
        plt.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5, label='Moderate conservation (>0.5)')
        
        plt.xlabel('Alignment Position')
        plt.ylabel('Conservation Score')
        plt.title(f'Amino Acid Conservation - {gene_name}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"  Created conservation plot: {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Error creating plot for {gene_name}: {e}")
        return False

def main():
    """Main function for Snakemake integration"""
    
    # Get parameters from Snakemake
    input_dir = Path(snakemake.input[0])
    output_dir = Path(snakemake.output[0])
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    group = snakemake.params.group
    
    logging.info(f"Amino Acid Conservation Analysis")
    logging.info(f"="*50)
    logging.info(f"Analysis: {analysis}")
    logging.info(f"Paramset: {paramset}")
    logging.info(f"Group: {group}")
    logging.info(f"Input directory: {input_dir}")
    logging.info(f"Output directory: {output_dir}")
    
    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Find trimmed alignment files
    alignment_files = list(input_dir.glob("*_trimmed.fasta"))
    
    if not alignment_files:
        logging.warning(f"No trimmed alignment files found in {input_dir}")
        return
    
    logging.info(f"Found {len(alignment_files)} trimmed alignments to analyze")
    
    # Analyze each alignment
    all_summaries = []
    all_position_data = {}
    
    for align_file in alignment_files:
        gene_name = align_file.stem.replace('_trimmed', '')
        logging.info(f"Analyzing conservation for {gene_name}...")
        
        # Analyze conservation
        result = analyze_alignment_conservation(align_file)
        
        if result:
            summary = result['summary']
            position_data = result['position_data']
            
            all_summaries.append(summary)
            all_position_data[gene_name] = position_data
            
            # Create conservation plot
            plot_file = plots_dir / f"{gene_name}_conservation.png"
            create_conservation_plot(position_data, plot_file, gene_name)
            
            # Save detailed position data
            pos_df = pd.DataFrame(position_data)
            pos_file = output_dir / f"{gene_name}_conservation_details.tsv"
            pos_df.to_csv(pos_file, sep='\t', index=False)
            
            logging.info(f"  {gene_name}: {summary['highly_conserved_positions']} highly conserved positions")
    
    # Create summary report
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        summary_file = output_dir / "conservation_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        
        # Save as JSON for programmatic access
        json_file = output_dir / "conservation_summary.json"
        with open(json_file, 'w') as f:
            json.dump({
                'summary': all_summaries,
                'position_data': all_position_data
            }, f, indent=2)
        
        logging.info(f"Saved conservation summary: {summary_file}")
        logging.info(f"Conservation analysis complete!")
        
        # Print summary stats
        logging.info(f"Summary statistics:")
        logging.info(f"  Total genes analyzed: {len(all_summaries)}")
        logging.info(f"  Mean conservation across all genes: {summary_df['mean_conservation'].mean():.3f}")
        logging.info(f"  Total highly conserved positions: {summary_df['highly_conserved_positions'].sum()}")
    
    else:
        logging.warning("No successful conservation analyses completed")

if __name__ == "__main__":
    main()