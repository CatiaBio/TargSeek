#!/usr/bin/env python3
"""
Extract Conserved Sequences from MSAs
====================================

This script identifies conserved regions (≥5 amino acids) from multiple sequence alignments
and outputs them in TSV format with:
- consensus_sequence: sequence with ties shown as G/A/L
- positions: start-end positions in alignment  
- list_aa_conservation: per-amino acid conservation percentages

Based on the existing analyze_conservation_adaptive.py pipeline.
"""

from Bio import AlignIO
import numpy as np
import pandas as pd
from pathlib import Path
import logging
from collections import Counter
import sys
import json

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_consensus_with_ties_and_percentages(column, min_frequency=0.3):
    """Get consensus with ties and per-AA percentages"""
    # Exclude gaps from consensus calculation (standard practice)
    aa_counts = Counter([aa for aa in column if aa != '-'])
    total_sequences = len(column)
    
    if not aa_counts:
        return '-', {}  # All gaps
    
    total_non_gap = sum(aa_counts.values())
    
    # Calculate percentages for all amino acids (including gaps)
    all_counts = Counter(column)  # Include gaps for percentage calculation
    aa_percentages = {aa: (count / total_sequences * 100) for aa, count in all_counts.items()}
    
    # Find maximum frequency among non-gap amino acids
    max_count = max(aa_counts.values())
    max_freq = max_count / total_non_gap
    
    # Get all amino acids with maximum frequency (ties)
    max_aas = [aa for aa, count in aa_counts.items() if count == max_count]
    
    if max_freq >= min_frequency:
        # Sort for consistent output: G/A/L not A/G/L
        consensus = '/'.join(sorted(max_aas))
        return consensus, aa_percentages
    else:
        return '-', aa_percentages  # No amino acid meets threshold

def analyze_alignment_for_conserved_sequences(alignment_file, gene_name, min_conservation=0.5, min_length=5):
    """Analyze alignment and extract conserved sequences ≥5 amino acids"""
    try:
        # Read alignment
        alignment = AlignIO.read(alignment_file, "fasta")
        
        if len(alignment) < 2:
            logging.warning(f"Alignment {alignment_file} has fewer than 2 sequences")
            return []
        
        alignment_length = alignment.get_alignment_length()
        num_sequences = len(alignment)
        
        logging.info(f"  Analyzing {gene_name}: {num_sequences} sequences, {alignment_length} positions")
        
        # Calculate conservation and consensus for each position
        position_data = []
        
        for i in range(alignment_length):
            column = alignment[:, i]  # Get column (all AAs at position i)
            
            # Calculate Shannon entropy for conservation score
            aa_counts = Counter([aa for aa in column if aa != '-'])
            if aa_counts:
                total = sum(aa_counts.values())
                entropy = 0
                for count in aa_counts.values():
                    if count > 0:
                        freq = count / total
                        entropy -= freq * np.log2(freq)
                
                max_entropy = np.log2(20)  # Maximum possible entropy for 20 amino acids
                conservation_score = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
            else:
                conservation_score = 0  # All gaps
            
            # Get consensus with ties and percentages
            consensus_aa, aa_percentages = get_consensus_with_ties_and_percentages(column)
            
            position_data.append({
                'position': i + 1,  # 1-indexed
                'conservation_score': conservation_score,
                'consensus_aa': consensus_aa,
                'aa_percentages': aa_percentages
            })
        
        # Identify conserved regions (≥min_length amino acids)
        conserved_sequences = []
        current_region = []
        
        for pos_data in position_data:
            if pos_data['conservation_score'] >= min_conservation and pos_data['consensus_aa'] != '-':
                current_region.append(pos_data)
            else:
                # End of a conserved region
                if len(current_region) >= min_length:
                    conserved_sequences.append(extract_conserved_sequence_info(current_region, gene_name))
                current_region = []
        
        # Don't forget the last region
        if len(current_region) >= min_length:
            conserved_sequences.append(extract_conserved_sequence_info(current_region, gene_name))
        
        logging.info(f"    Found {len(conserved_sequences)} conserved sequences (≥{min_length} AA)")
        return conserved_sequences
        
    except Exception as e:
        logging.error(f"Error analyzing {alignment_file}: {e}")
        return []

def extract_conserved_sequence_info(region_data, gene_name):
    """Extract conserved sequence information for TSV output"""
    # Build consensus sequence
    consensus_sequence = ''.join([pos['consensus_aa'] for pos in region_data])
    
    # Get positions (start-end)
    start_pos = region_data[0]['position']
    end_pos = region_data[-1]['position']
    positions = f"{start_pos}-{end_pos}"
    
    # Build conservation list per amino acid
    conservation_list = []
    for pos in region_data:
        aa = pos['consensus_aa']
        percentages = pos['aa_percentages']
        
        if '/' in aa:  # Handle ties: G/A -> "G/A:45%"
            tied_aas = aa.split('/')
            # Get the percentage of the tied amino acids (they should be the same)
            tied_percentages = [f"{tied_aa}:{percentages.get(tied_aa, 0):.1f}%" for tied_aa in tied_aas]
            conservation_list.append('/'.join(tied_percentages))
        else:
            # Single amino acid
            percentage = percentages.get(aa, 0)
            conservation_list.append(f"{aa}:{percentage:.1f}%")
    
    return {
        'gene': gene_name,
        'consensus_sequence': consensus_sequence,
        'positions': positions,
        'list_aa_conservation': ','.join(conservation_list),
        'length': len(consensus_sequence),
        'avg_conservation': np.mean([pos['conservation_score'] for pos in region_data])
    }

def process_alignment_directory(alignment_dir, output_file, min_conservation=0.5, min_length=5):
    """Process all alignments in a directory"""
    if not alignment_dir.exists():
        logging.warning(f"Alignment directory not found: {alignment_dir}")
        return
    
    logging.info(f"Processing alignments from: {alignment_dir}")
    
    # Find all FASTA alignment files
    alignment_files = list(alignment_dir.glob("*.fasta"))
    
    if not alignment_files:
        logging.warning(f"No FASTA files found in {alignment_dir}")
        return
    
    logging.info(f"Found {len(alignment_files)} alignment files")
    
    all_conserved_sequences = []
    
    for alignment_file in alignment_files:
        gene_name = alignment_file.stem.replace('_trimmed', '').replace('_aligned', '')
        
        try:
            conserved_seqs = analyze_alignment_for_conserved_sequences(
                alignment_file, gene_name, min_conservation, min_length
            )
            all_conserved_sequences.extend(conserved_seqs)
            
        except Exception as e:
            logging.error(f"Error processing {alignment_file}: {e}")
            continue
    
    # Save results to TSV
    if all_conserved_sequences:
        df = pd.DataFrame(all_conserved_sequences)
        
        # Sort by gene name and average conservation (descending)
        df = df.sort_values(['gene', 'avg_conservation'], ascending=[True, False])
        
        # Select columns for output (exclude avg_conservation and length - they were for sorting)
        output_columns = ['gene', 'consensus_sequence', 'positions', 'list_aa_conservation']
        df[output_columns].to_csv(output_file, sep='\t', index=False)
        
        logging.info(f"Saved {len(all_conserved_sequences)} conserved sequences to {output_file}")
        
        # Log summary statistics
        logging.info(f"Summary:")
        logging.info(f"  Total conserved sequences: {len(all_conserved_sequences)}")
        logging.info(f"  Genes with conserved sequences: {df['gene'].nunique()}")
        logging.info(f"  Average sequence length: {df['length'].mean():.1f} AA")
        logging.info(f"  Length range: {df['length'].min()}-{df['length'].max()} AA")
        
    else:
        logging.warning("No conserved sequences found")

def main():
    """Main function for Snakemake integration or standalone use"""
    
    try:
        # Get processing options
        enable_no_3d = getattr(snakemake.params, 'enable_no_3d', True)
        enable_with_3d = getattr(snakemake.params, 'enable_with_3d', True)
        
        min_conservation = getattr(snakemake.params, 'min_conservation', 0.5)
        min_length = getattr(snakemake.params, 'min_length', 5)
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        # Get inputs and outputs based on enabled processing
        inputs = {}
        outputs = {}
        
        if enable_no_3d:
            if hasattr(snakemake.input, 'raw_main'):
                inputs['raw_main'] = Path(snakemake.input.raw_main)
            if hasattr(snakemake.input, 'trimmed_main'):
                inputs['trimmed_main'] = Path(snakemake.input.trimmed_main)
            if hasattr(snakemake.output, 'raw_sequences_main'):
                outputs['raw_sequences_main'] = Path(snakemake.output.raw_sequences_main)
            if hasattr(snakemake.output, 'trimmed_sequences_main'):
                outputs['trimmed_sequences_main'] = Path(snakemake.output.trimmed_sequences_main)
        
        if enable_with_3d:
            if hasattr(snakemake.input, 'raw_structure'):
                inputs['raw_structure'] = Path(snakemake.input.raw_structure)
            if hasattr(snakemake.input, 'trimmed_structure'):
                inputs['trimmed_structure'] = Path(snakemake.input.trimmed_structure)
            if hasattr(snakemake.output, 'raw_sequences_structure'):
                outputs['raw_sequences_structure'] = Path(snakemake.output.raw_sequences_structure)
            if hasattr(snakemake.output, 'trimmed_sequences_structure'):
                outputs['trimmed_sequences_structure'] = Path(snakemake.output.trimmed_sequences_structure)
        
    except NameError:
        # Standalone mode
        if len(sys.argv) < 3:
            print("Usage: python extract_conserved_sequences.py <alignment_dir> <output_file> [min_conservation] [min_length]")
            print("Example: python extract_conserved_sequences.py results/msa_trimmed conserved_sequences.tsv 0.5 5")
            sys.exit(1)
        
        alignment_dir = Path(sys.argv[1])
        output_file = Path(sys.argv[2])
        min_conservation = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
        min_length = int(sys.argv[4]) if len(sys.argv) > 4 else 5
        
        analysis = "standalone"
        paramset = "default"
        group = "all"
    
    logging.info(f"Conserved Sequence Extraction")
    logging.info(f"="*50)
    logging.info(f"Analysis: {analysis}")
    logging.info(f"Paramset: {paramset}")
    logging.info(f"Group: {group}")
    logging.info(f"Minimum conservation score: {min_conservation}")
    logging.info(f"Minimum sequence length: {min_length} amino acids")
    
    # Process alignments based on enabled processing options
    if 'snakemake' in globals():
        # Snakemake mode - process only enabled combinations
        datasets = []
        
        if enable_no_3d:
            logging.info(f"No-3D processing enabled")
            if 'raw_main' in inputs and 'raw_sequences_main' in outputs:
                datasets.append((inputs['raw_main'], outputs['raw_sequences_main'], "raw main"))
            if 'trimmed_main' in inputs and 'trimmed_sequences_main' in outputs:
                datasets.append((inputs['trimmed_main'], outputs['trimmed_sequences_main'], "trimmed main"))
        else:
            logging.info(f"No-3D processing disabled")
        
        if enable_with_3d:
            logging.info(f"With-3D processing enabled")
            if 'raw_structure' in inputs and 'raw_sequences_structure' in outputs:
                datasets.append((inputs['raw_structure'], outputs['raw_sequences_structure'], "raw structure"))
            if 'trimmed_structure' in inputs and 'trimmed_sequences_structure' in outputs:
                datasets.append((inputs['trimmed_structure'], outputs['trimmed_sequences_structure'], "trimmed structure"))
        else:
            logging.info(f"With-3D processing disabled")
        
        if not datasets:
            logging.warning("No processing enabled - both enable_no_3d and enable_with_3d are False")
            return
        
        for input_dir, output_file, description in datasets:
            logging.info(f"Processing {description} alignments: {input_dir}")
            logging.info(f"Output: {output_file}")
            output_file.parent.mkdir(parents=True, exist_ok=True)
            process_alignment_directory(input_dir, output_file, min_conservation, min_length)
    else:
        # Standalone mode - process single directory
        logging.info(f"Input directory: {alignment_dir}")
        logging.info(f"Output file: {output_file}")
        output_file.parent.mkdir(parents=True, exist_ok=True)
        process_alignment_directory(alignment_dir, output_file, min_conservation, min_length)
    
    logging.info("Conserved sequence extraction complete!")

if __name__ == "__main__":
    main()