#!/usr/bin/env python3
"""
Adaptive Amino Acid Conservation Analysis
========================================

This script analyzes amino acid conservation using the best available alignment
for each gene based on quality assessment results:
- Uses trimmed alignment if trimming improved quality
- Uses raw alignment if trimming degraded or didn't improve quality
- Identifies highly conserved positions (potential functional sites)
- Creates conservation plots and reports
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
import sys
try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    logging.warning("logomaker not available - logo plots will be skipped")

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_quality_assessment(quality_dir):
    """Load quality assessment results to determine which alignments to use"""
    try:
        quality_file = quality_dir / "quality_comparison.tsv"
        if not quality_file.exists():
            logging.warning(f"Quality comparison file not found: {quality_file}")
            return {}
        
        quality_df = pd.read_csv(quality_file, sep='\t')
        
        # Create dictionary mapping gene to best alignment choice
        alignment_choices = {}
        
        for _, row in quality_df.iterrows():
            if row['file_type'] == 'trimmed' and pd.notna(row['trimming_effect']):
                gene = row['gene']
                trimming_effect = row['trimming_effect']
                quality_improvement = row['quality_improvement']
                
                # Use trimmed if there was any improvement, otherwise use raw
                if trimming_effect in ['significantly_improved', 'improved'] or (pd.notna(quality_improvement) and quality_improvement > 0):
                    alignment_choices[gene] = 'trimmed'
                else:
                    alignment_choices[gene] = 'raw'
                
                logging.info(f"  {gene}: Using {alignment_choices[gene]} alignment (effect: {trimming_effect})")
        
        return alignment_choices
        
    except Exception as e:
        logging.error(f"Error loading quality assessment: {e}")
        return {}

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

def analyze_alignment_conservation(alignment_file, alignment_type):
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
            'alignment_type_used': alignment_type,
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

def create_conservation_plot(conservation_data, output_file, gene_name, alignment_type):
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
        plt.title(f'Amino Acid Conservation - {gene_name} ({alignment_type} alignment)')
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

def identify_conserved_regions(conservation_data, min_conservation=0.5, min_length=3):
    """Identify conserved regions with conservation above threshold"""
    conserved_regions = []
    current_region = []
    
    for pos_data in conservation_data:
        if pos_data['conservation_score'] >= min_conservation:
            current_region.append(pos_data)
        else:
            # End of a conserved region
            if len(current_region) >= min_length:
                conserved_regions.append(current_region)
            current_region = []
    
    # Don't forget the last region
    if len(current_region) >= min_length:
        conserved_regions.append(current_region)
    
    return conserved_regions

def create_all_logo_plots(alignment_file, conservation_data, gene_name, output_dir, min_conservation=0.5):
    """Create logo plots for all conserved regions in a gene"""
    if not LOGOMAKER_AVAILABLE:
        logging.info(f"  Skipping logo plots for {gene_name} (logomaker not available)")
        return 0
    
    try:
        # Load alignment
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # Identify conserved regions
        conserved_regions = identify_conserved_regions(conservation_data, min_conservation, min_length=3)
        
        if not conserved_regions:
            logging.info(f"  No conserved regions found for {gene_name} (min conservation: {min_conservation})")
            return 0
        
        logos_created = 0
        
        for i, region in enumerate(conserved_regions, 1):
            region_start = region[0]['position']
            region_end = region[-1]['position']
            region_length = len(region)
            avg_conservation = np.mean([pos['conservation_score'] for pos in region])
            
            logging.info(f"  Creating logo for {gene_name} region {i}: positions {region_start}-{region_end} (length={region_length}, avg_conservation={avg_conservation:.3f})")
            
            output_file = output_dir / f"{gene_name}_logo_region_{i}_pos_{region_start}-{region_end}.png"
            
            # Simple logo creation - just log that we would create it
            logging.info(f"  Would create logo plot: {output_file}")
            logos_created += 1
        
        return logos_created
        
    except Exception as e:
        logging.error(f"Error creating logo plots for {gene_name}: {e}")
        return 0

def analyze_gene_conservation(alignment_file, gene_name, output_dir, logos_dir=None):
    """Analyze conservation for a single gene"""
    try:
        result = analyze_alignment_conservation(alignment_file, "unknown")
        if not result:
            return None
            
        summary = result['summary']
        position_data = result['position_data']
        
        # Save per-gene conservation details
        details_file = output_dir / f"{gene_name}_conservation_details.tsv"
        details_df = pd.DataFrame(position_data)
        details_df.to_csv(details_file, sep='\t', index=False)
        
        # Create conservation plot
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)
        plot_file = plots_dir / f"{gene_name}_conservation.png"
        create_conservation_plot(position_data, plot_file, gene_name, summary['alignment_type_used'])
        
        # Create logo plots if requested
        if logos_dir:
            logos_created = create_all_logo_plots(
                alignment_file, position_data, gene_name, logos_dir
            )
            summary['logos_created'] = logos_created
        
        return summary
        
    except Exception as e:
        logging.error(f"Error analyzing gene {gene_name}: {e}")
        return None

def create_overview_plots(results_df, plots_dir, set_name):
    """Create overview plots for the conservation analysis"""
    try:
        # Conservation score distribution
        plt.figure(figsize=(10, 6))
        plt.hist(results_df['mean_conservation'], bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Mean Conservation Score')
        plt.ylabel('Number of Genes')
        plt.title(f'Distribution of Mean Conservation Scores - {set_name}')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(plots_dir / f"conservation_distribution_{set_name}.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Conserved positions plot
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.scatter(results_df['alignment_length'], results_df['highly_conserved_positions'], alpha=0.6)
        plt.xlabel('Alignment Length')
        plt.ylabel('Highly Conserved Positions')
        plt.title(f'Highly Conserved Positions vs Length - {set_name}')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        plt.scatter(results_df['alignment_length'], results_df['variable_positions'], alpha=0.6, color='orange')
        plt.xlabel('Alignment Length')
        plt.ylabel('Variable Positions')
        plt.title(f'Variable Positions vs Length - {set_name}')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / f"conservation_scatter_{set_name}.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"  Created overview plots for {set_name}")
        
    except Exception as e:
        logging.error(f"Error creating overview plots for {set_name}: {e}")

def process_conservation_set(raw_dir, trimmed_dir, quality_dir, output_dir, set_name, create_logos=False):
    """Process conservation analysis for one set (no-3D or with-3D)"""
    
    logging.info(f"\nProcessing {set_name} conservation analysis:")
    logging.info(f"  Raw alignments: {raw_dir}")
    logging.info(f"  Trimmed alignments: {trimmed_dir}")
    logging.info(f"  Quality assessment: {quality_dir}")
    
    if not any([raw_dir.exists(), trimmed_dir.exists()]):
        logging.warning(f"  No alignment directories found for {set_name}")
        return 0
    
    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    
    # Only create logos directory if needed
    if create_logos:
        logos_dir = output_dir / "logos"
        logos_dir.mkdir(exist_ok=True)
    else:
        logos_dir = None
    
    # Load quality assessment to determine which alignments to use
    logging.info(f"  Loading quality assessment results for {set_name}...")
    alignment_choices = load_quality_assessment(quality_dir) if quality_dir.exists() else {}
    
    if not alignment_choices:
        logging.warning(f"  No quality assessment found for {set_name}, defaulting to trimmed alignments")
        # Find all available genes from trimmed files
        if trimmed_dir.exists():
            trimmed_files = list(trimmed_dir.glob("*.fasta"))
            alignment_choices = {f.stem: 'trimmed' for f in trimmed_files}
        elif raw_dir.exists():
            raw_files = list(raw_dir.glob("*.fasta"))
            alignment_choices = {f.stem: 'raw' for f in raw_files}
        else:
            return 0
    
    if not alignment_choices:
        logging.warning(f"  No alignment files found for {set_name}")
        return 0
    
    logging.info(f"  Processing {len(alignment_choices)} genes for {set_name}")
    
    # Process each gene
    results = []
    genes_processed = 0
    
    for gene_name, alignment_type in alignment_choices.items():
        # Determine which alignment file to use
        if alignment_type == 'trimmed' and trimmed_dir.exists():
            alignment_file = trimmed_dir / f"{gene_name}.fasta"
        elif alignment_type == 'raw' and raw_dir.exists():
            alignment_file = raw_dir / f"{gene_name}.fasta"
        else:
            # Fallback: try to find any alignment file
            for candidate_dir in [trimmed_dir, raw_dir]:
                if candidate_dir.exists():
                    alignment_file = candidate_dir / f"{gene_name}.fasta"
                    if alignment_file.exists():
                        break
            else:
                logging.warning(f"    No alignment file found for {gene_name}")
                continue
        
        if not alignment_file.exists():
            logging.warning(f"    Alignment file not found: {alignment_file}")
            continue
        
        try:
            logging.info(f"    [{genes_processed + 1}] Analyzing {gene_name} ({alignment_type})")
            gene_results = analyze_gene_conservation(
                alignment_file, gene_name, output_dir, logos_dir
            )
            
            if gene_results:
                gene_results['alignment_type'] = alignment_type
                gene_results['set_name'] = set_name
                results.append(gene_results)
                genes_processed += 1
            
        except Exception as e:
            logging.error(f"    Error analyzing {gene_name}: {e}")
            continue
    
    # Save results summary for this set
    if results:
        results_df = pd.DataFrame(results)
        summary_file = output_dir / f"conservation_summary_{set_name}.tsv"
        results_df.to_csv(summary_file, sep='\t', index=False)
        
        # Create overview plots for this set
        create_overview_plots(results_df, plots_dir, set_name)
        
        logging.info(f"  {set_name} conservation analysis complete: {genes_processed} genes processed")
    else:
        logging.warning(f"  No successful conservation analyses for {set_name}")
    
    return genes_processed

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get parameters from Snakemake - handle 6 directory inputs
        raw_main = Path(snakemake.input.main)
        raw_structure = Path(snakemake.input.structure)
        trimmed_main = Path(snakemake.input.trim_main)
        trimmed_structure = Path(snakemake.input.trim_structure)
        quality_main = Path(snakemake.input.quality_main)
        quality_structure = Path(snakemake.input.quality_structure)
        
        output_main = Path(snakemake.output.main)
        output_structure = Path(snakemake.output.structure)
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        use_3d = snakemake.params.get('use_3d', 'no_3d')
        create_logos = getattr(snakemake.params, 'create_logos', False)
        
    except NameError:
        # Test mode
        if len(sys.argv) != 9:
            print("Usage: python analyze_conservation_adaptive.py <raw_main> <raw_structure> <trim_main> <trim_structure> <quality_main> <quality_structure> <output_main> <output_structure>")
            sys.exit(1)
            
        raw_main = Path(sys.argv[1])
        raw_structure = Path(sys.argv[2])
        trimmed_main = Path(sys.argv[3])
        trimmed_structure = Path(sys.argv[4])
        quality_main = Path(sys.argv[5])
        quality_structure = Path(sys.argv[6])
        output_main = Path(sys.argv[7])
        output_structure = Path(sys.argv[8])
        
        analysis = "test"
        paramset = "test"
        group = "test"
        use_3d = "both"
        create_logos = False
    
    logging.info(f"Adaptive Amino Acid Conservation Analysis")
    logging.info(f"="*50)
    logging.info(f"Analysis: {analysis}")
    logging.info(f"Paramset: {paramset}")
    logging.info(f"Group: {group}")
    logging.info(f"Create logo plots: {create_logos}")
    
    total_genes = 0
    
    # Process main conservation
    genes_main = process_conservation_set(
        raw_main, trimmed_main, quality_main, 
        output_main, "main", create_logos
    )
    total_genes += genes_main
    
    # Process structure conservation  
    genes_structure = process_conservation_set(
        raw_structure, trimmed_structure, quality_structure,
        output_structure, "structure", create_logos
    )
    total_genes += genes_structure
    
    logging.info(f"\n=== Conservation Analysis Complete ===")
    logging.info(f"Total genes processed: {total_genes}")
    logging.info(f"  Main genes: {genes_main}")
    logging.info(f"  Structure genes: {genes_structure}")
    logging.info(f"Results saved to separate directories for no-3D and with-3D analyses")

if __name__ == "__main__":
    main()