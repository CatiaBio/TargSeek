#!/usr/bin/env python3
"""
Alignment Quality Assessment - Before/After Comparison
====================================================

This script compares alignment quality before and after trimAl processing:
- Basic statistics (length, gaps, sequences)
- Conservation analysis
- Gap distribution analysis
- Quality scoring and improvement metrics
- Generates comparison plots and reports
"""

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import logging
import json

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_basic_stats(alignment):
    """Calculate basic alignment statistics"""
    stats = {
        'num_sequences': len(alignment),
        'alignment_length': alignment.get_alignment_length(),
        'total_positions': len(alignment) * alignment.get_alignment_length()
    }
    
    # Calculate gap statistics
    gap_count = 0
    total_chars = 0
    
    for record in alignment:
        sequence = str(record.seq)
        gap_count += sequence.count('-')
        total_chars += len(sequence)
    
    stats['gap_percentage'] = (gap_count / total_chars) * 100 if total_chars > 0 else 0
    stats['non_gap_positions'] = total_chars - gap_count
    
    return stats

def calculate_conservation_metrics(alignment):
    """Calculate conservation metrics for alignment"""
    if len(alignment) == 0:
        return {'mean_conservation': 0, 'conservation_positions': 0}
    
    length = alignment.get_alignment_length()
    conservation_scores = []
    
    for i in range(length):
        column = alignment[:, i]
        
        # Count non-gap characters
        non_gap_chars = [aa for aa in column if aa != '-']
        if len(non_gap_chars) == 0:
            continue
            
        # Calculate conservation as frequency of most common amino acid
        aa_counts = Counter(non_gap_chars)
        most_common_count = aa_counts.most_common(1)[0][1]
        conservation = most_common_count / len(non_gap_chars)
        conservation_scores.append(conservation)
    
    return {
        'mean_conservation': np.mean(conservation_scores) if conservation_scores else 0,
        'conservation_positions': len([s for s in conservation_scores if s > 0.8])
    }

def calculate_quality_score(stats, conservation):
    """Calculate overall quality score for alignment"""
    # Base score from conservation
    conservation_score = conservation['mean_conservation'] * 100
    
    # Penalty for excessive gaps
    gap_penalty = max(0, stats['gap_percentage'] - 20) * 2  # Penalty starts at 20% gaps
    
    # Bonus for having sequences
    sequence_bonus = min(stats['num_sequences'] * 2, 20)  # Max 20 points
    
    quality_score = conservation_score + sequence_bonus - gap_penalty
    return max(0, min(100, quality_score))  # Clamp to 0-100

def analyze_single_alignment(alignment_file, file_type):
    """Analyze a single alignment file"""
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        
        basic_stats = calculate_basic_stats(alignment)
        conservation = calculate_conservation_metrics(alignment)
        quality_score = calculate_quality_score(basic_stats, conservation)
        
        return {
            'gene': alignment_file.stem.replace('_aligned', '').replace('_trimmed', ''),
            'file_type': file_type,
            'quality_score': quality_score,
            **basic_stats,
            **conservation
        }
        
    except Exception as e:
        logging.error(f"Error analyzing {alignment_file}: {e}")
        return None

def create_comparison_plots(comparison_df, output_dir):
    """Create before/after comparison plots"""
    
    # Set up the plot style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. Quality Score Comparison
    plt.figure(figsize=(12, 8))
    
    # Subplot 1: Quality scores
    plt.subplot(2, 2, 1)
    raw_scores = comparison_df[comparison_df['file_type'] == 'raw']['quality_score']
    trimmed_scores = comparison_df[comparison_df['file_type'] == 'trimmed']['quality_score']
    
    plt.boxplot([raw_scores, trimmed_scores], labels=['Raw MAFFT', 'After trimAl'])
    plt.ylabel('Quality Score')
    plt.title('Alignment Quality Comparison')
    plt.grid(True, alpha=0.3)
    
    # Subplot 2: Gap percentage
    plt.subplot(2, 2, 2)
    raw_gaps = comparison_df[comparison_df['file_type'] == 'raw']['gap_percentage']
    trimmed_gaps = comparison_df[comparison_df['file_type'] == 'trimmed']['gap_percentage']
    
    plt.boxplot([raw_gaps, trimmed_gaps], labels=['Raw MAFFT', 'After trimAl'])
    plt.ylabel('Gap Percentage (%)')
    plt.title('Gap Content Comparison')
    plt.grid(True, alpha=0.3)
    
    # Subplot 3: Conservation
    plt.subplot(2, 2, 3)
    raw_cons = comparison_df[comparison_df['file_type'] == 'raw']['mean_conservation']
    trimmed_cons = comparison_df[comparison_df['file_type'] == 'trimmed']['mean_conservation']
    
    plt.boxplot([raw_cons, trimmed_cons], labels=['Raw MAFFT', 'After trimAl'])
    plt.ylabel('Mean Conservation')
    plt.title('Conservation Comparison')
    plt.grid(True, alpha=0.3)
    
    # Subplot 4: Alignment length
    plt.subplot(2, 2, 4)
    raw_length = comparison_df[comparison_df['file_type'] == 'raw']['alignment_length']
    trimmed_length = comparison_df[comparison_df['file_type'] == 'trimmed']['alignment_length']
    
    plt.boxplot([raw_length, trimmed_length], labels=['Raw MAFFT', 'After trimAl'])
    plt.ylabel('Alignment Length')
    plt.title('Length Comparison')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "quality_comparison_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Gene-by-gene improvement plot
    plt.figure(figsize=(14, 8))
    
    # Merge raw and trimmed data for each gene
    raw_df = comparison_df[comparison_df['file_type'] == 'raw'].set_index('gene')
    trimmed_df = comparison_df[comparison_df['file_type'] == 'trimmed'].set_index('gene')
    
    common_genes = set(raw_df.index) & set(trimmed_df.index)
    if common_genes:
        improvements = []
        gene_names = []
        
        for gene in common_genes:
            raw_score = raw_df.loc[gene, 'quality_score']
            trimmed_score = trimmed_df.loc[gene, 'quality_score']
            improvement = trimmed_score - raw_score
            improvements.append(improvement)
            gene_names.append(gene)
        
        # Sort by improvement
        sorted_data = sorted(zip(gene_names, improvements), key=lambda x: x[1], reverse=True)
        gene_names, improvements = zip(*sorted_data)
        
        # Color bars based on improvement
        colors = ['green' if imp > 0 else 'red' if imp < -5 else 'orange' for imp in improvements]
        
        plt.bar(range(len(gene_names)), improvements, color=colors, alpha=0.7)
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        plt.xlabel('Genes')
        plt.ylabel('Quality Score Improvement')
        plt.title('trimAl Improvement per Gene (Green=Improved, Red=Degraded, Orange=Minimal change)')
        plt.xticks(range(len(gene_names)), gene_names, rotation=45)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_dir / "gene_improvement_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()

def add_improvement_classification(comparison_df):
    """Add improvement classification column to comparison DataFrame"""
    
    # Separate raw and trimmed data
    raw_df = comparison_df[comparison_df['file_type'] == 'raw'].set_index('gene')
    trimmed_df = comparison_df[comparison_df['file_type'] == 'trimmed'].set_index('gene')
    
    # Add improvement columns to the original DataFrame
    improvement_classification = []
    quality_improvement = []
    
    for _, row in comparison_df.iterrows():
        gene = row['gene']
        file_type = row['file_type']
        
        # Only calculate improvement for trimmed sequences
        if file_type == "trimmed" and gene in raw_df.index and gene in trimmed_df.index:
            raw_quality = raw_df.loc[gene, 'quality_score']
            trimmed_quality = trimmed_df.loc[gene, 'quality_score']
            improvement = trimmed_quality - raw_quality
            
            # Classify improvement
            if improvement > 5:
                classification = "significantly_improved"
            elif improvement > 1:
                classification = "improved"
            elif improvement > -1:
                classification = "no_change"
            elif improvement > -5:
                classification = "slightly_degraded"
            else:
                classification = "significantly_degraded"
            
            quality_improvement.append(improvement)
            improvement_classification.append(classification)
        else:
            # For raw sequences or unpaired sequences, leave empty
            quality_improvement.append(None)
            improvement_classification.append(None)
    
    # Add new columns
    comparison_df['trimming_effect'] = improvement_classification
    comparison_df['quality_improvement'] = quality_improvement
    
    return comparison_df

def process_alignment_set(raw_dir, trimmed_dir, output_dir, set_name):
    """Process a set of alignments (no-3D or with-3D)"""
    
    logging.info(f"\nProcessing {set_name} alignments:")
    logging.info(f"  Raw alignments: {raw_dir}")
    logging.info(f"  Trimmed alignments: {trimmed_dir}")
    
    if not raw_dir.exists():
        logging.warning(f"  Raw directory does not exist: {raw_dir}")
        return [], []
    
    if not trimmed_dir.exists():
        logging.warning(f"  Trimmed directory does not exist: {trimmed_dir}")
        return [], []
    
    # Find alignment files - look for *.fasta files (not specific patterns)
    raw_files = list(raw_dir.glob("*.fasta"))
    trimmed_files = list(trimmed_dir.glob("*.fasta"))
    
    logging.info(f"  Found {len(raw_files)} raw alignments and {len(trimmed_files)} trimmed alignments")
    
    # Analyze all alignments
    all_results = []
    
    # Analyze raw alignments
    for raw_file in raw_files:
        result = analyze_single_alignment(raw_file, f"raw_{set_name}")
        if result:
            all_results.append(result)
    
    # Analyze trimmed alignments
    for trimmed_file in trimmed_files:
        result = analyze_single_alignment(trimmed_file, f"trimmed_{set_name}")
        if result:
            all_results.append(result)
    
    return all_results, raw_files + trimmed_files

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get parameters from Snakemake - handle 4 directory inputs
        raw_main = Path(snakemake.input.main)
        raw_structure = Path(snakemake.input.structure)
        trimmed_main = Path(snakemake.input.trim_main)
        trimmed_structure = Path(snakemake.input.trim_structure)
        
        output_main = Path(snakemake.output.main)
        output_structure = Path(snakemake.output.structure)
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        use_3d = snakemake.params.get('use_3d', 'no_3d')
        
    except NameError:
        # Test mode
        import sys
        if len(sys.argv) != 7:
            print("Usage: python assess_alignment_quality_comparison.py <raw_main> <raw_structure> <trim_main> <trim_structure> <output_main> <output_structure>")
            sys.exit(1)
            
        raw_main = Path(sys.argv[1])
        raw_structure = Path(sys.argv[2])
        trimmed_main = Path(sys.argv[3])
        trimmed_structure = Path(sys.argv[4])
        output_main = Path(sys.argv[5])
        output_structure = Path(sys.argv[6])
        
        analysis = "test"
        paramset = "test"
        group = "test"
        use_3d = "both"
    
    logging.info(f"Alignment Quality Comparison Analysis")
    logging.info(f"="*50)
    logging.info(f"Analysis: {analysis}")
    logging.info(f"Paramset: {paramset}")
    logging.info(f"Group: {group}")
    
    # Create output directories
    output_main.mkdir(parents=True, exist_ok=True)
    output_structure.mkdir(parents=True, exist_ok=True)
    
    # Process main alignments
    main_results, main_files = process_alignment_set(raw_main, trimmed_main, output_main, "main")
    
    # Process structure alignments
    structure_results, structure_files = process_alignment_set(raw_structure, trimmed_structure, output_structure, "structure")
    
    # Generate reports for each set separately
    all_sets = [
        (main_results, output_main, "main"),
        (structure_results, output_structure, "structure")
    ]
    
    for results, output_dir, set_name in all_sets:
        if not results:
            logging.warning(f"No successful analyses completed for {set_name}")
            continue
            
        logging.info(f"\nGenerating reports for {set_name}...")
        
        # Create comparison DataFrame
        comparison_df = pd.DataFrame(results)
        
        # Add improvement classification column
        comparison_df = add_improvement_classification(comparison_df)
        
        # Save detailed results
        comparison_file = output_dir / "quality_comparison.tsv"
        comparison_df.to_csv(comparison_file, sep='\t', index=False)
        
        # Create summary statistics for this set
        summary_stats = {}
        for file_type_prefix in [f'raw_{set_name}', f'trimmed_{set_name}']:
            subset = comparison_df[comparison_df['file_type'] == file_type_prefix]
            if len(subset) > 0:
                summary_stats[file_type_prefix] = {
                    'mean_quality_score': subset['quality_score'].mean(),
                    'mean_gap_percentage': subset['gap_percentage'].mean(),
                    'mean_conservation': subset['mean_conservation'].mean(),
                    'mean_alignment_length': subset['alignment_length'].mean(),
                    'num_alignments': len(subset)
                }
        
        # Calculate improvements for this set
        raw_key = f'raw_{set_name}'
        trimmed_key = f'trimmed_{set_name}'
        
        if raw_key in summary_stats and trimmed_key in summary_stats:
            improvement = {
                'quality_improvement': summary_stats[trimmed_key]['mean_quality_score'] - summary_stats[raw_key]['mean_quality_score'],
                'gap_reduction': summary_stats[raw_key]['mean_gap_percentage'] - summary_stats[trimmed_key]['mean_gap_percentage'],
                'conservation_change': summary_stats[trimmed_key]['mean_conservation'] - summary_stats[raw_key]['mean_conservation'],
                'length_reduction': summary_stats[raw_key]['mean_alignment_length'] - summary_stats[trimmed_key]['mean_alignment_length']
            }
            summary_stats['improvement'] = improvement
        
        # Save summary for this set
        summary_file = output_dir / "quality_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary_stats, f, indent=2)
        
        # Create comparison plots for this set
        create_comparison_plots(comparison_df, output_dir)
        
        # Log summary for this set
        logging.info(f"Quality comparison complete for {set_name}!")
        if 'improvement' in summary_stats:
            imp = summary_stats['improvement']
            logging.info(f"trimAl Results for {set_name}:")
            logging.info(f"  Quality improvement: {imp['quality_improvement']:+.1f} points")
            logging.info(f"  Gap reduction: {imp['gap_reduction']:+.1f}%")
            logging.info(f"  Conservation change: {imp['conservation_change']:+.3f}")
            logging.info(f"  Length reduction: {imp['length_reduction']:+.0f} positions")
        
        # Log improvement classification counts for this set
        if 'trimming_effect' in comparison_df.columns:
            effect_counts = comparison_df['trimming_effect'].value_counts()
            logging.info(f"Improvement classification for {set_name}:")
            for effect, count in effect_counts.items():
                logging.info(f"  {effect}: {count} genes")
    
    logging.info(f"\n=== Quality Assessment Analysis Complete ===")
    logging.info(f"Reports generated for both no-3D and with-3D alignments")

if __name__ == "__main__":
    main()