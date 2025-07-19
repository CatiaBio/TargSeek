#!/usr/bin/env python3
"""
Multiple Sequence Alignment Quality Assessment
=============================================

This script assesses the quality of MAFFT alignments using multiple metrics:
- Basic statistics (length, gaps, sequences)
- Conservation analysis
- Gap distribution analysis
- Quality scoring and recommendations
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
    gap_counts = []
    seq_lengths = []
    
    for record in alignment:
        seq_str = str(record.seq)
        gap_count = seq_str.count('-')
        gap_counts.append(gap_count)
        seq_lengths.append(len(seq_str) - gap_count)
    
    stats.update({
        'total_gaps': sum(gap_counts),
        'gap_percentage': (sum(gap_counts) / stats['total_positions']) * 100,
        'avg_gaps_per_seq': np.mean(gap_counts),
        'max_gaps_per_seq': max(gap_counts),
        'avg_sequence_length': np.mean(seq_lengths),
        'min_sequence_length': min(seq_lengths),
        'max_sequence_length': max(seq_lengths)
    })
    
    return stats

def calculate_conservation_scores(alignment):
    """Calculate per-position conservation scores"""
    summary_align = AlignInfo.SummaryInfo(alignment)
    
    # Calculate conservation using Shannon entropy
    conservation_scores = []
    consensus_seq = ""
    
    for i in range(alignment.get_alignment_length()):
        column = [record.seq[i] for record in alignment]
        
        # Remove gaps for conservation calculation
        column_no_gaps = [aa for aa in column if aa != '-']
        
        if len(column_no_gaps) == 0:
            conservation_scores.append(0)
            consensus_seq += "-"
            continue
        
        # Calculate Shannon entropy
        aa_counts = Counter(column_no_gaps)
        total = len(column_no_gaps)
        entropy = 0
        
        for count in aa_counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)
        
        # Convert entropy to conservation score (0-1, higher = more conserved)
        max_entropy = np.log2(min(20, len(column_no_gaps)))  # 20 amino acids
        conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
        conservation_scores.append(conservation)
        
        # Get most common amino acid for consensus
        most_common = aa_counts.most_common(1)[0][0]
        consensus_seq += most_common
    
    return conservation_scores, consensus_seq

def analyze_gap_distribution(alignment):
    """Analyze gap distribution patterns"""
    gap_analysis = {
        'gap_columns': 0,
        'fully_conserved_columns': 0,
        'partially_gapped_columns': 0,
        'gap_runs': []
    }
    
    for i in range(alignment.get_alignment_length()):
        column = [record.seq[i] for record in alignment]
        gap_count = column.count('-')
        
        if gap_count == len(column):
            gap_analysis['gap_columns'] += 1
        elif gap_count == 0:
            gap_analysis['fully_conserved_columns'] += 1
        else:
            gap_analysis['partially_gapped_columns'] += 1
    
    # Analyze gap runs (consecutive gaps) for each sequence
    for record in alignment:
        seq_str = str(record.seq)
        gap_runs = []
        current_run = 0
        
        for char in seq_str:
            if char == '-':
                current_run += 1
            else:
                if current_run > 0:
                    gap_runs.append(current_run)
                    current_run = 0
        
        if current_run > 0:
            gap_runs.append(current_run)
        
        gap_analysis['gap_runs'].extend(gap_runs)
    
    return gap_analysis

def assess_alignment_quality(alignment_file, output_dir):
    """Comprehensive alignment quality assessment"""
    
    try:
        # Read alignment
        alignment = AlignIO.read(alignment_file, "fasta")
        gene_name = Path(alignment_file).stem.replace('_aligned', '')
        
        logging.info(f"Assessing alignment quality for gene: {gene_name}")
        
        # Calculate all metrics
        basic_stats = calculate_basic_stats(alignment)
        conservation_scores, consensus_seq = calculate_conservation_scores(alignment)
        gap_analysis = analyze_gap_distribution(alignment)
        
        # Overall quality assessment
        quality_score = calculate_overall_quality(basic_stats, conservation_scores, gap_analysis)
        
        # Compile results
        results = {
            'gene': gene_name,
            'basic_stats': basic_stats,
            'conservation': {
                'mean_conservation': np.mean(conservation_scores),
                'min_conservation': np.min(conservation_scores),
                'max_conservation': np.max(conservation_scores),
                'highly_conserved_positions': sum(1 for score in conservation_scores if score > 0.8),
                'poorly_conserved_positions': sum(1 for score in conservation_scores if score < 0.3)
            },
            'gap_analysis': gap_analysis,
            'quality_score': quality_score,
            'recommendations': generate_recommendations(basic_stats, conservation_scores, gap_analysis, quality_score)
        }
        
        # Save detailed results
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save JSON report
        json_file = output_path / f"{gene_name}_quality_report.json"
        with open(json_file, 'w') as f:
            # Convert numpy types for JSON serialization
            json_results = convert_numpy_types(results)
            json.dump(json_results, f, indent=2)
        
        # Save summary TSV
        summary = create_summary_row(results)
        
        logging.info(f"Quality assessment completed for {gene_name}")
        logging.info(f"Quality score: {quality_score:.2f}/10")
        
        return summary
        
    except Exception as e:
        logging.error(f"Error assessing {alignment_file}: {e}")
        return None

def calculate_overall_quality(basic_stats, conservation_scores, gap_analysis):
    """Calculate overall quality score (0-10)"""
    
    # Sequence count score (2-8 sequences = good, >8 = excellent)
    seq_score = min(2.0, basic_stats['num_sequences'] / 4.0)
    
    # Gap score (penalize high gap percentage)
    gap_penalty = min(2.0, basic_stats['gap_percentage'] / 25.0)  # 25% gaps = 2 point penalty
    gap_score = 2.0 - gap_penalty
    
    # Conservation score
    mean_conservation = np.mean(conservation_scores)
    conservation_score = mean_conservation * 3.0  # 0-3 points
    
    # Length consistency score
    length_range = basic_stats['max_sequence_length'] - basic_stats['min_sequence_length']
    avg_length = basic_stats['avg_sequence_length']
    length_consistency = 1.0 - min(1.0, length_range / (avg_length * 0.5))  # Penalize >50% length variation
    length_score = length_consistency * 2.0
    
    # Alignment completeness (penalize too many gap columns)
    total_columns = basic_stats['alignment_length']
    gap_columns = gap_analysis['gap_columns']
    completeness = 1.0 - (gap_columns / total_columns)
    completeness_score = completeness * 1.0
    
    overall_score = seq_score + gap_score + conservation_score + length_score + completeness_score
    return min(10.0, overall_score)

def generate_recommendations(basic_stats, conservation_scores, gap_analysis, quality_score):
    """Generate quality improvement recommendations"""
    recommendations = []
    
    if quality_score < 5.0:
        recommendations.append("LOW QUALITY: Consider manual inspection and potential re-alignment")
    
    if basic_stats['num_sequences'] < 3:
        recommendations.append("Too few sequences for reliable alignment")
    
    if basic_stats['gap_percentage'] > 30:
        recommendations.append("High gap percentage - consider trimming or different alignment parameters")
    
    if np.mean(conservation_scores) < 0.4:
        recommendations.append("Low conservation - verify that sequences are homologous")
    
    if gap_analysis['gap_columns'] > basic_stats['alignment_length'] * 0.2:
        recommendations.append("Many gap-only columns - consider trimming alignment")
    
    if basic_stats['max_sequence_length'] / basic_stats['min_sequence_length'] > 2.0:
        recommendations.append("Large sequence length variation - check for partial sequences")
    
    if not recommendations:
        recommendations.append("Alignment quality appears good")
    
    return recommendations

def convert_numpy_types(obj):
    """Convert numpy types to Python native types for JSON serialization"""
    if isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

def create_summary_row(results):
    """Create summary row for TSV output"""
    return {
        'gene': results['gene'],
        'num_sequences': results['basic_stats']['num_sequences'],
        'alignment_length': results['basic_stats']['alignment_length'],
        'gap_percentage': round(results['basic_stats']['gap_percentage'], 2),
        'mean_conservation': round(results['conservation']['mean_conservation'], 3),
        'quality_score': round(results['quality_score'], 2),
        'status': 'GOOD' if results['quality_score'] >= 7 else 'FAIR' if results['quality_score'] >= 5 else 'POOR'
    }

def process_alignment_directory(alignment_dir, output_dir):
    """Process all alignments in a directory"""
    
    alignment_path = Path(alignment_dir)
    output_path = Path(output_dir)
    
    if not alignment_path.exists():
        logging.error(f"Alignment directory not found: {alignment_dir}")
        return
    
    # Find all aligned FASTA files
    alignment_files = list(alignment_path.glob("*_aligned.fasta"))
    
    if not alignment_files:
        logging.warning(f"No aligned FASTA files found in {alignment_dir}")
        return
    
    logging.info(f"Found {len(alignment_files)} alignment files to assess")
    
    # Process each alignment
    summaries = []
    for alignment_file in alignment_files:
        summary = assess_alignment_quality(alignment_file, output_dir)
        if summary:
            summaries.append(summary)
    
    # Create combined summary report
    if summaries:
        summary_df = pd.DataFrame(summaries)
        summary_file = output_path / "alignment_quality_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        
        logging.info(f"Summary report saved to: {summary_file}")
        
        # Print quick statistics
        good_alignments = len(summary_df[summary_df['quality_score'] >= 7])
        fair_alignments = len(summary_df[(summary_df['quality_score'] >= 5) & (summary_df['quality_score'] < 7)])
        poor_alignments = len(summary_df[summary_df['quality_score'] < 5])
        
        logging.info(f"\nQuality Summary:")
        logging.info(f"Good alignments (≥7.0): {good_alignments}")
        logging.info(f"Fair alignments (5.0-6.9): {fair_alignments}")
        logging.info(f"Poor alignments (<5.0): {poor_alignments}")

def main():
    """Main function for Snakemake integration"""
    logging.info("MSA Quality Assessment")
    logging.info("="*50)
    
    try:
        # Get inputs from Snakemake
        alignment_dir = snakemake.input[0]
        output_dir = snakemake.output[0]
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        
        # Process all alignments in the directory
        process_alignment_directory(alignment_dir, output_dir)
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        
        alignment_dir = "results/msa_alignments/analysis_1_params_1_gram_positive"
        output_dir = "results/msa_quality/test"
        
        if Path(alignment_dir).exists():
            process_alignment_directory(alignment_dir, output_dir)
        else:
            logging.error(f"Test alignment directory not found: {alignment_dir}")

if __name__ == "__main__":
    main()