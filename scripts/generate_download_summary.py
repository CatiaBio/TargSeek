#!/usr/bin/env python3
"""
Download Summary Generator
=========================

This script generates a comprehensive summary table showing:
- Expected species count per gene (from coverage analysis)
- Actual downloaded species count per gene
- Real coverage calculations based on downloaded sequences
- Download success rates and missing species analysis
"""

import pandas as pd
import json
from pathlib import Path
import logging
from datetime import datetime
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def count_downloaded_sequences(protein_fasta_dir, gene):
    """Count actual downloaded FASTA files for a gene"""
    gene_dir = Path(protein_fasta_dir) / gene
    
    if not gene_dir.exists():
        return 0, []
    
    fasta_files = list(gene_dir.glob("*.fasta"))
    species_list = []
    
    for fasta_file in fasta_files:
        # Extract species name from filename
        species_name = fasta_file.stem.replace('_', ' ')
        species_list.append(species_name)
    
    return len(fasta_files), species_list

def load_download_summaries(protein_fasta_dir):
    """Load individual gene download summaries if available"""
    summaries = {}
    protein_dir = Path(protein_fasta_dir)
    
    if not protein_dir.exists():
        return summaries
    
    # Look for download summary files in each gene directory
    for gene_dir in protein_dir.iterdir():
        if gene_dir.is_dir():
            gene_name = gene_dir.name
            
            # Check for different types of summary files
            summary_files = [
                gene_dir / "download_summary.json",
                gene_dir / "final_download_summary.json"
            ]
            
            for summary_file in summary_files:
                if summary_file.exists():
                    try:
                        with open(summary_file, 'r') as f:
                            summary_data = json.load(f)
                            summaries[gene_name] = summary_data
                            break
                    except Exception as e:
                        logging.warning(f"Error reading {summary_file}: {e}")
    
    return summaries

def create_download_summary_table(analysis, paramset, group, 
                                 coverage_file, protein_fasta_dir, 
                                 proteins_to_study_file, output_dir):
    """Create comprehensive download summary table"""
    
    logging.info(f"Creating download summary for {analysis}_{paramset}_gram_{group}")
    
    # Load coverage data (expected counts)
    try:
        coverage_df = pd.read_csv(coverage_file, sep='\t')
        logging.info(f"Loaded coverage data: {len(coverage_df)} genes")
    except Exception as e:
        logging.error(f"Error loading coverage file {coverage_file}: {e}")
        return
    
    # Load proteins_to_study data
    try:
        proteins_df = pd.read_csv(proteins_to_study_file, sep='\t')
        selected_genes = set(proteins_df['gene'].unique())
        logging.info(f"Loaded proteins to study: {len(selected_genes)} genes selected")
    except Exception as e:
        logging.warning(f"Error loading proteins_to_study file {proteins_to_study_file}: {e}")
        selected_genes = set()
    
    # Load download summaries
    download_summaries = load_download_summaries(protein_fasta_dir)
    
    # Create comprehensive summary
    summary_data = []
    
    for _, row in coverage_df.iterrows():
        gene = row['gene']
        expected_count = row['count']
        expected_coverage = row['coverage_percentage']
        
        # Count actual downloaded sequences
        actual_count, downloaded_species = count_downloaded_sequences(protein_fasta_dir, gene)
        
        # Calculate real coverage based on total species in the group
        # We need to determine total species count for this group
        total_species_file = Path(f"data/bacdive/{analysis}/gram_{group}.txt")
        total_species = 0
        if total_species_file.exists():
            try:
                with open(total_species_file, 'r') as f:
                    total_species = len([line.strip() for line in f if line.strip()])
            except:
                total_species = expected_count  # Fallback
        
        if total_species == 0:
            total_species = expected_count
        
        actual_coverage = (actual_count / total_species) * 100 if total_species > 0 else 0
        
        # Get download details if available
        download_details = download_summaries.get(gene, {})
        
        # Determine if gene was selected for analysis
        selected_for_analysis = gene in selected_genes
        
        # Calculate download success rate
        download_success_rate = (actual_count / expected_count) * 100 if expected_count > 0 else 0
        
        summary_data.append({
            'gene': gene,
            'selected_for_analysis': selected_for_analysis,
            'expected_species_count': expected_count,
            'actual_species_count': actual_count,
            'download_success_rate': download_success_rate,
            'expected_coverage_percent': expected_coverage,
            'actual_coverage_percent': actual_coverage,
            'coverage_difference': actual_coverage - expected_coverage,
            'total_species_in_group': total_species,
            'uniprot_batch_found': download_details.get('uniprot_batch_found', 0),
            'uniprot_individual_found': download_details.get('uniprot_individual_found', 0),
            'ncbi_found': download_details.get('ncbi_found', 0),
            'already_existed': download_details.get('already_existed', 0),
            'final_missing': download_details.get('final_missing', expected_count - actual_count),
            'downloaded_species_list': '; '.join(downloaded_species[:10]) + ('...' if len(downloaded_species) > 10 else '')
        })
    
    # Create DataFrame and sort by actual coverage (descending)
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('actual_coverage_percent', ascending=False)
    
    # Save comprehensive summary
    output_file = Path(output_dir) / f"{analysis}_{paramset}_gram_{group}_download_summary.tsv"
    summary_df.to_csv(output_file, sep='\t', index=False)
    
    # Create a focused summary for selected genes only
    selected_summary_df = summary_df[summary_df['selected_for_analysis'] == True].copy()
    selected_output_file = Path(output_dir) / f"{analysis}_{paramset}_gram_{group}_selected_genes_download_summary.tsv"
    selected_summary_df.to_csv(selected_output_file, sep='\t', index=False)
    
    # Generate statistics summary
    stats = {
        'analysis': analysis,
        'paramset': paramset,
        'group': group,
        'total_genes_analyzed': int(len(coverage_df)),
        'genes_selected_for_analysis': int(len(selected_summary_df)),
        'total_species_in_group': int(total_species),
        'overall_stats': {
            'mean_expected_coverage': float(summary_df['expected_coverage_percent'].mean()),
            'mean_actual_coverage': float(summary_df['actual_coverage_percent'].mean()),
            'mean_download_success_rate': float(summary_df['download_success_rate'].mean()),
            'total_expected_downloads': int(summary_df['expected_species_count'].sum()),
            'total_actual_downloads': int(summary_df['actual_species_count'].sum()),
        },
        'selected_genes_stats': {
            'mean_expected_coverage': float(selected_summary_df['expected_coverage_percent'].mean()) if not selected_summary_df.empty else 0.0,
            'mean_actual_coverage': float(selected_summary_df['actual_coverage_percent'].mean()) if not selected_summary_df.empty else 0.0,
            'mean_download_success_rate': float(selected_summary_df['download_success_rate'].mean()) if not selected_summary_df.empty else 0.0,
            'total_expected_downloads': int(selected_summary_df['expected_species_count'].sum()) if not selected_summary_df.empty else 0,
            'total_actual_downloads': int(selected_summary_df['actual_species_count'].sum()) if not selected_summary_df.empty else 0,
        }
    }
    
    # Save statistics
    stats_file = Path(output_dir) / f"{analysis}_{paramset}_gram_{group}_download_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Generate human-readable summary report
    report_file = Path(output_dir) / f"{analysis}_{paramset}_gram_{group}_download_report.txt"
    with open(report_file, 'w') as f:
        f.write(f"Download Summary Report - {group.capitalize()} Bacteria\n")
        f.write("=" * 60 + "\n")
        f.write(f"Analysis: {analysis} | Parameter Set: {paramset}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"OVERVIEW:\n")
        f.write(f"Total species in {group} group: {total_species}\n")
        f.write(f"Total genes analyzed: {len(coverage_df)}\n")
        f.write(f"Genes selected for MSA: {len(selected_summary_df)}\n\n")
        
        f.write(f"DOWNLOAD PERFORMANCE (All Genes):\n")
        f.write(f"Expected downloads: {stats['overall_stats']['total_expected_downloads']}\n")
        f.write(f"Actual downloads: {stats['overall_stats']['total_actual_downloads']}\n")
        f.write(f"Download success rate: {stats['overall_stats']['mean_download_success_rate']:.1f}%\n")
        f.write(f"Mean expected coverage: {stats['overall_stats']['mean_expected_coverage']:.1f}%\n")
        f.write(f"Mean actual coverage: {stats['overall_stats']['mean_actual_coverage']:.1f}%\n\n")
        
        if not selected_summary_df.empty:
            f.write(f"SELECTED GENES PERFORMANCE:\n")
            f.write(f"Expected downloads: {stats['selected_genes_stats']['total_expected_downloads']}\n")
            f.write(f"Actual downloads: {stats['selected_genes_stats']['total_actual_downloads']}\n")
            f.write(f"Download success rate: {stats['selected_genes_stats']['mean_download_success_rate']:.1f}%\n")
            f.write(f"Mean expected coverage: {stats['selected_genes_stats']['mean_expected_coverage']:.1f}%\n")
            f.write(f"Mean actual coverage: {stats['selected_genes_stats']['mean_actual_coverage']:.1f}%\n\n")
            
            f.write(f"TOP 10 SELECTED GENES BY ACTUAL COVERAGE:\n")
            f.write("-" * 60 + "\n")
            top_genes = selected_summary_df.head(10)
            for _, row in top_genes.iterrows():
                f.write(f"{row['gene']}: {row['actual_species_count']}/{row['expected_species_count']} ")
                f.write(f"({row['actual_coverage_percent']:.1f}% coverage, ")
                f.write(f"{row['download_success_rate']:.1f}% success)\n")
            
            if len(selected_summary_df) > 0:
                f.write(f"\nGENES WITH LOW DOWNLOAD SUCCESS (<80%):\n")
                f.write("-" * 60 + "\n")
                low_success = selected_summary_df[selected_summary_df['download_success_rate'] < 80]
                if not low_success.empty:
                    for _, row in low_success.iterrows():
                        f.write(f"{row['gene']}: {row['download_success_rate']:.1f}% success ")
                        f.write(f"({row['actual_species_count']}/{row['expected_species_count']} downloaded)\n")
                else:
                    f.write("All selected genes have >80% download success rate!\n")
    
    logging.info(f"Generated download summaries:")
    logging.info(f"  - Full summary: {output_file}")
    logging.info(f"  - Selected genes: {selected_output_file}")
    logging.info(f"  - Statistics: {stats_file}")
    logging.info(f"  - Report: {report_file}")
    
    return summary_df, selected_summary_df, stats

def main():
    """Main function for Snakemake integration"""
    
    # Check if running under Snakemake
    try:
        # Get parameters from Snakemake
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        coverage_file = snakemake.input.coverage
        protein_fasta_dir = snakemake.input.protein_fasta
        proteins_to_study_file = snakemake.input.proteins_to_study
        output_dir = Path(snakemake.output[0])
        
    except NameError:
        # Running standalone - use default values for testing
        analysis = "analysis_1"
        paramset = "params_1"
        group = "positive"
        
        coverage_file = f"results/coverage/{analysis}_{paramset}_gram_{group}_coverage_count.tsv"
        protein_fasta_dir = f"results/protein_fasta/{analysis}_{paramset}_gram_{group}"
        proteins_to_study_file = f"results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv"
        output_dir = Path(f"results/download_summary")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    logging.info(f"Generating download summary for {analysis}_{paramset}_gram_{group}")
    logging.info(f"Coverage file: {coverage_file}")
    logging.info(f"Protein FASTA dir: {protein_fasta_dir}")
    logging.info(f"Proteins to study: {proteins_to_study_file}")
    logging.info(f"Output directory: {output_dir}")
    
    # Generate comprehensive summary
    summary_df, selected_df, stats = create_download_summary_table(
        analysis, paramset, group,
        coverage_file, protein_fasta_dir, proteins_to_study_file, output_dir
    )
    
    logging.info("Download summary generation complete!")
    
    # Print quick summary to console
    if not selected_df.empty:
        print(f"\nQuick Summary - {group.capitalize()} Bacteria:")
        print(f"Selected genes: {len(selected_df)}")
        print(f"Mean actual coverage: {selected_df['actual_coverage_percent'].mean():.1f}%")
        print(f"Mean download success: {selected_df['download_success_rate'].mean():.1f}%")
        print(f"Total downloads: {selected_df['actual_species_count'].sum()}")

if __name__ == "__main__":
    main()