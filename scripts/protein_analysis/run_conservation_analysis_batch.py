#!/usr/bin/env python3
"""
Batch Epitope Conservation Analysis
===================================

This script runs epitope conservation analysis for all genes and creates
a comprehensive HTML report with the same styling as the download report.

Usage:
    python scripts/protein_analysis/run_conservation_analysis_batch.py
"""

import subprocess
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
import glob
import logging

# Configure logging for Snakemake compatibility
if 'snakemake' in globals():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
else:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def find_all_genes(epitope_predictions_dir):
    """Find all genes with epitope predictions"""
    genes = []
    epitope_dir = Path(epitope_predictions_dir)
    
    if epitope_dir.exists():
        for gene_dir in epitope_dir.iterdir():
            if gene_dir.is_dir() and any(gene_dir.glob("*_linear_epitopes.tsv")):
                genes.append(gene_dir.name)
    
    return sorted(genes)

def run_conservation_analysis(gene, gram_type, analysis, paramset, output_dir):
    """Run conservation analysis for a single gene"""
    try:
        cmd = [
            "python", "scripts/protein_analysis/analyze_epitope_conservation.py",
            gene,
            "--gram-type", gram_type,
            "--analysis", analysis,
            "--paramset", paramset,
            "--output-dir", str(output_dir)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logging.info(f"✓ Completed analysis for {gene} ({gram_type})")
            return True
        else:
            logging.warning(f"✗ Failed analysis for {gene} ({gram_type}): {result.stderr}")
            return False
            
    except Exception as e:
        logging.error(f"Error running analysis for {gene}: {e}")
        return False

def load_gene_structure_mapping(pdb_numbering_file):
    """Load gene to structure mapping from PDB numbering file"""
    mapping = {}
    
    if Path(pdb_numbering_file).exists():
        try:
            df = pd.read_csv(pdb_numbering_file, sep='\t')
            for _, row in df.iterrows():
                gene = row['Gene']
                structure = row['Structure_ID']
                if gene not in mapping:
                    mapping[gene] = structure
        except Exception as e:
            logging.warning(f"Error loading PDB mapping: {e}")
    
    return mapping


def main():
    """Main function for both Snakemake and command line usage"""
    
    if 'snakemake' in globals():
        # Running from Snakemake
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        output_dir = Path(snakemake.params.conservation_analysis_dir)
        sentinel_file = snakemake.output.conservation_analysis_sentinel
        
        # Input files
        epitope_tables_sentinel = snakemake.input.epitope_tables_sentinel
        pdb_numbering_mapping = snakemake.input.pdb_numbering_mapping
        
        logging.info(f"Running epitope conservation analysis for {analysis}_{paramset}")
        
        # Derive epitope predictions directory from sentinel file
        epitope_predictions_dir = Path(epitope_tables_sentinel).parent
        
    else:
        # Command line usage (fallback)
        analysis = "analysis1"
        paramset = "params1"
        output_dir = Path("epitope_conservation_analysis")
        sentinel_file = output_dir / "conservation_analysis_complete.sentinel"
        epitope_predictions_dir = Path(f"results/{analysis}_{paramset}/protein_analysis/sequences_with_structure/epitope_predictions_bepipred")
        pdb_numbering_mapping = f"data/protein_structures/{analysis}_{paramset}_pdb_numbering_mapping.tsv"
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all genes with epitope predictions
    genes = find_all_genes(epitope_predictions_dir)
    
    if not genes:
        logging.error(f"No genes with epitope predictions found in {epitope_predictions_dir}")
        return
    
    logging.info(f"Found {len(genes)} genes for conservation analysis")
    
    # Load gene to structure mapping
    gene_structure_mapping = load_gene_structure_mapping(pdb_numbering_mapping)
    
    # Run conservation analysis for each gene
    successful_analyses = 0
    failed_analyses = 0
    
    for gene in genes:
        # Determine gram type based on gene presence in selected paths
        gram_types = []
        
        # Check both gram types
        for gram_type in ['positive', 'negative']:
            selected_paths_file = f"results/{analysis}_{paramset}/protein_analysis/sequences_with_structure/selected_3d_paths_gram_{gram_type}.txt"
            
            if Path(selected_paths_file).exists():
                with open(selected_paths_file, 'r') as f:
                    content = f.read()
                    if f"/{gene}/" in content:
                        gram_types.append(gram_type)
        
        # Run analysis for each applicable gram type
        for gram_type in gram_types:
            success = run_conservation_analysis(gene, gram_type, analysis, paramset, output_dir)
            if success:
                successful_analyses += 1
            else:
                failed_analyses += 1
    
    # Create sentinel file to mark analysis completion
    with open(sentinel_file, 'w') as f:
        f.write(f"Epitope conservation analysis completed for {analysis}_{paramset}\\n")
        f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        f.write(f"Genes processed: {len(genes)}\\n")
        f.write(f"Successful analyses: {successful_analyses}\\n")
        f.write(f"Failed analyses: {failed_analyses}\\n")
        f.write(f"Output directory: {output_dir}\\n")
    
    logging.info(f"Analysis complete!")
    logging.info(f"Results: {successful_analyses} successful, {failed_analyses} failed")
    logging.info(f"JSON results saved to: {output_dir}")
    logging.info(f"Sentinel file: {sentinel_file}")

if __name__ == "__main__":
    main()