#!/usr/bin/env python3
"""
MAFFT Multiple Sequence Alignment Runner
=======================================

This script runs MAFFT alignments on all genes in a group, reading from
MSA-ready sequences and producing aligned FASTA files.
"""

import pandas as pd
from pathlib import Path
import subprocess
import os
import sys
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_mafft_alignment(input_fasta, output_fasta, threads=8):
    """
    Run MAFFT alignment on a single gene
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output aligned FASTA file
        threads: Number of threads for MAFFT
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Create output directory if needed
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        
        # Run MAFFT command - fast but better than auto
        cmd = [
            "mafft",
            "--thread", str(threads),
            "--retree", "2",  # Fast tree estimation, better than --auto
            "--maxiterate", "0",  # No refinement for speed
            str(input_fasta)
        ]
        
        logging.info(f"Running MAFFT command: {' '.join(cmd)}")
        
        with open(output_fasta, 'w') as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        # Check if output file was created and has content
        if output_fasta.exists() and output_fasta.stat().st_size > 0:
            logging.info(f"✓ Successfully aligned: {output_fasta}")
            return True
        else:
            logging.error(f"✗ Output file is empty or missing: {output_fasta}")
            return False
            
    except subprocess.CalledProcessError as e:
        logging.error(f"✗ MAFFT failed for {input_fasta}: {e}")
        if e.stderr:
            logging.error(f"MAFFT stderr: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"✗ Unexpected error for {input_fasta}: {e}")
        return False

def process_gene_alignments(msa_dir, protein_list_file, output_dir, threads=8):
    """
    Process all genes in a group for MAFFT alignment
    
    Args:
        msa_dir: Directory containing MSA-ready sequences
        protein_list_file: TSV file with list of proteins to process
        output_dir: Output directory for aligned sequences
        threads: Number of threads for MAFFT
    """
    
    logging.info(f"Processing MAFFT alignments")
    logging.info(f"Input MSA directory: {msa_dir}")
    logging.info(f"Protein list: {protein_list_file}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"MAFFT threads: {threads}")
    
    # Read the protein list to get genes
    try:
        proteins_df = pd.read_csv(protein_list_file, sep='\t')
        genes = proteins_df['gene'].unique().tolist()
        logging.info(f"Found {len(genes)} genes to process")
    except Exception as e:
        logging.error(f"Error reading protein list: {e}")
        return
    
    # Convert paths to Path objects
    msa_path = Path(msa_dir)
    output_path = Path(output_dir)
    
    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Track statistics
    successful_alignments = 0
    failed_alignments = 0
    skipped_genes = 0
    
    # Process each gene
    for i, gene in enumerate(genes, 1):
        logging.info(f"[{i}/{len(genes)}] Processing gene: {gene}")
        
        input_fasta = msa_path / f"{gene}.fasta"
        output_fasta = output_path / f"{gene}_aligned.fasta"
        
        if not input_fasta.exists():
            logging.warning(f"  ⚠ Skipping {gene}: No input FASTA found at {input_fasta}")
            skipped_genes += 1
            continue
        
        # Check if sequences exist in the file
        try:
            with open(input_fasta, 'r') as f:
                content = f.read().strip()
                if not content or content.count('>') < 2:
                    logging.warning(f"  ⚠ Skipping {gene}: Less than 2 sequences in {input_fasta}")
                    skipped_genes += 1
                    continue
        except Exception as e:
            logging.error(f"  ✗ Error reading {input_fasta}: {e}")
            failed_alignments += 1
            continue
        
        # Run MAFFT alignment
        if run_mafft_alignment(input_fasta, output_fasta, threads):
            successful_alignments += 1
        else:
            failed_alignments += 1
    
    # Summary statistics
    logging.info(f"\n=== MAFFT Alignment Summary ===")
    logging.info(f"Total genes requested: {len(genes)}")
    logging.info(f"Successful alignments: {successful_alignments}")
    logging.info(f"Failed alignments: {failed_alignments}")
    logging.info(f"Skipped genes: {skipped_genes}")
    logging.info(f"Success rate: {successful_alignments/(len(genes)-skipped_genes)*100:.1f}%" if len(genes) > skipped_genes else "No alignments attempted")

def main():
    """Main function for Snakemake integration"""
    logging.info("MAFFT Alignment Runner")
    logging.info("="*50)
    
    try:
        # Get inputs from Snakemake
        msa_dir = snakemake.input.msa_dir
        protein_list_file = snakemake.input.protein_list
        output_dir = snakemake.output[0]
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        threads = snakemake.params.get('threads', 8)
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        
        # Run the alignment process
        process_gene_alignments(msa_dir, protein_list_file, output_dir, threads)
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        
        # Example test parameters
        msa_dir = "results/msa_sequences/analysis_1_params_1_gram_positive"
        protein_list_file = "results/proteins_to_study/analysis_1_params_1_gram_positive.tsv"
        output_dir = "results/msa_alignments/test"
        threads = 4
        
        if Path(msa_dir).exists() and Path(protein_list_file).exists():
            process_gene_alignments(msa_dir, protein_list_file, output_dir, threads)
        else:
            logging.error("Test files not found")
            logging.info(f"Expected MSA directory: {msa_dir}")
            logging.info(f"Expected protein list: {protein_list_file}")

if __name__ == "__main__":
    main()