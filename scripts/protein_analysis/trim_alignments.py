#!/usr/bin/env python3
"""
MSA Alignment Trimming (ClipKIT)
================================

This script runs ClipKIT on all alignment files to remove 
poorly aligned regions, gaps, and spurious sequences, improving alignment quality 
for downstream analysis.

ClipKIT retains phylogenetically informative sites while removing spurious gaps.
"""

import subprocess
from pathlib import Path
import logging
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_clipkit(input_file, output_file, mode="smart-gap"):
    """
    Run ClipKIT on a single alignment file (preferred method)
    
    Args:
        input_file: Path to input aligned FASTA file
        output_file: Path to output trimmed FASTA file
        mode: ClipKIT mode (smart-gap, kpic, kpic-smart-gap, kpi, etc.)
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Create output directory if needed
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Run ClipKIT command
        cmd = [
            "clipkit",
            str(input_file),
            "-o", str(output_file),
            "-m", mode
        ]
        
        logging.info(f"Running ClipKIT: {input_file.name} -> {output_file.name}")
        logging.debug(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        logging.debug(f"ClipKIT stdout: {result.stdout}")
        return True
        
    except subprocess.CalledProcessError as e:
        logging.error(f"ClipKIT failed for {input_file.name}: {e}")
        logging.error(f"ClipKIT stderr: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error running ClipKIT on {input_file.name}: {e}")
        return False


def process_alignment_directory(input_dir, output_dir, dir_name):
    """Process all alignments in a directory"""
    
    logging.info(f"\nProcessing {dir_name} alignments:")
    logging.info(f"  Input: {input_dir}")
    logging.info(f"  Output: {output_dir}")
    
    if not input_dir.exists():
        logging.warning(f"  Directory does not exist: {input_dir}")
        return 0, 0
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all FASTA files (alignments)
    fasta_files = list(input_dir.glob("*.fasta"))
    
    if not fasta_files:
        logging.warning(f"  No FASTA files found in {input_dir}")
        return 0, 0
    
    logging.info(f"  Found {len(fasta_files)} alignment files to trim")
    
    # Process each alignment file
    successful = 0
    failed = 0
    
    for fasta_file in fasta_files:
        # Create output filename
        gene_name = fasta_file.stem
        output_file = output_dir / f"{gene_name}.fasta"
        
        # Run ClipKIT on alignment
        clipkit_mode = getattr(snakemake.params, 'clipkit_mode', 'smart-gap') if 'snakemake' in globals() else 'smart-gap'
        if run_clipkit(fasta_file, output_file, mode=clipkit_mode):
            successful += 1
        else:
            failed += 1
    
    return successful, failed

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        input_main = Path(snakemake.input.main)
        input_structure = Path(snakemake.input.structure)
        output_main = Path(snakemake.output.trim_main)
        output_structure = Path(snakemake.output.trim_structure)
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
    except NameError:
        # Test mode
        if len(sys.argv) != 5:
            print("Usage: python trim_alignments.py <input_main> <input_structure> <output_main> <output_structure>")
            sys.exit(1)
            
        input_main = Path(sys.argv[1])
        input_structure = Path(sys.argv[2])
        output_main = Path(sys.argv[3])
        output_structure = Path(sys.argv[4])
        
        analysis = "test"
        paramset = "test"
        group = "test"
    
    logging.info(f"ClipKIT Alignment Trimming")
    logging.info(f"="*50)
    logging.info(f"Analysis: {analysis}")
    logging.info(f"Paramset: {paramset}")
    logging.info(f"Group: {group}")
    
    total_successful = 0
    total_failed = 0
    
    # Process main alignments
    successful, failed = process_alignment_directory(input_main, output_main, "main")
    total_successful += successful
    total_failed += failed
    
    # Process structure alignments
    successful, failed = process_alignment_directory(input_structure, output_structure, "structure")
    total_successful += successful
    total_failed += failed
    
    # Summary
    logging.info(f"\nClipKIT processing complete:")
    logging.info(f"  ✓ Total successful: {total_successful}")
    logging.info(f"  ✗ Total failed: {total_failed}")
    
    if total_failed > 0:
        logging.warning(f"Some files failed trimming. Check logs for details.")

if __name__ == "__main__":
    main()