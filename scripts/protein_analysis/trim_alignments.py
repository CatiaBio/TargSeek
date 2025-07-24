#!/usr/bin/env python3
"""
trimAl Alignment Trimming
========================

This script runs trimAl on all alignment files to remove poorly aligned regions,
gaps, and spurious sequences, improving alignment quality for downstream analysis.
"""

import subprocess
from pathlib import Path
import logging
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_trimal(input_file, output_file, method="automated1"):
    """
    Run trimAl on a single alignment file
    
    Args:
        input_file: Path to input aligned FASTA file
        output_file: Path to output trimmed FASTA file
        method: trimAl method to use (automated1, strict, etc.)
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Create output directory if needed
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Run trimAl command
        cmd = [
            "trimal",
            "-in", str(input_file),
            "-out", str(output_file),
            f"-{method}"  # e.g., -automated1
        ]
        
        logging.info(f"Running trimAl: {input_file.name} -> {output_file.name}")
        logging.debug(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Check if output file was created and has content
        if output_file.exists() and output_file.stat().st_size > 0:
            logging.info(f"  ✓ Successfully trimmed {input_file.name}")
            return True
        else:
            logging.warning(f"  ✗ trimAl produced empty output for {input_file.name}")
            return False
            
    except subprocess.CalledProcessError as e:
        logging.error(f"  ✗ trimAl failed for {input_file.name}: {e}")
        logging.error(f"  stderr: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"  ✗ Error processing {input_file.name}: {e}")
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
        
        # Run trimAl
        if run_trimal(fasta_file, output_file, method="automated1"):
            successful += 1
        else:
            failed += 1
    
    return successful, failed

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        input_no_3d = Path(snakemake.input.no_3d)
        input_with_3d = Path(snakemake.input.with_3d)
        output_no_3d = Path(snakemake.output.trim_no_3d)
        output_with_3d = Path(snakemake.output.trim_with_3d)
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
    except NameError:
        # Test mode
        if len(sys.argv) != 5:
            print("Usage: python trim_alignments.py <input_no_3d> <input_with_3d> <output_no_3d> <output_with_3d>")
            sys.exit(1)
            
        input_no_3d = Path(sys.argv[1])
        input_with_3d = Path(sys.argv[2])
        output_no_3d = Path(sys.argv[3])
        output_with_3d = Path(sys.argv[4])
        
        analysis = "test"
        paramset = "test"
        group = "test"
    
    logging.info(f"trimAl Alignment Trimming")
    logging.info(f"="*50)
    logging.info(f"Analysis: {analysis}")
    logging.info(f"Paramset: {paramset}")
    logging.info(f"Group: {group}")
    
    total_successful = 0
    total_failed = 0
    
    # Process no-3D alignments
    successful, failed = process_alignment_directory(input_no_3d, output_no_3d, "no-3D")
    total_successful += successful
    total_failed += failed
    
    # Process with-3D alignments
    successful, failed = process_alignment_directory(input_with_3d, output_with_3d, "with-3D")
    total_successful += successful
    total_failed += failed
    
    # Summary
    logging.info(f"\ntrimAl processing complete:")
    logging.info(f"  ✓ Total successful: {total_successful}")
    logging.info(f"  ✗ Total failed: {total_failed}")
    
    if total_failed > 0:
        logging.warning(f"Some files failed trimming. Check logs for details.")

if __name__ == "__main__":
    main()