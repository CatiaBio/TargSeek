#!/usr/bin/env python3
"""
MAFFT Multiple Sequence Alignment Runner
=======================================

This script runs MAFFT alignments on all genes in a group.
It reads input FASTA files from:
  - main_msa_sequences/gram_{group}      → main MSA sequences only
  - structure_msa_sequences/gram_{group} → main MSA + 3D structures
And writes aligned outputs to:
  - main_alignments/gram_{group}      → alignments from main MSA
  - structure_alignments/gram_{group} → alignments from structure MSA
"""

import sys
import logging
from pathlib import Path
import subprocess
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_mafft(input_fasta: Path, output_fasta: Path, threads: int = 8):
    try:
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            "mafft", "--thread", str(threads), "--retree", "2", "--maxiterate", "0", str(input_fasta)
        ]
        with open(output_fasta, 'w') as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=True)
        logging.info(f"✓ Aligned: {input_fasta.name} → {output_fasta.name}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"✗ MAFFT failed for {input_fasta.name}: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"✗ Error for {input_fasta.name}: {e}")
        return False

def main():
    try:
        msa_main = Path(snakemake.input.msa_main)
        msa_structure = Path(snakemake.input.msa_structure)
        output_main = Path(snakemake.output.main)
        output_structure = Path(snakemake.output.structure)
        threads = snakemake.params.threads
        enable_no_3d = snakemake.params.enable_no_3d
        enable_with_3d = snakemake.params.enable_with_3d
    except NameError:
        msa_main = Path(sys.argv[1])
        msa_structure = Path(sys.argv[2])
        output_main = Path(sys.argv[3])
        output_structure = Path(sys.argv[4])
        threads = int(sys.argv[5]) if len(sys.argv) > 5 else 4
        enable_no_3d = True
        enable_with_3d = True

    output_main.mkdir(parents=True, exist_ok=True)
    output_structure.mkdir(parents=True, exist_ok=True)

    logging.info("=== MAFFT Multiple Sequence Alignment ===")
    logging.info(f"Main MSA input: {msa_main}")
    logging.info(f"Structure MSA input: {msa_structure}")
    logging.info(f"Main alignments output: {output_main}")
    logging.info(f"Structure alignments output: {output_structure}")
    logging.info(f"Threads: {threads}")

    total_genes = 0
    successful_main = 0
    successful_structure = 0

    # Get all FASTA files from structure MSA (this should be the complete set)
    fasta_files = list(msa_structure.glob("*.fasta"))
    total_genes = len(fasta_files)
    
    logging.info(f"Found {total_genes} genes to align")

    # Align main sequences (if enabled)
    if enable_no_3d and msa_main.exists():
        logging.info("\n--- Aligning Main MSA Sequences ---")
        for fasta_path in sorted(msa_main.glob("*.fasta")):
            gene = fasta_path.stem
            output_fasta = output_main / f"{gene}.fasta"
            if run_mafft(fasta_path, output_fasta, threads):
                successful_main += 1
    else:
        logging.info("--- Skipping Main MSA Alignment (disabled or missing) ---")

    # Align structure sequences (if enabled)
    if enable_with_3d and msa_structure.exists():
        logging.info("\n--- Aligning Structure MSA Sequences ---")
        for fasta_path in sorted(msa_structure.glob("*.fasta")):
            gene = fasta_path.stem
            output_fasta = output_structure / f"{gene}.fasta"
            if run_mafft(fasta_path, output_fasta, threads):
                successful_structure += 1
    else:
        logging.info("--- Skipping Structure MSA Alignment (disabled or missing) ---")

    # Summary
    logging.info(f"\n=== MAFFT Alignment Summary ===")
    logging.info(f"Total genes found: {total_genes}")
    if enable_no_3d:
        logging.info(f"Main alignments successful: {successful_main}")
    if enable_with_3d:
        logging.info(f"Structure alignments successful: {successful_structure}")
    logging.info(f"=== MAFFT alignment completed! ===")

if __name__ == "__main__":
    main()