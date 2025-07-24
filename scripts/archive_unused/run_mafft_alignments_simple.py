#!/usr/bin/env python3
"""
MAFFT Multiple Sequence Alignment Runner - Updated Version
===========================================================

This script runs MAFFT alignments directly on per-gene MSA FASTA files:
- Input: results/{analysis}_{paramset}/msa_sequences/gram_{group}/{gene}.fasta
- Output: results/{analysis}_{paramset}/msa_alignments/gram_{group}/{gene}.fasta
"""

import sys
import logging
from pathlib import Path
import subprocess

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_mafft_alignment(input_fasta, output_fasta, threads=8):
    try:
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        cmd = ["mafft", "--auto", "--thread", str(threads), str(input_fasta)]
        with open(output_fasta, 'w') as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, timeout=300)
        if result.returncode == 0:
            logging.info(f"✓ Aligned {input_fasta.name} -> {output_fasta.name}")
            return True
        else:
            logging.error(f"✗ MAFFT failed for {input_fasta.name}: {result.stderr}")
            return False
    except Exception as e:
        logging.error(f"✗ Exception for {input_fasta.name}: {e}")
        return False

def process_all_alignments(msa_dir, align_dir, threads):
    msa_dir = Path(msa_dir)
    align_dir = Path(align_dir)
    align_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(msa_dir.glob("*.fasta"))
    if not fasta_files:
        logging.warning(f"No FASTA files found in {msa_dir}")
        return

    for fasta in fasta_files:
        output_fasta = align_dir / fasta.name
        run_mafft_alignment(fasta, output_fasta, threads)

def main():
    logging.info("Running MAFFT alignment from flat MSA FASTA files")
    try:
        msa_dir = Path(snakemake.input.msa_dir)
        align_dir = Path(snakemake.output.no_3d)
        threads = snakemake.params.threads
    except NameError:
        msa_dir = Path(sys.argv[1])
        align_dir = Path(sys.argv[2])
        threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4

    process_all_alignments(msa_dir, align_dir, threads)

if __name__ == "__main__":
    main()
