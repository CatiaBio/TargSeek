#!/usr/bin/env python3
"""
MAFFT Multiple Sequence Alignment Runner
=======================================

This script runs MAFFT alignments on all genes in a group.
It reads input FASTA files from:
  - msa_sequences/gram_{group}         → no 3D structure
  - msa_sequences_with_3d_fasta/gram_{group} → with 3D structure
And writes aligned outputs to:
  - msa_alignments/gram_{group}
  - msa_alignments_with_3d_fasta/gram_{group}
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
        msa_no_3d = Path(snakemake.input.msa_no_3d)
        msa_with_3d = Path(snakemake.input.msa_with_3d)
        output_no_3d = Path(snakemake.output.no_3d)
        output_with_3d = Path(snakemake.output.with_3d)
        threads = snakemake.params.threads
    except NameError:
        msa_no_3d = Path(sys.argv[1])
        msa_with_3d = Path(sys.argv[2])
        output_no_3d = Path(sys.argv[3])
        output_with_3d = Path(sys.argv[4])
        threads = int(sys.argv[5]) if len(sys.argv) > 5 else 4

    output_no_3d.mkdir(parents=True, exist_ok=True)
    output_with_3d.mkdir(parents=True, exist_ok=True)

    # Align no-3D sequences
    for fasta_path in sorted(msa_no_3d.glob("*.fasta")):
        gene = fasta_path.stem
        output_fasta = output_no_3d / f"{gene}.fasta"
        run_mafft(fasta_path, output_fasta, threads)

    # Align with-3D sequences
    for fasta_path in sorted(msa_with_3d.glob("*.fasta")):
        gene = fasta_path.stem
        output_fasta = output_with_3d / f"{gene}.fasta"
        run_mafft(fasta_path, output_fasta, threads)

if __name__ == "__main__":
    main()