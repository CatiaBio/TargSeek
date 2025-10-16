#!/usr/bin/env python3
"""
MAFFT alignment runner with WSL2 compatibility
Handles large sequence sets with adaptive timeouts
"""

import os
import sys
import subprocess
import time
from pathlib import Path

def count_sequences(fasta_file):
    """Count sequences in FASTA file"""
    try:
        with open(fasta_file, 'r') as f:
            return sum(1 for line in f if line.startswith('>'))
    except:
        return 0

def run_mafft(input_fasta, output_fasta, seq_count):
    """Run MAFFT with adaptive parameters"""

    # Choose parameters based on sequence count
    if seq_count > 100:
        timeout = 1200  # 20 minutes
        args = ["mafft", "--maxiterate", "2", "--retree", "1", input_fasta]
        mode = "fast mode (large dataset)"
    elif seq_count > 50:
        timeout = 600   # 10 minutes
        args = ["mafft", "--auto", input_fasta]
        mode = "standard mode (medium dataset)"
    else:
        timeout = 300   # 5 minutes
        args = ["mafft", "--auto", input_fasta]
        mode = "standard mode (small dataset)"

    print(f"    🚀 Using {mode}, timeout: {timeout}s")

    try:
        # Run MAFFT with timeout - separate stdout (alignment) from stderr (diagnostics)
        with open(output_fasta, 'w') as out_file:
            result = subprocess.run(
                args,
                stdout=out_file,
                stderr=subprocess.DEVNULL,  # Discard diagnostic output
                timeout=timeout,
                check=True
            )
        return True

    except subprocess.TimeoutExpired:
        print(f"    ⏰ Timeout after {timeout}s")
        if os.path.exists(output_fasta):
            os.remove(output_fasta)
        return False

    except subprocess.CalledProcessError as e:
        print(f"    ❌ MAFFT failed with exit code {e.returncode}")
        if os.path.exists(output_fasta):
            os.remove(output_fasta)
        return False

    except Exception as e:
        print(f"    💥 Error: {e}")
        if os.path.exists(output_fasta):
            os.remove(output_fasta)
        return False

def main():
    if len(sys.argv) != 4:
        print("Usage: run_mafft_python.py <input_dir> <gene_list> <output_dir>")
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    gene_list_file = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    print("🧬 Starting MAFFT alignments with Python runner...")
    print(f"📁 Input: {input_dir}")
    print(f"📋 Genes: {gene_list_file}")
    print(f"📁 Output: {output_dir}")

    # Read gene list
    try:
        with open(gene_list_file, 'r') as f:
            genes = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    except Exception as e:
        print(f"❌ Error reading gene list: {e}")
        sys.exit(1)

    total_genes = len(genes)
    print(f"🎯 Processing {total_genes} genes")
    print()

    successful = 0

    for i, gene in enumerate(genes, 1):
        print(f"[{i}/{total_genes}] Processing {gene}...")

        input_fasta = input_dir / f"{gene}.fasta"
        output_fasta = output_dir / f"{gene}_aligned.fasta"

        if not input_fasta.exists():
            print(f"    ⚠️  Input not found: {input_fasta}")
            continue

        # Count sequences
        seq_count = count_sequences(input_fasta)
        print(f"    📊 Sequences: {seq_count}")

        # Run MAFFT
        start_time = time.time()
        if run_mafft(str(input_fasta), str(output_fasta), seq_count):
            elapsed = time.time() - start_time
            lines = sum(1 for _ in open(output_fasta))
            print(f"    ✅ Completed in {elapsed:.1f}s ({lines} lines)")
            successful += 1
        else:
            print(f"    ❌ Failed")

        print()

    print(f"🎯 Summary: {successful}/{total_genes} genes aligned successfully")

    if successful > 0:
        print(f"📁 Results in: {output_dir}")
        aligned_files = list(output_dir.glob("*_aligned.fasta"))
        print(f"🔍 Created {len(aligned_files)} alignment files")

        # Show first few files
        for f in aligned_files[:5]:
            print(f"    • {f.name}")
        if len(aligned_files) > 5:
            print(f"    • ... and {len(aligned_files) - 5} more")

    # Exit code
    sys.exit(0 if successful > 0 else 1)

if __name__ == "__main__":
    main()