#!/usr/bin/env bash
set -euo pipefail

SEQUENCES_DIR="$1"
ALIGNMENTS_DIR="$2"

mkdir -p "$ALIGNMENTS_DIR"

for fasta_file in "$SEQUENCES_DIR"/*.fasta; do
    if [[ -f "$fasta_file" ]]; then
        gene_name=$(basename "$fasta_file" .fasta)
        echo "Running MAFFT on $gene_name..."
        mafft --localpair --maxiterate 1000 "$fasta_file" > "${ALIGNMENTS_DIR}/${gene_name}.fasta"
    fi
done
