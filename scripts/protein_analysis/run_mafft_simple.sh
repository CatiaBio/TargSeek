#!/bin/bash
# Simplified MAFFT script that works with WSL2
set -euo pipefail

INPUT_DIR="$1"
GENE_LIST="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

echo "Starting MAFFT alignments..."
echo "Input: $INPUT_DIR"
echo "Genes: $GENE_LIST"
echo "Output: $OUTPUT_DIR"

TOTAL_GENES=$(grep -c -v '^[[:space:]]*$' "$GENE_LIST" 2>/dev/null || echo "0")
echo "Processing $TOTAL_GENES genes"

PROCESSED=0
SUCCESSFUL=0

while IFS= read -r gene; do
    [[ -z "$gene" || "$gene" =~ ^#.* ]] && continue

    INPUT_FASTA="$INPUT_DIR/${gene}.fasta"
    OUTPUT_FASTA="$OUTPUT_DIR/${gene}_aligned.fasta"

    ((PROCESSED++))
    echo "[$PROCESSED/$TOTAL_GENES] Processing $gene..."

    if [[ ! -f "$INPUT_FASTA" ]]; then
        echo "  ⚠ Input not found: $INPUT_FASTA"
        continue
    fi

    # Run MAFFT with WSL2-compatible redirection
    if timeout 300 mafft --auto "$INPUT_FASTA" > "$OUTPUT_FASTA" 2>&1; then
        echo "  ✓ Completed $gene"
        ((SUCCESSFUL++))
    else
        echo "  ✗ Failed $gene"
        rm -f "$OUTPUT_FASTA"
    fi

done < "$GENE_LIST"

echo ""
echo "Summary: $SUCCESSFUL/$PROCESSED genes aligned successfully"
[[ $SUCCESSFUL -gt 0 ]] && echo "Results in: $OUTPUT_DIR"