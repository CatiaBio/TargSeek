#!/bin/bash
# Optimized MAFFT script for WSL2 with adaptive timeouts and speed optimization
set -euo pipefail

INPUT_DIR="$1"
GENE_LIST="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

echo "Starting optimized MAFFT alignments..."
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

    # Count sequences to determine optimal parameters
    SEQ_COUNT=$(grep -c "^>" "$INPUT_FASTA" 2>/dev/null || echo "0")
    echo "  📊 Sequences: $SEQ_COUNT"

    # Choose timeout and parameters based on sequence count
    if [[ $SEQ_COUNT -gt 100 ]]; then
        TIMEOUT=1200  # 20 minutes for large alignments
        MAFFT_ARGS="--maxiterate 2 --retree 1"  # Faster for large sets
        echo "  🚀 Using fast mode for large dataset"
    elif [[ $SEQ_COUNT -gt 50 ]]; then
        TIMEOUT=600   # 10 minutes for medium alignments
        MAFFT_ARGS="--auto"
        echo "  ⚡ Using standard mode"
    else
        TIMEOUT=300   # 5 minutes for small alignments
        MAFFT_ARGS="--auto"
        echo "  🔧 Using standard mode"
    fi

    # Run MAFFT with adaptive timeout and WSL2-compatible redirection
    echo "  ⏱️ Timeout: ${TIMEOUT}s, Starting alignment..."
    if timeout $TIMEOUT mafft $MAFFT_ARGS "$INPUT_FASTA" > "$OUTPUT_FASTA" 2>&1; then
        echo "  ✅ Completed $gene ($(wc -l < "$OUTPUT_FASTA") lines)"
        ((SUCCESSFUL++))
    else
        echo "  ❌ Failed/timeout $gene"
        rm -f "$OUTPUT_FASTA"
    fi

done < "$GENE_LIST"

echo ""
echo "🎯 Summary: $SUCCESSFUL/$PROCESSED genes aligned successfully"
if [[ $SUCCESSFUL -gt 0 ]]; then
    echo "📁 Results in: $OUTPUT_DIR"
    echo "🔍 Completed files:"
    ls -1 "$OUTPUT_DIR"/*.fasta 2>/dev/null | head -5
    [[ $(ls -1 "$OUTPUT_DIR"/*.fasta 2>/dev/null | wc -l) -gt 5 ]] && echo "   ... and more"
fi

# Exit with error if no alignments succeeded
[[ $SUCCESSFUL -eq 0 ]] && exit 1 || exit 0