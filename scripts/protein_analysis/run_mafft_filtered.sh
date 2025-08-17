#!/bin/bash
# Run MAFFT alignment only for genes in the filtered surface gene list
#
# Usage: run_mafft_filtered.sh <input_sequences_dir> <surface_gene_list> <output_dir>
#
# Arguments:
#   input_sequences_dir: Directory containing FASTA files for all genes
#   surface_gene_list: Text file with one gene name per line (filtered genes only)
#   output_dir: Directory to save MAFFT alignment results
#

set -euo pipefail

INPUT_DIR="$1"
GENE_LIST="$2"
OUTPUT_DIR="$3"

# Check arguments
if [[ $# -ne 3 ]]; then
    echo "Error: Wrong number of arguments"
    echo "Usage: $0 <input_sequences_dir> <surface_gene_list> <output_dir>"
    exit 1
fi

# Check if input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input sequences directory not found: $INPUT_DIR"
    exit 1
fi

# Check if gene list file exists
if [[ ! -f "$GENE_LIST" ]]; then
    echo "Error: Surface gene list file not found: $GENE_LIST"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Starting MAFFT alignment for filtered surface genes..."
echo "Input directory: $INPUT_DIR"
echo "Gene list: $GENE_LIST"
echo "Output directory: $OUTPUT_DIR"

# Count total genes in list
TOTAL_GENES=$(wc -l < "$GENE_LIST")
echo "Processing $TOTAL_GENES surface genes"

PROCESSED=0
SUCCESSFUL=0
FAILED=0

# Process each gene in the filtered list
while IFS= read -r gene; do
    # Skip empty lines and comments
    [[ -z "$gene" || "$gene" =~ ^#.* ]] && continue
    
    INPUT_FASTA="$INPUT_DIR/${gene}.fasta"
    OUTPUT_FASTA="$OUTPUT_DIR/${gene}_aligned.fasta"
    
    ((PROCESSED++))
    
    if [[ ! -f "$INPUT_FASTA" ]]; then
        echo "Warning: Input FASTA not found for gene $gene: $INPUT_FASTA"
        ((FAILED++))
        continue
    fi
    
    echo "Processing gene $gene ($PROCESSED/$TOTAL_GENES)..."
    
    # Run MAFFT alignment
    if mafft --auto --quiet "$INPUT_FASTA" > "$OUTPUT_FASTA" 2>/dev/null; then
        echo "  ✓ Successfully aligned $gene"
        ((SUCCESSFUL++))
    else
        echo "  ✗ Failed to align $gene"
        ((FAILED++))
        # Remove failed output file
        rm -f "$OUTPUT_FASTA"
    fi
    
done < "$GENE_LIST"

echo ""
echo "MAFFT alignment completed:"
echo "  Total genes processed: $PROCESSED"
echo "  Successful alignments: $SUCCESSFUL"
echo "  Failed alignments: $FAILED"

if [[ $SUCCESSFUL -eq 0 ]]; then
    echo "Error: No successful alignments produced"
    exit 1
fi

echo "Alignment results saved to: $OUTPUT_DIR"