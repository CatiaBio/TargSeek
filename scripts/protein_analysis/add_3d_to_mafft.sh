#!/usr/bin/env bash
set -euo pipefail

ALIGNMENTS_DIR="$1"
SELECTED_3D_PATHS="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

# Load all structure paths into an associative array keyed by gene
declare -A STRUCT_MAP
if [[ -f "$SELECTED_3D_PATHS" ]]; then
    while read -r struct_path; do
        [[ -z "$struct_path" ]] && continue
        gene=$(basename "$(dirname "$struct_path")")
        STRUCT_MAP[$gene]="$struct_path"
    done < "$SELECTED_3D_PATHS"
fi

# Process each original alignment
for aln_file in "$ALIGNMENTS_DIR"/*.fasta; do
    gene=$(basename "$aln_file" .fasta)
    out_file="${OUTPUT_DIR}/${gene}.fasta"

    if [[ -n "${STRUCT_MAP[$gene]:-}" && -f "${STRUCT_MAP[$gene]}" ]]; then
        echo "Adding 3D structure for $gene..."
        mafft --add "${STRUCT_MAP[$gene]}" --reorder "$aln_file" > "$out_file"
    else
        cp "$aln_file" "$out_file"
    fi
done
