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

# Check if there are any alignment files to process
if ls "$ALIGNMENTS_DIR"/*.fasta 1> /dev/null 2>&1; then
    # Process each original alignment
    for aln_file in "$ALIGNMENTS_DIR"/*.fasta; do
        gene=$(basename "$aln_file" .fasta)
        # Remove _aligned suffix if present
        gene=${gene%_aligned}
        out_file="${OUTPUT_DIR}/${gene}_aligned.fasta"

        if [[ -n "${STRUCT_MAP[$gene]:-}" && -f "${STRUCT_MAP[$gene]}" ]]; then
            echo "Adding 3D structure for $gene..."
            
            # Create temporary file with simplified header for ConSurf compatibility
            temp_struct_file="${OUTPUT_DIR}/${gene}_temp_struct.fasta"
            
            # Extract PDB ID from the structure file path and create simple header
            struct_file="${STRUCT_MAP[$gene]}"
            pdb_id=$(basename "$struct_file" .fasta)
            
            # Remove _1, _2 suffixes but keep names with - (like AF-C5D794)
            pdb_id=${pdb_id%_[0-9]*}
            
            # Create simplified FASTA file with PDB ID as header (e.g., >AF-C5D794)
            if head -n 1 "$struct_file" | grep -q "^>"; then
                # Keep the sequence but replace header with simple PDB ID
                echo ">${pdb_id}" > "$temp_struct_file"
                tail -n +2 "$struct_file" >> "$temp_struct_file"
                
                echo "  Using simplified header: >${pdb_id}"
                if mafft --add "$temp_struct_file" --reorder "$aln_file" > "$out_file" 2>/dev/null; then
                    echo "  ✓ Successfully added 3D structure to $gene"
                else
                    echo "  ✗ Failed to add 3D structure to $gene, copying original alignment"
                    cp "$aln_file" "$out_file"
                fi
                
                # Clean up temporary file
                rm -f "$temp_struct_file"
            else
                echo "  Warning: Invalid FASTA format in $struct_file"
                cp "$aln_file" "$out_file"
            fi
        else
            cp "$aln_file" "$out_file"
        fi
    done
else
    echo "No alignment files found in $ALIGNMENTS_DIR - skipping 3D structure addition"
fi
