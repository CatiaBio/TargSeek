#!/bin/bash

# Script to extract PDB amino acid numbering from structure files
# Usage: ./extract_pdb_numbering.sh [pdb_file] or ./extract_pdb_numbering.sh (to process all structures)

extract_pdb_numbering() {
    local pdb_file="$1"
    local base_name
    local output_file
    
    # Handle both .pdb and .pdb.gz extensions
    if [[ "$pdb_file" == *.pdb.gz ]]; then
        base_name="${pdb_file%.pdb.gz}"
        output_file="${base_name}_numbering.txt"
    else
        base_name="${pdb_file%.pdb}"
        output_file="${base_name}_numbering.txt"
    fi
    
    echo "Processing: $pdb_file"
    echo "Output: $output_file"
    
    # Extract chain, residue number, and residue type
    # Handle both compressed and uncompressed files
    if [[ "$pdb_file" == *.gz ]]; then
        zgrep "^ATOM" "$pdb_file" | awk '{print $5, $6, $4}' | uniq > "$output_file"
    else
        grep "^ATOM" "$pdb_file" | awk '{print $5, $6, $4}' | uniq > "$output_file"
    fi
    
    echo "Found $(wc -l < "$output_file") unique residues"
    echo "First 5 residues:"
    head -5 "$output_file"
    echo "---"
}

# If a specific file is provided, process only that file
if [ $# -eq 1 ]; then
    if [ -f "$1" ]; then
        extract_pdb_numbering "$1"
    else
        echo "Error: File $1 not found"
        exit 1
    fi
else
    # Process all PDB files in data/protein_structures/
    echo "Processing all selected 3D structures..."
    echo "========================================"
    
    # Find all PDB files (both .pdb and .pdb.gz) in the protein structures directory
    find data/protein_structures/ \( -name "*.pdb" -o -name "*.pdb.gz" \) -type f | while read -r pdb_file; do
        extract_pdb_numbering "$pdb_file"
    done
    
    echo "========================================"
    echo "Processing complete!"
    echo ""
    echo "Summary of numbering files created:"
    find data/protein_structures/ -name "*_numbering.txt" -type f | wc -l
    echo "numbering files created"
fi