#!/bin/bash

# Script to create structure range files showing where 3D structure starts and ends
# Format: structure_id/chain/start-end (e.g., 5D0Q/A/24-806)

create_structure_range_file() {
    local gene_name="$1"
    local structure_file="$2"
    
    echo "Processing: $structure_file"
    
    # Extract structure ID from filename
    local structure_id=$(basename "$structure_file" | sed 's/_numbering\.txt$//')
    
    # Output directory structure
    local base_output_dir="results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/${gene_name}"
    
    # Process each chain separately
    local current_chain=""
    local chain_start=""
    local chain_end=""
    local -a chain_data
    
    # Read the numbering file and group by chain
    while IFS=' ' read -r pdb_chain pdb_num pdb_aa; do
        if [[ "$pdb_chain" != "$current_chain" ]]; then
            # Save previous chain data if exists
            if [[ -n "$current_chain" && -n "$chain_start" && -n "$chain_end" ]]; then
                chain_data+=("${current_chain}:${chain_start}-${chain_end}")
            fi
            
            # Start new chain
            current_chain="$pdb_chain"
            chain_start="$pdb_num"
            chain_end="$pdb_num"
        else
            # Update end position for current chain
            chain_end="$pdb_num"
        fi
    done < "$structure_file"
    
    # Don't forget the last chain
    if [[ -n "$current_chain" && -n "$chain_start" && -n "$chain_end" ]]; then
        chain_data+=("${current_chain}:${chain_start}-${chain_end}")
    fi
    
    # Create range files for each chain
    for chain_info in "${chain_data[@]}"; do
        local chain="${chain_info%%:*}"
        local range="${chain_info##*:}"
        
        # Create directory structure: structure_id/tchain/tstart-end
        local range_dir="${base_output_dir}/${structure_id}/t${chain}"
        local range_filename="t${range}"
        local full_path="${range_dir}/${range_filename}"
        
        # Create directories
        mkdir -p "$range_dir"
        
        # Create the range file with metadata
        cat > "$full_path" << EOF
# 3D Structure Range File
# Gene: $gene_name
# Structure: $structure_id
# Chain: $chain
# Range: $range
# Generated: $(date)
#
# This file indicates that the 3D structure covers
# amino acid positions $range for chain $chain
EOF
        
        echo "Created: $full_path"
        echo "  Structure: $structure_id"
        echo "  Chain: $chain"
        echo "  Range: $range"
    done
    
    echo "---"
}

# Main execution
if [ $# -eq 0 ]; then
    # Process all numbering files
    echo "Creating structure range files for all 3D structures..."
    echo "======================================================"
    
    find data/protein_structures/ -name "*_numbering.txt" -type f | while read -r numbering_file; do
        # Extract gene name from path
        local gene_name=$(echo "$numbering_file" | sed -n 's|data/protein_structures/\([^/]*\)/.*|\1|p')
        
        create_structure_range_file "$gene_name" "$numbering_file"
    done
    
    echo "======================================================"
    echo "Processing complete!"
    echo ""
    echo "Range files created in results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/"
else
    echo "Usage: $0"
    echo "This script processes all *_numbering.txt files automatically"
fi