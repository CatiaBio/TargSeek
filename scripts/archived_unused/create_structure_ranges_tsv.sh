#!/bin/bash

# Script to create a single TSV file with all 3D structure ranges

output_file="data/protein_structures/structure_ranges.tsv"

echo "Creating structure ranges TSV file..."

# Create header
echo -e "Gene\tStructure_ID\tChain\tStart\tEnd\tRange" > "$output_file"

# Process all numbering files
find data/protein_structures/ -name "*_numbering.txt" -type f | sort | while read -r numbering_file; do
    # Extract gene name from path
    gene_name=$(echo "$numbering_file" | sed -n 's|data/protein_structures/\([^/]*\)/.*|\1|p')
    
    # Extract structure ID from filename
    structure_id=$(basename "$numbering_file" | sed 's/_numbering\.txt$//')
    
    echo "Processing: $gene_name - $structure_id"
    
    # Process each chain separately
    current_chain=""
    chain_start=""
    chain_end=""
    
    # Read the numbering file and group by chain
    while IFS=' ' read -r pdb_chain pdb_num pdb_aa; do
        # Skip non-numeric residue numbers (coordinates, etc.)
        if ! [[ "$pdb_num" =~ ^[0-9]+$ ]]; then
            continue
        fi
        
        if [[ "$pdb_chain" != "$current_chain" ]]; then
            # Save previous chain data if exists
            if [[ -n "$current_chain" && -n "$chain_start" && -n "$chain_end" ]]; then
                echo -e "${gene_name}\t${structure_id}\t${current_chain}\t${chain_start}\t${chain_end}\t${chain_start}-${chain_end}" >> "$output_file"
            fi
            
            # Start new chain
            current_chain="$pdb_chain"
            chain_start="$pdb_num"
            chain_end="$pdb_num"
        else
            # Update end position for current chain (keep track of min/max)
            if (( pdb_num < chain_start )); then
                chain_start="$pdb_num"
            fi
            if (( pdb_num > chain_end )); then
                chain_end="$pdb_num"
            fi
        fi
    done < "$numbering_file"
    
    # Don't forget the last chain
    if [[ -n "$current_chain" && -n "$chain_start" && -n "$chain_end" ]]; then
        echo -e "${gene_name}\t${structure_id}\t${current_chain}\t${chain_start}\t${chain_end}\t${chain_start}-${chain_end}" >> "$output_file"
    fi
done

echo "Structure ranges TSV created: $output_file"
echo ""
echo "Summary:"
wc -l "$output_file"
echo ""
echo "First 10 lines:"
head -10 "$output_file"