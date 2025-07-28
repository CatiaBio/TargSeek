#!/bin/bash

# Script to map FASTA sequences used for BepiPred to PDB numbering
# Usage: ./map_fasta_to_pdb_numbering.sh [fasta_file] [structure_id] or ./map_fasta_to_pdb_numbering.sh (to process all)

map_fasta_to_pdb() {
    local bepipred_fasta="$1"
    local structure_id="$2"
    # Determine gene name from the fasta file path
    local gene_name=$(echo "$bepipred_fasta" | sed -n 's|.*/epitope_predictions_bepipred/\([^/]*\)/.*|\1|p')
    local numbering_file="data/protein_structures/${gene_name}/${structure_id}_numbering.txt"
    local output_file="${bepipred_fasta%.fasta}_pdb_mapping.txt"
    
    echo "Processing: $bepipred_fasta"
    echo "Structure ID: $structure_id"
    echo "Numbering file: $numbering_file"
    echo "Output: $output_file"
    
    # Check if numbering file exists
    if [[ ! -f "$numbering_file" ]]; then
        echo "ERROR: Numbering file not found: $numbering_file"
        echo "Run ./extract_pdb_numbering.sh first to generate numbering files"
        echo "---"
        return 1
    fi
    
    # Check if FASTA file exists
    if [[ ! -f "$bepipred_fasta" ]]; then
        echo "ERROR: BepiPred FASTA file not found: $bepipred_fasta"
        echo "---"
        return 1
    fi
    
    # Extract header and determine chain
    local header=$(head -1 "$bepipred_fasta")
    local chain="A"  # Default to chain A
    
    # Try to extract chain from header (look for chain information)
    if [[ "$header" =~ [Cc]hain[s]?[[:space:]]*([A-Z]) ]]; then
        chain="${BASH_REMATCH[1]}"
    elif [[ "$header" =~ _chain_([0-9]+) ]]; then
        # For files like 5D0Q_chain_1.fasta, map chain numbers to letters
        local chain_num="${BASH_REMATCH[1]}"
        case $chain_num in
            1) chain="A" ;;
            2) chain="B" ;;
            3) chain="C" ;;
            4) chain="D" ;;
            5) chain="E" ;;
            6) chain="F" ;;
            7) chain="G" ;;
            *) chain="A" ;;
        esac
    fi
    
    echo "Detected chain: $chain"
    
    # Extract sequence (skip header line)
    local sequence=$(tail -n +2 "$bepipred_fasta" | tr -d '\n' | tr -d ' ')
    local seq_length=${#sequence}
    
    echo "Sequence length: $seq_length residues"
    
    # Filter numbering file for the detected chain
    local chain_numbering=$(grep "^$chain " "$numbering_file")
    local chain_count=$(echo "$chain_numbering" | wc -l)
    
    echo "Chain $chain has $chain_count residues in PDB"
    
    # Create mapping output
    echo "# Mapping between FASTA sequence position and PDB numbering" > "$output_file"
    echo "# BepiPred FASTA: $bepipred_fasta" >> "$output_file"
    echo "# Structure: $structure_id" >> "$output_file"
    echo "# Chain: $chain" >> "$output_file"
    echo "# FASTA_POS	PDB_CHAIN	PDB_NUM	PDB_AA	FASTA_AA	MATCH" >> "$output_file"
    
    # Create arrays for PDB numbering
    local -a pdb_nums
    local -a pdb_aas
    
    while IFS=' ' read -r pdb_chain pdb_num pdb_aa; do
        if [[ "$pdb_chain" == "$chain" ]]; then
            pdb_nums+=("$pdb_num")
            pdb_aas+=("$pdb_aa")
        fi
    done <<< "$chain_numbering"
    
    # Convert PDB sequence to 1-letter codes for alignment
    local pdb_sequence=""
    for pdb_aa in "${pdb_aas[@]}"; do
        pdb_sequence+=$(convert_aa_code "$pdb_aa")
    done
    
    echo "PDB sequence length: ${#pdb_sequence}"
    
    # Find alignment between FASTA and PDB sequences
    local best_offset=0
    local best_matches=0
    local uppercase_sequence=$(echo "$sequence" | tr '[:lower:]' '[:upper:]')
    
    # Try different offsets to find best alignment
    echo "Finding best alignment..."
    for ((offset=0; offset<=100 && offset<${#pdb_sequence}; offset++)); do
        local matches_at_offset=0
        local max_compare=$((seq_length < (${#pdb_sequence} - offset) ? seq_length : (${#pdb_sequence} - offset)))
        
        for ((j=0; j<max_compare && j<50; j++)); do  # Check first 50 positions
            local fasta_aa="${uppercase_sequence:$j:1}"
            local pdb_aa="${pdb_sequence:$((offset+j)):1}"
            if [[ "$fasta_aa" == "$pdb_aa" ]]; then
                ((matches_at_offset++))
            fi
        done
        
        if ((matches_at_offset > best_matches)); then
            best_matches=$matches_at_offset
            best_offset=$offset
        fi
    done
    
    echo "Best alignment: offset=$best_offset, matches=$best_matches"
    
    # Map each FASTA position to PDB position using best offset
    local matches=0
    local mismatches=0
    
    for ((i=0; i<seq_length; i++)); do
        local fasta_aa="${uppercase_sequence:$i:1}"
        local fasta_pos=$((i+1))
        local pdb_index=$((best_offset + i))
        
        if ((pdb_index < ${#pdb_nums[@]})); then
            local pdb_aa="${pdb_aas[$pdb_index]}"
            local pdb_num="${pdb_nums[$pdb_index]}"
            local pdb_aa_1letter=$(convert_aa_code "$pdb_aa")
            
            local match_status="MATCH"
            if [[ "$fasta_aa" != "$pdb_aa_1letter" ]]; then
                match_status="MISMATCH"
                ((mismatches++))
            else
                ((matches++))
            fi
            
            echo -e "${fasta_pos}\t${chain}\t${pdb_num}\t${pdb_aa}\t${fasta_aa}\t${match_status}" >> "$output_file"
        else
            # FASTA position beyond PDB structure
            echo -e "${fasta_pos}\t${chain}\t-\t-\t${fasta_aa}\tBEYOND_PDB" >> "$output_file"
        fi
    done
    
    echo "Mapping complete: $matches matches, $mismatches mismatches"
    
    # Add summary to output file
    echo "" >> "$output_file"
    echo "# SUMMARY" >> "$output_file"
    echo "# Total mapped positions: $((matches + mismatches))" >> "$output_file"
    echo "# Matches: $matches" >> "$output_file"
    echo "# Mismatches: $mismatches" >> "$output_file"
    echo "# Match rate: $(( matches * 100 / (matches + mismatches) ))%" >> "$output_file"
    
    echo "---"
}

# Function to convert 3-letter amino acid codes to 1-letter
convert_aa_code() {
    case "$1" in
        ALA) echo "A" ;;
        ARG) echo "R" ;;
        ASN) echo "N" ;;
        ASP) echo "D" ;;
        CYS) echo "C" ;;
        GLU) echo "E" ;;
        GLN) echo "Q" ;;
        GLY) echo "G" ;;
        HIS) echo "H" ;;
        ILE) echo "I" ;;
        LEU) echo "L" ;;
        LYS) echo "K" ;;
        MET) echo "M" ;;
        PHE) echo "F" ;;
        PRO) echo "P" ;;
        SER) echo "S" ;;
        THR) echo "T" ;;
        TRP) echo "W" ;;
        TYR) echo "Y" ;;
        VAL) echo "V" ;;
        *) echo "X" ;;
    esac
}

# Main execution
if [ $# -eq 2 ]; then
    # Process specific file
    map_fasta_to_pdb "$1" "$2"
elif [ $# -eq 0 ]; then
    # Process all BepiPred FASTA files
    echo "Processing all BepiPred FASTA files..."
    echo "========================================"
    
    find results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/ -name "*_Bcell_epitope_preds.fasta" -type f | while read -r fasta_file; do
        # Extract structure ID from filename
        local basename=$(basename "$fasta_file")
        local structure_id="${basename%_Bcell_epitope_preds.fasta}"
        structure_id="${structure_id#*_}"  # Remove gene name prefix
        
        map_fasta_to_pdb "$fasta_file" "$structure_id"
    done
    
    echo "========================================"
    echo "Processing complete!"
    echo ""
    echo "Summary of mapping files created:"
    find results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/ -name "*_pdb_mapping.txt" -type f | wc -l
    echo "mapping files created"
else
    echo "Usage:"
    echo "  $0                              # Process all BepiPred FASTA files"
    echo "  $0 [fasta_file] [structure_id]  # Process specific file"
    echo ""
    echo "Examples:"
    echo "  $0 results/.../bamA_5D0Q_Bcell_epitope_preds.fasta 5D0Q"
    echo "  $0  # Process all files automatically"
fi