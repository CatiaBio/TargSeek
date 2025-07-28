#!/bin/bash

# Script to create BepiPred PDB numbering mapping files
# This finds where the PDB structure sequence matches within the FASTA sequence
# and creates a mapping file for BepiPred analysis

create_bepipred_mapping() {
    local bepipred_fasta="$1"
    local structure_id="$2"
    
    # Determine gene name from the fasta file path
    local gene_name=$(echo "$bepipred_fasta" | sed -n 's|.*/epitope_predictions_bepipred/\([^/]*\)/.*|\1|p')
    local numbering_file="data/protein_structures/${gene_name}/${structure_id}_numbering.txt"
    local output_dir=$(dirname "$bepipred_fasta")
    local output_file="${output_dir}/${gene_name}_${structure_id}_pdb_numbering.txt"
    
    echo "Processing: $bepipred_fasta"
    echo "Gene: $gene_name"
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
    
    # Try to extract chain from header if present
    if [[ "$header" =~ [Cc]hain[s]?[[:space:]]*([A-Z]) ]]; then
        chain="${BASH_REMATCH[1]}"
    elif [[ "$structure_id" =~ _chain_([0-9]+) ]]; then
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
    
    # Extract FASTA sequence (skip header line, convert to uppercase)
    local fasta_sequence=$(tail -n +2 "$bepipred_fasta" | tr -d '\n' | tr -d ' ' | tr '[:lower:]' '[:upper:]')
    local fasta_length=${#fasta_sequence}
    
    echo "FASTA sequence length: $fasta_length residues"
    
    # Extract PDB sequence for the specific chain
    local pdb_sequence=""
    local -a pdb_nums
    local -a pdb_aas
    
    while IFS=' ' read -r pdb_chain pdb_num pdb_aa; do
        if [[ "$pdb_chain" == "$chain" ]]; then
            pdb_nums+=("$pdb_num")
            pdb_aas+=("$pdb_aa")
            pdb_sequence+=$(convert_aa_code "$pdb_aa")
        fi
    done < "$numbering_file"
    
    local pdb_length=${#pdb_sequence}
    echo "PDB sequence length: $pdb_length residues"
    
    # Find where PDB sequence matches in FASTA sequence
    echo "Searching for PDB sequence within FASTA sequence..."
    
    local best_start_pos=0
    local best_matches=0
    local found_match=false
    
    # Try all possible starting positions in FASTA
    for ((start=0; start<=(fasta_length-pdb_length+50); start++)); do
        local matches=0
        local max_check=$((pdb_length < 200 ? pdb_length : 200))  # Check first 200 residues for better accuracy
        local actual_check=$((start + max_check <= fasta_length ? max_check : fasta_length - start))
        
        if ((actual_check <= 0)); then
            continue
        fi
        
        for ((i=0; i<actual_check; i++)); do
            if ((start + i < fasta_length && i < pdb_length)); then
                local fasta_aa="${fasta_sequence:$((start+i)):1}"
                local pdb_aa="${pdb_sequence:$i:1}"
                
                if [[ "$fasta_aa" == "$pdb_aa" ]]; then
                    ((matches++))
                fi
            fi
        done
        
        # Calculate match percentage for this position
        if ((actual_check > 0)); then
            local match_percent=$((matches * 100 / actual_check))
            
            if ((matches > best_matches)); then
                best_matches=$matches
                best_start_pos=$start
                if ((match_percent >= 70)); then
                    found_match=true
                fi
            fi
        fi
    done
    
    if [[ "$found_match" == "false" ]]; then
        echo "WARNING: Could not find high-quality alignment (>70%)"
        echo "Best match: $best_matches matches at position $((best_start_pos+1))"
        echo "Using best available alignment"
    else
        echo "Found good alignment: FASTA position $((best_start_pos+1)) matches PDB start"
        echo "Match quality: $best_matches matches"
    fi
    
    # Create the mapping file
    echo "# BepiPred PDB Numbering Mapping" > "$output_file"
    echo "# Generated from: $bepipred_fasta" >> "$output_file"
    echo "# PDB Structure: $structure_id" >> "$output_file"
    echo "# Chain: $chain" >> "$output_file"
    echo "# FASTA length: $fasta_length" >> "$output_file"
    echo "# PDB length: $pdb_length" >> "$output_file"
    echo "# PDB sequence starts at FASTA position: $((best_start_pos+1))" >> "$output_file"
    echo "# PDB sequence ends at FASTA position: $((best_start_pos+pdb_length))" >> "$output_file"
    echo "#" >> "$output_file"
    echo "# FASTA_POS	PDB_NUM	PDB_AA	FASTA_AA	STATUS" >> "$output_file"
    
    # Create detailed mapping
    local total_matches=0
    local total_mismatches=0
    
    for ((i=0; i<pdb_length; i++)); do
        local fasta_pos=$((best_start_pos + i + 1))
        local pdb_num="${pdb_nums[$i]}"
        local pdb_aa="${pdb_aas[$i]}"
        local pdb_aa_1letter=$(convert_aa_code "$pdb_aa")
        
        if ((fasta_pos <= fasta_length)); then
            local fasta_aa="${fasta_sequence:$((fasta_pos-1)):1}"
            local status="MATCH"
            
            if [[ "$fasta_aa" != "$pdb_aa_1letter" ]]; then
                status="MISMATCH"
                ((total_mismatches++))
            else
                ((total_matches++))
            fi
            
            echo -e "${fasta_pos}\t${pdb_num}\t${pdb_aa}\t${fasta_aa}\t${status}" >> "$output_file"
        else
            echo -e "${fasta_pos}\t${pdb_num}\t${pdb_aa}\t-\tBEYOND_FASTA" >> "$output_file"
        fi
    done
    
    # Add summary
    echo "" >> "$output_file"
    echo "# SUMMARY" >> "$output_file"
    echo "# Total positions mapped: $pdb_length" >> "$output_file"
    echo "# Matches: $total_matches" >> "$output_file"
    echo "# Mismatches: $total_mismatches" >> "$output_file"
    if ((total_matches + total_mismatches > 0)); then
        echo "# Match rate: $(( total_matches * 100 / (total_matches + total_mismatches) ))%" >> "$output_file"
    fi
    
    echo "Mapping complete: $total_matches matches, $total_mismatches mismatches"
    
    # Show key information
    echo "PDB structure covers FASTA positions: $((best_start_pos+1))-$((best_start_pos+pdb_length))"
    echo "PDB residue range: ${pdb_nums[0]}-${pdb_nums[-1]}"
    echo "Mapping file created: $output_file"
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
    create_bepipred_mapping "$1" "$2"
elif [ $# -eq 0 ]; then
    # Process all BepiPred FASTA files
    echo "Creating BepiPred PDB numbering mappings for all structures..."
    echo "================================================================"
    
    find results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/ -name "*_Bcell_epitope_preds.fasta" -type f | while read -r fasta_file; do
        # Extract structure ID from filename
        local basename=$(basename "$fasta_file")
        local structure_id="${basename%_Bcell_epitope_preds.fasta}"
        structure_id="${structure_id#*_}"  # Remove gene name prefix
        
        create_bepipred_mapping "$fasta_file" "$structure_id"
    done
    
    echo "================================================================"
    echo "Processing complete!"
    echo ""
    echo "Summary of BepiPred PDB mapping files created:"
    find results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/ -name "*_pdb_numbering.txt" -type f | wc -l
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