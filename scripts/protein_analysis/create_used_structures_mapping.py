#!/usr/bin/env python3

"""
Script to create mapping of used structures from selected_3d_paths files
for Snakemake analysis pipeline. Creates used_structures_mapping.tsv in epitope_predictions_bepipred folder.
"""

import os
import sys
import re
from pathlib import Path

def parse_fasta_header(header_line):
    """
    Parse FASTA header to extract chain information
    Example: >5D0Q_5D0Q_1 PDB:5D0Q 5D0Q_1|Chains A, E[auth F]|Outer membrane protein assembly factor BamA|Escherichia coli K-12 (562)
    """
    # Extract structure ID from PDB field
    pdb_match = re.search(r'PDB:(\w+)', header_line)
    if not pdb_match:
        return None, None
    
    structure_id = pdb_match.group(1)
    
    # Extract chain information
    chain_match = re.search(r'Chains?\s+([A-Za-z0-9, \[\]]+)', header_line)
    if not chain_match:
        # Try alternative format or assume single chain
        return structure_id, "A"
    
    chain_text = chain_match.group(1)
    
    # Parse chain text - could be "A", "A, B", "A, E[auth F]", etc.
    # Extract the first primary chain (before any auth notation)
    primary_chain = chain_text.split(',')[0].strip()
    
    # Remove any bracketed authorization info
    primary_chain = re.sub(r'\[.*?\]', '', primary_chain).strip()
    
    return structure_id, primary_chain

def process_selected_paths(selected_paths_file):
    """Process selected 3D paths file and extract structure usage info"""
    
    used_structures = []
    
    if not os.path.exists(selected_paths_file):
        print(f"Warning: {selected_paths_file} does not exist")
        return used_structures
    
    with open(selected_paths_file, 'r') as f:
        for line in f:
            fasta_path = line.strip()
            if not fasta_path or fasta_path.startswith('#'):
                continue
            
            # Extract gene name from path
            # Example: data/protein_structures/bamA/5D0Q_chain_1.fasta
            path_parts = Path(fasta_path).parts
            if len(path_parts) < 3:
                print(f"Warning: Cannot parse path {fasta_path}")
                continue
            
            gene_name = path_parts[-2]  # bamA
            fasta_filename = path_parts[-1]  # 5D0Q_chain_1.fasta
            
            # Read FASTA header to get actual chain info
            try:
                with open(fasta_path, 'r') as fasta_file:
                    header_line = fasta_file.readline().strip()
                    
                    if not header_line.startswith('>'):
                        print(f"Warning: Invalid FASTA header in {fasta_path}")
                        continue
                    
                    structure_id, chain = parse_fasta_header(header_line)
                    
                    if structure_id and chain:
                        used_structures.append({
                            'gene': gene_name,
                            'structure_id': structure_id,
                            'chain': chain,
                            'fasta_path': fasta_path
                        })
                        print(f"Found: {gene_name} -> {structure_id} -> Chain {chain}")
                    else:
                        print(f"Warning: Could not parse structure/chain from {fasta_path}")
                        
            except Exception as e:
                print(f"Error reading {fasta_path}: {e}")
                continue
    
    return used_structures

def create_bepipred_structure_mapping(used_structures, output_file):
    """Create mapping file for BepiPred structure usage"""
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("Gene\tStructure_ID\tChain\tFASTA_Path\n")
        
        for struct in used_structures:
            f.write(f"{struct['gene']}\t{struct['structure_id']}\t{struct['chain']}\t{struct['fasta_path']}\n")
    
    print(f"Created mapping file: {output_file}")
    print(f"Total structures mapped: {len(used_structures)}")

def main():
    """Main function for Snakemake integration"""
    
    # Get parameters from snakemake object or command line
    if 'snakemake' in globals():
        # Running from Snakemake
        gram_negative_file = snakemake.input.selected_3d_paths_negative
        gram_positive_file = snakemake.input.selected_3d_paths_positive
        output_file = snakemake.output.mapping_tsv
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
    else:
        # Running from command line for testing
        if len(sys.argv) != 5:
            print("Usage: python create_used_structures_mapping.py <analysis> <paramset> <gram_negative_file> <gram_positive_file>")
            print("Example: python create_used_structures_mapping.py analysis1 params1 selected_3d_paths_gram_negative.txt selected_3d_paths_gram_positive.txt")
            sys.exit(1)
        
        analysis = sys.argv[1]
        paramset = sys.argv[2]
        gram_negative_file = sys.argv[3]
        gram_positive_file = sys.argv[4]
        
        # Define output path manually for testing
        output_dir = f"results/{analysis}_{paramset}/protein_analysis/epitope_predictions_bepipred"
        output_file = os.path.join(output_dir, "used_structures_mapping.tsv")
    
    print(f"Creating used structures mapping for {analysis}_{paramset}")
    print(f"Input files:")
    print(f"  Gram-negative: {gram_negative_file}")
    print(f"  Gram-positive: {gram_positive_file}")
    print(f"Output file: {output_file}")
    print()
    
    all_used_structures = []
    
    # Process gram-negative file
    if os.path.exists(gram_negative_file):
        print(f"Processing gram-negative structures: {gram_negative_file}")
        gram_neg_structures = process_selected_paths(gram_negative_file)
        all_used_structures.extend(gram_neg_structures)
        print(f"Found {len(gram_neg_structures)} gram-negative structures")
    else:
        print(f"Warning: {gram_negative_file} not found")
    
    # Process gram-positive file  
    if os.path.exists(gram_positive_file):
        print(f"Processing gram-positive structures: {gram_positive_file}")
        gram_pos_structures = process_selected_paths(gram_positive_file)
        all_used_structures.extend(gram_pos_structures)
        print(f"Found {len(gram_pos_structures)} gram-positive structures")
    else:
        print(f"Warning: {gram_positive_file} not found")
    
    if not all_used_structures:
        print("Warning: No structures found to process")
        # Create empty mapping file
        create_bepipred_structure_mapping([], output_file)
        return
    
    # Create mapping file
    create_bepipred_structure_mapping(all_used_structures, output_file)
    
    # Print summary by gene
    print("\nSummary by gene:")
    gene_counts = {}
    for struct in all_used_structures:
        gene = struct['gene']
        if gene not in gene_counts:
            gene_counts[gene] = []
        gene_counts[gene].append(f"{struct['structure_id']}:{struct['chain']}")
    
    for gene, structures in sorted(gene_counts.items()):
        print(f"  {gene}: {', '.join(structures)}")
    
    print(f"\nMapping completed successfully!")

if __name__ == "__main__":
    main()