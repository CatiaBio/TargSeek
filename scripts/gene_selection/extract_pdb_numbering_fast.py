#!/usr/bin/env python3

"""
Fast PDB numbering extraction script - optimized version
Only extracts chain start/end ranges from structure files without FASTA association
This version skips the expensive FASTA header parsing and complex mapping logic

Performance improvements:
- No FASTA file reading or header parsing
- No complex chain ID mapping between FASTA and PDB
- Simplified structure parsing focused on chain ranges only
- Reduced file I/O operations
"""

import os
import sys
import gzip
import pandas as pd
from pathlib import Path

def extract_structure_chain_info_fast(structure_file_path):
    """
    Fast extraction of chain information from PDB/mmCIF files
    Only extracts chain IDs with their start/end residue numbers
    
    Args:
        structure_file_path: Path to structure file
    
    Returns:
        dict: {chain_id: {'min_res': int, 'max_res': int}}
    """
    structure_path = Path(structure_file_path)
    
    if not structure_path.exists():
        return {}
    
    chains = {}
    
    try:
        # Handle compressed files
        if structure_path.name.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
        
        with open_func(structure_path, mode) as f:
            # Quick format detection
            first_line = f.readline().strip()
            f.seek(0)
            
            if structure_path.name.endswith(('.cif', '.cif.gz')) or first_line.startswith('data_'):
                chains = _parse_mmcif_chains_fast(f)
            else:
                chains = _parse_pdb_chains_fast(f)
                
    except Exception as e:
        print(f"Error reading {structure_file_path}: {e}")
        return {}
    
    return chains

def _parse_pdb_chains_fast(file_handle):
    """Fast PDB format parsing - only reads ATOM/HETATM lines"""
    chains = {}
    
    for line in file_handle:
        if line.startswith(('ATOM', 'HETATM')):
            chain_id = line[21].strip() or 'A'
            
            try:
                res_num = int(line[22:26].strip())
                
                if chain_id not in chains:
                    chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                else:
                    chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                    chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                    
            except ValueError:
                continue
    
    return chains

def _parse_mmcif_chains_fast(file_handle):
    """Fast mmCIF format parsing - simplified column detection"""
    chains = {}
    in_atom_section = False
    chain_col = None
    seq_id_col = None
    headers_found = 0
    
    for line in file_handle:
        line = line.strip()
        
        # Find column headers
        if line.startswith('_atom_site.'):
            in_atom_section = True
            if 'auth_asym_id' in line:
                chain_col = headers_found
            elif 'auth_seq_id' in line:
                seq_id_col = headers_found
            headers_found += 1
            
        # Process data lines
        elif in_atom_section and not line.startswith('_') and not line.startswith('#'):
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    parts = line.split()
                    if chain_col is not None and seq_id_col is not None and len(parts) > max(chain_col, seq_id_col):
                        chain_id = parts[chain_col]
                        res_num = int(parts[seq_id_col])
                        
                        if chain_id not in chains:
                            chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                        else:
                            chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                            chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                except (ValueError, IndexError):
                    continue
                    
        elif line.startswith('#') and in_atom_section:
            break
    
    return chains

def create_fast_pdb_numbering_mapping(structures_dir, mapping_file, output_tsv, analysis, paramset):
    """
    Create simplified TSV with structure chain ranges only
    No FASTA association - just pure structure data
    
    Args:
        structures_dir: Directory containing PDB files
        mapping_file: Original mapping file (used only for structure list)
        output_tsv: Output TSV file
        analysis: Analysis identifier
        paramset: Parameter set identifier
    
    Returns:
        dict: Processing results
    """
    structures_path = Path(structures_dir)
    mapping_path = Path(mapping_file)
    
    if not structures_path.exists():
        print(f"Error: Structures directory {structures_dir} does not exist")
        return {}
        
    if not mapping_path.exists():
        print(f"Error: Mapping file {mapping_file} does not exist")
        return {}
    
    results = {
        'processed': 0,
        'failed': 0,
        'total_structures': 0,
        'unique_structures': 0
    }
    
    # Load mapping to get structure list
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        print(f"Loaded {len(mapping_df)} entries from mapping file")
    except Exception as e:
        print(f"Error loading mapping file: {e}")
        return {}
    
    # Get unique structures to process
    unique_structures = mapping_df[['structure_id', 'structure_path', 'structure_type', 'gene_name']].drop_duplicates()
    results['unique_structures'] = len(unique_structures)
    
    tsv_data = []
    
    print(f"Processing {len(unique_structures)} unique structures")
    print("=" * 50)
    
    for _, row in unique_structures.iterrows():
        structure_path = Path(row['structure_path'])
        structure_id = row['structure_id']
        structure_type = row['structure_type']
        gene_name = row['gene_name']
        
        if not structure_path.exists():
            print(f"Warning: Structure file not found: {structure_path}")
            continue
            
        results['total_structures'] += 1
        
        print(f"Processing {structure_id} ({structure_type}) for {gene_name}")
        
        # Extract chain information
        chain_info = extract_structure_chain_info_fast(str(structure_path))
        
        if chain_info:
            results['processed'] += 1
            
            # Add entry for each chain
            for chain_id, residue_info in chain_info.items():
                tsv_data.append({
                    'Gene': gene_name,
                    'Structure_ID': structure_id,
                    'Chain': chain_id,
                    'Start': residue_info['min_res'],
                    'End': residue_info['max_res'],
                    'Range': f"{residue_info['min_res']}-{residue_info['max_res']}",
                    'Structure_Path': str(structure_path),
                    'Structure_Type': structure_type
                })
        else:
            results['failed'] += 1
            print(f"  Failed to extract chain info for {structure_path}")
    
    # Save results
    if tsv_data:
        df = pd.DataFrame(tsv_data)
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"\nCreated fast numbering TSV: {output_tsv}")
        print(f"Total chain entries: {len(tsv_data)}")
    else:
        print("No data to write to TSV file")
        # Create empty TSV with headers
        df = pd.DataFrame(columns=['Gene', 'Structure_ID', 'Chain', 'Start', 'End', 'Range', 'Structure_Path', 'Structure_Type'])
        df.to_csv(output_tsv, sep='\t', index=False)
    
    print("=" * 50)
    print(f"Fast processing complete!")
    print(f"Unique structures processed: {results['processed']}/{results['unique_structures']}")
    print(f"Failed: {results['failed']}")
    print(f"Total chain entries created: {len(tsv_data)}")
    
    return results

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        mapping_file = snakemake.input.mapping_file
        numbering_tsv = snakemake.output.numbering_tsv
        sentinel_file = snakemake.output.numbering_sentinel
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        structures_dir = snakemake.params.structures_dir
        
        print(f"Fast PDB numbering extraction for {analysis}_{paramset}")
        print(f"Structures directory: {structures_dir}")
        print(f"Mapping file: {mapping_file}")
        print(f"Output TSV: {numbering_tsv}")
        
    except NameError:
        # Test mode - command line arguments
        if len(sys.argv) != 4:
            print("Usage: python extract_pdb_numbering_fast.py <analysis> <paramset> <structures_dir>")
            print("Example: python extract_pdb_numbering_fast.py analysis1 params1 data/protein_structures/pdb_files")
            sys.exit(1)
        
        analysis = sys.argv[1]
        paramset = sys.argv[2]
        structures_dir = sys.argv[3]
        
        mapping_file = f"data/protein_structures/{analysis}_{paramset}_fasta_structure_mapping.tsv"
        numbering_tsv = f"data/protein_structures/pdb_files/{analysis}_{paramset}_pdb_numbering_mapping_fast.tsv"
        sentinel_file = f"data/protein_structures/pdb_files/{analysis}_{paramset}_pdb_numbering_fast_complete.sentinel"
    
    # Create the fast mapping
    results = create_fast_pdb_numbering_mapping(structures_dir, mapping_file, numbering_tsv, analysis, paramset)
    
    # Create sentinel file
    with open(sentinel_file, 'w') as f:
        f.write(f"Fast PDB numbering extraction completed for {analysis}_{paramset}\n")
        f.write(f"Unique structures processed: {results.get('processed', 0)}\n")
        f.write(f"Total chain entries: {results.get('total_structures', 0)}\n")
        f.write(f"Output TSV: {numbering_tsv}\n")
    
    print(f"Fast PDB numbering extraction completed. Sentinel: {sentinel_file}")

if __name__ == "__main__":
    main()