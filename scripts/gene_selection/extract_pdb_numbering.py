#!/usr/bin/env python3

"""
Script to extract PDB amino acid numbering from 3D structure files
Creates a consolidated TSV file mapping genes to PDB structures with numbering information

Updated to work with new architecture:
- FASTA files in gene directories: data/protein_structures/{gene_name}/
- 3D structures in shared directory: data/protein_structures/pdb_files/
- Uses mapping file to associate FASTA files with structures
"""

import os
import sys
import gzip
import pandas as pd
import re
from pathlib import Path

def extract_standard_chains_from_fasta(fasta_file_path):
    """
    Extract standard chain IDs from FASTA header
    
    Args:
        fasta_file_path: Path to FASTA file
    
    Returns:
        list: List of standard chain IDs (e.g., ['NA', 'OA', 'PA'] from multi-chain or ['B'] from single chain)
    """
    try:
        with open(fasta_file_path, 'r') as f:
            header_line = f.readline().strip()
            
        # Handle multi-chain case: "Chains AI[auth B1], AJ[auth B2], AK[auth B7], ..." or "Chains A, E[auth F]" or "Chains A, B"
        chains_match = re.search(r'Chains\s+([^|]+)', header_line)
        if chains_match:
            chains_text = chains_match.group(1).strip()
            # Extract all chain IDs (both with and without [auth ...] notation)
            # This regex captures: A, E[auth F], NA[auth k], B, etc.
            # Updated to handle 1-2 letter chain IDs including numbers
            chain_matches = re.findall(r'([A-Z]{1,2}[A-Z0-9]*)(?:\[auth\s+[^]]+\])?', chains_text)
            if chain_matches:
                return chain_matches
        
        # Handle single chain case: "Chain A[auth Y]" or just "Chain A" or "Chain B"
        single_chain_match = re.search(r'Chain\s+([A-Z]{1,2}[A-Z0-9]*)(?:\[auth\s+[^]]+\])?', header_line)
        if single_chain_match:
            return [single_chain_match.group(1)]
            
    except Exception as e:
        print(f"Error reading FASTA header from {fasta_file_path}: {e}")
    
    return []

def extract_structure_chain_info(structure_file_path):
    """
    Extract chain information (amino acid residue numbering) from PDB or mmCIF structure files
    
    Args:
        structure_file_path: Path to structure file (PDB or mmCIF format, compressed or uncompressed)
    
    Returns:
        dict: Dictionary mapping chain IDs to residue information
              Format: {chain_id: {'min_res': int, 'max_res': int}}
    """
    
    structure_path = Path(structure_file_path)
    
    # Check if file exists
    if not structure_path.exists():
        print(f"Structure file not found: {structure_file_path}")
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
            # Determine file format
            first_line = f.readline().strip()
            f.seek(0)  # Reset to beginning
            
            if structure_path.name.endswith(('.cif', '.cif.gz')) or first_line.startswith('data_'):
                # mmCIF format
                chains = _parse_mmcif_chains(f)
            else:
                # PDB format
                chains = _parse_pdb_chains(f)
                
    except Exception as e:
        print(f"Error reading structure file {structure_file_path}: {e}")
        return {}
    
    return chains

def _parse_pdb_chains(file_handle):
    """Parse PDB format file for chain information"""
    chains = {}
    
    for line in file_handle:
        if line.startswith(('ATOM', 'HETATM')):
            # Extract chain ID and residue number
            chain_id = line[21].strip()
            if not chain_id:
                chain_id = 'A'  # Default chain if missing
            
            try:
                res_num = int(line[22:26].strip())
                
                if chain_id not in chains:
                    chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                else:
                    chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                    chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                    
            except ValueError:
                continue  # Skip invalid residue numbers
    
    return chains

def _parse_mmcif_chains(file_handle):
    """Parse mmCIF format file for chain information"""
    chains = {}
    in_atom_site = False
    chain_col = None
    seq_id_col = None
    
    for line in file_handle:
        line = line.strip()
        
        if line.startswith('_atom_site.'):
            in_atom_site = True
            # Find column indices
            if 'auth_asym_id' in line:
                chain_col = len([x for x in line.split('.')[0].split('_') if x]) - 1
            elif 'auth_seq_id' in line:
                seq_id_col = len([x for x in line.split('.')[0].split('_') if x]) - 1
        elif in_atom_site and not line.startswith('_') and not line.startswith('#'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    parts = line.split()
                    if chain_col is not None and seq_id_col is not None:
                        chain_id = parts[chain_col] if chain_col < len(parts) else 'A'
                        res_num = int(parts[seq_id_col]) if seq_id_col < len(parts) else 0
                        
                        if chain_id not in chains:
                            chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                        else:
                            chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                            chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                except (ValueError, IndexError):
                    continue
        elif line.startswith('#') and in_atom_site:
            break  # End of atom_site loop
    
    return chains

def create_pdb_numbering_mapping(base_dir, structures_dir, mapping_file, output_tsv, analysis, paramset):
    """
    Create consolidated TSV file mapping genes to PDB structures with numbering information
    Uses the new architecture: FASTA files in gene directories, structures in shared pdb_files/
    
    Args:
        base_dir: Base directory containing gene folders with FASTA files
        structures_dir: Directory containing shared 3D structure files
        mapping_file: TSV file mapping FASTA files to structure files
        output_tsv: Output TSV file path
        analysis: Analysis identifier
        paramset: Parameter set identifier
    
    Returns:
        dict: Summary of processing results
    """
    base_path = Path(base_dir)
    structures_path = Path(structures_dir)
    mapping_path = Path(mapping_file)
    
    if not base_path.exists():
        print(f"Error: Base directory {base_dir} does not exist")
        return {}
    
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
        'genes_found': 0
    }
    
    # Load the mapping file
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        print(f"Loaded {len(mapping_df)} entries from mapping file")
    except Exception as e:
        print(f"Error loading mapping file {mapping_file}: {e}")
        return {}
    
    # Collect data for TSV
    tsv_data = []
    
    # Process entries from mapping file
    print(f"Processing {len(mapping_df)} FASTA-structure mappings")
    print("=" * 50)
    
    # Group by gene for processing
    genes_processed = set()
    
    for _, mapping_row in mapping_df.iterrows():
        fasta_path = Path(mapping_row['fasta_path'])
        structure_path = Path(mapping_row['structure_path'])
        gene_name = mapping_row['gene_name']
        structure_id = mapping_row['structure_id']
        structure_type = mapping_row['structure_type']
        
        # Track genes
        if gene_name not in genes_processed:
            results['genes_found'] += 1
            genes_processed.add(gene_name)
        
        # Check if FASTA file exists
        if not fasta_path.exists():
            print(f"Warning: FASTA file not found: {fasta_path}")
            continue
            
        # Check if structure file exists
        if not structure_path.exists():
            print(f"Warning: Structure file not found: {structure_path}")
            continue
        
        results['total_structures'] += 1
        
        print(f"Processing {gene_name}/{structure_id} ({structure_type})")
        
        # Extract chain information using the unified function
        chain_info = extract_structure_chain_info(str(structure_path))
        
        if chain_info:
            results['processed'] += 1
            
            # Extract standard chain IDs from FASTA header
            standard_chain_ids = extract_standard_chains_from_fasta(fasta_path)
            
            # Get all author chain IDs from structure info
            author_chain_ids = list(chain_info.keys())
            
            # Use the chains directly from the FASTA header
            if standard_chain_ids:
                # Find the first matching chain in the structure for residue numbering
                author_chain_id = None
                matched_chain = None
                
                for standard_chain_id in standard_chain_ids:
                    if standard_chain_id in author_chain_ids:
                        author_chain_id = standard_chain_id
                        matched_chain = standard_chain_id
                        break
                
                # If no exact match found, use first available chain from structure
                if not author_chain_id and author_chain_ids:
                    author_chain_id = author_chain_ids[0]
                    matched_chain = standard_chain_ids[0]  # Use first chain from header for display
                    print(f"  Warning: No exact chain match found, using structure chain {author_chain_id} for header chains {standard_chain_ids}")
                
                if author_chain_id:
                    # Create single entry with comma-separated chain list for display
                    display_chains = ','.join(standard_chain_ids)
                    
                    tsv_data.append({
                        'Gene': gene_name,
                        'Structure_ID': structure_id,
                        'Chain': display_chains,
                        'Start': chain_info[author_chain_id]['min_res'],
                        'End': chain_info[author_chain_id]['max_res'],
                        'Range': f"{chain_info[author_chain_id]['min_res']}-{chain_info[author_chain_id]['max_res']}",
                        'FASTA_Path': str(fasta_path),
                        'Structure_Path': str(structure_path),
                        'Structure_Type': structure_type
                    })
            else:
                # Fallback: no chains found in header, use first available chain from structure
                if author_chain_ids:
                    author_chain_id = author_chain_ids[0]
                    
                    tsv_data.append({
                        'Gene': gene_name,
                        'Structure_ID': structure_id,
                        'Chain': author_chain_id,
                        'Start': chain_info[author_chain_id]['min_res'],
                        'End': chain_info[author_chain_id]['max_res'],
                        'Range': f"{chain_info[author_chain_id]['min_res']}-{chain_info[author_chain_id]['max_res']}",
                        'FASTA_Path': str(fasta_path),
                        'Structure_Path': str(structure_path),
                        'Structure_Type': structure_type
                    })
                    print(f"  Warning: No chains found in FASTA header, using {author_chain_id}")
        else:
            results['failed'] += 1
            print(f"  Failed to extract chain info for {structure_path}")
    
    # Create DataFrame and save to TSV
    if tsv_data:
        df = pd.DataFrame(tsv_data)
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"\nCreated mapping TSV: {output_tsv}")
        print(f"Total entries: {len(tsv_data)}")
    else:
        print(f"No data to write to TSV file")
        # Create empty TSV with headers
        df = pd.DataFrame(columns=['Gene', 'Structure_ID', 'Chain', 'Start', 'End', 'Range', 'FASTA_Path', 'Structure_Path', 'Structure_Type'])
        df.to_csv(output_tsv, sep='\t', index=False)
    
    print("=" * 50)
    print(f"Processing complete!")
    print(f"Genes with structures: {results['genes_found']}")
    print(f"Successfully processed: {results['processed']} PDB files")
    print(f"Failed: {results['failed']} PDB files")
    print(f"Total structures: {results['total_structures']}")
    
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
        base_dir = snakemake.params.base_dir
        structures_dir = snakemake.params.structures_dir
        
        print(f"PDB numbering extraction for {analysis}_{paramset}")
        print(f"Base directory: {base_dir}")
        print(f"Structures directory: {structures_dir}")
        print(f"Mapping file: {mapping_file}")
        print(f"Output TSV: {numbering_tsv}")
        
    except NameError:
        # Test mode - command line arguments
        if len(sys.argv) != 5:
            print("Usage: python extract_pdb_numbering.py <analysis> <paramset> <base_dir> <structures_dir>")
            print("Example: python extract_pdb_numbering.py analysis1 params1 data/protein_structures data/protein_structures/pdb_files")
            sys.exit(1)
        
        analysis = sys.argv[1]
        paramset = sys.argv[2]
        base_dir = sys.argv[3]
        structures_dir = sys.argv[4]
        
        mapping_file = f"data/protein_structures/{analysis}_{paramset}_fasta_structure_mapping.tsv"
        numbering_tsv = f"data/protein_structures/pdb_files/{analysis}_{paramset}_pdb_numbering_mapping.tsv"
        sentinel_file = f"data/protein_structures/pdb_files/{analysis}_{paramset}_pdb_numbering_complete.sentinel"
    
    # Create the consolidated mapping
    results = create_pdb_numbering_mapping(base_dir, structures_dir, mapping_file, numbering_tsv, analysis, paramset)
    
    # Create sentinel file
    with open(sentinel_file, 'w') as f:
        f.write(f"PDB numbering extraction completed for {analysis}_{paramset}\n")
        f.write(f"Genes processed: {results.get('genes_found', 0)}\n")
        f.write(f"PDB files processed: {results.get('processed', 0)}\n")
        f.write(f"Total structures: {results.get('total_structures', 0)}\n")
        f.write(f"Output TSV: {numbering_tsv}\n")
    
    print(f"PDB numbering extraction completed. Sentinel: {sentinel_file}")

if __name__ == "__main__":
    main()