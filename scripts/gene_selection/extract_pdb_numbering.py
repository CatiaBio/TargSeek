#!/usr/bin/env python3

"""
Script to extract PDB amino acid numbering from 3D structure files
Creates a consolidated TSV file mapping genes to PDB structures with numbering information
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
        list: List of standard chain IDs (e.g., ['AI', 'AJ', 'AK'] from multi-chain or ['A'] from single chain)
    """
    try:
        with open(fasta_file_path, 'r') as f:
            header_line = f.readline().strip()
            
        # Handle multi-chain case: "Chains AI[auth B1], AJ[auth B2], AK[auth B7], ..." or "Chains A, E[auth F]" or "Chains A, B"
        chains_match = re.search(r'Chains\s+([^|]+)', header_line)
        if chains_match:
            chains_text = chains_match.group(1)
            # Extract all chain IDs (both with and without [auth ...] notation)
            # This regex captures: A, E[auth F], AI[auth B1], B, etc.
            chain_matches = re.findall(r'([A-Z]{1,2})(?:\[auth\s+[^]]+\])?', chains_text)
            if chain_matches:
                return chain_matches
        
        # Handle single chain case: "Chain A[auth Y]" or just "Chain A"
        single_chain_match = re.search(r'Chain\s+([A-Z]{1,2})(?:\[auth\s+[^]]+\])?', header_line)
        if single_chain_match:
            return [single_chain_match.group(1)]
            
    except Exception as e:
        print(f"Error reading FASTA header from {fasta_file_path}: {e}")
    
    return []

def extract_structure_chain_info(structure_file_path):
    """
    Extract chain information and residue numbering from a structure file (PDB or mmCIF)
    
    Args:
        structure_file_path: Path to structure file (.pdb, .pdb.gz, .cif, .cif.gz)
    
    Returns:
        dict: Chain information with residue ranges
    """
    # Determine file type and delegate to appropriate function
    if '.cif' in str(structure_file_path):
        return extract_mmcif_chain_info(structure_file_path)
    else:
        return extract_pdb_chain_info(structure_file_path)

def extract_pdb_chain_info(pdb_file_path):
    """
    Extract chain information and residue numbering from a PDB file using ATOM records
    
    Args:
        pdb_file_path: Path to PDB file (.pdb or .pdb.gz)
    
    Returns:
        dict: Chain information with residue ranges
    """
    chains = {}
    
    try:
        import subprocess
        
        # Use the approach similar to: grep "^ATOM" structure.pdb | awk '{print $5, $6, $4}' | uniq
        if pdb_file_path.endswith('.gz'):
            # For compressed files, use zcat
            cmd = f"zcat '{pdb_file_path}' | grep '^ATOM' | awk '{{print $5, $6, $4}}' | uniq"
        else:
            # For uncompressed files, use grep directly
            cmd = f"grep '^ATOM' '{pdb_file_path}' | awk '{{print $5, $6, $4}}' | uniq"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running command for {pdb_file_path}: {result.stderr}")
            return {}
        
        # Parse the output: Chain ResidueNumber ResidueName
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                parts = line.strip().split()
                if len(parts) == 3:
                    chain, res_num_str, res_name = parts
                    
                    # Handle empty chain (use 'A' as default)
                    if not chain or chain == ' ':
                        chain = 'A'
                    
                    try:
                        res_num = int(res_num_str)
                        
                        if chain not in chains:
                            chains[chain] = {
                                'residues': [],
                                'min_res': res_num,
                                'max_res': res_num,
                                'count': 0
                            }
                        
                        # Track unique residues
                        if res_num not in chains[chain]['residues']:
                            chains[chain]['residues'].append(res_num)
                            chains[chain]['min_res'] = min(chains[chain]['min_res'], res_num)
                            chains[chain]['max_res'] = max(chains[chain]['max_res'], res_num)
                            chains[chain]['count'] += 1
                            
                    except (ValueError, IndexError):
                        continue
        
        # Sort residues for each chain
        for chain in chains:
            chains[chain]['residues'].sort()
            
    except Exception as e:
        print(f"Error processing {pdb_file_path}: {e}")
        return {}
    
    return chains

def extract_mmcif_chain_info(mmcif_file_path):
    """
    Extract chain information and residue numbering from an mmCIF file
    
    Args:
        mmcif_file_path: Path to mmCIF file (.cif or .cif.gz)
    
    Returns:
        dict: Chain information with residue ranges
    """
    chains = {}
    
    try:
        # Open file with appropriate method
        if mmcif_file_path.endswith('.gz'):
            import gzip
            file_handle = gzip.open(mmcif_file_path, 'rt')
        else:
            file_handle = open(mmcif_file_path, 'r')
        
        with file_handle as f:
            in_atom_site_loop = False
            atom_site_columns = {}
            
            for line in f:
                line = line.strip()
                
                if not line or line.startswith('#'):
                    continue
                
                # Look for atom_site loop
                if line == 'loop_':
                    in_atom_site_loop = False
                    atom_site_columns = {}
                    continue
                elif line.startswith('_atom_site.'):
                    in_atom_site_loop = True
                    # Map column names to indices
                    col_name = line[11:]  # Remove '_atom_site.'
                    atom_site_columns[col_name] = len(atom_site_columns)
                    continue
                elif in_atom_site_loop and not line.startswith('_'):
                    # This is an atom record
                    fields = line.split()
                    if len(fields) >= len(atom_site_columns):
                        try:
                            # Extract relevant fields
                            atom_data = {}
                            for col_name, col_index in atom_site_columns.items():
                                if col_index < len(fields):
                                    atom_data[col_name] = fields[col_index]
                            
                            # Only process ATOM records (skip HETATM)
                            if atom_data.get('group_PDB', 'ATOM') == 'ATOM':
                                chain_id = atom_data.get('label_asym_id', 'A')[:1]
                                res_seq = atom_data.get('label_seq_id', '1')
                                
                                try:
                                    res_num = int(res_seq)
                                    
                                    if chain_id not in chains:
                                        chains[chain_id] = {
                                            'residues': [],
                                            'min_res': res_num,
                                            'max_res': res_num,
                                            'count': 0
                                        }
                                    
                                    # Track unique residues
                                    if res_num not in chains[chain_id]['residues']:
                                        chains[chain_id]['residues'].append(res_num)
                                        chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                                        chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                                        chains[chain_id]['count'] += 1
                                        
                                except (ValueError, KeyError):
                                    continue
                                    
                        except (IndexError, ValueError):
                            continue
                    else:
                        # End of atom_site loop
                        in_atom_site_loop = False
                elif not line.startswith('_') and in_atom_site_loop:
                    in_atom_site_loop = False
        
        # Sort residues for each chain
        for chain in chains:
            chains[chain]['residues'].sort()
            
    except Exception as e:
        print(f"Error processing {mmcif_file_path}: {e}")
        return {}
    
    return chains

def create_pdb_numbering_mapping(base_dir, output_tsv, analysis, paramset):
    """
    Create consolidated TSV file mapping genes to PDB structures with numbering information
    
    Args:
        base_dir: Base directory containing gene folders with PDB files
        output_tsv: Output TSV file path
        analysis: Analysis identifier
        paramset: Parameter set identifier
    
    Returns:
        dict: Summary of processing results
    """
    base_path = Path(base_dir)
    
    if not base_path.exists():
        print(f"Error: Base directory {base_dir} does not exist")
        return {}
    
    results = {
        'processed': 0,
        'failed': 0,
        'total_structures': 0,
        'genes_found': 0
    }
    
    # Collect data for TSV
    tsv_data = []
    # Track seen FASTA paths to avoid duplicates
    seen_fasta_paths = set()
    
    # Find all gene directories
    gene_dirs = [d for d in base_path.iterdir() if d.is_dir()]
    
    print(f"Found {len(gene_dirs)} gene directories to process")
    print("=" * 50)
    
    for gene_dir in gene_dirs:
        gene_name = gene_dir.name
        
        # Find structure files in this gene directory (both PDB and mmCIF)
        structure_files = (list(gene_dir.glob("*.pdb")) + 
                          list(gene_dir.glob("*.pdb.gz")) + 
                          list(gene_dir.glob("*.cif")) + 
                          list(gene_dir.glob("*.cif.gz")))
        
        if not structure_files:
            continue
            
        results['genes_found'] += 1
        gene_has_structures = False
        
        for structure_file in structure_files:
            # Extract structure ID from filename
            if structure_file.suffix == '.gz':
                if '.pdb.gz' in structure_file.name:
                    structure_id = structure_file.name.replace('.pdb.gz', '')
                elif '.cif.gz' in structure_file.name:
                    structure_id = structure_file.name.replace('.cif.gz', '')
                else:
                    structure_id = structure_file.stem
            else:
                structure_id = structure_file.stem
            
            print(f"Processing {gene_name}/{structure_id} ({'mmCIF' if '.cif' in structure_file.name else 'PDB'})")
            
            # Extract chain information using the unified function
            chain_info = extract_structure_chain_info(str(structure_file))
            
            if chain_info:
                results['processed'] += 1
                gene_has_structures = True
                
                # Find corresponding FASTA file(s)
                fasta_files = list(gene_dir.glob(f"{structure_id}*.fasta"))
                
                if not fasta_files:
                    # Create a generic entry without FASTA
                    for chain_id in chain_info.keys():
                        tsv_data.append({
                            'Gene': gene_name,
                            'Structure_ID': structure_id,
                            'Chain': chain_id,
                            'Start': chain_info[chain_id]['min_res'],
                            'End': chain_info[chain_id]['max_res'],
                            'Range': f"{chain_info[chain_id]['min_res']}-{chain_info[chain_id]['max_res']}",
                            'FASTA_Path': f"data/protein_structures/{gene_name}/{structure_id}.fasta"
                        })
                else:
                    # Match chains to FASTA files
                    for fasta_file in fasta_files:
                        # Extract standard chain IDs from FASTA header
                        standard_chain_ids = extract_standard_chains_from_fasta(fasta_file)
                        
                        if not standard_chain_ids:
                            # Fallback: try filename-based extraction for backwards compatibility
                            fasta_name = fasta_file.stem
                            if '_chain_' in fasta_name.lower():
                                chain_part = fasta_name.split('_chain_')[-1]
                                if chain_part.isdigit():
                                    chain_num = int(chain_part)
                                    if chain_num >= 1:
                                        standard_chain_ids = [chr(ord('A') + chain_num - 1)]
                                    else:
                                        standard_chain_ids = ['A']
                                else:
                                    standard_chain_ids = ['A']
                            else:
                                standard_chain_ids = ['A']
                        
                        # Get sorted author chain IDs from structure
                        author_chain_ids = sorted(chain_info.keys())
                        
                        # For multi-chain FASTA files, only keep the first chain to avoid duplicates
                        # Create FASTA path for deduplication check
                        fasta_path = f"data/protein_structures/{gene_name}/{fasta_file.name}"
                        
                        if fasta_path not in seen_fasta_paths:
                            seen_fasta_paths.add(fasta_path)
                            
                            # Use the first standard chain ID
                            standard_chain_id = standard_chain_ids[0] if standard_chain_ids else 'A'
                            
                            # Map to first author chain by position
                            fasta_name = fasta_file.stem
                            if '_chain_' in fasta_name.lower():
                                # Extract chain number from filename (e.g., "1KMI_chain_1" -> chain 1 -> author_chain_ids[0])
                                chain_part = fasta_name.split('_chain_')[-1]
                                if chain_part.isdigit():
                                    chain_num = int(chain_part)
                                    if chain_num <= len(author_chain_ids):
                                        author_chain_id = author_chain_ids[chain_num - 1]
                                    else:
                                        author_chain_id = author_chain_ids[0] if author_chain_ids else 'A'
                                else:
                                    author_chain_id = author_chain_ids[0] if author_chain_ids else 'A'
                            else:
                                # For multi-chain FASTA without filename indicators, use first author chain
                                author_chain_id = author_chain_ids[0] if author_chain_ids else 'A'
                            
                            if author_chain_id in chain_info:
                                tsv_data.append({
                                    'Gene': gene_name,
                                    'Structure_ID': structure_id,
                                    'Chain': standard_chain_id,  # Use first standard chain ID from FASTA header
                                    'Start': chain_info[author_chain_id]['min_res'],
                                    'End': chain_info[author_chain_id]['max_res'],
                                    'Range': f"{chain_info[author_chain_id]['min_res']}-{chain_info[author_chain_id]['max_res']}",
                                    'FASTA_Path': fasta_path
                                })
            else:
                results['failed'] += 1
                print(f"  Failed to extract chain info")
        
        if gene_has_structures:
            results['total_structures'] += len(structure_files)
    
    # Create DataFrame and save to TSV
    if tsv_data:
        df = pd.DataFrame(tsv_data)
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"\nCreated mapping TSV: {output_tsv}")
        print(f"Total entries: {len(tsv_data)}")
    else:
        print(f"No data to write to TSV file")
        # Create empty TSV with headers
        df = pd.DataFrame(columns=['Gene', 'Structure_ID', 'Chain', 'Start', 'End', 'Range', 'FASTA_Path'])
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
        metadata_tsv = snakemake.input.metadata_tsv
        numbering_tsv = snakemake.output.numbering_tsv
        sentinel_file = snakemake.output.numbering_sentinel
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        base_dir = snakemake.params.base_dir
        
        print(f"PDB numbering extraction for {analysis}_{paramset}")
        print(f"Base directory: {base_dir}")
        print(f"Output TSV: {numbering_tsv}")
        
    except NameError:
        # Test mode - command line arguments
        if len(sys.argv) != 4:
            print("Usage: python extract_pdb_numbering.py <analysis> <paramset> <base_dir>")
            print("Example: python extract_pdb_numbering.py analysis1 params1 data/protein_structures")
            sys.exit(1)
        
        analysis = sys.argv[1]
        paramset = sys.argv[2]
        base_dir = sys.argv[3]
        
        numbering_tsv = f"data/protein_structures/{analysis}_{paramset}_pdb_numbering_mapping.tsv"
        sentinel_file = f"data/protein_structures/{analysis}_{paramset}_pdb_numbering_complete.sentinel"
    
    # Create the consolidated mapping
    results = create_pdb_numbering_mapping(base_dir, numbering_tsv, analysis, paramset)
    
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