#!/usr/bin/env python3

"""
Script to extract PDB amino acid numbering from 3D structure files
Processes both .pdb and .pdb.gz files to create numbering files for each structure
"""

import os
import sys
import gzip
from pathlib import Path

def extract_pdb_numbering(pdb_file_path, output_file_path):
    """
    Extract amino acid numbering from a PDB file
    
    Args:
        pdb_file_path: Path to PDB file (.pdb or .pdb.gz)
        output_file_path: Path for output numbering file
    
    Returns:
        int: Number of unique residues found
    """
    residues = []
    
    try:
        # Handle both compressed and uncompressed files
        if pdb_file_path.endswith('.gz'):
            with gzip.open(pdb_file_path, 'rt') as f:
                lines = f.readlines()
        else:
            with open(pdb_file_path, 'r') as f:
                lines = f.readlines()
        
        # Extract ATOM lines and parse chain, residue number, residue type
        seen_residues = set()
        
        for line in lines:
            if line.startswith('ATOM'):
                # Parse PDB ATOM record
                # Columns: chain_id (22), res_seq_num (23-26), res_name (18-20)
                chain_id = line[21].strip()
                res_seq_num = line[22:26].strip()
                res_name = line[17:20].strip()
                
                # Create unique identifier for residue
                residue_key = f"{chain_id}_{res_seq_num}_{res_name}"
                
                if residue_key not in seen_residues:
                    residues.append((chain_id, res_seq_num, res_name))
                    seen_residues.add(residue_key)
        
        # Write output file
        with open(output_file_path, 'w') as f:
            for chain_id, res_seq_num, res_name in residues:
                f.write(f"{chain_id} {res_seq_num} {res_name}\n")
        
        print(f"Processed {pdb_file_path} -> {output_file_path} ({len(residues)} residues)")
        return len(residues)
        
    except Exception as e:
        print(f"Error processing {pdb_file_path}: {e}")
        return 0

def process_all_structures(base_dir):
    """
    Process all PDB files in the protein structures directory
    
    Args:
        base_dir: Base directory containing gene folders with PDB files
    
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
        'total_residues': 0,
        'files_created': []
    }
    
    # Find all PDB files recursively
    pdb_files = list(base_path.rglob("*.pdb")) + list(base_path.rglob("*.pdb.gz"))
    
    print(f"Found {len(pdb_files)} PDB files to process")
    print("=" * 50)
    
    for pdb_file in pdb_files:
        # Create output filename
        if pdb_file.suffix == '.gz':
            # Remove .pdb.gz and add _numbering.txt
            output_file = pdb_file.with_suffix('').with_suffix('') / f"{pdb_file.stem.replace('.pdb', '')}_numbering.txt"
            output_file = pdb_file.parent / f"{pdb_file.name.replace('.pdb.gz', '')}_numbering.txt"
        else:
            # Remove .pdb and add _numbering.txt
            output_file = pdb_file.with_suffix('') / f"{pdb_file.stem}_numbering.txt"
            output_file = pdb_file.parent / f"{pdb_file.stem}_numbering.txt"
        
        # Extract numbering
        residue_count = extract_pdb_numbering(str(pdb_file), str(output_file))
        
        if residue_count > 0:
            results['processed'] += 1
            results['total_residues'] += residue_count
            results['files_created'].append(str(output_file))
        else:
            results['failed'] += 1
    
    print("=" * 50)
    print(f"Processing complete!")
    print(f"Successfully processed: {results['processed']} files")
    print(f"Failed: {results['failed']} files")
    print(f"Total residues extracted: {results['total_residues']}")
    print(f"Numbering files created: {len(results['files_created'])}")
    
    return results

def process_specific_structures(analysis, paramset, base_dir):
    """
    Process structures for specific analysis and parameter set
    
    Args:
        analysis: Analysis identifier
        paramset: Parameter set identifier  
        base_dir: Base directory for protein structures
    
    Returns:
        str: Path to completion sentinel file
    """
    print(f"Processing PDB numbering extraction for {analysis}_{paramset}")
    
    # Process all structures in the base directory
    results = process_all_structures(base_dir)
    
    # Create sentinel file to mark completion
    sentinel_file = os.path.join(base_dir, f"{analysis}_{paramset}_pdb_numbering_complete.sentinel")
    
    with open(sentinel_file, 'w') as f:
        f.write(f"PDB numbering extraction completed for {analysis}_{paramset}\n")
        f.write(f"Files processed: {results['processed']}\n")
        f.write(f"Files failed: {results['failed']}\n")
        f.write(f"Total residues: {results['total_residues']}\n")
        f.write(f"Numbering files created: {len(results['files_created'])}\n")
        f.write("\nGenerated files:\n")
        for file_path in results['files_created']:
            f.write(f"  {file_path}\n")
    
    print(f"Sentinel file created: {sentinel_file}")
    return sentinel_file

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_pdb_numbering.py <analysis> <paramset> <base_dir>")
        print("Example: python extract_pdb_numbering.py analysis1 params1 data/protein_structures")
        sys.exit(1)
    
    analysis = sys.argv[1]
    paramset = sys.argv[2]
    base_dir = sys.argv[3]
    
    sentinel_file = process_specific_structures(analysis, paramset, base_dir)
    print(f"PDB numbering extraction completed. Sentinel: {sentinel_file}")

if __name__ == "__main__":
    main()