#!/usr/bin/env python3
"""
Extract chain information (chain ID, start residue, end residue) from PDB structures.
Handles both .pdb.gz and .cif.gz files.

Output format: TSV files with columns: chain, start, end
Example: A 24 808
"""

import os
import gzip
import glob
import sys
from pathlib import Path
import re
from collections import defaultdict


def extract_chain_info_pdb(file_path):
    """Extract chain information from PDB format files."""
    chain_ranges = defaultdict(lambda: {'start': None, 'end': None})
    
    try:
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                lines = f.readlines()
        else:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        
        for line in lines:
            if line.startswith('ATOM'):
                # PDB format: ATOM records contain chain and residue info
                # Column positions (1-indexed): chain=22, resSeq=23-26
                chain = line[21]  # Chain identifier (0-indexed: position 21)
                try:
                    res_num = int(line[22:26].strip())  # Residue number (0-indexed: positions 22-25)
                    
                    if chain_ranges[chain]['start'] is None or res_num < chain_ranges[chain]['start']:
                        chain_ranges[chain]['start'] = res_num
                    if chain_ranges[chain]['end'] is None or res_num > chain_ranges[chain]['end']:
                        chain_ranges[chain]['end'] = res_num
                        
                except (ValueError, IndexError):
                    continue
                    
    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)
        return {}
    
    return dict(chain_ranges)


def extract_chain_info_cif(file_path):
    """Extract chain information from CIF format files."""
    chain_ranges = defaultdict(lambda: {'start': None, 'end': None})
    
    try:
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                content = f.read()
        else:
            with open(file_path, 'r') as f:
                content = f.read()
        
        # Look for atom_site loop in CIF format
        in_atom_site = False
        headers = []
        
        for line in content.split('\n'):
            line = line.strip()
            
            if line.startswith('_atom_site.'):
                in_atom_site = True
                headers.append(line[11:])  # Remove '_atom_site.' prefix
                continue
            elif in_atom_site and line.startswith('_'):
                # End of atom_site loop
                break
            elif in_atom_site and line and not line.startswith('#'):
                # Process atom site data
                parts = line.split()
                if len(parts) >= len(headers):
                    try:
                        # Find indices for chain and residue number
                        chain_idx = headers.index('auth_asym_id') if 'auth_asym_id' in headers else headers.index('label_asym_id')
                        res_idx = headers.index('auth_seq_id') if 'auth_seq_id' in headers else headers.index('label_seq_id')
                        
                        chain = parts[chain_idx]
                        res_num = int(parts[res_idx])
                        
                        if chain_ranges[chain]['start'] is None or res_num < chain_ranges[chain]['start']:
                            chain_ranges[chain]['start'] = res_num
                        if chain_ranges[chain]['end'] is None or res_num > chain_ranges[chain]['end']:
                            chain_ranges[chain]['end'] = res_num
                            
                    except (ValueError, IndexError):
                        continue
                        
    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)
        return {}
    
    return dict(chain_ranges)


def process_structure_file(file_path):
    """Process a single structure file and extract chain information."""
    file_str = str(file_path)
    if file_str.endswith('.pdb.gz') or file_str.endswith('.pdb'):
        return extract_chain_info_pdb(file_str)
    elif file_str.endswith('.cif.gz') or file_str.endswith('.cif'):
        return extract_chain_info_cif(file_str)
    else:
        print(f"Unsupported file format: {file_path}", file=sys.stderr)
        return {}


def save_chain_info(chain_ranges, output_file):
    """Save chain information to TSV file."""
    with open(output_file, 'w') as f:
        f.write("chain\tstart\tend\n")
        for chain, ranges in sorted(chain_ranges.items()):
            if ranges['start'] is not None and ranges['end'] is not None:
                f.write(f"{chain}\t{ranges['start']}\t{ranges['end']}\n")


def main():
    """Main function to process all structure files."""
    # Define input directory
    pdb_dir = Path("data/protein_structures/pdb_files")
    
    if not pdb_dir.exists():
        print(f"Directory not found: {pdb_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Find all PDB and CIF files
    pdb_files = list(pdb_dir.glob("*.pdb.gz")) + list(pdb_dir.glob("*.pdb"))
    cif_files = list(pdb_dir.glob("*.cif.gz")) + list(pdb_dir.glob("*.cif"))
    
    all_files = pdb_files + cif_files
    
    if not all_files:
        print("No structure files found", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {len(all_files)} structure files...")
    
    processed_count = 0
    for file_path in all_files:
        # Extract structure ID from filename
        structure_id = file_path.stem
        if structure_id.endswith('.pdb') or structure_id.endswith('.cif'):
            structure_id = structure_id.rsplit('.', 1)[0]
        
        # Process the file
        chain_ranges = process_structure_file(file_path)
        
        if chain_ranges:
            # Save results
            output_file = pdb_dir / f"{structure_id}_chains.tsv"
            save_chain_info(chain_ranges, output_file)
            print(f"Processed {structure_id}: {len(chain_ranges)} chains -> {output_file}")
            processed_count += 1
        else:
            print(f"No chain information found for {structure_id}", file=sys.stderr)
    
    print(f"Successfully processed {processed_count}/{len(all_files)} files")


if __name__ == "__main__":
    main()