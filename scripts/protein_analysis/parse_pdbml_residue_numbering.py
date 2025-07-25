#!/usr/bin/env python3
"""
Parse PDBML/XML Residue Numbering
=================================

This script demonstrates how to parse PDBML/XML files to extract correct residue 
numbering information for accurate epitope mapping. PDBML files contain the 
auth_seq_id which represents the original sequence numbering.

Usage:
    python parse_pdbml_residue_numbering.py <pdbml_file> [--chain CHAIN]
"""

import xml.etree.ElementTree as ET
import gzip
import argparse
import pandas as pd
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_pdbml_residues(pdbml_file, target_chain=None):
    """
    Parse PDBML/XML file to extract residue numbering information
    
    Args:
        pdbml_file (Path): Path to PDBML/XML file (can be gzipped)
        target_chain (str): Specific chain to extract (default: all chains)
    
    Returns:
        dict: Dictionary mapping chain IDs to residue information
    """
    
    try:
        # Handle gzipped files
        if str(pdbml_file).endswith('.gz'):
            with gzip.open(pdbml_file, 'rt') as f:
                tree = ET.parse(f)
        else:
            tree = ET.parse(pdbml_file)
        
        root = tree.getroot()
        
        # Find the namespace
        namespace = {'pdbx': 'http://pdbml.pdb.org/schema/pdbx-v50.xsd'}
        if root.tag.startswith('{'):
            # Extract namespace from root tag
            ns_end = root.tag.find('}')
            if ns_end != -1:
                ns_uri = root.tag[1:ns_end]
                namespace = {'pdbx': ns_uri}
        
        # Find atom_site category which contains residue information
        atom_sites = root.findall('.//pdbx:atom_siteCategory/pdbx:atom_site', namespace)
        
        if not atom_sites:
            # Try without namespace
            atom_sites = root.findall('.//atom_site')
        
        chains_data = {}
        
        for atom_site in atom_sites:
            # Extract relevant information
            chain_id = None
            auth_seq_id = None
            auth_comp_id = None
            auth_atom_id = None
            
            # Try with namespace first
            for child in atom_site:
                tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                
                if tag == 'auth_asym_id':
                    chain_id = child.text
                elif tag == 'auth_seq_id':
                    auth_seq_id = child.text
                elif tag == 'auth_comp_id':
                    auth_comp_id = child.text
                elif tag == 'auth_atom_id':
                    auth_atom_id = child.text
            
            # Skip if we don't have essential information or if it's not a CA atom
            if not all([chain_id, auth_seq_id, auth_comp_id]) or auth_atom_id != 'CA':
                continue
            
            # Filter by target chain if specified
            if target_chain and chain_id != target_chain:
                continue
            
            # Skip non-amino acid residues
            amino_acids = {
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
            }
            
            if auth_comp_id not in amino_acids:
                continue
            
            # Initialize chain data if not exists
            if chain_id not in chains_data:
                chains_data[chain_id] = []
            
            # Convert auth_seq_id to integer for proper sorting
            try:
                seq_id_int = int(auth_seq_id)
            except ValueError:
                # Handle insertion codes or non-numeric sequences
                seq_id_int = auth_seq_id
            
            # Store residue information
            residue_info = {
                'auth_seq_id': auth_seq_id,
                'auth_seq_id_int': seq_id_int,
                'auth_comp_id': auth_comp_id,
                'chain_id': chain_id
            }
            
            # Avoid duplicates (multiple atoms per residue)
            if not any(r['auth_seq_id'] == auth_seq_id for r in chains_data[chain_id]):
                chains_data[chain_id].append(residue_info)
        
        # Sort residues by sequence ID within each chain
        for chain_id in chains_data:
            chains_data[chain_id].sort(key=lambda x: x['auth_seq_id_int'] if isinstance(x['auth_seq_id_int'], int) else float('inf'))
        
        return chains_data
        
    except Exception as e:
        logging.error(f"Error parsing PDBML file {pdbml_file}: {e}")
        return {}

def convert_aa_code(three_letter):
    """Convert 3-letter amino acid code to 1-letter"""
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_dict.get(three_letter, 'X')

def create_residue_mapping_table(pdbml_file, output_file=None, target_chain=None):
    """
    Create a detailed residue mapping table from PDBML file
    
    Args:
        pdbml_file (Path): Path to PDBML file
        output_file (Path): Output CSV file path
        target_chain (str): Specific chain to process
    """
    
    # Parse PDBML file
    chains_data = parse_pdbml_residues(pdbml_file, target_chain)
    
    if not chains_data:
        logging.error(f"No residue data found in {pdbml_file}")
        return
    
    # Create mapping table
    mapping_data = []
    
    for chain_id, residues in chains_data.items():
        for i, residue_info in enumerate(residues):
            mapping_data.append({
                'Chain': chain_id,
                'Sequential_Position': i + 1,  # 1-based sequential numbering
                'Auth_Seq_ID': residue_info['auth_seq_id'],
                'Auth_Comp_ID': residue_info['auth_comp_id'],
                'AA_Code': convert_aa_code(residue_info['auth_comp_id'])
            })
    
    # Create DataFrame
    df = pd.DataFrame(mapping_data)
    
    if output_file:
        df.to_csv(output_file, index=False)
        logging.info(f"Residue mapping table saved to: {output_file}")
    
    # Display summary
    pdb_id = pdbml_file.stem.replace('.xml', '').replace('.pdb', '')
    print(f"\n{'='*80}")
    print(f"PDBML Residue Analysis: {pdb_id}")
    print(f"{'='*80}")
    print(f"File: {pdbml_file}")
    print(f"Total chains: {len(chains_data)}")
    
    for chain_id, residues in chains_data.items():
        print(f"\nChain {chain_id}:")
        print(f"  Total residues: {len(residues)}")
        if residues:
            first_residue = residues[0]
            last_residue = residues[-1]
            print(f"  First residue: {first_residue['auth_seq_id']} ({convert_aa_code(first_residue['auth_comp_id'])})")
            print(f"  Last residue:  {last_residue['auth_seq_id']} ({convert_aa_code(last_residue['auth_comp_id'])})")
            
            # Build sequence
            sequence = ''.join([convert_aa_code(r['auth_comp_id']) for r in residues])
            print(f"  Sequence (first 50): {sequence[:50]}...")
            if len(sequence) > 50:
                print(f"  Sequence (last 50):  ...{sequence[-50:]}")
    
    if output_file:
        print(f"\nDetailed mapping saved to: {output_file}")
    
    # Show first few rows of mapping
    if not df.empty:
        print(f"\nFirst 10 residues from mapping table:")
        print(df.head(10).to_string(index=False))
        if len(df) > 10:
            print(f"\nLast 10 residues from mapping table:")
            print(df.tail(10).to_string(index=False))
    
    return df

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Parse PDBML/XML residue numbering information')
    parser.add_argument('pdbml_file', help='Path to PDBML/XML file (can be gzipped)')
    parser.add_argument('--chain', help='Specific chain to analyze (default: all chains)')
    parser.add_argument('--output', help='Output CSV file for detailed mapping table')
    
    args = parser.parse_args()
    
    pdbml_file = Path(args.pdbml_file)
    if not pdbml_file.exists():
        logging.error(f"PDBML file not found: {pdbml_file}")
        return
    
    # Generate default output filename if not provided
    output_file = None
    if args.output:
        output_file = Path(args.output)
    else:
        pdb_id = pdbml_file.stem.replace('.xml', '').replace('.pdb', '')
        output_file = pdbml_file.parent / f"{pdb_id}_pdbml_residue_mapping.csv"
    
    # Create residue mapping table
    create_residue_mapping_table(pdbml_file, output_file, args.chain)

if __name__ == "__main__":
    main()