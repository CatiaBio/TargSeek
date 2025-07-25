#!/usr/bin/env python3
"""
Enhanced FASTA to PDB Residue Mapping with PDBML Support
========================================================

This enhanced version uses PDBML/XML files when available to get the most
accurate residue mapping, as PDBML contains the auth_seq_id which represents
the original sequence numbering used by the authors.

Usage:
    python map_fasta_to_pdb_residues_enhanced.py <bepipred_dir> [--gene GENE] [--prefer-pdbml]
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from Bio.PDB import PDBParser
import argparse
import sys
import gzip
import xml.etree.ElementTree as ET

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class EnhancedPDBMapper:
    """Enhanced PDB residue mapper with PDBML support"""
    
    def __init__(self, prefer_pdbml=True):
        self.prefer_pdbml = prefer_pdbml
        self.parser = PDBParser(QUIET=True)
    
    def parse_pdbml_residues(self, pdbml_file, target_chain='A'):
        """
        Parse PDBML/XML file to extract accurate residue numbering
        
        Args:
            pdbml_file (Path): Path to PDBML/XML file
            target_chain (str): Chain to extract
        
        Returns:
            list: List of residue information with auth_seq_id
        """
        
        try:
            # Handle gzipped files
            if str(pdbml_file).endswith('.gz'):
                with gzip.open(pdbml_file, 'rt') as f:
                    tree = ET.parse(f)
            else:
                tree = ET.parse(pdbml_file)
            
            root = tree.getroot()
            
            # Find namespace
            namespace = {'pdbx': 'http://pdbml.pdb.org/schema/pdbx-v50.xsd'}
            if root.tag.startswith('{'):
                ns_end = root.tag.find('}')
                if ns_end != -1:
                    ns_uri = root.tag[1:ns_end]
                    namespace = {'pdbx': ns_uri}
            
            # Find atom_site category
            atom_sites = root.findall('.//pdbx:atom_siteCategory/pdbx:atom_site', namespace)
            
            if not atom_sites:
                # Try without namespace
                atom_sites = root.findall('.//atom_site')
            
            residues = []
            
            for atom_site in atom_sites:
                chain_id = None
                auth_seq_id = None
                auth_comp_id = None
                auth_atom_id = None
                
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
                
                # Only process CA atoms from target chain
                if (chain_id == target_chain and auth_atom_id == 'CA' and 
                    auth_seq_id and auth_comp_id):
                    
                    # Convert to 1-letter amino acid code
                    aa_code = self.three_to_one(auth_comp_id)
                    if aa_code:
                        try:
                            seq_id_int = int(auth_seq_id)
                        except ValueError:
                            seq_id_int = auth_seq_id
                        
                        residue_info = {
                            'auth_seq_id': auth_seq_id,
                            'auth_seq_id_int': seq_id_int,
                            'aa_code': aa_code,
                            'comp_id': auth_comp_id
                        }
                        
                        # Avoid duplicates
                        if not any(r['auth_seq_id'] == auth_seq_id for r in residues):
                            residues.append(residue_info)
            
            # Sort by sequence ID
            residues.sort(key=lambda x: x['auth_seq_id_int'] if isinstance(x['auth_seq_id_int'], int) else float('inf'))
            
            logging.info(f"Extracted {len(residues)} residues from PDBML chain {target_chain}")
            return residues
            
        except Exception as e:
            logging.warning(f"Failed to parse PDBML file {pdbml_file}: {e}")
            return []
    
    def extract_pdb_sequence_legacy(self, pdb_file, chain_id='A'):
        """Extract sequence from legacy PDB format (fallback method)"""
        try:
            if str(pdb_file).endswith('.gz'):
                with gzip.open(pdb_file, 'rt') as f:
                    structure = self.parser.get_structure('protein', f)
            else:
                structure = self.parser.get_structure('protein', pdb_file)
            
            model = structure[0]
            if chain_id not in model:
                available_chains = [c.id for c in model]
                logging.warning(f"Chain {chain_id} not found in {pdb_file}. Available: {available_chains}")
                if available_chains:
                    chain_id = available_chains[0]
                    logging.info(f"Using chain {chain_id} instead")
                else:
                    return "", []
            
            chain = model[chain_id]
            sequence = ""
            residue_numbers = []
            
            for residue in chain:
                if residue.id[0] == ' ':  # Standard amino acid
                    resnum = residue.id[1]
                    resname = residue.resname
                    aa_code = self.three_to_one(resname)
                    if aa_code:
                        sequence += aa_code
                        residue_numbers.append(resnum)
            
            return sequence, residue_numbers
            
        except Exception as e:
            logging.error(f"Error parsing PDB file {pdb_file}: {e}")
            return "", []
    
    def three_to_one(self, three_letter):
        """Convert 3-letter amino acid code to 1-letter"""
        aa_dict = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        return aa_dict.get(three_letter, None)
    
    def extract_bepipred_sequence(self, csv_file):
        """Extract sequence from BepiPred CSV file"""
        try:
            df = pd.read_csv(csv_file)
            residues = []
            for _, row in df.iterrows():
                residues.append(row['Residue'])
            sequence = ''.join(residues)
            return sequence
        except Exception as e:
            logging.error(f"Error reading BepiPred CSV {csv_file}: {e}")
            return ""
    
    def map_bepipred_to_pdb_enhanced(self, bepipred_csv, pdb_structures_dir, gene_name, structure_id):
        """
        Enhanced mapping using PDBML when available, fallback to legacy PDB
        
        Args:
            bepipred_csv (Path): BepiPred CSV file
            pdb_structures_dir (Path): Directory containing PDB structures
            gene_name (str): Gene name
            structure_id (str): Structure identifier
        
        Returns:
            dict: Mapping information
        """
        
        # Find structure files
        gene_dir = pdb_structures_dir / gene_name
        pdbml_file = gene_dir / f"{structure_id}.xml.gz"
        pdb_file = gene_dir / f"{structure_id}.pdb.gz"
        
        # Extract BepiPred sequence
        bepipred_seq = self.extract_bepipred_sequence(bepipred_csv)
        if not bepipred_seq:
            return {"error": "Could not extract BepiPred sequence"}
        
        # Try PDBML first if preferred and available
        pdb_residues = []
        mapping_method = "unknown"
        
        if self.prefer_pdbml and pdbml_file.exists():
            logging.info(f"Using PDBML file for accurate mapping: {pdbml_file}")
            pdb_residues = self.parse_pdbml_residues(pdbml_file, 'A')
            if pdb_residues:
                mapping_method = "pdbml_auth_seq_id"
                # Extract sequence from PDBML residues
                pdb_seq = ''.join([r['aa_code'] for r in pdb_residues])
                # Extract auth_seq_ids
                pdb_numbers = [r['auth_seq_id'] for r in pdb_residues]
            else:
                logging.warning(f"PDBML parsing failed, falling back to legacy PDB")
        
        # Fallback to legacy PDB if PDBML not available or failed
        if not pdb_residues and pdb_file.exists():
            logging.info(f"Using legacy PDB file: {pdb_file}")
            pdb_seq, pdb_numbers = self.extract_pdb_sequence_legacy(pdb_file, 'A')
            if pdb_seq:
                mapping_method = "legacy_pdb_residue_id"
            else:
                return {"error": "Could not extract PDB sequence"}
        
        if not pdb_residues and not pdb_seq:
            return {"error": f"No structure files found for {structure_id}"}
        
        # For PDBML, we already have the data
        if mapping_method == "pdbml_auth_seq_id":
            pass  # pdb_seq and pdb_numbers already set above
        
        # Perform sequence alignment
        logging.info(f"BepiPred sequence length: {len(bepipred_seq)}")
        logging.info(f"PDB sequence length: {len(pdb_seq)}")
        logging.info(f"Mapping method: {mapping_method}")
        
        # Find alignment offset
        bepipred_to_pdb_offset = 0
        alignment_found = False
        
        # Check if PDB sequence is a substring of BepiPred sequence
        if pdb_seq in bepipred_seq:
            bepipred_start_pos = bepipred_seq.find(pdb_seq)
            bepipred_to_pdb_offset = -bepipred_start_pos
            alignment_found = True
            logging.info(f"PDB sequence found at BepiPred position {bepipred_start_pos + 1}")
        
        # Check if BepiPred sequence is a substring of PDB sequence
        elif bepipred_seq in pdb_seq:
            pdb_start_pos = pdb_seq.find(bepipred_seq)
            bepipred_to_pdb_offset = pdb_start_pos
            alignment_found = True
            logging.info(f"BepiPred sequence found at PDB position {pdb_start_pos + 1}")
        
        # Try alignment from the beginning
        elif len(pdb_seq) >= 10 and len(bepipred_seq) >= 10:
            # Check if sequences start similarly
            if pdb_seq[:10] == bepipred_seq[:10]:
                bepipred_to_pdb_offset = 0
                alignment_found = True
                logging.info("Sequences appear to start at the same position")
        
        if not alignment_found:
            logging.warning("Could not find clear sequence alignment")
            # Try to estimate offset from sequence lengths
            if abs(len(pdb_seq) - len(bepipred_seq)) < len(bepipred_seq) * 0.2:
                bepipred_to_pdb_offset = 0
                logging.info("Using zero offset due to similar sequence lengths")
            else:
                return {"error": "Could not align sequences"}
        
        return {
            "bepipred_sequence": bepipred_seq,
            "pdb_sequence": pdb_seq,
            "pdb_residue_numbers": pdb_numbers,
            "bepipred_to_pdb_offset": bepipred_to_pdb_offset,
            "mapping_method": mapping_method,
            "structure_files": {
                "pdbml": str(pdbml_file) if pdbml_file.exists() else None,
                "pdb": str(pdb_file) if pdb_file.exists() else None
            }
        }
    
    def map_bepipred_position_to_pdb(self, bepipred_pos, mapping_info):
        """
        Map BepiPred position to PDB residue number
        
        Args:
            bepipred_pos (int): 1-based BepiPred position
            mapping_info (dict): Mapping information from map_bepipred_to_pdb_enhanced
        
        Returns:
            str or None: PDB residue number or None if out of range
        """
        
        if "error" in mapping_info:
            return None
        
        # Convert to 0-based index
        bepipred_idx = bepipred_pos - 1
        pdb_idx = bepipred_idx + mapping_info["bepipred_to_pdb_offset"]
        
        # Check bounds
        pdb_numbers = mapping_info["pdb_residue_numbers"]
        if 0 <= pdb_idx < len(pdb_numbers):
            return pdb_numbers[pdb_idx]
        else:
            return None

def process_bepipred_epitope_file(epitope_file, mapper, pdb_structures_dir):
    """Process a BepiPred epitope file and add enhanced PDB mapping"""
    
    try:
        # Extract information from filename
        file_parts = epitope_file.stem.split('_')
        gene_name = file_parts[0]
        structure_id = '_'.join(file_parts[1:-2])  # Remove gene name and file suffix
        
        # Find corresponding CSV file
        bepipred_dir = epitope_file.parent.parent
        csv_file = None
        for potential_csv in bepipred_dir.glob(f"{gene_name}_{structure_id}*.csv"):
            csv_file = potential_csv
            break
        
        if not csv_file or not csv_file.exists():
            logging.error(f"Could not find BepiPred CSV file for {epitope_file}")
            return
        
        logging.info(f"Processing {epitope_file.name}")
        
        # Get enhanced mapping
        mapping_info = mapper.map_bepipred_to_pdb_enhanced(
            csv_file, pdb_structures_dir, gene_name, structure_id
        )
        
        if "error" in mapping_info:
            logging.error(f"Mapping failed for {structure_id}: {mapping_info['error']}")
            return
        
        # Read existing epitope file
        with open(epitope_file, 'r') as f:
            lines = f.readlines()
        
        # Find data lines (skip header and comments)
        data_lines = []
        header_line = None
        
        for line in lines:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            elif line.startswith('No.\t'):
                header_line = line
            else:
                data_lines.append(line)
        
        if not header_line or not data_lines:
            logging.warning(f"No epitope data found in {epitope_file}")
            return
        
        # Create enhanced output file
        output_file = epitope_file.parent / f"{epitope_file.stem}_enhanced_pdbml_mapping.tsv"
        
        with open(output_file, 'w') as f:
            # Write header with mapping method info
            f.write(f"# Enhanced Linear Epitope Mapping using {mapping_info['mapping_method'].upper()}\n")
            f.write(f"# Structure: {structure_id}\n")
            f.write(f"# Mapping method: {mapping_info['mapping_method']}\n")
            f.write(f"# BepiPred sequence length: {len(mapping_info['bepipred_sequence'])}\n")
            f.write(f"# PDB sequence length: {len(mapping_info['pdb_sequence'])}\n")
            f.write(f"# Available files: PDB={mapping_info['structure_files']['pdb'] is not None}, ")
            f.write(f"PDBML={mapping_info['structure_files']['pdbml'] is not None}\n")
            
            # Enhanced header
            enhanced_header = header_line + "\tPDB_Start\tPDB_End\tPDB_Residues\tMapping_Method"
            f.write(enhanced_header + "\n")
            
            # Process each epitope
            for line in data_lines:
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                
                try:
                    start_pos = int(parts[2])  # Start position
                    end_pos = int(parts[3])    # End position
                    
                    # Map positions to PDB
                    pdb_start = mapper.map_bepipred_position_to_pdb(start_pos, mapping_info)
                    pdb_end = mapper.map_bepipred_position_to_pdb(end_pos, mapping_info)
                    
                    # Create PDB residue list
                    pdb_residues = []
                    if pdb_start is not None and pdb_end is not None:
                        for pos in range(start_pos, end_pos + 1):
                            pdb_pos = mapper.map_bepipred_position_to_pdb(pos, mapping_info)
                            if pdb_pos is not None:
                                pdb_residues.append(str(pdb_pos))
                    
                    # Format output
                    pdb_start_str = str(pdb_start) if pdb_start is not None else "?"
                    pdb_end_str = str(pdb_end) if pdb_end is not None else "?"
                    pdb_residues_str = ",".join(pdb_residues) if pdb_residues else "?"
                    
                    enhanced_line = f"{line}\t{pdb_start_str}\t{pdb_end_str}\t{pdb_residues_str}\t{mapping_info['mapping_method']}"
                    f.write(enhanced_line + "\n")
                    
                except (ValueError, IndexError) as e:
                    logging.warning(f"Error processing line: {line} - {e}")
                    enhanced_line = f"{line}\t?\t?\t?\t{mapping_info['mapping_method']}"
                    f.write(enhanced_line + "\n")
        
        logging.info(f"Enhanced mapping saved to: {output_file}")
        
    except Exception as e:
        logging.error(f"Error processing {epitope_file}: {e}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Enhanced FASTA to PDB residue mapping with PDBML support')
    parser.add_argument('bepipred_dir', help='Directory containing BepiPred results')
    parser.add_argument('--gene', help='Process specific gene only')
    parser.add_argument('--prefer-pdbml', action='store_true', default=True,
                       help='Prefer PDBML files for accurate mapping (default: True)')
    parser.add_argument('--pdb-structures-dir', default='data/protein_structures',
                       help='Directory containing PDB structure files')
    
    args = parser.parse_args()
    
    bepipred_dir = Path(args.bepipred_dir)
    pdb_structures_dir = Path(args.pdb_structures_dir)
    
    if not bepipred_dir.exists():
        logging.error(f"BepiPred directory not found: {bepipred_dir}")
        sys.exit(1)
    
    if not pdb_structures_dir.exists():
        logging.error(f"PDB structures directory not found: {pdb_structures_dir}")
        sys.exit(1)
    
    # Initialize enhanced mapper
    mapper = EnhancedPDBMapper(prefer_pdbml=args.prefer_pdbml)
    
    # Find epitope files
    pattern = f"{args.gene}_*_linear_epitopes_with_pdb_mapping.tsv" if args.gene else "*_linear_epitopes_with_pdb_mapping.tsv"
    epitope_files = list(bepipred_dir.glob(f"**/{pattern}"))
    
    if not epitope_files:
        logging.error(f"No epitope files found matching pattern: {pattern}")
        sys.exit(1)
    
    logging.info(f"Found {len(epitope_files)} epitope files to process")
    
    # Process each file
    for epitope_file in epitope_files:
        process_bepipred_epitope_file(epitope_file, mapper, pdb_structures_dir)
    
    logging.info("Enhanced mapping complete!")

if __name__ == "__main__":
    main()