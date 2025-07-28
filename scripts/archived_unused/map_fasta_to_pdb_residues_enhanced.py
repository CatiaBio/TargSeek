#!/usr/bin/env python3
"""
Enhanced FASTA to PDB Residue Mapping
=====================================

This script creates mappings between FASTA sequence positions and PDB residue numbers
using XML structure files. It handles both experimental structures (with XML files)
and computational models (AlphaFold, where FASTA positions = PDB residues).

For experimental structures, it parses the XML to extract the correspondence between:
- FASTA position (1-based indexing in the FASTA file)
- PDB residue number (auth_seq_id from the XML structure)

The mapping enables correct interpretation of epitope predictions relative to
the original protein sequence numbering.
"""

import pandas as pd
import gzip
import xml.etree.ElementTree as ET
import logging
from pathlib import Path
import json
from Bio import SeqIO
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class FASTAToPDBMapper:
    """Maps FASTA positions to PDB residue numbers using XML structure files"""
    
    def __init__(self):
        self.mappings = {}
        
    def is_computational_model(self, pdb_id):
        """Check if this is a computational model (AlphaFold, ChimeraX, etc.)"""
        computational_prefixes = ['AF-', 'CHM-', 'MODEL-']
        return any(pdb_id.startswith(prefix) for prefix in computational_prefixes)
    
    def parse_xml_for_chain_mapping(self, xml_path, chain_id):
        """
        Parse XML file to extract residue mapping for a specific chain
        
        Args:
            xml_path: Path to the XML.gz file
            chain_id: Chain identifier (e.g., 'A', '1', '2')
            
        Returns:
            dict: {fasta_position: pdb_residue_number}
        """
        try:
            # Handle both .gz and uncompressed files
            if xml_path.suffix == '.gz':
                with gzip.open(xml_path, 'rt') as f:
                    tree = ET.parse(f)
            else:
                tree = ET.parse(xml_path)
            
            root = tree.getroot()
            
            # Define namespace
            ns = {'PDBx': 'http://pdbml.pdb.org/schema/pdbx-v50.xsd'}
            
            # Find atom sites for the specific chain
            mapping = {}
            seen_label_seq_ids = set()
            
            for atom_site in root.findall('.//PDBx:atom_site', ns):
                # Get chain information
                auth_asym_id = atom_site.find('PDBx:auth_asym_id', ns)
                label_asym_id = atom_site.find('PDBx:label_asym_id', ns)
                
                # Check if this atom belongs to our chain
                chain_match = False
                if auth_asym_id is not None and auth_asym_id.text == chain_id:
                    chain_match = True
                elif label_asym_id is not None and label_asym_id.text == chain_id:
                    chain_match = True
                
                if not chain_match:
                    continue
                
                # Get sequence positions
                auth_seq_id = atom_site.find('PDBx:auth_seq_id', ns)
                label_seq_id = atom_site.find('PDBx:label_seq_id', ns)
                
                if auth_seq_id is not None and label_seq_id is not None:
                    try:
                        # Check if text is not None before converting
                        if auth_seq_id.text is not None and label_seq_id.text is not None:
                            pdb_residue = int(auth_seq_id.text)
                            fasta_entity_pos = int(label_seq_id.text)
                            
                            # Only record each label_seq_id once (avoid duplicate atoms)
                            if fasta_entity_pos not in seen_label_seq_ids:
                                mapping[fasta_entity_pos] = pdb_residue
                                seen_label_seq_ids.add(fasta_entity_pos)
                    except (ValueError, TypeError):
                        continue
            
            return mapping
            
        except Exception as e:
            logging.error(f"Error parsing XML {xml_path} for chain {chain_id}: {e}")
            return {}
    
    def get_chain_id_from_fasta_file(self, fasta_path):
        """
        Extract chain ID from FASTA file header information
        
        Examples:
            5D0Q_chain_1.fasta with header "Chains A, E[auth F]" -> 'A'
            2MPR.fasta -> 'A' (default)
            AF-Q9HTH0.fasta -> 'A' (computational model)
        """
        try:
            # First try to get chain info from FASTA header
            sequences = list(SeqIO.parse(fasta_path, "fasta"))
            if sequences:
                header = sequences[0].description
                
                # Look for chain information in header
                # Pattern: "Chains A, E[auth F]" or "Chain A"
                if "|Chains " in header:
                    chain_info = header.split("|Chains ")[1].split("|")[0]
                    # Take the first chain mentioned
                    first_chain = chain_info.split(",")[0].split("[")[0].strip()
                    return first_chain
                elif "|Chain " in header:
                    chain_info = header.split("|Chain ")[1].split("|")[0]
                    return chain_info.split("[")[0].strip()
        except Exception as e:
            logging.warning(f"Could not parse FASTA header for {fasta_path}: {e}")
        
        # Fallback to filename-based detection
        filename = Path(fasta_path).stem
        
        if '_chain_' in filename:
            return filename.split('_chain_')[-1]
        else:
            # Default to chain A for single-chain structures
            return 'A'
    
    def get_pdb_id_from_fasta_path(self, fasta_path):
        """Extract PDB ID from FASTA file path"""
        filename = Path(fasta_path).stem
        
        if '_chain_' in filename:
            return filename.split('_chain_')[0]
        else:
            return filename
    
    def get_xml_path_from_fasta_path(self, fasta_path):
        """Get corresponding XML path from FASTA path"""
        fasta_path = Path(fasta_path)
        pdb_id = self.get_pdb_id_from_fasta_path(fasta_path)
        
        # XML file should be in the same directory as FASTA
        xml_path = fasta_path.parent / f"{pdb_id}.xml.gz"
        
        if not xml_path.exists():
            # Try without .gz extension
            xml_path = fasta_path.parent / f"{pdb_id}.xml"
        
        return xml_path if xml_path.exists() else None
    
    def create_fasta_to_pdb_mapping(self, fasta_path):
        """
        Create mapping from FASTA position to PDB residue number
        
        Args:
            fasta_path: Path to the FASTA file
            
        Returns:
            dict: {
                'pdb_id': str,
                'chain_id': str, 
                'is_computational': bool,
                'mapping': {fasta_pos: pdb_residue, ...},
                'fasta_length': int,
                'mapping_coverage': float
            }
        """
        fasta_path = Path(fasta_path)
        pdb_id = self.get_pdb_id_from_fasta_path(fasta_path)
        chain_id = self.get_chain_id_from_fasta_file(fasta_path)
        
        logging.info(f"Processing {fasta_path.name}: PDB={pdb_id}, Chain={chain_id}")
        
        # Read FASTA to get sequence length
        try:
            sequences = list(SeqIO.parse(fasta_path, "fasta"))
            if not sequences:
                logging.error(f"No sequences found in {fasta_path}")
                return None
            
            fasta_length = len(sequences[0].seq)
            
        except Exception as e:
            logging.error(f"Error reading FASTA {fasta_path}: {e}")
            return None
        
        # Check if this is a computational model
        is_computational = self.is_computational_model(pdb_id)
        
        if is_computational:
            # For computational models, FASTA position = PDB residue number
            logging.info(f"Computational model detected: {pdb_id}")
            mapping = {i: i for i in range(1, fasta_length + 1)}
            coverage = 1.0
            
        else:
            # For experimental structures, parse XML
            xml_path = self.get_xml_path_from_fasta_path(fasta_path)
            
            if xml_path is None:
                logging.warning(f"No XML file found for {pdb_id}")
                # Fallback: assume 1:1 mapping
                mapping = {i: i for i in range(1, fasta_length + 1)}
                coverage = 0.0  # Unknown coverage
            else:
                logging.info(f"Parsing XML: {xml_path}")
                mapping = self.parse_xml_for_chain_mapping(xml_path, chain_id)
                
                if not mapping:
                    logging.warning(f"No mapping found in XML for chain {chain_id}")
                    # Fallback: assume 1:1 mapping
                    mapping = {i: i for i in range(1, fasta_length + 1)}
                    coverage = 0.0
                else:
                    coverage = len(mapping) / fasta_length
        
        result = {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'is_computational': is_computational,
            'mapping': mapping,
            'fasta_length': fasta_length,
            'mapping_coverage': coverage
        }
        
        logging.info(f"Mapping created: {len(mapping)} positions, {coverage:.2%} coverage")
        return result
    
    def create_mapping_for_structure_list(self, fasta_paths_file, output_file):
        """
        Create mappings for all structures in the list
        
        Args:
            fasta_paths_file: File containing list of FASTA paths
            output_file: Output JSON file for mappings
        """
        logging.info(f"Creating mappings from {fasta_paths_file}")
        
        # Read FASTA paths
        with open(fasta_paths_file, 'r') as f:
            fasta_paths = [line.strip() for line in f if line.strip()]
        
        logging.info(f"Found {len(fasta_paths)} FASTA files to process")
        
        all_mappings = {}
        
        for fasta_path in fasta_paths:
            if not os.path.exists(fasta_path):
                logging.warning(f"FASTA file not found: {fasta_path}")
                continue
            
            mapping_result = self.create_fasta_to_pdb_mapping(fasta_path)
            if mapping_result:
                # Use FASTA path as key for easy lookup
                all_mappings[fasta_path] = mapping_result
        
        # Save mappings to JSON
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert int keys to strings for JSON serialization
        json_mappings = {}
        for fasta_path, result in all_mappings.items():
            json_result = result.copy()
            json_result['mapping'] = {str(k): v for k, v in result['mapping'].items()}
            json_mappings[fasta_path] = json_result
        
        with open(output_path, 'w') as f:
            json.dump(json_mappings, f, indent=2)
        
        logging.info(f"Saved mappings for {len(all_mappings)} structures to {output_path}")
        
        # Create summary report
        self.create_mapping_summary_report(all_mappings, output_path.parent / "mapping_summary.txt")
        
        return all_mappings
    
    def create_mapping_summary_report(self, mappings, report_file):
        """Create a human-readable summary report"""
        
        total_structures = len(mappings)
        computational_models = sum(1 for m in mappings.values() if m['is_computational'])
        experimental_structures = total_structures - computational_models
        
        with open(report_file, 'w') as f:
            f.write("FASTA to PDB Residue Mapping Summary Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total structures processed: {total_structures}\n")
            f.write(f"Computational models: {computational_models}\n")
            f.write(f"Experimental structures: {experimental_structures}\n\n")
            
            f.write("STRUCTURE DETAILS:\n")
            f.write("-" * 30 + "\n")
            
            for fasta_path, result in mappings.items():
                filename = Path(fasta_path).name
                f.write(f"\n{filename}:\n")
                f.write(f"  PDB ID: {result['pdb_id']}\n")
                f.write(f"  Chain: {result['chain_id']}\n")
                f.write(f"  Type: {'Computational' if result['is_computational'] else 'Experimental'}\n")
                f.write(f"  FASTA length: {result['fasta_length']} AA\n")
                f.write(f"  Mapping coverage: {result['mapping_coverage']:.2%}\n")
                
                if result['mapping']:
                    first_pos = min(result['mapping'].keys())
                    last_pos = max(result['mapping'].keys())
                    first_pdb = result['mapping'][first_pos]
                    last_pdb = result['mapping'][last_pos]
                    f.write(f"  PDB residue range: {first_pdb}-{last_pdb}\n")
                    f.write(f"  FASTA position range: {first_pos}-{last_pos}\n")
        
        logging.info(f"Created summary report: {report_file}")

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get parameters from Snakemake
        selected_3d_paths_positive = snakemake.input.selected_3d_paths_positive
        selected_3d_paths_negative = snakemake.input.selected_3d_paths_negative
        output_mapping_json = snakemake.output.mapping_json
        
        # Get analysis parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        
    except NameError:
        # Test mode - use command line arguments
        import sys
        if len(sys.argv) < 4:
            print("Usage: python map_fasta_to_pdb_residues_enhanced.py <selected_3d_paths_positive> <selected_3d_paths_negative> <output_mapping_json>")
            sys.exit(1)
        
        selected_3d_paths_positive = sys.argv[1]
        selected_3d_paths_negative = sys.argv[2]
        output_mapping_json = sys.argv[3]
        analysis = "analysis1"
        paramset = "params1"
    
    logging.info("=== Enhanced FASTA to PDB Residue Mapping ===")
    logging.info(f"Positive paths: {selected_3d_paths_positive}")
    logging.info(f"Negative paths: {selected_3d_paths_negative}")
    logging.info(f"Output mapping: {output_mapping_json}")
    
    # Initialize mapper
    mapper = FASTAToPDBMapper()
    
    # Read all FASTA paths
    all_fasta_paths = []
    
    for paths_file in [selected_3d_paths_positive, selected_3d_paths_negative]:
        if os.path.exists(paths_file):
            with open(paths_file, 'r') as f:
                paths = [line.strip() for line in f if line.strip()]
                all_fasta_paths.extend(paths)
                logging.info(f"Found {len(paths)} paths in {paths_file}")
    
    if not all_fasta_paths:
        logging.error("No FASTA paths found")
        return
    
    logging.info(f"Total FASTA files to process: {len(all_fasta_paths)}")
    
    # Create mappings for all structures
    all_mappings = {}
    
    for fasta_path in all_fasta_paths:
        if not os.path.exists(fasta_path):
            logging.warning(f"FASTA file not found: {fasta_path}")
            continue
        
        mapping_result = mapper.create_fasta_to_pdb_mapping(fasta_path)
        if mapping_result:
            all_mappings[fasta_path] = mapping_result
    
    # Save results
    output_path = Path(output_mapping_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert int keys to strings for JSON serialization
    json_mappings = {}
    for fasta_path, result in all_mappings.items():
        json_result = result.copy()
        json_result['mapping'] = {str(k): v for k, v in result['mapping'].items()}
        json_mappings[fasta_path] = json_result
    
    with open(output_path, 'w') as f:
        json.dump(json_mappings, f, indent=2)
    
    logging.info(f"Saved mappings for {len(all_mappings)} structures to {output_path}")
    
    # Create summary report
    mapper.create_mapping_summary_report(all_mappings, output_path.parent / "fasta_pdb_mapping_summary.txt")
    
    logging.info("=== FASTA to PDB Mapping Complete ===")

if __name__ == "__main__":
    main()