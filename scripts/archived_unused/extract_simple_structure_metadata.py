#!/usr/bin/env python3
"""
Simple Structure Metadata Extractor
===================================

Extracts key metadata for protein structures with focus on the information requested:
- Method (e.g., X-RAY DIFFRACTION)
- Resolution (e.g., 2.60 Å)
- Entity ID, Molecule name, Chains, Sequence Length, Organism

This script provides essential structure information in a clean, accessible format.
"""

import requests
import pandas as pd
from pathlib import Path
import logging
import json
import time
from datetime import datetime
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_pdb_structure_info(pdb_id):
    """
    Get essential structure information for a PDB entry
    
    Args:
        pdb_id (str): PDB identifier
    
    Returns:
        dict: Structure information
    """
    try:
        # Get basic entry information
        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(entry_url, timeout=30)
        response.raise_for_status()
        entry_data = response.json()
        
        # Get FASTA sequences to determine chain information
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
        fasta_response = requests.get(fasta_url, timeout=30)
        
        # Parse basic information
        structure_info = {
            "pdb_id": pdb_id,
            "title": entry_data.get("struct", {}).get("title", "Unknown"),
            "method": "Unknown",
            "resolution": None,
            "resolution_unit": "Å",
            "entities": [],
            "total_chains": 0,
            "deposit_date": entry_data.get("rcsb_accession_info", {}).get("deposit_date", "Unknown"),
            "release_date": entry_data.get("rcsb_accession_info", {}).get("initial_release_date", "Unknown")
        }
        
        # Get experimental method
        exptl_methods = entry_data.get("exptl", [])
        if exptl_methods:
            structure_info["method"] = exptl_methods[0].get("method", "Unknown")
        
        # Get resolution
        refine_info = entry_data.get("refine", [])
        if refine_info:
            resolution = refine_info[0].get("ls_d_res_high")
            if resolution:
                structure_info["resolution"] = float(resolution)
        
        # Try alternative resolution sources
        if not structure_info["resolution"]:
            rcsb_entry_info = entry_data.get("rcsb_entry_info", {})
            resolution_combined = rcsb_entry_info.get("resolution_combined")
            if resolution_combined:
                structure_info["resolution"] = float(resolution_combined[0])
        
        # Parse FASTA to get chain and sequence information
        if fasta_response.status_code == 200:
            fasta_text = fasta_response.text
            chains = []
            current_chain = {}
            current_sequence = ""
            
            for line in fasta_text.split('\n'):
                if line.startswith('>'):
                    # Save previous chain if exists
                    if current_chain and current_sequence:
                        current_chain["sequence_length"] = len(current_sequence.replace(' ', '').replace('\n', ''))
                        current_chain["sequence"] = current_sequence.replace(' ', '').replace('\n', '')
                        chains.append(current_chain)
                    
                    # Parse header: >1OAP_1|Chain A|PEPTIDOGLYCAN-ASSOCIATED LIPOPROTEIN|Escherichia coli
                    parts = line.replace('>', '').split('|')
                    chain_id = parts[1].replace('Chain ', '') if len(parts) > 1 else "Unknown"
                    molecule_name = parts[2] if len(parts) > 2 else "Unknown"
                    organism = parts[3] if len(parts) > 3 else "Unknown"
                    
                    current_chain = {
                        "entity_id": len(chains) + 1,
                        "chain_id": chain_id,
                        "molecule_name": molecule_name,
                        "organism": organism,
                        "sequence_length": 0,
                        "sequence": ""
                    }
                    current_sequence = ""
                else:
                    current_sequence += line.strip()
            
            # Save last chain
            if current_chain and current_sequence:
                current_chain["sequence_length"] = len(current_sequence.replace(' ', '').replace('\n', ''))
                current_chain["sequence"] = current_sequence.replace(' ', '').replace('\n', '')
                chains.append(current_chain)
            
            structure_info["entities"] = chains
            structure_info["total_chains"] = len(chains)
        
        return structure_info
        
    except Exception as e:
        logging.warning(f"Failed to get structure info for PDB {pdb_id}: {e}")
        return {
            "pdb_id": pdb_id,
            "title": "Unknown",
            "method": "Unknown",
            "resolution": None,
            "entities": [],
            "total_chains": 0,
            "error": str(e)
        }

def get_alphafold_structure_info(alphafold_id):
    """
    Get structure information for an AlphaFold model
    
    Args:
        alphafold_id (str): AlphaFold identifier (e.g., AF-P12345)
    
    Returns:
        dict: Structure information
    """
    uniprot_accession = alphafold_id.replace('AF-', '')
    
    try:
        # Get AlphaFold metadata
        af_api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
        response = requests.get(af_api_url, timeout=30)
        response.raise_for_status()
        
        af_data = response.json()
        if not af_data:
            raise ValueError("No AlphaFold data found")
        
        af_entry = af_data[0] if isinstance(af_data, list) else af_data
        
        structure_info = {
            "alphafold_id": alphafold_id,
            "uniprot_accession": uniprot_accession,
            "title": f"AlphaFold model for {af_entry.get('gene', 'Unknown gene')}",
            "method": "Computational Prediction (AlphaFold)",
            "resolution": "N/A (Computed Model)",
            "confidence": af_entry.get('globalMetricValue', 'Unknown'),
            "model_version": af_entry.get('modelVersion', 'Unknown'),
            "created_date": af_entry.get('modelCreatedDate', 'Unknown'),
            "entities": [{
                "entity_id": 1,
                "chain_id": "A",
                "molecule_name": f"AlphaFold model of {af_entry.get('gene', 'Unknown')}",
                "organism": af_entry.get('organismScientificName', 'Unknown'),
                "sequence_length": af_entry.get('sequenceLength', 0),
                "gene_name": af_entry.get('gene', 'Unknown'),
                "confidence_score": af_entry.get('globalMetricValue', 'Unknown')
            }],
            "total_chains": 1
        }
        
        return structure_info
        
    except Exception as e:
        logging.warning(f"Failed to get AlphaFold info for {alphafold_id}: {e}")
        return {
            "alphafold_id": alphafold_id,
            "title": "Unknown",
            "method": "Computational Prediction (AlphaFold)",
            "resolution": "N/A (Computed Model)",
            "entities": [],
            "total_chains": 0,
            "error": str(e)
        }

def extract_structures_from_gene_directory(gene_dir):
    """
    Extract structure IDs from a gene directory
    
    Args:
        gene_dir (Path): Path to gene directory
    
    Returns:
        tuple: (pdb_ids, alphafold_ids)
    """
    pdb_ids = []
    alphafold_ids = []
    
    # Look for PDB structure files
    pdb_files = list(gene_dir.glob("*.pdb.gz"))
    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem.replace('.pdb', '')
        if not pdb_id.startswith('AF-'):
            pdb_ids.append(pdb_id)
    
    # Look for AlphaFold model files
    af_files = list(gene_dir.glob("AF-*.pdb"))
    for af_file in af_files:
        alphafold_id = af_file.stem
        alphafold_ids.append(alphafold_id)
    
    return pdb_ids, alphafold_ids

def process_gene_structures(gene_name, structures_base_dir):
    """
    Process all structures for a specific gene
    
    Args:
        gene_name (str): Gene name
        structures_base_dir (Path): Base directory containing gene structure folders
    
    Returns:
        dict: Gene structure metadata
    """
    gene_dir = structures_base_dir / gene_name
    
    if not gene_dir.exists():
        logging.warning(f"Gene directory not found: {gene_dir}")
        return {
            "gene_name": gene_name,
            "status": "directory_not_found",
            "structures": []
        }
    
    # Extract structure IDs
    pdb_ids, alphafold_ids = extract_structures_from_gene_directory(gene_dir)
    
    if not pdb_ids and not alphafold_ids:
        logging.info(f"No structures found for gene {gene_name}")
        return {
            "gene_name": gene_name,
            "status": "no_structures_found",
            "structures": []
        }
    
    logging.info(f"Processing {len(pdb_ids)} experimental + {len(alphafold_ids)} computed structures for gene {gene_name}")
    
    structures = []
    
    # Process experimental structures
    for pdb_id in pdb_ids:
        logging.info(f"Getting structure info for PDB {pdb_id}")
        structure_info = get_pdb_structure_info(pdb_id)
        structure_info["structure_type"] = "experimental"
        structures.append(structure_info)
        time.sleep(0.5)  # Rate limiting
    
    # Process computed models
    for alphafold_id in alphafold_ids:
        logging.info(f"Getting structure info for AlphaFold model {alphafold_id}")
        structure_info = get_alphafold_structure_info(alphafold_id)
        structure_info["structure_type"] = "computed_model"
        structures.append(structure_info)
        time.sleep(0.5)  # Rate limiting
    
    return {
        "gene_name": gene_name,
        "status": "completed",
        "experimental_structures_count": len(pdb_ids),
        "computed_models_count": len(alphafold_ids),
        "total_structures_count": len(structures),
        "structures": structures,
        "processing_date": datetime.now().isoformat()
    }

def main():
    """Main function for structure metadata extraction"""
    
    parser = argparse.ArgumentParser(description='Extract structure metadata for genes')
    parser.add_argument('--structures_dir', type=str, default='data/protein_structures',
                        help='Directory containing gene structure folders')
    parser.add_argument('--gene', type=str, help='Process specific gene only')
    parser.add_argument('--summary_file', type=str, help='Summary TSV file with gene list')
    
    args = parser.parse_args()
    
    structures_dir = Path(args.structures_dir)
    
    if not structures_dir.exists():
        logging.error(f"Structures directory not found: {structures_dir}")
        return
    
    # Determine which genes to process
    genes_to_process = []
    
    if args.gene:
        # Process single gene
        genes_to_process = [args.gene]
    elif args.summary_file:
        # Process genes from summary file
        summary_file = Path(args.summary_file)
        if summary_file.exists():
            df = pd.read_csv(summary_file, sep='\t')
            genes_to_process = df['gene_name'].unique().tolist()
            logging.info(f"Processing {len(genes_to_process)} genes from {summary_file}")
        else:
            logging.error(f"Summary file not found: {summary_file}")
            return
    else:
        # Process all gene directories
        gene_dirs = [d for d in structures_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
        genes_to_process = [d.name for d in gene_dirs]
        logging.info(f"Processing {len(genes_to_process)} genes from structure directories")
    
    if not genes_to_process:
        logging.error("No genes to process")
        return
    
    logging.info(f"Extracting structure metadata for {len(genes_to_process)} genes")
    
    # Process all genes
    for i, gene_name in enumerate(genes_to_process, 1):
        logging.info(f"Progress: {i}/{len(genes_to_process)} - Processing gene: {gene_name}")
        
        gene_data = process_gene_structures(gene_name, structures_dir)
        
        # Save individual gene metadata in the gene's folder
        gene_output_file = structures_dir / gene_name / "structure_metadata.json"
        if gene_output_file.parent.exists():
            with open(gene_output_file, 'w') as f:
                json.dump(gene_data, f, indent=2)
            logging.info(f"Saved structure metadata: {gene_output_file}")
    
    logging.info(f"=== Structure Metadata Extraction Complete ===")

if __name__ == "__main__":
    main()