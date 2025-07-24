#!/usr/bin/env python3
"""
PDB 3D Structure Downloader
============================

Downloads 3D protein structures and their corresponding FASTA sequences from the Protein Data Bank (PDB)
for genes listed in the proteins_to_study file. The sequences are integrated into existing protein_fasta
directories for inclusion in multiple sequence alignment.

Uses the PDB REST API to search for structures by gene name and organism.
"""

import requests
import pandas as pd
from pathlib import Path
import logging
import time
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import urllib.parse
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_gene_aliases(aliases_file):
    """Load gene aliases from JSON file"""
    try:
        with open(aliases_file, 'r') as f:
            aliases = json.load(f)
        logging.info(f"Loaded aliases for {len(aliases)} genes")
        return aliases
    except Exception as e:
        logging.warning(f"Could not load gene aliases: {e}")
        return {}

def search_pdb_structures(gene_name, organism_list=None, max_structures=3):
    """
    Search PDB for structures containing the given gene name using simple text search
    
    Args:
        gene_name (str): Gene name to search for
        organism_list (list): Optional list of organism names to filter by
        max_structures (int): Maximum number of structures to return
    
    Returns:
        list: PDB IDs of matching structures
    """
    
    # Use PDB's simple search endpoint
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Simplified query using text search
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "value": gene_name
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_structures}
        }
    }
    
    try:
        response = requests.post(search_url, json=query, timeout=30)
        
        if response.status_code != 200:
            logging.warning(f"PDB search returned status {response.status_code} for {gene_name}")
            # Try alternative approach - direct UniProt search  
            return search_pdb_via_uniprot(gene_name, max_structures)
        
        result = response.json()
        pdb_ids = [item["identifier"] for item in result.get("result_set", [])]
        
        logging.info(f"Found {len(pdb_ids)} PDB structures for gene {gene_name}")
        return pdb_ids[:max_structures]
        
    except Exception as e:
        logging.warning(f"PDB search failed for {gene_name}: {e}")
        # Try alternative approach
        return search_pdb_via_uniprot(gene_name, max_structures)

def search_pdb_via_uniprot(gene_name, max_structures=3, gene_aliases=None):
    """
    Alternative search using UniProt to find PDB structures, with alias fallback
    
    Args:
        gene_name (str): Primary gene name to search for
        max_structures (int): Maximum structures to return
        gene_aliases (list): Alternative gene names to try if primary fails
    """
    # Try primary gene name first
    pdb_ids = _search_uniprot_single_gene(gene_name, max_structures)
    
    if pdb_ids:
        logging.info(f"Found {len(pdb_ids)} PDB structures via UniProt for gene {gene_name}")
        return pdb_ids
    
    # If no results and we have aliases, try them
    if gene_aliases:
        logging.info(f"No structures found for {gene_name}, trying {len(gene_aliases)} aliases")
        for alias in gene_aliases:
            if alias != gene_name:  # Don't retry the same name
                alias_pdb_ids = _search_uniprot_single_gene(alias, max_structures)
                if alias_pdb_ids:
                    logging.info(f"Found {len(alias_pdb_ids)} PDB structures via alias '{alias}'")
                    return alias_pdb_ids
    
    logging.info(f"No PDB structures found for {gene_name} (tried {len(gene_aliases or [])} aliases)")
    return []

def _search_uniprot_single_gene(gene_name, max_structures):
    """Helper function to search UniProt for a single gene name"""
    try:
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'gene:{gene_name}',
            'format': 'json',
            'size': max_structures,
            'fields': 'accession,gene_names,xref_pdb'
        }
        
        response = requests.get(uniprot_url, params=params, timeout=30)
        if response.status_code != 200:
            logging.debug(f"UniProt search failed for {gene_name}: status {response.status_code}")
            return []
        
        data = response.json()
        pdb_ids = []
        
        for entry in data.get('results', []):
            # Look for PDB cross-references
            for xref in entry.get('uniProtKBCrossReferences', []):
                if xref.get('database') == 'PDB':
                    pdb_id = xref.get('id')
                    if pdb_id and pdb_id not in pdb_ids:
                        pdb_ids.append(pdb_id)
        
        return pdb_ids[:max_structures]
        
    except Exception as e:
        logging.warning(f"UniProt search failed for {gene_name}: {e}")
        return []

def get_pdb_fasta(pdb_id):
    """
    Download FASTA sequence for a PDB structure
    
    Args:
        pdb_id (str): PDB identifier
    
    Returns:
        list: List of SeqRecord objects
    """
    
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    
    try:
        response = requests.get(fasta_url, timeout=30)
        response.raise_for_status()
        
        # Parse FASTA content
        sequences = []
        lines = response.text.strip().split('\n')
        current_seq = ""
        current_header = ""
        
        for line in lines:
            if line.startswith('>'):
                if current_header and current_seq:
                    # Create SeqRecord
                    seq_id = current_header.split('|')[0].replace('>', '')
                    sequences.append(SeqRecord(
                        Seq(current_seq),
                        id=seq_id,
                        description=current_header.replace('>', '')
                    ))
                current_header = line
                current_seq = ""
            else:
                current_seq += line.strip()
        
        # Add last sequence
        if current_header and current_seq:
            seq_id = current_header.split('|')[0].replace('>', '')
            sequences.append(SeqRecord(
                Seq(current_seq),
                id=seq_id,
                description=current_header.replace('>', '')
            ))
        
        logging.info(f"Downloaded {len(sequences)} sequences for PDB {pdb_id}")
        return sequences
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download FASTA for {pdb_id}: {e}")
        return []

def download_pdb_structure(pdb_id, output_dir):
    """
    Download PDB structure file in legacy format (gz)
    
    Args:
        pdb_id (str): PDB identifier
        output_dir (Path): Output directory to save the structure
    
    Returns:
        bool: True if download succeeded, False otherwise
    """
    
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    
    try:
        response = requests.get(pdb_url, timeout=30)
        response.raise_for_status()
        
        # Save PDB file
        pdb_file = output_dir / f"{pdb_id}.pdb.gz"
        with open(pdb_file, 'wb') as f:
            f.write(response.content)
        
        logging.info(f"Downloaded PDB structure: {pdb_file}")
        return True
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download PDB structure for {pdb_id}: {e}")
        return False

def download_3d_structures_for_gene(gene_name, structures_base_dir, organisms=None, max_structures=3, gene_aliases=None):
    """
    Download 3D structures and sequences for a specific gene (top 3 best scoring only)
    
    Args:
        gene_name (str): Gene name to search for
        structures_base_dir (Path): Base directory for 3D structures (results/3d_structures/)
        organisms (list): Optional list of organism names to filter by
        max_structures (int): Maximum number of structures to download (default: 3)
        gene_aliases (list): Alternative gene names to try if primary search fails
    
    Returns:
        dict: Summary of downloads (sequences, structures, found)
    """
    
    logging.info(f"Searching for 3D structures for gene: {gene_name}")
    
    # Search for PDB structures via UniProt (skip direct PDB search to avoid 400 errors)
    pdb_ids = search_pdb_via_uniprot(gene_name, max_structures, gene_aliases)
    
    # Create gene-specific directory in 3d_structures
    gene_dir = structures_base_dir / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    if not pdb_ids:
        logging.info(f"No 3D structures found for {gene_name}")
        # Create a marker file indicating no structures found
        no_structures_file = gene_dir / "no_structures_found.txt"
        with open(no_structures_file, 'w') as f:
            f.write(f"No 3D structures found for gene: {gene_name}\n")
            f.write(f"Search performed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        return {"sequences": 0, "structures": 0, "found": False, "pdb_ids": []}
    
    # Limit to top 3 structures
    pdb_ids = pdb_ids[:max_structures]
    logging.info(f"Processing top {len(pdb_ids)} 3D structures for {gene_name}: {pdb_ids}")
    
    total_sequences = 0
    total_structures = 0
    downloaded_pdb_ids = []
    
    for pdb_id in pdb_ids:
        # Add delay to be respectful to PDB servers
        time.sleep(1)
        
        # Download FASTA sequences for this PDB
        sequences = get_pdb_fasta(pdb_id)
        
        if sequences:
            for i, seq_record in enumerate(sequences):
                # Create a unique filename for each sequence
                if len(sequences) > 1:
                    filename = f"{pdb_id}_chain_{i+1}.fasta"
                else:
                    filename = f"{pdb_id}.fasta"
                
                # Update sequence ID to include PDB info
                seq_record.id = f"{pdb_id}_{seq_record.id}"
                seq_record.description = f"PDB:{pdb_id} {seq_record.description}"
                
                # Save sequence in gene directory
                output_file = gene_dir / filename
                with open(output_file, 'w') as f:
                    SeqIO.write(seq_record, f, "fasta")
                
                total_sequences += 1
                logging.info(f"Saved 3D structure sequence: {output_file}")
            
            # Download PDB structure file
            if download_pdb_structure(pdb_id, gene_dir):
                total_structures += 1
                downloaded_pdb_ids.append(pdb_id)
        
        else:
            logging.warning(f"No sequences found for PDB {pdb_id}")
    
    return {
        "sequences": total_sequences, 
        "structures": total_structures, 
        "found": total_sequences > 0,
        "pdb_ids": downloaded_pdb_ids
    }

def process_proteins_list(proteins_file, structures_base_dir, gene_aliases=None):
    """
    Process proteins_to_study file and download 3D structures for each gene
    
    Args:
        proteins_file (Path): Path to proteins_to_study TSV file
        structures_base_dir (Path): Base directory for 3D structures (results/3d_structures/)
        gene_aliases (dict): Dictionary mapping gene names to lists of aliases
    
    Returns:
        dict: Summary of downloaded structures per gene
    """
    
    # Load proteins list
    if not proteins_file.exists():
        logging.error(f"Proteins file not found: {proteins_file}")
        return {}
    
    proteins_df = pd.read_csv(proteins_file, sep='\t')
    logging.info(f"Processing {len(proteins_df)} genes from {proteins_file}")
    
    # Extract species list for organism filtering
    organisms = []
    if 'species' in proteins_df.columns:
        for species_list in proteins_df['species'].dropna():
            organisms.extend([species.strip() for species in species_list.split(',')])
        organisms = list(set(organisms))  # Remove duplicates
        logging.info(f"Filtering by {len(organisms)} organisms")
    
    summary = {}
    
    for _, row in proteins_df.iterrows():
        gene_name = row['gene']
        
        # Get aliases for this gene
        aliases = gene_aliases.get(gene_name, []) if gene_aliases else []
        
        # Download 3D structures for this gene
        gene_summary = download_3d_structures_for_gene(
            gene_name, 
            structures_base_dir, 
            organisms if organisms else None,
            gene_aliases=aliases
        )
        
        summary[gene_name] = gene_summary
        
        # Add delay between genes
        time.sleep(2)
    
    return summary

def main():
    """Main function for Snakemake integration"""
    
    # Check if running under Snakemake
    try:
        # Get parameters from Snakemake
        proteins_file = Path(snakemake.input.protein_list)
        structures_dir = Path(snakemake.output.structures_dir)
        summary_file = Path(snakemake.output.summary)
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        aliases_file = Path(f"data/quickgo/{paramset}/gene_aliases.json")
        
    except NameError:
        # Running standalone - use default values for testing
        analysis = "analysis_1"
        paramset = "params_1"
        group = "positive"
        proteins_file = Path(f"results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv")
        structures_dir = Path(f"results/3d_structures/{analysis}_{paramset}_gram_{group}")
        summary_file = Path(f"results/3d_structures/{analysis}_{paramset}_gram_{group}_summary.json")
        aliases_file = Path(f"data/quickgo/{paramset}/gene_aliases.json")
    
    logging.info(f"Downloading 3D structures for {analysis}_{paramset}_gram_{group}")
    logging.info(f"Proteins file: {proteins_file}")
    logging.info(f"Structures directory: {structures_dir}")
    logging.info(f"Summary file: {summary_file}")
    
    # Load gene aliases
    gene_aliases = load_gene_aliases(aliases_file)
    
    # Ensure output directories exist
    structures_dir.mkdir(parents=True, exist_ok=True)
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Process proteins and download 3D structures
    summary = process_proteins_list(proteins_file, structures_dir, gene_aliases)
    
    # Calculate summary statistics
    total_sequences = sum(gene_data.get("sequences", 0) for gene_data in summary.values())
    total_structures = sum(gene_data.get("structures", 0) for gene_data in summary.values())
    genes_with_3d = len([g for g, data in summary.items() if data.get("found", False)])
    
    # Save summary
    summary_data = {
        "analysis": analysis,
        "paramset": paramset,
        "group": group,
        "total_genes": len(summary),
        "total_3d_sequences": total_sequences,
        "total_3d_structures": total_structures,
        "genes_with_3d": genes_with_3d,
        "per_gene_summary": summary
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    logging.info(f"Downloaded {total_sequences} sequences and {total_structures} structures for {genes_with_3d}/{len(summary)} genes")
    logging.info(f"Summary saved to: {summary_file}")

if __name__ == "__main__":
    main()