#!/usr/bin/env python3
"""
Simple protein download script for Snakemake
Directly downloads proteins without subprocess calls
"""

import warnings
warnings.filterwarnings("ignore", message=".*Signature.*longdouble.*")

import pandas as pd
from Bio import Entrez, SeqIO
import time
import requests
from datetime import datetime
import os
import json
from pathlib import Path
from typing import List, Dict, Any, Set, Tuple
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def download_from_uniprot_batch(gene_name: str, species_list: List[str]) -> Dict[str, str]:
    """Download proteins from UniProt using batch query"""
    logging.info(f"Downloading {gene_name} from UniProt (batch mode)")
    
    # Create query for all species
    organism_queries = [f'(organism:"{species}")' for species in species_list]
    query = f'(gene:{gene_name}) AND ({" OR ".join(organism_queries)})'
    
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': query,
        'format': 'fasta',
        'size': 500
    }
    
    sequences = {}
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        
        if response.text.strip():
            # Parse FASTA sequences
            from io import StringIO
            for record in SeqIO.parse(StringIO(response.text), "fasta"):
                # Extract species from description
                desc = record.description.lower()
                for species in species_list:
                    if species.lower() in desc:
                        sequences[species] = record
                        break
                        
    except Exception as e:
        logging.warning(f"UniProt batch download failed: {e}")
    
    logging.info(f"UniProt batch found {len(sequences)} sequences")
    return sequences

def download_from_ncbi(gene_name: str, species: str, email: str = "user@example.com") -> str:
    """Download single protein from NCBI"""
    try:
        Entrez.email = email
        
        # Search for protein
        search_term = f'"{gene_name}"[Gene Name] AND "{species}"[Organism]'
        search_handle = Entrez.esearch(db="protein", term=search_term, retmax=5)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        if search_results["IdList"]:
            # Get the first protein
            protein_id = search_results["IdList"][0]
            
            # Fetch the sequence
            fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
            fasta_seq = fetch_handle.read()
            fetch_handle.close()
            
            return fasta_seq
            
    except Exception as e:
        logging.debug(f"NCBI download failed for {species}: {e}")
    
    return None

def process_gene(gene_name: str, species_file: str, output_dir: str, ncbi_email: str = "user@example.com"):
    """Process downloads for a single gene"""
    logging.info(f"Processing gene: {gene_name}")
    
    # Read species list
    with open(species_file, 'r') as f:
        species_list = [line.strip() for line in f if line.strip()]
    
    logging.info(f"Found {len(species_list)} species to download")
    
    # Create gene directory
    gene_dir = Path(output_dir) / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Try UniProt batch download
    uniprot_sequences = download_from_uniprot_batch(gene_name, species_list)
    
    # Step 2: Download missing species from NCBI
    missing_species = [s for s in species_list if s not in uniprot_sequences]
    ncbi_sequences = {}
    
    if missing_species:
        logging.info(f"Downloading {len(missing_species)} missing species from NCBI")
        for i, species in enumerate(missing_species[:50]):  # Limit NCBI downloads
            if i % 10 == 0:
                logging.info(f"NCBI progress: {i+1}/{min(len(missing_species), 50)}")
            
            fasta_seq = download_from_ncbi(gene_name, species, ncbi_email)
            if fasta_seq:
                # Parse the sequence
                from io import StringIO
                try:
                    record = next(SeqIO.parse(StringIO(fasta_seq), "fasta"))
                    ncbi_sequences[species] = record
                except:
                    pass
            
            time.sleep(0.2)  # Rate limiting
    
    # Step 3: Save all sequences
    all_sequences = {**uniprot_sequences, **ncbi_sequences}
    
    for species, record in all_sequences.items():
        # Clean species name for filename
        clean_species = species.replace(' ', '_').replace('/', '_')
        output_file = gene_dir / f"{clean_species}.fasta"
        
        with open(output_file, 'w') as f:
            SeqIO.write(record, f, "fasta")
    
    # Save summary
    summary = {
        "gene": gene_name,
        "total_species": len(species_list),
        "uniprot_found": len(uniprot_sequences),
        "ncbi_found": len(ncbi_sequences),
        "total_found": len(all_sequences),
        "coverage_percent": (len(all_sequences) / len(species_list)) * 100
    }
    
    summary_file = gene_dir / "download_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logging.info(f"Gene {gene_name}: {len(all_sequences)}/{len(species_list)} sequences ({summary['coverage_percent']:.1f}%)")

def main():
    """Main function for Snakemake integration"""
    logging.info("Simple Protein Download")
    logging.info("=" * 50)
    
    try:
        # Get inputs from Snakemake
        protein_lists_dir = snakemake.input.protein_lists
        ncbi_config = snakemake.input.ncbi_info
        output_dir = snakemake.output.download_dir
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        logging.info(f"Protein lists directory: {protein_lists_dir}")
        logging.info(f"Output directory: {output_dir}")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        protein_lists_dir = "results/proteins_to_download/analysis_1_params_1_gram_positive"
        output_dir = "results/protein_fasta/test"
        ncbi_config = "config/login/ncbi_info.txt"
    
    # Read NCBI email if available
    ncbi_email = "user@example.com"
    try:
        if Path(ncbi_config).exists():
            with open(ncbi_config, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'email' in line.lower():
                        ncbi_email = line.split('=')[1].strip()
                        break
    except:
        pass
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process all gene files in the protein lists directory
    protein_lists_path = Path(protein_lists_dir)
    if not protein_lists_path.exists():
        logging.error(f"Protein lists directory not found: {protein_lists_dir}")
        return
    
    # Find all .txt files (gene species lists)
    gene_files = list(protein_lists_path.glob("*.txt"))
    logging.info(f"Found {len(gene_files)} gene files to process")
    
    if not gene_files:
        logging.warning(f"No gene files found in {protein_lists_dir}")
        return
    
    # Process each gene
    total_sequences = 0
    successful_genes = 0
    
    for gene_file in gene_files:
        gene = gene_file.stem  # filename without extension
        try:
            process_gene(gene, str(gene_file), output_dir, ncbi_email)
            successful_genes += 1
        except Exception as e:
            logging.error(f"Failed to process gene {gene}: {e}")
    
    # Create sentinel file to indicate completion
    sentinel_file = Path(output_dir) / ".download_complete"
    with open(sentinel_file, 'w') as f:
        f.write(f"Download completed at {datetime.now().isoformat()}\n")
        f.write(f"Processed {successful_genes}/{len(gene_files)} genes\n")
    
    logging.info(f"\nDownload Summary:")
    logging.info(f"Genes processed: {successful_genes}/{len(gene_files)}")
    logging.info(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()