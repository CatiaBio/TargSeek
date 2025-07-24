#!/usr/bin/env python3
"""
Download missing protein sequences from NCBI using not_found_species.txt file.
This script complements UniProt downloads by fetching species that weren't found.
"""

import argparse
from Bio import Entrez, SeqIO
import os
import time
import urllib.error
from pathlib import Path
from typing import List, Dict, Set
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def safe_filename(name: str) -> str:
    """Convert species name to safe filename."""
    return name.replace(" ", "_").replace("/", "_")


def search_gene_species_ncbi(gene: str, species: str, retries: int = 3) -> List[str]:
    """
    Search NCBI for a specific gene in a specific species.
    
    Args:
        gene: Gene name
        species: Species name
        retries: Number of retries for failed requests
        
    Returns:
        List of protein IDs
    """
    query = f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species}"[Organism]'
    
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.esearch(db="protein", term=query, retmax=100)
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        except urllib.error.HTTPError as e:
            logging.warning(f"HTTP error for {gene}-{species} attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(2 ** attempt)
        except Exception as e:
            logging.error(f"Error searching for {gene}-{species}: {e}")
            attempt += 1
            time.sleep(2 ** attempt)
    return []


def fetch_protein_sequences(id_list: List[str], retries: int = 3) -> List:
    """
    Fetch protein sequences from NCBI by ID.
    
    Args:
        id_list: List of protein IDs
        retries: Number of retries for failed requests
        
    Returns:
        List of SeqRecord objects
    """
    if not id_list:
        return []
    
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(db="protein", id=",".join(id_list), rettype="fasta", retmode="text")
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            return records
        except Exception as e:
            logging.warning(f"Fetch error attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(2 ** attempt)
    return []


def check_existing_files(output_dir: Path, species_list: List[str]) -> Set[str]:
    """Check which species already have FASTA files."""
    existing_species = set()
    for species in species_list:
        fasta_file = output_dir / f"{safe_filename(species)}.fasta"
        if fasta_file.exists() and fasta_file.stat().st_size > 0:
            existing_species.add(species)
    return existing_species


def update_not_found_file(not_found_file: Path, found_species: Set[str]) -> None:
    """Update the not_found_species.txt file by removing found species."""
    if not found_species:
        return
    
    # Read current not found species
    with open(not_found_file, 'r') as f:
        not_found_species = [line.strip() for line in f if line.strip()]
    
    # Remove found species
    remaining_not_found = [sp for sp in not_found_species if sp not in found_species]
    
    # Write updated list
    with open(not_found_file, 'w') as f:
        for species in remaining_not_found:
            f.write(f"{species}\n")
    
    logging.info(f"Updated not_found file: removed {len(found_species)} species, {len(remaining_not_found)} remain")


def download_missing_proteins(gene: str, not_found_file: Path, output_dir: Path, 
                            email: str, api_key: str = None) -> Dict[str, int]:
    """
    Download missing proteins for a gene based on not_found_species.txt.
    
    Args:
        gene: Gene name
        not_found_file: Path to not_found_species.txt
        output_dir: Directory to save FASTA files
        email: Email for NCBI
        api_key: Optional NCBI API key
        
    Returns:
        Dictionary with download statistics
    """
    # Set up NCBI
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    # Read not found species
    with open(not_found_file, 'r') as f:
        missing_species = [line.strip() for line in f if line.strip()]
    
    logging.info(f"Processing {len(missing_species)} missing species for gene {gene}")
    
    # Check for existing files
    existing_species = check_existing_files(output_dir, missing_species)
    if existing_species:
        logging.info(f"Found existing files for {len(existing_species)} species")
    
    # Filter out existing species
    species_to_download = [sp for sp in missing_species if sp not in existing_species]
    
    stats = {
        'total_missing': len(missing_species),
        'already_exist': len(existing_species),
        'attempted': len(species_to_download),
        'downloaded': 0,
        'still_missing': 0
    }
    
    found_species = set()
    
    # Download each species individually
    for i, species in enumerate(species_to_download):
        if i > 0 and i % 10 == 0:
            logging.info(f"Progress: {i}/{len(species_to_download)} species processed")
        
        # Search for the gene in this species
        ids = search_gene_species_ncbi(gene, species)
        
        if ids:
            # Fetch sequences
            records = fetch_protein_sequences(ids)
            
            if records:
                # Save to file
                fasta_file = output_dir / f"{safe_filename(species)}.fasta"
                with open(fasta_file, "w") as f:
                    SeqIO.write(records, f, "fasta")
                
                logging.info(f"Downloaded {len(records)} sequences for {species}")
                stats['downloaded'] += 1
                found_species.add(species)
            else:
                logging.warning(f"No sequences retrieved for {species}")
        else:
            logging.debug(f"No IDs found for {species}")
        
        # Rate limiting
        time.sleep(0.4)  # Stay well below NCBI rate limits
    
    stats['still_missing'] = len(species_to_download) - stats['downloaded']
    
    # Update not_found file
    if found_species:
        update_not_found_file(not_found_file, found_species)
    
    return stats


def main():
    parser = argparse.ArgumentParser(description='Download missing proteins from NCBI')
    parser.add_argument('gene', help='Gene name (e.g., bamA)')
    parser.add_argument('not_found_file', help='Path to not_found_species.txt from UniProt')
    parser.add_argument('output_dir', help='Output directory (same as UniProt output)')
    parser.add_argument('--ncbi-config', default='config/login/ncbi_info.txt', 
                       help='Path to NCBI credentials file (default: config/login/ncbi_info.txt)')
    
    args = parser.parse_args()
    
    not_found_file = Path(args.not_found_file)
    output_dir = Path(args.output_dir)
    ncbi_config = Path(args.ncbi_config)
    
    if not not_found_file.exists():
        logging.error(f"Not found file does not exist: {not_found_file}")
        return
    
    if not ncbi_config.exists():
        logging.error(f"NCBI config file does not exist: {ncbi_config}")
        return
    
    # Read NCBI credentials from config file
    with open(ncbi_config, 'r') as f:
        lines = f.readlines()
        email = lines[0].strip()
        api_key = lines[1].strip() if len(lines) > 1 else None
    
    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Download missing proteins
    stats = download_missing_proteins(
        args.gene,
        not_found_file,
        output_dir,
        email,
        api_key
    )
    
    # Print summary
    print(f"\n=== NCBI Download Summary for {args.gene} ===")
    print(f"Total missing species: {stats['total_missing']}")
    print(f"Already had files: {stats['already_exist']}")
    print(f"Download attempted: {stats['attempted']}")
    print(f"Successfully downloaded: {stats['downloaded']}")
    print(f"Still missing: {stats['still_missing']}")


if __name__ == '__main__':
    main()