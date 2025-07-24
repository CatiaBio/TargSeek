#!/usr/bin/env python3
"""
Download protein sequences from UniProt using an iterative approach:
1. First try batch query with all species
2. Then query remaining species individually
"""

import argparse
import requests
import time
from pathlib import Path
from typing import List, Dict, Any, Set
import json
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def search_uniprot(query: str, max_results: int = 500) -> List[Dict[str, Any]]:
    """
    Search UniProt and return list of entries.
    
    Args:
        query: UniProt query string
        max_results: Maximum number of results to retrieve
        
    Returns:
        List of UniProt entries
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    params = {
        'query': query,
        'format': 'json',
        'size': max_results,
        'fields': 'accession,organism_name,gene_names,protein_name,sequence'
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        
        data = response.json()
        results = data.get('results', [])
        
        return results
        
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 400:
            # Bad request - likely invalid query syntax
            logging.debug(f"Bad request for query: {params.get('query', '')[:100]}...")
        else:
            logging.error(f"HTTP error querying UniProt: {e}")
        return []
    except requests.exceptions.RequestException as e:
        logging.error(f"Error querying UniProt: {e}")
        return []


def search_gene_species_individual(gene_name: str, species: str, max_results: int = 50) -> List[Dict[str, Any]]:
    """
    Search UniProt for a specific gene in a specific species.
    
    Args:
        gene_name: Gene name
        species: Species name
        max_results: Maximum results per species
        
    Returns:
        List of UniProt entries
    """
    try:
        # Try exact species name first
        query = f'gene:{gene_name} AND organism_name:"{species}"'
        results = search_uniprot(query, max_results)
        
        if not results and not species.endswith(' sp.'):
            # Try without quotes for partial matching (but not for "sp." species)
            query = f'gene:{gene_name} AND organism:{species}'
            results = search_uniprot(query, max_results)
        
        return results
    except Exception as e:
        logging.debug(f"Error searching for {gene_name} in {species}: {e}")
        return []


def save_sequences(species_proteins: Dict[str, List[Dict]], gene_dir: Path, species: str) -> bool:
    """
    Save sequences for a species to FASTA file.
    
    Args:
        species_proteins: Dictionary of species to protein entries
        gene_dir: Directory to save files
        species: Species name
        
    Returns:
        True if sequences were saved, False otherwise
    """
    if species not in species_proteins or not species_proteins[species]:
        return False
    
    entries = species_proteins[species]
    fasta_file = gene_dir / f"{species.replace(' ', '_')}.fasta"
    
    with open(fasta_file, 'w') as f:
        for entry in entries:
            accession = entry.get('primaryAccession', 'Unknown')
            gene_names = entry.get('genes', [])
            gene_name_str = gene_names[0].get('geneName', {}).get('value', '') if gene_names else ''
            protein_name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown protein')
            sequence = entry.get('sequence', {}).get('value', '')
            
            # Write FASTA header and sequence
            header = f">{accession}|{gene_name_str}|{species}|{protein_name}"
            f.write(f"{header}\n")
            # Write sequence in 60-character lines
            for j in range(0, len(sequence), 60):
                f.write(f"{sequence[j:j+60]}\n")
    
    return True


def download_with_retry(gene_name: str, species_list: List[str], output_dir: Path, max_results: int = 500) -> Dict[str, Any]:
    """
    Download protein sequences using batch query first, then individual queries for missing species.
    
    Args:
        gene_name: Gene name
        species_list: List of species names
        output_dir: Directory to save FASTA files
        max_results: Maximum results for batch query
        
    Returns:
        Dictionary with download statistics
    """
    # Create output directory
    gene_dir = output_dir / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    # Check for existing files
    existing_species = set()
    for species in species_list:
        fasta_file = gene_dir / f"{species.replace(' ', '_')}.fasta"
        if fasta_file.exists() and fasta_file.stat().st_size > 0:
            existing_species.add(species)
    
    if existing_species:
        logging.info(f"Found existing files for {len(existing_species)} species")
    
    # Filter out existing species
    species_to_download = [s for s in species_list if s not in existing_species]
    
    stats = {
        'total_species': len(species_list),
        'already_existed': len(existing_species),
        'batch_found': 0,
        'individual_found': 0,
        'not_found': 0,
        'species_missing': []
    }
    
    if not species_to_download:
        logging.info("All species already downloaded")
        return stats
    
    # Step 1: Batch queries in groups of 10 species
    logging.info(f"Step 1: Batch queries for {len(species_to_download)} species (batches of 10)")
    
    batch_found = set()
    batch_size = 10
    
    for i in range(0, len(species_to_download), batch_size):
        batch_species = species_to_download[i:i+batch_size]
        batch_num = (i // batch_size) + 1
        total_batches = (len(species_to_download) + batch_size - 1) // batch_size
        
        logging.info(f"  Processing batch {batch_num}/{total_batches}: {len(batch_species)} species")
        
        # Create batch query for this subset
        organism_query = ' OR '.join([f'organism_name:"{species}"' for species in batch_species])
        batch_query = f'gene:{gene_name} AND ({organism_query})'
        
        batch_results = search_uniprot(batch_query, max_results)
        logging.info(f"  Batch {batch_num} returned {len(batch_results)} entries")
        
        # Process batch results
        species_proteins = {}
        for entry in batch_results:
            organism = entry.get('organism', {}).get('scientificName', 'Unknown')
            if organism not in species_proteins:
                species_proteins[organism] = []
            species_proteins[organism].append(entry)
        
        # Save sequences from this batch
        for species in batch_species:
            if save_sequences(species_proteins, gene_dir, species):
                batch_found.add(species)
                logging.info(f"    Batch {batch_num}: Found {len(species_proteins[species])} sequences for {species}")
        
        # Rate limiting between batches
        if i + batch_size < len(species_to_download):
            time.sleep(0.5)
    
    stats['batch_found'] = len(batch_found)
    logging.info(f"Batch queries completed: found {len(batch_found)} species total")
    
    # Step 2: Individual queries for missing species
    missing_species = [s for s in species_to_download if s not in batch_found]
    
    if missing_species:
        logging.info(f"\nStep 2: Individual queries for {len(missing_species)} missing species")
        
        for i, species in enumerate(missing_species):
            if i > 0 and i % 10 == 0:
                logging.info(f"  Progress: {i}/{len(missing_species)} species processed")
            
            # Search for this specific species
            results = search_gene_species_individual(gene_name, species)
            
            if results:
                # Process results
                species_proteins_individual = {species: results}
                if save_sequences(species_proteins_individual, gene_dir, species):
                    stats['individual_found'] += 1
                    logging.info(f"  Individual: Found {len(results)} sequences for {species}")
                else:
                    stats['species_missing'].append(species)
            else:
                stats['species_missing'].append(species)
                logging.debug(f"  Individual: No results for {species}")
            
            # Rate limiting
            time.sleep(0.2)
    
    stats['not_found'] = len(stats['species_missing'])
    
    # Save or update not found species list
    not_found_file = gene_dir / "not_found_species.txt"
    if stats['species_missing']:
        with open(not_found_file, 'w') as f:
            for species in stats['species_missing']:
                f.write(f"{species}\n")
        logging.info(f"Updated not_found_species.txt with {len(stats['species_missing'])} species")
    else:
        # Remove not_found file if all species were found
        if not_found_file.exists():
            not_found_file.unlink()
            logging.info("Removed not_found_species.txt - all species found!")
    
    # Save summary
    summary_file = gene_dir / "download_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    return stats


def main():
    parser = argparse.ArgumentParser(description='Download proteins from UniProt with iterative retry')
    parser.add_argument('gene', help='Gene name (e.g., bamA, eno)')
    parser.add_argument('species_file', help='File containing species names (one per line)')
    parser.add_argument('output_dir', help='Output directory for FASTA files')
    parser.add_argument('--max-results', type=int, default=500, help='Maximum results for batch query (default: 500)')
    parser.add_argument('--max-per-species', type=int, default=50, help='Maximum results per species for individual queries (default: 50)')
    
    args = parser.parse_args()
    
    # Read species list
    species_list = []
    with open(args.species_file, 'r') as f:
        for line in f:
            species = line.strip()
            if species:
                species_list.append(species)
    
    logging.info(f"Downloading {args.gene} sequences for {len(species_list)} species")
    
    # Download sequences
    output_dir = Path(args.output_dir)
    stats = download_with_retry(args.gene, species_list, output_dir, args.max_results)
    
    # Print summary
    print(f"\n=== Download Summary for {args.gene} ===")
    print(f"Total species requested: {stats['total_species']}")
    print(f"Already had files: {stats['already_existed']}")
    print(f"Found in batch query: {stats['batch_found']}")
    print(f"Found in individual queries: {stats['individual_found']}")
    print(f"Total downloaded: {stats['batch_found'] + stats['individual_found']}")
    print(f"Still missing: {stats['not_found']}")
    
    total_found = stats['already_existed'] + stats['batch_found'] + stats['individual_found']
    coverage = (total_found / stats['total_species']) * 100
    print(f"Coverage: {coverage:.1f}%")
    
    if stats['species_missing']:
        print(f"\nFirst 10 missing species:")
        for species in stats['species_missing'][:10]:
            print(f"  - {species}")
        if len(stats['species_missing']) > 10:
            print(f"  ... and {len(stats['species_missing']) - 10} more")


if __name__ == '__main__':
    main()