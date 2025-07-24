#!/usr/bin/env python3
"""
Download protein sequences from UniProt using gene name and species list.
This script queries UniProt with a single query containing all species for a given gene.
"""

import argparse
import requests
import time
from pathlib import Path
from typing import List, Dict, Any
import json


def build_uniprot_query(gene_name: str, species_list: List[str]) -> str:
    """
    Build UniProt query string for a gene across multiple species.
    
    Args:
        gene_name: Gene name (e.g., 'bamA', 'eno')
        species_list: List of species names
        
    Returns:
        UniProt query string
    """
    # Create organism query part
    organism_query = ' OR '.join([f'organism_name:"{species}"' for species in species_list])
    
    # Combine with gene name
    query = f'gene:{gene_name} AND ({organism_query})'
    
    return query


def search_uniprot(query: str, max_results: int = 500) -> List[Dict[str, Any]]:
    """
    Search UniProt and return list of accession IDs.
    
    Args:
        query: UniProt query string
        max_results: Maximum number of results to retrieve
        
    Returns:
        List of UniProt accession IDs
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    params = {
        'query': query,
        'format': 'json',
        'size': max_results,
        'fields': 'accession,organism_name,gene_names,protein_name,sequence'
    }
    
    print(f"Searching UniProt with query: {query[:100]}...")
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        
        data = response.json()
        results = data.get('results', [])
        
        print(f"Found {len(results)} entries")
        return results
        
    except requests.exceptions.RequestException as e:
        print(f"Error querying UniProt: {e}")
        return []


def download_sequences_batch(gene_name: str, species_list: List[str], output_dir: Path, max_results: int = 500) -> Dict[str, Any]:
    """
    Download all sequences for a gene across multiple species in one query.
    
    Args:
        gene_name: Gene name
        species_list: List of species names
        output_dir: Directory to save FASTA files
        
    Returns:
        Dictionary with download statistics
    """
    # Create output directory
    gene_dir = output_dir / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    # Build and execute query
    query = build_uniprot_query(gene_name, species_list)
    results = search_uniprot(query, max_results)
    
    # Process results
    species_proteins = {}
    stats = {
        'total_entries': len(results),
        'species_found': 0,
        'species_missing': [],
        'multiple_entries': {}
    }
    
    # Group by species
    for entry in results:
        organism = entry.get('organism', {}).get('scientificName', 'Unknown')
        
        if organism not in species_proteins:
            species_proteins[organism] = []
        
        species_proteins[organism].append(entry)
    
    # Save sequences for each species
    for species in species_list:
        if species in species_proteins:
            entries = species_proteins[species]
            stats['species_found'] += 1
            
            if len(entries) > 1:
                stats['multiple_entries'][species] = len(entries)
            
            # Save all sequences for this species
            fasta_file = gene_dir / f"{species.replace(' ', '_')}.fasta"
            with open(fasta_file, 'w') as f:
                for i, entry in enumerate(entries):
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
            
            print(f"  {species}: {len(entries)} sequence(s) saved")
        else:
            stats['species_missing'].append(species)
    
    # Save not found species list
    not_found_file = gene_dir / "not_found_species.txt"
    with open(not_found_file, 'w') as f:
        for species in stats['species_missing']:
            f.write(f"{species}\n")
    
    # Save summary
    summary_file = gene_dir / "download_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    return stats


def main():
    parser = argparse.ArgumentParser(description='Download proteins from UniProt using batch queries')
    parser.add_argument('gene', help='Gene name (e.g., bamA, eno)')
    parser.add_argument('species_file', help='File containing species names (one per line)')
    parser.add_argument('output_dir', help='Output directory for FASTA files')
    parser.add_argument('--max-results', type=int, default=500, help='Maximum results per query (default: 500)')
    
    args = parser.parse_args()
    
    # Read species list
    species_list = []
    with open(args.species_file, 'r') as f:
        for line in f:
            species = line.strip()
            if species:
                species_list.append(species)
    
    print(f"Downloading {args.gene} sequences for {len(species_list)} species")
    
    # Download sequences
    output_dir = Path(args.output_dir)
    stats = download_sequences_batch(args.gene, species_list, output_dir, args.max_results)
    
    # Print summary
    print(f"\nDownload complete!")
    print(f"Total entries retrieved: {stats['total_entries']}")
    print(f"Species with sequences: {stats['species_found']}/{len(species_list)}")
    print(f"Species missing: {len(stats['species_missing'])}")
    
    if stats['species_missing']:
        print("\nMissing species:")
        for species in stats['species_missing'][:10]:
            print(f"  - {species}")
        if len(stats['species_missing']) > 10:
            print(f"  ... and {len(stats['species_missing']) - 10} more")
    
    if stats['multiple_entries']:
        print("\nSpecies with multiple entries:")
        for species, count in list(stats['multiple_entries'].items())[:5]:
            print(f"  - {species}: {count} entries")


if __name__ == '__main__':
    main()