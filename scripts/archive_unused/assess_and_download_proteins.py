#!/usr/bin/env python3
"""
Combined coverage assessment and protein download script.

This script replaces the separate coverage assessment and download steps by:
1. Checking gene-species coverage (like gene_taxa_coverage_cached.py)
2. Immediately downloading sequences for found combinations
3. Using UniProt iterative approach (batch + individual) first
4. Using NCBI for remaining missing species
5. Creating coverage TSV and lists of found/not found species

Input:
- Species list file
- Gene list file  
- NCBI credentials

Output:
- Coverage TSV file (compatible with existing pipeline)
- Downloaded FASTA files organized by gene
- Found/not found species lists per gene
"""

import warnings
warnings.filterwarnings("ignore", message=".*Signature.*longdouble.*")

import pandas as pd
from Bio import Entrez, SeqIO
import time
import requests
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from urllib.error import HTTPError
import os
import json
from pathlib import Path
from typing import List, Dict, Any, Set, Tuple
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Cache file location
CACHE_DIR = Path("cache/gene_species")
CACHE_FILE = CACHE_DIR / "species_gene_cache.json"

print(f"=== Starting Combined Coverage Assessment and Download ===")
print(f"Script: {__file__}")
print(f"Coverage output: {snakemake.output.coverage}")
print(f"Download output: {snakemake.output.download_dir}")
print(f"Cache file: {CACHE_FILE}")

# Verify inputs exist
for input_name, input_path in [
    ("species_list", snakemake.input.species_list),
    ("gene_list", snakemake.input.gene_list), 
    ("ncbi_info", snakemake.input.ncbi_info)
]:
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_name} = {input_path}")

# Read input files
with open(snakemake.input.species_list, 'r') as f:
    species_list = [line.strip() for line in f if line.strip()]

with open(snakemake.input.gene_list, 'r') as f:
    gene_list = [line.strip() for line in f if line.strip()]

# Set up NCBI
with open(snakemake.input.ncbi_info, 'r') as f:
    lines = f.readlines()
    email = lines[0].strip()
    api_key = lines[1].strip() if len(lines) > 1 else None

Entrez.email = email
if api_key:
    Entrez.api_key = api_key

print(f"Loaded {len(species_list)} species")
print(f"Loaded {len(gene_list)} genes")
print(f"NCBI API key: {'Configured' if api_key else 'Not found'}")

# Create output directory
output_dir = Path(snakemake.output.download_dir)
try:
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Created output directory: {output_dir}")
except Exception as e:
    raise RuntimeError(f"Failed to create output directory {output_dir}: {e}")

# Load/save cache functions (same as original)
def load_cache():
    """Load existing species-gene cache"""
    if CACHE_FILE.exists():
        try:
            with open(CACHE_FILE, 'r') as f:
                cache = json.load(f)
            print(f"✓ Loaded cache with {len(cache)} entries")
            return cache
        except Exception as e:
            print(f"Warning: Could not load cache: {e}")
            return {}
    else:
        print("No existing cache found, starting fresh")
        return {}

def save_cache(cache):
    """Save species-gene cache"""
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with open(CACHE_FILE, 'w') as f:
            json.dump(cache, f, indent=2)
        print(f"✓ Saved cache with {len(cache)} entries to {CACHE_FILE}")
    except Exception as e:
        print(f"Warning: Could not save cache: {e}")

def get_cache_key(species, gene):
    """Generate cache key for species-gene pair"""
    return f"{species}||{gene}"

# UniProt functions (from uniprot_iterative script)
def search_uniprot(query: str, max_results: int = 500) -> List[Dict[str, Any]]:
    """Search UniProt and return list of entries"""
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
            logging.debug(f"Bad request for query: {params.get('query', '')[:100]}...")
        else:
            logging.error(f"HTTP error querying UniProt: {e}")
        return []
    except requests.exceptions.RequestException as e:
        logging.error(f"Error querying UniProt: {e}")
        return []

def search_gene_species_individual_uniprot(gene_name: str, species: str, max_results: int = 50) -> List[Dict[str, Any]]:
    """Search UniProt for a specific gene in a specific species"""
    try:
        query = f'gene:{gene_name} AND organism_name:"{species}"'
        results = search_uniprot(query, max_results)
        
        if not results and not species.endswith(' sp.'):
            query = f'gene:{gene_name} AND organism:{species}'
            results = search_uniprot(query, max_results)
        
        return results
    except Exception as e:
        logging.debug(f"Error searching UniProt for {gene_name} in {species}: {e}")
        return []

# NCBI functions
def search_gene_species_ncbi(gene: str, species: str, retries: int = 3) -> List[str]:
    """Search NCBI for a specific gene in a specific species"""
    query = f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species}"[Organism]'
    
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.esearch(db="protein", term=query, retmax=10)
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        except Exception as e:
            logging.debug(f"NCBI search error for {gene}-{species} attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(2 ** attempt)
    return []

def fetch_protein_sequences_ncbi(id_list: List[str], retries: int = 3) -> List:
    """Fetch protein sequences from NCBI by ID"""
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
            logging.debug(f"NCBI fetch error attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(2 ** attempt)
    return []

def save_sequences_to_fasta(sequences: List, gene_dir: Path, species: str, source: str = ""):
    """Save sequences to FASTA file"""
    if not sequences:
        return False
    
    fasta_file = gene_dir / f"{species.replace(' ', '_')}.fasta"
    
    with open(fasta_file, 'w') as f:
        if source == "uniprot":
            # UniProt format
            for entry in sequences:
                accession = entry.get('primaryAccession', 'Unknown')
                gene_names = entry.get('genes', [])
                gene_name_str = gene_names[0].get('geneName', {}).get('value', '') if gene_names else ''
                protein_name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown protein')
                sequence = entry.get('sequence', {}).get('value', '')
                
                header = f">{accession}|{gene_name_str}|{species}|{protein_name}"
                f.write(f"{header}\n")
                for j in range(0, len(sequence), 60):
                    f.write(f"{sequence[j:j+60]}\n")
        else:
            # NCBI format (SeqRecord objects)
            SeqIO.write(sequences, f, "fasta")
    
    return True

def download_gene_for_all_species(gene: str, species_list: List[str], output_dir: Path) -> Tuple[List[str], List[str], Dict[str, int]]:
    """
    Download a gene for all species using the iterative UniProt approach.
    This is similar to the uniprot_iterative script but integrated here.
    
    Returns:
        (found_species, not_found_species, source_stats)
    """
    gene_dir = output_dir / gene
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    # Check for existing files first
    existing_species = set()
    for species in species_list:
        fasta_file = gene_dir / f"{species.replace(' ', '_')}.fasta"
        if fasta_file.exists() and fasta_file.stat().st_size > 0:
            existing_species.add(species)
    
    if existing_species:
        logging.info(f"  Found existing files for {len(existing_species)} species")
    
    # Filter out existing species
    species_to_download = [s for s in species_list if s not in existing_species]
    
    source_stats = {"existing_file": len(existing_species), "uniprot_batch": 0, "uniprot_individual": 0, "ncbi": 0, "not_found": 0}
    
    if not species_to_download:
        logging.info(f"  All species already downloaded for {gene}")
        return list(existing_species), [], source_stats
    
    # Step 1: UniProt batch query
    logging.info(f"  Step 1: UniProt batch query for {len(species_to_download)} species")
    
    organism_query = ' OR '.join([f'organism_name:"{species}"' for species in species_to_download])
    batch_query = f'gene:{gene} AND ({organism_query})'
    
    batch_results = search_uniprot(batch_query, 500)
    logging.info(f"  Batch query returned {len(batch_results)} entries")
    
    # Process batch results
    species_proteins = {}
    for entry in batch_results:
        organism = entry.get('organism', {}).get('scientificName', 'Unknown')
        if organism not in species_proteins:
            species_proteins[organism] = []
        species_proteins[organism].append(entry)
    
    # Save sequences from batch query
    batch_found = set()
    for species in species_to_download:
        if species in species_proteins:
            if save_sequences_to_fasta(species_proteins[species], gene_dir, species, "uniprot"):
                batch_found.add(species)
                logging.info(f"    Batch: Found {len(species_proteins[species])} sequences for {species}")
    
    source_stats["uniprot_batch"] = len(batch_found)
    
    # Step 2: Individual UniProt queries for missing species
    missing_after_batch = [s for s in species_to_download if s not in batch_found]
    
    if missing_after_batch:
        logging.info(f"  Step 2: Individual UniProt queries for {len(missing_after_batch)} missing species")
        
        individual_found = set()
        for i, species in enumerate(missing_after_batch):
            if i > 0 and i % 10 == 0:
                logging.info(f"    Progress: {i}/{len(missing_after_batch)} species")
            
            results = search_gene_species_individual_uniprot(gene, species)
            if results:
                if save_sequences_to_fasta(results, gene_dir, species, "uniprot"):
                    individual_found.add(species)
                    logging.info(f"    Individual: Found {len(results)} sequences for {species}")
            
            time.sleep(0.2)  # Rate limiting
        
        source_stats["uniprot_individual"] = len(individual_found)
    else:
        individual_found = set()
    
    # Step 3: NCBI for remaining missing species
    missing_after_uniprot = [s for s in missing_after_batch if s not in individual_found]
    
    if missing_after_uniprot:
        logging.info(f"  Step 3: NCBI queries for {len(missing_after_uniprot)} remaining species")
        
        ncbi_found = set()
        for i, species in enumerate(missing_after_uniprot):
            if i > 0 and i % 10 == 0:
                logging.info(f"    Progress: {i}/{len(missing_after_uniprot)} species")
            
            ncbi_ids = search_gene_species_ncbi(gene, species)
            if ncbi_ids:
                ncbi_sequences = fetch_protein_sequences_ncbi(ncbi_ids)
                if ncbi_sequences:
                    if save_sequences_to_fasta(ncbi_sequences, gene_dir, species, "ncbi"):
                        ncbi_found.add(species)
                        logging.info(f"    NCBI: Found {len(ncbi_sequences)} sequences for {species}")
            
            time.sleep(0.4)  # Slower rate limiting for NCBI
        
        source_stats["ncbi"] = len(ncbi_found)
    else:
        ncbi_found = set()
    
    # Calculate final results
    all_found = existing_species | batch_found | individual_found | ncbi_found
    not_found = [s for s in species_list if s not in all_found]
    source_stats["not_found"] = len(not_found)
    
    return list(all_found), not_found, source_stats

# Main processing - Download FASTA first, then create statistics
try:
    print(f"Processing {len(gene_list)} genes for {len(species_list)} species")
    
    coverage_data = []
    gene_species_stats = {}
    
    for gene_idx, gene in enumerate(gene_list):
        print(f"\n=== Processing gene {gene_idx + 1}/{len(gene_list)}: {gene} ===")
        
        # Download sequences for this gene across all species
        found_species, not_found_species, source_stats = download_gene_for_all_species(gene, species_list, output_dir)
        
        # Save found/not found species lists
        gene_dir = output_dir / gene
        
        # Save found species list
        with open(gene_dir / "found_species.txt", 'w') as f:
            for species in found_species:
                f.write(f"{species}\n")
        
        # Save not found species list
        if not_found_species:
            with open(gene_dir / "not_found_species.txt", 'w') as f:
                for species in not_found_species:
                    f.write(f"{species}\n")
        else:
            # Create empty file if all species were found
            (gene_dir / "not_found_species.txt").touch()
        
        # Create coverage data for TSV
        for species in species_list:
            found = species in found_species
            coverage_data.append({
                "gene": gene,
                "species": species,
                "found": found
            })
        
        # Store statistics
        gene_species_stats[gene] = {
            "found": len(found_species),
            "not_found": len(not_found_species),
            "source_stats": source_stats
        }
        
        coverage_pct = (len(found_species) / len(species_list)) * 100
        print(f"  {gene}: {len(found_species)}/{len(species_list)} species ({coverage_pct:.1f}%)")
        print(f"    Sources: Existing={source_stats['existing_file']}, UniProt_batch={source_stats['uniprot_batch']}, UniProt_individual={source_stats['uniprot_individual']}, NCBI={source_stats['ncbi']}, Not_found={source_stats['not_found']}")

except Exception as e:
    print(f"Error in main processing: {e}")
    raise

# Note: We're not using cache in this approach since we download and save files directly

# Create coverage TSV (compatible with existing pipeline)
coverage_df = pd.DataFrame(coverage_data)
pivot_df = coverage_df.pivot(index='species', columns='gene', values='found').fillna(False)

# Convert to integer (1/0) for compatibility
pivot_df = pivot_df.astype(int)

# Save coverage TSV
pivot_df.to_csv(snakemake.output.coverage, sep='\t')

# Save summary statistics
summary_stats = {
    "total_combinations": len(coverage_data),
    "genes_processed": len(gene_list),
    "species_processed": len(species_list),
    "total_found": sum(stats["found"] for stats in gene_species_stats.values()),
    "total_downloaded": sum(
        stats["source_stats"]["uniprot_batch"] + 
        stats["source_stats"]["uniprot_individual"] + 
        stats["source_stats"]["ncbi"] 
        for stats in gene_species_stats.values()
    ),
    "gene_stats": gene_species_stats,
    "processing_date": datetime.now().isoformat()
}

summary_file = output_dir / "processing_summary.json"
with open(summary_file, 'w') as f:
    json.dump(summary_stats, f, indent=2)

print(f"\n=== Processing Complete ===")
print(f"Coverage TSV saved to: {snakemake.output.coverage}")
print(f"Downloads saved to: {snakemake.output.download_dir}")
print(f"Total combinations processed: {len(coverage_data)}")
print(f"Total sequences downloaded: {summary_stats['total_downloaded']}")