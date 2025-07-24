#!/usr/bin/env python3
"""
Optimized Protein Download Script
================================

Improvements:
- Larger batch sizes and better rate limiting
- Improved species name matching with fallbacks
- Progress tracking and estimated completion time
- Better error handling and recovery
- Parallel processing where possible
"""

from Bio import Entrez, SeqIO
import os
import time
import urllib.error
import logging
import pandas as pd
import re
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Set up logging with more detailed format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_protein_optimized.log'),
        logging.StreamHandler()
    ]
)

# Thread lock for NCBI API calls
api_lock = threading.Lock()

# ---------------- CONFIG FROM SNAKEMAKE ----------------

with open(snakemake.input.ncbi_info, "r") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

protein_file = snakemake.input.proteins
species_file = snakemake.input.species
output_base = snakemake.output.output_folder

os.makedirs(output_base, exist_ok=True)

# ---------------- LOAD INPUTS ----------------

try:
    proteins_df = pd.read_csv(protein_file, sep='\t')
    genes = proteins_df['gene'].unique().tolist()
    logging.info(f"Loaded {len(genes)} unique genes from {protein_file}")
except Exception as e:
    logging.error(f"Error reading protein file: {e}")
    raise

with open(species_file) as sf:
    species_list = [line.strip() for line in sf if line.strip()]
    logging.info(f"Loaded {len(species_list)} species from {species_file}")

# ---------------- IMPROVED FUNCTIONS ----------------

def safe_filename(name):
    """Create safe filename from species name"""
    return re.sub(r'[^\w\-_\.]', '_', name)

def normalize_species_name(species):
    """Normalize species name for better matching"""
    # Remove strain info and normalize
    normalized = re.sub(r'\s+(strain|subsp\.|spp\.|sp\.|var\.)\s+.*', '', species, flags=re.IGNORECASE)
    normalized = re.sub(r'\s+', ' ', normalized).strip()
    return normalized

def generate_search_queries(gene, species):
    """Generate multiple search query variants for better matching"""
    normalized_species = normalize_species_name(species)
    
    queries = [
        # Exact gene and species
        f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species}"[Organism]',
        # Normalized species
        f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{normalized_species}"[Organism]',
        # Genus only
        f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species.split()[0]}"[Organism]',
        # Gene symbol variations
        f'"{gene.upper()}"[Gene] AND "{normalized_species}"[Organism]',
        f'"{gene.lower()}"[Gene] AND "{normalized_species}"[Organism]'
    ]
    
    return queries

def search_gene_protein_optimized(gene, species_batch, retries=2):
    """Optimized search with better query construction and error handling"""
    
    # Try batch search first
    batch_queries = []
    for species in species_batch:
        queries = generate_search_queries(gene, species)
        batch_queries.extend(queries[:2])  # Use first 2 variants
    
    full_query = " OR ".join(f"({q})" for q in batch_queries)
    
    with api_lock:  # Ensure thread safety for NCBI API
        attempt = 0
        while attempt < retries:
            try:
                handle = Entrez.esearch(db="protein", term=full_query, retmax=2000)
                record = Entrez.read(handle)
                handle.close()
                return record["IdList"]
            except urllib.error.HTTPError as e:
                if e.code == 429:  # Rate limit
                    wait_time = 2 ** attempt
                    logging.warning(f"Rate limited, waiting {wait_time}s")
                    time.sleep(wait_time)
                else:
                    logging.warning(f"HTTP error for batch query {gene} attempt {attempt+1}: {e}")
                attempt += 1
                time.sleep(1)
            except Exception as e:
                logging.warning(f"Search error for {gene} attempt {attempt+1}: {e}")
                attempt += 1
                time.sleep(1)
    
    return []

def fetch_protein_sequences_optimized(id_list, retries=2):
    """Optimized fetch with better error handling"""
    if not id_list:
        return []
    
    with api_lock:  # Ensure thread safety for NCBI API
        attempt = 0
        while attempt < retries:
            try:
                # Split large requests
                chunk_size = 200
                all_records = []
                
                for i in range(0, len(id_list), chunk_size):
                    chunk = id_list[i:i + chunk_size]
                    handle = Entrez.efetch(
                        db="protein", 
                        id=",".join(chunk), 
                        rettype="fasta", 
                        retmode="text"
                    )
                    records = list(SeqIO.parse(handle, "fasta"))
                    handle.close()
                    all_records.extend(records)
                    
                    if i + chunk_size < len(id_list):
                        time.sleep(0.1)  # Brief pause between chunks
                
                return all_records
                
            except Exception as e:
                logging.warning(f"Fetch error attempt {attempt+1}: {e}")
                attempt += 1
                time.sleep(2)
    
    return []

def improved_species_matching(records, species_list):
    """Improved species matching with multiple strategies"""
    species_to_records = {species: [] for species in species_list}
    
    for record in records:
        description = record.description.lower()
        matched = False
        
        for species in species_list:
            if matched:
                break
                
            # Strategy 1: Exact species match
            if species.lower() in description:
                species_to_records[species].append(record)
                matched = True
                continue
            
            # Strategy 2: Normalized species match
            normalized = normalize_species_name(species).lower()
            if normalized in description:
                species_to_records[species].append(record)
                matched = True
                continue
            
            # Strategy 3: Genus match
            genus = species.split()[0].lower()
            if len(genus) > 3 and genus in description:
                # Additional check to avoid false positives
                species_words = species.lower().split()
                if len(species_words) >= 2:
                    if species_words[0] in description and species_words[1] in description:
                        species_to_records[species].append(record)
                        matched = True
                        continue
    
    return species_to_records

def check_existing_files(gene_dir, species_list):
    """Check which species already have FASTA files for this gene"""
    existing_species = set()
    if os.path.exists(gene_dir):
        for species in species_list:
            fasta_file = os.path.join(gene_dir, f"{safe_filename(species)}.fasta")
            if os.path.exists(fasta_file) and os.path.getsize(fasta_file) > 0:
                existing_species.add(species)
    return existing_species

# ---------------- MAIN PROCESSING ----------------

def process_gene(gene_idx, gene, species_list, output_base):
    """Process a single gene with progress tracking"""
    start_time = time.time()
    
    gene_dir = os.path.join(output_base, gene)
    os.makedirs(gene_dir, exist_ok=True)
    
    # Check existing files
    existing_species = check_existing_files(gene_dir, species_list)
    remaining_species = [s for s in species_list if s not in existing_species]
    
    if not remaining_species:
        logging.info(f"[{gene_idx+1}/{len(genes)}] Gene '{gene}': All {len(species_list)} species already downloaded")
        return {
            'gene': gene,
            'total_species': len(species_list),
            'existing_species': len(existing_species),
            'downloaded_species': 0,
            'success': True,
            'duration': time.time() - start_time
        }
    
    logging.info(f"[{gene_idx+1}/{len(genes)}] Processing gene '{gene}': {len(remaining_species)}/{len(species_list)} species to download")
    
    # Larger batch size for better efficiency
    batch_size = 25
    total_downloaded = 0
    
    for i in range(0, len(remaining_species), batch_size):
        batch = remaining_species[i:i + batch_size]
        batch_start = time.time()
        
        logging.info(f"  Processing batch {i//batch_size + 1}/{(len(remaining_species)-1)//batch_size + 1} ({len(batch)} species)")
        
        # Search for proteins
        ids = search_gene_protein_optimized(gene, batch)
        if not ids:
            logging.warning(f"  No IDs found for gene '{gene}' batch")
            continue
        
        # Fetch sequences
        records = fetch_protein_sequences_optimized(ids)
        if not records:
            logging.warning(f"  No sequences retrieved for gene '{gene}' batch")
            continue
        
        # Match to species
        species_to_records = improved_species_matching(records, batch)
        
        # Save files
        batch_downloaded = 0
        for species, recs in species_to_records.items():
            fasta_out = os.path.join(gene_dir, f"{safe_filename(species)}.fasta")
            if recs and not (os.path.exists(fasta_out) and os.path.getsize(fasta_out) > 0):
                with open(fasta_out, "w") as f:
                    SeqIO.write(recs, f, "fasta")
                batch_downloaded += 1
                total_downloaded += 1
        
        batch_time = time.time() - batch_start
        logging.info(f"  Batch completed: {batch_downloaded}/{len(batch)} species downloaded in {batch_time:.1f}s")
        
        # Adaptive delay based on performance
        if batch_time < 1.0:
            time.sleep(0.1)  # Very fast, minimal delay
        elif batch_time < 3.0:
            time.sleep(0.2)  # Normal speed
        else:
            time.sleep(0.5)  # Slow, give API more time
    
    gene_time = time.time() - start_time
    logging.info(f"[{gene_idx+1}/{len(genes)}] Gene '{gene}' completed: {total_downloaded} new downloads in {gene_time:.1f}s")
    
    return {
        'gene': gene,
        'total_species': len(species_list),
        'existing_species': len(existing_species),
        'downloaded_species': total_downloaded,
        'success': total_downloaded > 0 or len(existing_species) > 0,
        'duration': gene_time
    }

# ---------------- MAIN EXECUTION ----------------

# Log run parameters
try:
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    group = snakemake.params.group
    logging.info(f"Processing analysis={analysis}, paramset={paramset}, gram={group}")
except:
    logging.info("Running without Snakemake parameters")

logging.info(f"Starting optimized download for {len(genes)} genes across {len(species_list)} species")
logging.info(f"Output directory: {output_base}")

start_time = time.time()
results = []

# Process genes sequentially for now (can be parallelized if needed)
for i, gene in enumerate(genes):
    result = process_gene(i, gene, species_list, output_base)
    results.append(result)
    
    # Progress estimation
    elapsed = time.time() - start_time
    avg_time_per_gene = elapsed / (i + 1)
    remaining_genes = len(genes) - (i + 1)
    estimated_remaining = avg_time_per_gene * remaining_genes
    
    if remaining_genes > 0:
        eta = datetime.now() + timedelta(seconds=estimated_remaining)
        logging.info(f"Progress: {i+1}/{len(genes)} genes completed. ETA: {eta.strftime('%H:%M:%S')}")

# ---------------- FINAL SUMMARY ----------------
total_time = time.time() - start_time
successful_genes = sum(1 for r in results if r['success'])
total_downloaded = sum(r['downloaded_species'] for r in results)
total_existing = sum(r['existing_species'] for r in results)

logging.info(f"\n=== OPTIMIZED DOWNLOAD SUMMARY ===")
logging.info(f"Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
logging.info(f"Genes processed: {len(genes)}")
logging.info(f"Successful genes: {successful_genes}/{len(genes)}")
logging.info(f"Species files downloaded: {total_downloaded}")
logging.info(f"Species files already existed: {total_existing}")
logging.info(f"Average time per gene: {total_time/len(genes):.1f}s")

all_successful = successful_genes == len(genes)

# ---------------- WRITE COMPLETE FLAG ----------------
if all_successful:
    with open(snakemake.output.complete_flag, 'w') as f:
        f.write('Download complete\n')
        f.write(f'Total time: {total_time:.1f} seconds\n')
        f.write(f'Genes processed: {len(genes)}\n')
        f.write(f'Successful genes: {successful_genes}\n')
        f.write(f'Species files downloaded: {total_downloaded}\n')
        f.write(f'Species files already existed: {total_existing}\n')
    logging.info("Download complete flag written successfully.")
else:
    logging.warning("Some downloads failed; no complete flag written.")
    if os.path.exists(snakemake.output.complete_flag):
        os.remove(snakemake.output.complete_flag)