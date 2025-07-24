#!/usr/bin/env python3
"""
Download Proteins Using Gene-Specific Species Lists
==================================================

This script downloads proteins using the pre-generated gene-specific species lists
from the proteins_to_download directory. This approach is more organized and 
allows for better parallelization and debugging.

Input: 
- results/proteins_to_download/{analysis}_{paramset}_gram_{group}/*.txt
- results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv
Output: results/proteins/{analysis}_{paramset}_gram_{group}/{gene}/
"""

from Bio import Entrez, SeqIO
import os
import time
import urllib.error
import logging
import pandas as pd
import re
from datetime import datetime, timedelta
from pathlib import Path
import threading

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_protein_by_gene.log'),
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
gene_species_dir = snakemake.input.gene_species_lists
output_base = snakemake.output.output_folder

os.makedirs(output_base, exist_ok=True)

# ---------------- LOAD INPUTS ----------------

# Load genes from the protein file to ensure we only process relevant genes
try:
    proteins_df = pd.read_csv(protein_file, sep='\t')
    target_genes = set(proteins_df['gene'].unique())
    logging.info(f"Target genes from protein file: {len(target_genes)} genes")
except Exception as e:
    logging.error(f"Error reading protein file: {e}")
    raise

# Find gene-specific species list files
gene_species_path = Path(gene_species_dir)
gene_files = list(gene_species_path.glob("*.txt"))
gene_files = [f for f in gene_files if not f.name.startswith("_")]  # Skip summary files

logging.info(f"Found {len(gene_files)} gene-specific species list files")

# Filter to only target genes
available_genes = {}
for gene_file in gene_files:
    gene_name = gene_file.stem  # Remove .txt extension
    # Handle safe filename conversion back to original gene name
    original_gene = gene_name.replace("_", "/")  # Basic conversion, may need refinement
    
    if gene_name in target_genes or original_gene in target_genes:
        # Load species list for this gene
        with open(gene_file, 'r') as f:
            species_list = [line.strip() for line in f if line.strip()]
        available_genes[gene_name] = species_list
        logging.info(f"Loaded {len(species_list)} species for gene '{gene_name}'")

logging.info(f"Processing {len(available_genes)} genes with species lists")

# ---------------- IMPROVED FUNCTIONS ----------------

def safe_filename(name):
    """Create safe filename from species name"""
    return re.sub(r'[^\w\-_\.]', '_', name)

def normalize_species_name(species):
    """Normalize species name for better matching"""
    normalized = re.sub(r'\s+(strain|subsp\.|spp\.|sp\.|var\.)\s+.*', '', species, flags=re.IGNORECASE)
    normalized = re.sub(r'\s+', ' ', normalized).strip()
    return normalized

def generate_search_queries(gene, species):
    """Generate multiple search query variants for better matching"""
    normalized_species = normalize_species_name(species)
    
    queries = [
        f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species}"[Organism]',
        f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{normalized_species}"[Organism]',
        f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species.split()[0]}"[Organism]',
        f'"{gene.upper()}"[Gene] AND "{normalized_species}"[Organism]',
        f'"{gene.lower()}"[Gene] AND "{normalized_species}"[Organism]'
    ]
    
    return queries

def search_gene_protein_optimized(gene, species_batch, retries=2):
    """Search for proteins with improved query construction"""
    batch_queries = []
    for species in species_batch:
        queries = generate_search_queries(gene, species)
        batch_queries.extend(queries[:2])  # Use first 2 variants
    
    full_query = " OR ".join(f"({q})" for q in batch_queries)
    
    with api_lock:
        attempt = 0
        while attempt < retries:
            try:
                handle = Entrez.esearch(db="protein", term=full_query, retmax=2000)
                record = Entrez.read(handle)
                handle.close()
                return record["IdList"]
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    wait_time = 2 ** attempt
                    logging.warning(f"Rate limited, waiting {wait_time}s")
                    time.sleep(wait_time)
                else:
                    logging.warning(f"HTTP error for gene {gene} attempt {attempt+1}: {e}")
                attempt += 1
                time.sleep(1)
            except Exception as e:
                logging.warning(f"Search error for {gene} attempt {attempt+1}: {e}")
                attempt += 1
                time.sleep(1)
    
    return []

def fetch_protein_sequences_optimized(id_list, retries=2):
    """Fetch protein sequences with chunking for large requests"""
    if not id_list:
        return []
    
    with api_lock:
        attempt = 0
        while attempt < retries:
            try:
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
                        time.sleep(0.1)
                
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
            
            # Strategy 3: Genus match with additional verification
            genus = species.split()[0].lower()
            if len(genus) > 3 and genus in description:
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

def process_gene_with_species_list(gene_idx, gene, species_list, output_base):
    """Process a single gene using its specific species list"""
    start_time = time.time()
    
    gene_dir = os.path.join(output_base, gene)
    os.makedirs(gene_dir, exist_ok=True)
    
    # Check existing files
    existing_species = check_existing_files(gene_dir, species_list)
    remaining_species = [s for s in species_list if s not in existing_species]
    
    if not remaining_species:
        logging.info(f"[{gene_idx+1}/{len(available_genes)}] Gene '{gene}': All {len(species_list)} species already downloaded")
        return {
            'gene': gene,
            'total_species': len(species_list),
            'existing_species': len(existing_species),
            'downloaded_species': 0,
            'success': True,
            'duration': time.time() - start_time
        }
    
    logging.info(f"[{gene_idx+1}/{len(available_genes)}] Processing gene '{gene}': {len(remaining_species)}/{len(species_list)} species to download")
    
    # Process in batches
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
        
        # Adaptive delay
        if batch_time < 1.0:
            time.sleep(0.1)
        elif batch_time < 3.0:
            time.sleep(0.2)
        else:
            time.sleep(0.5)
    
    gene_time = time.time() - start_time
    logging.info(f"[{gene_idx+1}/{len(available_genes)}] Gene '{gene}' completed: {total_downloaded} new downloads in {gene_time:.1f}s")
    
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

logging.info(f"Starting gene-specific download for {len(available_genes)} genes")
logging.info(f"Output directory: {output_base}")

start_time = time.time()
results = []

# Process each gene with its specific species list
for i, (gene, species_list) in enumerate(available_genes.items()):
    result = process_gene_with_species_list(i, gene, species_list, output_base)
    results.append(result)
    
    # Progress estimation
    elapsed = time.time() - start_time
    avg_time_per_gene = elapsed / (i + 1)
    remaining_genes = len(available_genes) - (i + 1)
    estimated_remaining = avg_time_per_gene * remaining_genes
    
    if remaining_genes > 0:
        eta = datetime.now() + timedelta(seconds=estimated_remaining)
        logging.info(f"Progress: {i+1}/{len(available_genes)} genes completed. ETA: {eta.strftime('%H:%M:%S')}")

# ---------------- FINAL SUMMARY ----------------
total_time = time.time() - start_time
successful_genes = sum(1 for r in results if r['success'])
total_downloaded = sum(r['downloaded_species'] for r in results)
total_existing = sum(r['existing_species'] for r in results)

logging.info(f"\n=== GENE-SPECIFIC DOWNLOAD SUMMARY ===")
logging.info(f"Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
logging.info(f"Genes processed: {len(available_genes)}")
logging.info(f"Successful genes: {successful_genes}/{len(available_genes)}")
logging.info(f"Species files downloaded: {total_downloaded}")
logging.info(f"Species files already existed: {total_existing}")
logging.info(f"Average time per gene: {total_time/len(available_genes):.1f}s")

all_successful = successful_genes == len(available_genes)

# ---------------- WRITE COMPLETE FLAG ----------------
if all_successful:
    with open(snakemake.output.complete_flag, 'w') as f:
        f.write('Download complete\n')
        f.write(f'Total time: {total_time:.1f} seconds\n')
        f.write(f'Genes processed: {len(available_genes)}\n')
        f.write(f'Successful genes: {successful_genes}\n')
        f.write(f'Species files downloaded: {total_downloaded}\n')
        f.write(f'Species files already existed: {total_existing}\n')
    logging.info("Download complete flag written successfully.")
else:
    logging.warning("Some downloads failed; no complete flag written.")
    if os.path.exists(snakemake.output.complete_flag):
        os.remove(snakemake.output.complete_flag)