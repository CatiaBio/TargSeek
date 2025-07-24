#!/usr/bin/env python3
"""
Simple Sequential Protein Download
=================================

Downloads proteins one gene at a time with exact species matching.
Faster and more focused approach:
- Process one gene at a time
- Start with gram-negative, then gram-positive  
- Use exact species name matching only
- Track not-found species in separate files
- Clear progress reporting
"""

from Bio import Entrez, SeqIO
import os
import time
import urllib.error
import logging
import pandas as pd
import re
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_protein_simple.log'),
        logging.StreamHandler()
    ]
)

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

# Load genes from the protein file
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
gene_files = [f for f in gene_files if not f.name.startswith("_")]

# Filter and prepare gene list
available_genes = {}
for gene_file in gene_files:
    gene_name = gene_file.stem
    original_gene = gene_name.replace("_", "/")
    
    if gene_name in target_genes or original_gene in target_genes:
        with open(gene_file, 'r') as f:
            species_list = [line.strip() for line in f if line.strip()]
        available_genes[gene_name] = species_list

logging.info(f"Found {len(available_genes)} genes to process")

# ---------------- SIMPLE FUNCTIONS ----------------

def safe_filename(name):
    """Create safe filename from species name"""
    return re.sub(r'[^\w\-_\.]', '_', name)

def search_gene_protein_exact(gene, species, retries=2):
    """Search for protein with exact species name"""
    query = f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species}"[Organism]'
    
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.esearch(db="protein", term=query, retmax=50)
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait_time = 2 ** attempt
                logging.warning(f"Rate limited, waiting {wait_time}s")
                time.sleep(wait_time)
            else:
                logging.warning(f"HTTP error for {gene}-{species} attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(1)
        except Exception as e:
            logging.warning(f"Search error for {gene}-{species} attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(1)
    
    return []

def fetch_protein_sequences_simple(id_list, retries=2):
    """Fetch protein sequences with simple error handling"""
    if not id_list:
        return []
    
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(
                db="protein", 
                id=",".join(id_list), 
                rettype="fasta", 
                retmode="text"
            )
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            return records
        except Exception as e:
            logging.warning(f"Fetch error attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(2)
    
    return []

def process_single_gene(gene, species_list, output_base):
    """Process a single gene with all its species"""
    start_time = time.time()
    
    gene_dir = os.path.join(output_base, gene)
    os.makedirs(gene_dir, exist_ok=True)
    
    # Check existing files
    existing_species = set()
    for species in species_list:
        fasta_file = os.path.join(gene_dir, f"{safe_filename(species)}.fasta")
        if os.path.exists(fasta_file) and os.path.getsize(fasta_file) > 0:
            existing_species.add(species)
    
    remaining_species = [s for s in species_list if s not in existing_species]
    
    if not remaining_species:
        logging.info(f"Gene '{gene}': All {len(species_list)} species already downloaded")
        return {
            'gene': gene,
            'total_species': len(species_list),
            'downloaded': 0,
            'already_existed': len(existing_species),
            'not_found': [],
            'duration': time.time() - start_time
        }
    
    logging.info(f"Gene '{gene}': Processing {len(remaining_species)}/{len(species_list)} species")
    
    downloaded_count = 0
    not_found_species = []
    
    # Process each species individually
    for i, species in enumerate(remaining_species):
        logging.info(f"  [{i+1}/{len(remaining_species)}] Searching for {gene} in {species}")
        
        # Search for this specific gene-species combination
        ids = search_gene_protein_exact(gene, species)
        
        if not ids:
            logging.warning(f"    ✗ No IDs found for {gene} in {species}")
            not_found_species.append(species)
            time.sleep(0.5)  # Brief pause between failed searches
            continue
        
        # Fetch sequences
        records = fetch_protein_sequences_simple(ids)
        
        if not records:
            logging.warning(f"    ✗ No sequences retrieved for {gene} in {species}")
            not_found_species.append(species)
            time.sleep(0.5)
            continue
        
        # Save sequences
        fasta_out = os.path.join(gene_dir, f"{safe_filename(species)}.fasta")
        with open(fasta_out, "w") as f:
            SeqIO.write(records, f, "fasta")
        
        downloaded_count += 1
        logging.info(f"    ✓ Downloaded {len(records)} sequences for {species}")
        
        # Brief pause between successful downloads
        time.sleep(0.3)
    
    # Save not-found list
    if not_found_species:
        not_found_file = os.path.join(gene_dir, "not_found_species.txt")
        with open(not_found_file, 'w') as f:
            f.write(f"Species not found for gene '{gene}':\n")
            for species in not_found_species:
                f.write(f"{species}\n")
        logging.info(f"  Saved {len(not_found_species)} not-found species to: not_found_species.txt")
    
    gene_time = time.time() - start_time
    success_rate = (downloaded_count / len(remaining_species)) * 100 if remaining_species else 100
    
    logging.info(f"Gene '{gene}' completed: {downloaded_count}/{len(remaining_species)} downloaded ({success_rate:.1f}%) in {gene_time:.1f}s")
    
    return {
        'gene': gene,
        'total_species': len(species_list),
        'downloaded': downloaded_count,
        'already_existed': len(existing_species),
        'not_found': not_found_species,
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

# Sort genes for consistent processing order
# Process gram-positive first (usually fewer genes, good for testing)
gene_names = sorted(available_genes.keys())

logging.info(f"Starting simple sequential download for {len(gene_names)} genes")
logging.info(f"Output directory: {output_base}")
logging.info(f"Processing order: {', '.join(gene_names[:5])}{'...' if len(gene_names) > 5 else ''}")

start_time = time.time()
results = []

# Process each gene sequentially
for i, gene in enumerate(gene_names):
    species_list = available_genes[gene]
    
    logging.info(f"\n{'='*60}")
    logging.info(f"GENE {i+1}/{len(gene_names)}: {gene} ({len(species_list)} species)")
    logging.info(f"{'='*60}")
    
    result = process_single_gene(gene, species_list, output_base)
    results.append(result)
    
    # Progress summary
    elapsed = time.time() - start_time
    avg_time_per_gene = elapsed / (i + 1)
    remaining_genes = len(gene_names) - (i + 1)
    estimated_remaining = avg_time_per_gene * remaining_genes
    
    total_downloaded = sum(r['downloaded'] for r in results)
    total_not_found = sum(len(r['not_found']) for r in results)
    
    logging.info(f"\nProgress Summary after {i+1} genes:")
    logging.info(f"  Downloaded: {total_downloaded} species files")
    logging.info(f"  Not found: {total_not_found} species")
    logging.info(f"  Time elapsed: {elapsed/60:.1f} minutes")
    logging.info(f"  Estimated remaining: {estimated_remaining/60:.1f} minutes")

# ---------------- FINAL SUMMARY ----------------
total_time = time.time() - start_time
total_downloaded = sum(r['downloaded'] for r in results)
total_existed = sum(r['already_existed'] for r in results)
total_not_found = sum(len(r['not_found']) for r in results)
total_attempted = sum(r['total_species'] for r in results)

logging.info(f"\n{'='*60}")
logging.info(f"FINAL SUMMARY")
logging.info(f"{'='*60}")
logging.info(f"Total time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
logging.info(f"Genes processed: {len(gene_names)}")
logging.info(f"Species files downloaded: {total_downloaded}")
logging.info(f"Species files already existed: {total_existed}")
logging.info(f"Species not found: {total_not_found}")
logging.info(f"Total species attempted: {total_attempted}")
logging.info(f"Overall success rate: {(total_downloaded + total_existed)/total_attempted*100:.1f}%")
logging.info(f"Average time per gene: {total_time/len(gene_names):.1f}s")

# Write completion flag
all_successful = len(results) == len(gene_names)

if all_successful:
    with open(snakemake.output.complete_flag, 'w') as f:
        f.write('Download complete (simple sequential)\n')
        f.write(f'Total time: {total_time:.1f} seconds\n')
        f.write(f'Genes processed: {len(gene_names)}\n')
        f.write(f'Species files downloaded: {total_downloaded}\n')
        f.write(f'Species files already existed: {total_existed}\n')
        f.write(f'Species not found: {total_not_found}\n')
        f.write(f'Overall success rate: {(total_downloaded + total_existed)/total_attempted*100:.1f}%\n')
    logging.info("Download complete flag written successfully.")
else:
    logging.warning("Some genes failed; no complete flag written.")
    if os.path.exists(snakemake.output.complete_flag):
        os.remove(snakemake.output.complete_flag)