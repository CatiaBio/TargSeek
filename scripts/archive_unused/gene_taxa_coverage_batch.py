#!/usr/bin/env python3
"""
Batch gene-taxa coverage check using NCBI E-utilities efficiently
"""

import warnings
warnings.filterwarnings("ignore", message=".*Signature.*longdouble.*")

import pandas as pd
from Bio import Entrez
import time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from urllib.error import HTTPError

# Read input files
with open(snakemake.input.species_list, 'r') as f:
    species_list = [line.strip() for line in f if line.strip()]

with open(snakemake.input.gene_list, 'r') as f:
    gene_list = [line.strip() for line in f if line.strip()]

# Set up Entrez
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

# Rate limiting setup
rate_limiter = threading.Semaphore(10 if api_key else 3)
last_request_time = threading.Lock()
min_delay = 0.1 if api_key else 0.34

def rate_limited_search(query):
    """Perform rate-limited search"""
    with rate_limiter:
        with last_request_time:
            time.sleep(min_delay)
        
        try:
            handle = Entrez.esearch(db="protein", term=query, retmax=0)
            record = Entrez.read(handle)
            handle.close()
            return int(record["Count"])
        except HTTPError as e:
            if e.code == 429:
                time.sleep(2)  # Back off on rate limit
                return rate_limited_search(query)  # Retry
            else:
                return 0
        except Exception:
            return 0

def batch_search_species(species, genes):
    """Search all genes for a single species in batch"""
    results = {}
    
    # Build queries for all genes
    queries = []
    for gene in genes:
        query = f'("{gene}"[Gene Name] OR "{gene}"[All Fields]) AND "{species}"[Organism]'
        queries.append((gene, query))
    
    # Process in smaller batches to avoid timeouts
    batch_size = 20
    for i in range(0, len(queries), batch_size):
        batch = queries[i:i+batch_size]
        
        for gene, query in batch:
            count = rate_limited_search(query)
            results[gene] = count
    
    return species, results

# Use parallel processing for different species
print(f"\nStarting batch analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Processing {len(gene_list)} genes across {len(species_list)} species...")

# Determine number of workers based on API key
max_workers = 5 if api_key else 2

results_dict = {gene: {} for gene in gene_list}
start_time = time.time()

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    # Submit all species for processing
    future_to_species = {
        executor.submit(batch_search_species, species, gene_list): species 
        for species in species_list
    }
    
    completed = 0
    for future in as_completed(future_to_species):
        species, gene_counts = future.result()
        completed += 1
        
        # Store results
        for gene, count in gene_counts.items():
            results_dict[gene][species] = count
        
        # Progress update
        elapsed = time.time() - start_time
        rate = completed / elapsed if elapsed > 0 else 0
        eta = (len(species_list) - completed) / rate if rate > 0 else 0
        print(f"Progress: {completed}/{len(species_list)} species completed ({completed/len(species_list)*100:.1f}%) - ETA: {eta/60:.1f} min")

# Create DataFrame
data_for_df = []
for gene in gene_list:
    row = {'gene': gene}
    row.update(results_dict[gene])
    data_for_df.append(row)

df = pd.DataFrame(data_for_df)

# Calculate coverage statistics
coverage_stats = []
for idx, row in df.iterrows():
    gene = row['gene']
    # Count non-zero entries (excluding the 'gene' column)
    species_with_gene = sum(1 for species in species_list if row.get(species, 0) > 0)
    coverage_percentage = (species_with_gene / len(species_list)) * 100
    
    coverage_stats.append({
        'gene': gene,
        'species_with_gene': species_with_gene,
        'total_species': len(species_list),
        'coverage_percentage': coverage_percentage
    })

# Create coverage summary DataFrame
coverage_df = pd.DataFrame(coverage_stats)

# Sort by coverage
coverage_df = coverage_df.sort_values('coverage_percentage', ascending=False)

# Save results
coverage_df.to_csv(snakemake.output[0], sep='\t', index=False)

# Print summary
print(f"\n{'='*60}")
print(f"Analysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total time: {(time.time() - start_time)/60:.1f} minutes")
print(f"\nResults saved to: {snakemake.output[0]}")
print(f"\nTop 10 genes by coverage:")
print(coverage_df.head(10).to_string(index=False))