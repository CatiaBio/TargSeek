#!/usr/bin/env python3
"""
Check gene-taxa coverage with proper NCBI rate limiting
"""

import pandas as pd
from Bio import Entrez
import time
from datetime import datetime
import sys
import os
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
print(f"\nStarting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\nProcessing {len(gene_list)} genes across {len(species_list)} species...")

# Determine rate limit based on API key
if api_key:
    requests_per_second = 10  # With API key
    delay = 0.11  # Slightly more than 0.1 to be safe
else:
    requests_per_second = 3   # Without API key
    delay = 0.34  # Slightly more than 0.33 to be safe

print(f"Rate limit: {requests_per_second} requests/second (delay: {delay}s)")

# Initialize results
results = []
total_queries = len(gene_list) * len(species_list)
completed = 0

print(f"Total queries to perform: {total_queries}")

def search_protein_with_retry(gene, species, max_retries=3):
    """Search for protein with retry logic for rate limiting"""
    for attempt in range(max_retries):
        try:
            # Add delay before each request
            time.sleep(delay)
            
            # Search query
            query = f'("{gene}"[Gene Name] OR "{gene}"[All Fields]) AND "{species}"[Organism]'
            
            # Perform search
            handle = Entrez.esearch(db="protein", term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            
            return int(record["Count"])
            
        except HTTPError as e:
            if e.code == 429:  # Too Many Requests
                wait_time = 2 ** attempt  # Exponential backoff: 1, 2, 4 seconds
                print(f"\nRate limit hit for {gene}/{species}. Waiting {wait_time}s before retry {attempt+1}/{max_retries}")
                time.sleep(wait_time)
            else:
                print(f"\nHTTP Error {e.code} for {gene}/{species}: {e}")
                return 0
        except Exception as e:
            print(f"\nError querying {gene} for {species}: {str(e)}")
            return 0
    
    print(f"\nFailed after {max_retries} attempts for {gene}/{species}")
    return 0

# Process each gene-species combination
start_time = time.time()
for gene in gene_list:
    gene_results = {'gene': gene}
    
    for species in species_list:
        count = search_protein_with_retry(gene, species)
        gene_results[species] = count
        completed += 1
        
        # Progress update every 100 queries
        if completed % 100 == 0:
            elapsed = time.time() - start_time
            rate = completed / elapsed if elapsed > 0 else 0
            eta = (total_queries - completed) / rate if rate > 0 else 0
            print(f"Progress: {completed}/{total_queries} queries completed ({completed/total_queries*100:.1f}%) - Rate: {rate:.1f} q/s - ETA: {eta/60:.1f} min")
    
    results.append(gene_results)

# Create DataFrame
df = pd.DataFrame(results)

# Calculate coverage statistics
coverage_stats = []
for idx, row in df.iterrows():
    gene = row['gene']
    # Count non-zero entries (excluding the 'gene' column)
    species_with_gene = sum(1 for species in species_list if row[species] > 0)
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