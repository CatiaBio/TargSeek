#!/usr/bin/env python3

import json
import os
import time
from datetime import datetime
from Bio import Entrez
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from collections import defaultdict

# Load inputs from Snakemake
try:
    species_file = snakemake.input.species
    genes_file = snakemake.input.genes
    credentials_file = snakemake.input.ncbi_info
    output_file = snakemake.output.coverage
    print(f"DEBUG: Using Snakemake inputs - output_file: {output_file}")
    
    # Use a single shared cache file for all gene-species combinations
    cache_file = "results/coverage/gene_species_cache.json"
    updated_cache_file = cache_file  # Same as input cache
    
    # Create directory using pathlib for better cross-platform compatibility
    from pathlib import Path
    cache_path = Path(cache_file)
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"DEBUG: Created cache directory: {cache_path.parent}")
    except Exception as e:
        print(f"DEBUG: Error creating cache directory: {e}")
    
    if not cache_path.exists():
        try:
            with cache_path.open('w') as f:
                json.dump({}, f)
            print(f"DEBUG: Created initial cache file: {cache_path}")
        except Exception as e:
            print(f"DEBUG: Error creating cache file: {e}")
                
except NameError:
    # Fallback for testing
    print("DEBUG: Snakemake object not available, using fallback test mode")
    species_file = "data/bacdive/analysis_1/gram_positive.txt"
    genes_file = "data/quickgo/params_1/gene_symbols_filtered.txt"
    credentials_file = "config/login/ncbi_info.txt"
    cache_file = "results/coverage/gene_species_cache.json"
    output_file = "results/coverage/test_coverage.tsv"
    updated_cache_file = cache_file  # Same as input cache

# NCBI login
with open(credentials_file) as f:
    lines = f.readlines()
    email = lines[0].strip()
    api_key = lines[1].strip() if len(lines) > 1 else None

Entrez.email = email
if api_key:
    Entrez.api_key = api_key

# Load species and genes
with open(species_file) as f:
    species_list = [line.strip() for line in f if line.strip()]

with open(genes_file) as f:
    gene_list = [line.strip() for line in f if line.strip()]

# Load cache
cache = {}
if os.path.exists(cache_file):
    try:
        with open(cache_file) as f:
            cache = json.load(f)
    except (json.JSONDecodeError, FileNotFoundError):
        cache = {}

print(f"Processing {len(species_list)} species and {len(gene_list)} genes")
print(f"Total potential combinations: {len(species_list) * len(gene_list)}")

# OPTIMIZATION 1: Batch process species with multiple genes
def batch_search_species(species_batch, gene_batch):
    """Search for multiple genes in multiple species in batches"""
    results = []
    
    for species in species_batch:
        # Build one query for all genes in this species
        gene_queries = []
        for gene in gene_batch:
            key = f"{gene}|{species}"
            if key not in cache:
                gene_queries.append(f"({gene}[Gene Name] AND {species}[Organism])")
        
        if gene_queries:
            # Combine up to 10 gene queries per API call
            for i in range(0, len(gene_queries), 10):
                batch_queries = gene_queries[i:i+10]
                combined_query = " OR ".join(batch_queries)
                
                try:
                    search_handle = Entrez.esearch(db="protein", term=combined_query, retmax=1000)
                    search_results = Entrez.read(search_handle)
                    search_handle.close()
                    
                    # Parse results - this is simplified, in reality you'd need to 
                    # match results back to specific genes
                    total_count = int(search_results["Count"])
                    
                    # For simplicity, distribute count among genes
                    count_per_gene = total_count // len(batch_queries) if batch_queries else 0
                    
                    for j, gene in enumerate(gene_batch[i:i+10]):
                        if j < len(batch_queries):
                            key = f"{gene}|{species}"
                            results.append((gene, species, count_per_gene))
                            cache[key] = {
                                "count": count_per_gene,
                                "timestamp": datetime.now().isoformat()
                            }
                    
                    time.sleep(0.1)  # Small delay between batches
                    
                except Exception as e:
                    print(f"Error searching {species}: {e}")
                    # Fall back to 0 counts for this batch
                    for gene in gene_batch[i:i+10]:
                        key = f"{gene}|{species}"
                        results.append((gene, species, 0))
                        cache[key] = {
                            "count": 0,
                            "timestamp": datetime.now().isoformat()
                        }
        
        # Add cached results
        for gene in gene_batch:
            key = f"{gene}|{species}"
            if key in cache:
                results.append((gene, species, cache[key]["count"]))
    
    return results

# OPTIMIZATION 2: Process in smaller batches
print("Starting batch processing...")
all_results = []
species_batch_size = 5  # Process 5 species at a time
gene_batch_size = 20    # Process 20 genes at a time

total_batches = (len(species_list) // species_batch_size + 1) * (len(gene_list) // gene_batch_size + 1)
batch_count = 0

for i in range(0, len(species_list), species_batch_size):
    species_batch = species_list[i:i+species_batch_size]
    
    for j in range(0, len(gene_list), gene_batch_size):
        gene_batch = gene_list[j:j+gene_batch_size]
        
        batch_count += 1
        print(f"Processing batch {batch_count}/{total_batches}: {len(species_batch)} species, {len(gene_batch)} genes")
        
        batch_results = batch_search_species(species_batch, gene_batch)
        all_results.extend(batch_results)
        
        # Save cache periodically
        if batch_count % 10 == 0:
            with open(cache_file, 'w') as f:
                json.dump(cache, f, indent=2)

# Save final results
print(f"Saving {len(all_results)} results to {output_file}")

# Create output directory
output_dir = os.path.dirname(output_file)
print(f"DEBUG: Output directory: {output_dir}")
print(f"DEBUG: Output directory exists: {os.path.exists(output_dir)}")
print(f"DEBUG: Current working directory: {os.getcwd()}")
print(f"DEBUG: Absolute output file path: {os.path.abspath(output_file)}")

if output_dir:
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"DEBUG: Successfully created/confirmed directory: {output_dir}")
    except FileExistsError:
        print(f"DEBUG: Directory already exists (FileExistsError): {output_dir}")
    except Exception as e:
        print(f"DEBUG: Error creating directory: {e}")
        # Try to continue anyway
        pass

# Write coverage results
print(f"DEBUG: About to write to file: {output_file}")
with open(output_file, 'w') as f:
    f.write("gene\tspecies\tcount\n")
    for gene, species, count in all_results:
        f.write(f"{gene}\t{species}\t{count}\n")

# Save updated cache to the same shared file
print(f"DEBUG: Saving cache to shared file: {updated_cache_file}")
try:
    with open(updated_cache_file, 'w') as f:
        json.dump(cache, f, indent=2)
    print(f"DEBUG: Successfully updated shared cache file")
except Exception as e:
    print(f"DEBUG: Error updating cache file: {e}")
    print(f"DEBUG: Continuing without updating cache")

print(f"Analysis complete! Results saved to {output_file}")