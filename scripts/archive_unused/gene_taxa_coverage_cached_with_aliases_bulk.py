#!/usr/bin/env python3
"""
Enhanced cached gene-taxa coverage check with bulk alias queries
================================================================

This script performs gene-taxa coverage analysis with optimized bulk queries:
1. Primary search using original gene names
2. Bulk alias-based recovery for missing species using all aliases in single queries
3. Final consolidated results

Key optimization: Instead of checking aliases one by one, this version creates
bulk queries with all aliases for a gene, making it much faster.

Cache format:
- Key: "species||gene" (e.g., "Achromobacter xylosoxidans||eno")
- Value: boolean (True if species has gene, False if not)
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
import os
import json
from pathlib import Path

# Cache file location
CACHE_DIR = Path("cache/gene_species")
CACHE_FILE = CACHE_DIR / "species_gene_cache.json"

print(f"=== Starting Optimized Bulk Query Gene-Taxa Coverage Analysis ===")
print(f"Script: {__file__}")
print(f"Output file: {snakemake.output[0]}")
print(f"Cache file: {CACHE_FILE}")

# Read input files
with open(snakemake.input.species_list, 'r') as f:
    species_list = [line.strip() for line in f if line.strip()]

with open(snakemake.input.gene_list, 'r') as f:
    gene_list = [line.strip() for line in f if line.strip()]

# Load gene aliases (filter to only current gene list for efficiency)
gene_aliases = {}
aliases_file = Path(snakemake.input.aliases)
aliases_json_file = aliases_file.with_suffix('.json')

if aliases_json_file.exists():
    with open(aliases_json_file, 'r') as f:
        all_aliases = json.load(f)
    
    # Filter aliases to only include genes in current analysis
    gene_aliases = {gene: all_aliases[gene] for gene in gene_list if gene in all_aliases}
    print(f"Loaded aliases for {len(gene_aliases)} genes (filtered from {len(all_aliases)} total genes)")
else:
    print("Warning: gene_aliases.json not found, proceeding without alias-based recovery")

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

# Load existing cache
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

# Load existing cache
cache = load_cache()

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
                print(f"HTTP Error {e.code}: {e}")
                return 0
        except Exception as e:
            print(f"Search error: {e}")
            return 0

def search_species_gene_cached(species, gene, cache):
    """Search for species-gene pair with caching"""
    cache_key = get_cache_key(species, gene)
    
    # Check cache first
    if cache_key in cache:
        cached_value = cache[cache_key]
        return cached_value if isinstance(cached_value, int) else (1 if cached_value else 0)
    
    # If not in cache, query NCBI
    query = f'("{gene}"[Gene Name] OR "{gene}"[All Fields]) AND "{species}"[Organism]'
    count = rate_limited_search(query)
    
    # Store in cache as boolean
    cache[cache_key] = count > 0
    
    return count

def search_species_gene_bulk_aliases(species, gene, aliases, cache, max_query_length=2000, max_aliases_per_query=20):
    """
    Two-stage search with cache-first approach:
    
    Stage 1: Check cache for primary gene name
    Stage 2: If cache says false, check ALL aliases in bulk queries (with length protection)
    
    Args:
        species: Species name
        gene: Primary gene name  
        aliases: Dictionary of gene aliases
        cache: Query cache
        max_query_length: Maximum query length to avoid URL limits (default: 2000 chars)
        max_aliases_per_query: Maximum aliases per query to avoid complexity issues (default: 20)
    
    Returns:
        tuple: (count, found_via_alias_name or None)
    """
    # STAGE 1: Check cache first for primary gene
    cache_key = get_cache_key(species, gene)
    
    if cache_key in cache:
        cached_value = cache[cache_key]
        cached_count = cached_value if isinstance(cached_value, int) else (1 if cached_value else 0)
        
        # If cache says gene is found (True), return immediately
        if cached_count > 0:
            return cached_count, None
        
        # If cache says gene is NOT found (False), try aliases
        # (Don't query NCBI again for the primary gene - we trust the cache)
    else:
        # Not in cache - query NCBI for primary gene
        query = f'("{gene}"[Gene Name] OR "{gene}"[All Fields]) AND "{species}"[Organism]'
        primary_count = rate_limited_search(query)
        cache[cache_key] = primary_count > 0
        
        if primary_count > 0:
            return primary_count, None  # Found with primary name
    
    # STAGE 2: Primary gene not found (either cached as false or queried as 0) - check aliases
    # Only proceed if we actually have aliases for this gene
    if gene not in aliases:
        return 0, None  # No aliases available
    
    aliases_list = aliases[gene]
    # Filter out the primary gene name and create unique list
    unique_aliases = list(set([alias for alias in aliases_list if alias != gene]))
    
    if not unique_aliases:
        return 0, None  # No unique aliases to check
    
    # Create one comprehensive bulk query with ALL aliases
    # First try to fit everything in one query
    all_alias_terms = [f'("{alias}"[Gene Name] OR "{alias}"[All Fields])' for alias in unique_aliases]
    full_bulk_query = f'({" OR ".join(all_alias_terms)}) AND "{species}"[Organism]'
    
    # Check if the full query is within limits
    if len(full_bulk_query) <= max_query_length and len(unique_aliases) <= max_aliases_per_query:
        # Single bulk query with all aliases
        cache_key = get_cache_key(species, f"{gene}_all_aliases_bulk")
        
        if cache_key in cache:
            cached_result = cache[cache_key]
            bulk_count = cached_result if isinstance(cached_result, int) else (1 if cached_result else 0)
        else:
            bulk_count = rate_limited_search(full_bulk_query)
            cache[cache_key] = bulk_count > 0
        
        if bulk_count > 0:
            return bulk_count, unique_aliases[0]  # Return first alias as "found via"
        else:
            return 0, None
    
    # Query too long - need to split into chunks
    else:
        alias_chunks = []
        current_chunk = []
        current_query_length = len(f' AND "{species}"[Organism]')
        
        for alias in unique_aliases:
            alias_term = f'("{alias}"[Gene Name] OR "{alias}"[All Fields])'
            chunk_query_length = current_query_length + len(alias_term) + 4  # +4 for " OR "
            
            # Start new chunk if we exceed limits
            if (len(current_chunk) >= max_aliases_per_query or 
                chunk_query_length > max_query_length):
                if current_chunk:  # Don't add empty chunks
                    alias_chunks.append(current_chunk)
                current_chunk = [alias]
                current_query_length = len(f' AND "{species}"[Organism]') + len(alias_term)
            else:
                current_chunk.append(alias)
                current_query_length = chunk_query_length
        
        # Add final chunk
        if current_chunk:
            alias_chunks.append(current_chunk)
        
        # Query each chunk - return on first hit
        for chunk_idx, alias_chunk in enumerate(alias_chunks):
            alias_terms = [f'("{alias}"[Gene Name] OR "{alias}"[All Fields])' for alias in alias_chunk]
            chunk_query = f'({" OR ".join(alias_terms)}) AND "{species}"[Organism]'
            
            chunk_cache_key = get_cache_key(species, f"{gene}_aliases_chunk_{chunk_idx}")
            
            if chunk_cache_key in cache:
                cached_result = cache[chunk_cache_key]
                chunk_count = cached_result if isinstance(cached_result, int) else (1 if cached_result else 0)
            else:
                chunk_count = rate_limited_search(chunk_query)
                cache[chunk_cache_key] = chunk_count > 0
            
            if chunk_count > 0:
                return chunk_count, alias_chunk[0]  # Found in this chunk
        
        # All chunks checked, nothing found
        return 0, None

def batch_search_species_with_bulk_aliases(species, genes, aliases, cache):
    """Enhanced search for species with bulk alias queries"""
    results = {}
    alias_recoveries = {}
    
    # Check each gene
    for gene in genes:
        count, found_via_alias = search_species_gene_bulk_aliases(species, gene, aliases, cache)
        results[gene] = count
        if found_via_alias:
            alias_recoveries[gene] = found_via_alias
    
    return species, results, alias_recoveries

# Calculate cache statistics
total_queries_needed = len(species_list) * len(gene_list)
cached_queries = 0
for species in species_list:
    for gene in gene_list:
        cache_key = get_cache_key(species, gene)
        if cache_key in cache:
            cached_queries += 1

queries_to_make = total_queries_needed - cached_queries
cache_hit_rate = (cached_queries / total_queries_needed) * 100

print(f"\nCache Statistics:")
print(f"Total queries needed: {total_queries_needed}")
print(f"Already cached: {cached_queries}")
print(f"New queries to make: {queries_to_make}")
print(f"Cache hit rate: {cache_hit_rate:.1f}%")

if queries_to_make == 0:
    print("🎉 All queries cached! No NCBI requests needed.")
else:
    estimated_time = queries_to_make * min_delay / 60
    print(f"Estimated time for new queries: {estimated_time:.1f} minutes")

# Use parallel processing for different species
print(f"\nStarting bulk alias analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Determine number of workers based on cache hit rate
if cache_hit_rate > 90:
    max_workers = 10  # High cache hit rate, can go faster
elif cache_hit_rate > 50:
    max_workers = 5 if api_key else 2  # Medium cache hit rate
else:
    max_workers = 3 if api_key else 1  # Low cache hit rate, be conservative

results_dict = {gene: {} for gene in gene_list}
all_alias_recoveries = {}
start_time = time.time()

# Single stage: Primary gene search with bulk alias recovery
print("\n=== Integrated Primary + Bulk Alias Search ===")
try:
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all species for processing with bulk alias support
        future_to_species = {
            executor.submit(batch_search_species_with_bulk_aliases, species, gene_list, gene_aliases, cache): species 
            for species in species_list
        }
        
        completed = 0
        total_recoveries = 0
        
        for future in as_completed(future_to_species):
            try:
                species, gene_counts, alias_recoveries = future.result()
                completed += 1
                
                # Store results
                for gene, count in gene_counts.items():
                    results_dict[gene][species] = count
                
                # Track alias recoveries
                if alias_recoveries:
                    for gene, alias_used in alias_recoveries.items():
                        if gene not in all_alias_recoveries:
                            all_alias_recoveries[gene] = {}
                        all_alias_recoveries[gene][species] = alias_used
                        total_recoveries += 1
                
                # Progress update
                elapsed = time.time() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                eta = (len(species_list) - completed) / rate if rate > 0 else 0
                
                recovery_info = f", {len(alias_recoveries)} recoveries" if alias_recoveries else ""
                print(f"Progress: {completed}/{len(species_list)} species completed ({completed/len(species_list)*100:.1f}%){recovery_info} - ETA: {eta/60:.1f} min")
                
                # Save cache periodically (every 10 species)
                if completed % 10 == 0:
                    save_cache(cache)
                
            except Exception as e:
                print(f"Error processing species: {e}")
                completed += 1

except Exception as e:
    print(f"Error in thread executor: {e}")

# Save final cache
save_cache(cache)

# Create DataFrame
print("\n=== Final Results Processing ===")
data_for_df = []
for gene in gene_list:
    row = {'gene': gene}
    row.update(results_dict[gene])
    data_for_df.append(row)

df = pd.DataFrame(data_for_df)

# Calculate coverage statistics
print("Calculating bulk alias coverage statistics...")
coverage_stats = []
for idx, row in df.iterrows():
    gene = row['gene']
    # Count non-zero entries and collect species names
    species_with_gene_list = []
    species_with_gene_count = 0
    
    for species in species_list:
        if row.get(species, 0) > 0:
            species_with_gene_count += 1
            species_with_gene_list.append(species)
    
    coverage_percentage = (species_with_gene_count / len(species_list)) * 100
    species_names_with_gene = ','.join(species_with_gene_list)
    
    # Add recovery information
    recovered_via_aliases = len(all_alias_recoveries.get(gene, {}))
    
    coverage_stats.append({
        'gene': gene,
        'species_with_gene': species_with_gene_count,
        'total_species': len(species_list),
        'coverage_percentage': coverage_percentage,
        'species_names_with_gene': species_names_with_gene,
        'recovered_via_aliases': recovered_via_aliases
    })

# Create coverage summary DataFrame
coverage_df = pd.DataFrame(coverage_stats)

# Sort by coverage
coverage_df = coverage_df.sort_values('coverage_percentage', ascending=False)

# Ensure output directory exists
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

# Save results
print(f"Saving bulk alias results to: {snakemake.output[0]}")
try:
    coverage_df.to_csv(snakemake.output[0], sep='\t', index=False)
    print(f"✓ Successfully saved results to: {snakemake.output[0]}")
    
    # Save alias recovery details as supplementary file
    if all_alias_recoveries:
        recovery_file = output_path.with_suffix('.recovery_details.json')
        with open(recovery_file, 'w') as f:
            json.dump(all_alias_recoveries, f, indent=2)
        print(f"✓ Saved bulk alias recovery details to: {recovery_file}")
    
    # Verify file was created and has content
    if output_path.exists():
        file_size = output_path.stat().st_size
        print(f"✓ File exists with size: {file_size} bytes")
    else:
        print(f"✗ ERROR: Output file was not created!")
        
except Exception as e:
    print(f"✗ ERROR saving results: {e}")

# Print enhanced summary
print(f"\n{'='*60}")
print(f"Bulk alias analysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total time: {(time.time() - start_time)/60:.1f} minutes")
print(f"Genes processed: {len(gene_list)}")
print(f"Species processed: {len(species_list)}")
print(f"Cache entries: {len(cache)}")
print(f"Cache hit rate: {cache_hit_rate:.1f}%")

if all_alias_recoveries:
    total_recoveries = sum(len(species_dict) for species_dict in all_alias_recoveries.values())
    print(f"Bulk alias recoveries: {total_recoveries} species recovered for {len(all_alias_recoveries)} genes")

print(f"Output file: {snakemake.output[0]}")

if len(coverage_df) > 0:
    print(f"\nTop 10 genes by coverage (with bulk alias recoveries):")
    top_genes = coverage_df.head(10).copy()
    
    # Create display dataframe with recovery info
    display_cols = ['gene', 'coverage_percentage', 'species_with_gene', 'total_species', 'recovered_via_aliases']
    print(top_genes[display_cols].to_string(index=False))
    
    if 'recovered_via_aliases' in coverage_df.columns:
        total_recovered = coverage_df['recovered_via_aliases'].sum()
        if total_recovered > 0:
            print(f"\n🚀 Bulk alias recovery improved coverage for {len(coverage_df[coverage_df['recovered_via_aliases'] > 0])} genes")
            print(f"   Total species recovered: {total_recovered}")
else:
    print(f"\nWARNING: No results generated!")
    
print(f"=== Bulk Alias Analysis Complete ===")