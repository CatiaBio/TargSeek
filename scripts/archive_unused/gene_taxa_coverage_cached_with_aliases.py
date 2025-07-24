#!/usr/bin/env python3
"""
Enhanced cached gene-taxa coverage check with alias-based missing species recovery

This script performs gene-taxa coverage analysis with the following stages:
1. Primary search using original gene names
2. Alias-based recovery for missing species using gene aliases
3. Final consolidated results

Cache format:
- Key: "species||gene" (e.g., "Achromobacter xylosoxidans||eno")
- Value: boolean (True if species has gene, False if not)
- Backward compatible with integer counts from previous cache versions
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

print(f"=== Starting Enhanced Cached Gene-Taxa Coverage Analysis ===")
print(f"Script: {__file__}")
print(f"Output file: {snakemake.output[0]}")
print(f"Cache file: {CACHE_FILE}")

# Read input files
with open(snakemake.input.species_list, 'r') as f:
    species_list = [line.strip() for line in f if line.strip()]

with open(snakemake.input.gene_list, 'r') as f:
    gene_list = [line.strip() for line in f if line.strip()]

# Load gene aliases
gene_aliases = {}
aliases_file = Path(snakemake.input.aliases)
aliases_json_file = aliases_file.with_suffix('.json')

if aliases_json_file.exists():
    with open(aliases_json_file, 'r') as f:
        gene_aliases = json.load(f)
    print(f"Loaded aliases for {len(gene_aliases)} genes")
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

def parse_cache_key(key):
    """Parse cache key back to species, gene"""
    return key.split("||", 1)

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
        # Return the cached count (for backward compatibility with existing cache)
        cached_value = cache[cache_key]
        return cached_value if isinstance(cached_value, int) else (1 if cached_value else 0)
    
    # If not in cache, query NCBI
    query = f'("{gene}"[Gene Name] OR "{gene}"[All Fields]) AND "{species}"[Organism]'
    count = rate_limited_search(query)
    
    # Store in cache as boolean (True if gene found, False if not)
    cache[cache_key] = count > 0
    
    return count

def search_species_gene_with_aliases(species, gene, aliases, cache):
    """
    Enhanced search that tries aliases if primary gene search fails
    
    Returns:
        tuple: (count, found_via_alias_name or None)
    """
    # Try primary gene first
    primary_count = search_species_gene_cached(species, gene, cache)
    if primary_count > 0:
        return primary_count, None
    
    # If primary failed and we have aliases, try them
    if gene in aliases:
        gene_aliases_list = aliases[gene]
        for alias in gene_aliases_list:
            if alias != gene:  # Don't retry the same name
                alias_count = search_species_gene_cached(species, alias, cache)
                if alias_count > 0:
                    return alias_count, alias
    
    # No results found
    return 0, None

def batch_search_species_cached(species, genes, cache):
    """Search all genes for a single species with caching"""
    results = {}
    
    # Check each gene
    for gene in genes:
        count = search_species_gene_cached(species, gene, cache)
        results[gene] = count
    
    return species, results

def batch_search_species_with_aliases(species, genes, aliases, cache):
    """Enhanced search for species with alias-based recovery"""
    results = {}
    alias_recoveries = {}
    
    # Check each gene
    for gene in genes:
        count, found_via_alias = search_species_gene_with_aliases(species, gene, aliases, cache)
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
print(f"\nStarting enhanced analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

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

# Stage 1: Primary gene search
print("\n=== Stage 1: Primary Gene Search ===")
try:
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all species for processing
        future_to_species = {
            executor.submit(batch_search_species_cached, species, gene_list, cache): species 
            for species in species_list
        }
        
        completed = 0
        for future in as_completed(future_to_species):
            try:
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
                
                # Save cache periodically (every 10 species)
                if completed % 10 == 0:
                    save_cache(cache)
                
            except Exception as e:
                print(f"Error processing species: {e}")
                completed += 1

except Exception as e:
    print(f"Error in thread executor: {e}")

# Stage 2: Alias-based recovery for missing species
if gene_aliases:
    print(f"\n=== Stage 2: Alias-based Recovery ===")
    print(f"Looking for missing species using gene aliases...")
    
    # Identify genes with missing species (coverage < 100%)
    recovery_candidates = {}
    for gene in gene_list:
        species_with_gene = sum(1 for species in species_list if results_dict[gene].get(species, 0) > 0)
        coverage_pct = (species_with_gene / len(species_list)) * 100
        
        if coverage_pct < 100 and gene in gene_aliases and gene_aliases[gene]:
            missing_species = [s for s in species_list if results_dict[gene].get(s, 0) == 0]
            recovery_candidates[gene] = missing_species
    
    print(f"Found {len(recovery_candidates)} genes with missing species that have aliases")
    
    if recovery_candidates:
        recoveries_found = 0
        # Calculate total attempts (species × max 2 aliases per gene)
        total_recovery_attempts = sum(len(missing) * min(2, len(gene_aliases.get(gene, []))) 
                                     for gene, missing in recovery_candidates.items())
        recovery_attempt = 0
        
        for gene, missing_species in recovery_candidates.items():
            aliases_list = gene_aliases[gene]
            aliases_to_show = aliases_list[:2]  # Show only first 2 aliases in log
            print(f"Attempting recovery for {gene} ({len(missing_species)} missing species) using first 2 aliases: {', '.join(aliases_to_show)}")
            
            for species in missing_species:
                recovery_attempt += 1
                # Try up to first 2 aliases for this species
                aliases_to_try = [alias for alias in aliases_list[:2] if alias != gene]
                
                for alias in aliases_to_try:
                    count = search_species_gene_cached(species, alias, cache)
                    if count > 0:
                        # Found via alias! Update results
                        results_dict[gene][species] = count
                        recoveries_found += 1
                        
                        # Track the recovery
                        if gene not in all_alias_recoveries:
                            all_alias_recoveries[gene] = {}
                        all_alias_recoveries[gene][species] = alias
                        
                        print(f"  ✓ Found {species} for {gene} via alias '{alias}'")
                        break  # Stop trying other aliases for this species
                
                # Progress update for recovery
                if recovery_attempt % 50 == 0:
                    print(f"Recovery progress: {recovery_attempt}/{total_recovery_attempts} attempts ({recovery_attempt/total_recovery_attempts*100:.1f}%)")
        
        print(f"\nRecovery Summary:")
        print(f"  Species recovered: {recoveries_found}")
        print(f"  Genes with recoveries: {len(all_alias_recoveries)}")
        
        if all_alias_recoveries:
            print(f"  Recovery details:")
            for gene, species_aliases in all_alias_recoveries.items():
                print(f"    {gene}: {len(species_aliases)} species recovered")
                for species, alias in list(species_aliases.items())[:3]:  # Show first 3
                    print(f"      - {species} via '{alias}'")
                if len(species_aliases) > 3:
                    print(f"      ... and {len(species_aliases) - 3} more")

# Save final cache
save_cache(cache)

# Create DataFrame
print("\n=== Stage 3: Final Results Processing ===")
data_for_df = []
for gene in gene_list:
    row = {'gene': gene}
    row.update(results_dict[gene])
    data_for_df.append(row)

df = pd.DataFrame(data_for_df)

# Calculate coverage statistics
print("Calculating enhanced coverage statistics...")
coverage_stats = []
for idx, row in df.iterrows():
    gene = row['gene']
    # Count non-zero entries and collect species names (excluding the 'gene' column)
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

# Save results with error handling
print(f"Saving enhanced results to: {snakemake.output[0]}")
try:
    coverage_df.to_csv(snakemake.output[0], sep='\t', index=False)
    print(f"✓ Successfully saved results to: {snakemake.output[0]}")
    
    # Save alias recovery details as supplementary file
    if all_alias_recoveries:
        recovery_file = output_path.with_suffix('.recovery_details.json')
        with open(recovery_file, 'w') as f:
            json.dump(all_alias_recoveries, f, indent=2)
        print(f"✓ Saved alias recovery details to: {recovery_file}")
    
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
print(f"Enhanced analysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total time: {(time.time() - start_time)/60:.1f} minutes")
print(f"Genes processed: {len(gene_list)}")
print(f"Species processed: {len(species_list)}")
print(f"Cache entries: {len(cache)}")
print(f"Cache hit rate: {cache_hit_rate:.1f}%")

if all_alias_recoveries:
    total_recoveries = sum(len(species_dict) for species_dict in all_alias_recoveries.values())
    print(f"Alias-based recoveries: {total_recoveries} species recovered for {len(all_alias_recoveries)} genes")

print(f"Output file: {snakemake.output[0]}")

if len(coverage_df) > 0:
    print(f"\nTop 10 genes by coverage (including alias recoveries):")
    top_genes = coverage_df.head(10).copy()
    
    # Create display dataframe with recovery info
    display_cols = ['gene', 'coverage_percentage', 'species_with_gene', 'total_species', 'recovered_via_aliases']
    print(top_genes[display_cols].to_string(index=False))
    
    if 'recovered_via_aliases' in coverage_df.columns:
        total_recovered = coverage_df['recovered_via_aliases'].sum()
        if total_recovered > 0:
            print(f"\n🎯 Alias recovery improved coverage for {len(coverage_df[coverage_df['recovered_via_aliases'] > 0])} genes")
            print(f"   Total species recovered: {total_recovered}")
else:
    print(f"\nWARNING: No results generated!")
    
print(f"=== Enhanced Analysis Complete ===")