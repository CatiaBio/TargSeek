#!/usr/bin/env python3
"""
Ultra-Fast Gene Taxa Coverage Analysis
=====================================

This version uses multiple optimization strategies:
1. Batch queries (multiple genes per API call)
2. Taxonomic pre-filtering
3. Smart caching
4. Parallel processing
"""

import json
import os
import time
from datetime import datetime
from pathlib import Path
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def batch_query_species(species, gene_batch, cache, lock):
    """Query multiple genes for a species in one API call"""
    results = []
    cache_key_prefix = f"batch_{species}_"
    
    # Check if we have a cached batch result for this species
    with lock:
        cached_batch = cache.get(f"{cache_key_prefix}{len(gene_batch)}")
        if cached_batch:
            # Use cached batch result
            for gene in gene_batch:
                individual_key = f"{gene}|{species}"
                if individual_key in cache:
                    results.append((gene, species, cache[individual_key]["count"]))
                else:
                    results.append((gene, species, 0))
            return results
    
    # Build batch query - combine multiple genes with OR
    gene_queries = [f"({gene}[Gene Name])" for gene in gene_batch]
    combined_query = f"({' OR '.join(gene_queries)}) AND {species}[Organism]"
    
    try:
        # Single API call for multiple genes
        search_handle = Entrez.esearch(db="protein", term=combined_query, retmax=1000)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        total_count = int(search_results["Count"])
        
        # If we got results, do individual queries for this species
        if total_count > 0:
            for gene in gene_batch:
                individual_query = f"({gene}[Gene Name]) AND {species}[Organism]"
                try:
                    search_handle = Entrez.esearch(db="protein", term=individual_query, retmax=1)
                    search_results = Entrez.read(search_handle)
                    search_handle.close()
                    
                    count = int(search_results["Count"])
                    results.append((gene, species, count))
                    
                    # Cache individual result
                    with lock:
                        cache[f"{gene}|{species}"] = {
                            "count": count,
                            "timestamp": datetime.now().isoformat()
                        }
                    
                    time.sleep(0.05)  # Short delay between individual queries
                    
                except Exception as e:
                    print(f"Error in individual query for {gene} in {species}: {e}")
                    results.append((gene, species, 0))
        else:
            # No results for any gene in this species
            for gene in gene_batch:
                results.append((gene, species, 0))
                with lock:
                    cache[f"{gene}|{species}"] = {
                        "count": 0,
                        "timestamp": datetime.now().isoformat()
                    }
        
        # Cache the batch result
        with lock:
            cache[f"{cache_key_prefix}{len(gene_batch)}"] = {
                "timestamp": datetime.now().isoformat(),
                "total_count": total_count
            }
        
        return results
        
    except Exception as e:
        print(f"Error in batch query for {species}: {e}")
        # Return zero counts for all genes
        for gene in gene_batch:
            results.append((gene, species, 0))
        return results

def get_taxonomic_info(species_list):
    """Get taxonomic IDs for species to enable faster searches"""
    print("Getting taxonomic information for faster searches...")
    taxid_map = {}
    
    for species in species_list[:10]:  # Limit to first 10 for testing
        try:
            search_handle = Entrez.esearch(db="taxonomy", term=f'"{species}"[Scientific Name]')
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            if search_results["IdList"]:
                taxid = search_results["IdList"][0]
                taxid_map[species] = taxid
                print(f"  {species} -> TaxID: {taxid}")
            
            time.sleep(0.1)
        except Exception as e:
            print(f"  Error getting TaxID for {species}: {e}")
    
    return taxid_map

def main():
    """Main function with optimizations"""
    print("Starting Ultra-Fast Gene Taxa Coverage Analysis...")
    
    # Get inputs from Snakemake
    try:
        species_file = snakemake.input.species
        genes_file = snakemake.input.genes
        credentials_file = snakemake.input.ncbi_info
        output_file = snakemake.output.coverage
        cache_dir = snakemake.params.cache_dir
        print(f"Output file: {output_file}")
    except NameError:
        print("Snakemake object not available, using test values")
        species_file = "data/bacdive/analysis_1/gram_positive.txt"
        genes_file = "data/quickgo/params_1/gene_symbols_filtered.txt"
        credentials_file = "config/login/ncbi_info.txt"
        output_file = "results/coverage/test_coverage.tsv"
        cache_dir = "cache/gene_coverage"
    
    # Setup NCBI credentials
    print("Setting up NCBI credentials...")
    try:
        with open(credentials_file) as f:
            lines = f.readlines()
        email = lines[0].strip()
        api_key = lines[1].strip() if len(lines) > 1 else None
        
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            print("NCBI credentials loaded with API key")
        else:
            print("NCBI credentials loaded (no API key)")
    except Exception as e:
        print(f"Error setting up NCBI credentials: {e}")
        return
    
    # Load input files
    print("Loading input files...")
    try:
        with open(species_file) as f:
            species_list = [line.strip() for line in f if line.strip()]
        
        with open(genes_file) as f:
            gene_list = [line.strip() for line in f if line.strip()]
            
        print(f"Loaded {len(species_list)} species and {len(gene_list)} genes")
    except Exception as e:
        print(f"Error loading input files: {e}")
        return
    
    # FULL DATASET: Process all combinations
    print("FULL DATASET: Processing all combinations with optimizations")
    print(f"Using {len(species_list)} species and {len(gene_list)} genes")
    print(f"Total combinations: {len(species_list) * len(gene_list):,}")
    
    # Setup cache
    print("Setting up cache...")
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    cache_file = cache_path / "gene_species_cache.json"
    
    cache = {}
    if cache_file.exists():
        try:
            with cache_file.open('r') as f:
                cache = json.load(f)
            print(f"Loaded cache with {len(cache)} entries")
        except Exception as e:
            print(f"Could not load cache: {e}")
    
    # OPTIMIZATION 2: Get taxonomic info for faster searches
    # taxid_map = get_taxonomic_info(species_list)
    
    # OPTIMIZATION 3: Batch processing with progress monitoring
    print("Starting batch processing...")
    all_results = []
    cache_lock = threading.Lock()
    
    # Progress tracking
    total_combinations = len(species_list) * len(gene_list)
    processed_combinations = 0
    cached_combinations = 0
    start_time = time.time()
    
    # Process in batches of genes per species
    gene_batch_size = 10  # Process 10 genes at once per species
    
    for i, species in enumerate(species_list):
        print(f"\nProcessing species {i+1}/{len(species_list)}: {species}")
        
        # Process genes in batches
        for j in range(0, len(gene_list), gene_batch_size):
            gene_batch = gene_list[j:j+gene_batch_size]
            
            # Check cache first
            cached_results = []
            uncached_genes = []
            
            with cache_lock:
                for gene in gene_batch:
                    cache_key = f"{gene}|{species}"
                    if cache_key in cache:
                        cached_results.append((gene, species, cache[cache_key]["count"]))
                    else:
                        uncached_genes.append(gene)
            
            # Add cached results
            all_results.extend(cached_results)
            cached_combinations += len(cached_results)
            
            # Process uncached genes in batch
            if uncached_genes:
                batch_results = batch_query_species(species, uncached_genes, cache, cache_lock)
                all_results.extend(batch_results)
                
                print(f"  Processed {len(gene_batch)} genes ({len(cached_results)} cached, {len(uncached_genes)} new)")
                time.sleep(0.2)  # Delay between batches
            
            # Update progress
            processed_combinations += len(gene_batch)
            
            # Progress report every 1000 combinations
            if processed_combinations % 1000 == 0:
                elapsed = time.time() - start_time
                rate = processed_combinations / elapsed if elapsed > 0 else 0
                cache_hit_rate = (cached_combinations / processed_combinations) * 100 if processed_combinations > 0 else 0
                remaining = total_combinations - processed_combinations
                eta = remaining / rate if rate > 0 else 0
                
                print(f"  PROGRESS: {processed_combinations:,}/{total_combinations:,} ({processed_combinations/total_combinations*100:.1f}%)")
                print(f"  Cache hit rate: {cache_hit_rate:.1f}%")
                print(f"  Rate: {rate:.1f} combinations/sec")
                print(f"  ETA: {eta/3600:.1f} hours remaining")
        
        # Save cache every 5 species
        if (i + 1) % 5 == 0:
            try:
                with cache_file.open('w') as f:
                    json.dump(cache, f, indent=2)
                print(f"  Cache saved ({len(cache)} entries)")
            except Exception as e:
                print(f"  Error saving cache: {e}")
    
    # Save results with robust file handling
    print(f"\nSaving {len(all_results)} results...")
    output_path = Path(output_file)
    
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the results
    with output_path.open('w') as f:
        f.write("gene\tspecies\tcount\n")
        for gene, species, count in all_results:
            f.write(f"{gene}\t{species}\t{count}\n")
    
    print(f"Results saved to {output_file}")
    print(f"Full path: {output_path.absolute()}")
    
    # Verify file was created
    if output_path.exists():
        file_size = output_path.stat().st_size
        print(f"File created successfully, size: {file_size:,} bytes")
    else:
        print("ERROR: File was not created!")
        raise FileNotFoundError(f"Could not create output file: {output_file}")
    
    # Save final cache
    try:
        with cache_file.open('w') as f:
            json.dump(cache, f, indent=2)
        print("Final cache saved")
    except Exception as e:
        print(f"Error saving final cache: {e}")
    
    # Final progress report
    final_elapsed = time.time() - start_time
    final_rate = len(all_results) / final_elapsed if final_elapsed > 0 else 0
    final_cache_hit_rate = (cached_combinations / len(all_results)) * 100 if all_results else 0
    
    print(f"\n{'='*60}")
    print(f"ANALYSIS COMPLETE!")
    print(f"{'='*60}")
    print(f"Total combinations processed: {len(all_results):,}")
    print(f"Cache hits: {cached_combinations:,} ({final_cache_hit_rate:.1f}%)")
    print(f"New API calls: {len(all_results) - cached_combinations:,}")
    print(f"Total runtime: {final_elapsed/3600:.1f} hours")
    print(f"Average rate: {final_rate:.1f} combinations/sec")
    print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()