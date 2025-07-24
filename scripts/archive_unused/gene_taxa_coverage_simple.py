#!/usr/bin/env python3
"""
Simple Gene Taxa Coverage Analysis
=================================

A simplified version that prioritizes reliability over advanced features.
Provides basic progress monitoring and robust file handling.
"""

import json
import os
import time
from datetime import datetime
from pathlib import Path
from Bio import Entrez

def main():
    """Main function - simplified and reliable"""
    print("Starting Simple Gene Taxa Coverage Analysis...")
    
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
        cache_dir = "results/cache"
    
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
        print(f"Total combinations: {len(species_list) * len(gene_list):,}")
    except Exception as e:
        print(f"Error loading input files: {e}")
        return
    
    # Setup cache
    print("Setting up cache...")
    cache_path = Path(cache_dir)
    try:
        cache_path.mkdir(parents=True, exist_ok=True)
        cache_file = cache_path / "gene_species_cache.json"
        print(f"Cache directory: {cache_path}")
    except Exception as e:
        print(f"Error creating cache directory: {e}")
        # Fallback to current directory
        cache_file = Path("gene_species_cache.json")
        print(f"Using fallback cache: {cache_file}")
    
    # Load existing cache
    cache = {}
    if cache_file.exists():
        try:
            with cache_file.open('r') as f:
                cache = json.load(f)
            print(f"Loaded cache with {len(cache)} entries")
        except Exception as e:
            print(f"Could not load cache: {e}")
            cache = {}
    else:
        print("Starting with empty cache")
    
    # Process combinations
    print("Starting processing...")
    results = []
    processed = 0
    cache_hits = 0
    api_calls = 0
    start_time = time.time()
    
    for i, species in enumerate(species_list):
        print(f"\nProcessing species {i+1}/{len(species_list)}: {species}")
        
        for j, gene in enumerate(gene_list):
            cache_key = f"{gene}|{species}"
            
            # Check cache first
            if cache_key in cache:
                count = cache[cache_key]["count"]
                results.append((gene, species, count))
                cache_hits += 1
                processed += 1
                continue
            
            # Query NCBI
            try:
                query = f"({gene}[Gene Name] AND {species}[Organism])"
                search_handle = Entrez.esearch(db="protein", term=query, retmax=1)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                count = int(search_results["Count"])
                
                # Cache the result
                cache[cache_key] = {
                    "count": count,
                    "timestamp": datetime.now().isoformat()
                }
                
                results.append((gene, species, count))
                api_calls += 1
                processed += 1
                
                # Progress update every 50 items
                if processed % 50 == 0:
                    elapsed = time.time() - start_time
                    rate = processed / elapsed if elapsed > 0 else 0
                    remaining = len(species_list) * len(gene_list) - processed
                    eta = remaining / rate if rate > 0 else 0
                    print(f"  Progress: {processed:,} processed, {cache_hits} cached, {api_calls} API calls, Rate: {rate:.1f}/sec, ETA: {eta/60:.1f}min")
                
                # Small delay
                time.sleep(0.1)
                
            except Exception as e:
                print(f"  Error querying {gene} in {species}: {e}")
                results.append((gene, species, 0))
                processed += 1
                time.sleep(1)  # Longer delay after error
        
        # Save cache every 10 species
        if (i + 1) % 10 == 0:
            try:
                with cache_file.open('w') as f:
                    json.dump(cache, f, indent=2)
                print(f"  Cache saved ({len(cache)} entries)")
            except Exception as e:
                print(f"  Error saving cache: {e}")
    
    # Create output directory and save results
    print(f"\nSaving results to {output_file}")
    output_path = Path(output_file)
    
    try:
        # Create directory
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write results
        with output_path.open('w') as f:
            f.write("gene\tspecies\tcount\n")
            for gene, species, count in results:
                f.write(f"{gene}\t{species}\t{count}\n")
        
        print(f"Results saved successfully!")
        
    except Exception as e:
        print(f"Error saving results: {e}")
        # Fallback to current directory
        fallback_file = Path(output_path.name)
        try:
            with fallback_file.open('w') as f:
                f.write("gene\tspecies\tcount\n")
                for gene, species, count in results:
                    f.write(f"{gene}\t{species}\t{count}\n")
            print(f"Results saved to fallback location: {fallback_file}")
        except Exception as e2:
            print(f"All save attempts failed: {e2}")
    
    # Save final cache
    try:
        with cache_file.open('w') as f:
            json.dump(cache, f, indent=2)
        print("Final cache saved")
    except Exception as e:
        print(f"Error saving final cache: {e}")
    
    # Final report
    elapsed = time.time() - start_time
    print(f"\n=== FINAL REPORT ===")
    print(f"Processed: {processed:,} combinations")
    print(f"Cache hits: {cache_hits:,}")
    print(f"API calls: {api_calls:,}")
    print(f"Runtime: {elapsed/60:.1f} minutes")
    print(f"Rate: {processed/elapsed:.1f} items/second")
    print("Analysis complete!")

if __name__ == "__main__":
    main()