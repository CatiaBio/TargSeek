#!/usr/bin/env python3
"""
Gene Taxa Coverage Analysis - No Cache Version
=============================================

This version performs direct API queries without any caching.
Simplified for testing with small datasets.
"""

import os
import time
from datetime import datetime
from pathlib import Path
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed

def query_gene_species(gene, species):
    """Query a single gene for a species"""
    query = f"{gene}[Gene Name] AND {species}[Organism]"
    
    try:
        # Direct API call without caching
        search_handle = Entrez.esearch(db="protein", term=query, retmax=0)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        count = int(search_results["Count"])
        return (gene, species, count)
        
    except Exception as e:
        print(f"Error querying {gene} for {species}: {e}")
        return (gene, species, 0)

def process_species_batch(species_list, genes, max_workers=3):
    """Process all gene-species combinations"""
    results = []
    total_queries = len(species_list) * len(genes)
    completed = 0
    
    print(f"Processing {len(genes)} genes across {len(species_list)} species...")
    print(f"Total queries to perform: {total_queries}")
    
    # Create all gene-species pairs
    tasks = [(gene, species) for species in species_list for gene in genes]
    
    # Process with thread pool
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(query_gene_species, gene, species): (gene, species)
            for gene, species in tasks
        }
        
        # Process completed tasks
        for future in as_completed(future_to_task):
            result = future.result()
            results.append(result)
            completed += 1
            
            # Progress update
            if completed % 10 == 0 or completed == total_queries:
                print(f"Progress: {completed}/{total_queries} queries completed")
    
    return results

def main():
    """Main execution function"""
    # Set up Entrez
    with open(snakemake.input.ncbi_info, 'r') as f:
        lines = f.readlines()
        Entrez.email = lines[0].strip()
        Entrez.api_key = lines[1].strip() if len(lines) > 1 else None
    
    # Read input files
    print("Reading input files...")
    
    # Read species list
    with open(snakemake.input.species_list, 'r') as f:
        species_list = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(species_list)} species")
    
    # Read gene list
    with open(snakemake.input.gene_list, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(genes)} genes")
    
    # Get output file
    output_file = snakemake.output[0]
    
    # Start processing
    start_time = time.time()
    print(f"\nStarting analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Process all combinations
    all_results = process_species_batch(species_list, genes, max_workers=3)
    
    # Sort results
    all_results.sort(key=lambda x: (x[0], x[1]))
    
    # Save results
    print(f"\nSaving {len(all_results)} results...")
    output_path = Path(output_file)
    
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the results
    with output_path.open('w') as f:
        f.write("gene\tspecies\tcount\n")
        for gene, species, count in all_results:
            f.write(f"{gene}\t{species}\t{count}\n")
    
    elapsed_time = time.time() - start_time
    print(f"\nAnalysis completed in {elapsed_time:.2f} seconds")
    print(f"Results saved to {output_file}")
    
    # Summary statistics
    genes_with_hits = len(set(gene for gene, _, count in all_results if count > 0))
    total_hits = sum(count for _, _, count in all_results)
    print(f"\nSummary:")
    print(f"- Genes with at least one hit: {genes_with_hits}/{len(genes)}")
    print(f"- Total protein records found: {total_hits}")

if __name__ == "__main__":
    main()