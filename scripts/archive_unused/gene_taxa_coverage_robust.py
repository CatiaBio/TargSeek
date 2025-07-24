#!/usr/bin/env python3
"""
Robust Gene Taxa Coverage Analysis Script
=========================================

This script assesses how many species have each gene annotated using the NCBI Protein database.
It includes comprehensive error handling, cross-platform compatibility, and reliable caching.
"""

import json
import os
import time
import threading
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

# Thread-safe cache operations
cache_lock = threading.Lock()

class RobustCacheManager:
    """Thread-safe cache manager with proper error handling"""
    
    def __init__(self, cache_dir: str):
        self.cache_dir = Path(cache_dir)
        self.cache_file = self.cache_dir / "gene_species_cache.json"
        self.cache = {}
        self._ensure_cache_dir()
        self._load_cache()
    
    def _ensure_cache_dir(self):
        """Ensure cache directory exists"""
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            print(f"Cache directory ready: {self.cache_dir}")
        except Exception as e:
            print(f"Warning: Could not create cache directory: {e}")
            # Fallback to current directory
            self.cache_dir = Path(".")
            self.cache_file = self.cache_dir / "gene_species_cache.json"
    
    def _load_cache(self):
        """Load cache from file"""
        with cache_lock:
            if self.cache_file.exists():
                try:
                    with self.cache_file.open('r') as f:
                        self.cache = json.load(f)
                    print(f"Loaded cache with {len(self.cache)} entries")
                except Exception as e:
                    print(f"Warning: Could not load cache: {e}")
                    self.cache = {}
            else:
                self.cache = {}
    
    def save_cache(self):
        """Save cache to file"""
        with cache_lock:
            try:
                with self.cache_file.open('w') as f:
                    json.dump(self.cache, f, indent=2)
                print(f"Cache saved with {len(self.cache)} entries")
            except Exception as e:
                print(f"Warning: Could not save cache: {e}")
    
    def get(self, key: str) -> Optional[Dict]:
        """Get cached value"""
        with cache_lock:
            return self.cache.get(key)
    
    def set(self, key: str, value: Dict):
        """Set cached value"""
        with cache_lock:
            self.cache[key] = value

def setup_ncbi_credentials(credentials_file: str) -> Tuple[str, Optional[str]]:
    """Setup NCBI credentials from file"""
    try:
        with open(credentials_file) as f:
            lines = f.readlines()
        email = lines[0].strip()
        api_key = lines[1].strip() if len(lines) > 1 else None
        return email, api_key
    except Exception as e:
        raise ValueError(f"Could not read NCBI credentials from {credentials_file}: {e}")

def load_input_files(species_file: str, genes_file: str) -> Tuple[List[str], List[str]]:
    """Load species and genes from input files"""
    try:
        with open(species_file) as f:
            species_list = [line.strip() for line in f if line.strip()]
        
        with open(genes_file) as f:
            gene_list = [line.strip() for line in f if line.strip()]
        
        return species_list, gene_list
    except Exception as e:
        raise ValueError(f"Could not read input files: {e}")

def search_gene_species_combination(gene: str, species: str, cache_manager: RobustCacheManager) -> int:
    """Search for a gene-species combination in NCBI"""
    cache_key = f"{gene}|{species}"
    
    # Check cache first
    cached_result = cache_manager.get(cache_key)
    if cached_result:
        return cached_result["count"]
    
    # Search NCBI
    try:
        query = f"({gene}[Gene Name] AND {species}[Organism])"
        search_handle = Entrez.esearch(db="protein", term=query, retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        count = int(search_results["Count"])
        
        # Cache the result
        cache_manager.set(cache_key, {
            "count": count,
            "timestamp": datetime.now().isoformat()
        })
        
        return count
        
    except Exception as e:
        print(f"Error searching {gene} in {species}: {e}")
        return 0

def process_batch(species_batch: List[str], gene_batch: List[str], 
                 cache_manager: RobustCacheManager) -> List[Tuple[str, str, int]]:
    """Process a batch of species-gene combinations"""
    results = []
    
    for species in species_batch:
        for gene in gene_batch:
            count = search_gene_species_combination(gene, species, cache_manager)
            results.append((gene, species, count))
            time.sleep(0.1)  # Rate limiting
    
    return results

def ensure_output_directory(output_file: str):
    """Ensure output directory exists"""
    output_path = Path(output_file)
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Output directory ready: {output_path.parent}")
    except Exception as e:
        raise ValueError(f"Could not create output directory: {e}")

def save_results(results: List[Tuple[str, str, int]], output_file: str):
    """Save results to TSV file"""
    ensure_output_directory(output_file)
    
    try:
        with open(output_file, 'w') as f:
            f.write("gene\tspecies\tcount\n")
            for gene, species, count in results:
                f.write(f"{gene}\t{species}\t{count}\n")
        print(f"Results saved to {output_file}")
    except Exception as e:
        raise ValueError(f"Could not save results to {output_file}: {e}")

def main():
    """Main function"""
    print("Starting Gene Taxa Coverage Analysis...")
    
    # Get inputs from Snakemake or fallback to test values
    try:
        species_file = snakemake.input.species
        genes_file = snakemake.input.genes
        credentials_file = snakemake.input.ncbi_info
        output_file = snakemake.output.coverage
        cache_dir = snakemake.params.cache_dir
        print(f"Using Snakemake inputs - output: {output_file}")
    except NameError:
        # Fallback for testing
        print("Snakemake object not available, using test values")
        species_file = "data/bacdive/analysis_1/gram_positive.txt"
        genes_file = "data/quickgo/params_1/gene_symbols_filtered.txt"
        credentials_file = "config/login/ncbi_info.txt"
        output_file = "results/coverage/test_coverage.tsv"
        cache_dir = "results/cache"
    
    # Setup NCBI credentials
    email, api_key = setup_ncbi_credentials(credentials_file)
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    # Load input files
    species_list, gene_list = load_input_files(species_file, genes_file)
    print(f"Processing {len(species_list)} species and {len(gene_list)} genes")
    
    # Initialize cache manager
    cache_manager = RobustCacheManager(cache_dir)
    
    # Process in batches
    all_results = []
    batch_size = 5
    
    for i in range(0, len(species_list), batch_size):
        species_batch = species_list[i:i+batch_size]
        batch_results = process_batch(species_batch, gene_list, cache_manager)
        all_results.extend(batch_results)
        
        # Save cache periodically
        if i % (batch_size * 10) == 0:
            cache_manager.save_cache()
    
    # Save final results and cache
    save_results(all_results, output_file)
    cache_manager.save_cache()
    
    print(f"Analysis complete! Processed {len(all_results)} combinations.")

if __name__ == "__main__":
    main()