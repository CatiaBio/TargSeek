#!/usr/bin/env python3
"""
Enhanced Gene Taxa Coverage Analysis with Progress Monitoring
===========================================================

This script provides comprehensive progress tracking, intermediate saves,
recovery capabilities, and detailed monitoring for long-running NCBI queries.
"""

import json
import os
import time
import threading
import signal
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import Entrez
from collections import defaultdict

class ProgressMonitor:
    """Comprehensive progress monitoring and reporting"""
    
    def __init__(self, total_combinations: int, species_count: int, gene_count: int):
        self.total = total_combinations
        self.species_count = species_count
        self.gene_count = gene_count
        self.processed = 0
        self.cached_hits = 0
        self.api_calls = 0
        self.errors = 0
        self.start_time = time.time()
        self.last_report = time.time()
        self.lock = threading.Lock()
        
        # Rate limiting
        self.api_call_times = []
        self.max_calls_per_minute = 300  # Conservative rate limit
        
        print(f"=== PROGRESS MONITOR INITIALIZED ===")
        print(f"Total combinations: {self.total:,}")
        print(f"Species: {self.species_count}, Genes: {self.gene_count}")
        print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"==========================================")
    
    def update(self, cached: bool = False, error: bool = False):
        """Update progress counters"""
        with self.lock:
            self.processed += 1
            if cached:
                self.cached_hits += 1
            else:
                self.api_calls += 1
                self.api_call_times.append(time.time())
                # Keep only last minute of API calls
                cutoff = time.time() - 60
                self.api_call_times = [t for t in self.api_call_times if t > cutoff]
            
            if error:
                self.errors += 1
            
            # Report progress every 100 items or every 60 seconds
            now = time.time()
            if (self.processed % 100 == 0) or (now - self.last_report > 60):
                self.report_progress()
                self.last_report = now
    
    def report_progress(self):
        """Generate detailed progress report"""
        with self.lock:
            elapsed = time.time() - self.start_time
            if elapsed == 0:
                return
            
            # Calculate rates and estimates
            overall_rate = self.processed / elapsed
            api_rate = len(self.api_call_times)  # API calls in last minute
            
            remaining = self.total - self.processed
            eta_seconds = remaining / overall_rate if overall_rate > 0 else 0
            eta_time = datetime.now() + timedelta(seconds=eta_seconds)
            
            # Progress percentage
            progress_pct = (self.processed / self.total) * 100
            
            # Cache hit rate
            cache_hit_rate = (self.cached_hits / self.processed) * 100 if self.processed > 0 else 0
            
            print(f"\n{'='*60}")
            print(f"PROGRESS REPORT - {datetime.now().strftime('%H:%M:%S')}")
            print(f"{'='*60}")
            print(f"Processed: {self.processed:,} / {self.total:,} ({progress_pct:.1f}%)")
            print(f"Cache hits: {self.cached_hits:,} ({cache_hit_rate:.1f}%)")
            print(f"API calls: {self.api_calls:,}")
            print(f"Errors: {self.errors:,}")
            print(f"Runtime: {elapsed/3600:.1f} hours")
            print(f"Rate: {overall_rate:.1f} items/sec")
            print(f"API rate: {api_rate} calls/minute")
            print(f"ETA: {eta_time.strftime('%H:%M:%S')} ({eta_seconds/3600:.1f} hours remaining)")
            print(f"{'='*60}\n")
    
    def should_throttle(self) -> bool:
        """Check if we should throttle API calls"""
        return len(self.api_call_times) >= self.max_calls_per_minute
    
    def final_report(self):
        """Generate final completion report"""
        elapsed = time.time() - self.start_time
        print(f"\n{'='*60}")
        print(f"FINAL REPORT - COMPLETED")
        print(f"{'='*60}")
        print(f"Total processed: {self.processed:,}")
        print(f"Cache hits: {self.cached_hits:,} ({self.cached_hits/self.processed*100:.1f}%)")
        print(f"API calls: {self.api_calls:,}")
        print(f"Errors: {self.errors:,}")
        print(f"Total runtime: {elapsed/3600:.1f} hours")
        print(f"Average rate: {self.processed/elapsed:.1f} items/sec")
        print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*60}\n")

class RobustCacheManager:
    """Thread-safe cache manager with recovery capabilities"""
    
    def __init__(self, cache_dir: str):
        self.cache_dir = Path(cache_dir)
        self.cache_file = self.cache_dir / "gene_species_cache.json"
        self.backup_file = self.cache_dir / "gene_species_cache_backup.json"
        self.cache = {}
        self.lock = threading.Lock()
        self.save_counter = 0
        self._ensure_cache_dir()
        self._load_cache()
    
    def _ensure_cache_dir(self):
        """Ensure cache directory exists"""
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            print(f"Cache directory ready: {self.cache_dir}")
        except Exception as e:
            print(f"Warning: Could not create cache directory: {e}")
            self.cache_dir = Path(".")
            self.cache_file = self.cache_dir / "gene_species_cache.json"
            self.backup_file = self.cache_dir / "gene_species_cache_backup.json"
    
    def _load_cache(self):
        """Load cache from file with backup recovery"""
        with self.lock:
            # Try main cache file first
            if self.cache_file.exists():
                try:
                    with self.cache_file.open('r') as f:
                        self.cache = json.load(f)
                    print(f"Loaded cache with {len(self.cache)} entries")
                    return
                except Exception as e:
                    print(f"Warning: Could not load main cache: {e}")
            
            # Try backup file
            if self.backup_file.exists():
                try:
                    with self.backup_file.open('r') as f:
                        self.cache = json.load(f)
                    print(f"Loaded cache from backup with {len(self.cache)} entries")
                    return
                except Exception as e:
                    print(f"Warning: Could not load backup cache: {e}")
            
            # Start with empty cache
            self.cache = {}
            print("Starting with empty cache")
    
    def save_cache(self, force: bool = False):
        """Save cache to file with backup"""
        with self.lock:
            self.save_counter += 1
            
            # Save every 50 updates or when forced
            if not force and self.save_counter % 50 != 0:
                return
            
            try:
                # Create backup of current cache
                if self.cache_file.exists():
                    self.cache_file.rename(self.backup_file)
                
                # Save new cache
                with self.cache_file.open('w') as f:
                    json.dump(self.cache, f, indent=2)
                
                print(f"Cache saved with {len(self.cache)} entries")
                
            except Exception as e:
                print(f"Warning: Could not save cache: {e}")
    
    def get(self, key: str) -> Optional[Dict]:
        """Get cached value"""
        with self.lock:
            return self.cache.get(key)
    
    def set(self, key: str, value: Dict):
        """Set cached value"""
        with self.lock:
            self.cache[key] = value

class IntermediateResultsSaver:
    """Saves intermediate results to prevent data loss"""
    
    def __init__(self, output_file: str):
        self.output_file = Path(output_file)
        self.temp_file = self.output_file.with_suffix('.tmp')
        self.results = []
        self.lock = threading.Lock()
        
        # Ensure output directory exists with robust creation
        try:
            self.output_file.parent.mkdir(parents=True, exist_ok=True)
            print(f"Output directory created: {self.output_file.parent}")
        except Exception as e:
            print(f"Error creating output directory: {e}")
            # Try alternative approach
            import os
            try:
                os.makedirs(str(self.output_file.parent), exist_ok=True)
                print(f"Output directory created with os.makedirs: {self.output_file.parent}")
            except Exception as e2:
                print(f"Both directory creation methods failed: {e2}")
                raise ValueError(f"Cannot create output directory: {self.output_file.parent}")
        
        # Initialize temp file with header
        try:
            with self.temp_file.open('w') as f:
                f.write("gene\tspecies\tcount\n")
            print(f"Temp file initialized: {self.temp_file}")
        except Exception as e:
            print(f"Error creating temp file: {e}")
            # Try fallback without temp file
            self.temp_file = None
            print("Disabled temp file - will save directly to final file")
    
    def add_result(self, gene: str, species: str, count: int):
        """Add a result and save to temp file"""
        with self.lock:
            self.results.append((gene, species, count))
            
            # Append to temp file immediately if available
            if self.temp_file:
                try:
                    with self.temp_file.open('a') as f:
                        f.write(f"{gene}\t{species}\t{count}\n")
                except Exception as e:
                    print(f"Error writing to temp file: {e}")
                    # Disable temp file if it's causing issues
                    self.temp_file = None
    
    def finalize(self):
        """Move temp file to final location or write results directly"""
        with self.lock:
            try:
                if self.temp_file and self.temp_file.exists():
                    # Move temp file to final location
                    self.temp_file.rename(self.output_file)
                    print(f"Results finalized from temp file: {self.output_file}")
                else:
                    # Write results directly to final file
                    with self.output_file.open('w') as f:
                        f.write("gene\tspecies\tcount\n")
                        for gene, species, count in self.results:
                            f.write(f"{gene}\t{species}\t{count}\n")
                    print(f"Results written directly to final file: {self.output_file}")
            except Exception as e:
                print(f"Error finalizing results: {e}")
                # Last resort - try to write to current directory
                try:
                    fallback_file = Path(self.output_file.name)
                    with fallback_file.open('w') as f:
                        f.write("gene\tspecies\tcount\n")
                        for gene, species, count in self.results:
                            f.write(f"{gene}\t{species}\t{count}\n")
                    print(f"Results written to fallback location: {fallback_file}")
                except Exception as e2:
                    print(f"All file writing methods failed: {e2}")

def setup_signal_handlers(cache_manager, results_saver):
    """Setup signal handlers for graceful shutdown"""
    def signal_handler(signum, frame):
        print("\n\nReceived interrupt signal. Saving progress...")
        cache_manager.save_cache(force=True)
        results_saver.finalize()
        print("Progress saved. Exiting.")
        sys.exit(0)
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

def main():
    """Main function with comprehensive monitoring"""
    print("Starting Enhanced Gene Taxa Coverage Analysis...")
    
    # Get inputs from Snakemake
    try:
        species_file = snakemake.input.species
        genes_file = snakemake.input.genes
        credentials_file = snakemake.input.ncbi_info
        output_file = snakemake.output.coverage
        cache_dir = snakemake.params.cache_dir
        print(f"Using Snakemake inputs - output: {output_file}")
    except NameError:
        print("Snakemake object not available, using test values")
        species_file = "data/bacdive/analysis_1/gram_positive.txt"
        genes_file = "data/quickgo/params_1/gene_symbols_filtered.txt"
        credentials_file = "config/login/ncbi_info.txt"
        output_file = "results/coverage/test_coverage.tsv"
        cache_dir = "results/cache"
    
    # Setup NCBI credentials
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
    try:
        with open(species_file) as f:
            species_list = [line.strip() for line in f if line.strip()]
        
        with open(genes_file) as f:
            gene_list = [line.strip() for line in f if line.strip()]
            
        print(f"Loaded {len(species_list)} species and {len(gene_list)} genes")
    except Exception as e:
        print(f"Error loading input files: {e}")
        return
    
    # Initialize monitoring and management systems
    total_combinations = len(species_list) * len(gene_list)
    progress = ProgressMonitor(total_combinations, len(species_list), len(gene_list))
    cache_manager = RobustCacheManager(cache_dir)
    results_saver = IntermediateResultsSaver(output_file)
    
    # Setup signal handlers for graceful shutdown
    setup_signal_handlers(cache_manager, results_saver)
    
    print("Starting NCBI queries with progress monitoring...")
    
    # Process combinations with monitoring
    for i, species in enumerate(species_list):
        print(f"\nProcessing species {i+1}/{len(species_list)}: {species}")
        
        for j, gene in enumerate(gene_list):
            cache_key = f"{gene}|{species}"
            
            # Check cache first
            cached_result = cache_manager.get(cache_key)
            if cached_result:
                count = cached_result["count"]
                results_saver.add_result(gene, species, count)
                progress.update(cached=True)
                continue
            
            # Throttle if necessary
            if progress.should_throttle():
                print("Rate limiting: waiting 60 seconds...")
                time.sleep(60)
            
            # Query NCBI
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
                
                results_saver.add_result(gene, species, count)
                progress.update(cached=False)
                
                # Respectful delay
                time.sleep(0.1)
                
            except Exception as e:
                print(f"Error querying {gene} in {species}: {e}")
                results_saver.add_result(gene, species, 0)
                progress.update(cached=False, error=True)
                time.sleep(1)  # Longer delay after error
        
        # Save cache and progress after each species
        cache_manager.save_cache()
        print(f"Completed species {i+1}/{len(species_list)}")
    
    # Finalize results
    results_saver.finalize()
    cache_manager.save_cache(force=True)
    progress.final_report()
    
    print(f"Analysis complete! Results saved to {output_file}")

if __name__ == "__main__":
    main()