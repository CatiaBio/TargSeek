#!/usr/bin/env python3
"""
Cached Protein FASTA Download System
====================================

This script provides a caching system for protein FASTA downloads to avoid
re-downloading the same sequences. Uses a similar approach to the gene-taxa
coverage cache.

Cache structure:
- Key: "gene||species||database" (e.g., "bamA||Escherichia coli||uniprot")  
- Value: {"sequence": "MKLVL...", "accession": "P0A937", "found": true/false}
"""

import json
import requests
import time
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import hashlib
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Cache configuration
CACHE_DIR = Path("cache/protein_sequences")
CACHE_FILE = CACHE_DIR / "protein_fasta_cache.json"

class ProteinSequenceCache:
    """Manages caching of protein sequences"""
    
    def __init__(self):
        self.cache = {}
        self.cache_updates = 0
        self.load_cache()
    
    def load_cache(self):
        """Load existing protein sequence cache"""
        if CACHE_FILE.exists():
            try:
                with open(CACHE_FILE, 'r') as f:
                    self.cache = json.load(f)
                logging.debug(f"✓ Loaded protein cache with {len(self.cache)} entries")
            except Exception as e:
                logging.warning(f"Could not load protein cache: {e}")
                self.cache = {}
        else:
            logging.info("No existing protein cache found, starting fresh")
            self.cache = {}
    
    def save_cache(self):
        """Save protein sequence cache"""
        try:
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with open(CACHE_FILE, 'w') as f:
                json.dump(self.cache, f, indent=2)
            logging.info(f"✓ Saved protein cache with {len(self.cache)} entries")
        except Exception as e:
            logging.warning(f"Could not save protein cache: {e}")
    
    def get_cache_key(self, gene: str, species: str, database: str = "uniprot"):
        """Generate cache key for gene-species-database combination"""
        return f"{gene}||{species}||{database}"
    
    def get_cached_sequence(self, gene: str, species: str, database: str = "uniprot"):
        """Get cached sequence if available"""
        cache_key = self.get_cache_key(gene, species, database)
        cached_data = self.cache.get(cache_key)
        
        # If not found, return None
        if cached_data is None:
            return None
            
        # If found and it's a positive result, return it
        if cached_data.get("found", False):
            return cached_data
            
        # If it's a negative result (not found), check if it's older than 6 months
        if not cached_data.get("found", False):
            cached_time = cached_data.get("cached_at")
            if cached_time:
                try:
                    from datetime import datetime
                    cache_date = datetime.fromisoformat(cached_time)
                    age_days = (datetime.now() - cache_date).days
                    
                    # If older than 6 months (180 days), treat as cache miss
                    if age_days > 180:
                        logging.info(f"Cache entry for {gene}||{species} is {age_days} days old (expired), will retry")
                        # Remove expired entry
                        del self.cache[cache_key]
                        self.cache_updates += 1
                        return None
                except:
                    pass
                    
        return cached_data
    
    def cache_sequence(self, gene: str, species: str, sequence_data: Dict, database: str = "uniprot"):
        """Cache a sequence result (including negative results)"""
        cache_key = self.get_cache_key(gene, species, database)
        self.cache[cache_key] = {
            **sequence_data,
            "cached_at": datetime.now().isoformat(),
            "gene": gene,
            "species": species,
            "database": database
        }
        self.cache_updates += 1
        
        # Save cache every 50 updates
        if self.cache_updates % 50 == 0:
            self.save_cache()

class CachedProteinDownloader:
    """Enhanced protein downloader with caching"""
    
    def __init__(self):
        self.cache = ProteinSequenceCache()
        self.stats = {
            "cache_hits": 0,
            "cache_misses": 0, 
            "sequences_found": 0,
            "sequences_not_found": 0,
            "api_calls": 0
        }
    
    def search_uniprot_cached(self, gene: str, species: str) -> Optional[Dict[str, Any]]:
        """
        Search UniProt for gene-species combination with caching
        
        Returns:
            Dictionary with sequence data or None if not found
        """
        
        # Check cache first
        cached_result = self.cache.get_cached_sequence(gene, species, "uniprot")
        if cached_result is not None:
            self.stats["cache_hits"] += 1
            if cached_result.get("found", False):
                self.stats["sequences_found"] += 1
                return cached_result
            else:
                # Negative result was not expired (otherwise get_cached_sequence would return None)
                self.stats["sequences_not_found"] += 1
                return None
        
        # Not in cache, query UniProt
        self.stats["cache_misses"] += 1
        self.stats["api_calls"] += 1
        
        try:
            # Build UniProt query
            query = f'gene:{gene} AND organism_name:"{species}"'
            
            base_url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                'query': query,
                'format': 'json',
                'size': 1,  # Just need the best match
                'fields': 'accession,organism_name,gene_names,protein_name,sequence'
            }
            
            response = requests.get(base_url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            results = data.get('results', [])
            
            if results:
                # Found sequence
                entry = results[0]
                sequence_data = {
                    "found": True,
                    "accession": entry.get('primaryAccession', ''),
                    "organism": entry.get('organism', {}).get('scientificName', species),
                    "protein_name": entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    "sequence": entry.get('sequence', {}).get('value', ''),
                    "length": entry.get('sequence', {}).get('length', 0)
                }
                
                self.cache.cache_sequence(gene, species, sequence_data, "uniprot")
                self.stats["sequences_found"] += 1
                return sequence_data
            else:
                # Not found
                negative_result = {
                    "found": False,
                    "reason": "No UniProt entry found"
                }
                
                self.cache.cache_sequence(gene, species, negative_result, "uniprot")
                self.stats["sequences_not_found"] += 1
                return None
                
        except requests.exceptions.RequestException as e:
            logging.warning(f"UniProt search failed for {gene} in {species}: {e}")
            # Don't cache failures due to network issues
            self.stats["sequences_not_found"] += 1
            return None
        except Exception as e:
            logging.error(f"Error searching UniProt for {gene} in {species}: {e}")
            self.stats["sequences_not_found"] += 1
            return None
    
    def download_gene_sequences(self, gene: str, species_list: List[str], output_dir: Path) -> Dict[str, Any]:
        """
        Download sequences for a gene across multiple species with caching
        
        Args:
            gene: Gene name to search for
            species_list: List of species names
            output_dir: Output directory for FASTA files
            
        Returns:
            Summary dictionary with download statistics
        """
        
        logging.info(f"Downloading sequences for {gene} across {len(species_list)} species")
        
        # Create gene-specific directory
        gene_dir = output_dir / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        
        sequences = []
        found_species = []
        not_found_species = []
        
        cache_hits_in_batch = 0
        new_lookups_in_batch = 0
        
        for i, species in enumerate(species_list, 1):
            # Check if this will be a cache hit
            cache_key = self.cache.get_cache_key(gene, species, "uniprot")
            is_cached = cache_key in self.cache.cache
            
            if is_cached:
                cache_hits_in_batch += 1
            else:
                new_lookups_in_batch += 1
                logging.info(f"  Processing {species} ({i}/{len(species_list)}) - NEW LOOKUP")
            
            # Search for sequence (with caching)
            seq_data = self.search_uniprot_cached(gene, species)
            
            if seq_data:
                # Create FASTA record
                safe_species = species.replace(' ', '_').replace('/', '_')
                seq_record = SeqRecord(
                    Seq(seq_data['sequence']),
                    id=f"{seq_data['accession']}_{safe_species}",
                    description=f"{seq_data['protein_name']} [{species}] {gene}"
                )
                
                sequences.append(seq_record)
                found_species.append(species)
                
                # Save individual FASTA file
                individual_file = gene_dir / f"{safe_species}.fasta"
                with open(individual_file, 'w') as f:
                    SeqIO.write(seq_record, f, "fasta")
                    
            else:
                not_found_species.append(species)
            
            # Rate limiting only for new lookups
            if not is_cached:
                time.sleep(0.2)
            
            # Progress update every 50 species when using cache
            if is_cached and i % 50 == 0:
                logging.info(f"  Progress: {i}/{len(species_list)} species ({cache_hits_in_batch} from cache, {new_lookups_in_batch} new lookups)")
        
        # Save combined FASTA file
        if sequences:
            combined_file = gene_dir / f"{gene}_all_sequences.fasta"
            with open(combined_file, 'w') as f:
                SeqIO.write(sequences, f, "fasta")
            logging.info(f"✓ Saved {len(sequences)} sequences to {combined_file}")
        
        # Save not found species
        if not_found_species:
            not_found_file = gene_dir / f"{gene}_not_found.txt"
            with open(not_found_file, 'w') as f:
                f.write(f"# Species without {gene} sequences\n")
                f.write(f"# Generated: {datetime.now().isoformat()}\n")
                for species in not_found_species:
                    f.write(f"{species}\n")
        
        # Create summary
        summary = {
            "gene": gene,
            "total_species": len(species_list),
            "sequences_found": len(sequences),
            "sequences_not_found": len(not_found_species),
            "success_rate": len(sequences) / len(species_list) * 100,
            "found_species": found_species,
            "not_found_species": not_found_species,
            "output_directory": str(gene_dir)
        }
        
        # Save summary
        summary_file = gene_dir / f"{gene}_download_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return summary
    
    def get_cache_stats(self):
        """Get cache performance statistics"""
        total_requests = self.stats["cache_hits"] + self.stats["cache_misses"]
        if total_requests > 0:
            cache_hit_rate = (self.stats["cache_hits"] / total_requests) * 100
        else:
            cache_hit_rate = 0
            
        return {
            **self.stats,
            "cache_entries": len(self.cache.cache),
            "cache_hit_rate": cache_hit_rate
        }
    
    def finalize(self):
        """Save cache and print final statistics"""
        self.cache.save_cache()
        
        stats = self.get_cache_stats()
        logging.info("=== Download Statistics ===")
        logging.info(f"Cache entries: {stats['cache_entries']}")
        logging.info(f"Cache hits: {stats['cache_hits']}")
        logging.info(f"Cache misses: {stats['cache_misses']}")
        logging.info(f"Cache hit rate: {stats['cache_hit_rate']:.1f}%")
        logging.info(f"API calls made: {stats['api_calls']}")
        logging.info(f"Sequences found: {stats['sequences_found']}")
        logging.info(f"Sequences not found: {stats['sequences_not_found']}")

def main():
    """Test the cached downloader"""
    if len(sys.argv) != 4:
        print("Usage: python download_proteins_cached.py <gene> <species_file> <output_dir>")
        sys.exit(1)
    
    gene = sys.argv[1]
    species_file = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    
    # Load species list
    with open(species_file, 'r') as f:
        species_list = [line.strip() for line in f if line.strip()]
    
    # Download with caching
    downloader = CachedProteinDownloader()
    summary = downloader.download_gene_sequences(gene, species_list, output_dir)
    downloader.finalize()
    
    print(f"\\nDownload completed: {summary['sequences_found']}/{summary['total_species']} sequences found")

if __name__ == "__main__":
    import sys
    main()