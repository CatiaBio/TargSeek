#!/usr/bin/env python3
"""
Initialize 3D Structure Cache from Existing Files
=================================================

This script scans existing 3D structure files and populates the cache
to avoid re-downloading structures that already exist.

Usage: python initialize_3d_structure_cache.py
"""

import json
import logging
from pathlib import Path
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Cache configuration
CACHE_DIR = Path("cache/protein_structures")
CACHE_FILE = CACHE_DIR / "protein_3d_structure_cache.json"

def initialize_3d_structure_cache():
    """
    Scan existing 3D structure files and populate cache
    """
    
    logging.info("=== Initializing 3D Structure Cache from Existing Files ===")
    
    # Load existing cache if it exists
    cache = {}
    if CACHE_FILE.exists():
        try:
            with open(CACHE_FILE, 'r') as f:
                cache = json.load(f)
            logging.info(f"Loaded existing cache with {len(cache)} entries")
        except Exception as e:
            logging.warning(f"Could not load existing cache: {e}")
            cache = {}
    else:
        logging.info("No existing cache found, creating new cache")
    
    # Look for existing structure files in the shared directory
    shared_structures_dir = Path("data/proteins_3d_structure")
    if not shared_structures_dir.exists():
        logging.warning(f"3D structures directory not found: {shared_structures_dir}")
        return 0
    
    logging.info(f"Scanning 3D structure files in: {shared_structures_dir}")
    cached_count = 0
    updated_count = 0
    
    # Scan each gene directory
    for gene_dir in shared_structures_dir.iterdir():
        if not gene_dir.is_dir() or gene_dir.name.startswith('.'):
            continue
            
        gene_name = gene_dir.name
        logging.info(f"Processing gene directory: {gene_name}")
        
        # Count PDB files in this gene directory (both .pdb and .pdb.gz)
        pdb_files = list(gene_dir.glob("*.pdb")) + list(gene_dir.glob("*.pdb.gz"))
        if not pdb_files:
            logging.info(f"  No PDB files found in {gene_name}")
            continue
            
        logging.info(f"  Found {len(pdb_files)} PDB files in {gene_name}")
        
        # Check if there are any PDB files for this gene
        pdb_ids_found = []
        for pdb_file in pdb_files:
            if pdb_file.stat().st_size == 0:
                logging.warning(f"  Skipping empty file: {pdb_file.name}")
                continue
                
            # Extract PDB ID from filename (remove .pdb or .pdb.gz extension)
            filename = pdb_file.name
            if filename.endswith('.pdb.gz'):
                pdb_id = filename[:-7].upper()  # Remove .pdb.gz
            elif filename.endswith('.pdb'):
                pdb_id = filename[:-4].upper()  # Remove .pdb
            else:
                continue
                
            pdb_ids_found.append(pdb_id)
            
            # Create PDB download cache entry
            pdb_cache_key = f"pdb||{pdb_id}||structure"
            if pdb_cache_key not in cache:
                cache[pdb_cache_key] = {
                    "success": True,
                    "file_path": str(pdb_file),
                    "sequences_count": 1,  # Estimate
                    "download_date": datetime.now().isoformat(),
                    "cached_from_existing": True
                }
                cached_count += 1
                logging.info(f"  ✓ Cached PDB: {pdb_id} ({pdb_file.name})")
        
        # Create gene search cache entry
        if pdb_ids_found:
            gene_cache_key = f"{gene_name}||uniprot||bacterial"
            if gene_cache_key not in cache:
                cache[gene_cache_key] = {
                    "pdb_ids": pdb_ids_found,
                    "found": True,
                    "search_date": datetime.now().isoformat(),
                    "database": "uniprot",
                    "cached_from_existing": True
                }
                cached_count += 1
                logging.info(f"  ✓ Cached gene search: {gene_name} with {len(pdb_ids_found)} structures")
    
    # Save updated cache
    if cached_count > 0 or updated_count > 0:
        try:
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with open(CACHE_FILE, 'w') as f:
                json.dump(cache, f, indent=2)
            logging.info(f"✅ Successfully saved cache to: {CACHE_FILE}")
        except Exception as e:
            logging.error(f"Failed to save cache: {e}")
            return 0
    
    # Print summary
    logging.info(f"\n=== Cache Initialization Complete ===")
    logging.info(f"Total cache entries: {len(cache)}")
    logging.info(f"New entries added: {cached_count}")
    logging.info(f"Cache file: {CACHE_FILE}")
    
    if cached_count > 0:
        logging.info(f"\n✅ Cache successfully initialized with {cached_count} existing 3D structure files")
        logging.info("Future pipeline runs will skip downloading these structures")
    else:
        logging.info("No new structure files found to cache")
    
    return cached_count

def main():
    """Main function"""
    try:
        cached_count = initialize_3d_structure_cache()
        
        if cached_count > 0:
            print(f"\n🎉 Success! Initialized cache with {cached_count} existing 3D structure files")
            print("You can now run the download_3d_structures rule and it will skip these files")
        else:
            print("\nNo 3D structure files found to cache")
            print("Cache is ready for future downloads")
            
    except Exception as e:
        logging.error(f"Error during cache initialization: {e}")
        print(f"\n❌ Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())