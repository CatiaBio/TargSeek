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
CACHE_FILE = CACHE_DIR / "protein_structures_cache.json"

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
        
        # Count PDB files in this gene directory
        pdb_files = list(gene_dir.glob("*.pdb"))
        if not pdb_files:
            logging.info(f"  No PDB files found in {gene_name}")
            continue
            
        logging.info(f"  Found {len(pdb_files)} PDB files in {gene_name}")
        
        # Process each PDB file
        for pdb_file in pdb_files:
            if pdb_file.stat().st_size == 0:
                logging.warning(f"  Skipping empty file: {pdb_file.name}")
                continue
                
            # Extract species name from filename (convert back from safe filename)
            species_name = pdb_file.stem.replace('_', ' ')
            
            # Generate cache key
            cache_key = f"{gene_name}||{species_name}||pdb_search"
            
            # Check if already in cache
            if cache_key in cache:
                logging.debug(f"  Already cached: {gene_name} - {species_name}")
                continue
            
            # Try to extract PDB ID from file content
            pdb_id = None
            uniprot_acc = "unknown"
            try:
                with open(pdb_file, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith('HEADER'):
                        # PDB files have the ID in the header line
                        parts = first_line.split()
                        if len(parts) >= 4:
                            pdb_id = parts[-1].lower()
                        logging.debug(f"    Extracted PDB ID from header: {pdb_id}")
                    else:
                        # Try to find PDB ID in first few lines
                        f.seek(0)
                        for line_num, line in enumerate(f):
                            if line_num > 10:  # Don't search too far
                                break
                            if 'PDB' in line.upper():
                                # Try to extract 4-character PDB ID pattern
                                import re
                                pdb_match = re.search(r'\b[A-Za-z0-9]{4}\b', line)
                                if pdb_match:
                                    pdb_id = pdb_match.group().lower()
                                    logging.debug(f"    Found potential PDB ID: {pdb_id}")
                                    break
            except Exception as e:
                logging.debug(f"    Could not read file {pdb_file}: {e}")
            
            # If no PDB ID found, create a placeholder
            if not pdb_id:
                pdb_id = f"existing_{pdb_file.stem}"
                logging.debug(f"    Using placeholder PDB ID: {pdb_id}")
            
            # Create cache entry
            cache_entry = {
                "found": True,
                "pdb_ids": [pdb_id],
                "uniprot_acc": uniprot_acc,
                "found_via_gene": gene_name,
                "cached_from_existing": True,
                "file_path": str(pdb_file),
                "file_size": pdb_file.stat().st_size,
                "cached_at": datetime.now().isoformat(),
                "gene": gene_name,
                "species": species_name,
                "initialization_run": True
            }
            
            # Add to cache
            cache[cache_key] = cache_entry
            cached_count += 1
            logging.info(f"  ✓ Cached: {gene_name} - {species_name} ({pdb_file.name}, {cache_entry['file_size']} bytes)")
    
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