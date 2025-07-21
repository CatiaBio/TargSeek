#!/usr/bin/env python3
"""
Migrate existing results to new shared data structure
===================================================

This script moves existing files from the old results/ structure to the new shared data structure:

Old structure:
- results/protein_fasta/{analysis}_{paramset}_gram_{group}/{gene}/{species}.fasta
- results/3d_structures/{analysis}_{paramset}_gram_{group}/{gene}/...
- results/proteins_to_download/{analysis}_{paramset}_gram_{group}/{gene}.txt

New structure:
- data/proteins_fasta/{gene}/{species}.fasta + metadata.json
- data/proteins_3d_structure/{gene}/{species}.pdb + metadata.json  
- data/{analysis}_{paramset}/genes_species/gram_{group}/{gene}.txt
"""

import sys
import json
import shutil
import logging
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def create_directories():
    """Create the new directory structure"""
    dirs_to_create = [
        "data/proteins_fasta",
        "data/proteins_3d_structure",
        "data/analysis_1_params_1/genes_species/gram_positive",
        "data/analysis_1_params_1/genes_species/gram_negative",
        "cache/protein_sequences",
        "cache/gene_species"
    ]
    
    for dir_path in dirs_to_create:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        logging.info(f"Created directory: {dir_path}")

def migrate_protein_fasta():
    """Migrate protein FASTA files to shared structure"""
    logging.info("=== Migrating Protein FASTA Files ===")
    
    old_protein_dir = Path("results/protein_fasta")
    new_protein_dir = Path("data/proteins_fasta")
    
    if not old_protein_dir.exists():
        logging.info("No old protein_fasta directory found, skipping")
        return
    
    # Track gene metadata
    gene_metadata = defaultdict(lambda: {
        "available_species": set(),
        "failed_species": set(),
        "sources": []
    })
    
    # Process each analysis directory
    for analysis_dir in old_protein_dir.iterdir():
        if not analysis_dir.is_dir():
            continue
            
        logging.info(f"Processing analysis directory: {analysis_dir.name}")
        
        # Process each gene directory
        for gene_dir in analysis_dir.iterdir():
            if not gene_dir.is_dir() or gene_dir.name.startswith('_'):
                continue
                
            gene_name = gene_dir.name
            logging.info(f"  Migrating gene: {gene_name}")
            
            # Create gene directory in shared structure
            shared_gene_dir = new_protein_dir / gene_name
            shared_gene_dir.mkdir(exist_ok=True)
            
            # Count files moved
            files_moved = 0
            
            # Move all FASTA files
            for fasta_file in gene_dir.glob("*.fasta"):
                if fasta_file.name in ["not_found_species.txt", "download_summary.json", "final_download_summary.json"]:
                    continue
                    
                # Extract species name from filename
                species_name = fasta_file.stem.replace('_', ' ')
                target_file = shared_gene_dir / fasta_file.name
                
                # Only move if target doesn't exist or is smaller
                if not target_file.exists() or target_file.stat().st_size < fasta_file.stat().st_size:
                    shutil.copy2(fasta_file, target_file)
                    files_moved += 1
                    gene_metadata[gene_name]["available_species"].add(species_name)
                    gene_metadata[gene_name]["sources"].append(f"{analysis_dir.name}/{gene_name}")
            
            # Copy failure information
            not_found_file = gene_dir / "not_found_species.txt"
            if not_found_file.exists():
                with open(not_found_file, 'r') as f:
                    failed_species = [line.strip() for line in f if line.strip() and not line.startswith('#')]
                    gene_metadata[gene_name]["failed_species"].update(failed_species)
            
            logging.info(f"    Moved {files_moved} FASTA files for {gene_name}")
    
    # Create metadata files for each gene
    logging.info("Creating gene metadata files...")
    for gene_name, metadata in gene_metadata.items():
        gene_dir = new_protein_dir / gene_name
        metadata_file = gene_dir / f"{gene_name}_metadata.json"
        
        # Convert sets to sorted lists
        available_species = sorted(list(metadata["available_species"]))
        failed_species = sorted(list(metadata["failed_species"]))
        
        # Create metadata
        gene_meta = {
            "gene": gene_name,
            "last_updated": datetime.now().isoformat(),
            "available_species": available_species,
            "failed_species": failed_species,
            "total_species": len(available_species) + len(failed_species),
            "success_rate": len(available_species) / (len(available_species) + len(failed_species)) * 100 if (len(available_species) + len(failed_species)) > 0 else 0,
            "migration_sources": list(set(metadata["sources"])),
            "migrated_at": datetime.now().isoformat()
        }
        
        with open(metadata_file, 'w') as f:
            json.dump(gene_meta, f, indent=2)
        
        logging.info(f"  Created metadata for {gene_name}: {len(available_species)} available, {len(failed_species)} failed")

def migrate_3d_structures():
    """Migrate 3D structure files to shared structure"""
    logging.info("=== Migrating 3D Structure Files ===")
    
    old_structures_dir = Path("results/3d_structures")
    new_structures_dir = Path("data/proteins_3d_structure")
    
    if not old_structures_dir.exists():
        logging.info("No old 3d_structures directory found, skipping")
        return
    
    # Track gene structure metadata
    gene_metadata = defaultdict(lambda: {
        "available_species": set(),
        "failed_species": set(),
        "sources": []
    })
    
    # Process each analysis directory
    for analysis_dir in old_structures_dir.iterdir():
        if not analysis_dir.is_dir():
            continue
            
        logging.info(f"Processing structures from: {analysis_dir.name}")
        
        # Process each gene directory  
        for gene_dir in analysis_dir.iterdir():
            if not gene_dir.is_dir():
                continue
                
            gene_name = gene_dir.name
            logging.info(f"  Migrating structures for gene: {gene_name}")
            
            # Create gene directory in shared structure
            shared_gene_dir = new_structures_dir / gene_name
            shared_gene_dir.mkdir(exist_ok=True)
            
            files_moved = 0
            
            # Move PDB and FASTA files, but rename appropriately
            for structure_file in gene_dir.iterdir():
                if structure_file.is_file():
                    if structure_file.suffix in ['.pdb', '.gz', '.fasta']:
                        # Handle different file types
                        if structure_file.name == "no_structures_found.txt":
                            continue
                            
                        # For complex structure files like "7BIN_chain_1.fasta", we need to be more careful
                        if '_chain_' in structure_file.name:
                            # Skip chain-specific files for now, might need different handling
                            continue
                            
                        # Simple case: files like "2QCZ.pdb.gz" or "2QCZ.fasta"  
                        if structure_file.suffix == '.gz' and structure_file.name.endswith('.pdb.gz'):
                            # Convert to .pdb file
                            target_name = structure_file.name.replace('.pdb.gz', '.pdb')
                            target_file = shared_gene_dir / target_name
                            
                            if not target_file.exists():
                                # Copy and potentially decompress
                                shutil.copy2(structure_file, target_file)
                                files_moved += 1
                        elif structure_file.suffix == '.fasta':
                            # Structure FASTA file - we might want to handle differently
                            target_file = shared_gene_dir / structure_file.name
                            if not target_file.exists():
                                shutil.copy2(structure_file, target_file)
                                files_moved += 1
            
            logging.info(f"    Moved {files_moved} structure files for {gene_name}")
    
    logging.info("3D structures migration completed")

def migrate_gene_species_lists():
    """Migrate gene species lists to analysis-specific directories"""
    logging.info("=== Migrating Gene Species Lists ===")
    
    old_lists_dir = Path("results/proteins_to_download")
    
    if not old_lists_dir.exists():
        logging.info("No old proteins_to_download directory found, skipping")
        return
    
    # Process each analysis directory
    for analysis_dir in old_lists_dir.iterdir():
        if not analysis_dir.is_dir():
            continue
            
        # Parse analysis name to determine target directory
        # Expected format: analysis_1_params_1_gram_positive
        parts = analysis_dir.name.split('_')
        if len(parts) >= 5 and parts[2] == 'params':
            analysis = f"{parts[0]}_{parts[1]}"  # analysis_1
            paramset = f"{parts[2]}_{parts[3]}"  # params_1
            gram_type = f"gram_{'_'.join(parts[5:])}"  # gram_positive or gram_negative
            
            # Create target directory
            target_dir = Path(f"data/{analysis}_{paramset}/genes_species/{gram_type}")
            target_dir.mkdir(parents=True, exist_ok=True)
            
            logging.info(f"Migrating {analysis_dir.name} -> {target_dir}")
            
            # Copy all .txt files (gene species lists)
            files_moved = 0
            for gene_file in analysis_dir.glob("*.txt"):
                if not gene_file.name.startswith('_'):  # Skip summary files
                    target_file = target_dir / gene_file.name
                    shutil.copy2(gene_file, target_file)
                    files_moved += 1
            
            # Copy summary file if it exists
            summary_file = analysis_dir / "_gene_summary.txt"
            if summary_file.exists():
                target_summary = target_dir / "_gene_summary.txt"
                shutil.copy2(summary_file, target_summary)
            
            logging.info(f"  Moved {files_moved} gene species list files")

def create_initial_caches():
    """Create initial cache files from migrated data"""
    logging.info("=== Creating Initial Cache Files ===")
    
    # Create protein sequence cache from migrated FASTA files
    protein_cache = {}
    protein_dir = Path("data/proteins_fasta")
    
    if protein_dir.exists():
        for gene_dir in protein_dir.iterdir():
            if gene_dir.is_dir():
                gene_name = gene_dir.name
                
                # Read each species FASTA file and create cache entries
                for fasta_file in gene_dir.glob("*.fasta"):
                    if not fasta_file.name.endswith('_metadata.json'):
                        species_name = fasta_file.stem.replace('_', ' ')
                        
                        try:
                            # Read FASTA file to extract sequence info
                            with open(fasta_file, 'r') as f:
                                lines = f.readlines()
                            
                            if len(lines) >= 2:
                                header = lines[0].strip()
                                sequence = ''.join(line.strip() for line in lines[1:])
                                
                                # Parse header to get accession
                                accession = ""
                                if header.startswith('>'):
                                    parts = header[1:].split('|')
                                    if len(parts) > 0:
                                        accession = parts[0]
                                
                                # Create cache key and entry
                                cache_key = f"{gene_name}||{species_name}||uniprot"
                                cache_entry = {
                                    "found": True,
                                    "accession": accession,
                                    "organism": species_name,
                                    "sequence": sequence,
                                    "length": len(sequence),
                                    "cached_at": datetime.now().isoformat(),
                                    "gene": gene_name,
                                    "species": species_name,
                                    "database": "uniprot",
                                    "migrated": True
                                }
                                
                                protein_cache[cache_key] = cache_entry
                        except Exception as e:
                            logging.warning(f"Could not process {fasta_file}: {e}")
    
    # Save protein cache
    cache_dir = Path("cache/protein_sequences")
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "protein_fasta_cache.json"
    
    with open(cache_file, 'w') as f:
        json.dump(protein_cache, f, indent=2)
    
    logging.info(f"Created protein sequence cache with {len(protein_cache)} entries")
    
    # Create gene-species coverage cache (simpler version)
    coverage_cache = {}
    
    # For each gene-species combination we have FASTA for, mark as found
    for cache_key, cache_entry in protein_cache.items():
        gene, species, database = cache_key.split('||')
        coverage_key = f"{species}||{gene}"
        coverage_cache[coverage_key] = True  # Found
    
    # Save coverage cache
    coverage_cache_file = Path("cache/gene_species/species_gene_cache.json")
    coverage_cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(coverage_cache_file, 'w') as f:
        json.dump(coverage_cache, f, indent=2)
    
    logging.info(f"Created gene-species coverage cache with {len(coverage_cache)} entries")

def main():
    """Main migration function"""
    logging.info("=== Starting Migration to Shared Data Structure ===")
    
    try:
        # Create new directory structure
        create_directories()
        
        # Migrate protein FASTA files
        migrate_protein_fasta()
        
        # Migrate 3D structure files
        migrate_3d_structures()
        
        # Migrate gene species lists
        migrate_gene_species_lists()
        
        # Create initial cache files
        create_initial_caches()
        
        logging.info("=== Migration Completed Successfully ===")
        logging.info("New structure created:")
        logging.info("  - data/proteins_fasta/{gene}/ - Shared protein sequences")
        logging.info("  - data/proteins_3d_structure/{gene}/ - Shared 3D structures")
        logging.info("  - data/analysis_1_params_1/genes_species/gram_{group}/ - Analysis-specific gene lists")
        logging.info("  - cache/protein_sequences/ - Protein sequence cache")
        logging.info("  - cache/gene_species/ - Gene-species coverage cache")
        
    except Exception as e:
        logging.error(f"Migration failed: {e}")
        raise

if __name__ == "__main__":
    main()