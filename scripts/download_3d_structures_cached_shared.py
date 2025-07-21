#!/usr/bin/env python3
"""
Download 3D structures to shared data directory with caching
==========================================================

This script downloads 3D structures to a shared data/proteins_3d_structure/ directory
organized by gene, regardless of Gram classification. Uses caching and alias support.

Directory structure:
- data/proteins_3d_structure/{gene}/{species}.pdb
- data/proteins_3d_structure/{gene}/{gene}_structures_metadata.json

Based on the original download_3d_structures.py but adapted for shared directory structure.
"""

import sys
import json
import logging
from pathlib import Path
from datetime import datetime

# Add scripts directory to path for imports
sys.path.append(str(Path(__file__).parent))

# Import from the original 3D structures script
import requests
import time
from Bio.PDB import PDBParser
from Bio import PDB
import warnings
warnings.filterwarnings("ignore", category=PDB.PDBExceptions.PDBConstructionWarning)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Cache configuration
CACHE_DIR = Path("cache/protein_structures")
CACHE_FILE = CACHE_DIR / "protein_structures_cache.json"

class ProteinStructureCache:
    """Manages caching of protein 3D structure searches"""
    
    def __init__(self):
        self.cache = {}
        self.cache_updates = 0
        self.load_cache()
    
    def load_cache(self):
        """Load existing protein structure cache"""
        if CACHE_FILE.exists():
            try:
                with open(CACHE_FILE, 'r') as f:
                    self.cache = json.load(f)
                logging.info(f"✓ Loaded 3D structure cache with {len(self.cache)} entries")
            except Exception as e:
                logging.warning(f"Could not load 3D structure cache: {e}")
                self.cache = {}
        else:
            logging.info("No existing 3D structure cache found, starting fresh")
            self.cache = {}
    
    def save_cache(self):
        """Save protein structure cache"""
        try:
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with open(CACHE_FILE, 'w') as f:
                json.dump(self.cache, f, indent=2)
            logging.info(f"✓ Saved 3D structure cache with {len(self.cache)} entries")
        except Exception as e:
            logging.warning(f"Could not save 3D structure cache: {e}")
    
    def get_cache_key(self, gene: str, species: str):
        """Generate cache key for gene-species combination"""
        return f"{gene}||{species}||pdb_search"
    
    def get_cached_result(self, gene: str, species: str):
        """Get cached search result if available"""
        cache_key = self.get_cache_key(gene, species)
        return self.cache.get(cache_key)
    
    def cache_result(self, gene: str, species: str, search_result: dict):
        """Cache a search result (including negative results)"""
        cache_key = self.get_cache_key(gene, species)
        self.cache[cache_key] = {
            **search_result,
            "cached_at": datetime.now().isoformat(),
            "gene": gene,
            "species": species
        }
        self.cache_updates += 1
        
        # Save cache every 25 updates
        if self.cache_updates % 25 == 0:
            self.save_cache()

class CachedStructureDownloader:
    """Enhanced 3D structure downloader with caching"""
    
    def __init__(self):
        self.cache = ProteinStructureCache()
        self.stats = {
            "cache_hits": 0,
            "cache_misses": 0, 
            "structures_found": 0,
            "structures_not_found": 0,
            "api_calls": 0
        }
    
    def finalize(self):
        """Finalize downloads and save cache"""
        self.cache.save_cache()
    
    def get_cache_stats(self):
        """Get download and cache statistics"""
        total_requests = self.stats["cache_hits"] + self.stats["cache_misses"]
        cache_hit_rate = (self.stats["cache_hits"] / total_requests * 100) if total_requests > 0 else 0
        
        return {
            "cache_hit_rate": cache_hit_rate,
            "api_calls": self.stats["api_calls"],
            "structures_found": self.stats["structures_found"],
            "structures_not_found": self.stats["structures_not_found"]
        }
    
    def initialize_cache_from_existing_files(self):
        """
        Initialize cache with existing PDB structure files to avoid re-downloading
        This scans the shared data/proteins_3d_structure/ directory
        """
        shared_structures_dir = Path("data/proteins_3d_structure")
        if not shared_structures_dir.exists():
            logging.info("No existing 3D structures directory found, starting fresh")
            return 0
        
        logging.info("Scanning existing 3D structure files to populate cache...")
        cached_count = 0
        
        for gene_dir in shared_structures_dir.iterdir():
            if not gene_dir.is_dir() or gene_dir.name.startswith('.'):
                continue
                
            gene_name = gene_dir.name
            
            # Scan for PDB files in this gene directory
            for pdb_file in gene_dir.glob("*.pdb"):
                if pdb_file.stat().st_size > 0:  # Only cache non-empty files
                    # Extract species name from filename (convert back from safe filename)
                    species_name = pdb_file.stem.replace('_', ' ')
                    
                    # Check if already in cache
                    cached_result = self.cache.get_cached_result(gene_name, species_name)
                    if cached_result is None:
                        # Try to extract PDB ID from file content for more accurate caching
                        pdb_id = None
                        try:
                            with open(pdb_file, 'r') as f:
                                first_line = f.readline().strip()
                                if first_line.startswith('HEADER'):
                                    # PDB files often have the ID in the header
                                    parts = first_line.split()
                                    if len(parts) >= 4:
                                        pdb_id = parts[-1].lower()  # Usually the last part
                        except:
                            pass
                        
                        if not pdb_id:
                            pdb_id = f"existing_{pdb_file.stem}"
                        
                        # Create cache entry for this existing file
                        cache_entry = {
                            "found": True,
                            "pdb_ids": [pdb_id],
                            "uniprot_acc": "unknown",
                            "found_via_gene": gene_name,
                            "cached_from_existing": True,
                            "file_path": str(pdb_file)
                        }
                        
                        self.cache.cache_result(gene_name, species_name, cache_entry)
                        cached_count += 1
                        logging.debug(f"  Cached existing: {gene_name} - {species_name} ({pdb_file.name})")
        
        if cached_count > 0:
            self.cache.save_cache()
            logging.info(f"✓ Initialized cache with {cached_count} existing 3D structure files")
        else:
            logging.info("No existing 3D structure files found to cache")
        
        return cached_count

def load_gene_aliases():
    """Load gene aliases from the standard location"""
    aliases_file = None
    
    # Find the most recent gene_aliases.json file
    for paramset_dir in Path("data/quickgo").glob("params_*"):
        alias_file = paramset_dir / "gene_aliases.json"
        if alias_file.exists():
            aliases_file = alias_file
            break
    
    if not aliases_file:
        logging.warning("No gene aliases file found")
        return {}
    
    try:
        with open(aliases_file, 'r') as f:
            aliases = json.load(f)
        logging.info(f"Loaded aliases for {len(aliases)} genes from {aliases_file}")
        return aliases
    except Exception as e:
        logging.warning(f"Could not load gene aliases: {e}")
        return {}

def search_pdb_via_uniprot_cached(downloader, gene, species, gene_aliases=None):
    """
    Search for PDB structures via UniProt API with caching and alias support
    Enhanced with proper caching for alias-based results
    """
    
    # Check cache first for primary gene name
    cached_result = downloader.cache.get_cached_result(gene, species)
    if cached_result is not None:
        downloader.stats["cache_hits"] += 1
        if cached_result.get("found", False):
            downloader.stats["structures_found"] += 1
            # If this was cached from existing file, log that we're skipping download
            if cached_result.get("cached_from_existing", False):
                logging.debug(f"  ✓ Found existing structure file for {gene} - {species}, skipping download")
            return cached_result
        else:
            downloader.stats["structures_not_found"] += 1
            return None
    
    # Not in cache, perform search
    downloader.stats["cache_misses"] += 1
    
    def try_search(search_gene, attempt_num):
        """Try searching with a specific gene name"""
        try:
            downloader.stats["api_calls"] += 1
            logging.debug(f"  Attempt {attempt_num}: Searching UniProt for {search_gene} in {species}")
            
            # Search UniProt for the gene in this species
            query = f'gene:{search_gene} AND organism_name:"{species}"'
            uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                'query': query,
                'format': 'json',
                'size': 10,
                'fields': 'accession,organism_name,gene_names,xref_pdb'
            }
            
            response = requests.get(uniprot_url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            results = data.get('results', [])
            
            if not results:
                logging.debug(f"    No UniProt entries found for {search_gene} in {species}")
                return None
            
            # Look for PDB cross-references
            for entry in results:
                pdb_refs = entry.get('uniProtKBCrossReferences', [])
                pdb_structures = [ref for ref in pdb_refs if ref.get('database') == 'PDB']
                
                if pdb_structures:
                    pdb_ids = [ref['id'] for ref in pdb_structures]
                    logging.debug(f"    Found PDB structures for {search_gene}: {pdb_ids}")
                    return {
                        'found': True,
                        'pdb_ids': pdb_ids,
                        'uniprot_acc': entry.get('primaryAccession', ''),
                        'found_via_gene': search_gene
                    }
            
            logging.debug(f"    UniProt entries found for {search_gene} but no PDB structures")
            return None
            
        except Exception as e:
            logging.debug(f"    Search failed for {search_gene}: {e}")
            return None
    
    # Try primary gene name first
    result = try_search(gene, 1)
    found_via_alias = None
    
    # If primary failed and we have aliases, try them
    if not result and gene_aliases and gene in gene_aliases:
        aliases_list = gene_aliases[gene][:3]  # Limit to first 3 aliases
        for i, alias in enumerate(aliases_list, 2):
            if alias != gene:  # Don't retry the same name
                logging.debug(f"  Trying alias '{alias}' for {gene}")
                result = try_search(alias, i)
                if result:
                    found_via_alias = alias
                    logging.info(f"  Found structure via alias '{alias}' for {gene}")
                    break
    
    # Cache the result (positive or negative)
    if result:
        downloader.stats["structures_found"] += 1
        downloader.cache.cache_result(gene, species, result)
        
        # If found via alias, also cache under original gene name for future searches
        if found_via_alias:
            logging.debug(f"  ✓ Cached structure result under original gene name '{gene}' for future searches")
        
        return result
    else:
        downloader.stats["structures_not_found"] += 1
        # Cache negative result
        negative_result = {"found": False, "pdb_ids": [], "searched_aliases": gene_aliases.get(gene, []) if gene_aliases else []}
        downloader.cache.cache_result(gene, species, negative_result)
        return None

def download_pdb_structure(pdb_id, output_file):
    """Download PDB structure file"""
    try:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(pdb_url, timeout=30)
        response.raise_for_status()
        
        with open(output_file, 'w') as f:
            f.write(response.text)
        
        return True
    except Exception as e:
        logging.debug(f"    Failed to download PDB {pdb_id}: {e}")
        return False

def update_structures_metadata(gene_dir, gene_name, successful_species, failed_species):
    """
    Update the metadata file for gene structures with download results
    """
    metadata_file = gene_dir / f"{gene_name}_structures_metadata.json"
    
    # Load existing metadata if it exists
    metadata = {
        "gene": gene_name,
        "last_updated": datetime.now().isoformat(),
        "available_species": [],
        "failed_species": [],
        "total_species": 0,
        "success_rate": 0.0
    }
    
    if metadata_file.exists():
        try:
            with open(metadata_file, 'r') as f:
                existing_metadata = json.load(f)
                # Keep existing successful species and add new ones
                existing_successful = set(existing_metadata.get("available_species", []))
                existing_failed = set(existing_metadata.get("failed_species", []))
            
            # Update with new results
            all_successful = existing_successful.union(set(successful_species))
            all_failed = existing_failed.union(set(failed_species))
            
            # Remove species from failed if they're now successful
            all_failed = all_failed - all_successful
            
            metadata["available_species"] = sorted(list(all_successful))
            metadata["failed_species"] = sorted(list(all_failed))
            
        except Exception as e:
            logging.warning(f"Could not load existing metadata for {gene_name}: {e}")
            metadata["available_species"] = successful_species
            metadata["failed_species"] = failed_species
    else:
        metadata["available_species"] = successful_species
        metadata["failed_species"] = failed_species
    
    # Update statistics
    metadata["total_species"] = len(metadata["available_species"]) + len(metadata["failed_species"])
    if metadata["total_species"] > 0:
        metadata["success_rate"] = len(metadata["available_species"]) / metadata["total_species"] * 100
    
    # Save updated metadata
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logging.info(f"Updated structures metadata for {gene_name}: {len(metadata['available_species'])} available, {len(metadata['failed_species'])} failed")

def download_structures_shared(protein_list_file, analysis, paramset, group):
    """
    Download 3D structures to shared directory structure with caching
    """
    
    logging.info(f"=== Downloading 3D Structures to Shared Directory ===")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    
    # Read protein list to get target genes
    import pandas as pd
    try:
        proteins_df = pd.read_csv(protein_list_file, sep='\t')
        # Handle both 'gene' and 'protein' column names
        if 'protein' in proteins_df.columns:
            target_genes = proteins_df['protein'].unique().tolist()
        elif 'gene' in proteins_df.columns:
            target_genes = proteins_df['gene'].unique().tolist()
        else:
            raise ValueError("Protein list file must have 'protein' or 'gene' column")
        logging.info(f"Target genes: {len(target_genes)}")
    except Exception as e:
        logging.error(f"Error reading protein list: {e}")
        return
    
    # Create shared structures data directory
    shared_structures_dir = Path("data/proteins_3d_structure")
    shared_structures_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize cached downloader
    downloader = CachedStructureDownloader()
    
    # Initialize cache with existing structure files to avoid re-downloading
    existing_cached = downloader.initialize_cache_from_existing_files()
    if existing_cached > 0:
        logging.info(f"Cache initialized with {existing_cached} existing structure files - these won't be re-downloaded")
    
    # Load gene aliases
    gene_aliases = load_gene_aliases()
    
    # Load all available species from the shared protein directory
    protein_fasta_dir = Path("data/proteins_fasta")
    
    total_genes = len(target_genes)
    total_structures_downloaded = 0
    
    # Process each gene
    for gene_idx, gene_name in enumerate(target_genes, 1):
        logging.info(f"\n=== Processing Gene {gene_idx}/{total_genes}: {gene_name} ===")
        
        # Create gene directory for structures
        gene_structures_dir = shared_structures_dir / gene_name
        gene_structures_dir.mkdir(exist_ok=True)
        
        # Find available species for this gene from FASTA directory
        gene_fasta_dir = protein_fasta_dir / gene_name
        available_species = []
        
        if gene_fasta_dir.exists():
            for fasta_file in gene_fasta_dir.glob("*.fasta"):
                if not fasta_file.name.startswith("_"):  # Skip metadata files
                    # Convert filename back to species name
                    species_name = fasta_file.stem.replace('_', ' ')
                    available_species.append(species_name)
        
        if not available_species:
            logging.warning(f"No FASTA sequences found for {gene_name}, skipping 3D structure search")
            continue
        
        logging.info(f"Found {len(available_species)} species with FASTA sequences for {gene_name}")
        
        # Check which species already have structures
        existing_species = []
        for species in available_species:
            safe_species = species.replace(' ', '_').replace('/', '_')
            pdb_file = gene_structures_dir / f"{safe_species}.pdb"
            if pdb_file.exists() and pdb_file.stat().st_size > 0:
                existing_species.append(species)
        
        # Filter out existing species
        species_to_search = [s for s in available_species if s not in existing_species]
        
        logging.info(f"Already have structures: {len(existing_species)}, Need to search: {len(species_to_search)}")
        
        if not species_to_search:
            logging.info(f"All available species already have structures for {gene_name}")
            continue
        
        # Search for structures
        successful_species = list(existing_species)  # Start with existing
        failed_species = []
        
        for species_idx, species in enumerate(species_to_search, 1):
            logging.info(f"  Searching structures for {species} ({species_idx}/{len(species_to_search)})")
            
            # Search for PDB structures with caching and alias support
            pdb_info = search_pdb_via_uniprot_cached(downloader, gene_name, species, gene_aliases)
            
            if pdb_info and pdb_info['pdb_ids']:
                # Try to download the first available PDB structure
                downloaded = False
                for pdb_id in pdb_info['pdb_ids'][:3]:  # Try up to 3 structures
                    safe_species = species.replace(' ', '_').replace('/', '_')
                    pdb_file = gene_structures_dir / f"{safe_species}.pdb"
                    
                    logging.info(f"    Downloading PDB {pdb_id}...")
                    if download_pdb_structure(pdb_id, pdb_file):
                        successful_species.append(species)
                        total_structures_downloaded += 1
                        gene_used = pdb_info.get('found_via_gene', gene_name)
                        if gene_used != gene_name:
                            logging.info(f"    ✓ Downloaded {pdb_id} for {species} (found via alias '{gene_used}')")
                        else:
                            logging.info(f"    ✓ Downloaded {pdb_id} for {species}")
                        downloaded = True
                        break
                    else:
                        time.sleep(0.5)  # Brief delay between download attempts
                
                if not downloaded:
                    failed_species.append(species)
                    logging.info(f"    ✗ Failed to download any structures for {species}")
            else:
                failed_species.append(species)
                logging.info(f"    ✗ No PDB structures found for {species}")
            
            # Rate limiting
            time.sleep(0.3)
        
        # Update gene structures metadata
        update_structures_metadata(gene_structures_dir, gene_name, successful_species, failed_species)
        
        # Gene summary
        success_rate = len(successful_species) / len(available_species) * 100 if available_species else 0
        logging.info(f"Gene {gene_name} structures summary: {len(successful_species)}/{len(available_species)} species ({success_rate:.1f}%)")
    
    # Final cache save and statistics
    downloader.finalize()
    
    # Print final summary
    cache_stats = downloader.get_cache_stats()
    logging.info(f"\n=== 3D Structures Download Summary ===")
    logging.info(f"Total genes processed: {total_genes}")
    logging.info(f"Structures downloaded: {total_structures_downloaded}")
    logging.info(f"Cache hit rate: {cache_stats['cache_hit_rate']:.1f}%")
    logging.info(f"API calls made: {cache_stats['api_calls']}")
    logging.info(f"Structures found: {cache_stats['structures_found']}")
    logging.info(f"Structures not found: {cache_stats['structures_not_found']}")
    logging.info(f"Shared structures directory: {shared_structures_dir}")

def main():
    """Main function for Snakemake integration"""
    logging.info("=== Cached Shared 3D Structures Downloader ===")
    
    try:
        # Get inputs from Snakemake
        protein_list_file = snakemake.input.protein_list
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        logging.info(f"Protein list file: {protein_list_file}")
        
        # Download structures to shared directory
        download_structures_shared(protein_list_file, analysis, paramset, group)
        
        # Create sentinel file to indicate completion
        sentinel_file = Path(snakemake.output.sentinel)
        sentinel_file.parent.mkdir(parents=True, exist_ok=True)
        sentinel_file.touch()
        
        logging.info(f"3D structures download completed. Sentinel file created: {sentinel_file}")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        print("Usage: This script should be run from Snakemake")
    except Exception as e:
        logging.error(f"Error during 3D structures download: {e}")
        raise

if __name__ == "__main__":
    main()