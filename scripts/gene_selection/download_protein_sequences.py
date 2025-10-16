#!/usr/bin/env python3
"""
Download proteins to shared data directory with caching
======================================================

This script downloads proteins to a shared data/protein_sequences/ directory
organized by gene, regardless of Gram classification. Uses caching to avoid
re-downloading the same sequences.

Directory structure:
- data/protein_sequences/{gene}/{species}.fasta
- data/protein_sequences/{gene}/{gene}_metadata.json (tracks available species)

Cache structure:
- cache/protein_sequences/protein_fasta_cache.json
"""

import sys
import json
import logging
from pathlib import Path
from datetime import datetime

# Add scripts directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

from protein_cache import CachedProteinDownloader
from Bio import Entrez
import time
from urllib.error import HTTPError

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def search_ncbi_bulk_aliases(species, gene_name, aliases, ncbi_credentials, max_query_length=2000):
    """
    Search NCBI with bulk query for all gene aliases at once
    
    Args:
        species: Species name to search for
        gene_name: Primary gene name 
        aliases: List of alias names to search
        ncbi_credentials: Dict with email and api_key
        max_query_length: Maximum query length to avoid URL limits
    
    Returns:
        bool: True if any alias found sequences, False otherwise
    """
    if not aliases:
        return False
    
    # Set up Entrez
    Entrez.email = ncbi_credentials['email']
    if ncbi_credentials.get('api_key'):
        Entrez.api_key = ncbi_credentials['api_key']
    
    try:
        # Create bulk query with all aliases
        unique_aliases = list(set([alias for alias in aliases if alias != gene_name]))
        alias_terms = [f'("{alias}"[Gene Name] OR "{alias}"[All Fields])' for alias in unique_aliases]
        bulk_query = f'({" OR ".join(alias_terms)}) AND "{species}"[Organism]'
        
        # Check query length
        if len(bulk_query) > max_query_length:
            logging.warning(f"Bulk query too long ({len(bulk_query)} chars), falling back to individual searches")
            return False
        
        # Perform bulk search
        logging.debug(f"    Bulk NCBI search: {len(unique_aliases)} aliases for '{species}'")
        handle = Entrez.esearch(db="protein", term=bulk_query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        
        count = int(record["Count"])
        if count > 0:
            logging.info(f"    ✓ Bulk search found {count} sequences for '{species}' with aliases")
            return True
        else:
            logging.info(f"    ✗ Bulk search found no sequences for '{species}' with aliases")
            return False
            
    except HTTPError as e:
        if e.code == 429:
            time.sleep(2)  # Rate limit backoff
            logging.warning(f"Rate limited, retrying bulk search for '{species}'")
            return search_ncbi_bulk_aliases(species, gene_name, aliases, ncbi_credentials, max_query_length)
        else:
            logging.error(f"HTTP Error {e.code} in bulk search for '{species}': {e}")
            return False
    except Exception as e:
        logging.error(f"Error in bulk alias search for '{species}': {e}")
        return False

def load_gene_species_lists(gene_lists_dir):
    """
    Load all gene-specific species lists from the analysis directory
    
    Returns:
        dict: {gene_name: [species_list]}
    """
    gene_lists_path = Path(gene_lists_dir)
    gene_species_mapping = {}
    
    if not gene_lists_path.exists():
        logging.error(f"Gene lists directory not found: {gene_lists_dir}")
        return {}
    
    # Find all .txt files (gene species lists)
    gene_files = list(gene_lists_path.glob("*.txt"))
    logging.info(f"Found {len(gene_files)} gene files to process")
    
    for gene_file in gene_files:
        if gene_file.name.startswith("_"):
            continue  # Skip summary files
            
        gene_name = gene_file.stem
        
        # Load species list for this gene
        with open(gene_file, 'r') as f:
            species_list = [line.strip() for line in f if line.strip()]
        
        gene_species_mapping[gene_name] = species_list
        logging.info(f"Loaded {len(species_list)} species for gene '{gene_name}'")
    
    return gene_species_mapping

def update_gene_metadata(gene_dir, gene_name, successful_species, failed_species):
    """
    Update the metadata file for a gene with download results
    
    Args:
        gene_dir: Path to gene directory
        gene_name: Name of the gene
        successful_species: List of species successfully downloaded
        failed_species: List of species that failed to download
    """
    metadata_file = gene_dir / f"{gene_name}_metadata.json"
    
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
    
    logging.info(f"Updated metadata for {gene_name}: {len(metadata['available_species'])} available, {len(metadata['failed_species'])} failed")

def find_existing_sequences_in_other_analyses(gene_name, species, base_protein_dir, current_analysis_paramset):
    """
    Check if sequence exists in other analysis directories

    Args:
        gene_name: Name of the gene
        species: Species name (may include "spp." entries)
        base_protein_dir: Base protein sequences directory (data/protein_sequences)
        current_analysis_paramset: Current analysis_paramset to exclude

    Returns:
        Path or None: Path to existing sequence file if found
    """

    # Handle "spp." entries by looking for any species in that genus
    if "spp." in species:
        genus = species.replace(" spp.", "").strip()
        logging.debug(f"Looking for any {genus} species in previous analyses")

        # Check all analysis directories except current one
        for analysis_dir in base_protein_dir.iterdir():
            if analysis_dir.is_dir() and analysis_dir.name != current_analysis_paramset:
                gene_dir = analysis_dir / gene_name
                if gene_dir.exists():
                    # Look for any FASTA file that starts with the genus
                    for fasta_file in gene_dir.glob("*.fasta"):
                        # Check if filename starts with genus
                        if fasta_file.stem.lower().startswith(genus.lower().replace(' ', '_')):
                            if fasta_file.stat().st_size > 0:
                                logging.debug(f"Found {genus} species {fasta_file.name} in {analysis_dir.name}")
                                return fasta_file
    else:
        # Exact species match
        safe_species = species.replace(' ', '_').replace('/', '_')
        filename = f"{safe_species}.fasta"

        # Check all analysis directories except current one
        for analysis_dir in base_protein_dir.iterdir():
            if analysis_dir.is_dir() and analysis_dir.name != current_analysis_paramset:
                gene_dir = analysis_dir / gene_name
                sequence_file = gene_dir / filename
                if sequence_file.exists() and sequence_file.stat().st_size > 0:
                    logging.debug(f"Found existing sequence for {species} in {analysis_dir.name}")
                    return sequence_file

    return None

def download_proteins_shared(gene_species_mapping, gene_aliases, ncbi_credentials, analysis, paramset, group):
    """
    Download proteins to analysis-specific directory structure with cross-analysis copying and caching

    Args:
        gene_species_mapping: Dict of {gene: [species_list]}
        gene_aliases: Dict of {gene: [alias_list]} for fallback searches
        ncbi_credentials: Dict with NCBI credentials for bulk alias searches
        analysis: Analysis name for directory organization
        paramset: Parameter set name for directory organization
        group: Gram group name for logging
    """

    analysis_paramset = f"{analysis}_{paramset}"
    logging.info(f"=== Downloading Proteins to Analysis-Specific Directory ===")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    logging.info(f"Analysis-Paramset: {analysis_paramset}")
    logging.info(f"Genes to process: {len(gene_species_mapping)}")

    # Create analysis-specific protein data directory
    sentinel_path = Path(snakemake.output.sentinel)
    base_protein_dir = sentinel_path.parent.parent  # data/protein_sequences
    analysis_protein_dir = base_protein_dir / analysis_paramset
    analysis_protein_dir.mkdir(parents=True, exist_ok=True)

    # Initialize cached downloader with analysis-specific cache
    downloader = CachedProteinDownloader(analysis=analysis, paramset=paramset)

    # Show cache status once
    cache_size = len(downloader.cache.cache)
    logging.info(f"Initialized analysis-specific protein cache with {cache_size} entries")
    logging.info(f"Cache file: {downloader.cache.cache_file}")
    logging.info(f"Protein directory: {analysis_protein_dir}")
    
    total_genes = len(gene_species_mapping)
    total_successful_downloads = 0
    total_species_attempted = 0

    # Track overall cache performance and cross-analysis copying
    total_cache_hits = 0
    total_new_lookups = 0
    total_copied_from_other_analyses = 0

    # Process each gene
    for gene_idx, (gene_name, species_list) in enumerate(gene_species_mapping.items(), 1):
        # Create gene directory in analysis-specific location
        gene_dir = analysis_protein_dir / gene_name
        gene_dir.mkdir(exist_ok=True)

        # Check which species already exist in current analysis
        existing_species = []
        species_to_process = []

        for species in species_list:
            safe_species = species.replace(' ', '_').replace('/', '_')
            fasta_file = gene_dir / f"{safe_species}.fasta"

            if fasta_file.exists() and fasta_file.stat().st_size > 0:
                existing_species.append(species)
            else:
                # Check if sequence exists in other analyses
                existing_file = find_existing_sequences_in_other_analyses(
                    gene_name, species, base_protein_dir, analysis_paramset
                )

                if existing_file:
                    # Copy from other analysis
                    import shutil
                    shutil.copy2(existing_file, fasta_file)
                    existing_species.append(species)
                    total_copied_from_other_analyses += 1
                    logging.debug(f"Copied {species} from {existing_file.parent.parent.name}")
                else:
                    species_to_process.append(species)

        # Filter species that need downloading or cross-analysis lookup
        species_to_download = species_to_process
        
        if species_to_download:
            logging.info(f"\n=== Processing Gene {gene_idx}/{total_genes}: {gene_name} ===")
            logging.info(f"Species to download: {len(species_list)} (Already have: {len(existing_species)}, Need: {len(species_to_download)})")
        elif gene_idx % 10 == 0:
            # Progress update every 10 genes when all species exist
            logging.info(f"Progress: {gene_idx}/{total_genes} genes checked")
            
        if not species_to_download:
            continue
        
        # Download missing species
        successful_species = list(existing_species)  # Start with existing
        failed_species = []
        
        cache_hits_gene = 0
        new_lookups_gene = 0
        
        for species_idx, species in enumerate(species_to_download, 1):
            total_species_attempted += 1
            
            # Check if this will be a cache hit (and whether it should be skipped)
            cached_result = downloader.cache.get_cached_sequence(gene_name, species, "uniprot")
            
            if cached_result is not None:
                # We have a cache entry
                cache_hits_gene += 1
                total_cache_hits += 1
                is_cached = True
                
                # If it's a negative result (failure) that hasn't expired, skip this species entirely
                if not cached_result.get("found", False):
                    failed_species.append(species)
                    continue  # Skip to next species
                else:
                    # Use the cached successful result directly
                    seq_data = cached_result
                    
            else:
                # No cache entry, need to make new lookup
                new_lookups_gene += 1
                total_new_lookups += 1
                is_cached = False
                logging.info(f"  Downloading {species} ({species_idx}/{len(species_to_download)}) - NEW LOOKUP")
                
                # Handle "spp." entries by searching for representative species
                if "spp." in species:
                    genus = species.replace(" spp.", "").strip()
                    logging.info(f"    Searching for representative {genus} species...")

                    # Try common representative species for this genus
                    representative_species = [
                        f"{genus} sp.",
                        f"{genus} species",
                    ]

                    # Add genus-specific representatives
                    genus_representatives = {
                        "Acinetobacter": ["Acinetobacter baumannii", "Acinetobacter pittii"],
                        "Klebsiella": ["Klebsiella pneumoniae", "Klebsiella oxytoca"],
                        "Pseudomonas": ["Pseudomonas aeruginosa", "Pseudomonas fluorescens"],
                        "Serratia": ["Serratia marcescens"],
                        "Lactobacillus": ["Lactobacillus plantarum", "Lactobacillus acidophilus"],
                        "Microbacterium": ["Microbacterium lacticum"],
                        "Arthrobacter": ["Arthrobacter sp."],
                    }

                    if genus in genus_representatives:
                        representative_species.extend(genus_representatives[genus])

                    # Try each representative species
                    for rep_species in representative_species:
                        seq_data = downloader.search_uniprot_cached(gene_name, rep_species)
                        if seq_data:
                            logging.info(f"    ✓ Found via representative species: {rep_species}")
                            break
                else:
                    # Search for exact species name with caching
                    seq_data = downloader.search_uniprot_cached(gene_name, species)

                # If found, also update cache in other analysis caches to maintain cross-analysis consistency
                if seq_data:
                    # Update caches in other analyses for future cross-analysis copying
                    for other_analysis_dir in base_protein_dir.iterdir():
                        if other_analysis_dir.is_dir() and other_analysis_dir.name != analysis_paramset:
                            other_cache_file = f"cache/protein_sequences/{other_analysis_dir.name}/protein_fasta_cache.json"
                            if Path(other_cache_file).exists():
                                # Load other cache and update with this sequence
                                try:
                                    other_cache = CachedProteinDownloader(analysis=other_analysis_dir.name.split('_')[0],
                                                                        paramset=other_analysis_dir.name.split('_')[1])
                                    other_cache.cache.cache_sequence(gene_name, species, seq_data, "uniprot")
                                    other_cache.finalize()
                                except Exception as e:
                                    logging.debug(f"Could not update cache for {other_analysis_dir.name}: {e}")

            found_via_alias = None

            # If primary search failed for new lookups, try aliases with bulk NCBI check first
            if not seq_data and not is_cached and gene_name in gene_aliases:
                aliases = gene_aliases[gene_name]
                unique_aliases = [alias for alias in aliases if alias != gene_name]
                
                if unique_aliases:
                    logging.info(f"    Primary '{gene_name}' failed for '{species}', checking {len(unique_aliases)} aliases...")
                    
                    # First, do bulk NCBI search to see if any aliases have sequences for this species
                    bulk_found = search_ncbi_bulk_aliases(species, gene_name, aliases, ncbi_credentials)
                    
                    if bulk_found:
                        # Bulk search found sequences, now try each alias to get the actual sequence
                        for alias_name in unique_aliases:
                            seq_data = downloader.search_uniprot_cached(alias_name, species)
                            if seq_data:
                                found_via_alias = alias_name
                                logging.info(f"      ✓ Downloaded '{species}' via alias '{alias_name}'")
                                
                                # Cache this result with the ORIGINAL gene name for future searches
                                downloader.cache.cache_sequence(gene_name, species, seq_data, "uniprot")
                                break
                        
                        if not seq_data:
                            logging.warning(f"      Alias downloads failed for '{species}'")
                    else:
                        logging.debug(f"    No alias sequences available for '{species}'")
            
            if seq_data:
                # Save FASTA file
                safe_species = species.replace(' ', '_').replace('/', '_')
                fasta_file = gene_dir / f"{safe_species}.fasta"
                
                # Create FASTA content with alias info if used
                gene_for_header = found_via_alias if found_via_alias else gene_name
                alias_info = f" (found via alias: {found_via_alias})" if found_via_alias else ""
                header = f">{seq_data['accession']}|{gene_for_header}|{species}|{seq_data.get('protein_name', 'Unknown protein')}{alias_info}"
                sequence = seq_data['sequence']
                
                with open(fasta_file, 'w') as f:
                    f.write(f"{header}\n")
                    # Write sequence in 60-character lines
                    for i in range(0, len(sequence), 60):
                        f.write(f"{sequence[i:i+60]}\n")
                
                successful_species.append(species)
                total_successful_downloads += 1
                alias_msg = f" via alias '{found_via_alias}'" if found_via_alias else ""
                if not is_cached:
                    logging.info(f"    ✓ Downloaded {species} ({seq_data['accession']}, {len(sequence)} aa{alias_msg})")
            else:
                failed_species.append(species)
                alias_info = f" (tried {len(gene_aliases.get(gene_name, []))} aliases)" if gene_name in gene_aliases else ""
                if not is_cached:
                    logging.info(f"    ✗ Failed to download {species}{alias_info}")
            
            # Progress update every 20 species when using cache
            if is_cached and species_idx % 20 == 0:
                logging.info(f"  Progress: {species_idx}/{len(species_to_download)} species ({cache_hits_gene} from cache, {new_lookups_gene} new lookups)")
        
        # Update gene metadata
        update_gene_metadata(gene_dir, gene_name, successful_species, failed_species)
        
        # Gene summary
        success_rate = len(successful_species) / len(species_list) * 100 if species_list else 0
        logging.info(f"Gene {gene_name} summary: {len(successful_species)}/{len(species_list)} species ({success_rate:.1f}%)")
    
    # Final cache save and statistics
    downloader.finalize()
    
    # Print final summary
    cache_stats = downloader.get_cache_stats()
    logging.info(f"\n=== Download Summary ===")
    logging.info(f"Analysis-Paramset: {analysis_paramset}")
    logging.info(f"Total genes processed: {total_genes}")
    logging.info(f"Total species attempted: {total_species_attempted}")
    logging.info(f"Successful downloads: {total_successful_downloads}")
    logging.info(f"Copied from other analyses: {total_copied_from_other_analyses}")
    logging.info(f"Cache performance:")
    logging.info(f"  Cache hits: {total_cache_hits}")
    logging.info(f"  New lookups: {total_new_lookups}")
    logging.info(f"  Cache hit rate: {(total_cache_hits / (total_cache_hits + total_new_lookups) * 100) if (total_cache_hits + total_new_lookups) > 0 else 0:.1f}%")
    logging.info(f"Analysis-specific protein directory: {analysis_protein_dir}")
    logging.info(f"Analysis-specific cache file: {downloader.cache.cache_file}")

def main():
    """Main function for Snakemake integration"""
    logging.info("=== Cached Shared Protein Downloader with Alias Fallback ===")
    
    try:
        # Get inputs from Snakemake
        gene_lists_dir = snakemake.input.protein_lists
        aliases_file = snakemake.input.aliases
        ncbi_info_file = snakemake.input.ncbi_info
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        logging.info(f"Gene lists directory: {gene_lists_dir}")
        logging.info(f"Aliases file: {aliases_file}")
        
        # Load NCBI credentials
        with open(ncbi_info_file, 'r') as f:
            lines = f.readlines()
            ncbi_credentials = {
                'email': lines[0].strip(),
                'api_key': lines[1].strip() if len(lines) > 1 else None
            }
        logging.info(f"NCBI credentials loaded: {ncbi_credentials['email']}, API key: {'Yes' if ncbi_credentials['api_key'] else 'No'}")
        
        # Load gene aliases
        gene_aliases = {}
        aliases_path = Path(aliases_file)
        aliases_json_file = aliases_path.with_suffix('.json')
        
        if aliases_json_file.exists():
            with open(aliases_json_file, 'r') as f:
                all_aliases = json.load(f)
            logging.info(f"Loaded aliases for {len(all_aliases)} genes")
        else:
            logging.warning("Gene aliases JSON file not found, proceeding without alias fallback")
            all_aliases = {}
        
        # Load gene-species mappings
        gene_species_mapping = load_gene_species_lists(gene_lists_dir)
        
        # Filter aliases to only include genes in current analysis
        gene_aliases = {gene: all_aliases[gene] for gene in gene_species_mapping.keys() if gene in all_aliases}
        logging.info(f"Using aliases for {len(gene_aliases)} genes in current analysis")
        
        if not gene_species_mapping:
            logging.error("No gene-species mappings found!")
            return
        
        # Download proteins to shared directory
        download_proteins_shared(gene_species_mapping, gene_aliases, ncbi_credentials, analysis, paramset, group)
        
        # Create sentinel file to indicate completion
        sentinel_file = Path(snakemake.output.sentinel)
        sentinel_file.parent.mkdir(parents=True, exist_ok=True)
        sentinel_file.touch()
        
        logging.info(f"Download completed successfully. Sentinel file created: {sentinel_file}")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        print("Usage: This script should be run from Snakemake")
        print("For testing, modify the script to include test parameters")
    except Exception as e:
        logging.error(f"Error during download: {e}")
        raise

if __name__ == "__main__":
    main()