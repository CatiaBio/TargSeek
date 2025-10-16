#!/usr/bin/env python3
"""
Download proteins using NCBI directly with cross-analysis copying
================================================================

This script downloads proteins using NCBI Entrez directly, with support for:
- Cross-analysis copying (check previous analyses before downloading)
- "spp." entries (search genus + gene combinations)
- Analysis-specific caching and directory structure
- Proper handling of species names from BacDive gram classification
"""

import sys
import json
import logging
import time
from pathlib import Path
from datetime import datetime
from urllib.error import HTTPError
import random

# Add scripts directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

from Bio import Entrez, SeqIO
from io import StringIO

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class NCBIProteinDownloader:
    """NCBI-based protein downloader with caching and cross-analysis support"""

    def __init__(self, ncbi_credentials, analysis=None, paramset=None):
        self.analysis = analysis
        self.paramset = paramset
        self.analysis_paramset = f"{analysis}_{paramset}" if analysis and paramset else None

        # Set up NCBI credentials
        Entrez.email = ncbi_credentials['email']
        if ncbi_credentials.get('api_key'):
            Entrez.api_key = ncbi_credentials['api_key']

        # Cache for this session
        self.session_cache = {}

        # Analysis-specific cache file
        if self.analysis_paramset:
            cache_dir = Path(f"cache/protein_sequences/{self.analysis_paramset}")
            cache_dir.mkdir(parents=True, exist_ok=True)
            self.cache_file = cache_dir / "ncbi_protein_cache.json"
        else:
            self.cache_file = Path("cache/protein_sequences/ncbi_protein_cache.json")

        # Load existing cache
        self.persistent_cache = self._load_cache()

    def _load_cache(self):
        """Load persistent cache from file"""
        try:
            if self.cache_file.exists():
                with open(self.cache_file, 'r') as f:
                    return json.load(f)
        except Exception as e:
            logging.warning(f"Could not load cache: {e}")
        return {}

    def _save_cache(self):
        """Save persistent cache to file"""
        try:
            self.cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(self.cache_file, 'w') as f:
                json.dump(self.persistent_cache, f, indent=2)
        except Exception as e:
            logging.error(f"Could not save cache: {e}")

    def _cache_key(self, gene, species):
        """Generate cache key for gene-species combination"""
        return f"{gene}:{species.lower()}:ncbi"

    def search_ncbi_protein(self, gene, species, max_retries=3):
        """Search NCBI protein database for gene in species"""

        cache_key = self._cache_key(gene, species)

        # Check session cache first
        if cache_key in self.session_cache:
            return self.session_cache[cache_key]

        # Check persistent cache
        if cache_key in self.persistent_cache:
            cached_data = self.persistent_cache[cache_key]
            self.session_cache[cache_key] = cached_data
            return cached_data

        # Perform NCBI search
        for attempt in range(max_retries):
            try:
                # Create search query
                if "spp." in species:
                    # For "spp." entries, search by genus
                    genus = species.replace(" spp.", "").strip()
                    query = f'({gene}[Gene Name] OR {gene}[All Fields]) AND {genus}[Organism]'
                else:
                    # Exact species search
                    query = f'({gene}[Gene Name] OR {gene}[All Fields]) AND "{species}"[Organism]'

                logging.debug(f"NCBI query: {query}")

                # Search for protein IDs
                handle = Entrez.esearch(db="protein", term=query, retmax=5)
                search_results = Entrez.read(handle)
                handle.close()

                if search_results["IdList"]:
                    # Get the first protein sequence
                    protein_id = search_results["IdList"][0]

                    # Fetch sequence
                    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                    fasta_data = handle.read()
                    handle.close()

                    # Parse FASTA
                    if fasta_data.strip():
                        fasta_io = StringIO(fasta_data)
                        records = list(SeqIO.parse(fasta_io, "fasta"))
                        if records:
                            record = records[0]

                            # Extract organism from description if possible
                            description_parts = record.description.split('[')
                            organism = species
                            if len(description_parts) > 1:
                                org_part = description_parts[-1].replace(']', '')
                                if org_part and len(org_part) > 3:
                                    organism = org_part

                            result = {
                                'found': True,
                                'accession': record.id,
                                'sequence': str(record.seq),
                                'description': record.description,
                                'organism': organism,
                                'source': 'ncbi',
                                'cached_at': datetime.now().isoformat(),
                                'gene': gene,
                                'species': species
                            }

                            # Cache result
                            self.session_cache[cache_key] = result
                            self.persistent_cache[cache_key] = result

                            return result

                # No results found
                result = {
                    'found': False,
                    'cached_at': datetime.now().isoformat(),
                    'gene': gene,
                    'species': species,
                    'source': 'ncbi'
                }

                # Cache negative result
                self.session_cache[cache_key] = result
                self.persistent_cache[cache_key] = result

                return result

            except HTTPError as e:
                if e.code == 429:
                    wait_time = (attempt + 1) * 2 + random.uniform(1, 3)
                    logging.warning(f"Rate limited, waiting {wait_time:.1f}s before retry {attempt + 1}/{max_retries}")
                    time.sleep(wait_time)
                    continue
                else:
                    logging.error(f"HTTP Error {e.code} for {gene} in {species}: {e}")
                    break
            except Exception as e:
                logging.error(f"Error searching {gene} in {species}: {e}")
                break

        # Return failure result
        result = {
            'found': False,
            'cached_at': datetime.now().isoformat(),
            'gene': gene,
            'species': species,
            'source': 'ncbi',
            'error': 'Search failed after retries'
        }

        self.session_cache[cache_key] = result
        self.persistent_cache[cache_key] = result
        return result

    def finalize(self):
        """Save cache and cleanup"""
        self._save_cache()
        logging.info(f"Cache saved to {self.cache_file}")

def find_existing_sequences_in_other_analyses(gene_name, species, base_protein_dir, current_analysis_paramset):
    """
    Check if sequence exists in other analysis directories
    Handles both exact species matches and genus matches for "spp." entries
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
                                logging.info(f"Found {genus} species {fasta_file.name} in {analysis_dir.name} for copying")
                                return fasta_file
    else:
        # Exact species match
        safe_species = species.replace(' ', '_').replace('/', '_').replace(':', '_')
        filename = f"{safe_species}.fasta"

        # Check all analysis directories except current one
        for analysis_dir in base_protein_dir.iterdir():
            if analysis_dir.is_dir() and analysis_dir.name != current_analysis_paramset:
                gene_dir = analysis_dir / gene_name
                sequence_file = gene_dir / filename
                if sequence_file.exists() and sequence_file.stat().st_size > 0:
                    logging.info(f"Found existing sequence for {species} in {analysis_dir.name} for copying")
                    return sequence_file

    return None

def download_proteins_with_ncbi(gene_species_mapping, gene_aliases, ncbi_credentials, analysis, paramset, group):
    """
    Download proteins using NCBI with cross-analysis copying and proper caching
    """

    analysis_paramset = f"{analysis}_{paramset}"
    logging.info(f"=== NCBI Protein Download with Cross-Analysis Copying ===")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    logging.info(f"Analysis-Paramset: {analysis_paramset}")
    logging.info(f"Genes to process: {len(gene_species_mapping)}")

    # Create analysis-specific protein data directory
    sentinel_path = Path(snakemake.output.sentinel)
    base_protein_dir = sentinel_path.parent.parent  # data/protein_sequences
    analysis_protein_dir = base_protein_dir / analysis_paramset
    analysis_protein_dir.mkdir(parents=True, exist_ok=True)

    # Initialize NCBI downloader
    downloader = NCBIProteinDownloader(ncbi_credentials, analysis, paramset)

    logging.info(f"Analysis-specific protein directory: {analysis_protein_dir}")
    logging.info(f"Analysis-specific cache file: {downloader.cache_file}")

    total_genes = len(gene_species_mapping)
    total_successful_downloads = 0
    total_species_attempted = 0
    total_copied_from_other_analyses = 0
    total_cache_hits = 0
    total_new_lookups = 0

    # Process each gene
    for gene_idx, (gene_name, species_list) in enumerate(gene_species_mapping.items(), 1):
        # Create gene directory in analysis-specific location
        gene_dir = analysis_protein_dir / gene_name
        gene_dir.mkdir(exist_ok=True)

        logging.info(f"\n=== Processing Gene {gene_idx}/{total_genes}: {gene_name} ===")
        logging.info(f"Species to process: {len(species_list)}")

        successful_species = []
        failed_species = []

        # Process each species for this gene
        for species_idx, species in enumerate(species_list, 1):
            total_species_attempted += 1

            # Create safe filename
            safe_species = species.replace(' ', '_').replace('/', '_').replace(':', '_')
            fasta_file = gene_dir / f"{safe_species}.fasta"

            # Check if file already exists in current analysis
            if fasta_file.exists() and fasta_file.stat().st_size > 0:
                successful_species.append(species)
                logging.info(f"  ✓ {species} (already exists)")
                continue

            # Check if sequence exists in other analyses
            existing_file = find_existing_sequences_in_other_analyses(
                gene_name, species, base_protein_dir, analysis_paramset
            )

            if existing_file:
                # Copy from other analysis
                import shutil
                shutil.copy2(existing_file, fasta_file)
                successful_species.append(species)
                total_copied_from_other_analyses += 1
                logging.info(f"  ✓ {species} (copied from {existing_file.parent.parent.name})")
                continue

            # Need to download from NCBI
            logging.info(f"  Downloading {species} ({species_idx}/{len(species_list)}) from NCBI...")

            # Search NCBI
            result = downloader.search_ncbi_protein(gene_name, species)

            if result.get('found'):
                # Save FASTA file
                header = f">{result['accession']}|{gene_name}|{result['organism']}|{result.get('description', 'Unknown protein')}"
                sequence = result['sequence']

                with open(fasta_file, 'w') as f:
                    f.write(f"{header}\n")
                    # Write sequence in 60-character lines
                    for i in range(0, len(sequence), 60):
                        f.write(f"{sequence[i:i+60]}\n")

                successful_species.append(species)
                total_successful_downloads += 1
                total_new_lookups += 1
                logging.info(f"    ✓ Downloaded {species} ({result['accession']}, {len(sequence)} aa)")
            else:
                failed_species.append(species)
                total_new_lookups += 1
                logging.info(f"    ✗ Failed to download {species}")

        # Gene summary
        success_rate = len(successful_species) / len(species_list) * 100 if species_list else 0
        logging.info(f"Gene {gene_name} summary: {len(successful_species)}/{len(species_list)} species ({success_rate:.1f}%)")

    # Final cleanup and statistics
    downloader.finalize()

    # Print final summary
    logging.info(f"\n=== Download Summary ===")
    logging.info(f"Analysis-Paramset: {analysis_paramset}")
    logging.info(f"Total genes processed: {total_genes}")
    logging.info(f"Total species attempted: {total_species_attempted}")
    logging.info(f"Successful downloads from NCBI: {total_successful_downloads}")
    logging.info(f"Copied from other analyses: {total_copied_from_other_analyses}")
    logging.info(f"New NCBI lookups: {total_new_lookups}")
    logging.info(f"Analysis-specific protein directory: {analysis_protein_dir}")
    logging.info(f"Analysis-specific cache file: {downloader.cache_file}")

def load_gene_species_lists(gene_lists_dir):
    """Load all gene-specific species lists from the analysis directory"""
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

def main():
    """Main function for Snakemake integration"""
    logging.info("=== NCBI Protein Downloader with Cross-Analysis Copying ===")

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

        # Load gene aliases (not used in NCBI version but kept for compatibility)
        gene_aliases = {}

        # Load gene-species mappings
        gene_species_mapping = load_gene_species_lists(gene_lists_dir)

        if not gene_species_mapping:
            logging.error("No gene-species mappings found!")
            return

        # Download proteins using NCBI
        download_proteins_with_ncbi(gene_species_mapping, gene_aliases, ncbi_credentials, analysis, paramset, group)

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