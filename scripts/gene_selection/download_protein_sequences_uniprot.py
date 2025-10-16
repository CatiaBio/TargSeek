#!/usr/bin/env python3
"""
Download proteins using UniProt filters with cross-analysis copying
==================================================================

This script downloads proteins using UniProt REST API with specific filters:
- For exact species: (protein_name:gene) AND (taxonomy_name:"Escherichia coli O157:H7")
- For spp. entries: (protein_name:gene) AND (taxonomy_name:"Alcaligenes")

Features:
- Cross-analysis copying (check previous analyses before downloading)
- Analysis-specific caching and directory structure
- Combines all sequences per species into single FASTA files
- Proper handling of species names from BacDive gram classification
"""

import sys
import json
import logging
import time
import requests
from pathlib import Path
from datetime import datetime

# Add scripts directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class UniProtDownloader:
    """UniProt-based protein downloader with caching and cross-analysis support"""

    def __init__(self, analysis=None, paramset=None):
        self.analysis = analysis
        self.paramset = paramset
        self.analysis_paramset = f"{analysis}_{paramset}" if analysis and paramset else None

        # Session for requests
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'TargSeek/1.0 (https://github.com/user/targseek; contact@example.com)'
        })

        # Cache for this session
        self.session_cache = {}

        # Analysis-specific cache file
        if self.analysis_paramset:
            cache_dir = Path(f"cache/protein_sequences/{self.analysis_paramset}")
            cache_dir.mkdir(parents=True, exist_ok=True)
            self.cache_file = cache_dir / "uniprot_protein_cache.json"
        else:
            self.cache_file = Path("cache/protein_sequences/uniprot_protein_cache.json")

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
        return f"{gene}:{species.lower()}:uniprot"

    def clean_species_for_search(self, species):
        """Clean species name for UniProt search"""
        # Remove 'spp.' and get genus for genus-level searches
        if "spp." in species:
            genus = species.replace(" spp.", "").strip()
            return genus, True  # genus, is_genus_search
        else:
            return species, False  # exact species, not genus search

    def search_uniprot_with_filters(self, protein_name, species, max_results=50):
        """
        Search UniProt using specific protein_name and taxonomy_name filters
        """

        cache_key = self._cache_key(protein_name, species)

        # Check session cache first
        if cache_key in self.session_cache:
            return self.session_cache[cache_key]

        # Check persistent cache
        if cache_key in self.persistent_cache:
            cached_data = self.persistent_cache[cache_key]
            self.session_cache[cache_key] = cached_data
            return cached_data

        # Clean species name for search
        search_term, is_genus = self.clean_species_for_search(species)

        # Build query based on whether it's genus or exact species
        if is_genus:
            query = f'(protein_name:{protein_name}) AND (taxonomy_name:"{search_term}")'
            logging.debug(f"Genus search: {query}")
        else:
            query = f'(protein_name:{protein_name}) AND (taxonomy_name:"{search_term}")'
            logging.debug(f"Exact species search: {query}")

        # UniProt REST API parameters
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': query,
            'format': 'json',
            'size': max_results,
            'fields': 'accession,id,protein_name,sequence,organism_name,gene_names,length,reviewed'
        }

        try:
            logging.debug(f"Making UniProt API request for {protein_name} in {species}")
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()
            results = data.get('results', [])

            logging.debug(f"Found {len(results)} results")

            # Process results
            sequences = []
            for i, entry in enumerate(results):
                try:
                    sequence_data = {
                        'accession': entry.get('primaryAccession', ''),
                        'entry_name': entry.get('uniProtkbId', ''),
                        'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown'),
                        'sequence': entry.get('sequence', {}).get('value', ''),
                        'length': entry.get('sequence', {}).get('length', 0),
                        'organism': entry.get('organism', {}).get('scientificName', search_term),
                        'gene_names': entry.get('genes', []),
                        'reviewed': entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)',
                        'search_species': species,
                        'search_term': search_term
                    }

                    if sequence_data['sequence']:  # Only include if we have a sequence
                        sequences.append(sequence_data)

                except Exception as e:
                    logging.debug(f"Error processing entry {i+1}: {e}")
                    continue

            # Cache result
            result = {
                'found': len(sequences) > 0,
                'sequences': sequences,
                'cached_at': datetime.now().isoformat(),
                'gene': protein_name,
                'species': species,
                'source': 'uniprot'
            }

            self.session_cache[cache_key] = result
            self.persistent_cache[cache_key] = result

            return result

        except requests.exceptions.RequestException as e:
            logging.error(f"HTTP request failed for {protein_name} in {species}: {e}")
        except json.JSONDecodeError as e:
            logging.error(f"Failed to parse JSON response for {protein_name} in {species}: {e}")
        except Exception as e:
            logging.error(f"Unexpected error for {protein_name} in {species}: {e}")

        # Return failure result
        result = {
            'found': False,
            'sequences': [],
            'cached_at': datetime.now().isoformat(),
            'gene': protein_name,
            'species': species,
            'source': 'uniprot',
            'error': 'Search failed'
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

def save_sequences_to_fasta(sequences, gene_name, species, output_dir):
    """Save all sequences for a species to a single FASTA file"""

    if not sequences:
        return 0

    # Create gene directory
    gene_dir = Path(output_dir) / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)

    # Clean species name for filename
    safe_species = species.replace(' ', '_').replace('/', '_').replace(':', '_')
    filename = f"{safe_species}.fasta"
    fasta_file = gene_dir / filename

    # Write all sequences to single FASTA file
    with open(fasta_file, 'w') as f:
        for seq_data in sequences:
            # Create FASTA header
            gene_names = []
            if seq_data['gene_names']:
                for gene_group in seq_data['gene_names']:
                    if 'geneName' in gene_group:
                        gene_names.append(gene_group['geneName']['value'])

            gene_name_str = ','.join(gene_names) if gene_names else gene_name

            header = f">{seq_data['accession']}|{gene_name_str}|{seq_data['organism']}|{seq_data['protein_name']}"
            sequence = seq_data['sequence']

            # Write header
            f.write(f"{header}\n")

            # Write sequence in 60-character lines
            for j in range(0, len(sequence), 60):
                f.write(f"{sequence[j:j+60]}\n")

    return len(sequences)

def download_proteins_with_uniprot(gene_species_mapping, gene_aliases, ncbi_credentials, analysis, paramset, group):
    """
    Download proteins using UniProt with cross-analysis copying and proper caching
    """

    analysis_paramset = f"{analysis}_{paramset}"
    logging.info(f"=== UniProt Protein Download with Cross-Analysis Copying ===")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    logging.info(f"Analysis-Paramset: {analysis_paramset}")
    logging.info(f"Genes to process: {len(gene_species_mapping)}")

    # Create analysis-specific protein data directory
    sentinel_path = Path(snakemake.output.sentinel)
    base_protein_dir = sentinel_path.parent.parent  # data/protein_sequences
    analysis_protein_dir = base_protein_dir / analysis_paramset
    analysis_protein_dir.mkdir(parents=True, exist_ok=True)

    # Initialize UniProt downloader
    downloader = UniProtDownloader(analysis, paramset)

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

            # Need to download from UniProt
            logging.info(f"  Downloading {species} ({species_idx}/{len(species_list)}) from UniProt...")

            # Search UniProt
            result = downloader.search_uniprot_with_filters(gene_name, species)

            if result.get('found') and result.get('sequences'):
                sequences = result['sequences']

                # Save sequences to FASTA file
                saved_count = save_sequences_to_fasta(sequences, gene_name, species, analysis_protein_dir)

                if saved_count > 0:
                    successful_species.append(species)
                    total_successful_downloads += 1
                    total_new_lookups += 1
                    logging.info(f"    ✓ Downloaded {species} ({saved_count} sequences)")
                else:
                    failed_species.append(species)
                    total_new_lookups += 1
                    logging.info(f"    ✗ Failed to save sequences for {species}")
            else:
                failed_species.append(species)
                total_new_lookups += 1
                logging.info(f"    ✗ No sequences found for {species}")

            # Small delay to be nice to UniProt API
            time.sleep(0.5)

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
    logging.info(f"Successful downloads from UniProt: {total_successful_downloads}")
    logging.info(f"Copied from other analyses: {total_copied_from_other_analyses}")
    logging.info(f"New UniProt lookups: {total_new_lookups}")
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
    logging.info("=== UniProt Protein Downloader with Cross-Analysis Copying ===")

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

        # NCBI credentials not used for UniProt but kept for compatibility
        ncbi_credentials = {}

        # Load gene aliases (not used in UniProt version but kept for compatibility)
        gene_aliases = {}

        # Load gene-species mappings
        gene_species_mapping = load_gene_species_lists(gene_lists_dir)

        if not gene_species_mapping:
            logging.error("No gene-species mappings found!")
            return

        # Download proteins using UniProt
        download_proteins_with_uniprot(gene_species_mapping, gene_aliases, ncbi_credentials, analysis, paramset, group)

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