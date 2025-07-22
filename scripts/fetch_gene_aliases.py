#!/usr/bin/env python3
"""
Fetch Gene Aliases from NCBI
=============================

This script fetches gene aliases/synonyms from NCBI for a list of gene symbols,
helping to identify alternative names like bamA/yaeT that refer to the same protein.
"""

import requests
import time
import pandas as pd
from pathlib import Path
import logging
import json
from Bio import Entrez
import xml.etree.ElementTree as ET
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Cache configuration
CACHE_DIR = Path("cache/gene_aliases")
CACHE_FILE = CACHE_DIR / "gene_aliases_cache.json"

class GeneAliasCache:
    """Manages caching of gene alias lookups"""
    
    def __init__(self, max_age_days=180):
        """Initialize cache with 6-month expiration by default"""
        self.cache = {}
        self.cache_updates = 0
        self.max_age_days = max_age_days
        self.load_cache()
    
    def load_cache(self):
        """Load existing gene alias cache"""
        if CACHE_FILE.exists():
            try:
                with open(CACHE_FILE, 'r') as f:
                    self.cache = json.load(f)
                logging.info(f"✓ Loaded gene alias cache with {len(self.cache)} entries")
            except Exception as e:
                logging.warning(f"Could not load gene alias cache: {e}")
                self.cache = {}
        else:
            logging.info("No existing gene alias cache found, starting fresh")
            self.cache = {}
    
    def save_cache(self):
        """Save gene alias cache"""
        try:
            CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with open(CACHE_FILE, 'w') as f:
                json.dump(self.cache, f, indent=2)
            logging.info(f"✓ Saved gene alias cache with {len(self.cache)} entries")
        except Exception as e:
            logging.warning(f"Could not save gene alias cache: {e}")
    
    def is_cache_valid(self, cache_entry):
        """Check if cache entry is still valid"""
        if "timestamp" not in cache_entry:
            return False
        
        try:
            cache_time = datetime.fromisoformat(cache_entry["timestamp"])
            age_days = (datetime.now() - cache_time).days
            return age_days <= self.max_age_days
        except (ValueError, TypeError):
            return False
    
    def get_cache_key(self, gene_symbol, source="ncbi"):
        """Generate cache key for gene alias lookup"""
        return f"{gene_symbol}||{source}"
    
    def get_cached_aliases(self, gene_symbol, source="ncbi"):
        """Get cached aliases if available and valid"""
        cache_key = self.get_cache_key(gene_symbol, source)
        cache_entry = self.cache.get(cache_key)
        
        if cache_entry and self.is_cache_valid(cache_entry):
            return cache_entry.get("aliases", [])
        return None
    
    def cache_aliases(self, gene_symbol, aliases, source="ncbi"):
        """Cache gene aliases"""
        cache_key = self.get_cache_key(gene_symbol, source)
        self.cache[cache_key] = {
            "aliases": aliases,
            "timestamp": datetime.now().isoformat(),
            "source": source
        }
        self.cache_updates += 1
        
        # Save cache every 20 updates
        if self.cache_updates % 20 == 0:
            self.save_cache()

def setup_ncbi_credentials():
    """Set up NCBI credentials from config file"""
    try:
        # Try to read NCBI info from config
        ncbi_info_file = Path("config/login/ncbi_info.txt")
        if ncbi_info_file.exists():
            with open(ncbi_info_file, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]
                
                # Handle different formats
                if len(lines) >= 2:
                    # Format: email on line 1, api_key on line 2
                    Entrez.email = lines[0]
                    Entrez.api_key = lines[1]
                    logging.info(f"Loaded NCBI credentials for: {Entrez.email}")
                else:
                    # Try key:value format
                    for line in lines:
                        if ':' in line:
                            key, value = line.split(':', 1)
                            key = key.strip().lower()
                            value = value.strip()
                            if key == 'email':
                                Entrez.email = value
                            elif key in ['api_key', 'apikey']:
                                Entrez.api_key = value
        
        if not hasattr(Entrez, 'email') or not Entrez.email:
            Entrez.email = "user@example.com"  # Default fallback
            logging.warning("No email found in NCBI config, using default")
        else:
            logging.info(f"Using NCBI email: {Entrez.email}")
            
        if hasattr(Entrez, 'api_key') and Entrez.api_key:
            logging.info("NCBI API key loaded successfully")
        else:
            logging.info("No NCBI API key found, using default rate limits")
            
    except Exception as e:
        logging.warning(f"Could not read NCBI credentials: {e}")
        Entrez.email = "user@example.com"

def fetch_gene_aliases_ncbi(gene_symbol, taxon_id="511145", cache=None):  # E. coli as default
    """
    Fetch gene aliases from NCBI Gene database using E-utilities with API key
    
    Args:
        gene_symbol (str): Gene symbol to search for
        taxon_id (str): NCBI taxonomy ID (default: E. coli)
        cache (GeneAliasCache): Cache instance to use
    
    Returns:
        list: List of gene aliases/synonyms
    """
    # Check cache first
    if cache:
        cached_aliases = cache.get_cached_aliases(gene_symbol, "ncbi")
        if cached_aliases is not None:
            logging.debug(f"Using cached NCBI aliases for {gene_symbol}: {cached_aliases}")
            return cached_aliases
    
    aliases = set()
    
    try:
        # Search for gene in NCBI Gene database
        search_term = f"{gene_symbol}[Gene Name] AND {taxon_id}[Taxonomy ID]"
        
        # Search for gene IDs using API
        search_results = Entrez.read(Entrez.esearch(
            db="gene", 
            term=search_term, 
            retmax=5,
            api_key=getattr(Entrez, 'api_key', None),
            email=getattr(Entrez, 'email', 'user@example.com')
        ))
        
        gene_ids = search_results["IdList"]
        
        if not gene_ids:
            logging.debug(f"No NCBI Gene entries found for {gene_symbol}")
            return list(aliases)
        
        # Fetch detailed information for found genes
        for gene_id in gene_ids[:3]:  # Limit to first 3 results
            try:
                # Use efetch with API key
                handle = Entrez.efetch(
                    db="gene", 
                    id=gene_id, 
                    rettype="xml",
                    api_key=getattr(Entrez, 'api_key', None),
                    email=getattr(Entrez, 'email', 'user@example.com')
                )
                gene_data = handle.read()
                handle.close()
                
                # Parse XML to extract aliases
                root = ET.fromstring(gene_data)
                
                # Look for gene names and aliases in various XML elements
                for elem in root.iter():
                    if elem.tag == "Gene-ref_locus":
                        if elem.text and elem.text.strip():
                            aliases.add(elem.text.strip())
                    elif elem.tag == "Gene-ref_syn":
                        for syn_elem in elem:
                            if syn_elem.text and syn_elem.text.strip():
                                aliases.add(syn_elem.text.strip())
                    elif elem.tag == "Gene-ref_desc":
                        # Sometimes aliases are mentioned in descriptions
                        if elem.text and "also known as" in elem.text.lower():
                            desc = elem.text.lower()
                            if "also known as" in desc:
                                alias_part = desc.split("also known as")[1].split(".")[0]
                                potential_aliases = [a.strip() for a in alias_part.replace("or", ",").split(",")]
                                for alias in potential_aliases:
                                    clean_alias = alias.strip().replace("(", "").replace(")", "")
                                    if clean_alias and len(clean_alias) < 20:
                                        aliases.add(clean_alias)
                
                # Rate limiting (with API key: 10 requests/second, without: 3/second)
                sleep_time = 0.1 if hasattr(Entrez, 'api_key') and Entrez.api_key else 0.34
                time.sleep(sleep_time)
                
            except Exception as e:
                logging.warning(f"Error fetching details for gene ID {gene_id}: {e}")
                continue
        
        # Remove the original gene symbol and common non-gene terms
        aliases.discard(gene_symbol.lower())
        aliases.discard(gene_symbol.upper())
        aliases.discard(gene_symbol)
        
        # Filter out non-gene-like terms
        filtered_aliases = []
        for alias in aliases:
            if alias and len(alias) >= 2 and len(alias) <= 10 and alias.replace("_", "").replace("-", "").isalnum():
                filtered_aliases.append(alias)
        
        result = sorted(filtered_aliases)
        
        # Cache the result
        if cache:
            cache.cache_aliases(gene_symbol, result, "ncbi")
        
        return result
        
    except Exception as e:
        logging.warning(f"Error fetching aliases for {gene_symbol}: {e}")
        # Cache empty result to avoid repeated failed lookups
        if cache:
            cache.cache_aliases(gene_symbol, [], "ncbi")
        return []

def fetch_gene_aliases_uniprot(gene_symbol, cache=None):
    """
    Alternative method: Fetch gene aliases from UniProt
    
    Args:
        gene_symbol (str): Gene symbol to search for
        cache (GeneAliasCache): Cache instance to use
    
    Returns:
        list: List of gene aliases/synonyms
    """
    # Check cache first
    if cache:
        cached_aliases = cache.get_cached_aliases(gene_symbol, "uniprot")
        if cached_aliases is not None:
            logging.debug(f"Using cached UniProt aliases for {gene_symbol}: {cached_aliases}")
            return cached_aliases
    
    aliases = set()
    
    try:
        # Search UniProt for the gene
        uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'gene:{gene_symbol} AND organism_id:83333',  # E. coli
            'format': 'json',
            'fields': 'gene_names,gene_synonym',
            'size': 5
        }
        
        response = requests.get(uniprot_url, params=params, timeout=30)
        if response.status_code != 200:
            logging.warning(f"UniProt search failed for {gene_symbol}: status {response.status_code}")
            return []
        
        data = response.json()
        
        for entry in data.get('results', []):
            # Extract gene names
            gene_info = entry.get('genes', [])
            for gene in gene_info:
                # Primary name
                if 'geneName' in gene:
                    primary_name = gene['geneName'].get('value', '')
                    if primary_name and primary_name != gene_symbol:
                        aliases.add(primary_name)
                
                # Synonyms
                if 'synonyms' in gene:
                    for synonym in gene['synonyms']:
                        syn_name = synonym.get('value', '')
                        if syn_name and syn_name != gene_symbol:
                            aliases.add(syn_name)
        
        result = sorted(list(aliases))
        
        # Cache the result
        if cache:
            cache.cache_aliases(gene_symbol, result, "uniprot")
        
        return result
        
    except Exception as e:
        logging.warning(f"Error fetching UniProt aliases for {gene_symbol}: {e}")
        # Cache empty result to avoid repeated failed lookups
        if cache:
            cache.cache_aliases(gene_symbol, [], "uniprot")
        return []

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        genes_file = snakemake.input.genes
        output_file = snakemake.output.aliases_file
        
        logging.info(f"Fetching gene aliases for genes in: {genes_file}")
        
    except NameError:
        # Test mode
        genes_file = "data/quickgo/params_1/gene_symbols.txt"
        output_file = "data/quickgo/params_1/gene_aliases.txt"
        logging.info("Running in test mode")
    
    # Set up NCBI credentials
    setup_ncbi_credentials()
    
    # Initialize cache
    cache = GeneAliasCache(max_age_days=180)  # 6-month cache expiration
    
    # Read and filter genes
    try:
        with open(genes_file, 'r') as f:
            all_genes = [line.strip() for line in f if line.strip()]
        
        logging.info(f"Found {len(all_genes)} total genes in input file")
        
        # Filter out database-specific identifiers
        filtered_genes = []
        excluded_genes = []
        
        for gene in all_genes:
            # Skip genes that start with a capital letter
            if gene and gene[0].isupper():
                excluded_genes.append(gene)
            else:
                filtered_genes.append(gene)
        
        genes = filtered_genes
        
        logging.info(f"After filtering:")
        logging.info(f"  Genes to process (lowercase): {len(genes)}")
        logging.info(f"  Excluded (start with capital): {len(excluded_genes)}")
        
        if excluded_genes:
            logging.info(f"  Example excluded: {excluded_genes[:5]}")
        
        if len(genes) == 0:
            logging.error("No valid gene names found after filtering!")
            return
        
    except Exception as e:
        logging.error(f"Error reading filtered genes file: {e}")
        return
    
    # Collect aliases for all genes with consolidation
    all_aliases = {}
    gene_to_primary = {}  # Map any name (gene or alias) to primary gene name
    processed = 0
    cache_hits = 0
    new_lookups = 0
    
    # Check if we're mostly using cache
    total_cache_size = len(cache.cache)
    logging.info(f"Using cache with {total_cache_size} entries")
    
    for gene in genes:
        # Check cache directly to avoid function call overhead
        ncbi_aliases = cache.get_cached_aliases(gene, "ncbi")
        uniprot_aliases = cache.get_cached_aliases(gene, "uniprot")
        
        if ncbi_aliases is not None or uniprot_aliases is not None:
            # Use cached results
            cache_hits += 1
            is_cached = True
            aliases_ncbi = ncbi_aliases if ncbi_aliases is not None else []
            aliases_uniprot = uniprot_aliases if uniprot_aliases is not None else []
        else:
            # Need to fetch
            new_lookups += 1
            is_cached = False
            logging.info(f"Processing gene {gene} ({processed + 1}/{len(genes)}) - NEW LOOKUP")
            
            # Only call fetch functions for new lookups
            aliases_ncbi = fetch_gene_aliases_ncbi(gene, cache=cache)
            aliases_uniprot = fetch_gene_aliases_uniprot(gene, cache=cache)
        
        # Combine and deduplicate
        new_aliases = list(set(aliases_ncbi + aliases_uniprot))
        
        # For now, just store the raw aliases without consolidation
        # We'll do consolidation in a single pass at the end
        if new_aliases:
            all_aliases[gene] = sorted(new_aliases)
            # Only log details for new lookups
            if not is_cached:
                logging.info(f"  Aliases for {gene}: {', '.join(new_aliases)}")
        else:
            all_aliases[gene] = []
            if not is_cached:
                logging.info(f"  No aliases found for {gene}")
        
        processed += 1
        
        # No rate limiting needed - API calls have their own rate limiting
        
        # Progress update every 100 genes when using cache
        if is_cached and processed % 100 == 0:
            logging.info(f"Progress: {processed}/{len(genes)} genes processed ({cache_hits} from cache, {new_lookups} new lookups)")
    
    # Add summary before consolidation
    logging.info(f"\n=== Consolidating aliases ===")
    logging.info(f"Starting consolidation for {len(all_aliases)} genes...")
    
    # Now do consolidation in a single efficient pass
    gene_to_primary = {}
    consolidated_aliases = {}
    
    # First pass: build gene_to_primary mapping
    for gene, aliases in all_aliases.items():
        # Add the gene itself to the mapping
        if gene not in gene_to_primary:
            gene_to_primary[gene] = gene
        
        # Add aliases
        for alias in aliases:
            if alias not in gene_to_primary:
                gene_to_primary[alias] = gene
            elif gene_to_primary[alias] != gene:
                # This alias belongs to multiple genes - keep the first one
                logging.debug(f"Alias conflict: {alias} claimed by both {gene_to_primary[alias]} and {gene}")
    
    # Second pass: consolidate aliases under primary genes
    for gene in all_aliases:
        primary = gene_to_primary.get(gene, gene)
        aliases = all_aliases[gene]
        
        if primary not in consolidated_aliases:
            consolidated_aliases[primary] = set()
        
        # Add all aliases
        consolidated_aliases[primary].update(aliases)
        
        # Also add the gene itself if it's different from primary
        if gene != primary:
            consolidated_aliases[primary].add(gene)
    
    # Clean up: convert sets to sorted lists and remove self-references
    for primary in consolidated_aliases:
        alias_set = consolidated_aliases[primary]
        alias_set.discard(primary)  # Remove self-reference
        consolidated_aliases[primary] = sorted(list(alias_set))
    
    # Replace all_aliases with consolidated version
    all_aliases = consolidated_aliases
    
    logging.info(f"Consolidation complete: {len(all_aliases)} primary genes")
    
    # Final cleanup: only include genes that are NOT aliases of other genes
    # Create reverse mapping: alias -> primary gene
    alias_to_primary = {}
    for primary, aliases in all_aliases.items():
        for alias in aliases:
            alias_to_primary[alias] = primary
    
    final_aliases = {}
    excluded_as_aliases = []
    
    for gene in genes:
        # Check if this gene appears as an alias of another gene
        if gene in alias_to_primary and alias_to_primary[gene] != gene:
            # This gene is an alias of another gene, exclude it
            excluded_as_aliases.append(gene)
            logging.info(f"Excluding {gene} - it's an alias of {alias_to_primary[gene]}")
        else:
            # Include this gene with its aliases
            primary = gene_to_primary.get(gene, gene)
            if primary in all_aliases:
                # Remove the gene itself from its alias list
                aliases = [alias for alias in all_aliases[primary] if alias != gene]
                final_aliases[gene] = aliases
            else:
                final_aliases[gene] = []
    
    # Write results to file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write("# Gene Aliases/Synonyms (Consolidated)\n")
        f.write("# Format: gene_symbol:alias1,alias2,alias3\n")
        f.write(f"# Generated on: {pd.Timestamp.now()}\n")
        f.write(f"# Total input genes: {len(all_genes)}\n")
        f.write(f"# Filtered genes processed: {len(genes)}\n")
        f.write(f"# Excluded (start with capital): {len(excluded_genes)}\n")
        f.write(f"# Excluded (are aliases of other genes): {len(excluded_as_aliases)}\n")
        f.write(f"# Final genes in output: {len(final_aliases)}\n")
        f.write(f"# Genes with aliases: {len([g for g in final_aliases.values() if g])}\n\n")
        
        for gene in sorted(final_aliases.keys()):
            if final_aliases[gene]:
                aliases_str = ','.join(final_aliases[gene])
                f.write(f"{gene}:{aliases_str}\n")
            else:
                f.write(f"{gene}:\n")
    
    # Save excluded genes for reference
    excluded_file = output_path.with_name(output_path.stem + '_excluded.txt')
    with open(excluded_file, 'w') as f:
        f.write("# Genes excluded from alias search\n")
        f.write(f"# Generated on: {pd.Timestamp.now()}\n")
        f.write(f"# Excluded (capital letters): {len(excluded_genes)}\n")
        f.write(f"# Excluded (are aliases): {len(excluded_as_aliases)}\n")
        f.write(f"# Total excluded: {len(excluded_genes) + len(excluded_as_aliases)}\n\n")
        
        if excluded_genes:
            f.write("# Genes starting with capital letters:\n")
            for gene in sorted(excluded_genes):
                f.write(f"{gene}\n")
        
        if excluded_as_aliases:
            f.write(f"\n# Genes that are aliases of other genes:\n")
            for gene in sorted(excluded_as_aliases):
                primary = alias_to_primary[gene]
                f.write(f"{gene} -> alias of {primary}\n")
    
    # Also save as JSON for easier programmatic access
    json_file = output_path.with_suffix('.json')
    with open(json_file, 'w') as f:
        json.dump(final_aliases, f, indent=2)
    
    # Save final cache
    cache.save_cache()
    
    # Add final summary of cache usage
    logging.info(f"\n=== Gene Alias Processing Complete ===")
    logging.info(f"Cache performance:")
    logging.info(f"  Cache hits: {cache_hits}")
    logging.info(f"  New lookups: {new_lookups}")
    if processed > 0:
        cache_hit_rate = (cache_hits / processed) * 100
        logging.info(f"  Cache hit rate: {cache_hit_rate:.1f}%")
    
    logging.info(f"\nResults saved to:")
    logging.info(f"  Aliases: {output_path}")
    logging.info(f"  JSON: {json_file}")
    logging.info(f"  Excluded: {excluded_file}")
    logging.info(f"  Cache: {CACHE_FILE}")
    
    genes_with_aliases = len([g for g in final_aliases.values() if g])
    logging.info(f"\nProcessing summary:")
    logging.info(f"  Input genes: {len(all_genes)}")
    logging.info(f"  Filtered for processing: {len(genes)}")
    logging.info(f"  Excluded (capital letters): {len(excluded_genes)}")
    logging.info(f"  Excluded (are aliases): {len(excluded_as_aliases)}")
    logging.info(f"  Final genes in output: {len(final_aliases)}")
    logging.info(f"  Found aliases for: {genes_with_aliases} genes")
    
    # Only print detailed alias list if there were new lookups
    if new_lookups > 0:
        # Print summary of interesting aliases
        interesting_aliases = {gene: aliases for gene, aliases in final_aliases.items() if aliases}
        if interesting_aliases:
            logging.info("\nNew genes with aliases found:")
            new_genes_with_aliases = 0
            for gene, aliases in interesting_aliases.items():
                # Only show if it was a new lookup
                if cache.get_cached_aliases(gene, "ncbi") is None or cache.get_cached_aliases(gene, "uniprot") is None:
                    logging.info(f"  {gene} -> {', '.join(aliases)}")
                    new_genes_with_aliases += 1
                    if new_genes_with_aliases >= 10:
                        logging.info("  ... (showing first 10 new genes with aliases)")
                        break

if __name__ == "__main__":
    main()