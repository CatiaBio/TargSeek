#!/usr/bin/env python3
"""
PDB 3D Structure Downloader - Bacterial Only
=============================================

Downloads 3D protein structures and their corresponding FASTA sequences from the Protein Data Bank (PDB)
for genes listed in the proteins_to_study file. Only searches for bacterial proteins (taxonomy ID 2).

Uses UniProt to find bacterial proteins with 3D structures, then downloads from PDB.
"""

import requests
import pandas as pd
from pathlib import Path
import logging
import time
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import urllib.parse
from datetime import datetime
import concurrent.futures
from functools import partial

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Cache configuration - use absolute paths to avoid working directory issues
# Get the project root directory (assuming script is in scripts/gene_selection/)
PROJECT_ROOT = Path(__file__).parent.parent.parent
CACHE_DIR = PROJECT_ROOT / "cache" / "protein_structures"
CACHE_FILE = CACHE_DIR / "protein_3d_structure_cache.json"
GENE_DOWNLOAD_CACHE_FILE = CACHE_DIR / "gene_downloads_cache.json"

class GeneDownloadCache:
    """Manages caching of gene-level download completion status"""
    
    def __init__(self):
        self.cache = {}
        self.cache_updates = 0
        self.load_cache()
    
    def load_cache(self):
        """Load existing gene download cache"""
        if GENE_DOWNLOAD_CACHE_FILE.exists():
            try:
                with open(GENE_DOWNLOAD_CACHE_FILE, 'r') as f:
                    self.cache = json.load(f)
                logging.info(f"✓ Loaded gene download cache with {len(self.cache)} entries")
            except Exception as e:
                logging.warning(f"Could not load gene download cache: {e}")
                self.cache = {}
        else:
            logging.info("Creating new gene download cache")
            self.cache = {}
    
    def save_cache(self):
        """Save cache to disk if there were updates"""
        if self.cache_updates > 0:
            try:
                CACHE_DIR.mkdir(parents=True, exist_ok=True)
                with open(GENE_DOWNLOAD_CACHE_FILE, 'w') as f:
                    json.dump(self.cache, f, indent=2)
                logging.info(f"✓ Saved gene download cache with {self.cache_updates} new/updated entries")
                self.cache_updates = 0
            except Exception as e:
                logging.error(f"Failed to save gene download cache: {e}")
    
    def is_gene_completed(self, gene_name):
        """Check if gene download is already completed"""
        return self.cache.get(gene_name, {}).get('completed', False)
    
    def get_gene_structures(self, gene_name):
        """Get list of downloaded structures for a gene"""
        return self.cache.get(gene_name, {}).get('downloaded_structures', [])
    
    def mark_gene_completed(self, gene_name, structures_info):
        """Mark gene as completed with structure information"""
        self.cache[gene_name] = {
            'completed': True,
            'completion_date': datetime.now().isoformat(),
            'downloaded_structures': structures_info.get('structure_ids', []),
            'experimental_count': structures_info.get('experimental_count', 0),
            'computed_count': structures_info.get('computed_count', 0),
            'sequences_count': structures_info.get('sequences', 0),
            'structures_count': structures_info.get('structures', 0)
        }
        self.cache_updates += 1
        logging.info(f"Marked gene {gene_name} as completed in cache with {len(structures_info.get('structure_ids', []))} structures")

class Protein3DStructureCache:
    """Manages caching of protein 3D structure searches and downloads"""
    
    def __init__(self):
        self.cache = {}
        self.cache_updates = 0
        self.load_cache()
    
    def load_cache(self):
        """Load existing 3D structure cache"""
        if CACHE_FILE.exists():
            try:
                with open(CACHE_FILE, 'r') as f:
                    self.cache = json.load(f)
                logging.info(f"✓ Loaded 3D structure cache with {len(self.cache)} entries")
            except Exception as e:
                logging.warning(f"Could not load 3D structure cache: {e}")
                self.cache = {}
        else:
            logging.info("Creating new 3D structure cache")
            self.cache = {}
    
    def save_cache(self):
        """Save cache to disk if there were updates"""
        if self.cache_updates > 0:
            try:
                CACHE_DIR.mkdir(parents=True, exist_ok=True)
                with open(CACHE_FILE, 'w') as f:
                    json.dump(self.cache, f, indent=2)
                logging.info(f"✓ Saved 3D structure cache with {self.cache_updates} new entries")
                self.cache_updates = 0
            except Exception as e:
                logging.error(f"Failed to save 3D structure cache: {e}")
    
    def get_cache_key(self, gene_name, database="uniprot"):
        """Generate cache key for 3D structure search"""
        return f"{gene_name}||{database}||bacterial"
    
    def get_pdb_cache_key(self, pdb_id):
        """Generate cache key for PDB structure download"""
        return f"pdb||{pdb_id}||structure"
    
    def get_cached_search(self, gene_name, database="uniprot"):
        """Get cached 3D structure search results"""
        key = self.get_cache_key(gene_name, database)
        return self.cache.get(key)
    
    def cache_search_result(self, gene_name, pdb_ids, database="uniprot"):
        """Cache 3D structure search results"""
        key = self.get_cache_key(gene_name, database)
        self.cache[key] = {
            "pdb_ids": pdb_ids,
            "found": len(pdb_ids) > 0,
            "search_date": datetime.now().isoformat(),
            "database": database
        }
        self.cache_updates += 1
    
    def get_cached_pdb(self, pdb_id):
        """Get cached PDB download info"""
        key = self.get_pdb_cache_key(pdb_id)
        return self.cache.get(key)
    
    def cache_pdb_result(self, pdb_id, success, file_path=None, sequences_count=0):
        """Cache PDB download results"""
        key = self.get_pdb_cache_key(pdb_id)
        self.cache[key] = {
            "success": success,
            "file_path": str(file_path) if file_path else None,
            "sequences_count": sequences_count,
            "download_date": datetime.now().isoformat()
        }
        self.cache_updates += 1

def load_gene_aliases(aliases_file):
    """Load gene aliases from JSON or TSV file"""
    try:
        if aliases_file.suffix == '.json':
            with open(aliases_file, 'r') as f:
                aliases = json.load(f)
            logging.info(f"Loaded aliases for {len(aliases)} genes from JSON")
            return aliases
        elif aliases_file.suffix == '.tsv':
            # Load from TSV format
            df = pd.read_csv(aliases_file, sep='\t')
            aliases = {}
            for _, row in df.iterrows():
                gene = row['gene']
                alias_list = row.get('alias_list', '')
                if pd.isna(alias_list) or alias_list == '':
                    aliases[gene] = []
                else:
                    aliases[gene] = [alias.strip() for alias in alias_list.split(',')]
            logging.info(f"Loaded aliases for {len(aliases)} genes from TSV")
            return aliases
        else:
            logging.warning(f"Unsupported alias file format: {aliases_file}")
            return {}
    except Exception as e:
        logging.warning(f"Could not load gene aliases: {e}")
        return {}

def search_pdb_via_uniprot_bacterial(gene_name, max_structures=10, gene_aliases=None, cache=None, prioritize_extracellular=True, extracellular_keywords=None, intracellular_keywords=None):
    """
    Search UniProt for bacterial proteins with PDB structures, with alias fallback
    
    Args:
        gene_name (str): Primary gene name to search for
        max_structures (int): Maximum structures to return
        gene_aliases (list): Alternative gene names to try if primary fails
        cache (Protein3DStructureCache): Cache instance to use
        prioritize_extracellular (bool): Whether to prioritize extracellular proteins
        extracellular_keywords (list): Keywords indicating extracellular localization
        intracellular_keywords (list): Keywords indicating intracellular localization
    """
    # Check cache first
    if cache:
        cached_result = cache.get_cached_search(gene_name)
        if cached_result:
            logging.info(f"Found cached 3D structure search for {gene_name}: {len(cached_result['pdb_ids'])} structures")
            return cached_result['pdb_ids'][:max_structures]
    
    # Try primary gene name first
    pdb_ids = _search_uniprot_bacterial_single_gene(gene_name, max_structures, prioritize_extracellular, extracellular_keywords, intracellular_keywords)
    
    if pdb_ids:
        logging.info(f"Found {len(pdb_ids)} bacterial PDB structures via UniProt for gene {gene_name}")
        if cache:
            cache.cache_search_result(gene_name, pdb_ids)
        return pdb_ids
    
    # If no results and we have aliases, try them
    if gene_aliases:
        logging.info(f"No bacterial structures found for {gene_name}, trying {len(gene_aliases)} aliases")
        for alias in gene_aliases:
            if alias != gene_name:  # Don't retry the same name
                alias_pdb_ids = _search_uniprot_bacterial_single_gene(alias, max_structures, prioritize_extracellular, extracellular_keywords, intracellular_keywords)
                if alias_pdb_ids:
                    logging.info(f"Found {len(alias_pdb_ids)} bacterial PDB structures via alias '{alias}'")
                    if cache:
                        cache.cache_search_result(gene_name, alias_pdb_ids)
                    return alias_pdb_ids
    
    logging.info(f"No bacterial PDB structures found for {gene_name} (tried {len(gene_aliases or [])} aliases)")
    if cache:
        cache.cache_search_result(gene_name, [])  # Cache empty result
    return []

def score_protein_localization(entry, extracellular_keywords, intracellular_keywords):
    """Score a protein entry based on cellular localization for extracellular prioritization"""
    score = 0
    
    # Get relevant fields to search for localization info
    fields_to_search = []
    
    # Protein name and description
    if 'proteinDescription' in entry:
        protein_desc = entry['proteinDescription']
        if 'recommendedName' in protein_desc:
            if 'fullName' in protein_desc['recommendedName']:
                fields_to_search.append(protein_desc['recommendedName']['fullName']['value'].lower())
        if 'alternativeNames' in protein_desc:
            for alt_name in protein_desc['alternativeNames']:
                if 'fullName' in alt_name:
                    fields_to_search.append(alt_name['fullName']['value'].lower())
    
    # Comments (includes subcellular location)
    if 'comments' in entry:
        for comment in entry['comments']:
            if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                if 'texts' in comment:
                    for text in comment['texts']:
                        if 'value' in text:
                            fields_to_search.append(text['value'].lower())
            elif comment.get('commentType') == 'FUNCTION':
                if 'texts' in comment:
                    for text in comment['texts']:
                        if 'value' in text:
                            fields_to_search.append(text['value'].lower())
    
    # Keywords
    if 'keywords' in entry:
        for keyword in entry['keywords']:
            if 'value' in keyword:
                fields_to_search.append(keyword['value'].lower())
    
    # Gene names
    if 'genes' in entry:
        for gene in entry['genes']:
            if 'geneName' in gene and 'value' in gene['geneName']:
                fields_to_search.append(gene['geneName']['value'].lower())
    
    # Search for extracellular keywords (positive score)
    search_text = ' '.join(fields_to_search)
    for i, keyword in enumerate(extracellular_keywords):
        if keyword.lower() in search_text:
            # Higher priority keywords get higher scores
            keyword_score = len(extracellular_keywords) - i
            score += keyword_score * 10
            logging.debug(f"Extracellular keyword '{keyword}' found, +{keyword_score * 10} points")
    
    # Search for intracellular keywords (negative score)
    for keyword in intracellular_keywords:
        if keyword.lower() in search_text:
            score -= 5  # Penalty for intracellular keywords
            logging.debug(f"Intracellular keyword '{keyword}' found, -5 points")
    
    return score

def get_pdb_release_dates(pdb_ids):
    """
    Get release dates for a list of PDB IDs
    
    Args:
        pdb_ids (list): List of PDB identifiers
    
    Returns:
        dict: Dictionary mapping PDB ID to release date (or None if failed)
    """
    release_dates = {}
    
    for pdb_id in pdb_ids:
        try:
            # Get PDB metadata including release date
            api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            response = requests.get(api_url, timeout=15)
            response.raise_for_status()
            
            data = response.json()
            
            # Try different possible date fields
            release_date = None
            
            # Method 1: revision_date (most recent update)
            if 'rcsb_accession_info' in data:
                if 'revision_date' in data['rcsb_accession_info']:
                    release_date = data['rcsb_accession_info']['revision_date']
                elif 'initial_release_date' in data['rcsb_accession_info']:
                    release_date = data['rcsb_accession_info']['initial_release_date']
            
            # Method 2: struct entry revision date
            if not release_date and 'pdbx_database_status' in data:
                if 'recvd_initial_deposition_date' in data['pdbx_database_status']:
                    release_date = data['pdbx_database_status']['recvd_initial_deposition_date']
            
            release_dates[pdb_id] = release_date
            logging.debug(f"PDB {pdb_id} release date: {release_date}")
                
        except Exception as e:
            logging.debug(f"Could not get release date for {pdb_id}: {e}")
            release_dates[pdb_id] = None
    
    return release_dates

def filter_unique_pdb_prefixes(pdb_with_dates, max_structures):
    """
    Filter PDB structures to keep only those with unique first digit prefixes
    
    Args:
        pdb_with_dates (list): List of (pdb_id, release_date) tuples, sorted by date
        max_structures (int): Maximum number of structures to return
    
    Returns:
        list: Filtered list of PDB IDs with unique prefixes (first digit)
    """
    if not pdb_with_dates:
        return []
    
    seen_prefixes = set()
    unique_structures = []
    
    for pdb_id, release_date in pdb_with_dates:
        # Get the first character of PDB ID (the digit)
        first_digit = pdb_id[0] if pdb_id else ''
        
        if first_digit not in seen_prefixes:
            # This is a new prefix, include it
            seen_prefixes.add(first_digit)
            unique_structures.append(pdb_id)
            logging.debug(f"Including {pdb_id} (new prefix: {first_digit})")
        else:
            # This prefix was already seen, skip it
            logging.debug(f"Excluding {pdb_id} (duplicate prefix: {first_digit})")
        
        # Stop when we have enough unique structures
        if len(unique_structures) >= max_structures:
            break
    
    return unique_structures

def _search_uniprot_bacterial_single_gene(gene_name, max_structures, prioritize_extracellular=True, extracellular_keywords=None, intracellular_keywords=None):
    """Helper function to search UniProt for bacterial proteins with structures for a single gene name"""
    
    # Default keywords if not provided
    if extracellular_keywords is None:
        extracellular_keywords = ["extracellular", "outer membrane", "cell surface", "secreted", "periplasm", "cell wall", "membrane", "surface"]
    if intracellular_keywords is None:
        intracellular_keywords = ["cytoplasm", "cytosol", "intracellular", "ribosome", "nucleus", "nucleoid"]
    
    try:
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/search"
        
        # Build basic query without protein name
        query = f'gene:{gene_name} AND structure_3d:true'
        
        params = {
            'query': query,
            'format': 'json',
            'size': 500,  # Get many results to find all available bacterial structures
            'fields': 'accession,gene_names,organism_name,xref_pdb,protein_name,cc_subcellular_location,cc_function,keyword'
        }
        
        logging.info(f"UniProt search for {gene_name} (gene name only)")
        logging.debug(f"Query: {query}")
        
        response = requests.get(uniprot_url, params=params, timeout=30)
        if response.status_code != 200:
            logging.debug(f"UniProt bacterial search failed for {gene_name}: status {response.status_code}")
            return []
        
        data = response.json()
        bacterial_entries = []
        other_entries = []
        
        # Common bacterial genera for identification
        bacterial_genera = ['escherichia', 'bacillus', 'staphylococcus', 'streptococcus',
                          'pseudomonas', 'salmonella', 'listeria', 'mycobacterium',
                          'clostridium', 'enterococcus', 'lactococcus', 'vibrio',
                          'corynebacterium', 'actinomyces', 'bifidobacterium']
        
        for entry in data.get('results', []):
            organism = entry.get('organism', {}).get('scientificName', 'Unknown')
            organism_lower = organism.lower()
            
            # Check if this looks like a bacterial organism
            is_bacterial = any(genus in organism_lower for genus in bacterial_genera)
            
            # Look for PDB cross-references
            pdb_refs = []
            for xref in entry.get('uniProtKBCrossReferences', []):
                if xref.get('database') == 'PDB':
                    pdb_id = xref.get('id')
                    if pdb_id:
                        pdb_refs.append(pdb_id)
            
            if pdb_refs:
                # Score for extracellular prioritization
                localization_score = 0
                if prioritize_extracellular:
                    localization_score = score_protein_localization(entry, extracellular_keywords, intracellular_keywords)
                
                entry_info = {
                    'entry': entry,
                    'pdb_ids': pdb_refs,
                    'organism': organism,
                    'is_bacterial': is_bacterial,
                    'localization_score': localization_score
                }
                
                if is_bacterial:
                    bacterial_entries.append(entry_info)
                    logging.debug(f"Found bacterial protein from {organism} with PDB: {pdb_refs[0]} (score: {localization_score})")
                else:
                    other_entries.append(entry_info)
                    logging.debug(f"Found non-bacterial protein from {organism}")
        
        # Sort and select bacterial structures
        if bacterial_entries:
            # Extract all PDB IDs first
            all_bacterial_pdb_ids = []
            for entry_info in bacterial_entries:
                for pdb_id in entry_info['pdb_ids']:
                    if pdb_id not in all_bacterial_pdb_ids:
                        all_bacterial_pdb_ids.append(pdb_id)
            
            # Get release dates for all PDB IDs
            logging.info(f"Getting release dates for {len(all_bacterial_pdb_ids)} bacterial structures for {gene_name}")
            release_dates = get_pdb_release_dates(all_bacterial_pdb_ids)
            
            # Create list of (pdb_id, release_date) tuples for sorting
            pdb_with_dates = []
            for pdb_id in all_bacterial_pdb_ids:
                release_date = release_dates.get(pdb_id)
                pdb_with_dates.append((pdb_id, release_date))
            
            # Sort by release date (most recent first), putting None dates at the end
            def sort_by_date(item):
                pdb_id, release_date = item
                if release_date is None:
                    return ('9999-12-31', pdb_id)  # Put entries without dates at the end
                return (release_date, pdb_id)
            
            pdb_with_dates.sort(key=sort_by_date, reverse=True)
            
            # Log top entries by date (before filtering)
            logging.info(f"Structures for {gene_name} sorted by release date (most recent first):")
            for i, (pdb_id, release_date) in enumerate(pdb_with_dates[:5]):
                date_str = release_date if release_date else "Unknown"
                logging.info(f"  Rank {i+1}: {pdb_id} (released: {date_str})")
            
            # Filter for unique PDB prefixes (fast filtering by first digit)
            logging.info(f"Filtering for unique PDB prefixes among {len(pdb_with_dates)} bacterial structures for {gene_name}")
            unique_bacterial_pdb_ids = filter_unique_pdb_prefixes(pdb_with_dates, max_structures)
            
            logging.info(f"Found {len(pdb_with_dates)} bacterial structures for {gene_name}, "
                        f"filtered to {len(unique_bacterial_pdb_ids)} unique prefixes")
            return unique_bacterial_pdb_ids
        elif other_entries:
            # Apply same date-based prioritization to non-bacterial entries
            all_other_pdb_ids = []
            for entry_info in other_entries:
                for pdb_id in entry_info['pdb_ids']:
                    if pdb_id not in all_other_pdb_ids:
                        all_other_pdb_ids.append(pdb_id)
            
            # Get release dates for non-bacterial structures
            logging.info(f"Getting release dates for {len(all_other_pdb_ids)} non-bacterial structures for {gene_name}")
            release_dates = get_pdb_release_dates(all_other_pdb_ids)
            
            # Create and sort by date
            pdb_with_dates = []
            for pdb_id in all_other_pdb_ids:
                release_date = release_dates.get(pdb_id)
                pdb_with_dates.append((pdb_id, release_date))
            
            def sort_by_date(item):
                pdb_id, release_date = item
                if release_date is None:
                    return ('9999-12-31', pdb_id)
                return (release_date, pdb_id)
            
            pdb_with_dates.sort(key=sort_by_date, reverse=True)
            
            # Filter for unique PDB prefixes among non-bacterial structures too
            logging.info(f"Filtering for unique PDB prefixes among {len(pdb_with_dates)} non-bacterial structures for {gene_name}")
            unique_other_pdb_ids = filter_unique_pdb_prefixes(pdb_with_dates, max_structures)
            
            logging.info(f"No bacterial structures found for {gene_name}, using {len(unique_other_pdb_ids)} unique non-bacterial structures by release date")
            return unique_other_pdb_ids
        else:
            return []
        
    except Exception as e:
        logging.warning(f"UniProt bacterial search failed for {gene_name}: {e}")
        return []

def search_alphafold_models(gene_name, gene_aliases=None, cache=None):
    """
    Search for AlphaFold computed models for a gene
    
    Args:
        gene_name (str): Gene name to search for
        gene_aliases (list): Alternative gene names to try
        cache (Protein3DStructureCache): Cache instance
    
    Returns:
        list: List of AlphaFold model identifiers
    """
    # Check cache first
    if cache:
        cache_key = f"{gene_name}||alphafold||bacterial"
        cached_result = cache.cache.get(cache_key)
        if cached_result:
            logging.info(f"Found cached AlphaFold search for {gene_name}: {len(cached_result.get('pdb_ids', []))} models")
            return cached_result.get('pdb_ids', [])
    
    alphafold_ids = _search_alphafold_single_gene(gene_name)
    
    if alphafold_ids:
        logging.info(f"Found {len(alphafold_ids)} AlphaFold models for gene {gene_name}")
        if cache:
            cache.cache[f"{gene_name}||alphafold||bacterial"] = {
                "pdb_ids": alphafold_ids,
                "found": len(alphafold_ids) > 0,
                "search_date": datetime.now().isoformat(),
                "database": "alphafold"
            }
            cache.cache_updates += 1
        return alphafold_ids
    
    # Try aliases if available
    if gene_aliases:
        logging.info(f"No AlphaFold models found for {gene_name}, trying {len(gene_aliases)} aliases")
        for alias in gene_aliases:
            if alias != gene_name:
                alias_models = _search_alphafold_single_gene(alias)
                if alias_models:
                    logging.info(f"Found {len(alias_models)} AlphaFold models via alias '{alias}'")
                    if cache:
                        cache.cache[f"{gene_name}||alphafold||bacterial"] = {
                            "pdb_ids": alias_models,
                            "found": len(alias_models) > 0,
                            "search_date": datetime.now().isoformat(),
                            "database": "alphafold"
                        }
                        cache.cache_updates += 1
                    return alias_models
    
    logging.info(f"No AlphaFold models found for {gene_name}")
    if cache:
        cache.cache[f"{gene_name}||alphafold||bacterial"] = {
            "pdb_ids": [],
            "found": False,
            "search_date": datetime.now().isoformat(),
            "database": "alphafold"
        }
        cache.cache_updates += 1
    return []

def _search_alphafold_single_gene(gene_name):
    """
    Search AlphaFold database for a specific gene
    """
    try:
        # Search UniProt for the gene to get UniProt accessions
        uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'gene:{gene_name}',
            'format': 'json',
            'size': 50,
            'fields': 'accession,gene_names,organism_name'
        }
        
        response = requests.get(uniprot_url, params=params, timeout=30)
        if response.status_code != 200:
            return []
        
        data = response.json()
        alphafold_ids = []
        
        # Check each UniProt entry for AlphaFold models
        for entry in data.get('results', []):
            accession = entry.get('primaryAccession')
            organism = entry.get('organism', {}).get('scientificName', '')
            
            # Focus on bacterial organisms
            bacterial_genera = ['escherichia', 'bacillus', 'staphylococcus', 'streptococcus',
                              'pseudomonas', 'salmonella', 'listeria', 'mycobacterium',
                              'clostridium', 'enterococcus', 'lactococcus', 'vibrio']
            
            is_bacterial = any(genus in organism.lower() for genus in bacterial_genera)
            
            if accession and is_bacterial:
                # Check if AlphaFold model exists
                alphafold_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
                try:
                    af_response = requests.get(alphafold_url, timeout=10)
                    if af_response.status_code == 200:
                        af_data = af_response.json()
                        if af_data:
                            # Use the UniProt accession as the identifier
                            alphafold_ids.append(f"AF-{accession}")
                            logging.debug(f"Found AlphaFold model for {accession} ({organism})")
                except:
                    continue
        
        return alphafold_ids[:10]  # Limit to 10 models max
        
    except Exception as e:
        logging.warning(f"AlphaFold search failed for {gene_name}: {e}")
        return []

def get_pdb_fasta(pdb_id):
    """
    Download FASTA sequence for a PDB structure
    
    Args:
        pdb_id (str): PDB identifier
    
    Returns:
        list: List of SeqRecord objects
    """
    
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    
    try:
        response = requests.get(fasta_url, timeout=30)
        response.raise_for_status()
        
        # Parse FASTA content
        sequences = []
        lines = response.text.strip().split('\n')
        current_seq = ""
        current_header = ""
        
        for line in lines:
            if line.startswith('>'):
                if current_header and current_seq:
                    # Create SeqRecord
                    seq_id = current_header.split('|')[0].replace('>', '')
                    sequences.append(SeqRecord(
                        Seq(current_seq),
                        id=seq_id,
                        description=current_header.replace('>', '')
                    ))
                current_header = line
                current_seq = ""
            else:
                current_seq += line.strip()
        
        # Add last sequence
        if current_header and current_seq:
            seq_id = current_header.split('|')[0].replace('>', '')
            sequences.append(SeqRecord(
                Seq(current_seq),
                id=seq_id,
                description=current_header.replace('>', '')
            ))
        
        logging.info(f"Downloaded {len(sequences)} sequences for PDB {pdb_id}")
        return sequences
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download FASTA for {pdb_id}: {e}")
        return []

def get_alphafold_fasta(alphafold_id):
    """
    Get FASTA sequence for an AlphaFold model
    
    Args:
        alphafold_id (str): AlphaFold identifier (e.g., AF-P12345)
    
    Returns:
        list: List of SeqRecord objects
    """
    uniprot_accession = alphafold_id.replace('AF-', '')
    
    try:
        # Get sequence from UniProt
        uniprot_fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.fasta"
        response = requests.get(uniprot_fasta_url, timeout=30)
        response.raise_for_status()
        
        # Parse FASTA content
        sequences = []
        lines = response.text.strip().split('\n')
        current_seq = ""
        current_header = ""
        
        for line in lines:
            if line.startswith('>'):
                if current_header and current_seq:
                    seq_id = current_header.split('|')[1] if '|' in current_header else alphafold_id
                    sequences.append(SeqRecord(
                        Seq(current_seq),
                        id=f"{alphafold_id}_{seq_id}",
                        description=f"AlphaFold:{alphafold_id} {current_header.replace('>', '')}"
                    ))
                current_header = line
                current_seq = ""
            else:
                current_seq += line.strip()
        
        # Add last sequence
        if current_header and current_seq:
            seq_id = current_header.split('|')[1] if '|' in current_header else alphafold_id
            sequences.append(SeqRecord(
                Seq(current_seq),
                id=f"{alphafold_id}_{seq_id}",
                description=f"AlphaFold:{alphafold_id} {current_header.replace('>', '')}"
            ))
        
        logging.info(f"Downloaded {len(sequences)} sequences for AlphaFold model {alphafold_id}")
        return sequences
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download FASTA for AlphaFold model {alphafold_id}: {e}")
        return []

def download_single_pdb_structure(pdb_id, output_dir, cache=None):
    """
    Download a single PDB structure - helper function for concurrent downloads
    """
    return download_pdb_structure(pdb_id, output_dir, cache)

def download_alphafold_model(alphafold_id, output_dir, cache=None):
    """
    Download AlphaFold computed model
    
    Args:
        alphafold_id (str): AlphaFold identifier (e.g., AF-P12345)
        output_dir (Path): Output directory
        cache (Protein3DStructureCache): Cache instance
    
    Returns:
        bool: True if download succeeded
    """
    # Extract UniProt accession from AlphaFold ID
    uniprot_accession = alphafold_id.replace('AF-', '')
    
    # Check if file already exists
    pdb_file = output_dir / f"{alphafold_id}.pdb"
    if pdb_file.exists() and pdb_file.stat().st_size > 0:
        if cache:
            cache.cache_pdb_result(alphafold_id, True, pdb_file)
        return True
    
    # Check cache
    if cache:
        cached_af = cache.get_cached_pdb(alphafold_id)
        if cached_af and cached_af.get('success'):
            cached_path = Path(cached_af.get('file_path', ''))
            if cached_path.exists() and cached_path.stat().st_size > 0:
                return True
    
    try:
        # Download AlphaFold model
        af_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-model_v4.pdb"
        response = requests.get(af_url, timeout=30)
        response.raise_for_status()
        
        # Save model file (uncompressed for AlphaFold)
        with open(pdb_file, 'w') as f:
            f.write(response.text)
        
        logging.info(f"Downloaded AlphaFold model: {pdb_file}")
        if cache:
            cache.cache_pdb_result(alphafold_id, True, pdb_file)
        return True
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download AlphaFold model for {alphafold_id}: {e}")
        if cache:
            cache.cache_pdb_result(alphafold_id, False)
        return False

def download_pdb_structure(pdb_id, output_dir, cache=None):
    """
    Download PDB structure files in both legacy format (gz) and PDBML/XML format (gz) with caching
    
    Args:
        pdb_id (str): PDB identifier
        output_dir (Path): Output directory to save the structure
        cache (Protein3DStructureCache): Cache instance to use
    
    Returns:
        bool: True if download succeeded, False otherwise
    """
    
    # Define file paths for both formats
    pdb_file = output_dir / f"{pdb_id}.pdb.gz"
    pdbml_file = output_dir / f"{pdb_id}.xml.gz"
    
    # Check if both files already exist (quick check first)
    pdb_exists = pdb_file.exists() and pdb_file.stat().st_size > 0
    pdbml_exists = pdbml_file.exists() and pdbml_file.stat().st_size > 0
    
    if pdb_exists and pdbml_exists:
        # Both files exist and are not empty, no need to re-download
        if cache:
            cache.cache_pdb_result(pdb_id, True, pdb_file)
        return True
    
    # Check cache for download status only if files don't exist
    if cache:
        cached_pdb = cache.get_cached_pdb(pdb_id)
        if cached_pdb and cached_pdb.get('success'):
            cached_path = Path(cached_pdb.get('file_path', ''))
            cached_pdbml_path = cached_path.parent / f"{pdb_id}.xml.gz"
            if (cached_path.exists() and cached_path.stat().st_size > 0 and 
                cached_pdbml_path.exists() and cached_pdbml_path.stat().st_size > 0):
                return True
    
    # URLs for both formats
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    pdbml_url = f"https://files.rcsb.org/download/{pdb_id}.xml.gz"
    
    success_count = 0
    
    try:
        # Download PDB format (legacy format)
        if not pdb_exists:
            response = requests.get(pdb_url, timeout=30)
            response.raise_for_status()
            
            with open(pdb_file, 'wb') as f:
                f.write(response.content)
            
            logging.info(f"Downloaded PDB structure: {pdb_file}")
            success_count += 1
        else:
            logging.info(f"PDB file already exists: {pdb_file}")
            success_count += 1
        
        # Download PDBML/XML format (contains correct sequence numbering)
        if not pdbml_exists:
            response = requests.get(pdbml_url, timeout=30)
            response.raise_for_status()
            
            with open(pdbml_file, 'wb') as f:
                f.write(response.content)
            
            logging.info(f"Downloaded PDBML/XML structure: {pdbml_file}")
            success_count += 1
        else:
            logging.info(f"PDBML file already exists: {pdbml_file}")
            success_count += 1
        
        # Consider successful if we have both files
        if success_count >= 2:
            if cache:
                cache.cache_pdb_result(pdb_id, True, pdb_file)
            return True
        else:
            if cache:
                cache.cache_pdb_result(pdb_id, False)
            return False
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download structure files for {pdb_id}: {e}")
        if cache:
            cache.cache_pdb_result(pdb_id, False)
        return False

def download_3d_structures_for_gene(gene_name, structures_base_dir, organisms=None, max_structures=10, gene_aliases=None, cache=None, gene_cache=None, include_computed_models=True, prioritize_extracellular=True, extracellular_keywords=None, intracellular_keywords=None):
    """
    Download 3D structures and sequences for a specific bacterial gene (up to max_structures)
    
    Args:
        gene_name (str): Gene name to search for
        structures_base_dir (Path): Base directory for 3D structures
        organisms (list): Optional list of organism names to filter by (ignored for bacterial search)
        max_structures (int): Maximum number of structures to download (default: 10)
        gene_aliases (list): Alternative gene names to try if primary search fails
        cache (Protein3DStructureCache): Cache instance to use
        gene_cache (GeneDownloadCache): Gene-level cache to track completed downloads
        include_computed_models (bool): Include AlphaFold models if no experimental structures
    
    Returns:
        dict: Summary of downloads (sequences, structures, found)
    """
    
    # Create gene-specific directory in 3d_structures
    gene_dir = structures_base_dir / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    # Check gene cache first for completion status
    if gene_cache and gene_cache.is_gene_completed(gene_name):
        cached_structures = gene_cache.get_gene_structures(gene_name)
        logging.info(f"Gene {gene_name} already completed in cache with {len(cached_structures)} structures - SKIPPING")
        
        # Verify files still exist
        existing_pdb_files = list(gene_dir.glob("*.pdb.gz"))
        existing_pdbml_files = list(gene_dir.glob("*.xml.gz"))
        existing_fasta_files = list(gene_dir.glob("*.fasta"))
        
        if existing_pdb_files or existing_pdbml_files or existing_fasta_files:
            return {
                "sequences": len(existing_fasta_files),
                "structures": len(existing_pdb_files), 
                "found": True,
                "skipped": True,
                "pdb_ids": cached_structures
            }
        else:
            logging.warning(f"Gene {gene_name} marked as completed in cache but files missing - will retry download")
    
    # Fallback to file-based completion checking if cache doesn't have info
    existing_pdb_files = list(gene_dir.glob("*.pdb.gz"))
    existing_pdbml_files = list(gene_dir.glob("*.xml.gz"))
    existing_fasta_files = list(gene_dir.glob("*.fasta"))
    has_no_structures_marker = (gene_dir / "no_structures_found.txt").exists()
    
    # If directory has any structure files (PDB, PDBML, or FASTA), consider it complete
    if existing_pdb_files or existing_pdbml_files or existing_fasta_files:
        logging.info(f"Gene {gene_name} already has structures: {len(existing_pdb_files)} PDB files, {len(existing_pdbml_files)} PDBML files, {len(existing_fasta_files)} FASTA files - SKIPPING")
        
        # Update cache with existing structures
        if gene_cache:
            existing_structure_ids = [f.stem.replace('.pdb', '') for f in existing_pdb_files]
            structures_info = {
                "sequences": len(existing_fasta_files),
                "structures": len(existing_pdb_files),
                "structure_ids": existing_structure_ids,
                "experimental_count": len([sid for sid in existing_structure_ids if not sid.startswith('AF-')]),
                "computed_count": len([sid for sid in existing_structure_ids if sid.startswith('AF-')])
            }
            gene_cache.mark_gene_completed(gene_name, structures_info)
        
        return {
            "sequences": len(existing_fasta_files),
            "structures": len(existing_pdb_files), 
            "found": True,
            "skipped": True,
            "pdb_ids": [f.stem.replace('.pdb', '') for f in existing_pdb_files]
        }
    
    # Check if marked as no structures found (and directory is still empty)
    if has_no_structures_marker:
        logging.info(f"Gene {gene_name} already marked as no structures found - SKIPPING")
        return {"sequences": 0, "structures": 0, "found": False, "skipped": True, "completed": True, "pdb_ids": []}
    
    logging.info(f"Processing bacterial 3D structures for gene: {gene_name}")
    
    # Search for bacterial PDB structures via UniProt
    pdb_ids = search_pdb_via_uniprot_bacterial(gene_name, max_structures, gene_aliases, cache, prioritize_extracellular, extracellular_keywords, intracellular_keywords)
    alphafold_ids = []
    
    # If no experimental structures found and computed models are enabled, search AlphaFold
    if not pdb_ids and include_computed_models:
        logging.info(f"No experimental structures found for {gene_name}, searching AlphaFold models...")
        alphafold_ids = search_alphafold_models(gene_name, gene_aliases, cache)
    
    # Combine structure IDs
    all_structure_ids = pdb_ids + alphafold_ids
    
    if not all_structure_ids:
        logging.info(f"No 3D structures or computed models found for {gene_name}")
        # Create a marker file indicating no structures found
        no_structures_file = gene_dir / "no_structures_found.txt"
        with open(no_structures_file, 'w') as f:
            f.write(f"No 3D structures or computed models found for gene: {gene_name}\n")
            f.write(f"Search performed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Searched experimental structures: {'Yes' if max_structures > 0 else 'No'}\n")
            f.write(f"Searched computed models: {'Yes' if include_computed_models else 'No'}\n")
            if gene_aliases:
                f.write(f"Aliases tried: {', '.join(gene_aliases)}\n")
            else:
                f.write("No aliases available\n")
            f.write(f"Status: COMPLETED (no structures available)\n")
        
        return {"sequences": 0, "structures": 0, "found": False, "skipped": False, "completed": True, "structure_ids": [], "pdb_ids": [], "alphafold_ids": []}
    
    # Limit total structures to max_structures
    if len(all_structure_ids) > max_structures:
        # Prioritize experimental structures over computed models within the limit
        if len(pdb_ids) >= max_structures:
            all_structure_ids = pdb_ids[:max_structures]
            alphafold_ids = []
        else:
            remaining_slots = max_structures - len(pdb_ids)
            all_structure_ids = pdb_ids + alphafold_ids[:remaining_slots]
            alphafold_ids = alphafold_ids[:remaining_slots]
    
    logging.info(f"Processing {len(all_structure_ids)} structures for {gene_name}: {len(pdb_ids)} experimental + {len(alphafold_ids)} computed")
    logging.info(f"  PDB IDs: {pdb_ids}")
    logging.info(f"  AlphaFold IDs: {alphafold_ids}")
    
    total_sequences = 0
    total_structures = 0
    downloaded_pdb_ids = []
    
    # Process FASTA downloads for both experimental and computed structures
    structures_with_sequences = []
    
    # Process PDB structures
    for pdb_id in pdb_ids:
        sequences = get_pdb_fasta(pdb_id)
        if sequences:
            structures_with_sequences.append((pdb_id, sequences, 'experimental'))
            for i, seq_record in enumerate(sequences):
                filename = f"{pdb_id}_chain_{i+1}.fasta" if len(sequences) > 1 else f"{pdb_id}.fasta"
                seq_record.id = f"{pdb_id}_{seq_record.id}"
                seq_record.description = f"PDB:{pdb_id} {seq_record.description}"
                
                output_file = gene_dir / filename
                with open(output_file, 'w') as f:
                    SeqIO.write(seq_record, f, "fasta")
                
                total_sequences += 1
                logging.info(f"Saved experimental structure sequence: {output_file}")
        else:
            logging.warning(f"No sequences found for PDB {pdb_id}")
    
    # Process AlphaFold models
    for alphafold_id in alphafold_ids:
        sequences = get_alphafold_fasta(alphafold_id)
        if sequences:
            structures_with_sequences.append((alphafold_id, sequences, 'computed'))
            for i, seq_record in enumerate(sequences):
                filename = f"{alphafold_id}_chain_{i+1}.fasta" if len(sequences) > 1 else f"{alphafold_id}.fasta"
                
                output_file = gene_dir / filename
                with open(output_file, 'w') as f:
                    SeqIO.write(seq_record, f, "fasta")
                
                total_sequences += 1
                logging.info(f"Saved computed model sequence: {output_file}")
        else:
            logging.warning(f"No sequences found for AlphaFold model {alphafold_id}")
    
    # Download structures concurrently
    if structures_with_sequences:
        structures_to_download = [(struct_id, struct_type) for struct_id, _, struct_type in structures_with_sequences]
        logging.info(f"Downloading {len(structures_to_download)} structures concurrently...")
        
        # Use ThreadPoolExecutor for concurrent downloads (limited by max_structures or available structures)
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(max_structures, len(structures_to_download))) as executor:
            future_to_structure = {}
            
            # Submit download tasks for each structure type
            for struct_id, struct_type in structures_to_download:
                if struct_type == 'experimental':
                    download_func = partial(download_pdb_structure, output_dir=gene_dir, cache=cache)
                    future = executor.submit(download_func, struct_id)
                else:  # computed model
                    download_func = partial(download_alphafold_model, output_dir=gene_dir, cache=cache)
                    future = executor.submit(download_func, struct_id)
                
                future_to_structure[future] = (struct_id, struct_type)
            
            # Collect results
            for future in concurrent.futures.as_completed(future_to_structure):
                struct_id, struct_type = future_to_structure[future]
                try:
                    success = future.result()
                    if success:
                        total_structures += 1
                        downloaded_pdb_ids.append(struct_id)
                        structure_type_name = "experimental structure" if struct_type == 'experimental' else "computed model"
                        logging.info(f"✓ Downloaded {structure_type_name}: {struct_id}")
                    else:
                        structure_type_name = "experimental structure" if struct_type == 'experimental' else "computed model"
                        logging.warning(f"✗ Failed to download {structure_type_name}: {struct_id}")
                except Exception as e:
                    logging.error(f"Error downloading {struct_id}: {e}")
    
    # Mark gene as completed in cache if we successfully processed it
    if total_structures > 0 and gene_cache:
        structures_info = {
            "sequences": total_sequences,
            "structures": total_structures,
            "structure_ids": downloaded_pdb_ids,
            "experimental_count": len([sid for sid in downloaded_pdb_ids if not sid.startswith('AF-')]),
            "computed_count": len([sid for sid in downloaded_pdb_ids if sid.startswith('AF-')])
        }
        gene_cache.mark_gene_completed(gene_name, structures_info)
    
    # Create completion marker file for backward compatibility
    if total_structures > 0:
        completion_file = gene_dir / ".processing_complete.txt"
        with open(completion_file, 'w') as f:
            f.write(f"Gene: {gene_name}\n")
            f.write(f"Completed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Downloaded structures: {total_structures}\n")
            f.write(f"Downloaded sequences: {total_sequences}\n")
            f.write(f"Structure IDs: {', '.join(downloaded_pdb_ids)}\n")
            f.write(f"Experimental structures: {len([sid for sid in downloaded_pdb_ids if not sid.startswith('AF-')])}\n")
            f.write(f"Computed models: {len([sid for sid in downloaded_pdb_ids if sid.startswith('AF-')])}\n")
            f.write(f"Status: COMPLETED (structures downloaded)\n")
    
    # Separate experimental and computed structures for return
    downloaded_experimental = [sid for sid in downloaded_pdb_ids if not sid.startswith('AF-')]
    downloaded_computed = [sid for sid in downloaded_pdb_ids if sid.startswith('AF-')]
    
    return {
        "sequences": total_sequences, 
        "structures": total_structures, 
        "found": total_sequences > 0,
        "skipped": False,
        "completed": True,
        "structure_ids": downloaded_pdb_ids,
        "pdb_ids": downloaded_experimental,
        "alphafold_ids": downloaded_computed,
        "experimental_count": len(downloaded_experimental),
        "computed_count": len(downloaded_computed)
    }


def get_alphafold_metadata(alphafold_id):
    """
    Get metadata for an AlphaFold model
    
    Args:
        alphafold_id (str): AlphaFold identifier (e.g., AF-P12345)
    
    Returns:
        dict: Metadata dictionary
    """
    uniprot_accession = alphafold_id.replace('AF-', '')
    
    try:
        # Get AlphaFold metadata
        af_api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
        response = requests.get(af_api_url, timeout=30)
        response.raise_for_status()
        
        af_data = response.json()
        if not af_data:
            raise ValueError("No AlphaFold data found")
        
        af_entry = af_data[0] if isinstance(af_data, list) else af_data
        
        # Get additional UniProt metadata
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.json"
        uniprot_response = requests.get(uniprot_url, timeout=30)
        uniprot_data = {}
        if uniprot_response.status_code == 200:
            uniprot_data = uniprot_response.json()
        
        metadata = {
            "id": alphafold_id,
            "type": "computed_model",
            "database": "AlphaFold",
            "uniprot_accession": uniprot_accession,
            "gene_name": af_entry.get('gene', 'Unknown'),
            "organism": af_entry.get('organismScientificName', 'Unknown'),
            "model_confidence": af_entry.get('globalMetricValue', 'Unknown'),
            "model_version": af_entry.get('modelVersion', 'Unknown'),
            "sequence_length": af_entry.get('sequenceLength', 'Unknown'),
            "created_date": af_entry.get('modelCreatedDate', 'Unknown'),
            "protein_name": "Unknown",
            "method": "AlphaFold prediction",
            "resolution": "N/A (Computed Model)"
        }
        
        # Add protein name from UniProt if available
        if uniprot_data:
            protein_name = uniprot_data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
            metadata["protein_name"] = protein_name
        
        return metadata
        
    except Exception as e:
        logging.warning(f"Failed to get AlphaFold metadata for {alphafold_id}: {e}")
        return {
            "id": alphafold_id,
            "type": "computed_model",
            "database": "AlphaFold",
            "uniprot_accession": uniprot_accession,
            "gene_name": "Unknown",
            "organism": "Unknown",
            "model_confidence": "Unknown",
            "protein_name": "Unknown",
            "method": "AlphaFold prediction",
            "resolution": "N/A (Computed Model)",
            "error": str(e)
        }

def get_pdb_metadata(pdb_id):
    """Get metadata for a PDB structure from RCSB API"""
    try:
        api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        
        # Extract key metadata
        metadata = {
            "id": pdb_id,
            "type": "experimental_structure",
            "database": "PDB",
            "title": data.get("struct", {}).get("title", "Unknown"),
            "resolution": None,
            "method": "Unknown",
            "organism": "Unknown",
            "deposit_date": data.get("rcsb_accession_info", {}).get("deposit_date", "Unknown"),
            "authors": [],
            "journal": "Unknown"
        }
        
        # Get experimental method and resolution
        exptl_methods = data.get("exptl", [])
        if exptl_methods:
            metadata["method"] = exptl_methods[0].get("method", "Unknown")
        
        # Get resolution for X-ray structures
        refine_info = data.get("refine", [])
        if refine_info:
            resolution = refine_info[0].get("ls_d_res_high")
            if resolution:
                metadata["resolution"] = float(resolution)
        
        # Get organism information
        entity_src_gen = data.get("entity_src_gen", [])
        if entity_src_gen:
            organism = entity_src_gen[0].get("pdbx_gene_src_scientific_name")
            if organism:
                metadata["organism"] = organism
        
        # Get authors
        audit_author = data.get("audit_author", [])
        if audit_author:
            metadata["authors"] = [author.get("name", "") for author in audit_author[:5]]  # Top 5 authors
        
        # Get journal information
        citation_info = data.get("citation", [])
        if citation_info:
            citation = citation_info[0]
            journal = citation.get("journal_abbrev", "")
            if journal:
                metadata["journal"] = journal
        
        return metadata
        
    except Exception as e:
        logging.warning(f"Failed to get metadata for PDB {pdb_id}: {e}")
        return {
            "id": pdb_id,
            "type": "experimental_structure",
            "database": "PDB",
            "title": "Unknown",
            "resolution": None,
            "method": "Unknown",
            "organism": "Unknown",
            "deposit_date": "Unknown",
            "authors": [],
            "journal": "Unknown",
            "error": str(e)
        }

def load_protein_summary(summary_file):
    """Load protein summary with gene names and protein names"""
    try:
        df = pd.read_csv(summary_file, sep='\t')
        protein_info = {}
        for _, row in df.iterrows():
            gene_name = row['gene_name']
            protein_name = row.get('protein_name', '')
            protein_info[gene_name] = protein_name
        logging.info(f"Loaded protein information for {len(protein_info)} genes from {summary_file}")
        return protein_info
    except Exception as e:
        logging.warning(f"Could not load protein summary: {e}")
        return {}

def main():
    """Main function for Snakemake integration"""
    
    # Check if running under Snakemake
    try:
        # Get parameters from Snakemake
        proteins_file = Path(snakemake.input.protein_list)
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        max_structures = snakemake.params.max_structures  # Get from config
        include_computed_models = snakemake.params.include_computed_models  # Get from config
        prioritize_extracellular = snakemake.params.prioritize_extracellular  # Get from config
        extracellular_keywords = snakemake.params.extracellular_keywords  # Get from config
        intracellular_keywords = snakemake.params.intracellular_keywords  # Get from config
        
        # Create sentinel file path and derive structures directory from it
        sentinel_file = Path(snakemake.output.sentinel)
        structures_dir = sentinel_file.parent
        aliases_file = Path(f"data/quickgo/{paramset}/gene_aliases.json")
        
        # Look for summary file with protein names
        summary_file = Path(f"results/{analysis}_{paramset}/gene_selection/summary.tsv")
        
    except NameError:
        # Running standalone - use default values for testing
        analysis = "analysis1"
        paramset = "params1" 
        group = "positive"
        max_structures = 10  # Default value
        include_computed_models = True  # Default value
        prioritize_extracellular = True  # Default value
        extracellular_keywords = ["extracellular", "outer membrane", "cell surface", "secreted", "periplasm", "cell wall", "membrane", "surface"]
        intracellular_keywords = ["cytoplasm", "cytosol", "intracellular", "ribosome", "nucleus", "nucleoid"]
        proteins_file = Path(f"results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv")
        structures_dir = Path("data/protein_structures")
        aliases_file = Path(f"data/quickgo/{paramset}/gene_aliases.json")
        sentinel_file = Path(f"data/protein_structures/.{analysis}_{paramset}_{group}_structures_complete")
        summary_file = Path(f"results/{analysis}_{paramset}/gene_selection/summary.tsv")
    
    logging.info(f"Downloading bacterial 3D structures for {analysis}_{paramset} (unified Gram-independent)")
    logging.info(f"Proteins file: {proteins_file}")
    logging.info(f"Structures directory: {structures_dir}")
    logging.info(f"Summary file: {summary_file}")
    
    # Load gene aliases
    gene_aliases = load_gene_aliases(aliases_file)
    
    # Load protein summary with protein names
    protein_summary = load_protein_summary(summary_file)
    
    # Initialize caches
    cache = Protein3DStructureCache()
    gene_cache = GeneDownloadCache()
    logging.info(f"Structure cache initialized with {len(cache.cache)} entries from {CACHE_FILE}")
    logging.info(f"Gene download cache initialized with {len(gene_cache.cache)} entries from {GENE_DOWNLOAD_CACHE_FILE}")
    
    # Ensure output directories exist
    structures_dir.mkdir(parents=True, exist_ok=True)
    sentinel_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Process proteins and download bacterial 3D structures
    logging.info(f"Using max_structures = {max_structures}, include_computed_models = {include_computed_models}")
    
    # Process all unique genes from summary.tsv regardless of Gram classification
    def process_proteins_with_config(proteins_file, structures_base_dir, gene_aliases=None, cache=None):
        """Process proteins list with max_structures from config"""
        # Load all unique genes from summary.tsv instead of gram-specific file
        if not summary_file.exists():
            logging.error(f"Summary file not found: {summary_file}")
            return {}
        
        # Load summary.tsv and get unique genes
        summary_df = pd.read_csv(summary_file, sep='\t')
        unique_genes = summary_df['gene_name'].unique()
        total_genes = len(unique_genes)
        logging.info(f"Processing {total_genes} unique genes from {summary_file} (independent of Gram classification)")
        
        summary = {}
        genes_skipped = 0
        genes_processed = 0
        
        for gene_idx, gene_name in enumerate(unique_genes, 1):
            
            # Progress logging every 10 genes
            if gene_idx % 10 == 0 or gene_idx == 1:
                logging.info(f"Progress: {gene_idx}/{total_genes} genes ({gene_idx/total_genes*100:.1f}%) - Processed: {genes_processed}, Skipped: {genes_skipped}")
            
            # Get aliases for this gene
            aliases = gene_aliases.get(gene_name, []) if gene_aliases else []
            
            # Download bacterial 3D structures for this gene (without protein name)
            gene_summary = download_3d_structures_for_gene(
                gene_name, 
                structures_base_dir, 
                organisms=None,  # Ignored for bacterial search
                max_structures=max_structures,
                gene_aliases=aliases,
                cache=cache,
                gene_cache=gene_cache,
                include_computed_models=include_computed_models,
                prioritize_extracellular=prioritize_extracellular,
                extracellular_keywords=extracellular_keywords,
                intracellular_keywords=intracellular_keywords
            )
            
            # Collect metadata for all downloaded structures
            metadata_list = []
            
            # Get metadata for experimental structures
            for pdb_id in gene_summary.get('pdb_ids', []):
                metadata = get_pdb_metadata(pdb_id)
                metadata_list.append(metadata)
            
            # Get metadata for computed models
            for alphafold_id in gene_summary.get('alphafold_ids', []):
                metadata = get_alphafold_metadata(alphafold_id)
                metadata_list.append(metadata)
            
            if metadata_list:
                gene_summary['metadata'] = metadata_list
                
                # Save metadata to file
                gene_dir = structures_base_dir / gene_name
                metadata_file = gene_dir / "metadata.json"
                with open(metadata_file, 'w') as f:
                    json.dump(metadata_list, f, indent=2)
                logging.info(f"Saved metadata for {len(metadata_list)} structures ({gene_summary.get('experimental_count', 0)} experimental + {gene_summary.get('computed_count', 0)} computed): {metadata_file}")
            
            summary[gene_name] = gene_summary
            
            # Update counters
            if gene_summary.get('skipped', False):
                genes_skipped += 1
            else:
                genes_processed += 1
            
            # Minimal delay only every 5th gene for performance
            if gene_idx % 5 == 0:
                time.sleep(0.1)
        
        logging.info(f"Final processing summary: {genes_processed} processed, {genes_skipped} skipped out of {total_genes} total genes")
        return summary
    
    summary = process_proteins_with_config(proteins_file, structures_dir, gene_aliases, cache)
    
    # Calculate summary statistics
    total_sequences = sum(gene_data.get("sequences", 0) for gene_data in summary.values())
    total_structures = sum(gene_data.get("structures", 0) for gene_data in summary.values())
    total_experimental = sum(gene_data.get("experimental_count", 0) for gene_data in summary.values())
    total_computed = sum(gene_data.get("computed_count", 0) for gene_data in summary.values())
    genes_with_3d = len([g for g, data in summary.items() if data.get("found", False)])
    genes_with_experimental = len([g for g, data in summary.items() if data.get("experimental_count", 0) > 0])
    genes_with_computed = len([g for g, data in summary.items() if data.get("computed_count", 0) > 0])
    genes_skipped = len([g for g, data in summary.items() if data.get("skipped", False)])
    genes_completed = len([g for g, data in summary.items() if data.get("completed", False)])
    genes_processed = len(summary) - genes_skipped
    
    # Save summary
    summary_file = structures_dir / f"{analysis}_{paramset}_{group}_summary.json"
    summary_data = {
        "analysis": analysis,
        "paramset": paramset,
        "group": group,
        "total_genes": len(summary),
        "genes_processed": genes_processed,
        "genes_skipped": genes_skipped,
        "genes_completed": genes_completed,
        "total_3d_sequences": total_sequences,
        "total_3d_structures": total_structures,
        "total_experimental_structures": total_experimental,
        "total_computed_models": total_computed,
        "genes_with_3d": genes_with_3d,
        "genes_with_experimental": genes_with_experimental,
        "genes_with_computed": genes_with_computed,
        "max_structures_per_gene": max_structures,
        "computed_models_enabled": include_computed_models,
        "per_gene_summary": summary
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    # Save caches
    cache.save_cache()
    gene_cache.save_cache()
    
    # Create sentinel file to indicate completion
    sentinel_file.touch()
    
    logging.info(f"=== Final Summary ===")
    logging.info(f"Total genes: {len(summary)}")
    logging.info(f"  - Newly processed: {genes_processed}")
    logging.info(f"  - Skipped (already had structures): {genes_skipped}")
    logging.info(f"  - Total completed: {genes_completed}")
    logging.info(f"Downloaded {total_sequences} sequences and {total_structures} structures")
    logging.info(f"  - Experimental structures: {total_experimental}")
    logging.info(f"  - Computed models: {total_computed}")
    logging.info(f"Genes with any 3D data: {genes_with_3d}/{len(summary)} ({genes_with_3d/len(summary)*100:.1f}%)")
    logging.info(f"  - With experimental structures: {genes_with_experimental}")
    logging.info(f"  - With computed models: {genes_with_computed}")
    logging.info(f"Configuration: max_structures={max_structures}, include_computed_models={include_computed_models}")
    logging.info(f"Summary saved to: {summary_file}")
    logging.info(f"Sentinel file created: {sentinel_file}")
    
    # Report details on skipped genes if any
    if genes_skipped > 0:
        skipped_with_structures = [g for g, data in summary.items() if data.get("skipped", False) and data.get("found", False)]
        skipped_no_structures = [g for g, data in summary.items() if data.get("skipped", False) and not data.get("found", False)]
        if skipped_with_structures:
            logging.info(f"Skipped {len(skipped_with_structures)} genes that already had structures")
        if skipped_no_structures:
            logging.info(f"Skipped {len(skipped_no_structures)} genes previously marked as having no structures")

if __name__ == "__main__":
    main()