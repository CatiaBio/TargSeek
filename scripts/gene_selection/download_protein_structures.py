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
import re
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
STRUCTURE_DOWNLOAD_CACHE_FILE = CACHE_DIR / "structure_downloads_cache.json"

class StructureDownloadCache:
    """Manages caching of individual gene|structure|fasta combinations"""
    
    def __init__(self):
        self.cache = {}
        self.cache_updates = 0
        self.load_cache()
    
    def load_cache(self):
        """Load existing structure download cache"""
        if STRUCTURE_DOWNLOAD_CACHE_FILE.exists():
            try:
                with open(STRUCTURE_DOWNLOAD_CACHE_FILE, 'r') as f:
                    self.cache = json.load(f)
                logging.info(f"✓ Loaded structure download cache with {len(self.cache)} entries")
            except Exception as e:
                logging.warning(f"Could not load structure download cache: {e}")
                self.cache = {}
        else:
            logging.info("Creating new structure download cache")
            self.cache = {}
    
    def save_cache(self):
        """Save cache to disk if there were updates"""
        if self.cache_updates > 0:
            try:
                CACHE_DIR.mkdir(parents=True, exist_ok=True)
                with open(STRUCTURE_DOWNLOAD_CACHE_FILE, 'w') as f:
                    json.dump(self.cache, f, indent=2)
                logging.info(f"✓ Saved structure download cache with {self.cache_updates} new/updated entries")
                self.cache_updates = 0
            except Exception as e:
                logging.error(f"Failed to save structure download cache: {e}")
    
    def _make_key(self, gene_name: str, structure_id: str, fasta_type: str = "main") -> str:
        """Create cache key for gene|structure|fasta combination"""
        return f"{gene_name}|{structure_id}|{fasta_type}"
    
    def is_structure_downloaded(self, gene_name: str, structure_id: str, fasta_type: str = "main") -> bool:
        """Check if specific gene|structure|fasta combination is downloaded"""
        key = self._make_key(gene_name, structure_id, fasta_type)
        entry = self.cache.get(key, {})
        return entry.get('downloaded', False)
    
    def get_downloaded_structures_for_gene(self, gene_name: str) -> set:
        """Get all downloaded structure IDs for a gene"""
        structures = set()
        for key in self.cache:
            if key.startswith(f"{gene_name}|") and self.cache[key].get('downloaded', False):
                _, structure_id, _ = key.split('|', 2)
                structures.add(structure_id)
        return structures
    
    def mark_structure_downloaded(self, gene_name: str, structure_id: str, fasta_type: str = "main", 
                                structure_file_path: str = None, fasta_file_paths: list = None):
        """Mark specific gene|structure|fasta combination as downloaded"""
        key = self._make_key(gene_name, structure_id, fasta_type)
        self.cache[key] = {
            'downloaded': True,
            'download_date': datetime.now().isoformat(),
            'structure_file': structure_file_path,
            'fasta_files': fasta_file_paths or [],
            'gene': gene_name,
            'structure_id': structure_id,
            'fasta_type': fasta_type
        }
        self.cache_updates += 1
        logging.info(f"Marked {key} as downloaded in cache")
    
    def get_missing_structures_for_gene(self, gene_name: str, required_structures: set) -> set:
        """Get structure IDs that are missing for a gene"""
        downloaded = self.get_downloaded_structures_for_gene(gene_name)
        return required_structures - downloaded

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

def search_pdb_via_uniprot_bacterial(gene_name, max_structures=10, gene_aliases=None, cache=None, prioritize_extracellular=True, extracellular_keywords=None, intracellular_keywords=None, structure_filters=None, uniprot_protein_info=None):
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
        structure_filters (dict): Quality filtering criteria
        uniprot_protein_info (dict): UniProt protein information for keyword filtering
    """
    # Load UniProt protein info if not provided
    if uniprot_protein_info is None:
        uniprot_protein_info = load_uniprot_protein_info()
    
    # Check cache first
    if cache:
        cached_result = cache.get_cached_search(gene_name)
        if cached_result:
            logging.info(f"Found cached 3D structure search for {gene_name}: {len(cached_result['pdb_ids'])} structures")
            return cached_result['pdb_ids'][:max_structures]
    
    # Try primary gene name first
    pdb_ids = _search_uniprot_bacterial_single_gene(gene_name, max_structures, prioritize_extracellular, extracellular_keywords, intracellular_keywords, structure_filters, uniprot_protein_info)
    
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
                alias_pdb_ids = _search_uniprot_bacterial_single_gene(alias, max_structures, prioritize_extracellular, extracellular_keywords, intracellular_keywords, structure_filters, uniprot_protein_info)
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

def filter_structures_by_quality(pdb_ids, structure_filters=None):
    """
    Filter PDB structures based on quality criteria
    
    Args:
        pdb_ids (list): List of PDB identifiers
        structure_filters (dict): Quality filtering criteria
    
    Returns:
        list: Filtered list of PDB IDs that pass quality checks
    """
    if not structure_filters:
        return pdb_ids
    
    # Default filters if not provided
    default_filters = {
        "max_resolution": 3.5,
        "min_year": 2000,
        "exclude_obsolete": True,
        "preferred_methods": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", "SOLUTION NMR"],
        "min_sequence_length": 50,
        "max_sequence_length": 2000,
        "min_validation_confidence": 0.3
    }
    
    filters = {**default_filters, **structure_filters}
    filtered_pdb_ids = []
    
    logging.info(f"Applying quality filters to {len(pdb_ids)} structures:")
    logging.info(f"  - Max resolution: {filters['max_resolution']}Å")
    logging.info(f"  - Min year: {filters['min_year']}")
    logging.info(f"  - Exclude obsolete: {filters['exclude_obsolete']}")
    logging.info(f"  - Preferred methods: {filters['preferred_methods']}")
    
    for pdb_id in pdb_ids:
        try:
            metadata = get_pdb_metadata(pdb_id)
            
            # Check if structure is obsolete
            if filters['exclude_obsolete'] and metadata.get('is_obsolete', False):
                logging.debug(f"Excluding {pdb_id}: obsolete structure")
                continue
            
            # Check resolution
            resolution = metadata.get('resolution')
            if resolution and resolution > filters['max_resolution']:
                logging.debug(f"Excluding {pdb_id}: resolution {resolution}Å > {filters['max_resolution']}Å")
                continue
            
            # Check year
            year = metadata.get('year')
            if year and year < filters['min_year']:
                logging.debug(f"Excluding {pdb_id}: year {year} < {filters['min_year']}")
                continue
            
            # Check experimental method
            method = metadata.get('method', '').upper()
            if method != "UNKNOWN" and method not in filters['preferred_methods']:
                logging.debug(f"Excluding {pdb_id}: method '{method}' not in preferred methods")
                continue
            
            # Check sequence length
            seq_length = metadata.get('sequence_length')
            if seq_length:
                if seq_length < filters['min_sequence_length']:
                    logging.debug(f"Excluding {pdb_id}: sequence length {seq_length} < {filters['min_sequence_length']}")
                    continue
                if seq_length > filters['max_sequence_length']:
                    logging.debug(f"Excluding {pdb_id}: sequence length {seq_length} > {filters['max_sequence_length']}")
                    continue
            
            # Structure passed all filters
            filtered_pdb_ids.append(pdb_id)
            logging.debug(f"Including {pdb_id}: passed quality filters")
            
        except Exception as e:
            logging.warning(f"Error checking quality for {pdb_id}: {e}")
            # Include structure if we can't check quality (to avoid losing data)
            filtered_pdb_ids.append(pdb_id)
    
    logging.info(f"Quality filtering: {len(filtered_pdb_ids)}/{len(pdb_ids)} structures passed")
    return filtered_pdb_ids

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

def load_gene_annotations_from_data(analysis, paramset, base_dir=None):
    """
    Load gene annotations and GO terms from existing pipeline data
    Returns dict mapping gene_name -> expected keywords/functions
    """
    if base_dir is None:
        base_dir = Path(__file__).parent.parent.parent
    
    gene_annotations = {}
    
    try:
        # Load GO validation report with functional annotations
        go_report_path = Path(base_dir) / "data" / "quickgo" / paramset / "go_validation_report.tsv"
        if go_report_path.exists():
            go_df = pd.read_csv(go_report_path, sep='\t')
            for _, row in go_df.iterrows():
                gene = row.get('gene_symbol', '').lower()
                go_terms = row.get('go_terms', '')
                if gene and go_terms:
                    # Extract keywords from GO terms
                    keywords = [gene]  # Always include gene name
                    go_keywords = re.findall(r'\b\w+\b', go_terms.lower())
                    keywords.extend([kw for kw in go_keywords if len(kw) > 3])  # Filter short words
                    gene_annotations[gene] = keywords
        
        # Load protein information from UniProt data
        protein_info_path = Path(base_dir) / "data" / "uniprot" / "protein_info.json"
        if protein_info_path.exists():
            with open(protein_info_path, 'r') as f:
                protein_data = json.load(f)
            
            for gene, info in protein_data.items():
                gene_lower = gene.lower()
                if gene_lower not in gene_annotations:
                    gene_annotations[gene_lower] = []
                
                # Add gene name and aliases
                gene_annotations[gene_lower].append(gene_lower)
                aliases = info.get('aliases', [])
                gene_annotations[gene_lower].extend([alias.lower() for alias in aliases])
                
                # Add protein names and functions if available
                protein_names = info.get('protein_names', [])
                gene_annotations[gene_lower].extend([name.lower() for name in protein_names])
        
        # Load gene aliases
        gene_aliases_path = Path(base_dir) / "data" / "quickgo" / paramset / "gene_aliases.json"
        if gene_aliases_path.exists():
            with open(gene_aliases_path, 'r') as f:
                aliases_data = json.load(f)
            
            for gene, aliases in aliases_data.items():
                gene_lower = gene.lower()
                if gene_lower not in gene_annotations:
                    gene_annotations[gene_lower] = []
                gene_annotations[gene_lower].extend([alias.lower() for alias in aliases])
        
        logging.info(f"Loaded validation data for {len(gene_annotations)} genes from pipeline data")
        
    except Exception as e:
        logging.warning(f"Could not load gene annotations from pipeline data: {e}")
        # Fallback to minimal gene name matching
        return {}
    
    return gene_annotations

def load_uniprot_protein_info(base_dir=None):
    """
    Load protein information from data/uniprot/protein_info.json
    Returns dict mapping gene_name -> protein_name and keywords
    """
    if base_dir is None:
        base_dir = Path(__file__).parent.parent.parent
    
    protein_info = {}
    
    try:
        protein_info_path = Path(base_dir) / "data" / "uniprot" / "protein_info.json"
        if protein_info_path.exists():
            with open(protein_info_path, 'r') as f:
                uniprot_data = json.load(f)
            
            for gene, info in uniprot_data.items():
                gene_lower = gene.lower()
                protein_name = info.get('protein_name', '').strip()
                
                if protein_name:
                    # Extract keywords from protein name (split by spaces, remove short words)
                    protein_keywords = []
                    words = re.findall(r'\b\w+\b', protein_name.lower())
                    for word in words:
                        if len(word) >= 3 and word not in ['the', 'and', 'for', 'with', 'from']:
                            protein_keywords.append(word)
                    
                    protein_info[gene_lower] = {
                        'protein_name': protein_name,
                        'keywords': protein_keywords,
                        'gene_names': [name.lower() for name in info.get('gene_names', [])]
                    }
            
            logging.info(f"Loaded UniProt protein info for {len(protein_info)} genes")
            
    except Exception as e:
        logging.warning(f"Could not load UniProt protein info: {e}")
    
    return protein_info

def validate_protein_name_keywords(gene_name, structure_protein_names, uniprot_protein_info):
    """
    Validate if structure protein names contain keywords from the expected protein name in UniProt data
    
    Args:
        gene_name (str): Gene name to look up
        structure_protein_names (list): List of protein names from the structure
        uniprot_protein_info (dict): UniProt protein information loaded from JSON
    
    Returns:
        tuple: (passes_filter, match_score, matched_keywords)
    """
    gene_lower = gene_name.lower()
    
    # Get expected keywords from UniProt data
    if gene_lower not in uniprot_protein_info:
        logging.debug(f"No UniProt info found for gene {gene_name}, allowing structure")
        return True, 1.0, []  # Allow if no info available
    
    expected_keywords = uniprot_protein_info[gene_lower]['keywords']
    expected_protein_name = uniprot_protein_info[gene_lower]['protein_name']
    
    if not expected_keywords:
        logging.debug(f"No expected keywords for gene {gene_name}, allowing structure")
        return True, 1.0, []
    
    # Combine all structure protein names into one text
    combined_text = ' '.join(structure_protein_names).lower()
    
    # Count matched keywords
    matched_keywords = []
    for keyword in expected_keywords:
        if keyword in combined_text:
            matched_keywords.append(keyword)
    
    # Special handling for specific genes
    if gene_lower == 'cls':
        # For cardiolipin synthase, require "cardiolipin" or "synthase" keywords
        has_cardiolipin = 'cardiolipin' in combined_text
        has_synthase = 'synthase' in combined_text
        if has_cardiolipin or has_synthase:
            passes_filter = True
            match_score = 0.8 if (has_cardiolipin and has_synthase) else 0.6
        else:
            passes_filter = False
            match_score = 0.0
    else:
        # General case: require at least one keyword match
        passes_filter = len(matched_keywords) > 0
        match_score = len(matched_keywords) / len(expected_keywords) if expected_keywords else 0.0
    
    if passes_filter:
        logging.debug(f"✓ {gene_name}: Structure protein names contain expected keywords: {matched_keywords}")
        logging.debug(f"  Expected: {expected_protein_name}")
        logging.debug(f"  Found: {structure_protein_names}")
    else:
        logging.warning(f"✗ {gene_name}: Structure protein names missing expected keywords")
        logging.warning(f"  Expected keywords: {expected_keywords}")
        logging.warning(f"  Expected protein: {expected_protein_name}")
        logging.warning(f"  Found: {structure_protein_names}")
    
    return passes_filter, match_score, matched_keywords

def validate_protein_match(gene_name, protein_name, protein_function="", analysis="", paramset=""):
    """
    Validate if a protein matches the expected gene function using pipeline data
    Returns confidence score (0-1) and warnings
    """
    # Load expected functions from pipeline data
    gene_annotations = load_gene_annotations_from_data(analysis, paramset)
    
    expected_keywords = gene_annotations.get(gene_name.lower(), [])
    if not expected_keywords:
        # Fallback: use gene name itself
        expected_keywords = [gene_name.lower()]
    
    combined_text = f"{protein_name} {protein_function}".lower()
    warnings = []
    matches = 0
    
    # Check for expected keywords
    for keyword in expected_keywords:
        if keyword in combined_text:
            matches += 1
    
    # Special validation for known problematic patterns
    if gene_name.lower() == 'cls':
        if 'lsrb' in combined_text or 'lsr' in combined_text:
            return 0.1, ['LsrB protein found for CLS gene - likely wrong protein']
        if 'cardiolipin' not in combined_text and 'synthase' not in combined_text:
            warnings.append('CLS gene but no cardiolipin synthase keywords found')
    
    # Calculate confidence
    confidence = min(matches / max(len(expected_keywords), 1), 1.0) if expected_keywords else 0.5
    
    # Boost confidence if gene name appears in protein description
    if gene_name.lower() in combined_text:
        confidence = min(confidence + 0.3, 1.0)
    
    if confidence < 0.3:
        warnings.append(f'Low confidence match for {gene_name}: {matches}/{len(expected_keywords)} keywords found')
    
    return confidence, warnings

def _search_uniprot_bacterial_single_gene(gene_name, max_structures, prioritize_extracellular=True, extracellular_keywords=None, intracellular_keywords=None, structure_filters=None, uniprot_protein_info=None):
    """Helper function to search UniProt for bacterial proteins with structures for a single gene name"""
    
    # Default keywords if not provided
    if extracellular_keywords is None:
        extracellular_keywords = ["extracellular", "outer membrane", "cell surface", "secreted", "periplasm", "cell wall", "membrane", "surface"]
    if intracellular_keywords is None:
        intracellular_keywords = ["cytoplasm", "cytosol", "intracellular", "ribosome", "nucleus", "nucleoid"]
    
    # Load UniProt protein info if not provided
    if uniprot_protein_info is None:
        uniprot_protein_info = load_uniprot_protein_info()
    
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
                
                # Validate protein match
                protein_name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')
                protein_function = entry.get('comments', [{}])[0].get('texts', [{}])[0].get('value', '') if entry.get('comments') else ''
                # Note: analysis and paramset would need to be passed from the calling function
                confidence, warnings = validate_protein_match(gene_name, protein_name, protein_function)
                
                # Validate protein name keywords from UniProt data
                structure_protein_names = [protein_name] if protein_name else []
                # Also get alternative names if available
                alt_names = entry.get('proteinDescription', {}).get('alternativeName', [])
                if isinstance(alt_names, list):
                    for alt_name in alt_names:
                        if isinstance(alt_name, dict) and 'fullName' in alt_name:
                            alt_name_value = alt_name['fullName'].get('value', '')
                            if alt_name_value:
                                structure_protein_names.append(alt_name_value)
                elif isinstance(alt_names, dict) and 'fullName' in alt_names:
                    alt_name_value = alt_names['fullName'].get('value', '')
                    if alt_name_value:
                        structure_protein_names.append(alt_name_value)
                
                passes_keyword_filter, keyword_match_score, matched_keywords = validate_protein_name_keywords(
                    gene_name, structure_protein_names, uniprot_protein_info
                )
                
                # Log validation warnings
                if warnings:
                    for warning in warnings:
                        logging.warning(f"  {gene_name}/{pdb_refs[0]}: {warning}")
                
                entry_info = {
                    'entry': entry,
                    'pdb_ids': pdb_refs,
                    'organism': organism,
                    'is_bacterial': is_bacterial,
                    'localization_score': localization_score,
                    'validation_confidence': confidence,
                    'validation_warnings': warnings,
                    'passes_keyword_filter': passes_keyword_filter,
                    'keyword_match_score': keyword_match_score,
                    'matched_keywords': matched_keywords,
                    'structure_protein_names': structure_protein_names
                }
                
                if is_bacterial:
                    bacterial_entries.append(entry_info)
                    logging.debug(f"Found bacterial protein from {organism} with PDB: {pdb_refs[0]} (score: {localization_score})")
                else:
                    other_entries.append(entry_info)
                    logging.debug(f"Found non-bacterial protein from {organism}")
        
        # Sort and select bacterial structures
        if bacterial_entries:
            # Filter out very low confidence matches (< 0.3) for critical genes
            critical_genes = ['cls', 'cls2']  # Genes where accuracy is crucial
            if gene_name.lower() in critical_genes:
                high_confidence_entries = [e for e in bacterial_entries if e['validation_confidence'] >= 0.3]
                if high_confidence_entries:
                    bacterial_entries = high_confidence_entries
                    logging.info(f"Filtered to {len(bacterial_entries)} high-confidence matches for critical gene {gene_name}")
                else:
                    logging.warning(f"No high-confidence matches found for critical gene {gene_name}, using all matches")
            
            # Apply protein name keyword filtering
            original_count = len(bacterial_entries)
            keyword_filtered_entries = [e for e in bacterial_entries if e['passes_keyword_filter']]
            
            if keyword_filtered_entries:
                bacterial_entries = keyword_filtered_entries
                logging.info(f"Protein name keyword filtering: {len(bacterial_entries)}/{original_count} structures passed for {gene_name}")
                
                # Log filtered entries for debugging
                for entry_info in bacterial_entries:
                    logging.debug(f"  ✓ Kept: {entry_info['pdb_ids'][0]} - {entry_info['structure_protein_names']} (matched: {entry_info['matched_keywords']})")
                
                # Log rejected entries
                rejected_entries = [e for e in bacterial_entries if not e['passes_keyword_filter']]
                for entry_info in rejected_entries:
                    logging.debug(f"  ✗ Rejected: {entry_info['pdb_ids'][0]} - {entry_info['structure_protein_names']} (no keyword match)")
            else:
                logging.warning(f"No structures passed protein name keyword filter for {gene_name}, keeping all {original_count} structures")
                # Keep all entries if none pass the filter to avoid losing all data
                # This is a fallback to prevent empty results due to overly strict filtering
            
            # Sort by keyword match score first, then validation confidence, then localization score
            bacterial_entries.sort(key=lambda x: (-x['keyword_match_score'], -x['validation_confidence'], -x['localization_score']))
            
            # Extract all PDB IDs first
            all_bacterial_pdb_ids = []
            for entry_info in bacterial_entries:
                for pdb_id in entry_info['pdb_ids']:
                    if pdb_id not in all_bacterial_pdb_ids:
                        all_bacterial_pdb_ids.append(pdb_id)
            
            # Apply quality filtering first if enabled
            if structure_filters:
                logging.info(f"Applying quality filters to {len(all_bacterial_pdb_ids)} bacterial structures for {gene_name}")
                all_bacterial_pdb_ids = filter_structures_by_quality(all_bacterial_pdb_ids, structure_filters)
                if not all_bacterial_pdb_ids:
                    logging.info(f"No bacterial structures passed quality filters for {gene_name}")
                    return []
            
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
    Download PDB structure files with mmCIF fallback for large structures
    
    Args:
        pdb_id (str): PDB identifier
        output_dir (Path): Output directory to save the structure
        cache (Protein3DStructureCache): Cache instance to use
    
    Returns:
        bool: True if download succeeded, False otherwise
    """
    
    # Define file paths for both formats
    pdb_file = output_dir / f"{pdb_id}.pdb.gz"
    mmcif_file = output_dir / f"{pdb_id}.cif.gz"
    
    # Check if either file already exists
    pdb_exists = pdb_file.exists() and pdb_file.stat().st_size > 0
    mmcif_exists = mmcif_file.exists() and mmcif_file.stat().st_size > 0
    
    if pdb_exists or mmcif_exists:
        # File exists and is not empty, no need to re-download
        existing_file = pdb_file if pdb_exists else mmcif_file
        if cache:
            cache.cache_pdb_result(pdb_id, True, existing_file)
        return True
    
    # Check cache for download status only if no files exist
    if cache:
        cached_pdb = cache.get_cached_pdb(pdb_id)
        if cached_pdb and cached_pdb.get('success'):
            cached_path = Path(cached_pdb.get('file_path', ''))
            if cached_path.exists() and cached_path.stat().st_size > 0:
                return True
    
    # Try PDB format first
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    
    try:
        # Download PDB format
        response = requests.get(pdb_url, timeout=30)
        response.raise_for_status()
        
        with open(pdb_file, 'wb') as f:
            f.write(response.content)
        
        logging.info(f"Downloaded PDB structure: {pdb_file}")
        
        # Cache success and return
        if cache:
            cache.cache_pdb_result(pdb_id, True, pdb_file)
        return True
        
    except requests.exceptions.HTTPError as e:
        # If 404 error, try mmCIF format as fallback
        logging.debug(f"HTTPError caught for {pdb_id}: {e}")
        try:
            status_code = e.response.status_code
            logging.debug(f"Status code: {status_code}")
            if status_code == 404:
                logging.info(f"PDB format not available for {pdb_id} (404), trying mmCIF format...")
                return download_mmcif_structure(pdb_id, output_dir, cache)
        except (AttributeError, TypeError):
            logging.debug(f"Could not get status code from HTTPError for {pdb_id}")
        
        logging.warning(f"Failed to download PDB structure for {pdb_id}: {e}")
        if cache:
            cache.cache_pdb_result(pdb_id, False)
        return False
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download PDB structure for {pdb_id}: {e}")
        if cache:
            cache.cache_pdb_result(pdb_id, False)
        return False

def download_mmcif_structure(pdb_id, output_dir, cache=None):
    """
    Download mmCIF format structure file
    
    Args:
        pdb_id (str): PDB identifier
        output_dir (Path): Output directory to save the structure
        cache (Protein3DStructureCache): Cache instance to use
    
    Returns:
        bool: True if download succeeded, False otherwise
    """
    
    mmcif_file = output_dir / f"{pdb_id}.cif.gz"
    
    try:
        # Download mmCIF format
        mmcif_url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
        response = requests.get(mmcif_url, timeout=30)
        response.raise_for_status()
        
        with open(mmcif_file, 'wb') as f:
            f.write(response.content)
        
        logging.info(f"Downloaded mmCIF structure: {mmcif_file}")
        
        # Cache success and return
        if cache:
            cache.cache_pdb_result(pdb_id, True, mmcif_file)
        return True
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download mmCIF format for {pdb_id}: {e}")
        if cache:
            cache.cache_pdb_result(pdb_id, False)
        return False

def download_3d_structures_for_gene(gene_name, structures_base_dir, shared_structures_dir, organisms=None, max_structures=10, gene_aliases=None, cache=None, structure_cache=None, include_computed_models=True, prioritize_extracellular=True, extracellular_keywords=None, intracellular_keywords=None, structure_filters=None):
    """
    Download 3D structures and sequences for a specific bacterial gene (up to max_structures)
    
    Args:
        gene_name (str): Gene name to search for
        structures_base_dir (Path): Base directory for gene-specific FASTA files
        shared_structures_dir (Path): Shared directory for 3D structure files
        organisms (list): Optional list of organism names to filter by (ignored for bacterial search)
        max_structures (int): Maximum number of structures to download (default: 10)
        gene_aliases (list): Alternative gene names to try if primary search fails
        cache (Protein3DStructureCache): Cache instance to use
        structure_cache (StructureDownloadCache): Structure-level cache to track individual downloads
        include_computed_models (bool): Include AlphaFold models if no experimental structures
    
    Returns:
        dict: Summary of downloads (sequences, structures, found, fasta_structure_mapping)
    """
    
    # Create gene-specific directory for FASTA files
    gene_dir = structures_base_dir / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    # Create shared directory for 3D structure files
    shared_structures_dir.mkdir(parents=True, exist_ok=True)
    
    # Check existing files to identify what FASTA files need structure files
    existing_fasta_files = list(gene_dir.glob("*.fasta"))
    has_no_structures_marker = (gene_dir / "no_structures_found.txt").exists()
    
    # Check existing structures in shared directory
    existing_pdb_files = list(shared_structures_dir.glob("*.pdb.gz"))
    existing_mmcif_files = list(shared_structures_dir.glob("*.cif.gz"))
    existing_structure_files = existing_pdb_files + existing_mmcif_files
    
    # Extract structure IDs from existing FASTA files - these are our requirements
    required_structure_ids = set()
    for fasta_file in existing_fasta_files:
        filename = fasta_file.stem
        if '_chain_' in filename:
            structure_id = filename.split('_chain_')[0]
        elif '_' in filename and filename.split('_')[-1].isdigit():
            structure_id = '_'.join(filename.split('_')[:-1])
        else:
            structure_id = filename
        required_structure_ids.add(structure_id)
    
    # If no FASTA files exist, we need to search for new structures
    if not existing_fasta_files and not has_no_structures_marker:
        logging.info(f"No existing FASTA files for {gene_name}, will search for new structures")
        search_needed = True
        structures_to_download = set()
    else:
        # Check which structures are missing using the new cache system
        if structure_cache:
            missing_structures = structure_cache.get_missing_structures_for_gene(gene_name, required_structure_ids)
        else:
            # Fallback: check filesystem directly
            existing_structure_ids = set()
            for struct_file in existing_structure_files:
                filename = struct_file.name
                if filename.endswith('.pdb.gz'):
                    structure_id = filename[:-7]
                elif filename.endswith('.cif.gz'):
                    structure_id = filename[:-7]
                elif filename.endswith('.pdb'):
                    structure_id = filename[:-4]
                elif filename.endswith('.cif'):
                    structure_id = filename[:-4]
                else:
                    structure_id = filename.split('.')[0]
                existing_structure_ids.add(structure_id)
            missing_structures = required_structure_ids - existing_structure_ids
        
        if not missing_structures:
            logging.info(f"Gene {gene_name} already has all required structures: {len(existing_structure_files)} files for {len(required_structure_ids)} structure IDs - SKIPPING")
            return {
                "sequences": len(existing_fasta_files),
                "structures": len(existing_structure_files), 
                "found": True,
                "skipped": True,
                "pdb_ids": list(required_structure_ids)
            }
        else:
            logging.info(f"Gene {gene_name} missing structures: {missing_structures}")
            search_needed = False
            structures_to_download = missing_structures
    
    # Check if marked as no structures found (only if we need new structures)
    if has_no_structures_marker and search_needed:
        logging.info(f"Gene {gene_name} already marked as no structures found - SKIPPING")
        return {"sequences": 0, "structures": 0, "found": False, "skipped": True, "completed": True, "pdb_ids": []}
    
    logging.info(f"Processing bacterial 3D structures for gene: {gene_name}")
    
    # Determine which structures we need to find/download
    if search_needed:
        # Search for new bacterial PDB structures via UniProt
        # Load UniProt protein info for filtering
        uniprot_protein_info = load_uniprot_protein_info()
        pdb_ids = search_pdb_via_uniprot_bacterial(gene_name, max_structures, gene_aliases, cache, prioritize_extracellular, extracellular_keywords, intracellular_keywords, structure_filters, uniprot_protein_info)
        alphafold_ids = []
        
        # If no experimental structures found and computed models are enabled, search AlphaFold
        if not pdb_ids and include_computed_models:
            logging.info(f"No experimental structures found for {gene_name}, searching AlphaFold models...")
            alphafold_ids = search_alphafold_models(gene_name, gene_aliases, cache)
        
        # Combine structure IDs
        all_structure_ids = pdb_ids + alphafold_ids
        structures_to_download = set(all_structure_ids)
    else:
        # We only need to download the missing structures
        pdb_ids = []
        alphafold_ids = []
        all_structure_ids = list(structures_to_download)
        logging.info(f"Downloading missing structures: {structures_to_download}")
        
        # Split existing structure_to_download into PDB and AlphaFold based on naming
        for structure_id in structures_to_download:
            if structure_id.startswith('AF-'):
                alphafold_ids.append(structure_id)
            else:
                pdb_ids.append(structure_id)
    
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
    
    # Download structures concurrently (only missing ones)
    if structures_with_sequences:
        # Filter structures to download only missing ones
        structures_to_download = []
        for struct_id, _, struct_type in structures_with_sequences:
            # Check if structure file already exists in shared directory
            pdb_file = shared_structures_dir / f"{struct_id}.pdb.gz"
            mmcif_file = shared_structures_dir / f"{struct_id}.cif.gz"
            pdb_exists = pdb_file.exists() and pdb_file.stat().st_size > 0
            mmcif_exists = mmcif_file.exists() and mmcif_file.stat().st_size > 0
            
            if not (pdb_exists or mmcif_exists):
                structures_to_download.append((struct_id, struct_type))
                logging.info(f"Will download missing structure: {struct_id}")
            else:
                logging.info(f"Structure {struct_id} already exists in shared directory - skipping download")
    
    # Convert structures_to_download to the format expected by download logic
    download_candidates = []
    if not search_needed and structures_to_download:
        # These are missing structures we need to download
        logging.info(f"Preparing to download missing structures: {structures_to_download}")
        download_candidates = [(sid, 'experimental') for sid in structures_to_download]
    
    if structures_with_sequences or download_candidates:
        # Combine both types of structures to download
        all_downloads = structures_to_download + download_candidates
        
        if all_downloads:
            logging.info(f"Downloading {len(all_downloads)} structures concurrently...")
            
            # Use ThreadPoolExecutor for concurrent downloads
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(max_structures, len(all_downloads))) as executor:
                future_to_structure = {}
                
                # Submit download tasks for each structure type
                for struct_id, struct_type in all_downloads:
                    if struct_type == 'experimental':
                        download_func = partial(download_pdb_structure, output_dir=shared_structures_dir, cache=cache)
                        future = executor.submit(download_func, struct_id)
                    else:  # computed model
                        download_func = partial(download_alphafold_model, output_dir=shared_structures_dir, cache=cache)
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
                        
                        # Mark structure as downloaded in new cache system
                        if structure_cache:
                            # Find structure file path in shared directory
                            pdb_file = shared_structures_dir / f"{struct_id}.pdb.gz"
                            mmcif_file = shared_structures_dir / f"{struct_id}.cif.gz"
                            structure_file_path = pdb_file if pdb_file.exists() else mmcif_file
                            
                            # Find FASTA files for this structure in gene directory
                            fasta_files = list(gene_dir.glob(f"{struct_id}*.fasta"))
                            fasta_file_paths = [str(f) for f in fasta_files]
                            
                            structure_cache.mark_structure_downloaded(
                                gene_name, struct_id, "main", 
                                str(structure_file_path), fasta_file_paths
                            )
                    else:
                        structure_type_name = "experimental structure" if struct_type == 'experimental' else "computed model"
                        logging.warning(f"✗ Failed to download {structure_type_name}: {struct_id}")
                except Exception as e:
                    logging.error(f"Error downloading {struct_id}: {e}")
    
    # Save any structures that already existed but weren't tracked in cache
    if structure_cache:
        # Check for any existing structures not yet in cache in shared directory
        existing_structure_files = list(shared_structures_dir.glob("*.pdb.gz")) + list(shared_structures_dir.glob("*.cif.gz"))
        for struct_file in existing_structure_files:
            # Extract structure ID
            filename = struct_file.name
            if filename.endswith('.pdb.gz'):
                structure_id = filename[:-7]
            elif filename.endswith('.cif.gz'):
                structure_id = filename[:-7]
            else:
                structure_id = filename.split('.')[0]
            
            # Check if already tracked for this gene
            if not structure_cache.is_structure_downloaded(gene_name, structure_id):
                # Find FASTA files for this structure in gene directory
                fasta_files = list(gene_dir.glob(f"{structure_id}*.fasta"))
                fasta_file_paths = [str(f) for f in fasta_files]
                
                # Only mark if this gene has FASTA files for this structure
                if fasta_file_paths:
                    structure_cache.mark_structure_downloaded(
                        gene_name, structure_id, "main", 
                        str(struct_file), fasta_file_paths
                    )
    
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
            "year": None,
            "sequence_length": None,
            "is_obsolete": False,
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
        
        # Extract year from deposit date
        deposit_date = metadata.get("deposit_date", "")
        if deposit_date and deposit_date != "Unknown":
            try:
                year = int(deposit_date.split('-')[0])
                metadata["year"] = year
            except (ValueError, IndexError):
                pass
        
        # Get sequence length
        try:
            entity_info = data.get("entity", [])
            if entity_info:
                seq_length = entity_info[0].get("pdbx_number_of_molecules")
                if seq_length:
                    metadata["sequence_length"] = int(seq_length)
        except (ValueError, IndexError, TypeError):
            pass
        
        # Check if structure is obsolete
        pdbx_database_status = data.get("pdbx_database_status", {})
        status_code = pdbx_database_status.get("status_code", "")
        metadata["is_obsolete"] = status_code.upper() in ["OBS", "WDRN"]
        
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
            "year": None,
            "sequence_length": None,
            "is_obsolete": False,
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
        structure_filters = snakemake.params.structure_filters  # Get from config
        
        # Get output files from Snakemake
        sentinel_file = Path(snakemake.output.sentinel)
        mapping_file = Path(snakemake.output.mapping_file)
        structures_dir = sentinel_file.parent
        shared_structures_dir = Path(snakemake.config["paths"]["protein_structures"]["structures_dir"])
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
        structure_filters = {}  # Default empty filters
        proteins_file = Path(f"results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv")
        structures_dir = Path("data/protein_structures")
        shared_structures_dir = Path("data/protein_structures/pdb_files")
        aliases_file = Path(f"data/quickgo/{paramset}/gene_aliases.json")
        sentinel_file = Path(f"data/protein_structures/.{analysis}_{paramset}_{group}_structures_complete")
        mapping_file = Path(f"data/protein_structures/{analysis}_{paramset}_fasta_structure_mapping.tsv")
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
    structure_cache = StructureDownloadCache()
    logging.info(f"Structure cache initialized with {len(cache.cache)} entries from {CACHE_FILE}")
    logging.info(f"Structure download cache initialized with {len(structure_cache.cache)} entries from {STRUCTURE_DOWNLOAD_CACHE_FILE}")
    
    # Ensure output directories exist
    structures_dir.mkdir(parents=True, exist_ok=True)
    sentinel_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Process proteins and download bacterial 3D structures
    logging.info(f"Using max_structures = {max_structures}, include_computed_models = {include_computed_models}")
    
    # Process all unique genes from summary.tsv regardless of Gram classification
    def process_proteins_with_config(proteins_file, structures_base_dir, shared_structures_dir, gene_aliases=None, cache=None):
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
                shared_structures_dir,
                organisms=None,  # Ignored for bacterial search
                max_structures=max_structures,
                gene_aliases=aliases,
                cache=cache,
                structure_cache=structure_cache,
                include_computed_models=include_computed_models,
                prioritize_extracellular=prioritize_extracellular,
                extracellular_keywords=extracellular_keywords,
                intracellular_keywords=intracellular_keywords,
                structure_filters=structure_filters
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
    
    summary = process_proteins_with_config(proteins_file, structures_dir, shared_structures_dir, gene_aliases, cache)
    
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
    
    # Create FASTA to structure mapping file
    logging.info("Creating FASTA to structure mapping file...")
    mapping_entries = []
    
    # Iterate through all genes and their FASTA files
    for gene_name in summary.keys():
        gene_dir = structures_dir / gene_name
        if gene_dir.exists():
            fasta_files = list(gene_dir.glob("*.fasta"))
            for fasta_file in fasta_files:
                # Extract structure ID from FASTA filename
                filename = fasta_file.stem
                if '_chain_' in filename:
                    structure_id = filename.split('_chain_')[0]
                elif '_' in filename and filename.split('_')[-1].isdigit():
                    structure_id = '_'.join(filename.split('_')[:-1])
                else:
                    structure_id = filename
                
                # Look for corresponding structure file in shared directory
                structure_file = None
                pdb_file = shared_structures_dir / f"{structure_id}.pdb.gz"
                mmcif_file = shared_structures_dir / f"{structure_id}.cif.gz"
                alphafold_file = shared_structures_dir / f"{structure_id}.pdb"  # AlphaFold files are uncompressed
                
                if pdb_file.exists():
                    structure_file = pdb_file
                elif mmcif_file.exists():
                    structure_file = mmcif_file
                elif alphafold_file.exists():
                    structure_file = alphafold_file
                
                if structure_file:
                    mapping_entries.append({
                        'fasta_path': str(fasta_file),
                        'structure_path': str(structure_file),
                        'gene_name': gene_name,
                        'structure_id': structure_id,
                        'structure_type': 'computed' if structure_id.startswith('AF-') else 'experimental'
                    })
    
    # Write mapping file
    if mapping_entries:
        mapping_df = pd.DataFrame(mapping_entries)
        mapping_df.to_csv(mapping_file, sep='\t', index=False)
        logging.info(f"Created mapping file with {len(mapping_entries)} entries: {mapping_file}")
    else:
        # Create empty mapping file
        with open(mapping_file, 'w') as f:
            f.write("fasta_path\tstructure_path\tgene_name\tstructure_id\tstructure_type\n")
        logging.info(f"Created empty mapping file: {mapping_file}")
    
    # Save caches
    cache.save_cache()
    structure_cache.save_cache()
    
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