#!/usr/bin/env python3
"""
Simple Sequence Similarity-Based Structure Downloader
====================================================

Downloads protein structures based on sequence similarity to reference sequences.
Uses a reference sequence from standard bacteria (E. coli, Bacillus, Salmonella, Pseudomonas)
to find similar structures in the PDB via BLAST search.

Usage:
    python download_structures_by_similarity.py <gene_name>

Example:
    python download_structures_by_similarity.py bamA
"""

import requests
import json
import time
import logging
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import concurrent.futures
from functools import partial

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def fix_corrupted_alphafold_identifier(identifier):
    """
    Fix corrupted AlphaFold identifiers from PDB API
    
    Converts: AF_AFQ325W3F1 -> AF-Q325W3-F1
    Pattern: AF_AF{UNIPROT_ID}F1 -> AF-{UNIPROT_ID}-F1
    
    Args:
        identifier (str): Corrupted identifier like AF_AFQ325W3F1
        
    Returns:
        str: Fixed AlphaFold URL format like AF-Q325W3-F1
    """
    # Pattern: AF_AF followed by UniProt ID followed by F1
    pattern = r'^AF_AF([A-Za-z0-9_]+)F1$'
    match = re.match(pattern, identifier)
    
    if match:
        uniprot_id = match.group(1)
        fixed_id = f"AF-{uniprot_id}-F1"
        logging.debug(f"Fixed corrupted AlphaFold ID: {identifier} -> {fixed_id}")
        return fixed_id
    
    # If no match, return as-is
    logging.debug(f"No fix needed for identifier: {identifier}")
    return identifier


class SequenceSimilarityStructureDownloader:
    """Downloads structures based on sequence similarity"""
    
    def __init__(self, data_dir=None, output_dir=None, structures_dir=None):
        if data_dir is None:
            data_dir = Path(__file__).parent.parent.parent / "data"
        
        self.data_dir = Path(data_dir)
        self.sequences_dir = self.data_dir / "protein_sequences"
        
        # Main protein_structures directory
        self.protein_structures_base = Path(__file__).parent.parent.parent / "data" / "protein_structures"
        
        # Gene-specific FASTA directories under protein_structures/
        self.output_dir = self.protein_structures_base
        
        # All 3D structures go in protein_structures/pdb_files/
        self.structures_dir = self.protein_structures_base / "pdb_files"
        
        # Initialize caches
        self.cache_dir = Path(__file__).parent.parent.parent / "cache" / "protein_structures"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.structure_cache_file = self.cache_dir / "downloaded_structures.json"
        self.fasta_cache_file = self.cache_dir / "downloaded_fastas.json"
        
        # Load existing caches
        self.structure_cache = self._load_cache(self.structure_cache_file)
        self.fasta_cache = self._load_cache(self.fasta_cache_file)
        
        # Preferred reference organisms (in order of preference)
        # Include standard bacteria used for protein expression and research
        self.reference_organisms = [
            "Escherichia_coli",           # Most common expression host
            "Bacillus_subtilis",          # Gram-positive expression host
            "Salmonella_enterica",        # Enterobacteria family
            "Pseudomonas_aeruginosa",     # Gram-negative research model
            "Staphylococcus_aureus",      # Important Gram-positive pathogen
            "Streptococcus_pneumoniae",   # Gram-positive pathogen
            "Klebsiella_pneumoniae",      # Enterobacteria family
            "Enterococcus_faecalis",      # Gram-positive research model
            "Bacillus_cereus",            # Bacillus species
            "Pseudomonas_putida",         # Environmental Pseudomonas
            "Acinetobacter_baumannii",    # Gram-negative pathogen
            "Listeria_monocytogenes",     # Gram-positive pathogen
            "Streptococcus_pyogenes",     # Group A Streptococcus
            "Vibrio_cholerae",            # Marine pathogen
            "Campylobacter_jejuni"        # Foodborne pathogen
        ]
    
    def _load_cache(self, cache_file):
        """Load cache from JSON file"""
        try:
            if cache_file.exists():
                with open(cache_file, 'r') as f:
                    return json.load(f)
        except Exception as e:
            logging.debug(f"Error loading cache {cache_file}: {e}")
        return {}
    
    def _save_cache(self, cache_file, cache_data):
        """Save cache to JSON file"""
        try:
            with open(cache_file, 'w') as f:
                json.dump(cache_data, f, indent=2)
        except Exception as e:
            logging.debug(f"Error saving cache {cache_file}: {e}")
    
    def _is_structure_cached(self, structure_id, structure_type='experimental'):
        """Check if structure is already downloaded and cached"""
        cache_key = f"{structure_id}_{structure_type}"
        if cache_key in self.structure_cache:
            # Check if file actually exists
            cached_path = self.structure_cache[cache_key]
            if Path(cached_path).exists():
                return True
            else:
                # Remove invalid cache entry
                del self.structure_cache[cache_key]
                self._save_cache(self.structure_cache_file, self.structure_cache)
        return False
    
    def _cache_structure(self, structure_id, file_path, structure_type='experimental'):
        """Add structure to cache"""
        cache_key = f"{structure_id}_{structure_type}"
        self.structure_cache[cache_key] = str(file_path)
        self._save_cache(self.structure_cache_file, self.structure_cache)
    
    def _is_fasta_cached(self, gene_name, structure_id, structure_type='experimental'):
        """Check if FASTA is already downloaded and cached"""
        cache_key = f"{gene_name}_{structure_id}_{structure_type}"
        if cache_key in self.fasta_cache:
            # Check if file actually exists
            cached_path = self.fasta_cache[cache_key]
            if Path(cached_path).exists():
                return True
            else:
                # Remove invalid cache entry
                del self.fasta_cache[cache_key]
                self._save_cache(self.fasta_cache_file, self.fasta_cache)
        return False
    
    def _cache_fasta(self, gene_name, structure_id, file_path, structure_type='experimental'):
        """Add FASTA to cache"""
        cache_key = f"{gene_name}_{structure_id}_{structure_type}"
        self.fasta_cache[cache_key] = str(file_path)
        self._save_cache(self.fasta_cache_file, self.fasta_cache)
    
    def get_reference_sequence(self, gene_name):
        """Get reference sequence from preferred organisms"""
        gene_dir = self.sequences_dir / gene_name
        
        if not gene_dir.exists():
            logging.error(f"Gene directory not found: {gene_dir}")
            return None
        
        logging.debug(f"Searching for reference sequence for {gene_name}")
        
        # Try preferred organisms in order
        for organism in self.reference_organisms:
            fasta_file = gene_dir / f"{organism}.fasta"
            if fasta_file.exists():
                try:
                    sequences = list(SeqIO.parse(fasta_file, "fasta"))
                    if sequences:
                        # Use first sequence if multiple exist
                        ref_seq = sequences[0]
                        logging.info(f"✓ Using reference sequence from {organism}: {len(ref_seq.seq)} residues")
                        return ref_seq
                except Exception as e:
                    logging.warning(f"Error reading {fasta_file}: {e}")
                    continue
            else:
                logging.debug(f"No sequence file found for {organism}")
        
        logging.info(f"No preferred organism sequences found for {gene_name}, trying any available...")
        
        # If no preferred organism found, use the first available sequence
        fasta_files = list(gene_dir.glob("*.fasta"))
        if not fasta_files:
            logging.error(f"No FASTA files found in {gene_dir}")
            return None
            
        logging.info(f"Found {len(fasta_files)} alternative sequences, trying in order...")
        for fasta_file in fasta_files:
            if fasta_file.name.endswith("_metadata.json"):
                continue
            try:
                sequences = list(SeqIO.parse(fasta_file, "fasta"))
                if sequences:
                    ref_seq = sequences[0]
                    organism_name = fasta_file.stem.replace("_", " ")
                    logging.info(f"✓ Using reference sequence from {organism_name}: {len(ref_seq.seq)} residues")
                    return ref_seq
            except Exception as e:
                logging.warning(f"Error reading {fasta_file}: {e}")
                continue
        
        logging.error(f"No valid reference sequence found for {gene_name} - all files failed to parse")
        return None
    
    def search_similar_structures(self, gene_name, sequence, max_structures=20, identity_cutoff=0.3, evalue_cutoff=0.1, experimental_only=False):
        """Search for similar structures using PDB REST API"""
        search_type = "experimental only" if experimental_only else "experimental and computational"
        logging.info(f"Searching PDB for {search_type} structures similar to {gene_name} (identity >{identity_cutoff*100}%)")
        
        # PDB search API endpoint
        pdb_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        
        # Use the working PDB API query format (based on user-provided example)
        search_query = {
            "query": {
                "type": "group",
                "nodes": [
                    {
                        "type": "group",
                        "nodes": [
                            {
                                "type": "terminal",
                                "service": "full_text",
                                "parameters": {
                                    "value": gene_name
                                }
                            }
                        ],
                        "logical_operator": "and"
                    },
                    {
                        "type": "terminal",
                        "service": "sequence",
                        "parameters": {
                            "evalue_cutoff": evalue_cutoff,
                            "identity_cutoff": identity_cutoff,
                            "sequence_type": "protein",
                            "value": str(sequence.seq)
                        }
                    }
                ],
                "logical_operator": "and"
            },
            "return_type": "entry",
            "request_options": {
                "results_content_type": [
                    "experimental"
                ] if experimental_only else [
                    "computational",
                    "experimental"
                ],
                "paginate": {
                    "start": 0,
                    "rows": 25
                },
                "sort": [
                    {
                        "sort_by": "score",
                        "direction": "desc"
                    }
                ],
                "scoring_strategy": "combined"
            }
        }
        
        try:
            logging.info("Submitting PDB search...")
            response = requests.post(pdb_url, json=search_query, timeout=60)
            response.raise_for_status()
            
            # Check if response is empty or not JSON
            if not response.text.strip():
                logging.error("PDB search returned empty response")
                return []
            
            try:
                results = response.json()
            except json.JSONDecodeError as json_error:
                logging.error(f"PDB search returned invalid JSON. Response: {response.text[:200]}...")
                return []
                
            return self._parse_pdb_results(results, max_structures)
            
        except requests.exceptions.RequestException as req_error:
            logging.error(f"PDB search request failed: {req_error}")
            return []
        except Exception as e:
            logging.error(f"PDB search failed: {e}")
            return []
    
    def _is_likely_alphafold_structure(self, identifier, result):
        """Determine if a corrupted identifier likely represents an AlphaFold structure"""
        
        # Simple rule: if identifier starts with AF_, it's AlphaFold (corrupted or not)
        # If it doesn't start with AF_, it's experimental (format: PDB_ID_entity)
        return identifier.startswith('AF_')

    def _parse_pdb_results(self, pdb_data, max_structures):
        """Parse PDB search results and extract structure information"""
        similar_structures = []
        
        try:
            if 'result_set' not in pdb_data:
                logging.info("No PDB search results found")
                return []
            
            results = pdb_data['result_set']
            logging.info(f"Found {len(results)} PDB search results")
            
            experimental_structures = []
            alphafold_structures = []
            
            for result in results:
                identifier = result.get('identifier', '')
                score = result.get('score', 0)
                
                # Debug: Print raw identifiers to understand the corruption pattern
                logging.debug(f"Raw identifier: '{identifier}' -> is_AF_prefixed: {identifier.startswith('AF_')}")
                
                # Determine if this is experimental or AlphaFold structure
                # Due to PDB API corruption, all identifiers start with AF_ but we need to distinguish them
                is_alphafold = self._is_likely_alphafold_structure(identifier, result)
                
                if is_alphafold:
                    # AlphaFold structure - extract UniProt ID correctly
                    # Handle corrupted identifiers like AF_AFQ5L1S5F1 (double AF prefix)
                    # Expected format: AF_UNIPROTID_1 or AF_UNIPROTID_X
                    
                    # Fix corrupted AlphaFold identifier and extract UniProt ID
                    fixed_id = fix_corrupted_alphafold_identifier(identifier)
                    
                    # Extract UniProt ID from fixed identifier (AF-UNIPROTID-F1)
                    if fixed_id.startswith('AF-') and fixed_id.endswith('-F1'):
                        uniprot_id = fixed_id[3:-3]  # Remove 'AF-' prefix and '-F1' suffix
                    else:
                        # Fallback to original logic for non-standard formats
                        without_af = identifier[3:]  # Remove 'AF_'
                        parts = without_af.split('_') if '_' in without_af else [without_af]
                        uniprot_id = parts[0]
                    
                    # Validate that it looks like a proper UniProt ID (6+ characters, alphanumeric)
                    if len(uniprot_id) >= 6 and uniprot_id.replace('_', '').isalnum():
                        logging.debug(f"Parsed AlphaFold ID: {identifier} -> {uniprot_id}")
                        alphafold_structures.append({
                            'type': 'alphafold',
                            'alphafold_id': uniprot_id,
                            'identifier': identifier,
                            'score': score
                        })
                    else:
                        logging.debug(f"Invalid UniProt ID format, skipping: {identifier} -> {uniprot_id}")
                else:
                    # Experimental PDB structure
                    if '_' in identifier:
                        pdb_id = identifier.split('_')[0].upper()
                        entity_id = identifier.split('_')[1]
                    else:
                        pdb_id = identifier.upper()
                        entity_id = '1'
                    
                    # Extract sequence similarity information if available
                    services_results = result.get('services', [])
                    identity = None
                    evalue = None
                    for service in services_results:
                        if service.get('service_type') == 'sequence':
                            nodes = service.get('nodes', [])
                            for node in nodes:
                                if 'similarity_cutoff' in node:
                                    identity = node.get('similarity_cutoff', 0) * 100
                                if 'evalue_cutoff' in node:
                                    evalue = node.get('evalue_cutoff', 1.0)
                    
                    experimental_structures.append({
                        'type': 'experimental',
                        'pdb_id': pdb_id,
                        'entity_id': entity_id,
                        'identifier': identifier,
                        'score': score,
                        'identity': identity if identity else 'Unknown',
                        'evalue': evalue if evalue else 'Unknown'
                    })
                    
                    logging.debug(f"Found experimental structure: {pdb_id} (score: {score:.2f})")
            
            # Strict prioritization: experimental only, or AlphaFold only if no experimental
            logging.info(f"Found {len(experimental_structures)} experimental and {len(alphafold_structures)} AlphaFold structures")
            
            if experimental_structures:
                # Use only experimental structures if any are available
                similar_structures = experimental_structures[:max_structures]
                logging.info(f"Using {len(similar_structures)} experimental structures (ignoring {len(alphafold_structures)} AlphaFold)")
            else:
                # Use AlphaFold structures only if no experimental structures found
                similar_structures = alphafold_structures[:max_structures]
                logging.info(f"No experimental structures found, using {len(similar_structures)} AlphaFold structures")
            
            # Sort by score (descending)
            similar_structures.sort(key=lambda x: x['score'], reverse=True)
            
            logging.info(f"Processed {len(similar_structures)} similar structures")
            
        except Exception as e:
            logging.error(f"Error parsing PDB results: {e}")
        
        return similar_structures
    
    def download_pdb_structure(self, pdb_id, structures_dir):
        """Download PDB structure file to central structures directory"""
        structures_dir = Path(structures_dir)
        structures_dir.mkdir(parents=True, exist_ok=True)
        
        # Check cache first
        if self._is_structure_cached(pdb_id, 'experimental'):
            logging.debug(f"Structure {pdb_id} found in cache")
            return True
        
        # Try different formats
        formats = [
            ('pdb.gz', f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb.gz"),
            ('cif.gz', f"https://files.rcsb.org/download/{pdb_id.lower()}.cif.gz"),
            ('pdb', f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb")
        ]
        
        for ext, url in formats:
            output_file = structures_dir / f"{pdb_id}.{ext}"
            
            if output_file.exists() and output_file.stat().st_size > 0:
                logging.debug(f"Structure {pdb_id} already exists: {output_file}")
                return True
            
            try:
                logging.info(f"Downloading {pdb_id}.{ext}...")
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                
                if output_file.stat().st_size > 0:
                    logging.info(f"✓ Downloaded {pdb_id}.{ext}")
                    # Cache the downloaded structure
                    self._cache_structure(pdb_id, output_file, 'experimental')
                    return True
                else:
                    output_file.unlink()  # Remove empty file
                    
            except Exception as e:
                logging.debug(f"Failed to download {pdb_id}.{ext}: {e}")
                if output_file.exists():
                    output_file.unlink()
                continue
        
        logging.warning(f"✗ Failed to download structure: {pdb_id}")
        return False
    
    def download_alphafold_structure(self, uniprot_id, structures_dir):
        """Download AlphaFold structure if available to central structures directory"""
        structures_dir = Path(structures_dir)
        structures_dir.mkdir(parents=True, exist_ok=True)
        
        # Check cache first
        if self._is_structure_cached(f"AF-{uniprot_id}", 'computed'):
            logging.debug(f"AlphaFold structure AF-{uniprot_id} found in cache")
            return True
        
        # Try different formats: PDB first, then CIF
        formats = [
            ('pdb', f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"),
            ('cif', f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif")
        ]
        
        for format_ext, alphafold_url in formats:
            output_file = structures_dir / f"AF-{uniprot_id}.{format_ext}"
            
            # Check if file already exists
            if output_file.exists() and output_file.stat().st_size > 0:
                logging.debug(f"AlphaFold structure already exists: {output_file}")
                self._cache_structure(f"AF-{uniprot_id}", output_file, 'computed')
                return True
            
            try:
                logging.info(f"Downloading AlphaFold structure: AF-{uniprot_id} ({format_ext.upper()})")
                response = requests.get(alphafold_url, timeout=30)
                response.raise_for_status()
                
                with open(output_file, 'w') as f:
                    f.write(response.text)
                
                if output_file.stat().st_size > 0:
                    logging.info(f"✓ Downloaded AlphaFold structure: AF-{uniprot_id} ({format_ext.upper()})")
                    # Cache the downloaded structure
                    self._cache_structure(f"AF-{uniprot_id}", output_file, 'computed')
                    return True
                else:
                    logging.debug(f"Downloaded empty {format_ext.upper()} file for {uniprot_id}")
                    output_file.unlink()
                    
            except Exception as e:
                logging.debug(f"Failed to download AlphaFold {format_ext.upper()} {uniprot_id}: {e}")
                if output_file.exists():
                    output_file.unlink()
                continue
        
        logging.debug(f"No AlphaFold structure available for {uniprot_id} in any format")
        return False
    
    def download_pdb_fasta(self, pdb_id, entity_id, gene_name):
        """Download FASTA sequence from PDB for a specific structure"""
        gene_fasta_dir = self.protein_structures_base / gene_name
        gene_fasta_dir.mkdir(parents=True, exist_ok=True)
        
        # Check cache first
        if self._is_fasta_cached(gene_name, f"{pdb_id}_{entity_id}", 'experimental'):
            logging.debug(f"FASTA {pdb_id}_{entity_id} for {gene_name} found in cache")
            return True
        
        output_file = gene_fasta_dir / f"{pdb_id}_{entity_id}.fasta"
        
        if output_file.exists() and output_file.stat().st_size > 0:
            logging.debug(f"FASTA {pdb_id} already exists: {output_file}")
            return True
        
        # PDB GraphQL API for FASTA sequence with detailed metadata
        graphql_url = "https://data.rcsb.org/graphql"
        query = {
            "query": f'''
            {{
              entry(entry_id: "{pdb_id.upper()}") {{
                polymer_entities {{
                  entity_poly {{
                    pdbx_seq_one_letter_code_can
                  }}
                  rcsb_id
                  rcsb_polymer_entity_container_identifiers {{
                    entry_id
                    entity_id
                  }}
                  rcsb_entity_source_organism {{
                    ncbi_taxonomy_id
                    scientific_name
                  }}
                  rcsb_polymer_entity {{
                    pdbx_description
                  }}
                  polymer_entity_instances {{
                    rcsb_polymer_entity_instance_container_identifiers {{
                      auth_asym_id
                    }}
                  }}
                }}
                struct {{
                  title
                }}
              }}
            }}
            '''
        }
        
        try:
            logging.debug(f"Downloading FASTA for {pdb_id} entity {entity_id}...")
            response = requests.post(graphql_url, json=query, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            if 'data' in data and data['data']['entry']:
                polymer_entities = data['data']['entry']['polymer_entities']
                
                for entity in polymer_entities:
                    if entity['rcsb_polymer_entity_container_identifiers']['entity_id'] == entity_id:
                        sequence = entity['entity_poly']['pdbx_seq_one_letter_code_can']
                        if sequence:
                            # Extract additional metadata
                            chain_id = "A"  # Default
                            if (entity.get('polymer_entity_instances') and 
                                len(entity['polymer_entity_instances']) > 0 and
                                entity['polymer_entity_instances'][0].get('rcsb_polymer_entity_instance_container_identifiers') and
                                entity['polymer_entity_instances'][0]['rcsb_polymer_entity_instance_container_identifiers'].get('auth_asym_id')):
                                chain_id = entity['polymer_entity_instances'][0]['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id']
                            
                            # Extract organism information
                            organism_info = ""
                            if (entity.get('rcsb_entity_source_organism') and 
                                len(entity['rcsb_entity_source_organism']) > 0):
                                organism_data = entity['rcsb_entity_source_organism'][0]
                                scientific_name = organism_data.get('scientific_name', '')
                                taxonomy_id = organism_data.get('ncbi_taxonomy_id', '')
                                if scientific_name and taxonomy_id:
                                    organism_info = f"|{scientific_name} ({taxonomy_id})"
                                elif scientific_name:
                                    organism_info = f"|{scientific_name}"
                            
                            # Extract protein description
                            protein_description = ""
                            if (entity.get('rcsb_polymer_entity') and 
                                entity['rcsb_polymer_entity'].get('pdbx_description')):
                                protein_description = entity['rcsb_polymer_entity']['pdbx_description'][0] if isinstance(entity['rcsb_polymer_entity']['pdbx_description'], list) else entity['rcsb_polymer_entity']['pdbx_description']
                            
                            # Create detailed header: >PDB_ID_ENTITY|Chain X|Protein Description|Organism (TaxID)
                            if protein_description and organism_info:
                                header = f">{pdb_id.upper()}_{entity_id}|Chain {chain_id}|{protein_description}{organism_info}"
                            elif protein_description:
                                header = f">{pdb_id.upper()}_{entity_id}|Chain {chain_id}|{protein_description}"
                            else:
                                header = f">{pdb_id.upper()}_{entity_id}|Chain {chain_id}"
                            fasta_content = f"{header}\n{sequence}\n"
                            
                            with open(output_file, 'w') as f:
                                f.write(fasta_content)
                            
                            if output_file.stat().st_size > 0:
                                logging.debug(f"✓ Downloaded FASTA: {pdb_id}")
                                # Cache the downloaded FASTA
                                self._cache_fasta(gene_name, f"{pdb_id}_{entity_id}", output_file, 'experimental')
                                return True
                            else:
                                output_file.unlink()
                        break
            else:
                logging.debug(f"No sequence data found for {pdb_id} entity {entity_id}")
                
        except Exception as e:
            logging.debug(f"Failed to download FASTA {pdb_id}: {e}")
            if output_file.exists():
                output_file.unlink()
        
        return False
    
    def download_alphafold_fasta(self, uniprot_id, gene_name):
        """Download FASTA sequence from UniProt for AlphaFold structure"""
        gene_fasta_dir = self.protein_structures_base / gene_name
        gene_fasta_dir.mkdir(parents=True, exist_ok=True)
        
        # Check cache first
        if self._is_fasta_cached(gene_name, f"AF-{uniprot_id}", 'computed'):
            logging.debug(f"AlphaFold FASTA AF-{uniprot_id} for {gene_name} found in cache")
            return True
        
        output_file = gene_fasta_dir / f"AF-{uniprot_id}.fasta"
        
        if output_file.exists() and output_file.stat().st_size > 0:
            logging.debug(f"AlphaFold FASTA already exists: {output_file}")
            return True
        
        # UniProt FASTA URL
        fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        
        try:
            logging.debug(f"Downloading FASTA for AF-{uniprot_id}...")
            response = requests.get(fasta_url, timeout=30)
            response.raise_for_status()
            
            fasta_content = response.text.strip()
            if fasta_content and fasta_content.startswith('>'):
                # Enhance AlphaFold FASTA header with detailed metadata
                lines = fasta_content.split('\n')
                original_header = lines[0]
                
                # Extract detailed information from UniProt header
                # Format: >sp|UniProtID|NAME_ORGANISM Description OS=Organism OX=TaxID
                protein_description = ""
                organism_info = ""
                
                if ' OS=' in original_header:
                    # Extract protein description (everything before OS=)
                    desc_part = original_header.split(' OS=')[0]
                    if '|' in desc_part:
                        # Handle >sp|UniProtID|NAME format
                        parts = desc_part.split('|')
                        if len(parts) >= 3:
                            protein_description = ' '.join(parts[2:]).strip()
                    else:
                        # Direct description
                        protein_description = desc_part.replace('>', '').strip()
                    
                    # Extract organism and taxonomy ID
                    os_part = original_header.split(' OS=')[1]
                    if ' OX=' in os_part:
                        organism_name = os_part.split(' OX=')[0].strip()
                        taxonomy_id = os_part.split(' OX=')[1].split(' ')[0].strip()
                        organism_info = f"|{organism_name} ({taxonomy_id})"
                    else:
                        organism_name = os_part.split(' ')[0].strip()
                        organism_info = f"|{organism_name}"
                
                # Create detailed header: >AF-UniProtID|AlphaFold|Protein Description|Organism (TaxID)
                enhanced_header = f">AF-{uniprot_id}|AlphaFold|{protein_description}{organism_info}"
                
                # Reconstruct FASTA with enhanced header
                enhanced_fasta = enhanced_header + '\n' + '\n'.join(lines[1:])
                
                with open(output_file, 'w') as f:
                    f.write(enhanced_fasta)
                
                if output_file.stat().st_size > 0:
                    logging.debug(f"✓ Downloaded AlphaFold FASTA: AF-{uniprot_id}")
                    # Cache the downloaded FASTA
                    self._cache_fasta(gene_name, f"AF-{uniprot_id}", output_file, 'computed')
                    return True
                else:
                    output_file.unlink()
            else:
                logging.debug(f"Invalid FASTA content for AF-{uniprot_id}")
                
        except Exception as e:
            logging.debug(f"Failed to download AlphaFold FASTA {uniprot_id}: {e}")
            if output_file.exists():
                output_file.unlink()
        
        return False
    
    def search_alphafold_structures(self, gene_name, reference_seq):
        """Search for AlphaFold structures using UniProt"""
        logging.info(f"Searching AlphaFold structures for {gene_name}")
        
        uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'gene:{gene_name} AND reviewed:true AND taxonomy_id:2 AND structure_3d:true',
            'format': 'json',
            'limit': 10
        }
        
        try:
            response = requests.get(uniprot_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            alphafold_structures = []
            for entry in data.get('results', []):
                uniprot_id = entry.get('primaryAccession', '')
                organism = entry.get('organism', {}).get('scientificName', 'Unknown')
                
                if uniprot_id:
                    alphafold_structures.append({
                        'uniprot_id': uniprot_id,
                        'organism': organism,
                        'type': 'computed'
                    })
            
            logging.info(f"Found {len(alphafold_structures)} potential AlphaFold structures")
            return alphafold_structures
            
        except Exception as e:
            logging.error(f"Error searching AlphaFold structures: {e}")
            return []
    
    def get_reference_sequence_by_organism(self, gene_name, organism):
        """Get reference sequence from a specific organism"""
        gene_dir = self.sequences_dir / gene_name
        fasta_file = gene_dir / f"{organism}.fasta"
        
        if not fasta_file.exists():
            return None
            
        try:
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            if sequences:
                return sequences[0]
        except Exception as e:
            logging.debug(f"Error reading {fasta_file}: {e}")
        return None
    
    def download_structures_for_gene_with_fallback(self, gene_name, max_structures=20, include_alphafold=True, 
                                                 identity_cutoff=0.3, evalue_cutoff=0.1, min_experimental=5):
        """Download structures for a gene using two-phase approach: experimental first, then computational"""
        logging.info(f"Starting structure download for gene: {gene_name}")
        
        gene_dir = self.sequences_dir / gene_name
        if not gene_dir.exists():
            return {"success": False, "error": f"Gene directory not found: {gene_dir}"}
        
        # Create output directory
        gene_output_dir = self.output_dir / gene_name
        gene_output_dir.mkdir(parents=True, exist_ok=True)
        
        # PHASE 1: Search for experimental structures only
        logging.info("=== PHASE 1: Searching for experimental structures ===")
        experimental_result = self._search_structures_phase(
            gene_name, max_structures, identity_cutoff, evalue_cutoff, 
            min_experimental, experimental_only=True
        )
        
        if experimental_result and experimental_result.get("success"):
            logging.info(f"✓ Found {experimental_result.get('experimental_count', 0)} experimental structures")
            return experimental_result
        
        # PHASE 2: No experimental structures found, search for computational structures
        logging.info("=== PHASE 2: No experimental structures found, searching for computational structures ===")
        computational_result = self._search_structures_phase(
            gene_name, max_structures, identity_cutoff, evalue_cutoff, 
            min_experimental, experimental_only=False
        )
        
        if computational_result and computational_result.get("success"):
            logging.info(f"✓ Found {computational_result.get('total_count', 0)} computational structures")
            return computational_result
        
        return {"success": False, "error": "No structures found in either phase"}
    
    def _search_structures_phase(self, gene_name, max_structures, identity_cutoff, evalue_cutoff, 
                               min_experimental, experimental_only=False):
        """Search for structures in one phase (experimental or computational)"""
        best_result = None
        best_structure_count = 0
        
        for organism in self.reference_organisms:
            reference_seq = self.get_reference_sequence_by_organism(gene_name, organism)
            if not reference_seq:
                continue
                
            logging.info(f"Trying {organism} as reference ({len(reference_seq.seq)} residues)")
            
            # Search for structures using PDB API
            similar_structures = self.search_similar_structures(
                gene_name, reference_seq, max_structures, identity_cutoff, evalue_cutoff, experimental_only
            )
            
            structure_count = len(similar_structures)
            
            if structure_count >= min_experimental:
                logging.info(f"✓ Found {structure_count} structures with {organism}, proceeding with download")
                return self.download_structures_for_gene_with_sequence(
                    gene_name, reference_seq, similar_structures, max_structures, 
                    not experimental_only, min_experimental, identity_cutoff, evalue_cutoff  # include_alphafold = not experimental_only
                )
            elif structure_count > best_structure_count:
                best_result = (organism, reference_seq, similar_structures)
                best_structure_count = structure_count
                logging.info(f"Found {structure_count} structures with {organism} (saving as backup)")
            else:
                logging.info(f"Found only {structure_count} structures with {organism}")
        
        # If we didn't find enough structures with any preferred organism, use the best one
        if best_result:
            organism, reference_seq, similar_structures = best_result
            logging.info(f"Using best available option: {organism} with {len(similar_structures)} structures")
            return self.download_structures_for_gene_with_sequence(
                gene_name, reference_seq, similar_structures, max_structures, 
                not experimental_only, min_experimental, identity_cutoff, evalue_cutoff
            )
        
        # Last resort: try any available sequence
        reference_seq = self.get_reference_sequence(gene_name)
        if not reference_seq:
            return None
            
        similar_structures = self.search_similar_structures(
            gene_name, reference_seq, max_structures, identity_cutoff, evalue_cutoff, experimental_only
        )
        
        if similar_structures:
            return self.download_structures_for_gene_with_sequence(
                gene_name, reference_seq, similar_structures, max_structures, 
                not experimental_only, min_experimental, identity_cutoff, evalue_cutoff
            )
        
        return None
    
    def download_structures_for_gene(self, gene_name, max_structures=20, include_alphafold=True, 
                                   identity_cutoff=0.3, evalue_cutoff=0.1, min_experimental=5):
        """Main function to download structures for a gene (uses fallback method)"""
        return self.download_structures_for_gene_with_fallback(
            gene_name, max_structures, include_alphafold, 
            identity_cutoff, evalue_cutoff, min_experimental
        )
    
    def download_structures_for_gene_with_sequence(self, gene_name, reference_seq, similar_structures, 
                                                 max_structures=20, include_alphafold=True, min_experimental=5,
                                                 identity_cutoff=0.3, evalue_cutoff=0.1):
        """Download structures for a gene using a specific reference sequence"""
        experimental_count = len(similar_structures)
        downloaded_experimental = 0
        downloaded_computed = 0
        
        # Create gene output directory
        gene_output_dir = self.output_dir / gene_name
        gene_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Download structures (both experimental and AlphaFold from PDB search)
        if similar_structures:
            logging.info(f"Downloading {len(similar_structures)} structures from PDB search...")
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
                future_to_struct = {}
                
                for struct in similar_structures:
                    if struct['type'] == 'experimental':
                        # Download structure to central directory
                        download_func = partial(self.download_pdb_structure, structures_dir=self.structures_dir)
                        future = executor.submit(download_func, struct['pdb_id'])
                        future_to_struct[future] = struct
                    elif struct['type'] == 'alphafold':
                        # Download structure to central directory
                        download_func = partial(self.download_alphafold_structure, structures_dir=self.structures_dir)
                        future = executor.submit(download_func, struct['alphafold_id'])
                        future_to_struct[future] = struct
                
                for future in concurrent.futures.as_completed(future_to_struct):
                    struct = future_to_struct[future]
                    try:
                        success = future.result()
                        if success:
                            if struct['type'] == 'experimental':
                                downloaded_experimental += 1
                                logging.info(f"✓ Downloaded experimental: {struct['pdb_id']}")
                                # Download corresponding FASTA
                                self.download_pdb_fasta(struct['pdb_id'], struct['entity_id'], gene_name)
                            else:
                                downloaded_computed += 1
                                logging.info(f"✓ Downloaded AlphaFold: {struct['alphafold_id']}")
                                # Download corresponding FASTA
                                self.download_alphafold_fasta(struct['alphafold_id'], gene_name)
                    except Exception as e:
                        struct_id = struct.get('pdb_id') or struct.get('alphafold_id')
                        logging.error(f"Error downloading {struct_id}: {e}")
        
        # Download AlphaFold structures if experimental structures are insufficient
        alphafold_structures = []
        if include_alphafold and downloaded_experimental < min_experimental:
            logging.info(f"Only {downloaded_experimental} experimental structures found, searching AlphaFold models...")
            alphafold_structures = self.search_alphafold_structures(gene_name, reference_seq)
            
            if alphafold_structures:
                with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
                    download_func = partial(self.download_alphafold_structure, structures_dir=self.structures_dir)
                    future_to_alphafold = {
                        executor.submit(download_func, struct['uniprot_id']): struct 
                        for struct in alphafold_structures
                    }
                    
                    for future in concurrent.futures.as_completed(future_to_alphafold):
                        struct = future_to_alphafold[future]
                        try:
                            success = future.result()
                            if success:
                                downloaded_computed += 1
                                # Download corresponding FASTA
                                self.download_alphafold_fasta(struct['uniprot_id'], gene_name)
                        except Exception as e:
                            logging.error(f"Error downloading AlphaFold {struct['uniprot_id']}: {e}")
        
        # Save metadata
        metadata = {
            'gene_name': gene_name,
            'reference_organism': reference_seq.description if reference_seq else 'Unknown',
            'reference_sequence_length': len(reference_seq.seq) if reference_seq else 0,
            'identity_cutoff': identity_cutoff,
            'evalue_cutoff': evalue_cutoff,
            'experimental_structures_found': experimental_count,
            'experimental_structures_downloaded': downloaded_experimental,
            'computed_structures_downloaded': downloaded_computed,
            'total_structures_downloaded': downloaded_experimental + downloaded_computed,
            'similar_structures': similar_structures,
            'alphafold_structures': alphafold_structures
        }
        
        metadata_file = gene_output_dir / "download_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logging.info(f"=== Download Summary for {gene_name} ===")
        logging.info(f"Experimental structures: {downloaded_experimental}/{experimental_count}")
        logging.info(f"Computed structures: {downloaded_computed}")
        logging.info(f"Total downloaded: {downloaded_experimental + downloaded_computed}")
        
        return {
            "success": True,
            "gene_name": gene_name,
            "experimental_downloaded": downloaded_experimental,
            "computed_downloaded": downloaded_computed,
            "total_downloaded": downloaded_experimental + downloaded_computed,
            "metadata": metadata
        }

def download_all_genes_from_summary(summary_file, output_base_dir, structures_dir, 
                                  max_structures=20, identity_cutoff=0.3, evalue_cutoff=0.1, 
                                  min_experimental=5, include_alphafold=True):
    """Download structures for all genes in summary file (Snakemake integration)"""
    import pandas as pd
    
    # Read summary file to get list of genes
    try:
        if summary_file.suffix == '.tsv':
            df = pd.read_csv(summary_file, sep='\t')
        else:
            df = pd.read_csv(summary_file)
        
        # Get unique gene names from 'gene_name' column (this is the correct column)
        if 'gene_name' in df.columns:
            genes = df['gene_name'].unique().tolist()
        elif 'Gene' in df.columns:
            genes = df['Gene'].unique().tolist()
        elif 'gene' in df.columns:
            genes = df['gene'].unique().tolist()
        else:
            # Fallback: use second column (gene_name is typically the second column)
            genes = df.iloc[:, 1].unique().tolist()
        
        logging.info(f"Found {len(genes)} unique genes in summary file")
        
    except Exception as e:
        logging.error(f"Error reading summary file {summary_file}: {e}")
        return {"success": False, "error": str(e)}
    
    # Create downloader
    downloader = SequenceSimilarityStructureDownloader(
        output_dir=output_base_dir,
        structures_dir=structures_dir
    )
    
    # Download structures for each gene
    all_results = {}
    successful_downloads = 0
    total_experimental = 0
    total_computed = 0
    mapping_entries = []
    
    for i, gene_name in enumerate(genes, 1):
        logging.info(f"Processing gene {i}/{len(genes)}: {gene_name}")
        
        try:
            result = downloader.download_structures_for_gene(
                gene_name=gene_name,
                max_structures=max_structures,
                include_alphafold=include_alphafold,
                identity_cutoff=identity_cutoff,
                evalue_cutoff=evalue_cutoff,
                min_experimental=min_experimental
            )
            
            all_results[gene_name] = result
            
            if result["success"] and result["total_downloaded"] > 0:
                successful_downloads += 1
                total_experimental += result["experimental_downloaded"]
                total_computed += result["computed_downloaded"]
                
                # Create mapping entries for FASTA-structure mapping file
                gene_output_dir = output_base_dir / gene_name
                
                # Check for downloaded FASTA files to create mapping
                gene_fasta_dir = output_base_dir.parent / "protein_structures" / gene_name
                fasta_files = list(gene_fasta_dir.glob("*.fasta")) if gene_fasta_dir.exists() else []
                
                # Create mapping for structures that have corresponding FASTA files
                for fasta_file in fasta_files:
                    fasta_id = fasta_file.stem  # e.g., "8BNZ_1" or "AF-Q8EM16"
                    
                    # Parse FASTA header to extract metadata
                    try:
                        with open(fasta_file, 'r') as f:
                            header_line = f.readline().strip()
                            if not header_line.startswith('>'):
                                continue
                            
                            # Remove '>' and split by '|'
                            header_parts = header_line[1:].split('|')
                            
                            # Initialize metadata
                            structure_id = fasta_id
                            protein_name = ""
                            chain = ""
                            species = ""
                            
                            if len(header_parts) >= 2:
                                if fasta_id.startswith('AF-'):
                                    # AlphaFold format: >AF-UniProtID|AlphaFold|Protein Description|Organism (TaxID)
                                    structure_type = 'computed'
                                    structure_id = header_parts[0]  # AF-UniProtID
                                    if len(header_parts) >= 3:
                                        protein_name = header_parts[2] if header_parts[2] else ""
                                    if len(header_parts) >= 4:
                                        species = header_parts[3] if header_parts[3] else ""
                                    chain = "A"  # AlphaFold structures are single chain
                                else:
                                    # Experimental format: >PDB_ID_ENTITY|Chain X|Protein Description|Organism (TaxID)
                                    structure_type = 'experimental'
                                    structure_id = header_parts[0]  # PDB_ID_ENTITY
                                    if len(header_parts) >= 2 and header_parts[1].startswith('Chain '):
                                        chain = header_parts[1].replace('Chain ', '')
                                    if len(header_parts) >= 3:
                                        protein_name = header_parts[2] if header_parts[2] else ""
                                    if len(header_parts) >= 4:
                                        species = header_parts[3] if header_parts[3] else ""
                            
                    except Exception as e:
                        logging.warning(f"Error parsing FASTA header from {fasta_file}: {e}")
                        # Fallback to basic parsing
                        structure_type = 'computed' if fasta_id.startswith('AF-') else 'experimental'
                        structure_id = fasta_id
                        protein_name = ""
                        chain = ""
                        species = ""
                    
                    # Find corresponding structure file in central structures directory
                    struct_file = None
                    
                    if structure_type == 'computed':
                        # For AlphaFold, extract UniProt ID from FASTA filename
                        if fasta_id.startswith('AF-'):
                            uniprot_id = fasta_id  # AF-UniProtID
                        else:
                            uniprot_id = structure_id
                            
                        # Look for AlphaFold structure
                        potential_files = [
                            output_base_dir.parent / "protein_structures" / "pdb_files" / f"{uniprot_id}.pdb",
                            output_base_dir.parent / "protein_structures" / "pdb_files" / f"{uniprot_id}.pdb.gz"
                        ]
                    else:
                        # For experimental, extract PDB ID from FASTA filename (remove entity suffix)
                        if '_' in fasta_id:
                            pdb_id = fasta_id.split('_')[0]  # Extract PDB ID from "8BNZ_1"
                        else:
                            pdb_id = fasta_id
                            
                        # Look for experimental structure
                        potential_files = [
                            output_base_dir.parent / "protein_structures" / "pdb_files" / f"{pdb_id}.pdb.gz",
                            output_base_dir.parent / "protein_structures" / "pdb_files" / f"{pdb_id}.cif.gz",
                            output_base_dir.parent / "protein_structures" / "pdb_files" / f"{pdb_id}.pdb"
                        ]
                    
                    # Find the first existing structure file
                    for potential_file in potential_files:
                        if potential_file.exists():
                            struct_file = potential_file
                            break
                    
                    if struct_file:
                        mapping_entries.append({
                            'structure_type': structure_type,
                            'structure_id': structure_id,
                            'gene_name': gene_name,
                            'protein_name': protein_name,
                            'chain': chain,
                            'species': species,
                            'structure_path': str(struct_file),
                            'fasta_path': str(fasta_file)
                        })
            
        except Exception as e:
            logging.error(f"Error processing gene {gene_name}: {e}")
            all_results[gene_name] = {"success": False, "error": str(e)}
    
    # Create mapping file for Snakemake compatibility
    if mapping_entries:
        mapping_df = pd.DataFrame(mapping_entries)
        mapping_file = output_base_dir.parent / "fasta_structure_mapping.tsv"
        mapping_df.to_csv(mapping_file, sep='\t', index=False)
        logging.info(f"Created mapping file with {len(mapping_entries)} entries: {mapping_file}")
    
    # Summary
    logging.info(f"=== Final Summary ===")
    logging.info(f"Total genes processed: {len(genes)}")
    logging.info(f"Genes with structures: {successful_downloads}")
    logging.info(f"Total experimental structures: {total_experimental}")
    logging.info(f"Total computed structures: {total_computed}")
    logging.info(f"Total structures downloaded: {total_experimental + total_computed}")
    
    return {
        "success": True,
        "total_genes": len(genes),
        "successful_downloads": successful_downloads,
        "total_experimental": total_experimental,
        "total_computed": total_computed,
        "total_structures": total_experimental + total_computed,
        "results": all_results
    }

def main():
    """Command line interface"""
    # Check if running under Snakemake
    try:
        # Snakemake integration
        summary_file = Path(snakemake.input.protein_list)
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        max_structures = snakemake.params.max_structures
        include_alphafold = snakemake.params.include_computed_models
        
        # Get additional parameters with defaults
        identity_cutoff = getattr(snakemake.params, 'identity_cutoff', 0.3)
        evalue_cutoff = getattr(snakemake.params, 'evalue_cutoff', 0.1)
        min_experimental = getattr(snakemake.params, 'min_experimental', 5)
        
        # Create output directories
        sentinel_file = Path(snakemake.output.sentinel)
        mapping_file = Path(snakemake.output.mapping_file)
        
        # Use the protein_structures base directory directly
        protein_structures_base = Path(__file__).parent.parent.parent / "data" / "protein_structures"
        output_base_dir = protein_structures_base
        structures_dir = protein_structures_base / "pdb_files"
        
        logging.info(f"Running Snakemake integration for {analysis}_{paramset}")
        logging.info(f"Summary file: {summary_file}")
        logging.info(f"Output directory: {output_base_dir}")
        
        # Download all structures
        result = download_all_genes_from_summary(
            summary_file=summary_file,
            output_base_dir=output_base_dir,
            structures_dir=structures_dir,
            max_structures=max_structures,
            identity_cutoff=identity_cutoff,
            evalue_cutoff=evalue_cutoff,
            min_experimental=min_experimental,
            include_alphafold=include_alphafold
        )
        
        if result["success"]:
            # Create sentinel file
            sentinel_file.touch()
            
            # Copy mapping file to expected location if it exists
            temp_mapping = output_base_dir.parent / "fasta_structure_mapping.tsv"
            if temp_mapping.exists():
                import shutil
                shutil.move(str(temp_mapping), str(mapping_file))
            else:
                # Create empty mapping file with new column format
                with open(mapping_file, 'w') as f:
                    f.write("structure_type\tstructure_id\tgene_name\tprotein_name\tchain\tspecies\tstructure_path\tfasta_path\n")
            
            logging.info(f"Snakemake integration completed successfully")
            return 0
        else:
            logging.error(f"Snakemake integration failed: {result.get('error', 'Unknown error')}")
            return 1
            
    except NameError:
        # Command line interface
        parser = argparse.ArgumentParser(description="Download protein structures based on sequence similarity")
        parser.add_argument("gene_name", help="Gene name to search for (e.g., bamA)")
        parser.add_argument("--max-structures", type=int, default=20, 
                           help="Maximum number of structures to search for (default: 20)")
        parser.add_argument("--identity-cutoff", type=float, default=0.3,
                           help="Minimum sequence identity cutoff 0.0-1.0 (default: 0.3)")
        parser.add_argument("--evalue-cutoff", type=float, default=0.1,
                           help="Maximum E-value cutoff (default: 0.1)")
        parser.add_argument("--min-experimental", type=int, default=5,
                           help="Minimum experimental structures before including AlphaFold (default: 5)")
        parser.add_argument("--no-alphafold", action="store_true",
                           help="Don't include AlphaFold structures")
        parser.add_argument("--output-dir", type=str,
                           help="Output directory for structures")
        
        args = parser.parse_args()
        
        # Create downloader
        downloader = SequenceSimilarityStructureDownloader(output_dir=args.output_dir)
        
        # Download structures
        result = downloader.download_structures_for_gene(
            gene_name=args.gene_name,
            max_structures=args.max_structures,
            include_alphafold=not args.no_alphafold,
            identity_cutoff=args.identity_cutoff,
            evalue_cutoff=args.evalue_cutoff,
            min_experimental=args.min_experimental
        )
        
        if result["success"]:
            print(f"Successfully downloaded {result['total_downloaded']} structures for {args.gene_name}")
            print(f"  Experimental: {result['experimental_downloaded']}")
            print(f"  Computed: {result['computed_downloaded']}")
        else:
            print(f"Failed to download structures: {result.get('error', 'Unknown error')}")
            return 1
        
        return 0

if __name__ == "__main__":
    exit(main())