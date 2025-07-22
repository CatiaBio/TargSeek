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

# Cache configuration - use relative paths since Snakemake runs from project root
CACHE_DIR = Path("cache/protein_structures")
CACHE_FILE = CACHE_DIR / "protein_3d_structure_cache.json"

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
    """Load gene aliases from JSON file"""
    try:
        with open(aliases_file, 'r') as f:
            aliases = json.load(f)
        logging.info(f"Loaded aliases for {len(aliases)} genes")
        return aliases
    except Exception as e:
        logging.warning(f"Could not load gene aliases: {e}")
        return {}

def search_pdb_via_uniprot_bacterial(gene_name, max_structures=3, gene_aliases=None, cache=None):
    """
    Search UniProt for bacterial proteins with PDB structures, with alias fallback
    
    Args:
        gene_name (str): Primary gene name to search for
        max_structures (int): Maximum structures to return
        gene_aliases (list): Alternative gene names to try if primary fails
        cache (Protein3DStructureCache): Cache instance to use
    """
    # Check cache first
    if cache:
        cached_result = cache.get_cached_search(gene_name)
        if cached_result:
            logging.info(f"Found cached 3D structure search for {gene_name}: {len(cached_result['pdb_ids'])} structures")
            return cached_result['pdb_ids'][:max_structures]
    
    # Try primary gene name first
    pdb_ids = _search_uniprot_bacterial_single_gene(gene_name, max_structures)
    
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
                alias_pdb_ids = _search_uniprot_bacterial_single_gene(alias, max_structures)
                if alias_pdb_ids:
                    logging.info(f"Found {len(alias_pdb_ids)} bacterial PDB structures via alias '{alias}'")
                    if cache:
                        cache.cache_search_result(gene_name, alias_pdb_ids)
                    return alias_pdb_ids
    
    logging.info(f"No bacterial PDB structures found for {gene_name} (tried {len(gene_aliases or [])} aliases)")
    if cache:
        cache.cache_search_result(gene_name, [])  # Cache empty result
    return []

def _search_uniprot_bacterial_single_gene(gene_name, max_structures):
    """Helper function to search UniProt for bacterial proteins with structures for a single gene name"""
    try:
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'gene:{gene_name} AND structure_3d:true',
            'format': 'json',
            'size': max_structures * 3,  # Get more results to find bacterial ones (increased from 2x to 3x)
            'fields': 'accession,gene_names,organism_name,xref_pdb'
        }
        
        response = requests.get(uniprot_url, params=params, timeout=30)
        if response.status_code != 200:
            logging.debug(f"UniProt bacterial search failed for {gene_name}: status {response.status_code}")
            return []
        
        data = response.json()
        bacterial_pdb_ids = []
        other_pdb_ids = []
        
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
            for xref in entry.get('uniProtKBCrossReferences', []):
                if xref.get('database') == 'PDB':
                    pdb_id = xref.get('id')
                    if pdb_id:
                        if is_bacterial:
                            if pdb_id not in bacterial_pdb_ids:
                                bacterial_pdb_ids.append(pdb_id)
                                logging.debug(f"Found bacterial protein from {organism}")
                        else:
                            if pdb_id not in other_pdb_ids:
                                other_pdb_ids.append(pdb_id)
                                logging.debug(f"Found non-bacterial protein from {organism}")
        
        # Prefer bacterial structures, but use others if no bacterial ones found
        if bacterial_pdb_ids:
            logging.info(f"Found {len(bacterial_pdb_ids)} bacterial structures for {gene_name}")
            return bacterial_pdb_ids[:max_structures]
        elif other_pdb_ids:
            logging.info(f"No bacterial structures found for {gene_name}, using {len(other_pdb_ids)} non-bacterial structures")
            return other_pdb_ids[:max_structures]
        else:
            return []
        
    except Exception as e:
        logging.warning(f"UniProt bacterial search failed for {gene_name}: {e}")
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

def download_single_pdb_structure(pdb_id, output_dir, cache=None):
    """
    Download a single PDB structure - helper function for concurrent downloads
    """
    return download_pdb_structure(pdb_id, output_dir, cache)

def download_pdb_structure(pdb_id, output_dir, cache=None):
    """
    Download PDB structure file in legacy format (gz) with caching
    
    Args:
        pdb_id (str): PDB identifier
        output_dir (Path): Output directory to save the structure
        cache (Protein3DStructureCache): Cache instance to use
    
    Returns:
        bool: True if download succeeded, False otherwise
    """
    
    # Check if file already exists (quick check first)
    pdb_file = output_dir / f"{pdb_id}.pdb.gz"
    if pdb_file.exists() and pdb_file.stat().st_size > 0:
        # File exists and is not empty, no need to re-download
        if cache:
            cache.cache_pdb_result(pdb_id, True, pdb_file)
        return True
    
    # Check cache for download status only if file doesn't exist
    if cache:
        cached_pdb = cache.get_cached_pdb(pdb_id)
        if cached_pdb and cached_pdb.get('success'):
            cached_path = Path(cached_pdb.get('file_path', ''))
            if cached_path.exists() and cached_path.stat().st_size > 0:
                return True
    
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb.gz"
    
    try:
        response = requests.get(pdb_url, timeout=30)
        response.raise_for_status()
        
        # Save PDB file
        with open(pdb_file, 'wb') as f:
            f.write(response.content)
        
        logging.info(f"Downloaded PDB structure: {pdb_file}")
        if cache:
            cache.cache_pdb_result(pdb_id, True, pdb_file)
        return True
        
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to download PDB structure for {pdb_id}: {e}")
        if cache:
            cache.cache_pdb_result(pdb_id, False)
        return False

def download_3d_structures_for_gene(gene_name, structures_base_dir, organisms=None, max_structures=3, gene_aliases=None, cache=None):
    """
    Download 3D structures and sequences for a specific bacterial gene (top 3 best scoring only)
    
    Args:
        gene_name (str): Gene name to search for
        structures_base_dir (Path): Base directory for 3D structures
        organisms (list): Optional list of organism names to filter by (ignored for bacterial search)
        max_structures (int): Maximum number of structures to download (default: 3)
        gene_aliases (list): Alternative gene names to try if primary search fails
        cache (Protein3DStructureCache): Cache instance to use
    
    Returns:
        dict: Summary of downloads (sequences, structures, found)
    """
    
    # Create gene-specific directory in 3d_structures
    gene_dir = structures_base_dir / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if gene already has structures (file-based completion checking)
    existing_pdb_files = list(gene_dir.glob("*.pdb.gz"))
    existing_fasta_files = list(gene_dir.glob("*.fasta"))
    has_no_structures_marker = (gene_dir / "no_structures_found.txt").exists()
    has_completion_marker = (gene_dir / ".processing_complete.txt").exists()
    
    # If directory has any structure files (PDB or FASTA), consider it complete
    # Some genes might have less than 3 structures available
    if existing_pdb_files or existing_fasta_files:
        logging.info(f"Gene {gene_name} already has structures: {len(existing_pdb_files)} PDB files, {len(existing_fasta_files)} FASTA files - SKIPPING")
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
    
    # Check if marked as completed but somehow files were deleted
    if has_completion_marker and not existing_pdb_files and not existing_fasta_files:
        logging.warning(f"Gene {gene_name} marked as complete but no files found - will retry download")
        # Remove the completion marker to allow retry
        (gene_dir / ".processing_complete.txt").unlink(missing_ok=True)
    
    logging.info(f"Processing bacterial 3D structures for gene: {gene_name}")
    
    # Search for bacterial PDB structures via UniProt
    pdb_ids = search_pdb_via_uniprot_bacterial(gene_name, max_structures, gene_aliases, cache)
    
    if not pdb_ids:
        logging.info(f"No bacterial 3D structures found for {gene_name}")
        # Create a marker file indicating no structures found (completion detection)
        no_structures_file = gene_dir / "no_structures_found.txt"
        with open(no_structures_file, 'w') as f:
            f.write(f"No bacterial 3D structures found for gene: {gene_name}\n")
            f.write(f"Search performed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            # Include alias information for completeness
            if gene_aliases:
                f.write(f"Aliases tried: {', '.join(gene_aliases)}\n")
            else:
                f.write("No aliases available\n")
            f.write(f"Status: COMPLETED (no structures available)\n")
        
        return {"sequences": 0, "structures": 0, "found": False, "skipped": False, "completed": True, "pdb_ids": []}
    
    # Limit to top 3 structures
    pdb_ids = pdb_ids[:max_structures]
    logging.info(f"Processing top {len(pdb_ids)} bacterial 3D structures for {gene_name}: {pdb_ids}")
    
    total_sequences = 0
    total_structures = 0
    downloaded_pdb_ids = []
    
    # Process FASTA downloads first (sequential for logging clarity)
    pdb_with_sequences = []
    for pdb_id in pdb_ids:
        sequences = get_pdb_fasta(pdb_id)
        
        if sequences:
            pdb_with_sequences.append((pdb_id, sequences))
            for i, seq_record in enumerate(sequences):
                # Create a unique filename for each sequence
                if len(sequences) > 1:
                    filename = f"{pdb_id}_chain_{i+1}.fasta"
                else:
                    filename = f"{pdb_id}.fasta"
                
                # Update sequence ID to include PDB info
                seq_record.id = f"{pdb_id}_{seq_record.id}"
                seq_record.description = f"PDB:{pdb_id} {seq_record.description}"
                
                # Save sequence in gene directory
                output_file = gene_dir / filename
                with open(output_file, 'w') as f:
                    SeqIO.write(seq_record, f, "fasta")
                
                total_sequences += 1
                logging.info(f"Saved bacterial 3D structure sequence: {output_file}")
        else:
            logging.warning(f"No sequences found for PDB {pdb_id}")
    
    # Download PDB structures concurrently for better performance
    if pdb_with_sequences:
        pdb_ids_to_download = [pdb_id for pdb_id, _ in pdb_with_sequences]
        logging.info(f"Downloading {len(pdb_ids_to_download)} PDB structures concurrently...")
        
        # Use ThreadPoolExecutor for concurrent downloads (limited to 3 workers to be respectful)
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            # Create partial function with fixed arguments
            download_func = partial(download_pdb_structure, output_dir=gene_dir, cache=cache)
            
            # Submit all download tasks
            future_to_pdb = {executor.submit(download_func, pdb_id): pdb_id for pdb_id in pdb_ids_to_download}
            
            # Collect results
            for future in concurrent.futures.as_completed(future_to_pdb):
                pdb_id = future_to_pdb[future]
                try:
                    success = future.result()
                    if success:
                        total_structures += 1
                        downloaded_pdb_ids.append(pdb_id)
                        logging.info(f"✓ Downloaded PDB structure: {pdb_id}")
                    else:
                        logging.warning(f"✗ Failed to download PDB structure: {pdb_id}")
                except Exception as e:
                    logging.error(f"Error downloading PDB {pdb_id}: {e}")
    
    # Mark gene as completed by creating a completion marker if we successfully processed it
    if total_structures > 0:
        completion_file = gene_dir / ".processing_complete.txt"
        with open(completion_file, 'w') as f:
            f.write(f"Gene: {gene_name}\n")
            f.write(f"Completed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Downloaded structures: {total_structures}\n")
            f.write(f"Downloaded sequences: {total_sequences}\n")
            f.write(f"PDB IDs: {', '.join(downloaded_pdb_ids)}\n")
            f.write(f"Status: COMPLETED (structures downloaded)\n")
    
    return {
        "sequences": total_sequences, 
        "structures": total_structures, 
        "found": total_sequences > 0,
        "skipped": False,
        "completed": True,
        "pdb_ids": downloaded_pdb_ids
    }

def process_proteins_list(proteins_file, structures_base_dir, gene_aliases=None, cache=None):
    """
    Process proteins_to_study file and download bacterial 3D structures for each gene
    
    Args:
        proteins_file (Path): Path to proteins_to_study TSV file
        structures_base_dir (Path): Base directory for 3D structures
        gene_aliases (dict): Dictionary mapping gene names to lists of aliases
        cache (Protein3DStructureCache): Cache instance to use
    
    Returns:
        dict: Summary of downloaded structures per gene
    """
    
    # Load proteins list
    if not proteins_file.exists():
        logging.error(f"Proteins file not found: {proteins_file}")
        return {}
    
    proteins_df = pd.read_csv(proteins_file, sep='\t')
    total_genes = len(proteins_df)
    logging.info(f"Processing {total_genes} genes from {proteins_file}")
    
    summary = {}
    genes_skipped = 0
    genes_processed = 0
    
    for gene_idx, (_, row) in enumerate(proteins_df.iterrows(), 1):
        gene_name = row['protein']
        
        # Progress logging every 10 genes
        if gene_idx % 10 == 0 or gene_idx == 1:
            logging.info(f"Progress: {gene_idx}/{total_genes} genes ({gene_idx/total_genes*100:.1f}%) - Processed: {genes_processed}, Skipped: {genes_skipped}")
        
        # Get aliases for this gene
        aliases = gene_aliases.get(gene_name, []) if gene_aliases else []
        
        # Download bacterial 3D structures for this gene
        gene_summary = download_3d_structures_for_gene(
            gene_name, 
            structures_base_dir, 
            organisms=None,  # Ignored for bacterial search
            gene_aliases=aliases,
            cache=cache
        )
        
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

def main():
    """Main function for Snakemake integration"""
    
    # Check if running under Snakemake
    try:
        # Get parameters from Snakemake
        proteins_file = Path(snakemake.input[0])  # First input is the protein list
        structures_dir = Path("data/proteins_3d_structure")  # Shared directory
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        aliases_file = Path(f"data/quickgo/{paramset}/gene_aliases.json")
        
        # Create sentinel file path
        sentinel_file = Path(snakemake.output[0])
        
    except NameError:
        # Running standalone - use default values for testing
        analysis = "analysis_1"
        paramset = "params_1" 
        group = "positive"
        proteins_file = Path(f"results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv")
        structures_dir = Path("data/proteins_3d_structure")
        aliases_file = Path(f"data/quickgo/{paramset}/gene_aliases.json")
        sentinel_file = Path(f"data/proteins_3d_structure/.{analysis}_{paramset}_{group}_structures_complete")
    
    logging.info(f"Downloading bacterial 3D structures for {analysis}_{paramset}_gram_{group}")
    logging.info(f"Proteins file: {proteins_file}")
    logging.info(f"Structures directory: {structures_dir}")
    
    # Load gene aliases
    gene_aliases = load_gene_aliases(aliases_file)
    
    # Initialize cache
    cache = Protein3DStructureCache()
    logging.info(f"Cache initialized with {len(cache.cache)} entries from {CACHE_FILE}")
    
    # Ensure output directories exist
    structures_dir.mkdir(parents=True, exist_ok=True)
    sentinel_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Process proteins and download bacterial 3D structures
    summary = process_proteins_list(proteins_file, structures_dir, gene_aliases, cache)
    
    # Calculate summary statistics
    total_sequences = sum(gene_data.get("sequences", 0) for gene_data in summary.values())
    total_structures = sum(gene_data.get("structures", 0) for gene_data in summary.values())
    genes_with_3d = len([g for g, data in summary.items() if data.get("found", False)])
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
        "genes_with_3d": genes_with_3d,
        "per_gene_summary": summary
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    # Save cache
    cache.save_cache()
    
    # Create sentinel file to indicate completion
    sentinel_file.touch()
    
    logging.info(f"=== Final Summary ===")
    logging.info(f"Total genes: {len(summary)}")
    logging.info(f"  - Newly processed: {genes_processed}")
    logging.info(f"  - Skipped (already had structures): {genes_skipped}")
    logging.info(f"  - Total completed: {genes_completed}")
    logging.info(f"Downloaded {total_sequences} bacterial sequences and {total_structures} bacterial structures")
    logging.info(f"Genes with 3D structures: {genes_with_3d}/{len(summary)} ({genes_with_3d/len(summary)*100:.1f}%)")
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