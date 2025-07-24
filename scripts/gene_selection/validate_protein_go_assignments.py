#!/usr/bin/env python3
"""
Protein GO Term Assignment Validation & UniProt Data Collection
==============================================================

This script validates GO term assignments and collects comprehensive protein
information from UniProt. Uses gene aliases if primary gene search fails.
Includes persistent caching to avoid repeat downloads.

Creates:
- data/uniprot/{paramset}/protein_info.json: Comprehensive protein information
- cache/uniprot/{paramset}/{protein}.json: Individual protein cache files
- GO validation report TSV: Simple protein-GO mapping
"""

import requests
import json
import time
import pandas as pd
from pathlib import Path
from urllib.parse import urlencode
import sys
import os
from datetime import datetime


def load_protein_aliases_tsv(aliases_tsv_file):
    """Load protein aliases from the TSV file"""
    try:
        df = pd.read_csv(aliases_tsv_file, sep='\t')
        
        protein_aliases = {}
        for _, row in df.iterrows():
            gene = row['gene']
            alias_list = row['alias_list']
            
            if pd.isna(alias_list) or alias_list == '':
                protein_aliases[gene] = []
            else:
                # Split comma-separated aliases
                protein_aliases[gene] = [alias.strip() for alias in alias_list.split(',')]
        
        print(f"Loaded aliases for {len(protein_aliases)} proteins from TSV")
        return protein_aliases
        
    except Exception as e:
        print(f"Error: Could not load protein aliases from {aliases_tsv_file}: {e}")
        return {}


def get_comprehensive_protein_info(protein_name, protein_aliases=None, cache_dir=None, retries=3, delay=1):
    """
    Get comprehensive protein information from UniProt with caching
    
    Returns:
        dict with comprehensive protein information
    """
    
    # Check cache first
    if cache_dir:
        cache_file = os.path.join(cache_dir, f"{protein_name}.json")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'r') as f:
                    cached_data = json.load(f)
                return cached_data
            except Exception as e:
                print(f"Error reading cache for {protein_name}: {e}")
    
    # Try primary protein name first
    result = _try_comprehensive_uniprot_search(protein_name, retries, delay)
    if result['found_in_uniprot']:
        # Cache successful result
        if cache_dir:
            _cache_protein_data(cache_dir, protein_name, result)
        return result
    
    # Try aliases if available and primary search failed
    if protein_aliases:
        for alias in protein_aliases:
            if alias != protein_name:  # Don't retry the same name
                alias_result = _try_comprehensive_uniprot_search(alias, retries, delay)
                if alias_result['found_in_uniprot']:
                    alias_result['query_protein'] = protein_name  # Keep original protein name
                    alias_result['found_via_alias'] = alias
                    # Cache successful result
                    if cache_dir:
                        _cache_protein_data(cache_dir, protein_name, alias_result)
                    return alias_result
    
    # No protein found
    not_found_result = {
        'query_protein': protein_name,
        'found_in_uniprot': False,
        'accession': '',
        'protein_name': '',
        'gene_names': [],
        'organism': '',
        'taxonomy_id': '',
        'go_cellular_component': 'NA',
        'go_biological_process': 'NA',
        'go_molecular_function': 'NA',
        'keywords': [],
        'protein_families': [],
        'subcellular_location': [],
        'sequence_length': 0,
        'molecular_weight': 0,
        'tried_aliases': protein_aliases or [],
        'status': 'Not found in UniProt bacteria database',
        'search_date': datetime.now().isoformat()
    }
    
    # Cache negative result too
    if cache_dir:
        _cache_protein_data(cache_dir, protein_name, not_found_result)
    
    return not_found_result


def _cache_protein_data(cache_dir, protein_name, data):
    """Cache protein data to avoid repeat API calls"""
    try:
        os.makedirs(cache_dir, exist_ok=True)
        cache_file = os.path.join(cache_dir, f"{protein_name}.json")
        with open(cache_file, 'w') as f:
            json.dump(data, f, indent=2)
    except Exception as e:
        print(f"Warning: Could not cache data for {protein_name}: {e}")


def _try_comprehensive_uniprot_search(search_term, retries=3, delay=1):
    """Try a comprehensive UniProt search for all protein information"""
    
    # Build UniProt search query - gene name restricted to bacteria (taxid 2)
    query = f'gene_exact:{search_term} AND taxonomy_id:2'
    
    params = {
        'query': query,
        'format': 'json',
        'size': 1,
        'fields': 'accession,protein_name,gene_names,organism_name,organism_id,go,go_c,go_f,go_p,keyword,protein_families,cc_subcellular_location,sequence,mass'
    }
    
    url = f"https://rest.uniprot.org/uniprotkb/search?{urlencode(params)}"
    
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            if data['results']:
                entry = data['results'][0]
                
                # Extract comprehensive protein information
                # GO terms
                go_cc, go_bp, go_mf = _extract_go_terms(entry)
                
                # Gene names
                gene_names = _extract_gene_names(entry)
                
                # Keywords
                keywords = _extract_keywords(entry)
                
                # Protein families
                protein_families = _extract_protein_families(entry)
                
                # Subcellular location
                subcellular_location = _extract_subcellular_location(entry)
                
                # Sequence info
                sequence_info = _extract_sequence_info(entry)
                
                return {
                    'query_protein': search_term,
                    'found_in_uniprot': True,
                    'accession': entry.get('primaryAccession', ''),
                    'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'gene_names': gene_names,
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'taxonomy_id': entry.get('organism', {}).get('taxonId', ''),
                    'go_cellular_component': '; '.join(go_cc) if go_cc else 'NA',
                    'go_biological_process': '; '.join(go_bp) if go_bp else 'NA',
                    'go_molecular_function': '; '.join(go_mf) if go_mf else 'NA',
                    'keywords': keywords,
                    'protein_families': protein_families,
                    'subcellular_location': subcellular_location,
                    'sequence_length': sequence_info.get('length', 0),
                    'molecular_weight': sequence_info.get('mass', 0),
                    'status': 'Found in UniProt',
                    'search_date': datetime.now().isoformat()
                }
            else:
                return {
                    'query_protein': search_term,
                    'found_in_uniprot': False,
                    'status': 'Not found in UniProt bacteria database'
                }
                
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))
            else:
                return {
                    'query_protein': search_term,
                    'found_in_uniprot': False,
                    'status': 'Network error'
                }
        except Exception as e:
            return {
                'query_protein': search_term,
                'found_in_uniprot': False,
                'status': f'Error: {str(e)}'
            }
    
    return {
        'query_protein': search_term,
        'found_in_uniprot': False,
        'status': 'Query failed'
    }


def _extract_go_terms(entry):
    """Extract GO terms from UniProt entry"""
    go_cc, go_bp, go_mf = [], [], []
    
    references = entry.get('uniProtKBCrossReferences', [])
    for ref in references:
        if ref.get('database') == 'GO':
            properties = ref.get('properties', [])
            for prop in properties:
                if prop.get('key') == 'GoTerm':
                    go_term = prop.get('value', '')
                    if go_term.startswith('C:'):  # Cellular Component
                        go_cc.append(go_term.replace('C:', '').strip())
                    elif go_term.startswith('P:'):  # Biological Process
                        go_bp.append(go_term.replace('P:', '').strip())
                    elif go_term.startswith('F:'):  # Molecular Function
                        go_mf.append(go_term.replace('F:', '').strip())
    
    return go_cc, go_bp, go_mf


def _extract_gene_names(entry):
    """Extract gene names from UniProt entry"""
    gene_names = []
    genes = entry.get('genes', [])
    for gene in genes:
        if gene.get('geneName', {}).get('value'):
            gene_names.append(gene['geneName']['value'])
        for synonym in gene.get('synonyms', []):
            if synonym.get('value'):
                gene_names.append(synonym['value'])
    return gene_names


def _extract_keywords(entry):
    """Extract keywords from UniProt entry"""
    keywords = []
    for keyword in entry.get('keywords', []):
        if keyword.get('name'):
            keywords.append(keyword['name'])
    return keywords


def _extract_protein_families(entry):
    """Extract protein families from UniProt entry"""
    families = []
    comments = entry.get('comments', [])
    for comment in comments:
        if comment.get('commentType') == 'SIMILARITY':
            texts = comment.get('texts', [])
            for text in texts:
                if text.get('value'):
                    families.append(text['value'])
    return families


def _extract_subcellular_location(entry):
    """Extract subcellular location from UniProt entry"""
    locations = []
    comments = entry.get('comments', [])
    for comment in comments:
        if comment.get('commentType') == 'SUBCELLULAR_LOCATION':
            subcellular_locations = comment.get('subcellularLocations', [])
            for loc in subcellular_locations:
                location = loc.get('location', {})
                if location.get('value'):
                    locations.append(location['value'])
    return locations


def _extract_sequence_info(entry):
    """Extract sequence information from UniProt entry"""
    sequence = entry.get('sequence', {})
    return {
        'length': sequence.get('length', 0),
        'mass': sequence.get('molWeight', 0)
    }


def main():
    """Main function for Snakemake integration"""
    print("Protein GO Term Assignment Validation & UniProt Data Collection")
    print("=" * 65)
    
    try:
        # Get inputs from Snakemake
        aliases_file = snakemake.input.aliases
        output_file = snakemake.output.validation_report
        
        # Extract paramset from output path for UniProt folder
        output_path_parts = Path(output_file).parts
        paramset = None
        for part in output_path_parts:
            if part.startswith('params_'):
                paramset = part
                break
        
        if not paramset:
            paramset = 'params_1'  # fallback
        
        print(f"Aliases file: {aliases_file}")
        print(f"Output file: {output_file}")
        print(f"Paramset: {paramset}")
        
        # Use the TSV file directly
        aliases_tsv_file = aliases_file
        
        # Create UniProt cache directory in cache folder
        cache_dir = f"cache/uniprot/{paramset}"
        
        # Create UniProt data directory for comprehensive info
        uniprot_data_dir = f"data/uniprot/{paramset}"
        
    except NameError:
        # Test mode
        aliases_tsv_file = "data/quickgo/params_1/gene_aliases.txt"
        output_file = "data/quickgo/params_1/protein_go_validation_report.tsv"
        paramset = "params_1"
        cache_dir = f"cache/uniprot/{paramset}"
        uniprot_data_dir = f"data/uniprot/{paramset}"
        print("Running in test mode")
    
    # Load protein aliases from TSV
    protein_aliases_dict = load_protein_aliases_tsv(aliases_tsv_file)
    
    if not protein_aliases_dict:
        print("No protein aliases loaded. Exiting.")
        return
    
    # Get all proteins from the aliases dictionary
    all_proteins = list(protein_aliases_dict.keys())
    print(f"Found {len(all_proteins)} proteins to process from aliases JSON")
    
    # Create UniProt cache and data directories
    os.makedirs(cache_dir, exist_ok=True)
    os.makedirs(uniprot_data_dir, exist_ok=True)
    print(f"UniProt cache directory: {cache_dir}")
    print(f"UniProt data directory: {uniprot_data_dir}")
    
    # Collect comprehensive protein information
    all_protein_info = {}
    final_results = []
    proteins_with_go = 0
    proteins_found_via_alias = 0
    cache_hits = 0
    
    for i, protein in enumerate(all_proteins, 1):
        print(f"Processing protein {i}/{len(all_proteins)}: {protein}", end='... ')
        
        # Check if we'll use cache
        cached_file = os.path.join(cache_dir, f"{protein}.json")
        will_use_cache = os.path.exists(cached_file)
        
        # Get aliases for this protein
        protein_aliases = protein_aliases_dict.get(protein, [])
        
        # Get comprehensive protein information with caching
        result = get_comprehensive_protein_info(protein, protein_aliases, cache_dir)
        
        # Store comprehensive info
        all_protein_info[protein] = result
        
        # Count cache hits
        if will_use_cache:
            cache_hits += 1
        
        # Create simple two-column output for GO validation report
        final_results.append({
            'protein': protein,  # Always the original protein name
            'go_cellular_component': result['go_cellular_component']
        })
        
        if result.get('found_in_uniprot') and result['go_cellular_component'] != 'NA':
            proteins_with_go += 1
            if result.get('found_via_alias'):
                proteins_found_via_alias += 1
                print(f"✓ (via alias '{result['found_via_alias']}')")  
            else:
                print("✓")
        else:
            if protein_aliases:
                print(f"✗ (tried {len(protein_aliases)} aliases)")
            else:
                print("✗")
        
        # Rate limiting (skip for cached results)
        if not will_use_cache:
            time.sleep(0.2)
    
    # Save comprehensive protein information to JSON in data directory
    uniprot_json_file = os.path.join(uniprot_data_dir, "protein_info.json")
    with open(uniprot_json_file, 'w') as f:
        json.dump(all_protein_info, f, indent=2)
    
    # Create simple two-column dataframe for GO validation
    df = pd.DataFrame(final_results)
    
    # Save GO validation results
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    
    # Print summary
    print(f"\nProcessing Summary:")
    print(f"Total proteins processed: {len(all_proteins)}")
    print(f"Cache hits: {cache_hits} ({cache_hits/len(all_proteins)*100:.1f}%)")
    print(f"Proteins found in UniProt: {sum(1 for p in all_protein_info.values() if p.get('found_in_uniprot'))}")
    print(f"Proteins with GO cellular component: {proteins_with_go} ({proteins_with_go/len(all_proteins)*100:.1f}%)")
    print(f"Found via aliases: {proteins_found_via_alias}")
    print(f"No GO annotations (NA): {len(all_proteins) - proteins_with_go}")
    print(f"\nFiles created:")
    print(f"  GO validation report: {output_path}")
    print(f"  Comprehensive UniProt data: {uniprot_json_file}")
    print(f"  Individual protein caches: {cache_dir}/{{protein}}.json")
    print(f"\nFormat: protein\\tgo_cellular_component")


if __name__ == "__main__":
    main()