#!/usr/bin/env python3
"""
UniProt Information Fetcher with Coverage Enrichment
===================================================

This script:
1. Reads coverage count TSV files
2. Fetches subcellular location and GO cellular component for each gene from UniProt (bacteria only)
3. Creates enriched coverage files with location data
4. Saves detailed protein information
"""

import requests
import json
import time
import pandas as pd
from pathlib import Path
from urllib.parse import urlencode
import sys


def load_gene_aliases(aliases_file):
    """
    Load gene aliases from the aliases file
    
    Returns:
        dict: {gene_name: [list_of_aliases]}
    """
    gene_aliases = {}
    
    try:
        with open(aliases_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    if ':' in line:
                        gene, aliases_str = line.split(':', 1)
                        if aliases_str.strip():
                            aliases = [alias.strip() for alias in aliases_str.split(',')]
                            gene_aliases[gene] = aliases
                        else:
                            gene_aliases[gene] = []
                    
        print(f"Loaded aliases for {len(gene_aliases)} genes")
        return gene_aliases
        
    except Exception as e:
        print(f"Warning: Could not load gene aliases from {aliases_file}: {e}")
        return {}


def fetch_gene_location_info(gene_name, gene_aliases=None, retries=3, delay=1):
    """
    Fetch subcellular location and GO cellular component for a gene from UniProt
    Searches within bacteria (taxon ID 2) and returns first hit
    If no result found, tries aliases
    
    Args:
        gene_name: Primary gene name to search
        gene_aliases: List of aliases to try if primary search fails
        
    Returns:
        dict with location info and protein details
    """
    
    # Try primary gene name first
    result = _try_uniprot_search(gene_name, retries, delay)
    if result['found']:
        return result
    
    # Try aliases if available and primary search failed
    if gene_aliases:
        print(f" (trying aliases: {', '.join(gene_aliases)})", end='')
        for alias in gene_aliases:
            if alias != gene_name:  # Don't retry the same name
                alias_result = _try_uniprot_search(alias, retries, delay)
                if alias_result['found']:
                    alias_result['query_gene'] = gene_name  # Keep original gene name
                    alias_result['found_via_alias'] = alias
                    return alias_result
    
    # No results found even with aliases
    return {
        'query_gene': gene_name,
        'subcellular_location': 'NA',
        'go_cellular_component': 'NA',
        'found': False,
        'tried_aliases': gene_aliases or []
    }

def _try_uniprot_search(search_term, retries=3, delay=1):
    """
    Internal function to try a single UniProt search
    """
    # Build UniProt search query - gene name restricted to bacteria (taxid 2)
    query = f'gene_exact:{search_term} AND taxonomy_id:2'
    
    params = {
        'query': query,
        'format': 'json',
        'size': 1  # Just get the first hit for location info
    }
    
    url = f"https://rest.uniprot.org/uniprotkb/search?{urlencode(params)}"
    
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            if data['results']:
                entry = data['results'][0]
                
                # Extract subcellular location
                subcellular_locations = []
                comments = entry.get('comments', [])
                for comment in comments:
                    if comment.get('commentType') == 'SUBCELLULAR_LOCATION':
                        if 'subcellularLocations' in comment:
                            for loc in comment['subcellularLocations']:
                                if 'location' in loc:
                                    loc_name = loc['location'].get('value', '')
                                    if loc_name:
                                        subcellular_locations.append(loc_name)
                        else:
                            texts = comment.get('texts', [])
                            if texts:
                                subcellular_locations.append(texts[0].get('value', ''))
                
                # Extract GO cellular component (value only, no GO ID)
                go_cc = []
                cross_refs = entry.get('uniProtKBCrossReferences', [])
                for xref in cross_refs:
                    if xref.get('database') == 'GO':
                        properties = xref.get('properties', [])
                        for prop in properties:
                            if prop.get('key') == 'GoTerm':
                                go_term = prop.get('value', '')
                                if go_term.startswith('C:'):  # Cellular Component
                                    go_cc.append(go_term.replace('C:', '').strip())
                
                # Extract other useful information
                protein_info = {
                    'query_gene': search_term,
                    'subcellular_location': '; '.join(subcellular_locations) if subcellular_locations else 'Requires manual check',
                    'go_cellular_component': '; '.join(go_cc) if go_cc else 'Requires manual check',
                    'found': True,
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'organism_id': entry.get('organism', {}).get('taxonId', ''),
                    'accession': entry.get('primaryAccession', ''),
                    'entry_id': entry.get('uniProtkbId', ''),
                    'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'gene_names': ', '.join([gn.get('geneName', {}).get('value', '') for gn in entry.get('genes', [])]),
                    'length': entry.get('sequence', {}).get('length', ''),
                    'mass': entry.get('sequence', {}).get('molWeight', ''),
                    'reviewed': entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)',
                    'existence': entry.get('proteinExistence', ''),
                    'function': '',
                    'ec_numbers': ''
                }
                
                # Extract function
                for comment in comments:
                    if comment.get('commentType') == 'FUNCTION':
                        protein_info['function'] = comment.get('texts', [{}])[0].get('value', '')
                
                # Extract EC numbers
                ec_numbers = entry.get('ecNumbers', [])
                protein_info['ec_numbers'] = '; '.join([ec.get('value', '') for ec in ec_numbers])
                
                return protein_info
                
            else:
                return {
                    'query_gene': search_term,
                    'subcellular_location': 'Not found in UniProt bacteria database',
                    'go_cellular_component': 'Not found in UniProt bacteria database',
                    'found': False
                }
                
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for {search_term} (attempt {attempt + 1}): {e}")
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))
            else:
                return {
                    'query_gene': search_term,
                    'subcellular_location': 'Network error - requires manual check',
                    'go_cellular_component': 'Network error - requires manual check',
                    'found': False
                }
        except Exception as e:
            print(f"Unexpected error for {search_term}: {e}")
            return {
                'query_gene': search_term,
                'subcellular_location': 'Error - requires manual check',
                'go_cellular_component': 'Error - requires manual check',
                'found': False
            }
    
    return {
        'query_gene': search_term,
        'subcellular_location': 'Error - requires manual check',
        'go_cellular_component': 'Error - requires manual check',
        'found': False
    }


def fetch_multiple_hits_per_gene(gene_name, max_hits=3, retries=3, delay=1):
    """
    Fetch up to max_hits UniProt entries for a bacterial gene
    Used for detailed protein information collection
    """
    query = f'gene_exact:{gene_name} AND taxonomy_id:2'
    params = {
        'query': query,
        'format': 'json',
        'size': max_hits
    }
    url = f"https://rest.uniprot.org/uniprotkb/search?{urlencode(params)}"
    
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            results = []
            for entry in data.get('results', []):
                info = {
                    'query_gene': gene_name,
                    'accession': entry.get('primaryAccession', ''),
                    'entry_id': entry.get('uniProtkbId', ''),
                    'gene_names': ', '.join([gn.get('geneName', {}).get('value', '') for gn in entry.get('genes', [])]),
                    'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'organism_id': entry.get('organism', {}).get('taxonId', ''),
                    'length': entry.get('sequence', {}).get('length', ''),
                    'mass': entry.get('sequence', {}).get('molWeight', ''),
                    'reviewed': entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)',
                    'existence': entry.get('proteinExistence', ''),
                    'subcellular_location': '',
                    'go_cellular_component': '',
                    'go_biological_process': '',
                    'go_molecular_function': '',
                    'function': '',
                    'ec_numbers': '',
                    'keywords': '',
                    'retrieval_date': time.strftime('%Y-%m-%d %H:%M:%S')
                }
                
                # Extract subcellular location
                subcellular = []
                for comment in entry.get('comments', []):
                    if comment.get('commentType') == 'SUBCELLULAR_LOCATION':
                        if 'subcellularLocations' in comment:
                            for loc in comment['subcellularLocations']:
                                if 'location' in loc:
                                    subcellular.append(loc['location'].get('value', ''))
                        elif comment.get('texts'):
                            subcellular.append(comment['texts'][0].get('value', ''))
                    elif comment.get('commentType') == 'FUNCTION':
                        info['function'] = comment.get('texts', [{}])[0].get('value', '')
                
                info['subcellular_location'] = '; '.join(set(subcellular))
                
                # Extract GO terms (values only, no GO IDs)
                go_bp = []
                go_mf = []
                go_cc = []
                
                for xref in entry.get('uniProtKBCrossReferences', []):
                    if xref.get('database') == 'GO':
                        props = xref.get('properties', [])
                        for prop in props:
                            if prop.get('key') == 'GoTerm':
                                term = prop.get('value', '')
                                if term.startswith('C:'):
                                    go_cc.append(term.replace('C:', '').strip())
                                elif term.startswith('P:'):
                                    go_bp.append(term.replace('P:', '').strip())
                                elif term.startswith('F:'):
                                    go_mf.append(term.replace('F:', '').strip())
                
                info['go_cellular_component'] = '; '.join(go_cc)
                info['go_biological_process'] = '; '.join(go_bp)
                info['go_molecular_function'] = '; '.join(go_mf)
                
                # Extract EC numbers
                ec_numbers = entry.get('ecNumbers', [])
                info['ec_numbers'] = '; '.join([ec.get('value', '') for ec in ec_numbers])
                
                # Extract keywords
                keywords = entry.get('keywords', [])
                info['keywords'] = '; '.join([kw.get('name', '') for kw in keywords])
                
                results.append(info)
            
            return results
            
        except requests.exceptions.RequestException as e:
            print(f"Error fetching multiple hits for {gene_name} (attempt {attempt + 1}): {e}")
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))
            else:
                return []
        except Exception as e:
            print(f"Unexpected error for {gene_name}: {e}")
            return []
    
    return []


def process_coverage_file(coverage_file, output_dir, gene_aliases_dict=None):
    """
    Process a coverage count file:
    1. Add location columns to create enriched coverage file
    2. Collect detailed protein information
    3. Use gene aliases if primary gene search fails
    
    Args:
        coverage_file: Path to coverage TSV file
        output_dir: Output directory 
        gene_aliases_dict: Dictionary of gene aliases {gene: [aliases]}
    """
    
    print(f"\nProcessing {coverage_file}...")
    
    # Read coverage data
    df = pd.read_csv(coverage_file, sep='\t')
    print(f"Found {len(df)} genes to process")
    
    # Initialize new column
    df['go_cellular_component'] = ''
    
    # Store detailed protein info
    all_detailed_info = []
    
    # Process each gene
    for idx, row in df.iterrows():
        gene_name = row['gene']
        print(f"Processing gene {idx + 1}/{len(df)}: {gene_name}", end='... ')
        
        # Get aliases for this gene
        gene_aliases = gene_aliases_dict.get(gene_name, []) if gene_aliases_dict else []
        
        # Fetch location information (first hit only for enrichment)
        location_info = fetch_gene_location_info(gene_name, gene_aliases)
        
        # Update dataframe with GO cellular component info only
        df.at[idx, 'go_cellular_component'] = location_info['go_cellular_component']
        
        # Log if found via alias
        if location_info.get('found_via_alias'):
            print(f" -> found via alias '{location_info['found_via_alias']}'", end='')
        elif not location_info['found'] and gene_aliases:
            print(f" -> not found (tried {len(gene_aliases)} aliases)", end='')
        
        # Get multiple hits for detailed info
        detailed_hits = fetch_multiple_hits_per_gene(gene_name, max_hits=3)
        
        if detailed_hits:
            # Add species count info to each hit
            for hit in detailed_hits:
                hit['species_count_in_analysis'] = row['count']
                hit['species_list'] = row['species']
            all_detailed_info.extend(detailed_hits)
            print(f"Found {len(detailed_hits)} entries")
        else:
            print("Not found")
        
        # Rate limiting
        time.sleep(0.5)
    
    # Reorder columns to put GO cellular component after gene
    cols = df.columns.tolist()
    cols.remove('go_cellular_component')
    gene_idx = cols.index('gene')
    cols = cols[:gene_idx+1] + ['go_cellular_component'] + cols[gene_idx+1:]
    df = df[cols]
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save enriched coverage file
    input_stem = Path(coverage_file).stem
    enriched_file = output_path / f"{input_stem}_location.tsv"
    df.to_csv(enriched_file, sep='\t', index=False)
    print(f"\nEnriched coverage file saved to: {enriched_file}")
    
    # Save detailed protein information
    if all_detailed_info:
        detailed_df = pd.DataFrame(all_detailed_info)
        
        # Save as TSV
        detailed_tsv = output_path / f"{input_stem}_detailed_proteins.tsv"
        detailed_df.to_csv(detailed_tsv, sep='\t', index=False)
        print(f"Detailed protein info saved to: {detailed_tsv}")
        
        # Save as JSON
        detailed_json = output_path / f"{input_stem}_detailed_proteins.json"
        with open(detailed_json, 'w') as f:
            json.dump(all_detailed_info, f, indent=2)
        print(f"Detailed JSON saved to: {detailed_json}")
        
        # Print summary statistics
        print(f"\nSummary:")
        print(f"Total genes processed: {len(df)}")
        print(f"Genes with UniProt data: {sum(1 for _, row in df.iterrows() if row['go_cellular_component'] != 'Not found in UniProt bacteria database')}")
        print(f"Total protein entries retrieved: {len(all_detailed_info)}")
        
        # Count annotations
        with_subcellular = sum(1 for item in all_detailed_info if item['subcellular_location'])
        with_go_cc = sum(1 for item in all_detailed_info if item['go_cellular_component'])
        reviewed = sum(1 for item in all_detailed_info if item['reviewed'])
        
        print(f"Entries with subcellular location: {with_subcellular}")
        print(f"Entries with GO cellular component: {with_go_cc}")
        print(f"Reviewed (Swiss-Prot) entries: {reviewed}")


def main():
    """Main function for Snakemake integration"""
    print("UniProt Information Fetcher with Coverage Enrichment")
    print("="*60)
    
    # Get inputs from Snakemake
    try:
        input_file = snakemake.input.coverage
        aliases_file = snakemake.input.aliases
        output_dir = snakemake.output[0]
        
        print(f"Input file: {input_file}")
        print(f"Aliases file: {aliases_file}")
        print(f"Output directory: {output_dir}")
        
        # Load gene aliases
        gene_aliases_dict = load_gene_aliases(aliases_file)
        print(f"Loaded aliases for {len(gene_aliases_dict)} genes")
        
        # Process the coverage file with aliases
        process_coverage_file(input_file, output_dir, gene_aliases_dict)
        
    except NameError:
        print("Running in test mode (no Snakemake)")
        
        # Test with example files
        test_files = [
            "results/coverage/analysis_1_params_1_gram_positive_coverage_count.tsv",
            "results/coverage/analysis_1_params_1_gram_negative_coverage_count.tsv"
        ]
        
        # Try to load aliases for testing
        aliases_file = "data/quickgo/params_1/gene_aliases.txt"
        gene_aliases_dict = load_gene_aliases(aliases_file) if Path(aliases_file).exists() else {}
        
        for test_file in test_files:
            if Path(test_file).exists():
                output_dir = "results/uniprot_info/test"
                process_coverage_file(test_file, output_dir, gene_aliases_dict)
                break
        else:
            print("No test files found")
            print("Available coverage files:")
            for file in Path("results/coverage").glob("*coverage_count.tsv"):
                print(f"  - {file}")


if __name__ == "__main__":
    main()