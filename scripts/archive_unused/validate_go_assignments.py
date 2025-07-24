#!/usr/bin/env python3
"""
GO Term Assignment Validation
=============================

This script validates GO term assignments by checking if genes have 
GO cellular component annotations in UniProt. Uses gene aliases if
primary gene search fails. Returns "NA" for genes with no GO annotations.
"""

import requests
import json
import time
import pandas as pd
from pathlib import Path
from urllib.parse import urlencode
import sys


def load_gene_aliases_json(aliases_json_file):
    """Load gene aliases from the JSON file"""
    try:
        with open(aliases_json_file, 'r') as f:
            gene_aliases = json.load(f)
        
        print(f"Loaded aliases for {len(gene_aliases)} genes from JSON")
        return gene_aliases
        
    except Exception as e:
        print(f"Error: Could not load gene aliases from {aliases_json_file}: {e}")
        return {}


def validate_gene_go_assignment(gene_name, gene_aliases=None, retries=3, delay=1):
    """
    Check if a gene has GO cellular component annotations in UniProt
    
    Returns:
        dict with validation results
    """
    
    # Try primary gene name first
    result = _try_uniprot_go_search(gene_name, retries, delay)
    if result['has_go_cc']:
        return result
    
    # Try aliases if available and primary search failed
    if gene_aliases:
        for alias in gene_aliases:
            if alias != gene_name:  # Don't retry the same name
                alias_result = _try_uniprot_go_search(alias, retries, delay)
                if alias_result['has_go_cc']:
                    alias_result['query_gene'] = gene_name  # Keep original gene name
                    alias_result['found_via_alias'] = alias
                    return alias_result
    
    # No GO cellular component found
    return {
        'query_gene': gene_name,
        'go_cellular_component': 'NA',
        'has_go_cc': False,
        'tried_aliases': gene_aliases or [],
        'status': 'No GO cellular component found'
    }


def _try_uniprot_go_search(search_term, retries=3, delay=1):
    """Try a single UniProt search for GO cellular component"""
    
    # Build UniProt search query - gene name restricted to bacteria (taxid 2)
    query = f'gene_exact:{search_term} AND taxonomy_id:2'
    
    params = {
        'query': query,
        'format': 'json',
        'size': 1
    }
    
    url = f"https://rest.uniprot.org/uniprotkb/search?{urlencode(params)}"
    
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            if data['results']:
                entry = data['results'][0]
                
                # Extract GO cellular component terms
                go_cc = []
                references = entry.get('uniProtKBCrossReferences', [])
                for ref in references:
                    if ref.get('database') == 'GO':
                        properties = ref.get('properties', [])
                        for prop in properties:
                            if prop.get('key') == 'GoTerm':
                                go_term = prop.get('value', '')
                                if go_term.startswith('C:'):  # Cellular Component
                                    go_cc.append(go_term.replace('C:', '').strip())
                
                return {
                    'query_gene': search_term,
                    'go_cellular_component': '; '.join(go_cc) if go_cc else 'NA',
                    'has_go_cc': len(go_cc) > 0,
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'accession': entry.get('primaryAccession', ''),
                    'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'status': 'Found in UniProt'
                }
            else:
                return {
                    'query_gene': search_term,
                    'go_cellular_component': 'NA',
                    'has_go_cc': False,
                    'status': 'Not found in UniProt bacteria database'
                }
                
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))
            else:
                return {
                    'query_gene': search_term,
                    'go_cellular_component': 'NA',
                    'has_go_cc': False,
                    'status': 'Network error'
                }
        except Exception as e:
            return {
                'query_gene': search_term,
                'go_cellular_component': 'NA',
                'has_go_cc': False,
                'status': f'Error: {str(e)}'
            }
    
    return {
        'query_gene': search_term,
        'go_cellular_component': 'NA',
        'has_go_cc': False,
        'status': 'Query failed'
    }


def main():
    """Main function for Snakemake integration"""
    print("GO Term Assignment Validation")
    print("=" * 40)
    
    try:
        # Get inputs from Snakemake
        aliases_file = snakemake.input.aliases
        output_file = snakemake.output.validation_report
        
        print(f"Aliases file: {aliases_file}")
        print(f"Output file: {output_file}")
        
        # Use JSON file instead of text file
        aliases_json_file = Path(aliases_file).with_suffix('.json')
        
    except NameError:
        # Test mode
        aliases_json_file = "data/quickgo/params_1/gene_aliases.json"
        output_file = "data/quickgo/params_1/go_validation_report.tsv"
        print("Running in test mode")
    
    # Load gene aliases from JSON
    gene_aliases_dict = load_gene_aliases_json(aliases_json_file)
    
    if not gene_aliases_dict:
        print("No gene aliases loaded. Exiting.")
        return
    
    # Get all genes from the aliases dictionary
    all_genes = list(gene_aliases_dict.keys())
    print(f"Found {len(all_genes)} genes to validate from aliases JSON")
    
    # Validate GO assignments for each gene and create simple output
    final_results = []
    genes_with_go = 0
    genes_found_via_alias = 0
    
    for i, gene in enumerate(all_genes, 1):
        print(f"Validating gene {i}/{len(all_genes)}: {gene}", end='... ')
        
        # Get aliases for this gene
        gene_aliases = gene_aliases_dict.get(gene, [])
        
        # Validate GO assignment
        result = validate_gene_go_assignment(gene, gene_aliases)
        
        # Create simple two-column output (always use original gene name)
        final_results.append({
            'gene': gene,  # Always the original gene name
            'go_cellular_component': result['go_cellular_component']
        })
        
        if result['has_go_cc']:
            genes_with_go += 1
            if result.get('found_via_alias'):
                genes_found_via_alias += 1
                print(f"✓ (via alias '{result['found_via_alias']}')")
            else:
                print("✓")
        else:
            if gene_aliases:
                print(f"✗ (tried {len(gene_aliases)} aliases)")
            else:
                print("✗")
        
        # Rate limiting
        time.sleep(0.2)
    
    # Create simple two-column dataframe
    df = pd.DataFrame(final_results)
    
    # Save results in the exact format requested
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    
    # Print summary
    print(f"\nValidation Summary:")
    print(f"Total genes validated: {len(all_genes)}")
    print(f"Genes with GO cellular component: {genes_with_go} ({genes_with_go/len(all_genes)*100:.1f}%)")
    print(f"Found via aliases: {genes_found_via_alias}")
    print(f"No GO annotations (NA): {len(all_genes) - genes_with_go}")
    print(f"\nResults saved to: {output_path}")
    print(f"Format: gene\\tgo_cellular_component")


if __name__ == "__main__":
    main()