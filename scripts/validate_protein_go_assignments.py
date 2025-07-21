#!/usr/bin/env python3
"""
Protein GO Term Assignment Validation
====================================

This script validates GO term assignments by checking if proteins have 
GO cellular component annotations in UniProt. Uses gene aliases if
primary gene search fails. Returns "NA" for proteins with no GO annotations.
"""

import requests
import json
import time
import pandas as pd
from pathlib import Path
from urllib.parse import urlencode
import sys


def load_protein_aliases_json(aliases_json_file):
    """Load protein aliases from the JSON file"""
    try:
        with open(aliases_json_file, 'r') as f:
            protein_aliases = json.load(f)
        
        print(f"Loaded aliases for {len(protein_aliases)} proteins from JSON")
        return protein_aliases
        
    except Exception as e:
        print(f"Error: Could not load protein aliases from {aliases_json_file}: {e}")
        return {}


def validate_protein_go_assignment(protein_name, protein_aliases=None, retries=3, delay=1):
    """
    Check if a protein has GO cellular component annotations in UniProt
    
    Returns:
        dict with validation results
    """
    
    # Try primary protein name first
    result = _try_uniprot_go_search(protein_name, retries, delay)
    if result['has_go_cc']:
        return result
    
    # Try aliases if available and primary search failed
    if protein_aliases:
        for alias in protein_aliases:
            if alias != protein_name:  # Don't retry the same name
                alias_result = _try_uniprot_go_search(alias, retries, delay)
                if alias_result['has_go_cc']:
                    alias_result['query_protein'] = protein_name  # Keep original protein name
                    alias_result['found_via_alias'] = alias
                    return alias_result
    
    # No GO cellular component found
    return {
        'query_protein': protein_name,
        'go_cellular_component': 'NA',
        'has_go_cc': False,
        'tried_aliases': protein_aliases or [],
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
                    'query_protein': search_term,
                    'go_cellular_component': '; '.join(go_cc) if go_cc else 'NA',
                    'has_go_cc': len(go_cc) > 0,
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'accession': entry.get('primaryAccession', ''),
                    'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'status': 'Found in UniProt'
                }
            else:
                return {
                    'query_protein': search_term,
                    'go_cellular_component': 'NA',
                    'has_go_cc': False,
                    'status': 'Not found in UniProt bacteria database'
                }
                
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))
            else:
                return {
                    'query_protein': search_term,
                    'go_cellular_component': 'NA',
                    'has_go_cc': False,
                    'status': 'Network error'
                }
        except Exception as e:
            return {
                'query_protein': search_term,
                'go_cellular_component': 'NA',
                'has_go_cc': False,
                'status': f'Error: {str(e)}'
            }
    
    return {
        'query_protein': search_term,
        'go_cellular_component': 'NA',
        'has_go_cc': False,
        'status': 'Query failed'
    }


def main():
    """Main function for Snakemake integration"""
    print("Protein GO Term Assignment Validation")
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
        output_file = "data/quickgo/params_1/protein_go_validation_report.tsv"
        print("Running in test mode")
    
    # Load protein aliases from JSON
    protein_aliases_dict = load_protein_aliases_json(aliases_json_file)
    
    if not protein_aliases_dict:
        print("No protein aliases loaded. Exiting.")
        return
    
    # Get all proteins from the aliases dictionary
    all_proteins = list(protein_aliases_dict.keys())
    print(f"Found {len(all_proteins)} proteins to validate from aliases JSON")
    
    # Validate GO assignments for each protein and create simple output
    final_results = []
    proteins_with_go = 0
    proteins_found_via_alias = 0
    
    for i, protein in enumerate(all_proteins, 1):
        print(f"Validating protein {i}/{len(all_proteins)}: {protein}", end='... ')
        
        # Get aliases for this protein
        protein_aliases = protein_aliases_dict.get(protein, [])
        
        # Validate GO assignment
        result = validate_protein_go_assignment(protein, protein_aliases)
        
        # Create simple two-column output (always use original protein name)
        final_results.append({
            'protein': protein,  # Always the original protein name
            'go_cellular_component': result['go_cellular_component']
        })
        
        if result['has_go_cc']:
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
    print(f"Total proteins validated: {len(all_proteins)}")
    print(f"Proteins with GO cellular component: {proteins_with_go} ({proteins_with_go/len(all_proteins)*100:.1f}%)")
    print(f"Found via aliases: {proteins_found_via_alias}")
    print(f"No GO annotations (NA): {len(all_proteins) - proteins_with_go}")
    print(f"\nResults saved to: {output_path}")
    print(f"Format: protein\\tgo_cellular_component")


if __name__ == "__main__":
    main()