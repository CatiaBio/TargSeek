#!/usr/bin/env python3
"""
Enrich Coverage Files with Subcellular Location Data
===================================================

This script:
1. Reads coverage count TSV files
2. Fetches subcellular location and GO cellular component for each gene from UniProt
3. Adds these as new columns to the coverage files
4. Saves enriched files with "_location" suffix
5. Saves detailed protein information in separate files
"""

import requests
import json
import time
import pandas as pd
from pathlib import Path
from urllib.parse import urlencode
import sys


def fetch_gene_location_info(gene_name, retries=3, delay=1):
    """
    Fetch subcellular location and GO cellular component for a gene from UniProt
    Searches within bacteria (taxon ID 2) and returns first hit
    
    Returns:
        dict with 'subcellular_location' and 'go_cellular_component' keys
    """
    
    # Build UniProt search query - gene name restricted to bacteria (taxid 2)
    query = f'gene_exact:{gene_name} AND taxonomy_id:2'
    
    params = {
        'query': query,
        'format': 'json',
        'size': 1  # Just get the first hit
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
                
                # Extract GO cellular component
                go_cc = []
                cross_refs = entry.get('uniProtKBCrossReferences', [])
                for xref in cross_refs:
                    if xref.get('database') == 'GO':
                        go_id = xref.get('id', '')
                        properties = xref.get('properties', [])
                        for prop in properties:
                            if prop.get('key') == 'GoTerm':
                                go_term = prop.get('value', '')
                                if go_term.startswith('C:'):  # Cellular Component
                                    go_cc.append(f"{go_term.replace('C:', '')} ({go_id})")
                
                return {
                    'subcellular_location': '; '.join(subcellular_locations) if subcellular_locations else 'Requires manual check',
                    'go_cellular_component': '; '.join(go_cc) if go_cc else 'Requires manual check',
                    'found': True,
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'accession': entry.get('primaryAccession', ''),
                    'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'reviewed': entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)'
                }
            else:
                return {
                    'subcellular_location': 'Not found in UniProt bacteria database',
                    'go_cellular_component': 'Not found in UniProt bacteria database',
                    'found': False
                }
                
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for {gene_name} (attempt {attempt + 1}): {e}")
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))
            else:
                return {
                    'subcellular_location': 'Network error - requires manual check',
                    'go_cellular_component': 'Network error - requires manual check',
                    'found': False
                }
        except Exception as e:
            print(f"Unexpected error for {gene_name}: {e}")
            return {
                'subcellular_location': 'Error - requires manual check',
                'go_cellular_component': 'Error - requires manual check',
                'found': False
            }
    
    return {
        'subcellular_location': 'Error - requires manual check',
        'go_cellular_component': 'Error - requires manual check',
        'found': False
    }


def enrich_coverage_file(input_file, output_file, detailed_output_file):
    """
    Process a coverage count file and add location information
    
    Args:
        input_file: Path to input coverage count TSV
        output_file: Path to output enriched TSV
        detailed_output_file: Path to save detailed protein information
    """
    
    print(f"\nProcessing {input_file}...")
    
    # Read coverage data
    df = pd.read_csv(input_file, sep='\t')
    print(f"Found {len(df)} genes to enrich")
    
    # Initialize new columns
    df['subcellular_location'] = ''
    df['go_cellular_component'] = ''
    
    # Store detailed protein info
    detailed_info = []
    
    # Process each gene
    for idx, row in df.iterrows():
        gene_name = row['gene']
        print(f"Processing gene {idx + 1}/{len(df)}: {gene_name}", end='... ')
        
        # Fetch location information
        location_info = fetch_gene_location_info(gene_name)
        
        # Update dataframe
        df.at[idx, 'subcellular_location'] = location_info['subcellular_location']
        df.at[idx, 'go_cellular_component'] = location_info['go_cellular_component']
        
        # Store detailed info if found
        if location_info['found']:
            detailed_info.append({
                'gene': gene_name,
                'species_count': row['count'],
                'subcellular_location': location_info['subcellular_location'],
                'go_cellular_component': location_info['go_cellular_component'],
                'organism': location_info['organism'],
                'accession': location_info['accession'],
                'protein_name': location_info['protein_name'],
                'reviewed': location_info['reviewed']
            })
            print("Found")
        else:
            print("Not found")
        
        # Rate limiting
        time.sleep(0.5)
    
    # Reorder columns to put location info after gene
    cols = df.columns.tolist()
    # Remove the location columns from their current position
    cols.remove('subcellular_location')
    cols.remove('go_cellular_component')
    # Insert after 'gene' column
    gene_idx = cols.index('gene')
    cols = cols[:gene_idx+1] + ['subcellular_location', 'go_cellular_component'] + cols[gene_idx+1:]
    df = df[cols]
    
    # Save enriched coverage file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nEnriched coverage file saved to: {output_file}")
    
    # Save detailed protein information
    if detailed_info:
        detailed_df = pd.DataFrame(detailed_info)
        detailed_df.to_csv(detailed_output_file, sep='\t', index=False)
        print(f"Detailed protein info saved to: {detailed_output_file}")
        
        # Print summary statistics
        found_count = len(detailed_info)
        print(f"\nSummary:")
        print(f"Genes with UniProt data: {found_count}/{len(df)} ({found_count/len(df)*100:.1f}%)")
        
        # Count genes with location data
        with_subcellular = sum(1 for item in detailed_info if item['subcellular_location'] != 'Requires manual check')
        with_go_cc = sum(1 for item in detailed_info if item['go_cellular_component'] != 'Requires manual check')
        reviewed = sum(1 for item in detailed_info if item['reviewed'])
        
        print(f"Genes with subcellular location: {with_subcellular}")
        print(f"Genes with GO cellular component: {with_go_cc}")
        print(f"Reviewed (Swiss-Prot) entries: {reviewed}")


def main():
    """Main function"""
    print("Coverage File Location Enrichment")
    print("="*50)
    
    # Define input and output files
    files_to_process = [
        {
            'input': 'results/coverage/analysis_1_params_1_gram_positive_coverage_count.tsv',
            'output': 'results/coverage/analysis_1_params_1_gram_positive_coverage_count_location.tsv',
            'detailed': 'results/protein_info/gram_positive_protein_details.tsv'
        },
        {
            'input': 'results/coverage/analysis_1_params_1_gram_negative_coverage_count.tsv',
            'output': 'results/coverage/analysis_1_params_1_gram_negative_coverage_count_location.tsv',
            'detailed': 'results/protein_info/gram_negative_protein_details.tsv'
        }
    ]
    
    # Create output directory for detailed info
    Path('results/protein_info').mkdir(parents=True, exist_ok=True)
    
    # Process each file
    for file_info in files_to_process:
        if Path(file_info['input']).exists():
            enrich_coverage_file(
                file_info['input'],
                file_info['output'],
                file_info['detailed']
            )
        else:
            print(f"\nWarning: Input file not found: {file_info['input']}")
    
    print("\nProcessing complete!")


if __name__ == "__main__":
    main()