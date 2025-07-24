#!/usr/bin/env python3
"""
Bacteria-Specific Protein Information Fetcher for Snakemake
============================================================

This script fetches comprehensive protein information from UniProt
for selected genes, focusing on bacteria (taxid 2) and emphasizing
subcellular location and GO cellular component information.

Works with Snakemake input/output system.
"""

import requests
import json
import time
import pandas as pd
from pathlib import Path
from urllib.parse import urlencode
import sys

def fetch_bacteria_protein_info(gene_name, retries=3, delay=1):
    """
    Fetch comprehensive protein information from UniProt for a gene in bacteria
    
    Args:
        gene_name (str): Gene name (e.g., 'dnaK')
        retries (int): Number of retry attempts
        delay (float): Delay between requests in seconds
    
    Returns:
        dict: Complete protein information or None if not found
    """
    
    # Build UniProt search query - gene name restricted to bacteria (taxid 2)
    query = f'gene_exact:{gene_name} AND taxonomy_id:2'
    
    # Get all data without specifying fields (this works)
    params = {
        'query': query,
        'format': 'json',
        'size': 1  # Get only the first (best) hit
    }
    
    url = f"https://rest.uniprot.org/uniprotkb/search?{urlencode(params)}"
    
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            if data['results']:
                # Return the first (best) match
                entry = data['results'][0]
                
                # Extract comprehensive information
                protein_info = {
                    # Basic identifiers
                    'query_gene': gene_name,
                    'accession': entry.get('primaryAccession', ''),
                    'entry_id': entry.get('uniProtkbId', ''),
                    'reviewed': entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)',
                    'protein_existence': entry.get('proteinExistence', ''),
                    'annotation_score': entry.get('annotationScore', ''),
                    
                    # Gene and protein information
                    'gene_names': extract_gene_names(entry.get('genes', [])),
                    'protein_name': extract_protein_name(entry.get('proteinDescription', {})),
                    'protein_alternative_names': extract_alternative_names(entry.get('proteinDescription', {})),
                    
                    # Organism information
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'organism_id': entry.get('organism', {}).get('taxonId', ''),
                    'taxonomy_lineage': extract_lineage(entry.get('organism', {}).get('lineage', [])),
                    
                    # Sequence information
                    'sequence_length': entry.get('sequence', {}).get('length', ''),
                    'molecular_weight': entry.get('sequence', {}).get('molWeight', ''),
                    'sequence': entry.get('sequence', {}).get('value', ''),
                    
                    # Functional information
                    'function': '',
                    'catalytic_activity': '',
                    'pathway': '',
                    'subunit_structure': '',
                    'subcellular_location': '',
                    'ptm': '',
                    'similarity': '',
                    
                    # GO terms
                    'go_biological_process': '',
                    'go_molecular_function': '',
                    'go_cellular_component': '',
                    'go_terms_count': 0,
                    
                    # Structural features
                    'domains': '',
                    'regions': '',
                    'binding_sites': '',
                    'active_sites': '',
                    
                    # Biochemical properties
                    'ec_numbers': '',
                    'keywords': '',
                    
                    # Cross-references
                    'pdb_structures': '',
                    'kegg_pathways': '',
                    'interpro': '',
                    'pfam': '',
                    
                    # Metadata
                    'retrieval_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'api_url': url
                }
                
                # Extract detailed information
                extract_comments(entry.get('comments', []), protein_info)
                extract_go_terms(entry.get('uniProtKBCrossReferences', []), protein_info)
                extract_features(entry.get('features', []), protein_info)
                extract_keywords(entry.get('keywords', []), protein_info)
                extract_ec_numbers(entry.get('ecNumbers', []), protein_info)
                extract_cross_references(entry.get('uniProtKBCrossReferences', []), protein_info)
                
                return protein_info
                
            else:
                print(f"No UniProt entry found for gene: {gene_name}")
                return None
                
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for {gene_name} (attempt {attempt + 1}): {e}")
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))
            else:
                return None
        except Exception as e:
            print(f"Unexpected error for {gene_name}: {e}")
            return None
    
    return None

def extract_gene_names(genes):
    """Extract all gene names"""
    names = []
    for gene in genes:
        if 'geneName' in gene:
            names.append(gene['geneName'].get('value', ''))
    return '; '.join(filter(None, names))

def extract_protein_name(protein_desc):
    """Extract protein name"""
    if 'recommendedName' in protein_desc:
        return protein_desc['recommendedName'].get('fullName', {}).get('value', '')
    elif 'submissionNames' in protein_desc:
        return protein_desc['submissionNames'][0].get('fullName', {}).get('value', '')
    return ''

def extract_alternative_names(protein_desc):
    """Extract alternative protein names"""
    names = []
    for alt_name in protein_desc.get('alternativeNames', []):
        names.append(alt_name.get('fullName', {}).get('value', ''))
    return '; '.join(filter(None, names))

def extract_lineage(lineage):
    """Extract taxonomic lineage"""
    if isinstance(lineage, list):
        return ' > '.join([item.get('scientificName', '') if isinstance(item, dict) else str(item) for item in lineage])
    else:
        return str(lineage)

def extract_comments(comments, protein_info):
    """Extract functional information from comments"""
    for comment in comments:
        comment_type = comment.get('commentType', '')
        texts = comment.get('texts', [])
        text_value = texts[0].get('value', '') if texts else ''
        
        if comment_type == 'FUNCTION':
            protein_info['function'] = text_value
        elif comment_type == 'CATALYTIC_ACTIVITY':
            protein_info['catalytic_activity'] = text_value
        elif comment_type == 'PATHWAY':
            protein_info['pathway'] = text_value
        elif comment_type == 'SUBUNIT':
            protein_info['subunit_structure'] = text_value
        elif comment_type == 'SUBCELLULAR_LOCATION':
            # Extract subcellular location - can be complex structure
            if 'subcellularLocations' in comment:
                locations = []
                for loc in comment['subcellularLocations']:
                    if 'location' in loc:
                        loc_name = loc['location'].get('value', '')
                        if loc_name:
                            locations.append(loc_name)
                protein_info['subcellular_location'] = '; '.join(locations)
            else:
                protein_info['subcellular_location'] = text_value
        elif comment_type == 'PTM':
            protein_info['ptm'] = text_value
        elif comment_type == 'SIMILARITY':
            protein_info['similarity'] = text_value

def extract_go_terms(cross_refs, protein_info):
    """Extract GO terms from cross-references"""
    go_bp = []
    go_mf = []
    go_cc = []
    
    for xref in cross_refs:
        if xref.get('database') == 'GO':
            go_id = xref.get('id', '')
            properties = xref.get('properties', [])
            go_term = ''
            
            for prop in properties:
                if prop.get('key') == 'GoTerm':
                    go_term = prop.get('value', '')
                    break
            
            if go_term.startswith('P:'):  # Biological Process
                go_bp.append(f"{go_term.replace('P:', '')} ({go_id})")
            elif go_term.startswith('F:'):  # Molecular Function
                go_mf.append(f"{go_term.replace('F:', '')} ({go_id})")
            elif go_term.startswith('C:'):  # Cellular Component
                go_cc.append(f"{go_term.replace('C:', '')} ({go_id})")
    
    protein_info['go_biological_process'] = '; '.join(go_bp)
    protein_info['go_molecular_function'] = '; '.join(go_mf)
    protein_info['go_cellular_component'] = '; '.join(go_cc)
    protein_info['go_terms_count'] = len(go_bp) + len(go_mf) + len(go_cc)

def extract_features(features, protein_info):
    """Extract structural features"""
    domains = []
    regions = []
    binding_sites = []
    active_sites = []
    
    for feature in features:
        feature_type = feature.get('type', '')
        description = feature.get('description', '')
        
        if feature_type == 'DOMAIN':
            domains.append(description)
        elif feature_type == 'REGION':
            regions.append(description)
        elif feature_type == 'BINDING':
            binding_sites.append(description)
        elif feature_type == 'ACT_SITE':
            active_sites.append(description)
    
    protein_info['domains'] = '; '.join(filter(None, domains))
    protein_info['regions'] = '; '.join(filter(None, regions))
    protein_info['binding_sites'] = '; '.join(filter(None, binding_sites))
    protein_info['active_sites'] = '; '.join(filter(None, active_sites))

def extract_keywords(keywords, protein_info):
    """Extract keywords"""
    keyword_list = [kw.get('name', '') for kw in keywords]
    protein_info['keywords'] = '; '.join(keyword_list)

def extract_ec_numbers(ec_numbers, protein_info):
    """Extract EC numbers"""
    ec_list = [ec.get('value', '') for ec in ec_numbers]
    protein_info['ec_numbers'] = '; '.join(ec_list)

def extract_cross_references(cross_refs, protein_info):
    """Extract cross-references"""
    pdb_structures = []
    kegg_pathways = []
    interpro = []
    pfam = []
    
    for xref in cross_refs:
        database = xref.get('database', '')
        ref_id = xref.get('id', '')
        
        if database == 'PDB':
            pdb_structures.append(ref_id)
        elif database == 'KEGG':
            kegg_pathways.append(ref_id)
        elif database == 'InterPro':
            interpro.append(ref_id)
        elif database == 'Pfam':
            pfam.append(ref_id)
    
    protein_info['pdb_structures'] = '; '.join(pdb_structures)
    protein_info['kegg_pathways'] = '; '.join(kegg_pathways)
    protein_info['interpro'] = '; '.join(interpro)
    protein_info['pfam'] = '; '.join(pfam)

def main():
    """Main function for Snakemake"""
    print("Bacteria-Specific Protein Information Fetcher")
    print("=" * 50)
    
    # Get inputs from Snakemake
    try:
        coverage_file = snakemake.input[0]
        output_tsv = snakemake.output.tsv
        output_json = snakemake.output.json
        
        print(f"Input file: {coverage_file}")
        print(f"Output TSV: {output_tsv}")
        print(f"Output JSON: {output_json}")
        
        # Read coverage count file to get gene list
        df = pd.read_csv(coverage_file, sep='\t')
        gene_list = df['gene'].tolist()
        
        print(f"Processing {len(gene_list)} genes...")
        
        # Results storage
        all_results = []
        
        # Process each gene
        for idx, gene_name in enumerate(gene_list):
            print(f"\nProcessing gene {idx + 1}/{len(gene_list)}: {gene_name}")
            
            protein_info = fetch_bacteria_protein_info(gene_name)
            
            if protein_info:
                all_results.append(protein_info)
                print(f"  Found: {protein_info['protein_name']}")
                print(f"  Organism: {protein_info['organism']}")
                print(f"  Subcellular location: {protein_info['subcellular_location']}")
                print(f"  GO cellular component: {protein_info['go_cellular_component']}")
            else:
                print(f"  No information found for {gene_name}")
            
            # Rate limiting
            time.sleep(0.5)
        
        # Save results
        if all_results:
            # Create output directories
            Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
            Path(output_json).parent.mkdir(parents=True, exist_ok=True)
            
            # Save as TSV
            results_df = pd.DataFrame(all_results)
            results_df.to_csv(output_tsv, sep='\t', index=False)
            print(f"\nResults saved to {output_tsv}")
            
            # Save as JSON
            with open(output_json, 'w') as f:
                json.dump(all_results, f, indent=2)
            print(f"Detailed results saved to {output_json}")
            
            # Summary statistics
            print(f"\nSummary:")
            print(f"Total genes processed: {len(gene_list)}")
            print(f"Total protein entries retrieved: {len(all_results)}")
            print(f"Success rate: {len(all_results) / len(gene_list) * 100:.1f}%")
            
            # Show subcellular location and GO cellular component statistics
            subcellular_count = sum(1 for r in all_results if r.get('subcellular_location'))
            go_cc_count = sum(1 for r in all_results if r.get('go_cellular_component'))
            
            print(f"Proteins with subcellular location: {subcellular_count}")
            print(f"Proteins with GO cellular component: {go_cc_count}")
            
        else:
            print("No results retrieved!")
            
    except NameError:
        print("Error: This script is designed to run with Snakemake")
        sys.exit(1)

if __name__ == "__main__":
    main()