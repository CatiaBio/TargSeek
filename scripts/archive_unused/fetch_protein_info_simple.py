#!/usr/bin/env python3
"""
Comprehensive Protein Information Fetcher (Simplified)
======================================================

This script fetches comprehensive protein information from UniProt
for each unique gene, getting all available data from the first (best) hit.

Uses only verified API field names to avoid errors.
"""

import requests
import json
import time
import pandas as pd
from pathlib import Path
from urllib.parse import urlencode
import sys

def fetch_complete_protein_info(gene_name, retries=3, delay=1):
    """
    Fetch comprehensive protein information from UniProt for a gene
    
    Args:
        gene_name (str): Gene name (e.g., 'dnaK')
        retries (int): Number of retry attempts
        delay (float): Delay between requests in seconds
    
    Returns:
        dict: Complete protein information or None if not found
    """
    
    # Build UniProt search query - search for gene in any organism
    query = f'gene_exact:{gene_name}'
    
    # Use only working field names
    params = {
        'query': query,
        'format': 'json',
        'fields': 'accession,id,gene_names,protein_name,organism_name,organism_id,length,mass,sequence,cc_function,cc_pathway,cc_catalytic_activity,cc_subunit,cc_subcellular_location,cc_domain,cc_ptm,cc_disease,cc_disruption_phenotype,cc_similarity,cc_caution,go,go_p,go_f,go_c,ft_domain,ft_region,ft_site,ft_binding,ft_act_site,ft_metal,ft_carbohyd,ft_lipid,ft_mod_res,ec,reviewed,protein_existence',
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
                    'accession': entry.get('primaryAccession', ''),
                    'entry_id': entry.get('uniProtkbId', ''),
                    'gene_names': extract_gene_names(entry.get('genes', [])),
                    'protein_name': extract_protein_name(entry.get('proteinDescription', {})),
                    'protein_alternative_names': extract_alternative_names(entry.get('proteinDescription', {})),
                    
                    # Organism information
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'organism_common_name': entry.get('organism', {}).get('commonName', ''),
                    'organism_id': entry.get('organism', {}).get('taxonId', ''),
                    'taxonomy_lineage': extract_lineage(entry.get('organism', {}).get('lineage', [])),
                    
                    # Sequence information
                    'length': entry.get('sequence', {}).get('length', ''),
                    'mass': entry.get('sequence', {}).get('molWeight', ''),
                    'sequence': entry.get('sequence', {}).get('value', ''),
                    
                    # Quality indicators
                    'reviewed': entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)',
                    'protein_existence': entry.get('proteinExistence', ''),
                    'annotation_score': entry.get('annotationScore', ''),
                    
                    # Functional information
                    'function': '',
                    'catalytic_activity': '',
                    'pathway': '',
                    'subunit_structure': '',
                    'subcellular_location': '',
                    'domain_info': '',
                    'ptm': '',
                    'disease': '',
                    'disruption_phenotype': '',
                    'similarity': '',
                    'caution': '',
                    
                    # GO terms
                    'go_biological_process': '',
                    'go_molecular_function': '',
                    'go_cellular_component': '',
                    
                    # Structural features
                    'domains': '',
                    'regions': '',
                    'sites': '',
                    'binding_sites': '',
                    'active_sites': '',
                    'metal_binding': '',
                    'modifications': '',
                    
                    # Biochemical properties
                    'ec_numbers': '',
                    
                    # Cross-references
                    'pdb_structures': '',
                    'kegg_pathways': '',
                    'go_terms': '',
                    
                    # Publication info
                    'publication_count': 0,
                    'references': ''
                }
                
                # Extract functional information from comments
                extract_comments(entry.get('comments', []), protein_info)
                
                # Extract GO terms
                extract_go_terms(entry.get('goTerms', []), protein_info)
                
                # Extract structural features
                extract_features(entry.get('features', []), protein_info)
                
                # Extract EC numbers
                extract_ec_numbers(entry.get('ecNumbers', []), protein_info)
                
                # Extract cross-references
                extract_cross_references(entry.get('uniProtKBCrossReferences', []), protein_info)
                
                # Extract publications
                extract_publications(entry.get('references', []), protein_info)
                
                # Add query information
                protein_info['query_gene'] = gene_name
                
                return protein_info
                
            else:
                print(f"No UniProt entry found for gene: {gene_name}")
                return None
                
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for {gene_name} (attempt {attempt + 1}): {e}")
            if attempt < retries - 1:
                time.sleep(delay * (2 ** attempt))  # Exponential backoff
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
        for syn in gene.get('synonyms', []):
            names.append(syn.get('value', ''))
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
    return ' > '.join([item.get('scientificName', '') for item in lineage])

def extract_comments(comments, protein_info):
    """Extract functional information from comments"""
    for comment in comments:
        comment_type = comment.get('commentType', '')
        
        if comment_type == 'FUNCTION':
            protein_info['function'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'CATALYTIC_ACTIVITY':
            protein_info['catalytic_activity'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'PATHWAY':
            protein_info['pathway'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'SUBUNIT':
            protein_info['subunit_structure'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'SUBCELLULAR_LOCATION':
            protein_info['subcellular_location'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'DOMAIN':
            protein_info['domain_info'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'PTM':
            protein_info['ptm'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'DISEASE':
            protein_info['disease'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'DISRUPTION_PHENOTYPE':
            protein_info['disruption_phenotype'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'SIMILARITY':
            protein_info['similarity'] = comment.get('texts', [{}])[0].get('value', '')
        elif comment_type == 'CAUTION':
            protein_info['caution'] = comment.get('texts', [{}])[0].get('value', '')

def extract_go_terms(go_terms, protein_info):
    """Extract GO terms by category"""
    go_bp = []
    go_mf = []
    go_cc = []
    all_go = []
    
    for go_term in go_terms:
        aspect = go_term.get('aspect', '')
        term = f"{go_term.get('term', '')} ({go_term.get('id', '')})"
        all_go.append(term)
        
        if aspect == 'P':  # Biological Process
            go_bp.append(term)
        elif aspect == 'F':  # Molecular Function
            go_mf.append(term)
        elif aspect == 'C':  # Cellular Component
            go_cc.append(term)
    
    protein_info['go_biological_process'] = '; '.join(go_bp)
    protein_info['go_molecular_function'] = '; '.join(go_mf)
    protein_info['go_cellular_component'] = '; '.join(go_cc)
    protein_info['go_terms'] = '; '.join(all_go)

def extract_features(features, protein_info):
    """Extract structural features"""
    domains = []
    regions = []
    sites = []
    binding_sites = []
    active_sites = []
    metal_binding = []
    modifications = []
    
    for feature in features:
        feature_type = feature.get('type', '')
        description = feature.get('description', '')
        
        if feature_type == 'DOMAIN':
            domains.append(description)
        elif feature_type == 'REGION':
            regions.append(description)
        elif feature_type == 'SITE':
            sites.append(description)
        elif feature_type == 'BINDING':
            binding_sites.append(description)
        elif feature_type == 'ACT_SITE':
            active_sites.append(description)
        elif feature_type == 'METAL':
            metal_binding.append(description)
        elif feature_type in ['MOD_RES', 'LIPID', 'CARBOHYD']:
            modifications.append(description)
    
    protein_info['domains'] = '; '.join(filter(None, domains))
    protein_info['regions'] = '; '.join(filter(None, regions))
    protein_info['sites'] = '; '.join(filter(None, sites))
    protein_info['binding_sites'] = '; '.join(filter(None, binding_sites))
    protein_info['active_sites'] = '; '.join(filter(None, active_sites))
    protein_info['metal_binding'] = '; '.join(filter(None, metal_binding))
    protein_info['modifications'] = '; '.join(filter(None, modifications))

def extract_ec_numbers(ec_numbers, protein_info):
    """Extract EC numbers"""
    ec_list = [ec.get('value', '') for ec in ec_numbers]
    protein_info['ec_numbers'] = '; '.join(ec_list)

def extract_cross_references(cross_refs, protein_info):
    """Extract cross-references"""
    pdb_structures = []
    kegg_pathways = []
    
    for xref in cross_refs:
        database = xref.get('database', '')
        ref_id = xref.get('id', '')
        
        if database == 'PDB':
            pdb_structures.append(ref_id)
        elif database == 'KEGG':
            kegg_pathways.append(ref_id)
    
    protein_info['pdb_structures'] = '; '.join(pdb_structures)
    protein_info['kegg_pathways'] = '; '.join(kegg_pathways)

def extract_publications(references, protein_info):
    """Extract publication information"""
    if references:
        protein_info['publication_count'] = len(references)
        
        # Extract titles for first few publications
        titles = []
        for ref in references[:5]:  # First 5 publications
            if 'citation' in ref and 'title' in ref['citation']:
                titles.append(ref['citation']['title'])
        
        protein_info['references'] = '; '.join(titles)

def process_gene_list(gene_list, output_dir):
    """
    Process a list of genes and fetch UniProt information
    
    Args:
        gene_list (list): List of gene names
        output_dir (str): Output directory for results
    """
    
    print(f"Processing {len(gene_list)} genes...")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Results storage
    all_results = []
    
    # Process each gene
    for idx, gene_name in enumerate(gene_list):
        print(f"\nProcessing gene {idx + 1}/{len(gene_list)}: {gene_name}")
        
        protein_info = fetch_complete_protein_info(gene_name)
        
        if protein_info:
            all_results.append(protein_info)
            print(f"  Found: {protein_info['protein_name']} ({protein_info['organism']})")
            print(f"  Accession: {protein_info['accession']}")
            print(f"  Reviewed: {protein_info['reviewed']}")
            print(f"  Length: {protein_info['length']} aa")
            print(f"  Function: {protein_info['function'][:100]}..." if protein_info['function'] else "  Function: Not available")
        else:
            print(f"  No information found for {gene_name}")
        
        # Rate limiting
        time.sleep(0.5)
    
    # Save results
    if all_results:
        # Save as TSV
        results_df = pd.DataFrame(all_results)
        tsv_file = output_path / "protein_info_complete.tsv"
        results_df.to_csv(tsv_file, sep='\t', index=False)
        print(f"\nResults saved to {tsv_file}")
        
        # Save as JSON for detailed information
        json_file = output_path / "protein_info_complete.json"
        with open(json_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        print(f"Detailed results saved to {json_file}")
        
        # Summary statistics
        print(f"\nSummary:")
        print(f"Total genes processed: {len(gene_list)}")
        print(f"Total protein entries retrieved: {len(all_results)}")
        print(f"Success rate: {len(all_results) / len(gene_list) * 100:.1f}%")
        
        # Count reviewed vs unreviewed
        reviewed_count = sum(1 for r in all_results if r['reviewed'])
        print(f"Reviewed entries: {reviewed_count}")
        print(f"Unreviewed entries: {len(all_results) - reviewed_count}")
        
        # Show organisms represented
        organisms = set(r['organism'] for r in all_results)
        print(f"Organisms represented: {len(organisms)}")
        for org in sorted(organisms):
            print(f"  - {org}")
    
    else:
        print("No results retrieved!")

def main():
    """Main function"""
    print("Comprehensive Protein Information Fetcher")
    print("="*50)
    
    # Get inputs from Snakemake or use defaults
    try:
        coverage_file = snakemake.input[0]
        output_dir = snakemake.output[0]
        
        print(f"Input file: {coverage_file}")
        print(f"Output directory: {output_dir}")
        
        # Read coverage count file to get gene list
        df = pd.read_csv(coverage_file, sep='\t')
        gene_list = df['gene'].tolist()
        
        process_gene_list(gene_list, output_dir)
        
    except NameError:
        print("Snakemake object not available, using test values")
        
        # Test with gram_positive genes
        coverage_file = "results/coverage/analysis_1_params_1_gram_positive_coverage_count.tsv"
        output_dir = "results/protein_info"
        
        if Path(coverage_file).exists():
            df = pd.read_csv(coverage_file, sep='\t')
            gene_list = df['gene'].tolist()
            process_gene_list(gene_list, output_dir)
        else:
            print(f"Test file not found: {coverage_file}")
            # Use a small test set
            test_genes = ["dnaK", "eno", "era"]
            process_gene_list(test_genes, output_dir)

if __name__ == "__main__":
    main()