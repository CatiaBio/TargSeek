#!/usr/bin/env python3
"""
Final Comprehensive Protein Information Fetcher
===============================================

This script fetches comprehensive protein information from UniProt
for each unique gene, getting all available data from the first (best) hit.

Uses the working approach with comprehensive data extraction.
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
                    
                    # Gene information
                    'gene_names': extract_gene_names(entry.get('genes', [])),
                    'gene_primary': extract_primary_gene_name(entry.get('genes', [])),
                    'gene_synonyms': extract_gene_synonyms(entry.get('genes', [])),
                    
                    # Protein names
                    'protein_name': extract_protein_name(entry.get('proteinDescription', {})),
                    'protein_short_name': extract_protein_short_name(entry.get('proteinDescription', {})),
                    'protein_alternative_names': extract_alternative_names(entry.get('proteinDescription', {})),
                    
                    # Organism information
                    'organism': entry.get('organism', {}).get('scientificName', ''),
                    'organism_common_name': entry.get('organism', {}).get('commonName', ''),
                    'organism_id': entry.get('organism', {}).get('taxonId', ''),
                    'taxonomy_lineage': extract_lineage(entry.get('organism', {}).get('lineage', [])),
                    
                    # Sequence information
                    'sequence_length': entry.get('sequence', {}).get('length', ''),
                    'molecular_weight': entry.get('sequence', {}).get('molWeight', ''),
                    'sequence_checksum': entry.get('sequence', {}).get('crc64', ''),
                    'sequence': entry.get('sequence', {}).get('value', ''),
                    
                    # Dates
                    'created_date': entry.get('entryAudit', {}).get('firstPublicDate', ''),
                    'modified_date': entry.get('entryAudit', {}).get('lastAnnotationUpdateDate', ''),
                    'sequence_modified_date': entry.get('entryAudit', {}).get('lastSequenceUpdateDate', ''),
                    
                    # Comprehensive functional information
                    'function': '',
                    'catalytic_activity': '',
                    'cofactor': '',
                    'enzyme_regulation': '',
                    'biophysicochemical_properties': '',
                    'pathway': '',
                    'subunit_structure': '',
                    'interaction': '',
                    'subcellular_location': '',
                    'tissue_specificity': '',
                    'developmental_stage': '',
                    'induction': '',
                    'domain_info': '',
                    'ptm': '',
                    'disease': '',
                    'disruption_phenotype': '',
                    'allergen': '',
                    'toxic_dose': '',
                    'biotechnology': '',
                    'pharmaceutical': '',
                    'similarity': '',
                    'caution': '',
                    'miscellaneous': '',
                    
                    # GO terms
                    'go_biological_process': '',
                    'go_molecular_function': '',
                    'go_cellular_component': '',
                    'go_terms_count': 0,
                    
                    # Structural features
                    'domains': '',
                    'regions': '',
                    'sites': '',
                    'binding_sites': '',
                    'active_sites': '',
                    'metal_binding': '',
                    'modifications': '',
                    'secondary_structure': '',
                    'transmembrane_regions': '',
                    'topological_domains': '',
                    
                    # Biochemical properties
                    'ec_numbers': '',
                    'keywords': '',
                    'protein_families': '',
                    
                    # Cross-references
                    'pdb_structures': '',
                    'embl_nucleotide': '',
                    'refseq': '',
                    'kegg_pathways': '',
                    'drugbank_targets': '',
                    'interpro': '',
                    'pfam': '',
                    'prosite': '',
                    'reactome': '',
                    'string': '',
                    
                    # Publication and evidence
                    'publication_count': 0,
                    'reference_titles': '',
                    'evidence_codes': '',
                    
                    # Additional metadata
                    'entry_version': entry.get('entryAudit', {}).get('entryVersion', ''),
                    'data_source': 'UniProt',
                    'retrieval_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'api_url': url
                }
                
                # Extract detailed information
                extract_comments(entry.get('comments', []), protein_info)
                extract_go_terms(entry.get('goTerms', []), protein_info)
                extract_features(entry.get('features', []), protein_info)
                extract_keywords(entry.get('keywords', []), protein_info)
                extract_ec_numbers(entry.get('ecNumbers', []), protein_info)
                extract_cross_references(entry.get('uniProtKBCrossReferences', []), protein_info)
                extract_publications(entry.get('references', []), protein_info)
                extract_protein_families(entry.get('proteinDescription', {}), protein_info)
                
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
        for syn in gene.get('synonyms', []):
            names.append(syn.get('value', ''))
    return '; '.join(filter(None, names))

def extract_primary_gene_name(genes):
    """Extract primary gene name"""
    for gene in genes:
        if 'geneName' in gene:
            return gene['geneName'].get('value', '')
    return ''

def extract_gene_synonyms(genes):
    """Extract gene synonyms"""
    synonyms = []
    for gene in genes:
        for syn in gene.get('synonyms', []):
            synonyms.append(syn.get('value', ''))
    return '; '.join(synonyms)

def extract_protein_name(protein_desc):
    """Extract protein name"""
    if 'recommendedName' in protein_desc:
        return protein_desc['recommendedName'].get('fullName', {}).get('value', '')
    elif 'submissionNames' in protein_desc:
        return protein_desc['submissionNames'][0].get('fullName', {}).get('value', '')
    return ''

def extract_protein_short_name(protein_desc):
    """Extract protein short name"""
    if 'recommendedName' in protein_desc:
        short_names = protein_desc['recommendedName'].get('shortNames', [])
        if short_names:
            return short_names[0].get('value', '')
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
        
        # Extract text value
        texts = comment.get('texts', [])
        if texts:
            text_value = texts[0].get('value', '')
        else:
            text_value = ''
        
        # Map comment types to protein_info fields
        if comment_type == 'FUNCTION':
            protein_info['function'] = text_value
        elif comment_type == 'CATALYTIC_ACTIVITY':
            protein_info['catalytic_activity'] = text_value
        elif comment_type == 'COFACTOR':
            protein_info['cofactor'] = text_value
        elif comment_type == 'ENZYME_REGULATION':
            protein_info['enzyme_regulation'] = text_value
        elif comment_type == 'BIOPHYSICOCHEMICAL_PROPERTIES':
            protein_info['biophysicochemical_properties'] = text_value
        elif comment_type == 'PATHWAY':
            protein_info['pathway'] = text_value
        elif comment_type == 'SUBUNIT':
            protein_info['subunit_structure'] = text_value
        elif comment_type == 'INTERACTION':
            protein_info['interaction'] = text_value
        elif comment_type == 'SUBCELLULAR_LOCATION':
            protein_info['subcellular_location'] = text_value
        elif comment_type == 'TISSUE_SPECIFICITY':
            protein_info['tissue_specificity'] = text_value
        elif comment_type == 'DEVELOPMENTAL_STAGE':
            protein_info['developmental_stage'] = text_value
        elif comment_type == 'INDUCTION':
            protein_info['induction'] = text_value
        elif comment_type == 'DOMAIN':
            protein_info['domain_info'] = text_value
        elif comment_type == 'PTM':
            protein_info['ptm'] = text_value
        elif comment_type == 'DISEASE':
            protein_info['disease'] = text_value
        elif comment_type == 'DISRUPTION_PHENOTYPE':
            protein_info['disruption_phenotype'] = text_value
        elif comment_type == 'ALLERGEN':
            protein_info['allergen'] = text_value
        elif comment_type == 'TOXIC_DOSE':
            protein_info['toxic_dose'] = text_value
        elif comment_type == 'BIOTECHNOLOGY':
            protein_info['biotechnology'] = text_value
        elif comment_type == 'PHARMACEUTICAL':
            protein_info['pharmaceutical'] = text_value
        elif comment_type == 'SIMILARITY':
            protein_info['similarity'] = text_value
        elif comment_type == 'CAUTION':
            protein_info['caution'] = text_value
        elif comment_type == 'MISCELLANEOUS':
            protein_info['miscellaneous'] = text_value

def extract_go_terms(go_terms, protein_info):
    """Extract GO terms by category"""
    go_bp = []
    go_mf = []
    go_cc = []
    
    for go_term in go_terms:
        aspect = go_term.get('aspect', '')
        term = f"{go_term.get('term', '')} ({go_term.get('id', '')})"
        
        if aspect == 'P':  # Biological Process
            go_bp.append(term)
        elif aspect == 'F':  # Molecular Function
            go_mf.append(term)
        elif aspect == 'C':  # Cellular Component
            go_cc.append(term)
    
    protein_info['go_biological_process'] = '; '.join(go_bp)
    protein_info['go_molecular_function'] = '; '.join(go_mf)
    protein_info['go_cellular_component'] = '; '.join(go_cc)
    protein_info['go_terms_count'] = len(go_terms)

def extract_features(features, protein_info):
    """Extract structural features"""
    domains = []
    regions = []
    sites = []
    binding_sites = []
    active_sites = []
    metal_binding = []
    modifications = []
    secondary_structure = []
    transmembrane = []
    topological = []
    
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
        elif feature_type in ['HELIX', 'STRAND', 'TURN']:
            secondary_structure.append(f"{feature_type}: {description}")
        elif feature_type == 'TRANSMEM':
            transmembrane.append(description)
        elif feature_type == 'TOPO_DOM':
            topological.append(description)
    
    protein_info['domains'] = '; '.join(filter(None, domains))
    protein_info['regions'] = '; '.join(filter(None, regions))
    protein_info['sites'] = '; '.join(filter(None, sites))
    protein_info['binding_sites'] = '; '.join(filter(None, binding_sites))
    protein_info['active_sites'] = '; '.join(filter(None, active_sites))
    protein_info['metal_binding'] = '; '.join(filter(None, metal_binding))
    protein_info['modifications'] = '; '.join(filter(None, modifications))
    protein_info['secondary_structure'] = '; '.join(filter(None, secondary_structure))
    protein_info['transmembrane_regions'] = '; '.join(filter(None, transmembrane))
    protein_info['topological_domains'] = '; '.join(filter(None, topological))

def extract_keywords(keywords, protein_info):
    """Extract keywords"""
    keyword_list = [kw.get('value', '') for kw in keywords]
    protein_info['keywords'] = '; '.join(keyword_list)

def extract_ec_numbers(ec_numbers, protein_info):
    """Extract EC numbers"""
    ec_list = [ec.get('value', '') for ec in ec_numbers]
    protein_info['ec_numbers'] = '; '.join(ec_list)

def extract_cross_references(cross_refs, protein_info):
    """Extract cross-references"""
    pdb_structures = []
    embl_nucleotide = []
    refseq = []
    kegg_pathways = []
    drugbank_targets = []
    interpro = []
    pfam = []
    prosite = []
    reactome = []
    string_db = []
    
    for xref in cross_refs:
        database = xref.get('database', '')
        ref_id = xref.get('id', '')
        
        if database == 'PDB':
            pdb_structures.append(ref_id)
        elif database == 'EMBL':
            embl_nucleotide.append(ref_id)
        elif database == 'RefSeq':
            refseq.append(ref_id)
        elif database == 'KEGG':
            kegg_pathways.append(ref_id)
        elif database == 'DrugBank':
            drugbank_targets.append(ref_id)
        elif database == 'InterPro':
            interpro.append(ref_id)
        elif database == 'Pfam':
            pfam.append(ref_id)
        elif database == 'PROSITE':
            prosite.append(ref_id)
        elif database == 'Reactome':
            reactome.append(ref_id)
        elif database == 'STRING':
            string_db.append(ref_id)
    
    protein_info['pdb_structures'] = '; '.join(pdb_structures)
    protein_info['embl_nucleotide'] = '; '.join(embl_nucleotide)
    protein_info['refseq'] = '; '.join(refseq)
    protein_info['kegg_pathways'] = '; '.join(kegg_pathways)
    protein_info['drugbank_targets'] = '; '.join(drugbank_targets)
    protein_info['interpro'] = '; '.join(interpro)
    protein_info['pfam'] = '; '.join(pfam)
    protein_info['prosite'] = '; '.join(prosite)
    protein_info['reactome'] = '; '.join(reactome)
    protein_info['string'] = '; '.join(string_db)

def extract_publications(references, protein_info):
    """Extract publication information"""
    if references:
        protein_info['publication_count'] = len(references)
        
        # Extract titles for first few publications
        titles = []
        for ref in references[:3]:  # First 3 publications
            if 'citation' in ref and 'title' in ref['citation']:
                titles.append(ref['citation']['title'])
        
        protein_info['reference_titles'] = '; '.join(titles)

def extract_protein_families(protein_desc, protein_info):
    """Extract protein families"""
    families = []
    
    # Check for protein families in the description
    if 'proteinFamilies' in protein_desc:
        for family in protein_desc['proteinFamilies']:
            families.append(family.get('name', ''))
    
    protein_info['protein_families'] = '; '.join(families)

def process_gene_list(gene_list, output_dir):
    """
    Process a list of genes and fetch UniProt information
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
            print(f"  Found: {protein_info['protein_name']}")
            print(f"  Organism: {protein_info['organism']}")
            print(f"  Accession: {protein_info['accession']}")
            print(f"  Reviewed: {protein_info['reviewed']}")
            print(f"  Length: {protein_info['sequence_length']} aa")
            print(f"  GO Terms: {protein_info['go_terms_count']}")
            print(f"  Keywords: {len(protein_info['keywords'].split(';')) if protein_info['keywords'] else 0}")
            print(f"  PDB Structures: {len(protein_info['pdb_structures'].split(';')) if protein_info['pdb_structures'] else 0}")
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
        
        # Show functional annotation statistics
        function_count = sum(1 for r in all_results if r['function'])
        go_count = sum(1 for r in all_results if r['go_terms_count'] > 0)
        structure_count = sum(1 for r in all_results if r['pdb_structures'])
        
        print(f"\nAnnotation Statistics:")
        print(f"Proteins with function annotation: {function_count}")
        print(f"Proteins with GO terms: {go_count}")
        print(f"Proteins with PDB structures: {structure_count}")
    
    else:
        print("No results retrieved!")

def main():
    """Main function"""
    print("Final Comprehensive Protein Information Fetcher")
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