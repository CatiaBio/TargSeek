import requests
import json
import os
import time

# # Use Snakemake-provided file paths
# goid_file = snakemake.input.go_ids
# taxonid_file = snakemake.input.taxon_ids
# annotations_output = snakemake.output.annotations
# symbols_output = snakemake.output.genes

# Local input files
goid_file = "config/quickgo/go_ids.tsv"
taxonid_file = "config/quickgo/taxon_ids.tsv"

# Local output files
output_dir = "data/ebi_proteins"
os.makedirs(output_dir, exist_ok=True)
annotations_output = os.path.join(output_dir, "protein_annotations.tsv")
symbols_output = os.path.join(output_dir, "protein_symbols.txt")

# Read GO IDs from TSV file, ignoring lines that start with #
go_ids_list = []

with open(goid_file, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith("gene"):
            continue  # skip empty lines, comments, or headers
        go_id = line.split()[0]
        go_ids_list.append(go_id)

print(f"Found {len(go_ids_list)} GO IDs: {go_ids_list}")

# Read Taxon IDs from file, skipping headers and comments
taxon_ids_list = []

with open(taxonid_file, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith("taxon"):
            continue  # skip empty lines, comments, or headers
        taxon_id = line.split()[0]
        taxon_ids_list.append(taxon_id)

print(f"Found {len(taxon_ids_list)} Taxon IDs: {taxon_ids_list}")

# Ensure output directory exists
output_dir = os.path.dirname(annotations_output)
os.makedirs(output_dir, exist_ok=True)

# EBI Proteins API base URL
base_url = "https://www.ebi.ac.uk/proteins/api/proteins"

headers = {"Accept": "application/json"}

# Collect all protein data
all_proteins = []
page_size = 100
offset = 0

# The EBI Proteins API supports searching by taxonomy and GO terms
# We'll search by taxonomy first, then filter by GO terms

for taxon_id in taxon_ids_list:
    print(f"\nSearching for proteins in taxonomy: {taxon_id}")
    
    # Search for proteins in specific taxonomy
    params = {
        "offset": 0,
        "size": page_size,
        "taxid": taxon_id,
    }
    
    page_count = 0
    
    while True:
        params["offset"] = offset
        
        try:
            response = requests.get(base_url, params=params, headers=headers)
            
            if response.status_code == 200:
                data = response.json()
                
                # Handle different response formats
                if isinstance(data, list):
                    proteins = data
                elif isinstance(data, dict):
                    proteins = data.get('results', data.get('proteins', []))
                else:
                    proteins = []
                
                if not proteins:
                    print(f"No more results for taxon:{taxon_id}")
                    break
                
                # Filter proteins by GO terms
                filtered_proteins = []
                for protein in proteins:
                    # Check if protein has any of our target GO terms
                    protein_go_terms = []
                    db_refs = protein.get('dbReferences', [])
                    for ref in db_refs:
                        if ref.get('type') == 'GO':
                            protein_go_terms.append(ref.get('id'))
                    
                    # Check if any of the protein's GO terms match our target GO terms
                    if any(go_term in go_ids_list for go_term in protein_go_terms):
                        filtered_proteins.append(protein)
                
                all_proteins.extend(filtered_proteins)
                page_count += 1
                
                print(f"Retrieved {len(filtered_proteins)} matching proteins (page {page_count})")
                
                # Check if we've reached the end
                if len(proteins) < page_size:
                    break
                
                offset += page_size
                time.sleep(0.1)  # Rate limiting
                
            elif response.status_code == 400:
                print(f"Bad request for taxon:{taxon_id}. This may indicate invalid taxonomy ID.")
                break
            elif response.status_code == 404:
                print(f"No results found for taxon:{taxon_id}")
                break
            else:
                print(f"Error {response.status_code} for taxon:{taxon_id}: {response.text}")
                break
                
        except requests.exceptions.RequestException as e:
            print(f"Request failed for taxon:{taxon_id}: {str(e)}")
            break
    
    # Reset offset for next taxon
    offset = 0

print(f"\nTotal proteins found: {len(all_proteins)}")

# Save all protein data to JSON file
with open(annotations_output, "w", encoding="utf-8") as f:
    json.dump(all_proteins, f, indent=4)

print(f"All protein data saved to {annotations_output}")

# Extract gene symbols/names from the proteins
symbols = set()

for protein in all_proteins:
    # Try different possible fields for gene names/symbols
    gene_list = protein.get('gene', [])
    if gene_list:
        for gene in gene_list:
            gene_name = gene.get('name', {}).get('value', '')
            if gene_name:
                symbols.add(gene_name)
            
            # Also check for synonyms
            gene_synonyms = gene.get('synonyms', [])
            for syn in gene_synonyms:
                if syn.get('value'):
                    symbols.add(syn['value'])
    
    # Try protein names as fallback
    protein_info = protein.get('protein', {})
    if 'recommendedName' in protein_info:
        rec_name = protein_info['recommendedName'].get('fullName', {}).get('value', '')
        if rec_name:
            symbols.add(rec_name)
    
    # Also try submitted names
    submitted_names = protein_info.get('submittedName', [])
    for sub_name in submitted_names:
        name_value = sub_name.get('fullName', {}).get('value', '')
        if name_value:
            symbols.add(name_value)

symbols = sorted(list(symbols))

# Save gene symbols to file
with open(symbols_output, "w", encoding="utf-8") as f:
    for symbol in symbols:
        f.write(symbol + "\n")

print(f"Extracted {len(symbols)} unique gene symbols and saved to {symbols_output}")