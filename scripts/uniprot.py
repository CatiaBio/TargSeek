import requests
import time
import os

# ========== CONFIGURATION ========== #

file_path = "data/all_taxids_manual_check.txt"
GO_TERM = "GO:0005576"  # extracellular region
FORMAT = "tsv"
FIELDS = "accession,gene_names,organism_name,go_id"
OUTPUT_DIR = "results/"  # Save files in current directory
API_URL = "https://rest.uniprot.org/uniprotkb/search"
PAGE_SIZE = 500  # Number of entries per page

# ========== FUNCTIONS ========== #

def read_taxonomy_ids(file_path):
    """Read taxon IDs from file, ignoring comments and empty lines"""
    with open(file_path, "r") as f:
        return [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]

def build_query(go_term, tax_id):
    """Build UniProt query for a single taxon ID"""
    return f"(go:{go_term.split(':')[1]}) AND (taxonomy_id:{tax_id})"

def fetch_uniprot_results_paginated(query, format="tsv", fields=None, page_size=500):
    """Fetch all pages of UniProt results for a given query"""
    all_results = []
    cursor = "*"
    while True:
        params = {
            "query": query,
            "format": format,
            "size": page_size,
            "cursor": cursor
        }
        if fields and format in ["tsv", "json"]:
            params["fields"] = fields

        try:
            response = requests.get(API_URL, params=params)
            if response.status_code == 200:
                data = response.text
                if not data.strip():
                    break  # No more data
                all_results.append(data)
                # Check for 'Link' header for pagination
                link_header = response.headers.get("Link")
                if link_header and 'rel="next"' in link_header:
                    # Extract the next cursor from the Link header
                    next_cursor = link_header.split("cursor=")[-1].split(">")[0]
                    cursor = next_cursor
                    time.sleep(1)  # Be respectful to UniProt servers
                else:
                    break  # No more pages
            elif response.status_code == 500:
                print(f"❌ Internal Server Error for query: {query}")
                break
            else:
                print(f"❌ Error {response.status_code} for query: {query}")
                break
        except requests.exceptions.RequestException as e:
            print(f"❌ Request failed: {e}")
            break
    return "\n".join(all_results)

# ========== MAIN SCRIPT ========== #

taxonomy_ids = read_taxonomy_ids(file_path)

for tax_id in taxonomy_ids:
    query = build_query(GO_TERM, tax_id)
    print(f"🔍 Querying UniProt for taxon {tax_id}...")
    result = fetch_uniprot_results_paginated(query, FORMAT, FIELDS, PAGE_SIZE)

    if result and result.strip():  # Skip empty results
        output_path = os.path.join(OUTPUT_DIR, f"{tax_id}.tsv")
        with open(output_path, "w") as f:
            f.write(result)
        print(f"✅ Saved: {output_path}")
    else:
        print(f"⚠️ No results for taxon {tax_id}")
