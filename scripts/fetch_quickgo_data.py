import requests
import json
import os
import time

# Use Snakemake-provided file paths
goid_file = snakemake.input.go_ids
taxonid_file = snakemake.input.taxon_ids
annotations_output = snakemake.output.annotations
symbols_output = snakemake.output.symbols

# Read GO IDs from TSV file, ignoring lines that start with #
go_ids_list = []

with open(goid_file, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith("gene"):
            continue  # skip empty lines, comments, or headers
        go_id = line.split()[0]
        go_ids_list.append(go_id)

go_ids = ",".join(go_ids_list)

# Read Taxon IDs from file, skipping headers and comments
taxon_ids_list = []

with open(taxonid_file, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#') or line.lower().startswith("taxon"):
            continue  # skip empty lines, comments, or headers
        taxon_id = line.split()[0]
        taxon_ids_list.append(taxon_id)

taxon_ids = ",".join(taxon_ids_list)

# Ensure output directory exists
output_dir = os.path.dirname(annotations_output)
os.makedirs(output_dir, exist_ok=True)

base_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
params = {
    "proteome" : "complete",
    "includeFields": "taxonName",
    "selectedFields": "symbol",
    "geneProductType": "protein",
    "goId": go_ids,
    "goUsage": "exact",
    "goUsageRelationships": "is_a,part_of,occurs_in",
    "taxonId": taxon_ids,
    "taxonUsage": "descendants",
    "aspect": "cellular_component",
    "limit": 200,
    "page": 1
}

headers = {"Accept": "application/json"}

# Step 1: Get total pages using total results instead of relying on pageInfo.totalPages
response = requests.get(base_url, params=params, headers=headers)
if not response.ok:
    print(f"Error {response.status_code}: {response.text}")
    response.raise_for_status()

data = response.json()
total_results = data.get('pageInfo', {}).get('total', 0)
total_pages = (total_results // params["limit"]) + (1 if total_results % params["limit"] else 0)
print(f"Total results: {total_results}")
print(f"Total pages to fetch: {total_pages}")

# Step 2: Download each page and save temporarily with a delay
for page in range(1, total_pages + 1):
    params["page"] = page

    # Adjust limit for the last page
    if page == total_pages and total_results % 200 != 0:
        params["limit"] = total_results % 200
    else:
        params["limit"] = 200

    response = requests.get(base_url, params=params, headers=headers)
    if response.ok:
        data = response.json()
        filename = os.path.join(output_dir, f"annotations_page{page}.json")
        with open(filename, "w", encoding="utf-8") as f:
            json.dump(data.get("results", []), f, indent=4)
        print(f"Data successfully downloaded and saved (page {page}).")
        time.sleep(1)  # Pause to avoid overwhelming the API
    else:
        print(f"Error {response.status_code} on page {page}: {response.text}")
        response.raise_for_status()

# Step 3: Combine all pages into one list
all_results = []
for page in range(1, total_pages + 1):
    filename = os.path.join(output_dir, f"annotations_page{page}.json")
    with open(filename, "r", encoding="utf-8") as f:
        page_data = json.load(f)
        all_results.extend(page_data)
    os.remove(filename)

# Step 4: Save combined results to a single JSON file
with open(annotations_output, "w", encoding="utf-8") as f:
    json.dump(all_results, f, indent=4)

print(f"All data successfully combined and saved to JSON. Total items: {len(all_results)}")

# Step 5: Extract gene names (symbols) from the JSON and save to list file
symbols = sorted({entry["symbol"] for entry in all_results if "symbol" in entry and entry["symbol"]})

with open(symbols_output, "w", encoding="utf-8") as f:
    for symbol in symbols:
        f.write(symbol + "\n")

print(f"Extracted {len(symbols)} unique gene symbols and saved to quickgo_gene_symbols.txt")

# def is_valid_symbol(symbol):
#     return not symbol[0].isupper()

# symbols = sorted({
#     entry["symbol"]
#     for entry in all_results
#     if "symbol" in entry and entry["symbol"] and is_valid_symbol(entry["symbol"])
# })

# with open(symbols_output, "w", encoding="utf-8") as f:
#     for symbol in symbols:
#         f.write(symbol + "\n")

# print(f"Extracted {len(symbols)} unique gene symbols (filtered) and saved to quickgo/gene_symbols.txt")
