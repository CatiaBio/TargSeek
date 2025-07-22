import json
import requests
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# -------------------------
# Snakemake inputs/outputs
# -------------------------
params_file = snakemake.input.params_file
annotations_output = snakemake.output.annotations
symbols_output = snakemake.output.genes
output_dir = os.path.dirname(annotations_output)
os.makedirs(output_dir, exist_ok=True)

# -------------------------
# Load query parameters
# -------------------------
with open(params_file, "r", encoding="utf-8") as f:
    params = json.load(f)

# -------------------------
# Fetch total pages from QuickGO
# -------------------------
base_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
headers = {"Accept": "application/json"}

response = requests.get(base_url, params=params, headers=headers)
response.raise_for_status()
data = response.json()
total_results = data.get("pageInfo", {}).get("total", 0)
limit = int(params.get("limit", 200))
total_pages = (total_results // limit) + (1 if total_results % limit else 0)

# -------------------------
# Function to fetch a single page
# -------------------------
def fetch_page(page_num):
    page_params = params.copy()
    page_params["page"] = page_num
    page_params["limit"] = (
        total_results % limit if page_num == total_pages and total_results % limit != 0 else limit
    )
    
    response = requests.get(base_url, params=page_params, headers=headers)
    response.raise_for_status()
    return response.json().get("results", [])

# -------------------------
# Download annotations concurrently
# -------------------------
all_results = []
max_workers = min(5, total_pages)  # Limit concurrent requests to avoid overwhelming API

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    # Submit all page requests
    future_to_page = {executor.submit(fetch_page, page): page for page in range(1, total_pages + 1)}
    
    # Collect results as they complete
    for future in as_completed(future_to_page):
        page_num = future_to_page[future]
        try:
            page_data = future.result()
            all_results.extend(page_data)
            time.sleep(0.2)  # Reduced sleep to prevent API overload
        except Exception as e:
            print(f"Error fetching page {page_num}: {e}")
            raise

# -------------------------
# Save results and gene symbols
# -------------------------
with open(annotations_output, "w", encoding="utf-8") as f:
    json.dump(all_results, f, indent=4)

# Extract unique symbols more efficiently
symbols = set()
for entry in all_results:
    if "symbol" in entry:
        symbols.add(entry["symbol"])

# Write sorted symbols
with open(symbols_output, "w", encoding="utf-8") as f:
    for symbol in sorted(symbols):
        f.write(symbol + "\n")