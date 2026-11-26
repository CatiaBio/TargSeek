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

# QuickGO occasionally returns 500/502/503 for a few minutes at a time. A single
# transient error used to abort the whole download pipeline, so retry with
# exponential backoff before giving up.
RETRY_STATUS = {429, 500, 502, 503, 504}
MAX_ATTEMPTS = 5
BACKOFF_SECONDS = 5


def quickgo_get(query_params, what="request"):
    """GET from QuickGO, retrying transient server errors."""
    last_error = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            response = requests.get(base_url, params=query_params,
                                    headers=headers, timeout=120)
            if response.status_code in RETRY_STATUS and attempt < MAX_ATTEMPTS:
                wait = BACKOFF_SECONDS * (2 ** (attempt - 1))
                print(f"QuickGO {what}: HTTP {response.status_code}, "
                      f"retrying in {wait}s (attempt {attempt}/{MAX_ATTEMPTS})")
                time.sleep(wait)
                continue
            response.raise_for_status()
            return response.json()
        except (requests.exceptions.Timeout,
                requests.exceptions.ConnectionError) as e:
            last_error = e
            if attempt < MAX_ATTEMPTS:
                wait = BACKOFF_SECONDS * (2 ** (attempt - 1))
                print(f"QuickGO {what}: {type(e).__name__}, retrying in {wait}s "
                      f"(attempt {attempt}/{MAX_ATTEMPTS})")
                time.sleep(wait)
                continue
            raise
    if last_error:
        raise last_error
    raise RuntimeError(f"QuickGO {what} failed after {MAX_ATTEMPTS} attempts")


data = quickgo_get(params, what="page count")
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
    
    return quickgo_get(page_params, what=f"page {page_num}").get("results", [])

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