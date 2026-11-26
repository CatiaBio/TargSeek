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

# QuickGO's pageInfo.total is the number of PAGES, not the number of results.
# (Verified: the same query reports total=1903 at limit=1, 381 at limit=5 and
# 10 at limit=200 for ~1903 actual results.) This code used to treat it as a
# result count and divide by limit, so it fetched only total/limit pages -
# about 0.5% of the available annotations.
page_info = data.get("pageInfo", {})
limit = int(params.get("limit", 200))
available_pages = int(page_info.get("total", 0))

# Fetching every page of an all-bacteria query means millions of annotations,
# so cap it. The cap is recorded in the metadata file next to the output.
max_pages_cfg = params.pop("maxPages", None)
max_pages = int(max_pages_cfg) if max_pages_cfg else 0
total_pages = min(available_pages, max_pages) if max_pages else available_pages
truncated = total_pages < available_pages

print(f"QuickGO reports {available_pages:,} pages at limit={limit} "
      f"(~{available_pages * limit:,} annotations)")
if truncated:
    print(f"WARNING: fetching only the first {total_pages:,} pages "
          f"(~{total_pages * limit:,} annotations) because maxPages={max_pages} "
          f"is set in the parameter file. Results are a subset, not the full set.")
else:
    print(f"Fetching all {total_pages:,} pages")

# -------------------------
# Function to fetch a single page
# -------------------------
def fetch_page(page_num):
    page_params = params.copy()
    page_params["page"] = page_num
    page_params["limit"] = limit  # the last page simply returns fewer rows
    return quickgo_get(page_params, what=f"page {page_num}").get("results", [])

# -------------------------
# Download annotations concurrently
# -------------------------
all_results = []
max_workers = min(5, total_pages)
# 0.2s per page would add hours to a 12k-page fetch; concurrency is already
# capped at 5 in-flight requests, which is the real politeness control.
sleep_between = 0.2 if total_pages <= 200 else 0.02  # Limit concurrent requests to avoid overwhelming API

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    # Submit all page requests
    future_to_page = {executor.submit(fetch_page, page): page for page in range(1, total_pages + 1)}
    
    # Collect results as they complete
    done = 0
    for future in as_completed(future_to_page):
        page_num = future_to_page[future]
        try:
            page_data = future.result()
            all_results.extend(page_data)
            done += 1
            if done % 100 == 0 or done == total_pages:
                print(f"  fetched {done:,}/{total_pages:,} pages "
                      f"({len(all_results):,} annotations)", flush=True)
            time.sleep(sleep_between)
        except Exception as e:
            print(f"Error fetching page {page_num}: {e}")
            raise

# -------------------------
# Save results and gene symbols
# -------------------------
# indent=4 on millions of entries wastes a lot of disk; keep it compact
with open(annotations_output, "w", encoding="utf-8") as f:
    json.dump(all_results, f)

# Extract unique symbols more efficiently
symbols = set()
for entry in all_results:
    if "symbol" in entry:
        symbols.add(entry["symbol"])

# Write sorted symbols
with open(symbols_output, "w", encoding="utf-8") as f:
    for symbol in sorted(symbols):
        f.write(symbol + "\n")

# Record exactly what was fetched, so a truncated run is never mistaken for a
# complete one when the results are interpreted later.
taxa = {e.get("taxonId") for e in all_results if e.get("taxonId")}
metadata = {
    "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    "query": {k: v for k, v in params.items() if k != "page"},
    "pages_available": available_pages,
    "pages_fetched": total_pages,
    "max_pages_setting": max_pages or None,
    "truncated": truncated,
    "annotations_fetched": len(all_results),
    "distinct_symbols": len(symbols),
    "distinct_taxa": len(taxa),
}
with open(os.path.join(output_dir, "annotations_fetch_metadata.json"), "w",
          encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print(f"\nAnnotations: {len(all_results):,} | symbols: {len(symbols):,} "
      f"| taxa: {len(taxa):,}")
if truncated:
    print("NOTE: this is a TRUNCATED subset (see annotations_fetch_metadata.json)")