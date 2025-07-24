import requests
import os
import time
import csv
import re
from itertools import islice

# ── Config ─────────────────────────────────────────────────────────────────────
GO_FILE         = "config/quickgo/go_ids.tsv"
TAXON_FILE      = "config/quickgo/taxon_ids.tsv"
OUTPUT_DIR      = "data/uniprot_annotations"
ANNOTATIONS_TSV = os.path.join(OUTPUT_DIR, "annotations.tsv")
SYMBOLS_TXT     = os.path.join(OUTPUT_DIR, "gene_symbols.txt")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Helpers ────────────────────────────────────────────────────────────────────
def clean_go_ids(path):
    """Extract 7-digit GO IDs (no prefix) from lines like ‘GO:0008150’."""
    with open(path, 'r', encoding='utf-8') as f:
        return [
            m.group(1)
            for line in f
            if (m := re.search(r'GO:(\d{7})', line))
        ]

def clean_taxon_ids(path):
    """Read one numeric taxon ID per line."""
    with open(path, 'r', encoding='utf-8') as f:
        return [
            line.strip().split()[0]
            for line in f
            if line.strip().split()[0].isdigit()
        ]

def chunked(iterable, n):
    """Yield successive n-sized chunks from iterable."""
    it = iter(iterable)
    return iter(lambda: list(islice(it, n)), [])

def fetch_tsv(query):
    """
    Run a single UniProt search (with cursor-based pagination)
    and return all rows including a header.
    """
    url     = "https://rest.uniprot.org/uniprotkb/search"
    headers = {"Accept": "application/tsv"}
    params  = {
        "query":  query,
        "format": "tsv",
        "fields": "accession,gene_names,go_id,organism_name",
        "size":   500,
        "cursor": "*"
    }

    rows, first = [], True
    while params.get("cursor"):
        r = requests.get(url, params=params, headers=headers)
        if not r.ok:
            raise RuntimeError(f"{r.status_code} {r.text}")
        lines = r.text.strip().split("\n")
        if len(lines) <= 1:
            break
        if first:
            rows.append(lines[0].split("\t"))
            first = False
        for line in lines[1:]:
            rows.append(line.split("\t"))
        params["cursor"] = r.headers.get("x-next-page-cursor")
        time.sleep(1)
    return rows

# ── Load inputs ────────────────────────────────────────────────────────────────
go_ids  = clean_go_ids(GO_FILE)
tax_ids = clean_taxon_ids(TAXON_FILE)
print("GO IDs:",    go_ids)
print("Taxon IDs:", tax_ids)

# Build a single taxonomy clause
tax_clause = "(" + " OR ".join(f"organism_id:{tid}" for tid in tax_ids) + ")"

# ── Query in chunks & collect results ──────────────────────────────────────────
all_rows, header = [], None
CHUNK_SIZE = 10

for batch in chunked(go_ids, CHUNK_SIZE):
    go_clause = "(" + " OR ".join(f"go:{gid}" for gid in batch) + ")"
    # Now require GO AND reviewed AND taxonomy
    query     = f"{go_clause} AND reviewed:true AND {tax_clause}"
    print("Running query:", query)
    block = fetch_tsv(query)
    if not block:
        continue
    if header is None:
        header = block[0]
    # append data rows (skip header)
    all_rows.extend(tuple(r) for r in block[1:])

# ── Deduplicate & write TSV ───────────────────────────────────────────────────
seen, uniq = set(), []
for row in all_rows:
    if row not in seen:
        seen.add(row)
        uniq.append(row)

with open(ANNOTATIONS_TSV, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    writer.writerows(uniq)

# ── Extract & save unique gene symbols ────────────────────────────────────────
symbols = {sym for row in uniq for sym in row[1].split()}
with open(SYMBOLS_TXT, "w", encoding="utf-8") as f:
    for s in sorted(symbols):
        f.write(s + "\n")

print(f"Saved {len(uniq)} annotation records.")
print(f"Saved {len(symbols)} unique gene symbols.")
