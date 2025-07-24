import csv
import gzip
import os
import urllib.request

# --- Configuration ---
INPUT_FILE = "data/quickgo/gene_symbols.txt"
OUTPUT_FILE = "data/quickgo/gene_synonyms.tsv"
NCBI_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
LOCAL_GENE_INFO = "gene_info.gz"

# # --- Step 1: Download gene_info.gz if not already downloaded ---
# if not os.path.exists(LOCAL_GENE_INFO):
#     print("Downloading gene_info.gz from NCBI...")
#     urllib.request.urlretrieve(NCBI_URL, LOCAL_GENE_INFO)
#     print("Download complete.")

# --- Step 2: Load user-provided gene list ---
with open(INPUT_FILE, "r") as f:
    input_genes = [line.strip().lower() for line in f if line.strip()]

# --- Step 3: Parse gene_info.gz and build alias mapping ---
alias_map = {}  # key: alias -> value: (symbol, [synonyms])

with gzip.open(LOCAL_GENE_INFO, "rt") as f:
    reader = csv.reader(f, delimiter="\t")
    headers = next(reader)
    for row in reader:
        symbol = row[2]
        synonyms = row[4]
        if not symbol or symbol == "-":
            continue
        alias_list = [s for s in synonyms.split('|') if s and s != "-"]
        for name in [symbol] + alias_list:
            alias_map[name.lower()] = (symbol, alias_list)

# --- Step 4: Match input list to alias map ---
seen_symbols = set()
with open(OUTPUT_FILE, "w") as out:
    out.write("Gene\tAlso known as\n")
    for gene in input_genes:
        if gene in alias_map:
            symbol, aliases = alias_map[gene]
            if symbol not in seen_symbols:
                seen_symbols.add(symbol)
                aliases_clean = [a for a in aliases if a.lower() != symbol.lower()]
                out.write(f"{symbol}\t{', '.join(aliases_clean)}\n")
        else:
            print(f"Gene '{gene}' not found in NCBI gene_info")

print(f"\n✅ Finished! Output saved to: {OUTPUT_FILE}")
