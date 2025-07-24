import requests
import os
import time
import csv


# Input files
goid_file = "config/quickgo/go_ids.tsv"
taxonid_file = "config/quickgo/taxon_ids.tsv"

# Output files
output_dir = "data/uniprot_annotations"
os.makedirs(output_dir, exist_ok=True)
annotations_output = os.path.join(output_dir, "annotations.tsv")
symbols_output = os.path.join(output_dir, "gene_symbols.txt")

# Read GO and taxon IDs (first column only)
def clean_lines(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        return [line.strip().split()[0].lower()  # Convert to lowercase
                for line in f
                if line.strip() and not line.startswith('#') and not line.lower().startswith(('gene', 'taxon'))]

go_ids = clean_lines(goid_file)  # Convert GO IDs to lowercase
taxon_ids = clean_lines(taxonid_file)
print(go_ids)
print(taxon_ids)    

# Prepare request
base_url = "https://rest.uniprot.org/uniprotkb/search?"
fields = "accession,gene_names,go_id,organism_id,organism_name"
headers = {"Accept": "application/tsv"}
limit = 500
sleep_sec = 1

results = []
symbols = set()

# Loop over GO IDs and taxon IDs
for go in go_ids:
    for taxon in taxon_ids:
        print(f"Querying: GO {go}, Taxon {taxon}")
        query = f"({go}) AND (taxonomy_id:{taxon})"
        print(query)
        cursor = "*"

        while cursor:
            params = {
                "query": query,
                "format": "tsv",
                "fields": fields,
                "size": limit,
                "cursor": cursor
            }
            r = requests.get(base_url, params=params, headers=headers)
            print(r)
            if not r.ok:
                print(f"Error: {r.status_code} {r.text}")
                break

            lines = r.text.strip().split('\n')
            if len(lines) <= 1:
                break  # no data

            if cursor == "*":
                header = lines[0].split('\t')
                results.append(header)

            for line in lines[1:]:
                row = line.split('\t')
                results.append(row)
                if row[1]:  # gene_names
                    symbols.update(row[1].split())

            # Handle cursor pagination
            cursor = r.headers.get("x-next-page-cursor", None)
            time.sleep(sleep_sec)

# Save annotations
with open(annotations_output, "w", encoding="utf-8", newline="") as out:
    writer = csv.writer(out, delimiter='\t')
    for row in results:
        writer.writerow(row)

# Save gene symbols
with open(symbols_output, "w", encoding="utf-8") as out:
    for symbol in sorted(symbols):
        out.write(symbol + "\n")

print(f"Saved {len(results)-1} entries to {annotations_output}")
print(f"Saved {len(symbols)} gene symbols to {symbols_output}")
