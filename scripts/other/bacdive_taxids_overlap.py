#!/usr/bin/env python3
import bacdive
import time
import os

# ----------------------
# Snakemake I/O
# ----------------------
genera_file = snakemake.input.genera
reference_taxids_file = snakemake.input.taxid_list
credentials_file = snakemake.input.bacdive_info 
output_dir = snakemake.params.output_dir
output_common = snakemake.output.common_taxids

# ----------------------
# BacDive Auth
# ----------------------
with open(credentials_file) as bacdive_info:
    lines = bacdive_info.readlines()
    email = lines[0].strip()
    password = lines[1].strip() if len(lines) > 1 else None

    
# Initialize the client
client = bacdive.BacdiveClient(email, password)

# ----------------------
# Read inputs
# ----------------------
with open(genera_file) as f:
    genera = [line.strip() for line in f if line.strip()]

with open(reference_taxids_file) as f:
    known_taxids = {line.strip() for line in f if line.strip().isdigit()}

# ----------------------
# Prepare output
# ----------------------
os.makedirs(output_dir, exist_ok=True)
found_taxids_global = set()

def extract_ncbi_ids():
    ncbi_ids = []
    for strain in client.retrieve():
        general_info = strain.get("General", {})
        ncbi_info = general_info.get("NCBI tax id", {})

        if isinstance(ncbi_info, dict):
            taxid = ncbi_info.get("NCBI tax id")
            if taxid:
                ncbi_ids.append(str(taxid))
        elif isinstance(ncbi_info, list):
            for entry in ncbi_info:
                taxid = entry.get("NCBI tax id")
                if taxid:
                    ncbi_ids.append(str(taxid))
    return ncbi_ids

for genus in genera:
    print(f"🔍 Searching genus: {genus}")
    try:
        client.search(taxonomy=genus)
        genus_taxids = []

        while True:
            batch_taxids = extract_ncbi_ids()
            if not batch_taxids:
                break

            genus_taxids.extend(batch_taxids)

            if len(batch_taxids) < 100:
                break
            time.sleep(1.5)

        genus_taxids = sorted(set(genus_taxids))
        found_taxids_global.update(genus_taxids)

        with open(os.path.join(output_dir, f"{genus}_taxid.txt"), "w") as f:
            f.write("\n".join(genus_taxids))

        print(f"✅ {len(genus_taxids)} taxids saved for {genus}")

    except Exception as e:
        print(f"❌ Error processing genus {genus}: {e}")

# Save common taxids
common_taxids = sorted(found_taxids_global & known_taxids)
with open(output_common, "w") as f:
    f.write("\n".join(common_taxids))

print(f"\n📄 Saved {len(common_taxids)} overlapping taxids to {output_common}")
