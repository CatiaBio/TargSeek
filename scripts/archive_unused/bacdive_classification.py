#!/usr/bin/env python3

import os
import json
import time
import bacdive
import shutil

# ----------------------
# Snakemake I/O
# ----------------------
species_file = snakemake.input.species
credentials_file = snakemake.input.bacdive_info
output_json_path = snakemake.output[0]        # all_json
not_found_path = snakemake.output[1]          # not_found
output_errors_path = snakemake.output[2]      # errors
output_downloaded_path = snakemake.output[3]  # downloaded
output_gram_path = snakemake.output[4]        # gram_classification

# ----------------------
# BacDive Auth
# ----------------------
with open(credentials_file) as f:
    lines = f.readlines()
    email = lines[0].strip()
    password = lines[1].strip() if len(lines) > 1 else None

client = bacdive.BacdiveClient(email, password)

# ----------------------
# Load species list
# ----------------------
with open(species_file, "r") as f:
    species_list = [line.strip() for line in f if line.strip()]

# ----------------------
# Helper to flatten retrieve()
# ----------------------
def flatten_strain_entries(strains):
    for s in strains:
        if isinstance(s, dict):
            yield s
        elif isinstance(s, list):
            for item in flatten_strain_entries(s):
                yield item

# ----------------------
# Step 1: Save individual species JSONs
# ----------------------
not_found = []
errors = []
downloaded = []

cache_directory = "cache/bacdive/species_json"
os.makedirs(cache_directory, exist_ok=True)

for species in species_list:
    try:
        client.search(taxonomy=species)
        raw_data = list(client.retrieve())
        all_strains = list(flatten_strain_entries(raw_data))

        if not all_strains:
            not_found.append(species)
            continue

        species_filename = os.path.join(cache_directory, species.replace(" ", "_") + ".json")
        with open(species_filename, "w", encoding="utf-8") as f:
            json.dump(all_strains, f, indent=2, ensure_ascii=False)

        downloaded.append(species)
        time.sleep(1.5)

    except Exception as e:
        not_found.append(species)
        errors.append(f"{species}\t{str(e)}")

# ----------------------
# Step 2: Combine into one JSON file
# ----------------------
combined_data = {}
for file in os.listdir(cache_directory):
    if not file.endswith(".json"):
        continue
    species = file.replace("_", " ").replace(".json", "")
    with open(os.path.join(cache_directory, file), "r", encoding="utf-8") as f:
        combined_data[species] = json.load(f)

with open(output_json_path, "w", encoding="utf-8") as out_json:
    json.dump(combined_data, out_json, indent=2, ensure_ascii=False)

# ----------------------
# Save not found species
# ----------------------
with open(not_found_path, "w", encoding="utf-8") as nf:
    for name in sorted(not_found):
        nf.write(f"{name}\n")

# ----------------------
# Save species that raised errors
# ----------------------
with open(output_errors_path, "w", encoding="utf-8") as ef:
    for err in errors:
        ef.write(f"{err}\n")

# ----------------------
# Save successfully downloaded species
# ----------------------
with open(output_downloaded_path, "w", encoding="utf-8") as df:
    for species in sorted(downloaded):
        df.write(f"{species}\n")

# ----------------------
# Gram stain extraction helpers
# ----------------------
def extract_gram_from_morphology(entry):
    try:
        morph = entry.get("Morphology", {})
        cell_morph = morph.get("cell morphology", {})

        if isinstance(cell_morph, list):
            return [m.get("gram stain", "").lower() for m in cell_morph if isinstance(m, dict)]
        elif isinstance(cell_morph, dict):
            gram = cell_morph.get("gram stain", "").lower()
            return [gram] if gram else []
    except Exception:
        pass
    return []

def search_gram_stain_recursively(entry):
    grams = set()
    def recurse(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k.lower() == "gram stain":
                    if isinstance(v, str):
                        grams.add(v.lower())
                    elif isinstance(v, list):
                        grams.update(x.lower() for x in v if isinstance(x, str))
                else:
                    recurse(v)
        elif isinstance(obj, list):
            for item in obj:
                recurse(item)
    recurse(entry)
    return list(grams)

# ----------------------
# Classify Gram stain
# ----------------------
gram_classification = {}

for species, entries in combined_data.items():
    gram_result = "unknown"
    for entry in entries:
        grams = extract_gram_from_morphology(entry)
        if not grams:
            grams = search_gram_stain_recursively(entry)
        if "positive" in grams:
            gram_result = "positive"
            break
        elif "negative" in grams:
            gram_result = "negative"
            break
    gram_classification[species] = gram_result

# ----------------------
# Save Gram classification as TSV
# ----------------------
with open(output_gram_path, "w", encoding="utf-8") as gf:
    gf.write("species\tgram_classification\n")
    for species, gram in sorted(gram_classification.items()):
        gf.write(f"{species}\t{gram}\n")

# ----------------------
# Cleanup temporary files
# ----------------------
shutil.rmtree(cache_directory)
