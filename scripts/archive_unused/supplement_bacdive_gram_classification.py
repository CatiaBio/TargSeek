import sys
from pathlib import Path
import bacdive

# -------------------------
# USE FOR TESTING OUTSIDE SNAKEMAKE
# -------------------------
# Uncomment the following lines to run standalone:
# bacdive_gram_file = Path("data/bacdive/gram_classification.tsv")
# not_found_file = Path("data/bacdive/not_found.txt")
# fallback_file = Path("config/microbiome/fallback_gram_lookup.tsv")
# updated_gram_file = Path("data/bacdive/updated_gram_classification.tsv")
# updated_not_found_file = Path("data/bacdive/updated_not_found.txt")

# -------------------------
# USE INSIDE SNAKEMAKE
# -------------------------
bacdive_gram_file = Path(snakemake.input.bacdive_classification)
not_found_file = Path(snakemake.input.not_found)
updated_gram_file = Path(snakemake.output.updated_classification)
updated_not_found_file = Path(snakemake.output.updated_not_found)

# -----------------------------
# BacDive authentication
# -----------------------------
try:
    bacdive_info_path = snakemake.input.bacdive_info
except NameError:
    bacdive_info_path = "config/login/bacdive_info.txt"  # fallback for standalone

with open(bacdive_info_path, "r", encoding="utf-8") as f:
    lines = f.readlines()
    email = lines[0].strip()
    password = lines[1].strip() if len(lines) > 1 else None

client = bacdive.BacdiveClient(email, password)

# --- Load not-found species ---
with not_found_file.open("r", encoding="utf-8") as f:
    missing_species = [line.strip() for line in f if line.strip()]

# --- Load existing BacDive Gram classifications ---
existing_gram_data = {}
with bacdive_gram_file.open("r", encoding="utf-8") as f:
    next(f)  # skip header
    for line in f:
        species, gram = line.strip().split("\t")
        existing_gram_data[species] = gram.lower()

# --- Classify each not-found species via genus-based BacDive search ---
genus_cache = {}
updated_gram_data = existing_gram_data.copy()

for species in missing_species:
    genus = species.split()[0]
    try:
        if genus in genus_cache:
            records = genus_cache[genus]
        else:
            client.search(taxonomy=genus)
            records = list(client.retrieve())
            genus_cache[genus] = records

        strain_grams = set()

        def recurse(obj):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if k.lower() == "gram stain":
                        if isinstance(v, str):
                            strain_grams.add(v.lower())
                        elif isinstance(v, list):
                            strain_grams.update(x.lower() for x in v if isinstance(x, str))
                    else:
                        recurse(v)
            elif isinstance(obj, list):
                for item in obj:
                    recurse(item)

        for record in records:
            recurse(record)

        if "positive" in strain_grams:
            updated_gram_data[species] = "positive"
        elif "negative" in strain_grams:
            updated_gram_data[species] = "negative"
        else:
            updated_gram_data[species] = "unknown"

    except Exception:
        updated_gram_data[species] = "unknown"

# --- Save updated Gram classification with Unix line endings ---
with updated_gram_file.open("w", encoding="utf-8", newline="\n") as out:
    out.write("species\tgram_classification\n")
    for species in sorted(updated_gram_data):
        out.write(f"{species}\t{updated_gram_data[species]}\n")

print(f"✅ Updated Gram classification saved to: {updated_gram_file}")

# --- Update not_found file with still unknowns from original not_found list ---
still_unknown = [species for species in missing_species if updated_gram_data.get(species, "") == "unknown"]
with updated_not_found_file.open("w", encoding="utf-8") as nf:
    for species in sorted(still_unknown):
        nf.write(f"{species}\n")

print(f"📄 Updated list of unknowns saved to: {updated_not_found_file}")