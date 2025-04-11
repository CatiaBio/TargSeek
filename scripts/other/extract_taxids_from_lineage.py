#!/usr/bin/env python3

import csv
import os

# ----------------------
# Snakemake I/O
# ----------------------
lineage_file = snakemake.input.lineage
species_list_file = snakemake.input.genera
output_exact = snakemake.output[0]  # matched_taxids.tsv
output_2ndround = snakemake.params["second_round_file"]
log_1st = snakemake.params["not_found_log"]
log_2nd = snakemake.params["not_found_log_2nd"]

# ----------------------
# Load species list
# ----------------------
with open(species_list_file) as f:
    input_species_raw = [line.strip() for line in f if line.strip()]
    input_species = {name.lower(): name for name in input_species_raw}

# ----------------------
# Build first-two-words map
# ----------------------
input_species_base = {
    " ".join(name.split()[:2]).lower(): name
    for name in input_species_raw
}

# ----------------------
# Load lineage entries
# ----------------------
lineage_map = {}  # lowercased species -> taxid
with open(lineage_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        species = row.get("species", "").strip()
        taxid = row.get("taxid", "").strip()
        if species and taxid:
            lineage_map[species.lower()] = taxid

# ----------------------
# First round: exact matches
# ----------------------
matches_1 = {}
not_found_1 = []

for lower_name, original_name in input_species.items():
    if lower_name in lineage_map:
        matches_1[original_name] = lineage_map[lower_name]
    else:
        not_found_1.append(original_name)

# ----------------------
# Second round: match on first 2 words
# ----------------------
matches_2 = {}
not_found_2 = []

for name in not_found_1:
    key = " ".join(name.lower().split()[:2])
    if key in lineage_map:
        matches_2[name] = lineage_map[key]
    else:
        not_found_2.append(name)

# ----------------------
# Write output: exact matches
# ----------------------
os.makedirs(os.path.dirname(output_exact), exist_ok=True)
with open(output_exact, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["species", "taxid"], delimiter="\t")
    writer.writeheader()
    for species, taxid in matches_1.items():
        writer.writerow({"species": species, "taxid": taxid})

# ----------------------
# Write output: 2nd round matches
# ----------------------
with open(output_2ndround, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["species", "taxid"], delimiter="\t")
    writer.writeheader()
    for species, taxid in matches_2.items():
        writer.writerow({"species": species, "taxid": taxid})

# ----------------------
# Logs
# ----------------------
with open(log_1st, "w") as f:
    for name in not_found_1:
        f.write(name + "\n")

with open(log_2nd, "w") as f:
    for name in not_found_2:
        f.write(name + "\n")

# ----------------------
# Summary
# ----------------------
print(f" Round 1: {len(matches_1)} exact matches")
print(f" Round 2: {len(matches_2)} genus+species matches")
print(f" Still not found: {len(not_found_2)}")
