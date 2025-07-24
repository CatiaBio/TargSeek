#!/usr/bin/env python3

"""
For each genus in a list, find its matching rows in the lineage file and group by phylum.
Outputs one file per (Phylum, Genus) pair: e.g., Proteobacteria_Pseudomonas.tsv
"""

import csv
import argparse
import os

# ----------------------
# Parse arguments
# ----------------------
parser = argparse.ArgumentParser()
parser.add_argument("lineage_file", help="Input lineage file (e.g., txid1279_lineage.txt)")
parser.add_argument("genera_file", help="Text file with one genus name per line")
parser.add_argument("-o", "--output-dir", default="other", help="Output directory")
args = parser.parse_args()

# ----------------------
# Load genus list
# ----------------------
with open(args.genera_file) as f:
    target_genera = {line.strip() for line in f if line.strip()}

# ----------------------
# Create output dir
# ----------------------
os.makedirs(args.output_dir, exist_ok=True)

# ----------------------
# Read and group rows
# ----------------------
grouped = {}  # (phylum, genus) -> list of rows

with open(args.lineage_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        genus = row['species']
        phylum = row['genera']
        if genus in target_genera:
            key = (phylum, genus)
            grouped.setdefault(key, []).append(row)

# ----------------------
# Write one file per (phylum, genus)
# ----------------------
for (phylum, genus), rows in grouped.items():
    phylum_safe = phylum.replace(" ", "_").replace("/", "_")
    genus_safe = genus.replace(" ", "_").replace("/", "_")
    out_file = os.path.join(args.output_dir, f"{phylum_safe}_{genus_safe}.tsv")
    with open(out_file, "w", newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=rows[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

print(f"Generated {len(grouped)} files in {args.output_dir}")