#!/usr/bin/env python3

"""
Description:
This script processes NCBI taxonomy files (nodes.dmp and names.dmp) to create two files:
1. taxonomy.tsv – taxid, scientific name, rank, and full lineage
2. lineage.tsv – extracted ranks for each species (genus, family, etc.)
"""

import os
import csv
import argparse

# -----------------------------
# Parse command-line arguments
# -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--taxid", type=str, default="2", help="Root taxid to filter lineage under (default: 1279 = Staphylococcus)")
args = parser.parse_args()
selected_tax_id = args.taxid

# -----------------------------
# File paths
# -----------------------------
nodes_file = "other/nodes.dmp"
names_file = "other/names.dmp"
taxonomy_output_file = "other/taxonomy.tsv"
lineage_output_file = "other/lineage.tsv"

# -----------------------------
# Load names.dmp (scientific names)
# -----------------------------
print("Parsing names.dmp...")
name_dict = {}
with open(names_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        tax_name = split_line[1].strip()
        name_class = split_line[3].strip()
        if name_class == 'scientific name':
            name_dict[tax_id] = tax_name

# -----------------------------
# Load nodes.dmp (hierarchy + ranks)
# -----------------------------
print("Parsing nodes.dmp...")
parent_dict = {}
rank_dict = {}
with open(nodes_file, 'r') as f:
    for line in f:
        split_line = line.split('|')
        tax_id = split_line[0].strip()
        parent_tax_id = split_line[1].strip()
        rank = split_line[2].strip()
        parent_dict[tax_id] = parent_tax_id
        rank_dict[tax_id] = rank

# -----------------------------
# Helper: build full lineage
# -----------------------------
def get_lineage(tax_id):
    lineage = [tax_id]
    while tax_id in parent_dict and parent_dict[tax_id] != '1':
        tax_id = parent_dict[tax_id]
        lineage.append(tax_id)
    return list(reversed(lineage))

# -----------------------------
# Helper: check if tax_id is under selected_tax_id
# -----------------------------
def is_in_lineage_of(tax_id, root_taxid):
    while tax_id in parent_dict and tax_id != '1':
        if tax_id == root_taxid:
            return True
        tax_id = parent_dict[tax_id]
    return False

# -----------------------------
# Write taxonomy.tsv
# -----------------------------
print("Writing taxonomy.tsv...")
with open(taxonomy_output_file, 'w') as f:
    f.write("taxid\tname\trank\tlineage\n")
    for tax_id in name_dict:
        if is_in_lineage_of(tax_id, selected_tax_id):
            lineage = get_lineage(tax_id)
            lineage_str = ','.join(lineage)
            rank = rank_dict.get(tax_id, 'no rank')
            f.write(f"{tax_id}\t{name_dict[tax_id]}\t{rank}\t{lineage_str}\n")

print(f"Taxonomy written to {taxonomy_output_file}")

# -----------------------------
# Build lineage.tsv (species-level info)
# -----------------------------
def process_taxonomy_file(taxonomy_file, output_file):
    lineage_dict = {}

    with open(taxonomy_file, 'r') as f:
        for line in f:
            if line.startswith("taxid"):
                continue  # skip header
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                taxon_id, name, rank, lineage = parts
                lineage_dict[taxon_id] = {
                    'rank': rank,
                    'name': name,
                    'lineage': lineage.split(',')
                }

    with open(output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(['taxid', 'species', 'genus', 'family', 'order', 'class', 'phylum'])

        for taxon_id, info in lineage_dict.items():
            if info['rank'].lower() == 'species':
                row = {
                    'taxid': taxon_id,
                    'species': info['name'],
                    'genus': '',
                    'family': '',
                    'order': '',
                    'class': '',
                    'phylum': ''
                }
                for ancestor_id in info['lineage']:
                    ancestor_info = lineage_dict.get(ancestor_id)
                    if ancestor_info:
                        ancestor_rank = ancestor_info.get('rank', '').lower()
                        if ancestor_rank in row:
                            row[ancestor_rank] = ancestor_info.get('name', '')
                writer.writerow([
                    row['taxid'], row['species'], row['genus'],
                    row['family'], row['order'], row['class'], row['phylum']
                ])
    print(f"Lineage written to {output_file}")

# Run it
process_taxonomy_file(taxonomy_output_file, lineage_output_file)