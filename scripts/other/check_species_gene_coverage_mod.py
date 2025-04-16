#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import os
import time
import urllib.error

# ---------------- CONFIG ----------------

with open("config/login/ncbi_info.txt") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

species_file = "data/bacdive/gram_stain/gram_positive.txt"
gene_file = "data/quickgo/gene_symbols.txt"
output_file = "results/gene_species_coverage_gram_positive.tsv"

os.makedirs(os.path.dirname(output_file), exist_ok=True)

# ---------------- LOAD INPUT ----------------

with open(species_file) as sf:
    species_list = [line.strip() for line in sf if line.strip()]

with open(gene_file) as gf:
    genes = [line.strip() for line in gf if line.strip()]

# ---------------- MAIN ----------------

results = []

for gene in genes:
    print(f"\n🔍 Checking gene: {gene}")
    
    # Build one query for all species
    species_query = " OR ".join([f'"{species}"[Organism]' for species in species_list])
    full_query = f'"{gene}"[Gene] AND ({species_query})'
    print(f"\n🔍 Checking gene: {gene}")
    print(f"🔎 NCBI query: {full_query}")


    try:
        search_handle = Entrez.esearch(db="protein", term=full_query, retmax=0)  # Only get count
        search_result = Entrez.read(search_handle)
        search_handle.close()

        hit_count = int(search_result.get("Count", 0))
        print(f"✅ {gene}: {hit_count} total hits")

        results.append((gene, hit_count))

    except Exception as e:
        print(f"⚠️ Error while processing gene {gene}: {e}")
        results.append((gene, 0))

    time.sleep(0.4)

# ---------------- WRITE OUTPUT ----------------

with open(output_file, "w") as out:
    out.write("gene\ttotal_hits\n")
    for gene, count in results:
        out.write(f"{gene}\t{count}\n")

print(f"\n📄 Results saved to: {output_file}")
