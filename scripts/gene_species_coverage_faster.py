#!/usr/bin/env python3

from Bio import Entrez
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
output_file = "results/gene_species_coverage_gram_positive_faster.tsv"

os.makedirs(os.path.dirname(output_file), exist_ok=True)

# ---------------- LOAD INPUT ----------------

with open(species_file) as sf:
    species_list = [line.strip() for line in sf if line.strip()]

with open(gene_file) as gf:
    genes = [line.strip() for line in gf if line.strip()]

# ---------------- UTILS ----------------

def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

# ---------------- MAIN ----------------

results = []

for gene in genes:
    print(f"\n🔍 Checking gene: {gene}")
    
    # Build one big query combining gene and all species
    species_query = " OR ".join([f'"{species}"[Organism]' for species in species_list])
    full_query = f'"{gene}"[Gene] AND ({species_query})'
    #print(f"🔎 Query: {full_query}")

    try:
        # Step 1: Search for matching protein IDs
        search_handle = Entrez.esearch(db="protein", term=full_query, retmax=100000)
        search_result = Entrez.read(search_handle)
        search_handle.close()
        id_list = search_result["IdList"]

        if not id_list:
            print(f"❌ No hits found for gene {gene}")
            results.append((gene, 0, ""))
            continue

        # Step 2: Use esummary to fetch titles quickly
        matching_species = set()
        batch_size = 200

        for batch in chunk_list(id_list, batch_size):
            summary_handle = Entrez.esummary(db="protein", id=",".join(batch), retmode="xml")
            summaries = Entrez.read(summary_handle)
            summary_handle.close()

            for docsum in summaries:
                title = docsum.get("Title", "")
                for species in species_list:
                    if species.lower() in title.lower():
                        matching_species.add(species)

        print(f"✅ {gene}: found in {len(matching_species)} species")
        results.append((gene, len(matching_species), ",".join(sorted(matching_species))))

    except Exception as e:
        print(f"⚠️ Error while processing gene {gene}: {e}")
        results.append((gene, 0, ""))

    time.sleep(0.3)  # Safe for API key users

# ---------------- WRITE OUTPUT ----------------

with open(output_file, "w") as out:
    out.write("gene\tspecies_number\tspecies_with_gene\n")
    for gene, count, species_str in results:
        out.write(f"{gene}\t{count}\t{species_str}\n")

print(f"\n📄 Results saved to: {output_file}")
