#!/usr/bin/env python3

from Bio import Entrez
import os
import time

# ---------------- CONFIG ----------------

with open("config/login/ncbi_info.txt") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

species_file = "data/bacdive/gram_stain/gram_positive.txt"
gene_file = "data/quickgo/new_genes.txt"
output_file = "results/gene_species_coverage_gram_positive_bacteria.tsv"

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

with open(output_file, "w") as out:
    out.write("gene\tspecies_number\tspecies_with_gene\n")  # write header
    out.flush()

    for gene in genes:
        print(f"\n🔍 Checking gene: {gene}")

        species_query = " OR ".join([f'"{species}"[Organism]' for species in species_list])
        full_query = f'"{gene}"[Gene] AND ({species_query})'

        try:
            # Get total number of hits
            search_handle = Entrez.esearch(db="protein", term=full_query, retmax=0)
            search_result = Entrez.read(search_handle)
            search_handle.close()

            total_hits = int(search_result["Count"])
            print(f"🧬 Total hits: {total_hits}")

            if total_hits == 0:
                print(f"❌ No hits found for gene {gene}")
                out.write(f"{gene}\t0\t\n")
                out.flush()
                continue

            # Fetch all IDs with pagination
            id_list = []
            batch_size = 10000
            for start in range(0, total_hits, batch_size):
                print(f"📥 Fetching IDs {start + 1}–{min(start + batch_size, total_hits)}")
                search_handle = Entrez.esearch(
                    db="protein", term=full_query,
                    retstart=start, retmax=batch_size
                )
                search_result = Entrez.read(search_handle)
                search_handle.close()
                id_list.extend(search_result["IdList"])
                time.sleep(0.3)

            print(f"✅ Collected {len(id_list)} IDs")

            # Fetch titles and match species
            matching_species = set()
            for batch in chunk_list(id_list, 200):
                summary_handle = Entrez.esummary(db="protein", id=",".join(batch), retmode="xml")
                summaries = Entrez.read(summary_handle)
                summary_handle.close()

                for docsum in summaries:
                    title = docsum.get("Title", "")
                    for species in species_list:
                        if species.lower() in title.lower():
                            matching_species.add(species)

                time.sleep(0.3)

            print(f"✅ {gene}: found in {len(matching_species)} species")
            out.write(f"{gene}\t{len(matching_species)}\t{','.join(sorted(matching_species))}\n")
            out.flush()  # 💾 Force immediate write to disk

        except Exception as e:
            print(f"⚠️ Error while processing gene {gene}: {e}")
            out.write(f"{gene}\t0\t\n")
            out.flush()
            continue

print(f"\n📄 Results saved live to: {output_file}")
