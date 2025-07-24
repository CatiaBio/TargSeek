#!/usr/bin/env python3

from Bio import Entrez
import os
import time
import shutil

# ---------------- CONFIG ----------------

credentials_file = snakemake.input.ncbi_info
species_file = snakemake.input.species
gene_file = snakemake.input.genes  # this is the missing genes file
existing_coverage = snakemake.input.existing_coverage
output_file = snakemake.output.updated_coverage

# Setup Entrez
with open(credentials_file) as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

os.makedirs(os.path.dirname(output_file), exist_ok=True)

# ---------------- LOAD INPUT ----------------

with open(species_file) as sf:
    species_list = [line.strip() for line in sf if line.strip()]

with open(gene_file) as gf:
    genes = [line.strip() for line in gf if line.strip()]

# Copy existing results to new output file
shutil.copy(existing_coverage, output_file)

# ---------------- UTILS ----------------

def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

# ---------------- MAIN ----------------

with open(output_file, "a") as out:
    for gene in genes:
        print(f"\n🔍 Checking gene: {gene}")

        species_query = " OR ".join([f'"{species}"[Organism]' for species in species_list])
        full_query = f'"{gene}"[Gene] AND ({species_query})'

        try:
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

            id_list = []
            batch_size = 10000
            for start in range(0, total_hits, batch_size):
                search_handle = Entrez.esearch(
                    db="protein", term=full_query,
                    retstart=start, retmax=batch_size
                )
                search_result = Entrez.read(search_handle)
                search_handle.close()
                id_list.extend(search_result["IdList"])
                time.sleep(0.3)

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
            out.flush()

        except Exception as e:
            print(f"⚠️ Error while processing gene {gene}: {e}")
            out.write(f"{gene}\t0\t\n")
            out.flush()
            continue

print(f"\n✅ Final output saved to: {output_file}")
