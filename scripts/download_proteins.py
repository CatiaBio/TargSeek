#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import os
import time
import urllib.error
import logging

# Set up logging
logging.basicConfig(filename='download_protein.log', level=logging.INFO,
                    format='%(asctime)s %(levelname)s:%(message)s')

# ---------------- CONFIG FROM SNAKEMAKE ----------------

with open(snakemake.input.ncbi_info, "r") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

protein_file = snakemake.input.proteins
species_file = snakemake.input.species
output_base = snakemake.output.output_folder

os.makedirs(output_base, exist_ok=True)

# ---------------- LOAD INPUTS ----------------

with open(protein_file) as pf:
    genes = [line.strip() for line in pf if line.strip()]

with open(species_file) as sf:
    species_list = [line.strip() for line in sf if line.strip()]

# ---------------- FUNCTIONS ----------------

def safe_filename(name):
    return name.replace(" ", "_").replace("/", "_")

def search_gene_protein_batch(gene, species_batch, retries=3):
    queries = [f'("{gene}"[Gene] OR "{gene}"[Protein]) AND "{species}"[Organism]' for species in species_batch]
    full_query = " OR ".join(f"({q})" for q in queries)

    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.esearch(db="protein", term=full_query, retmax=1000)
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"]
        except urllib.error.HTTPError as e:
            logging.warning(f"HTTP error for batch query {gene} attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(2 ** attempt)
    return []

def fetch_protein_sequences(id_list, retries=3):
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.efetch(db="protein", id=",".join(id_list), rettype="fasta", retmode="text")
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            return records
        except Exception as e:
            logging.warning(f"Fetch error attempt {attempt+1}: {e}")
            attempt += 1
            time.sleep(2 ** attempt)
    return []

# ---------------- MAIN LOOP ----------------

batch_size = 10
all_successful = True

for gene in genes:
    gene_dir = os.path.join(output_base, gene)
    os.makedirs(gene_dir, exist_ok=True)

    for i in range(0, len(species_list), batch_size):
        species_batch = species_list[i:i + batch_size]
        logging.info(f"Processing gene '{gene}' for species batch {species_batch}")

        ids = search_gene_protein_batch(gene, species_batch)
        if not ids:
            logging.warning(f"No IDs found for batch {gene}-{species_batch}")
            all_successful = False
            continue

        records = fetch_protein_sequences(ids)
        if not records:
            logging.warning(f"No sequences retrieved for batch {gene}-{species_batch}")
            all_successful = False
            continue

        # Organize sequences per species
        species_to_records = {species: [] for species in species_batch}
        for record in records:
            for species in species_batch:
                if species.lower() in record.description.lower():
                    species_to_records[species].append(record)

        for species, recs in species_to_records.items():
            fasta_out = os.path.join(gene_dir, f"{safe_filename(species)}.fasta")
            if recs:
                if not (os.path.exists(fasta_out) and os.path.getsize(fasta_out) > 0):
                    with open(fasta_out, "w") as f:
                        SeqIO.write(recs, f, "fasta")
                    logging.info(f"Successfully wrote: {fasta_out}")
            else:
                logging.warning(f"No records matched precisely for {gene}-{species}")

        time.sleep(0.2)  # To maintain safely below 10 req/sec

# # ---------------- WRITE COMPLETE FLAG ----------------

# if all_successful:
#     with open(snakemake.output.complete_flag, 'w') as f:
#         f.write('Download complete\n')
#     logging.info("Download complete flag written successfully.")
# else:
#     logging.warning("Some downloads failed; no complete flag written.")