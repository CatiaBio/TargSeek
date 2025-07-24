#!/usr/bin/env python3

"""
Description:
This script searches NCBI for a given protein name or UniProt ID across bacterial species,
downloads matching protein sequences, and saves each one individually in compressed FASTA format.
"""

# Libraries
from Bio import Entrez, SeqIO
import os
import gzip
import re
import time

# Read input and output from Snakemake
ncbi_info = snakemake.input.ncbi_info       # Input: email and API key
protein_query = snakemake.params.protein    # Input: protein name or UniProt ID
organism_query = snakemake.params.organism  # Input: organism name (e.g., "Staphylococcus aureus")
output_base = snakemake.output[0]           # Output file used as completion marker


if not snakemake.params.protein or not snakemake.params.organism:
    raise ValueError(f"Could not resolve protein/organism for query: {snakemake.wildcards.query}")

# Load Entrez credentials
with open(ncbi_info, "r") as info:
    lines = info.readlines()
    email = lines[0].strip()
    api_key = lines[1].strip() if len(lines) > 1 else None

Entrez.email = email
if api_key:
    Entrez.api_key = api_key

# Create a safe folder name based on query
safe_query = re.sub(r'[^\w\s-]', '', f"{protein_query}_{organism_query}").strip().replace(' ', '_')
output_dir = os.path.join("data", safe_query)
os.makedirs(output_dir, exist_ok=True)

def search_bacterial_proteins(protein_query, organism_query, database="protein"):
    full_query = f'"{protein_query}"[Protein Name] AND "{organism_query}"[Organism]'
    id_list = []

    handle = Entrez.esearch(db=database, term=full_query, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    count = int(record['Count'])
    print(f"Found {count} records for query: '{full_query}'")

    if count == 0:
        return []

    batch_size = 500
    for start in range(0, count, batch_size):
        handle = Entrez.esearch(db=database, term=full_query, retstart=start, retmax=batch_size)
        record = Entrez.read(handle)
        handle.close()
        id_list.extend(record["IdList"])

    return id_list

def fetch_and_save_proteins(id_list, directory):
    """
    Fetch protein sequences and save them individually.
    """
    batch_size = 200
    for start in range(0, len(id_list), batch_size):
        end = min(start + batch_size, len(id_list))
        batch_ids = id_list[start:end]
        
        fetch_handle = Entrez.efetch(
            db="protein", id=",".join(batch_ids),
            rettype="fasta", retmode="text"
        )
        records = SeqIO.parse(fetch_handle, "fasta")
        
        for record in records:
            accession = record.id.split("|")[0]
            description = record.description
            org_match = re.search(r'\[(.*?)\]$', description)
            org_name = org_match.group(1) if org_match else "Unknown"
            org_name = org_name.replace(" ", "_").replace("/", "_").replace("\\", "_")
            filename = f"{org_name}_{accession}.fasta.gz"
            filepath = os.path.join(directory, filename)

            with gzip.open(filepath, "wt") as f:
                SeqIO.write(record, f, "fasta")

        fetch_handle.close()
        time.sleep(0.4)  # small pause between batches to reduce NCBI load

def main():
    ids = search_bacterial_proteins(protein_query, organism_query)
    if not ids:
        print("No matching bacterial proteins found.")
        with open(output_base, "w") as f:
            f.write("No sequences found.\n")
        return

    fetch_and_save_proteins(ids, output_dir)
    print(f"Saved protein sequences to: {output_dir}")

    # Write the output flag file after successful run
    with open(output_base, "w") as f:
        f.write(f"Downloaded {len(ids)} records for '{protein_query}' in '{organism_query}'\n")


if __name__ == "__main__":
    main()
