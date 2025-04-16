#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import time

# ---------- Config ----------
with open(snakemake.input.ncbi_info, "r") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

gene_file = snakemake.input.gene_file
species = snakemake.params.species
output_file = snakemake.output.protein_names

# ---------- Load genes ----------
with open(gene_file) as f:
    gene_names = [line.strip() for line in f if line.strip()]

# ---------- Open output file ----------
with open(output_file, "w") as out:
    out.write("Gene\tProtein\n")
    out.flush()

    for gene in gene_names:
        query = f'"{gene}"[Gene] AND "{species}"[Organism]'

        try:
            print(f"🔎 Searching for {gene}...")
            handle = Entrez.esearch(db="protein", term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()

            ids = record.get("IdList", [])
            if not ids:
                out.write(f"{gene}\tNO MATCH\n")
                out.flush()
                continue

            # Fetch the first GenBank record
            fetch_handle = Entrez.efetch(db="protein", id=ids[0], rettype="gb", retmode="text")
            record = SeqIO.read(fetch_handle, "genbank")
            fetch_handle.close()

            protein_name = "NO PRODUCT FOUND"
            for feature in record.features:
                if feature.type == "Protein" and "product" in feature.qualifiers:
                    protein_name = feature.qualifiers["product"][0]
                    break

            out.write(f"{gene}\t{protein_name}\n")
            out.flush()
            time.sleep(0.4)

        except Exception as e:
            print(f"❌ Error with {gene}: {e}")
            out.write(f"{gene}\tERROR: {e}\n")
            out.flush()

print(f"\n✅ Done! Results saved to: {output_file}")
