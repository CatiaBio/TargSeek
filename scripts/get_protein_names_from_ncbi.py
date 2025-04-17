#!/usr/bin/env python3

from Bio import Entrez
import time
import re

# ---------- Config ----------
with open(snakemake.input.ncbi_info, "r") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

gene_file = snakemake.input.gene_file
output_file = snakemake.output.protein_names

# ---------- Load genes ----------
with open(gene_file) as f:
    gene_names = [line.strip() for line in f if line.strip()]

# ---------- Open output file ----------
with open(output_file, "w") as out:
    out.write("Gene\tProtein\n")
    out.flush()

    for gene in gene_names:
        query = f'"{gene}"[Gene Name] OR {gene}[Gene] OR {gene}[All Fields]'

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

            # Use esummary to fetch the protein title
            summary_handle = Entrez.esummary(db="protein", id=ids[0])
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()

            if summary_record and "Title" in summary_record[0]:
                raw_title = summary_record[0]["Title"]
                protein_name = re.sub(r"\s*\[.*?\]$", "", raw_title).strip()
            else:
                protein_name = "NO TITLE FOUND"

            out.write(f"{gene}\t{protein_name}\n")
            out.flush()
            time.sleep(0.3)

        except Exception as e:
            print(f"❌ Error with {gene}: {e}")
            out.write(f"{gene}\tERROR: {e}\n")
            out.flush()

print(f"\n✅ Done! Results saved to: {output_file}")
