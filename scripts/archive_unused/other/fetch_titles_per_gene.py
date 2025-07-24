#!/usr/bin/env python3

from Bio import Entrez
import os
import time

# ---------------- CONFIG ----------------

with open("config/login/ncbi_info.txt") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

gene_file = "data/quickgo/gene_symbols.txt"
output_dir = "results/raw_titles"

os.makedirs(output_dir, exist_ok=True)

with open(gene_file) as gf:
    genes = [line.strip() for line in gf if line.strip()]

def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

# ---------------- MAIN LOOP ----------------

for gene in genes:
    print(f"\n🔍 Fetching titles for gene: {gene}")
    query = f'"{gene}"[Gene]'

    try:
        # Step 1: Get total number of hits
        search_handle = Entrez.esearch(db="protein", term=query, retmax=0)
        search_result = Entrez.read(search_handle)
        search_handle.close()

        total = int(search_result["Count"])
        print(f"🧬 Total hits for {gene}: {total}")

        if total == 0:
            print(f"❌ No results for {gene}")
            continue

        id_list = []
        batch_size = 10000

        # Step 2: Fetch all IDs using pagination
        for start in range(0, total, batch_size):
            print(f"📥 Fetching IDs {start + 1} to {min(start + batch_size, total)}...")
            search_handle = Entrez.esearch(
                db="protein",
                term=query,
                retstart=start,
                retmax=batch_size
            )
            search_result = Entrez.read(search_handle)
            search_handle.close()
            id_list.extend(search_result["IdList"])
            time.sleep(0.3)  # to stay within NCBI rate limits

        print(f"✅ Retrieved {len(id_list)} IDs for {gene}")

        # Step 3: Use esummary to fetch titles
        titles = []
        for batch in chunk_list(id_list, 200):  # esummary accepts up to 200 IDs
            summary_handle = Entrez.esummary(db="protein", id=",".join(batch), retmode="xml")
            summaries = Entrez.read(summary_handle)
            summary_handle.close()
            titles.extend(doc.get("Title", "") for doc in summaries)
            time.sleep(0.3)

        # Step 4: Save titles
        out_path = os.path.join(output_dir, f"{gene}_titles.txt")
        with open(out_path, "w") as out:
            for title in titles:
                out.write(title + "\n")

        print(f"📝 Saved {len(titles)} titles for {gene} → {out_path}")

    except Exception as e:
        print(f"⚠️ Error while processing gene {gene}: {e}")
