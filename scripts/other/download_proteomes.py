#!/usr/bin/env python3

import os
import re
import time
import urllib.request
from Bio import Entrez

# ----------------------
# Settings
# ----------------------

with open("config/ncbi_info.txt", "r") as f:
    lines = f.readlines()
    Entrez.email = lines[0].strip()
    Entrez.api_key = lines[1].strip() if len(lines) > 1 else None

species_file = "data/bacdive/gram_stain/gram_positive.txt"
output_dir = "data/proteomes"
os.makedirs(output_dir, exist_ok=True)

# ----------------------
# Load species list
# ----------------------
with open(species_file, "r") as f:
    species_list = [line.strip() for line in f if line.strip()]

# ----------------------
# Download proteome for each species
# ----------------------
for species_name in species_list:
    safe_name = species_name.replace(" ", "_")
    output_path = os.path.join(output_dir, f"{safe_name}_proteome.faa.gz")

    if os.path.exists(output_path):
        print(f"🟡 Skipping {species_name} (already exists)")
        continue

    print(f"🔍 Searching for: {species_name}")
    try:
        # Step 1: Search NCBI Assembly database
        handle = Entrez.esearch(db="assembly", term=f"{species_name}[Organism]", retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            print(f"❌ No assembly found for {species_name}")
            continue

        assembly_id = record['IdList'][0]

        # Step 2: Get FTP path from esummary
        summary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary = Entrez.read(summary_handle)
        summary_handle.close()

        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        ftp_path = docsum.get('FtpPath_RefSeq') or docsum.get('FtpPath_GenBank')

        if not ftp_path:
            print(f"❌ No FTP path for {species_name}")
            continue

        basename = os.path.basename(ftp_path)
        faa_url = f"{ftp_path}/{basename}_protein.faa.gz"

        print(f"⬇️  Downloading: {faa_url}")
        urllib.request.urlretrieve(faa_url, output_path)
        print(f"✅ Saved: {output_path}")
        time.sleep(1.0)  # Be kind to NCBI

    except Exception as e:
        print(f"❌ Error for {species_name}: {e}")
