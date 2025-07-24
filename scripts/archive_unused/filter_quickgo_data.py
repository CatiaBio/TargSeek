import json
import os
import re
from pathlib import Path

# --- Option A: Snakemake paths ---
# GENE_LIST_FILE = Path(snakemake.input.genes)
# ANNOTATION_FILE = Path(snakemake.input.annotations)
# SUMMARY_FILE = Path(snakemake.output.filtered_genes)

# --- Option B: For standalone testing ---
GENE_LIST_FILE = Path("data/quickgo_test2/gene_symbols.txt")
ANNOTATION_FILE = Path("data/quickgo_test2/annotations.tsv")
SUMMARY_FILE = Path("data/quickgo_test2/genes_filtered.txt")

# --- Create output folder based on annotation file name ---
stem = ANNOTATION_FILE.stem  # e.g., "annotations_all"
output_dir = Path(f"data/quickgo/{stem}_taxon_files")
output_dir.mkdir(parents=True, exist_ok=True)

# --- Load gene list ---
with GENE_LIST_FILE.open("r") as f:
    input_genes = [line.strip() for line in f if line.strip()]
    gene_names_lower = set(g.lower() for g in input_genes)

# --- Load JSON annotations ---
with ANNOTATION_FILE.open("r", encoding="utf-8") as f:
    annotations = json.load(f)

# --- Map genes to taxon names ---
gene_to_taxa = {gene: set() for gene in input_genes}

for entry in annotations:
    symbol = entry.get("symbol", "")
    taxon_name = entry.get("taxonName")
    if symbol and symbol.lower() in gene_names_lower and taxon_name:
        for gene in input_genes:
            if gene.lower() == symbol.lower():
                gene_to_taxa[gene].add(taxon_name)
                break

# --- Helper to sanitize file names ---
def safe_filename(name):
    return re.sub(r'[^a-zA-Z0-9_\-]', '_', name)

# --- Save per-gene taxonName.txt files (only if ≥5 taxa) ---
genes_with_5plus = []

for gene, taxon_names in gene_to_taxa.items():
    if len(taxon_names) >= 5:
        genes_with_5plus.append(gene)
        safe_name = safe_filename(gene)
        output_file = output_dir / f"{safe_name}_taxonName.txt"
        with output_file.open("w", encoding="utf-8") as f_out:
            for name in sorted(taxon_names):
                f_out.write(name + "\n")

# --- Save summary list ---
SUMMARY_FILE.parent.mkdir(parents=True, exist_ok=True)
with SUMMARY_FILE.open("w") as f:
    for gene in sorted(genes_with_5plus):
        f.write(gene + "\n")

print("✅ Done!")
print(f"🗂 Gene files in:      {output_dir}")
print(f"📄 Summary file saved: {SUMMARY_FILE}")
