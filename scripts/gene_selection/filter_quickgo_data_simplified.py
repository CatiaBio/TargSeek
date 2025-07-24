#!/usr/bin/env python3
"""
Simplified QuickGO Data Filtering for Taxa Coverage

This script filters genes based on taxa coverage (≥5 taxa per gene).
It starts with a pre-filtered list of surface accessible genes.

Input:
- surface_accessible_proteins.txt: Pre-filtered list of surface accessible genes
- annotations.json: QuickGO annotation data in JSON format
- gene_aliases.json: Gene aliases for enhanced matching

Output:
- proteins_to_be_tested.txt: Final filtered protein list (surface accessible + ≥5 taxa)
"""

from pathlib import Path
import json
import os
import re

GENE_LIST_FILE = Path(snakemake.input.genes)  # Now surface_accessible_proteins.txt
ANNOTATION_FILE = Path(snakemake.input.annotations)
ALIASES_FILE = Path(snakemake.input.aliases)
OUTPUT_FILE = Path(snakemake.output.proteins_to_test)

# Create a taxon file output directory based on the paramset
paramset = snakemake.wildcards.paramset
output_dir = Path(f"data/quickgo/{paramset}/annotations_taxon_files")
output_dir.mkdir(parents=True, exist_ok=True)

print("Simplified QuickGO Taxa Coverage Filtering")
print("=" * 50)

# Load gene aliases
gene_aliases = {}
if ALIASES_FILE.exists():
    aliases_json_file = ALIASES_FILE.with_suffix('.json')
    if aliases_json_file.exists():
        with aliases_json_file.open("r") as f:
            gene_aliases = json.load(f)
        print(f"Loaded aliases for {len(gene_aliases)} genes")
    else:
        print("Warning: gene_aliases.json not found, proceeding without aliases")

# Load surface accessible genes (our starting list)
with GENE_LIST_FILE.open("r") as f:
    input_genes = [line.strip() for line in f if line.strip()]

print(f"Processing {len(input_genes)} surface accessible genes")

# Load annotations (JSON format)
with ANNOTATION_FILE.open("r", encoding="utf-8") as f:
    annotations = json.load(f)

print(f"Processing {len(annotations)} annotations")

# Create mapping for gene to taxa
gene_to_taxa = {gene: set() for gene in input_genes}

# Function to check if a gene or its aliases match a symbol
def gene_matches_symbol(gene, symbol, aliases_dict):
    """Check if gene or any of its aliases match the symbol"""
    if gene.lower() == symbol.lower():
        return True
    
    # Check aliases
    gene_aliases_list = aliases_dict.get(gene, [])
    for alias in gene_aliases_list:
        if alias.lower() == symbol.lower():
            return True
    
    return False

# Process annotations with alias support
print("Matching genes to taxa using annotations...")
matches_found = 0

for entry in annotations:
    symbol = entry.get("symbol", "")
    taxon_name = entry.get("taxonName", "")
    
    if symbol and taxon_name:
        for gene in input_genes:
            if gene_matches_symbol(gene, symbol, gene_aliases):
                gene_to_taxa[gene].add(taxon_name)
                matches_found += 1
                break

print(f"Found {matches_found} gene-taxon matches")

def safe_filename(name):
    return re.sub(r'[^a-zA-Z0-9_\\-]', '_', name)

# Filter genes by taxa coverage (≥5 taxa)
genes_with_5plus_taxa = []
genes_insufficient_taxa = []

print(f"\nFiltering by taxa coverage (≥5 taxa per gene):")

for gene, taxon_names in gene_to_taxa.items():
    if len(taxon_names) >= 5:
        genes_with_5plus_taxa.append(gene)
        
        # Create taxon file for this gene
        safe_name = safe_filename(gene)
        output_file = output_dir / f"{safe_name}_taxonName.txt"
        with output_file.open("w", encoding="utf-8") as f_out:
            for name in sorted(taxon_names):
                f_out.write(name + "\n")
    else:
        genes_insufficient_taxa.append(gene)
        print(f"  Excluded {gene}: only {len(taxon_names)} taxa")

print(f"\nTaxa Coverage Filtering Results:")
print(f"  Surface accessible genes (input): {len(input_genes)}")
print(f"  Genes with ≥5 taxa: {len(genes_with_5plus_taxa)}")
print(f"  Excluded (<5 taxa): {len(genes_insufficient_taxa)}")

# Sort final genes
final_genes = sorted(genes_with_5plus_taxa)

# Save final filtered genes
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT_FILE.open("w") as f:
    for gene in final_genes:
        f.write(gene + "\n")

print(f"\n+ Final protein list saved to: {OUTPUT_FILE}")
print(f"+ Taxon files saved to: {output_dir}")

if final_genes:
    print(f"\nFirst 10 proteins to be tested:")
    for gene in final_genes[:10]:
        taxa_count = len(gene_to_taxa[gene])
        print(f"  {gene} ({taxa_count} taxa)")
    if len(final_genes) > 10:
        print(f"  ... and {len(final_genes) - 10} more")

print(f"\n=== Taxa Coverage Filtering Complete ===")