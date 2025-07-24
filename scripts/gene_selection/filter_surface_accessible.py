#!/usr/bin/env python3
"""
Simple Surface Accessibility Filter

This script filters proteins based on GO cellular component annotations
to identify only surface accessible proteins.

Input:
- protein_go_validation_report.tsv: Protein GO validation results (protein\tgo_cellular_component)
- surface_accessible.txt: List of surface accessible GO cellular component terms

Output:
- surface_accessible_proteins.txt: List of proteins with surface accessible GO annotations
"""

import pandas as pd
from pathlib import Path

# Get inputs from Snakemake
GO_VALIDATION_FILE = Path(snakemake.input.go_validation)
SURFACE_ACCESSIBLE_FILE = Path(snakemake.input.surface_accessible)
OUTPUT_FILE = Path(snakemake.output.surface_genes)

print("Surface Accessibility Filter")
print("=" * 40)
print(f"GO validation file: {GO_VALIDATION_FILE}")
print(f"Surface accessible terms: {SURFACE_ACCESSIBLE_FILE}")
print(f"Output file: {OUTPUT_FILE}")

# Create output directory
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# Load surface accessible GO terms
surface_accessible_terms = set()
if SURFACE_ACCESSIBLE_FILE.exists():
    with SURFACE_ACCESSIBLE_FILE.open("r") as f:
        surface_accessible_terms = {line.strip().lower() for line in f if line.strip()}
    print(f"Loaded {len(surface_accessible_terms)} surface accessible GO terms")
else:
    print("Error: surface_accessible.txt not found")
    exit(1)

# Load GO validation report
if not GO_VALIDATION_FILE.exists():
    print("Error: GO validation report not found")
    exit(1)

try:
    go_df = pd.read_csv(GO_VALIDATION_FILE, sep='\t')
    print(f"Loaded GO validation for {len(go_df)} genes")
except Exception as e:
    print(f"Error loading GO validation report: {e}")
    exit(1)

def is_surface_accessible_go(go_terms, surface_terms):
    """Check if any GO cellular component term indicates surface accessibility"""
    if pd.isna(go_terms) or go_terms == 'NA':
        return False
    
    # Split multiple GO terms (separated by semicolon)
    terms_list = [term.strip().lower() for term in str(go_terms).split(';')]
    
    # Check if any term matches surface accessible terms
    for term in terms_list:
        if term in surface_terms:
            return True
        
        # Also check if any surface term is contained in the GO term
        for surface_term in surface_terms:
            if surface_term in term or term in surface_term:
                return True
    
    return False

# Filter proteins by surface accessibility
surface_proteins = []
total_proteins = 0
proteins_with_go = 0

print(f"\nFiltering proteins by surface accessibility...")

for _, row in go_df.iterrows():
    protein = row['protein']
    go_cc = row['go_cellular_component']
    total_proteins += 1
    
    if pd.notna(go_cc) and go_cc != 'NA':
        proteins_with_go += 1
        
        if is_surface_accessible_go(go_cc, surface_accessible_terms):
            surface_proteins.append(protein)
            print(f"  ✓ {protein}: {go_cc}")

print(f"\nFiltering Results:")
print(f"  Total proteins in GO validation: {total_proteins}")
print(f"  Proteins with GO CC annotations: {proteins_with_go}")
print(f"  Surface accessible proteins: {len(surface_proteins)}")
print(f"  Excluded (no GO CC): {total_proteins - proteins_with_go}")
print(f"  Excluded (not surface accessible): {proteins_with_go - len(surface_proteins)}")

# Sort proteins alphabetically
surface_proteins.sort()

# Save surface accessible proteins
with OUTPUT_FILE.open("w") as f:
    for protein in surface_proteins:
        f.write(protein + "\n")

print(f"\n✅ Surface accessible proteins saved to: {OUTPUT_FILE}")

if surface_proteins:
    print(f"\nFirst 10 surface accessible proteins:")
    for protein in surface_proteins[:10]:
        print(f"  {protein}")
    if len(surface_proteins) > 10:
        print(f"  ... and {len(surface_proteins) - 10} more")
else:
    print("⚠️  No surface accessible proteins found!")

print("=" * 40)