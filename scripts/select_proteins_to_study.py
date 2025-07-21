#!/usr/bin/env python3
"""
Select Proteins to Study

This script selects the top proteins from coverage data for detailed analysis.
It filters the coverage results to include only proteins from the pre-filtered
proteins_to_be_tested.txt list and selects the top N proteins based on coverage.

Input:
- coverage_count.tsv: Coverage analysis results with protein coverage statistics
- proteins_to_be_tested.txt: Pre-filtered list of surface accessible proteins

Output:
- proteins_to_study.tsv: Selected proteins for detailed analysis
"""

import pandas as pd
from pathlib import Path

# Get inputs from Snakemake
COVERAGE_FILE = Path(snakemake.input.coverage)
TESTED_PROTEINS_FILE = Path(snakemake.input.tested_proteins)
OUTPUT_FILE = Path(snakemake.output[0])

# Get parameters
analysis = snakemake.params.analysis
paramset = snakemake.params.paramset
group = snakemake.params.group
max_proteins = snakemake.params.max_proteins

print(f"Selecting Proteins to Study")
print("=" * 40)
print(f"Analysis: {analysis}")
print(f"Paramset: {paramset}")
print(f"Group: {group}")
print(f"Max proteins to select: {max_proteins}")
print(f"Coverage file: {COVERAGE_FILE}")
print(f"Tested proteins file: {TESTED_PROTEINS_FILE}")
print(f"Output file: {OUTPUT_FILE}")

# Create output directory
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# Load tested proteins list
tested_proteins = set()
if TESTED_PROTEINS_FILE.exists():
    with TESTED_PROTEINS_FILE.open("r") as f:
        tested_proteins = {line.strip() for line in f if line.strip()}
    print(f"Loaded {len(tested_proteins)} tested proteins")
else:
    print("Error: proteins_to_be_tested.txt not found")
    exit(1)

# Load coverage data
if not COVERAGE_FILE.exists():
    print("Error: Coverage file not found")
    exit(1)

try:
    coverage_df = pd.read_csv(COVERAGE_FILE, sep='\t')
    print(f"Loaded unified coverage data for {len(coverage_df)} gene-gram combinations")
    
    # Filter by Gram type for this specific analysis
    if 'gram' in coverage_df.columns:
        # Use group parameter from Snakemake (positive/negative)
        gram_type = group  # This comes from snakemake params
        filtered_df = coverage_df[coverage_df['gram'] == gram_type].copy()
        print(f"Filtered to {len(filtered_df)} entries for Gram-{gram_type}")
        coverage_df = filtered_df
    else:
        print("Warning: No 'gram' column found, using all data")
        
except Exception as e:
    print(f"Error loading coverage file: {e}")
    exit(1)

# Filter coverage data to include only tested proteins
if 'gene' in coverage_df.columns:
    protein_col = 'gene'
elif 'protein' in coverage_df.columns:
    protein_col = 'protein'
else:
    print("Error: Coverage file must have 'gene' or 'protein' column")
    exit(1)

# Filter to include only proteins in tested_proteins list
filtered_coverage = coverage_df[coverage_df[protein_col].isin(tested_proteins)].copy()

print(f"Coverage filtering results:")
print(f"  Total proteins in coverage: {len(coverage_df)}")
print(f"  Proteins in tested list: {len(tested_proteins)}")
print(f"  Matching proteins: {len(filtered_coverage)}")

if len(filtered_coverage) == 0:
    print("Warning: No proteins found in both coverage data and tested proteins list")
    # Create empty output file
    with OUTPUT_FILE.open('w') as f:
        f.write("protein\tcoverage_percentage\tspecies_with_protein\ttotal_species\n")
    exit(0)

# Filter by coverage percentage threshold (≥50%)
coverage_threshold = 50.0
if 'coverage_percentage' in filtered_coverage.columns:
    # First filter by coverage threshold
    high_coverage = filtered_coverage[filtered_coverage['coverage_percentage'] >= coverage_threshold].copy()
    print(f"Filtered to {len(high_coverage)} proteins with ≥{coverage_threshold}% coverage (from {len(filtered_coverage)} total)")
    
    if len(high_coverage) == 0:
        print(f"Warning: No proteins found with ≥{coverage_threshold}% coverage")
        # Create empty output file
        with OUTPUT_FILE.open('w') as f:
            f.write("protein\tcoverage_percentage\tspecies_with_protein\ttotal_species\n")
        exit(0)
    
    # Sort by coverage percentage (descending) and select top N
    sorted_coverage = high_coverage.sort_values('coverage_percentage', ascending=False)
else:
    print("Error: Coverage file must have 'coverage_percentage' column")
    exit(1)

# Select ALL proteins that meet the coverage threshold (no limit)
selected_proteins = sorted_coverage.copy()
print(f"Selected {len(selected_proteins)} proteins with ≥{coverage_threshold}% coverage")

# Rename the protein column to 'protein' for consistency
if protein_col == 'gene':
    selected_proteins = selected_proteins.rename(columns={'gene': 'protein'})

print(f"\nSelected {len(selected_proteins)} proteins for study:")
for idx, row in selected_proteins.iterrows():
    protein = row['protein']
    coverage = row['coverage_percentage']
    species_count = row.get('species_with_protein', row.get('species_with_gene', 0))
    print(f"  {protein}: {coverage:.1f}% coverage ({species_count} species)")

# Save selected proteins
try:
    selected_proteins.to_csv(OUTPUT_FILE, sep='\t', index=False)
    print(f"\n✅ Selected proteins saved to: {OUTPUT_FILE}")
    
    # Verify file was created
    if OUTPUT_FILE.exists():
        file_size = OUTPUT_FILE.stat().st_size
        print(f"✓ File created with size: {file_size} bytes")
    else:
        print("✗ Error: Output file was not created!")
        
except Exception as e:
    print(f"✗ Error saving results: {e}")

print(f"\n=== Protein Selection Complete ===")