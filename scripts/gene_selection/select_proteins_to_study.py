#!/usr/bin/env python3
"""
Select Proteins to Study - Unified Version

This script selects proteins from coverage data for detailed analysis and creates:
1. A unified summary.tsv with both Gram-positive and Gram-negative data
2. Separate gene lists for each Gram type

Input:
- coverage_count.tsv: Unified coverage analysis results with Gram classification
- proteins_to_be_tested.txt: Pre-filtered list of surface accessible proteins
- protein_go_validation_report.tsv: Protein location information from GO validation
- uniprot/{paramset}/protein_info.json: UniProt protein information

Output:
- summary.tsv: Unified summary with both Gram types
- gram_positive.txt: List of selected Gram-positive genes
- gram_negative.txt: List of selected Gram-negative genes
"""

import pandas as pd
import json
from pathlib import Path
import yaml

# Get inputs from Snakemake
COVERAGE_FILE = Path(snakemake.input.coverage)
TESTED_PROTEINS_FILE = Path(snakemake.input.tested_proteins)
GO_VALIDATION_FILE = Path(snakemake.input.go_validation)
UNIPROT_INFO_FILE = Path(snakemake.input.uniprot_info)

# Get outputs from Snakemake
SUMMARY_TSV_FILE = Path(snakemake.output.summary_tsv)
GENE_LIST_POSITIVE_FILE = Path(snakemake.output.gene_list_positive)
GENE_LIST_NEGATIVE_FILE = Path(snakemake.output.gene_list_negative)

# Get parameters
analysis = snakemake.params.analysis
paramset = snakemake.params.paramset

print(f"Selecting Proteins to Study - Unified Version")
print("=" * 50)
print(f"Analysis: {analysis}")
print(f"Paramset: {paramset}")
print(f"Coverage file: {COVERAGE_FILE}")
print(f"Tested proteins file: {TESTED_PROTEINS_FILE}")
print(f"GO validation file: {GO_VALIDATION_FILE}")
print(f"UniProt info file: {UNIPROT_INFO_FILE}")
print(f"Summary TSV file: {SUMMARY_TSV_FILE}")
print(f"Positive gene list: {GENE_LIST_POSITIVE_FILE}")
print(f"Negative gene list: {GENE_LIST_NEGATIVE_FILE}")

# Create output directory
SUMMARY_TSV_FILE.parent.mkdir(parents=True, exist_ok=True)

# Load configuration for thresholds
try:
    with open('config/config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    
    coverage_thresholds = config.get('gram_thresholds', {'positive': 50, 'negative': 50})
    
    print(f"Coverage thresholds: {coverage_thresholds}")
except Exception as e:
    print(f"Warning: Could not load config: {e}")
    coverage_thresholds = {'positive': 50, 'negative': 50}

# Load UniProt protein information
uniprot_info = {}
if UNIPROT_INFO_FILE.exists():
    try:
        with UNIPROT_INFO_FILE.open('r') as f:
            uniprot_info = json.load(f)
        print(f"Loaded UniProt info for {len(uniprot_info)} proteins")
    except Exception as e:
        print(f"Warning: Could not load UniProt info: {e}")
else:
    print("Warning: UniProt info file not found, protein names will not be available")

# Load GO validation (location) data
location_data = {}
if GO_VALIDATION_FILE.exists():
    try:
        go_df = pd.read_csv(GO_VALIDATION_FILE, sep='\t')
        for _, row in go_df.iterrows():
            protein = row['protein']
            location = row.get('go_cellular_component', 'Unknown')
            # Use full location information
            if pd.notna(location) and location != 'NA':
                location = str(location).strip()
            else:
                location = 'Unknown'
            location_data[protein] = location
        print(f"Loaded location data for {len(location_data)} proteins")
    except Exception as e:
        print(f"Warning: Could not load location data: {e}")
else:
    print("Warning: GO validation file not found, locations will be 'Unknown'")

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
    print(f"Loaded coverage data for {len(coverage_df)} gene-gram combinations")
except Exception as e:
    print(f"Error loading coverage file: {e}")
    exit(1)

# Filter coverage data to include only tested proteins
filtered_coverage = coverage_df[coverage_df['gene'].isin(tested_proteins)].copy()

print(f"Coverage filtering results:")
print(f"  Total entries in coverage: {len(coverage_df)}")
print(f"  Proteins in tested list: {len(tested_proteins)}")
print(f"  Matching entries: {len(filtered_coverage)}")

if len(filtered_coverage) == 0:
    print("Warning: No proteins found in both coverage data and tested proteins list")
    # Create empty output files
    empty_df = pd.DataFrame(columns=['gram', 'gene_name', 'protein_name', 'location', 
                                   'coverage_percentage', 'species_with_gene', 'total_species', 
                                   'species_names_with_gene'])
    empty_df.to_csv(SUMMARY_TSV_FILE, sep='\t', index=False)
    GENE_LIST_POSITIVE_FILE.touch()
    GENE_LIST_NEGATIVE_FILE.touch()
    exit(0)

def get_protein_name(gene_symbol):
    """Get protein name from UniProt data"""
    if gene_symbol in uniprot_info:
        protein_data = uniprot_info[gene_symbol]
        if 'protein_name' in protein_data and protein_data['protein_name']:
            return protein_data['protein_name']
        elif 'gene_primary_name' in protein_data and protein_data['gene_primary_name']:
            return protein_data['gene_primary_name']
    return "Unknown protein"

def get_location(gene_symbol):
    """Get protein location from GO validation data"""
    return location_data.get(gene_symbol, 'Unknown')

# Process each Gram type separately
summary_data = []
selected_genes_positive = []
selected_genes_negative = []

for gram_type in ['positive', 'negative']:
    print(f"\n=== Processing Gram-{gram_type} proteins ===")
    
    # Filter by gram type
    gram_coverage = filtered_coverage[filtered_coverage['gram'] == gram_type].copy()
    print(f"Found {len(gram_coverage)} {gram_type} entries")
    
    if len(gram_coverage) == 0:
        print(f"No {gram_type} proteins found")
        continue
    
    # Apply coverage threshold
    threshold = coverage_thresholds[gram_type]
    high_coverage = gram_coverage[gram_coverage['coverage_percentage'] >= threshold].copy()
    print(f"After ≥{threshold}% coverage filter: {len(high_coverage)} proteins")
    
    if len(high_coverage) == 0:
        print(f"No {gram_type} proteins meet coverage threshold")
        continue
    
    # Sort by coverage (highest first)
    sorted_coverage = high_coverage.sort_values('coverage_percentage', ascending=False)
    
    # Select all proteins that meet the coverage threshold (no limits)
    selected_proteins = sorted_coverage.copy()
    print(f"Selected {len(selected_proteins)} {gram_type} proteins (all above threshold)")
    
    # Add to summary data
    for _, row in selected_proteins.iterrows():
        gene = row['gene']
        protein_name = get_protein_name(gene)
        location = get_location(gene)
        
        summary_data.append({
            'gram': gram_type,
            'gene_name': gene,
            'protein_name': protein_name,
            'location': location,
            'coverage_percentage': row['coverage_percentage'],
            'species_with_gene': row['species_with_gene'],
            'total_species': row['total_species'],
            'species_names_with_gene': row['species_names_with_gene']
        })
        
        # Add to gene lists
        if gram_type == 'positive':
            selected_genes_positive.append(gene)
        else:
            selected_genes_negative.append(gene)
    
    # Print top proteins for this gram type
    print(f"Top {gram_type} proteins:")
    for i, (_, row) in enumerate(selected_proteins.head(10).iterrows(), 1):
        gene = row['gene']
        protein_name = get_protein_name(gene)
        coverage = row['coverage_percentage']
        species_count = row['species_with_gene']
        print(f"  {i}. {gene} ({protein_name}): {coverage:.1f}% coverage ({species_count} species)")

# Create summary DataFrame and save
print(f"\n=== Creating unified summary ===")
summary_df = pd.DataFrame(summary_data)

# Sort by gram type and coverage
summary_df = summary_df.sort_values(['gram', 'coverage_percentage'], ascending=[True, False])

print(f"Final summary: {len(summary_df)} total proteins")
print(f"  Gram-positive: {len(selected_genes_positive)} proteins")
print(f"  Gram-negative: {len(selected_genes_negative)} proteins")

# Save summary TSV
summary_df.to_csv(SUMMARY_TSV_FILE, sep='\t', index=False)
print(f"✅ Summary saved to: {SUMMARY_TSV_FILE}")

# Save gene lists
with GENE_LIST_POSITIVE_FILE.open('w') as f:
    for gene in selected_genes_positive:
        f.write(f"{gene}\n")
print(f"✅ Positive gene list saved to: {GENE_LIST_POSITIVE_FILE}")

with GENE_LIST_NEGATIVE_FILE.open('w') as f:
    for gene in selected_genes_negative:
        f.write(f"{gene}\n")
print(f"✅ Negative gene list saved to: {GENE_LIST_NEGATIVE_FILE}")

print(f"\n=== Protein Selection Complete ===")
print(f"Summary file contains all selected proteins from both Gram types")
print(f"Selection criteria: Proteins meeting coverage thresholds (≥{coverage_thresholds['positive']}% positive, ≥{coverage_thresholds['negative']}% negative)")
print(f"Columns: gram, gene_name, protein_name, location, coverage_percentage, species_with_gene, total_species, species_names_with_gene")