#!/usr/bin/env python3
"""
Ultra-fast gene-taxa coverage using combined queries and caching
"""

import pandas as pd
from Bio import Entrez
import time
from datetime import datetime
import os
import json
from urllib.error import HTTPError

# Read input files
with open(snakemake.input.species_list, 'r') as f:
    species_list = [line.strip() for line in f if line.strip()]

with open(snakemake.input.gene_list, 'r') as f:
    gene_list = [line.strip() for line in f if line.strip()]

# Set up Entrez
with open(snakemake.input.ncbi_info, 'r') as f:
    lines = f.readlines()
    email = lines[0].strip()
    api_key = lines[1].strip() if len(lines) > 1 else None

Entrez.email = email
if api_key:
    Entrez.api_key = api_key

print(f"Loaded {len(species_list)} species")
print(f"Loaded {len(gene_list)} genes")
print(f"NCBI API key: {'Configured' if api_key else 'Not found'}")

# Strategy: Combine multiple genes in a single query using OR
def build_combined_query(genes, species, max_genes_per_query=10):
    """Build combined query for multiple genes and one species"""
    gene_queries = []
    for gene in genes:
        gene_queries.append(f'("{gene}"[Gene Name] OR "{gene}"[All Fields])')
    
    combined_gene_query = " OR ".join(gene_queries)
    full_query = f'({combined_gene_query}) AND "{species}"[Organism]'
    return full_query

def search_with_retry(query, retmax=500, max_retries=3):
    """Search with retry logic"""
    delay = 0.1 if api_key else 0.34
    
    for attempt in range(max_retries):
        try:
            time.sleep(delay)
            handle = Entrez.esearch(db="protein", term=query, retmax=retmax, retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            return record
        except HTTPError as e:
            if e.code == 429:
                wait_time = 2 ** attempt
                print(f"Rate limit hit, waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                return None
        except Exception as e:
            print(f"Error: {e}")
            return None
    return None

def get_gene_from_protein(protein_id):
    """Extract gene information from protein record"""
    delay = 0.1 if api_key else 0.34
    try:
        time.sleep(delay)
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        for record in records:
            for feature in record.get('GBSeq_feature-table', []):
                if feature.get('GBFeature_key') == 'gene':
                    for qualifier in feature.get('GBFeature_quals', []):
                        if qualifier.get('GBQualifier_name') == 'gene':
                            return qualifier.get('GBQualifier_value')
        return None
    except:
        return None

print(f"\nStarting ultra-fast analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Initialize results
results = {gene: {species: 0 for species in species_list} for gene in gene_list}

# Process species in batches
start_time = time.time()
total_queries = len(species_list)

for species_idx, species in enumerate(species_list):
    # Search for all genes at once for this species
    all_genes_query = build_combined_query(gene_list, species, max_genes_per_query=50)
    
    # Get more results to capture all genes
    record = search_with_retry(all_genes_query, retmax=1000)
    
    if record and 'IdList' in record:
        # Now we need to figure out which protein belongs to which gene
        # For efficiency, we'll use the title/summary to match genes
        if len(record['IdList']) > 0:
            # Fetch summaries for all proteins at once
            delay = 0.1 if api_key else 0.34
            time.sleep(delay)
            
            try:
                handle = Entrez.esummary(db="protein", id=",".join(record['IdList']))
                summaries = Entrez.read(handle)
                handle.close()
                
                # Match summaries to genes
                for summary in summaries:
                    title = summary.get('Title', '').lower()
                    # Check which gene this protein belongs to
                    for gene in gene_list:
                        if gene.lower() in title:
                            results[gene][species] += 1
                            
            except Exception as e:
                print(f"Error fetching summaries for {species}: {e}")
    
    # Progress update
    if (species_idx + 1) % 10 == 0:
        elapsed = time.time() - start_time
        rate = (species_idx + 1) / elapsed
        eta = (len(species_list) - species_idx - 1) / rate if rate > 0 else 0
        print(f"Progress: {species_idx + 1}/{len(species_list)} species ({(species_idx + 1)/len(species_list)*100:.1f}%) - ETA: {eta/60:.1f} min")

# Convert results to DataFrame format
coverage_stats = []
for gene in gene_list:
    species_with_gene = sum(1 for species in species_list if results[gene][species] > 0)
    coverage_percentage = (species_with_gene / len(species_list)) * 100
    
    coverage_stats.append({
        'gene': gene,
        'species_with_gene': species_with_gene,
        'total_species': len(species_list),
        'coverage_percentage': coverage_percentage
    })

# Create coverage summary DataFrame
coverage_df = pd.DataFrame(coverage_stats)
coverage_df = coverage_df.sort_values('coverage_percentage', ascending=False)

# Save results
coverage_df.to_csv(snakemake.output[0], sep='\t', index=False)

# Print summary
print(f"\n{'='*60}")
print(f"Analysis completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total time: {(time.time() - start_time)/60:.1f} minutes")
print(f"\nResults saved to: {snakemake.output[0]}")
print(f"\nTop 10 genes by coverage:")
print(coverage_df.head(10).to_string(index=False))