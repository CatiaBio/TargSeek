# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PureMilk is a **Snakemake-based protein discovery pipeline** that identifies conserved, functionally relevant proteins within specific microbial groups. The pipeline processes taxonomic data and Gene Ontology (GO) terms to extract candidate proteins for downstream analysis such as diagnostic marker discovery, target validation, or phylogenomic studies.

## Architecture

The pipeline follows a **linear workflow** with these key stages:

1. **Taxonomic Classification**: Uses BacDive API to classify species as Gram-positive or Gram-negative
2. **GO Annotation Fetching**: Retrieves gene symbols from QuickGO API based on provided GO terms
3. **Gene Coverage Assessment**: Checks how many taxa contain proteins for each gene using NCBI Protein database
4. **Protein Selection**: Filters and selects top proteins based on coverage thresholds
5. **Sequence Download**: Downloads protein sequences from NCBI for selected genes
6. **Multiple Sequence Alignment**: Uses MAFFT for sequence alignment
7. **Alignment Processing**: Quality checks with AliStat and trimming with trimAl
8. **Conserved Region Analysis**: Identifies conserved amino acid positions

## Core Components

- **Snakefile**: Main workflow definition with all pipeline rules
- **config/config.yaml**: Central configuration file containing all paths, thresholds, and parameters
- **scripts/**: Python scripts for each pipeline step
- **config/**: Configuration files including login credentials and input data
- **data/**: Raw and processed data storage
- **results/**: Pipeline outputs organized by Gram classification

## Common Commands

### Environment Setup
```bash
# Create conda environment
conda env create -f env.yml
conda activate puremilk
```

### Pipeline Execution
```bash
# Run entire pipeline
snakemake --cores 8

# Run specific rule
snakemake classify_taxa_by_gram --cores 4

# Run with specific target
snakemake results/gram_positive_proteins.txt --cores 4

# Dry run to check workflow
snakemake --dry-run

# Generate workflow visualization
snakemake --dag | dot -Tsvg > workflow.svg
```

### Configuration
- All pipeline parameters are centralized in `config/config.yaml`
- Modify thresholds, protein selection counts, and file paths without touching the Snakefile
- Login credentials for BacDive and NCBI APIs are stored in `config/login/`

## Key Scripts

- `scripts/bacdive_classification.py`: Gram stain classification via BacDive API
- `scripts/fetch_quickgo_data.py`: GO annotation retrieval from QuickGO
- `scripts/gene_taxa_coverage.py`: NCBI protein database coverage assessment
- `scripts/download_proteins.py`: Protein sequence download from NCBI
- `scripts/get_msa_sequences.py`: Representative sequence selection for MSA

## Dependencies

Core tools managed via conda (`env.yml`):
- **Snakemake**: Workflow management
- **MAFFT**: Multiple sequence alignment  
- **trimAl**: Alignment trimming
- **Biopython**: Sequence parsing and API interactions
- **pandas**: Data manipulation
- **matplotlib, logomaker**: Visualization

**Note**: AliStat (alignment quality assessment) requires manual installation from GitHub and compilation.

## Data Structure

- Input species lists: `config/microbiome/cow_milk/unique_species.txt`
- GO terms: `config/quickgo/go_ids.tsv`
- Taxonomic IDs: `config/quickgo/taxon_ids.tsv`
- Results organized by Gram classification: `results/gram_positive/` and `results/gram_negative/`
- FASTA sequences stored in hierarchical structure: `results/{group}/{gene}/species.fasta`

## Configuration Patterns

The pipeline uses a template-based approach for handling Gram-positive vs Gram-negative processing:
- `gram_species_template: "data/bacdive/gram_{group}.txt"` allows dynamic file resolution
- Thresholds and selection criteria are group-specific via `gram_thresholds` and `protein_selection`