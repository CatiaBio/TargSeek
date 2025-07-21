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

# Run with minimal output (suppress NumPy warnings)
PYTHONWARNINGS="ignore" snakemake all_downloaded_proteins --cores 4 --quiet

# Resume interrupted downloads (downloads use sentinel files for proper resumption)
snakemake all_downloaded_proteins --cores 4  # Will automatically detect and resume

# Run BepiPred 3.0 epitope predictions (requires Ubuntu setup)
snakemake all_epitope_predictions_bepipred --cores 4

# Run epitope predictions for specific group
snakemake results/epitope_predictions_bepipred/analysis_1_params_1_gram_positive --cores 4
```

### BepiPred 3.0 Setup (Ubuntu/Linux)
```bash
# 1. Run the setup script (on Ubuntu)
./setup_bepipred.sh

# 2. Test the installation
python test_bepipred.py

# 3. Run epitope predictions
snakemake all_epitope_predictions_bepipred --cores 4
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

## Cache Management

The pipeline uses persistent caching to avoid re-downloading data:

- **`cache/gene_species/`**: Gene-taxa coverage search results (NCBI API calls)
- **`cache/protein_sequences/`**: Protein sequence downloads (UniProt API calls) 
- **`cache/protein_structures/`**: 3D structure searches (UniProt/PDB API calls)

### Protecting Cache Data

**⚠️ IMPORTANT: The cache contains valuable data that can take hours to rebuild!**

```bash
# Backup cache before major changes
python utils/cache/backup_cache.py backup

# Restore cache if accidentally deleted  
python utils/cache/backup_cache.py restore cache_backup_20250721_174500
```

The cache directories are protected in git:
- Directory structure is preserved (`.gitkeep` files)
- Cache contents are gitignored but not deleted
- Reinstalling/cloning preserves cache structure

### Cache Initialization

To populate cache with existing files:

```bash
# Initialize protein sequence cache
python utils/cache/initialize_protein_sequence_cache.py

# Initialize 3D structure cache  
python utils/cache/initialize_3d_structure_cache.py
```

## Dependencies

Core tools managed via conda (`env.yml`):
- **Snakemake**: Workflow management
- **MAFFT**: Multiple sequence alignment  
- **trimAl**: Alignment trimming
- **Biopython**: Sequence parsing and API interactions
- **pandas**: Data manipulation
- **matplotlib, logomaker**: Visualization

**Note**: AliStat (alignment quality assessment) requires manual installation from GitHub and compilation.

## Project Organization

### Directory Structure

```
PureMilk/
├── Snakefile                    # Main workflow definition
├── env.yml                      # Conda environment specification
├── CLAUDE.md                    # Project guidance (this file)
├── README.md                    # Project overview
├── config/                      # Configuration files
│   ├── config.yaml             # Main pipeline configuration
│   ├── login/                  # API credentials
│   ├── microbiome/            # Species lists by analysis
│   └── quickgo/               # GO terms and parameters
├── scripts/                     # Active pipeline scripts
│   ├── [18 core pipeline scripts]
│   └── archive_unused/        # Deprecated scripts (preserved)
├── scripts_test/               # Testing and development scripts
├── utils/                      # Utility scripts
│   ├── cache/                 # Cache management utilities
│   ├── setup/                 # Installation and setup scripts
│   └── migration/             # Data migration utilities
├── docs/                       # Documentation files
├── archive/                    # Archived temporary files
├── cache/                      # Persistent API call caches
├── data/                       # Raw and intermediate data
├── results/                    # Pipeline analysis outputs
└── tools/                      # External tools (BepiPred, etc.)
```

### Key File Locations

**Pipeline Files:**
- Main workflow: `Snakefile`
- Environment: `env.yml` 
- Configuration: `config/config.yaml`

**Utility Scripts:**
- Cache management: `utils/cache/`
- BepiPred setup: `utils/setup/`
- Data migration: `utils/migration/`

**Documentation:**
- Project guidance: `CLAUDE.md` 
- Additional docs: `docs/`

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

## Technical Notes

### Download Resumption
The `download_proteins_to_analyse` rule uses a **sentinel file approach** (`.download_complete`) to ensure proper handling of interrupted downloads:

- **Problem**: Snakemake's `directory()` output considers jobs complete if the directory exists, even if downloads were interrupted
- **Solution**: Added `.download_complete` sentinel file that's only created when ALL downloads succeed
- **Benefits**: 
  - Interrupted downloads are automatically detected and resumed
  - Existing files are checked to prevent re-downloading
  - No manual cleanup required after cancellation
  - Complete data integrity and smart resumption

**Usage**: Simply re-run the same Snakemake command after interruption - it will automatically resume from where it left off.