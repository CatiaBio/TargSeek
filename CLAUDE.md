# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TargSeek is a **Snakemake-based protein discovery pipeline** that identifies conserved, functionally relevant proteins within specific microbial groups. The pipeline processes taxonomic data and Gene Ontology (GO) terms to extract candidate proteins for downstream analysis such as diagnostic marker discovery, target validation, or phylogenomic studies.

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
conda activate targseek
```

### Pipeline Execution

#### Download Pipeline (Snakefile_download)
```bash
# Run complete download pipeline
snakemake -s Snakefile_download all_download_data --cores 8

# Run individual download stages
snakemake -s Snakefile_download all_species_classification --cores 4
snakemake -s Snakefile_download all_gene_selection --cores 4
snakemake -s Snakefile_download all_downloads --cores 4

# Dry run to check download workflow
snakemake -s Snakefile_download --dry-run

# Generate download workflow visualization
snakemake -s Snakefile_download --dag | dot -Tsvg > download_workflow.svg

# Run with minimal output (suppress NumPy warnings)
PYTHONWARNINGS="ignore" snakemake -s Snakefile_download all_download_data --cores 4 --quiet

# Resume interrupted downloads (downloads use sentinel files for proper resumption)
snakemake -s Snakefile_download all_downloads --cores 4  # Will automatically detect and resume
```

#### Analysis Pipeline (Snakefile_analysis)
```bash
# Run complete analysis pipeline
snakemake -s Snakefile_analysis all_analysis --cores 8

# Run individual analysis stages
snakemake -s Snakefile_analysis all_sequence_preparation --cores 4
snakemake -s Snakefile_analysis all_alignments --cores 4
snakemake -s Snakefile_analysis all_quality_and_conservation --cores 4
snakemake -s Snakefile_analysis all_predictions_and_reports --cores 4

# Dry run to check analysis workflow
snakemake -s Snakefile_analysis --dry-run

# Generate analysis workflow visualization
snakemake -s Snakefile_analysis --dag | dot -Tsvg > analysis_workflow.svg

# Run BepiPred 3.0 epitope predictions (requires Ubuntu setup)
snakemake -s Snakefile_analysis all_predictions_and_reports --cores 4

# Run epitope predictions for specific group
snakemake -s Snakefile_analysis results/analysis1_params1/protein_analysis/epitope_predictions_bepipred/gram_positive --cores 4
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
- **Download pipeline**: All parameters centralized in `config/config_download.yaml`
- **Analysis pipeline**: All parameters centralized in `config/config_analysis.yaml`
- Modify thresholds, protein selection counts, and file paths without touching the Snakefiles
- Login credentials for BacDive and NCBI APIs are stored in `config/login/`
- Usage guides available: `USAGE_DOWNLOAD.txt` and `USAGE_ANALYSIS.txt`

## Key Scripts

### Download Pipeline Scripts (`scripts/gene_selection/`)
- `classify_gram.py`: Gram stain classification via BacDive API
- `fetch_quickgo_data.py`: GO annotation retrieval from QuickGO
- `gene_taxa_coverage_unified.py`: NCBI protein database coverage assessment
- `download_protein_sequences.py`: Protein sequence download with caching
- `download_3d_structures.py`: PDB structure download and integration
- `select_proteins_to_study.py`: Final protein selection based on thresholds

### Analysis Pipeline Scripts (`scripts/protein_analysis/`)
- `create_sequence_references.py`: MSA-ready sequence preparation
- `create_msa_fasta.py`: FASTA file creation for alignments
- `run_mafft_alignments.py`: Multiple sequence alignment with MAFFT
- `trim_alignments.py`: Alignment trimming with trimAl
- `assess_alignment_quality_comparison.py`: Quality assessment with AliStat
- `analyze_conservation_adaptive.py`: Conservation analysis and scoring
- `predict_epitopes_bepipred_3d_only.py`: BepiPred 3.0 epitope prediction

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

## Citation

When using TargSeek in your research, please cite:

> Baptista, C. (2024). TargSeek: Protein Discovery and Epitope Prediction Pipeline. GitHub. https://github.com/your-username/TargSeek

BibTeX format:
```bibtex
@software{baptista2024targseek,
  title={TargSeek: Protein Discovery and Epitope Prediction Pipeline},
  author={Cátia Baptista},
  year={2024},
  url={https://github.com/your-username/TargSeek}
}
```

## Project Organization

### Directory Structure

```
TargSeek/
├── Snakefile_download           # Download pipeline workflow
├── Snakefile_analysis           # Analysis pipeline workflow
├── env.yml                      # Conda environment specification
├── CLAUDE.md                    # Project guidance (this file)
├── README.md                    # Project overview
├── USAGE_DOWNLOAD.txt           # Download pipeline usage guide
├── USAGE_ANALYSIS.txt           # Analysis pipeline usage guide
├── config/                      # Configuration files
│   ├── config_download.yaml    # Download pipeline configuration
│   ├── config_analysis.yaml    # Analysis pipeline configuration
│   ├── login/                  # API credentials
│   ├── microbiome/            # Species lists by analysis
│   └── quickgo/               # GO terms and parameters
├── scripts/                     # Pipeline scripts (organized by function)
│   ├── gene_selection/        # Download pipeline scripts (14 scripts)
│   └── protein_analysis/      # Analysis pipeline scripts (8 scripts)
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
- Download workflow: `Snakefile_download`
- Analysis workflow: `Snakefile_analysis`
- Environment: `env.yml` 
- Download configuration: `config/config_download.yaml`
- Analysis configuration: `config/config_analysis.yaml`

**Script Organization:**
- Download pipeline: `scripts/gene_selection/` (14 scripts)
- Analysis pipeline: `scripts/protein_analysis/` (8 scripts)

**Utility Scripts:**
- Cache management: `utils/cache/`
- BepiPred setup: `utils/setup/`
- Data migration: `utils/migration/`

**Documentation:**
- Project guidance: `CLAUDE.md` 
- Download usage: `USAGE_DOWNLOAD.txt`
- Analysis usage: `USAGE_ANALYSIS.txt`
- Additional docs: `docs/`

## Data Structure

- Input species lists: `config/species/cow_milk/unique_species.txt`
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