# BLAST Database Setup for Epitope Analysis

This document explains how to download and configure BLAST databases for epitope cross-reactivity analysis in the TargSeek pipeline.

## Quick Start

### Option 1: Use the convenience script (recommended)
```bash
# Essential databases only (human, viral, fungal, cow, mouse, ecoli)
./download_blast_databases.sh --quick

# All available databases
./download_blast_databases.sh --comprehensive

# Show available databases and current settings
./download_blast_databases.sh --list
```

### Option 2: Use Snakemake directly
```bash
# Download databases using current configuration
snakemake -s Snakefile_download all_blast_databases --cores 4

# Download as part of full pipeline
snakemake -s Snakefile_download all_download_data --cores 4
```

## Configuration

Edit `config/config_download.yaml` to select which databases to download:

```yaml
blast_databases:
  # Standard databases
  human: true
  viral_proteins: true
  fungal_proteins: true
  
  # Farm animals & livestock (important for milk project)
  cow: true
  goat: false
  horse: false
  pig: false
  chicken: false
  buffalo: false
  
  # Research model organisms
  mouse: true
  rat: false
  
  # Aquatic food organisms
  salmon: false
  zebrafish: false
  
  # Plant crops (food allergens/contamination)
  arabidopsis: false
  rice: false
  soybean: false
  
  # Bacterial pathogens
  ecoli: true
  salmonella: false
  staph: false
```

## Available Databases

### Standard Databases
- **human** - Human proteome (~25 MB) - Essential for human cross-reactivity
- **viral_proteins** - NCBI viral proteins (~15 MB) - Viral contamination check
- **fungal_proteins** - NCBI fungal proteins (~80 MB) - Fungal contamination check

### Farm Animals & Livestock
- **cow** - Cow proteome (~22 MB) - Critical for milk project
- **goat** - Goat proteome (~20 MB)
- **horse** - Horse proteome (~20 MB)
- **pig** - Pig proteome (~25 MB)
- **chicken** - Chicken proteome (~17 MB)
- **buffalo** - Water buffalo proteome (~20 MB)
- **sheep** - Sheep proteome (~20 MB)

### Companion/Domestic Animals
- **dog** - Dog proteome (~25 MB) - Common farm companion
- **cat** - Cat proteome (~20 MB) - Common farm companion

### Research Model Organisms
- **mouse** - Mouse proteome (~55 MB) - Common research model
- **rat** - Rat proteome (~30 MB) - Research model
- **rabbit** - Rabbit proteome (~22 MB) - Research model and farm animal

### Aquatic Food Organisms
- **salmon** - Atlantic salmon proteome (~35 MB)
- **zebrafish** - Zebrafish proteome (~45 MB) - Model organism

### Plant Crops
- **arabidopsis** - Arabidopsis proteome (~35 MB) - Plant model
- **rice** - Rice proteome (~40 MB)
- **soybean** - Soybean proteome (~56 MB)

### Bacterial Pathogens
- **ecoli** - E. coli K-12 proteome (~1.3 MB)
- **salmonella** - Salmonella enterica proteome (~1.5 MB)
- **staph** - Staphylococcus aureus proteome (~0.8 MB)

## Output Structure

Downloaded databases are stored in:
```
data/blast_databases/
├── .blast_databases_complete.sentinel    # Completion marker
├── download.log                          # Download log (JSON format)
├── human_proteome.fasta                  # FASTA files
├── human_proteome.fasta.phr              # BLAST database files
├── human_proteome.fasta.pin
├── human_proteome.fasta.psq
├── cow_proteome.fasta
├── cow_proteome.fasta.phr
└── ...
```

## Usage in Analysis

The downloaded databases can be used for epitope cross-reactivity analysis:

```bash
# BLAST epitope against human proteome
blastp -query epitope.fasta -db data/blast_databases/human_proteome.fasta -out results.txt

# BLAST with specific parameters
blastp -query epitope.fasta \
       -db data/blast_databases/cow_proteome.fasta \
       -evalue 1e-3 \
       -max_target_seqs 100 \
       -outfmt 6 \
       -out cow_hits.tsv
```

## Recommended Configurations

### For Milk Project (Quick Setup)
Essential databases for milk-related epitope analysis:
- human (human cross-reactivity)
- cow (primary milk source)
- viral_proteins (contamination)
- fungal_proteins (contamination)
- mouse (research model)
- ecoli (pathogen control)

### For Farm/Dairy Analysis
Extended farm animal coverage:
- human, cow, goat, sheep (core milk animals)
- dog, cat (farm companions)
- mouse, rabbit (farm-associated models)
- viral_proteins, fungal_proteins (contamination)
- ecoli, salmonella (pathogens)

### For Comprehensive Analysis
All databases for thorough cross-reactivity screening across:
- Multiple livestock species
- Companion/domestic animals
- Common research models
- Known pathogens
- Potential food allergens

## Dependencies

Required tools (install via conda):
```bash
conda install -c bioconda blast snakemake
```

## Troubleshooting

### Download Issues
- Check internet connectivity
- Verify sufficient disk space (up to ~500 MB for all databases)
- Check that UniProt/NCBI FTP servers are accessible

### BLAST Formatting Issues
- Ensure `makeblastdb` is in PATH
- Check write permissions in `data/blast_databases/`

### Configuration Issues
- Validate YAML syntax in `config/config_download.yaml`
- Ensure database names match exactly (case-sensitive)

## Integration with Pipeline

The BLAST databases are automatically included in the full download pipeline:

```bash
# Downloads databases first, then proceeds with species classification, etc.
snakemake -s Snakefile_download all_download_data --cores 4
```

The databases can be used in subsequent analysis steps for epitope validation and cross-reactivity testing.