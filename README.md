# 🎯 TargSeek: Surface Protein Discovery & Epitope Prediction Pipeline

A **Snakemake-based protein discovery pipeline** that identifies conserved, functionally relevant proteins within specific microbial groups. The pipeline processes taxonomic data and Gene Ontology (GO) terms to extract candidate proteins for downstream analysis such as diagnostic marker discovery, target validation, or phylogenomic studies.

## 🚀 Quick Start

```bash
# Clone and setup
git clone https://github.com/CatiaBio/TargSeek.git
cd TargSeek
conda env create -f env.yml
conda activate targseek

# Configure API credentials
cp config/login/bacdive_info.example.txt config/login/bacdive_info.txt
cp config/login/ncbi_info.example.txt config/login/ncbi_info.txt
# Edit with your credentials

# Run complete pipelines
snakemake -s Snakefile_download all_download_data --cores 8
snakemake -s Snakefile_analysis all_analysis --cores 8
```

### Running individual stages

```bash
# Download pipeline
snakemake -s Snakefile_download stage1_species_classification --cores 4  # BacDive Gram classification
snakemake -s Snakefile_download stage2_gene_discovery --cores 4          # QuickGO genes + taxa coverage
snakemake -s Snakefile_download stage3_protein_selection --cores 4       # Protein selection
snakemake -s Snakefile_download stage4_downloads --cores 4               # Sequences + 3D structures

# Analysis pipeline
snakemake -s Snakefile_analysis stage1_sequence_preparation --cores 4
snakemake -s Snakefile_analysis stage2_alignments_topology --cores 4     # MAFFT + DeepTMHMM topology
snakemake -s Snakefile_analysis stage3_conservation_analysis --cores 4   # ConSurf
snakemake -s Snakefile_analysis stage4_epitope_predictions --cores 4     # BepiPred 3.0 + reports

# Dry run (show what would be executed)
snakemake -s Snakefile_analysis all_analysis --dry-run
```

## 🎯 Overview

TargSeek identifies candidate proteins for vaccine development, diagnostics, and therapeutic targets through:

- Taxonomic classification and gene coverage assessment
- Transmembrane topology prediction with DeepTMHMM
- 3D structure integration from PDB
- Conservation analysis with ConSurf and BepiPred 3.0 epitope prediction

## 🏗️ Pipeline Architecture

### Two-Stage System:

1. **Download Pipeline** (`Snakefile_download`): Data collection, classification, sequence/structure downloads
2. **Analysis Pipeline** (`Snakefile_analysis`): MSA, conservation analysis, topology prediction, epitope prediction

### Key Features:

- **Coverage threshold**: ≥50% species coverage for gene selection
- **Intelligent caching**: Persistent API call caching to avoid re-downloads
- **3D structure priority**: Structures selected first for MSA
- **Dual analysis**: With and without 3D structures
- **Epitope conservation**: Maps epitopes to MSA positions and scores conservation

## 📁 Project Structure

```
TargSeek/
├── Snakefile_download           # Download pipeline workflow
├── Snakefile_analysis           # Analysis pipeline workflow
├── env.yml                      # Conda environment specification
├── config/                      # Configuration files
│   ├── config_download.yaml    # Download pipeline configuration
│   ├── config_analysis.yaml    # Analysis pipeline configuration
│   ├── login/                  # API credentials
│   ├── species/               # Species lists per analysis batch
│   └── quickgo/               # GO terms and parameters
├── scripts/                     # Pipeline scripts (organized by function)
│   ├── gene_selection/        # Download pipeline scripts
│   └── protein_analysis/      # Analysis pipeline scripts
├── cache/                      # Persistent API call caches
├── data/                       # Raw and intermediate data
├── results/                    # Pipeline analysis outputs
└── tools/                      # External tools (not included; see below)
```

## 📊 Output Organization

```
results/{analysis}_{paramset}/
├── gene_selection/
│   ├── coverage_count.tsv
│   ├── selected_genes_gram_*.txt
│   └── protein_structures/
└── protein_analysis/
    ├── stage1_sequence_preparation_3d_selection/      # MSA input sequences, 3D structure selection
    ├── stage2_multiple_sequence_alignments_topology/  # DeepTMHMM topology, MAFFT alignments
    ├── stage3_quality_assessment_conservation_analysis/  # ConSurf conservation
    └── stage4_epitope_predictions_and_reports/        # BepiPred epitopes, final reports
```

## 🔧 Configuration

- **Download**: `config/config_download.yaml`
- **Analysis**: `config/config_analysis.yaml`
- Modify thresholds, protein counts, and paths without editing Snakefiles

## 🧬 Dependencies

**Core** (via conda):

- Snakemake, Biopython, pandas, MAFFT, ClipKIT

**APIs**:

- BacDive, QuickGO, UniProt, NCBI, PDB

**Optional**:

- BepiPred 3.0 (Ubuntu), AliStat (manual compilation)

## 📚 Documentation

- [USAGE\_DOWNLOAD.txt](USAGE_DOWNLOAD.txt) - Download pipeline guide
- [USAGE\_ANALYSIS.txt](USAGE_ANALYSIS.txt) - Analysis pipeline guide

## 📦 External Tools & Licenses

This pipeline orchestrates several third-party tools which are **not included** in this repository.\
Users must install them separately and comply with their licenses:

- **MAFFT** (Multiple sequence alignment) – [License](https://mafft.cbrc.jp/alignment/software/source.html)
- **ConSurf** (Conservation analysis) – academic/non-commercial use only
- **DeepTMHMM** (Transmembrane topology prediction) – [DTU License](https://dtu.biolib.com/DeepTMHMM/)
- **BepiPred 3.0** (Epitope prediction) – academic/non-commercial use only
- **AliStat, ClipKIT, Biopython, pandas** – see their respective repositories

⚠️ These tools remain under their **original licenses**. This repository only provides wrappers/workflow rules to call them.

---

## 🗄️ Databases & Data Usage

This pipeline fetches data from publicly available biological databases. Each database has its own citation/usage requirements:

- **NCBI** – free use, cite [NCBI Resource Coordinators, 2018](https://www.ncbi.nlm.nih.gov/)
- **UniProt** – cite [UniProt Consortium, 2023](https://www.uniprot.org/)
- **PDB** – cite [Berman et al., 2000](https://www.rcsb.org/)
- **BacDive** – cite [Reimer et al., 2019](https://bacdive.dsmz.de/)
- **QuickGO (EBI/GOA)** – cite [Gene Ontology Consortium, 2021](http://geneontology.org/)

Please follow the recommended citations when publishing results.

---

## 📝 License & Attribution

- This repository’s **original code** is released under the [MIT License](LICENSE).
- Third-party tools and databases used in this pipeline remain under their **own licenses** and terms of use.
- If you use this pipeline in research, please cite both this repository and the external resources/tools.

**Citation for TargSeek:**

```bibtex
@software{baptista2025targseek,
  title={TargSeek: Surface Protein Discovery and Epitope Prediction Pipeline},
  author={Cátia Baptista},
  year={2025},
  url={https://github.com/CatiaBio/TargSeek}
}
```

