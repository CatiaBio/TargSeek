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

# Run pipelines
snakemake -s Snakefile_download all_download_data --cores 8
snakemake -s Snakefile_analysis all_analysis --cores 8

# Run ConSurf conservation analysis
snakemake -s Snakefile_analysis all_consurf_analysis --cores 4

# Run topology prediction with DeepTMHMM
snakemake -s Snakefile_analysis all_topology_prediction --cores 4
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
├── CLAUDE.md                    # Project guidance
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
├── cache/                      # Persistent API call caches
├── data/                       # Raw and intermediate data
├── results/                    # Pipeline analysis outputs
└── tools/                      # External tools (BepiPred, ConSurf, etc.)
```

## 📊 Output Organization

```
results/{analysis}_{paramset}/
├── gene_selection/
│   ├── coverage_count.tsv
│   ├── selected_genes_gram_*.txt
│   └── protein_structures/
└── protein_analysis/
    ├── main_msa/              # Primary MSA directory
    ├── msa_3d_structures/     # 3D-selected sequences
    ├── conservation/
    ├── consurf_analysis/      # ConSurf conservation analysis
    ├── topology_predictions/  # DeepTMHMM topology predictions
    └── epitope_predictions_bepipred/
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

- [CLAUDE.md](CLAUDE.md) - Technical documentation
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
@software{baptista2024targseek,
  title={TargSeek: Surface Protein Discovery and Epitope Prediction Pipeline},
  author={Cátia Baptista},
  year={2025},
  url={https://github.com/CatiaBio/TargSeek}
}
```

