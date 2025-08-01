# 🎯 TargSeek: Surface Protein Discovery & Epitope Prediction Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A dual Snakemake pipeline for discovering conserved, functionally relevant proteins within microbial groups, with integrated 3D structure analysis and B-cell epitope prediction.

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
```

## 🎯 Overview

TargSeek identifies candidate proteins for vaccine development, diagnostics, and therapeutic targets through:
- Taxonomic classification and gene coverage assessment
- Surface accessibility filtering (membrane proteins)
- 3D structure integration from PDB
- Conservation analysis and BepiPred 3.0 epitope prediction

## 🏗️ Pipeline Architecture

### Two-Stage System:
1. **Download Pipeline** (`Snakefile_download`): Data collection, classification, sequence/structure downloads
2. **Analysis Pipeline** (`Snakefile_analysis`): MSA, conservation analysis, epitope prediction

### Key Features:
- **Coverage threshold**: ≥50% species coverage for gene selection
- **Intelligent caching**: Persistent API call caching to avoid re-downloads
- **3D structure priority**: Structures selected first for MSA
- **Dual analysis**: With and without 3D structures
- **Epitope conservation**: Maps epitopes to MSA positions and scores conservation

## 📁 Project Structure

```
TargSeek/
├── Snakefile_download      # Download workflow
├── Snakefile_analysis      # Analysis workflow
├── config/
│   ├── config_download.yaml
│   ├── config_analysis.yaml
│   └── species/           # Input species lists
├── scripts/
│   ├── gene_selection/    # Download pipeline scripts
│   └── protein_analysis/  # Analysis pipeline scripts
├── utils/                 # Utilities (cache, setup, conversion)
├── cache/                 # Persistent API caches
├── data/                  # Downloaded sequences/structures
└── results/               # Analysis outputs
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
- [USAGE_DOWNLOAD.txt](USAGE_DOWNLOAD.txt) - Download pipeline guide
- [USAGE_ANALYSIS.txt](USAGE_ANALYSIS.txt) - Analysis pipeline guide

## 📝 Citation

```bibtex
@software{baptista2024targseek,
  title={TargSeek: Protein Discovery and Epitope Prediction Pipeline},
  author={Cátia Baptista},
  year={2024},
  url={https://github.com/CatiaBio/TargSeek}
}
```

## 📝 License

MIT License - see [LICENSE](LICENSE) file.
