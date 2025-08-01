# 📁 TargSeek Folder Organization

## 🏗️ Project Root Structure

```
TargSeek/
├── 📄 Snakefile_download          # Download pipeline workflow
├── 📄 Snakefile_analysis          # Analysis pipeline workflow
├── 📄 env.yml                     # Conda environment specification
├── 📄 README.md                   # Main project documentation
├── 📄 CLAUDE.md                   # Technical documentation
├── 📄 USAGE_DOWNLOAD.txt          # Download pipeline guide
├── 📄 USAGE_ANALYSIS.txt          # Analysis pipeline guide
├── 📁 config/                     # Configuration and input files
├── 📁 scripts/                    # Pipeline scripts
├── 📁 utils/                      # Utility scripts
├── 📁 tools/                      # External tools (BepiPred, DSSP)
├── 📁 cache/                      # Persistent API caches (gitignored)
├── 📁 data/                       # Downloaded sequences/structures (gitignored)
├── 📁 results/                    # Pipeline outputs (gitignored)
└── 📁 logs/                       # Execution logs (gitignored)
```

## ⚙️ Configuration Directory

```
config/
├── 📄 config_download.yaml        # Download pipeline configuration
├── 📄 config_analysis.yaml        # Analysis pipeline configuration
├── 📁 login/                      # API credentials (actual files gitignored)
│   ├── 📄 bacdive_info.example.txt
│   └── 📄 ncbi_info.example.txt
├── 📁 species/                    # Species lists by analysis
│   ├── 📄 analysis1.txt           # Primary analysis species
│   └── 📁 cow_milk/               # Additional species data
└── 📁 quickgo/                    # GO terms and parameters
    ├── 📄 go_id_descriptionl.tsv
    ├── 📄 params1.json
    └── 📄 surface_accessible.txt
```

## 🔧 Scripts Directory

### **Gene Selection Scripts** (`scripts/gene_selection/`)
```
├── 📄 classify_gram.py                      # BacDive Gram classification
├── 📄 fetch_quickgo_data.py                 # GO annotation retrieval
├── 📄 gene_taxa_coverage_unified.py         # NCBI coverage assessment
├── 📄 download_protein_sequences.py         # Multi-source sequence download
├── 📄 download_protein_structures.py        # PDB structure integration
├── 📄 extract_pdb_numbering.py              # PDB residue numbering
├── 📄 select_proteins_to_study.py           # Final protein selection
└── 📁 archived_unused/                      # Deprecated scripts
```

### **Protein Analysis Scripts** (`scripts/protein_analysis/`)
```
├── 📄 create_main_msa_sequences.py          # Main MSA preparation
├── 📄 select_3d_from_main_msa.py            # 3D structure selection
├── 📄 run_mafft_alignments.py               # MAFFT alignment execution
├── 📄 trim_alignments.py                    # ClipKIT trimming
├── 📄 analyze_conservation_adaptive.py       # Conservation analysis
├── 📄 predict_epitopes_bepipred_3d_only.py  # BepiPred 3.0 predictions
├── 📄 analyze_epitope_conservation.py       # Epitope conservation mapping
└── 📄 visualize_epitopes_pymol.py           # PyMOL visualization
```

## 🛠️ Utilities Directory

```
utils/
├── 📁 cache/                      # Cache management utilities
├── 📁 setup/                      # Installation scripts
│   ├── 📄 setup_bepipred.sh      # BepiPred 3.0 setup
│   └── 📄 test_bepipred.py       # BepiPred validation
└── 📄 convert_cif_to_pdb.py      # CIF to PDB conversion
```

## 📊 Data Flow Structure

### **Download Pipeline Output** (`results/{analysis}_{paramset}/gene_selection/`)
```
gene_selection/
├── 📄 coverage_count.tsv          # Filtered coverage data
├── 📄 selected_genes_gram_*.txt   # Gene lists by Gram type
├── 📄 summary.tsv                 # Selected proteins summary
├── 📁 genes_species/              # Gene-specific species lists
└── 📁 protein_structures/         # Downloaded 3D structures
    └── 📁 {gene}/
        ├── 📄 *.pdb               # PDB structure files
        └── 📄 structure_info.json # Metadata
```

### **Analysis Pipeline Output** (`results/{analysis}_{paramset}/protein_analysis/`)
```
protein_analysis/
├── 📁 main_msa/                   # Primary MSA directory
│   ├── 📄 {gene}_main.fasta       # All sequences for MSA
│   └── 📄 sequence_counts.json
├── 📁 msa_3d_structures/          # 3D-selected sequences
│   ├── 📄 {gene}_3d.fasta         # Sequences with structures
│   └── 📄 structure_mapping.tsv
├── 📁 no_3d/                      # Analysis without structures
│   ├── 📁 msa_alignments/         # MAFFT outputs
│   ├── 📁 msa_trimmed/            # ClipKIT outputs
│   └── 📁 conservation/           # Conservation analysis
├── 📁 with_3d/                    # Analysis with structures
│   └── [same structure as no_3d]
└── 📁 epitope_predictions_bepipred/
    ├── 📁 gram_positive/
    │   ├── 📄 {gene}_epitopes.tsv
    │   └── 📄 {gene}_conservation_report.html
    └── 📄 epitope_summary.json
```

## 🔄 Cache Organization

```
cache/
├── 📁 gene_species/               # NCBI coverage search cache
│   └── 📄 {gene}_{species}.json
├── 📁 protein_sequences/          # UniProt/NCBI sequence cache
│   └── 📄 {gene}_{species}.json
└── 📁 protein_structures/         # PDB structure search cache
    └── 📄 {uniprot_id}_structures.json
```

## 📝 Key Configuration Parameters

### **Download Pipeline** (`config_download.yaml`)
- Coverage thresholds: 50% for both Gram types
- Max proteins per group: 10
- API credentials: BacDive, NCBI
- Cache directories

### **Analysis Pipeline** (`config_analysis.yaml`)
- Max sequences per gene: 100
- MAFFT parameters
- ClipKIT settings
- Conservation thresholds
- BepiPred 3.0 paths

## 🎯 File Naming Conventions

```
Pattern: {analysis}_{paramset}_{component}_{gram_type}_{suffix}

Examples:
- analysis1_params1_coverage_count.tsv
- analysis1_params1_gram_positive_selected.txt
- analysis1_params1_conservation_summary.json
```

## 🔒 Git Tracking Policy

### ✅ **Tracked**
- Pipeline files (Snakefile_*)
- Scripts (scripts/)
- Configuration templates
- Documentation files
- Environment specification

### ❌ **Ignored**
- Data files (data/, results/)
- Cache directories
- API credentials
- Log files
- Compiled tools