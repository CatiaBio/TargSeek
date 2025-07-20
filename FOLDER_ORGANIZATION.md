# 📁 PureMilk Folder Organization

This document provides a comprehensive overview of the PureMilk project structure, file organization, and data flow patterns.

---

## 🏗️ Project Root Structure

```
PureMilk/
├── 📄 README.md                    # Main project documentation
├── 📄 CLAUDE.md                    # Technical documentation for Claude Code
├── 📄 FOLDER_ORGANIZATION.md       # This file - project structure guide
├── 📄 Snakefile                    # Main Snakemake workflow definition
├── 📄 env.yml                      # Conda environment specification
├── 📄 .gitignore                   # Git ignore patterns
├── 📊 full_pipeline_performance.json # Pipeline performance metrics
├── 📁 config/                      # Configuration files and parameters
├── 📁 scripts/                     # Python scripts for pipeline execution
├── 📁 data/                        # Raw and processed data (gitignored)
├── 📁 results/                     # Pipeline outputs (gitignored)
├── 📁 logs/                        # Execution logs (gitignored)
├── 📁 .snakemake/                  # Snakemake metadata (gitignored)
└── 📁 backup_before_clean/         # Backup files (gitignored)
```

---

## ⚙️ Configuration Directory (`config/`)

### **Main Configuration**
```
config/
├── 📄 config.yaml                  # Central pipeline configuration
└── 📁 login/                       # API credentials (gitignored)
    ├── 📄 bacdive_info.txt          # BacDive API credentials
    ├── 📄 ncbi_info.txt             # NCBI email and API key
    ├── 📄 bacdive_info.example.txt  # Template for BacDive config
    └── 📄 ncbi_info.example.txt     # Template for NCBI config
```

### **Input Data Configuration**
```
config/
├── 📁 microbiome/                  # Species lists for analysis
│   └── 📁 cow_milk/
│       └── 📄 unique_species.txt    # Target species list (297 species)
└── 📁 quickgo/                     # GO term and parameter configurations
    ├── 📄 go_ids.tsv                # Gene Ontology terms of interest
    ├── 📄 taxon_ids.tsv             # Taxonomic restrictions
    ├── 📄 surface_accessible.txt    # Surface accessibility criteria
    ├── 📄 params_1.json             # Primary parameter set
    └── 📄 params_2.json             # Alternative parameter set
```

**Key Configuration Parameters:**
- **Coverage Thresholds**: 50% for both gram-positive and gram-negative
- **Species Counts**: 113 gram-positive, 184 gram-negative
- **API Rate Limits**: Built-in delays and retry mechanisms
- **Batch Sizes**: 10 species per UniProt batch query

---

## 🔧 Scripts Directory (`scripts/`)

### **Core Pipeline Scripts**
```
scripts/
├── 📄 classify_gram.py                    # BacDive taxonomic classification
├── 📄 supplement_gram_classification.py   # Genus-based classification inference
├── 📄 fetch_quickgo_data.py              # QuickGO API integration
├── 📄 filter_quickgo_data.py             # Gene symbol filtering
├── 📄 gene_taxa_coverage_cached.py       # NCBI coverage assessment with caching
├── 📄 filter_coverage_count.py           # Coverage threshold filtering
├── 📄 create_proteins_to_study_from_coverage.py # Coverage-based protein selection
├── 📄 fetch_uniprot_info.py              # UniProt data enrichment
├── 📄 filter_surface_accessible_proteins.py # Surface accessibility filtering
├── 📄 create_gene_species_lists_from_coverage.py # Gene-specific species lists
└── 📄 download_proteins_combined.py      # Multi-stage protein downloading
```

### **3D Structure & MSA Scripts**
```
scripts/
├── 📄 download_3d_structures.py          # PDB structure downloading
├── 📄 get_msa_sequences.py               # MSA sequence selection with 3D priority
├── 📄 run_mafft_alignments.py           # MAFFT alignment execution
├── 📄 trim_alignments.py                # trimAl alignment optimization
├── 📄 assess_alignment_quality_comparison.py # Quality assessment
└── 📄 analyze_conservation_adaptive.py   # Conservation analysis
```

### **Advanced Analysis Scripts**
```
scripts/
├── 📄 predict_epitopes.py               # IEDB epitope prediction integration
├── 📄 generate_download_summary.py      # Download performance analysis
├── 📄 generate_final_report.py          # HTML report generation
└── 📄 monitor_pipeline.py               # Performance monitoring
```

### **Deprecated/Archive Scripts** (gitignored)
```
scripts/
├── 📄 download_proteins_simple.py       # Simplified download (superseded)
├── 📄 download_proteins_uniprot_iterative.py # UniProt-specific download
├── 📄 download_proteins_ncbi_missing.py # NCBI fallback download
└── 📄 analyze_conservation.py           # Basic conservation (superseded)
```

---

## 📊 Data Directory (`data/`) - Runtime Generated

### **BacDive Classification Data**
```
data/
└── 📁 bacdive/
    └── 📁 analysis_1/                   # Analysis-specific results
        ├── 📄 gram_raw.json             # Raw BacDive API responses
        ├── 📄 gram.tsv                  # Initial classification results
        ├── 📄 updated_gram.tsv          # Supplemented with genus inference
        ├── 📄 gram_positive.txt         # Gram-positive species list (113)
        ├── 📄 gram_negative.txt         # Gram-negative species list (184)
        ├── 📄 not_found.txt             # Species not found in BacDive
        ├── 📄 updated_not_found.txt     # Final unclassified species
        ├── 📄 errors.txt                # API errors encountered
        └── 📄 downloaded.txt            # Successfully processed species
```

### **QuickGO Annotation Data**
```
data/
└── 📁 quickgo/
    └── 📁 params_1/                     # Parameter set specific data
        ├── 📄 annotations.tsv           # Raw GO annotations
        ├── 📄 gene_symbols.txt          # Extracted gene symbols
        └── 📄 gene_symbols_filtered.txt # Filtered gene list (~186 genes)
```

---

## 📈 Results Directory (`results/`) - Pipeline Outputs

### **Coverage Analysis Results**
```
results/
└── 📁 coverage/
    ├── 📄 analysis_1_params_1_gram_positive_coverage.tsv      # Raw coverage data
    ├── 📄 analysis_1_params_1_gram_positive_coverage_count.tsv # Filtered (50% threshold)
    ├── 📄 analysis_1_params_1_gram_negative_coverage.tsv
    └── 📄 analysis_1_params_1_gram_negative_coverage_count.tsv
```

**Coverage File Structure:**
- `gene`: Gene symbol
- `species_with_gene`: Count of species containing the gene
- `coverage_percentage`: Percentage coverage across target species
- `species_names_with_gene`: Semicolon-separated species list

### **UniProt Enrichment Results**
```
results/
└── 📁 uniprot_info/
    └── 📁 analysis_1_params_1_gram_positive_uniprot_info/
        ├── 📄 analysis_1_params_1_gram_positive_coverage_count_location.tsv # Main output
        ├── 📄 download_summary.txt                           # Processing summary
        └── 📁 cache/                                          # API response cache
            ├── 📄 uniprot_cache_bacteria_gene1.json
            └── 📄 uniprot_cache_bacteria_gene2.json
```

### **Protein Selection Results**
```
results/
└── 📁 proteins_to_study/
    ├── 📄 analysis_1_params_1_gram_positive.tsv              # Selected proteins (G+)
    └── 📄 analysis_1_params_1_gram_negative.tsv              # Selected proteins (G-)
```

**Selected Proteins Structure:**
- `gene`: Gene symbol
- `go_cellular_component`: GO cellular component annotation
- `count`: Number of species with this gene
- `coverage_percentage`: Coverage percentage

### **Gene-Specific Species Lists**
```
results/
└── 📁 proteins_to_download/
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 gene1.txt                 # Species list for gene1
        ├── 📄 gene2.txt                 # Species list for gene2
        └── 📄 ...                       # One file per selected gene
```

### **Downloaded Protein Sequences**
```
results/
└── 📁 protein_fasta/
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 .download_complete        # Sentinel file for successful completion
        └── 📁 gene1/                    # Gene-specific sequences
            ├── 📄 species1.fasta        # Protein sequences by species
            ├── 📄 species2.fasta
            ├── 📄 download_summary.json # Download statistics
            └── 📄 not_found_species.txt # Failed downloads
```

### **3D Structure Integration**
```
results/
└── 📁 3d_structures/
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 analysis_1_params_1_gram_positive_summary.json # Structure overview
        └── 📁 gene1/                    # Gene-specific structures
            ├── 📄 1ABC.pdb              # PDB structure files
            ├── 📄 1ABC.fasta            # Corresponding sequences
            ├── 📄 structure_info.json   # Structure metadata
            └── 📄 no_structures_found.txt # Marker for no structures
```

### **MSA Preparation**
```
results/
└── 📁 msa_sequences/
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 gene1.fasta               # MSA-ready sequences (3D structure priority)
        ├── 📄 gene2.fasta
        ├── 📄 3d_structure_selection_report.json # 3D structure tracking
        └── 📄 3d_structure_summary.txt  # Human-readable structure report
```

### **Multiple Sequence Alignments**
```
results/
├── 📁 msa_alignments/               # Raw MAFFT alignments
│   └── 📁 analysis_1_params_1_gram_positive/
│       ├── 📄 gene1_aligned.fasta
│       └── 📄 gene2_aligned.fasta
├── 📁 msa_trimmed/                  # trimAl optimized alignments
│   └── 📁 analysis_1_params_1_gram_positive/
│       ├── 📄 gene1_trimmed.fasta
│       └── 📄 gene2_trimmed.fasta
└── 📁 msa_quality/                  # Quality assessment
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 quality_comparison.tsv    # Raw vs trimmed comparison
        ├── 📄 quality_summary.json     # Aggregate statistics
        └── 📄 gene1_quality_plot.png   # Quality visualization
```

### **Conservation Analysis**
```
results/
└── 📁 conservation/
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 conservation_summary.tsv  # Per-gene conservation metrics
        ├── 📄 conservation_summary.json # Detailed statistics
        └── 📁 per_gene/                 # Gene-specific analysis
            ├── 📄 gene1_conservation.tsv # Position-specific scores
            ├── 📄 gene1_conservation.png # Conservation plot
            └── 📁 logos/                # Logo plots (if enabled)
                └── 📄 gene1_logo_region1.png
```

### **Epitope Predictions**
```
results/
└── 📁 epitope_predictions/
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 epitope_prediction_summary.json # Overall summary
        ├── 📄 epitope_prediction_summary.tsv  # Tabular summary
        ├── 📄 epitope_prediction_report.txt   # Human-readable report
        └── 📁 gene1/                          # Gene-specific predictions
            ├── 📄 gene1_epitopes.tsv          # Detailed predictions
            └── 📄 gene1_epitope_summary.json  # Gene summary
```

**Epitope Prediction Structure:**
- **MHC Class I**: 9-10 amino acid peptides with HLA binding predictions
- **MHC Class II**: Variable length peptides with HLA-DR/DQ/DP predictions
- **B-cell**: Linear epitopes based on hydrophilicity and surface accessibility
- **Conservation Scoring**: Integration with position-specific conservation data

### **Download Performance Analysis**
```
results/
└── 📁 download_summary/
    └── 📁 analysis_1_params_1_gram_positive/
        ├── 📄 analysis_1_params_1_gram_positive_download_summary.tsv # Complete summary
        ├── 📄 analysis_1_params_1_gram_positive_selected_genes_download_summary.tsv # Selected only
        ├── 📄 analysis_1_params_1_gram_positive_download_stats.json # Statistics
        └── 📄 analysis_1_params_1_gram_positive_download_report.txt # Human-readable
```

### **Final Reports**
```
results/
└── 📁 reports/
    ├── 📄 analysis_1_params_1_final_report.html # Comprehensive HTML report
    └── 📁 assets/                               # Report assets (if any)
```

---

## 📝 Logs Directory (`logs/`) - Execution Logs

```
logs/
└── 📁 coverage/                     # Coverage analysis logs
    ├── 📄 analysis_1_params_1_gram_positive_coverage.log
    └── 📄 analysis_1_params_1_gram_negative_coverage.log
```

---

## 🔄 Data Flow Patterns

### **File Naming Convention**
```
{analysis}_{paramset}_gram_{group}_{suffix}

Examples:
- analysis_1_params_1_gram_positive_coverage.tsv
- analysis_1_params_1_gram_negative_epitopes.json
```

### **Dependencies Chain**
```
Species List → BacDive Classification → QuickGO Genes → Coverage Assessment 
    ↓
UniProt Enrichment → Surface Filtering → Gene Lists → Protein Download
    ↓  
3D Structures → MSA Sequences → Alignments → Quality Assessment
    ↓
Conservation Analysis → Epitope Prediction → Final Reports
```

### **Parallel Processing Points**
- **Coverage Assessment**: Gram-positive and negative processed in parallel
- **Protein Downloads**: Multiple genes downloaded simultaneously
- **MSA Generation**: Independent alignment jobs per gene
- **Quality Assessment**: Parallel evaluation of alignment methods

---

## 🛠️ Maintenance & Cleanup

### **Cache Management**
- **UniProt Cache**: `results/uniprot_info/*/cache/` - Can be cleaned to force re-download
- **Snakemake Cache**: `.snakemake/` - Snakemake metadata, safe to remove
- **Backup Files**: `backup_before_clean/` - Manual backups, review before deletion

### **Space Optimization**
- **Large Files**: Raw alignment files in `msa_alignments/` can be compressed
- **Intermediate Data**: Coverage cache files can be removed after successful runs
- **Logs**: Old log files can be archived or removed

### **Resumable Components**
- **Downloads**: Sentinel files enable automatic resumption
- **Coverage Analysis**: Cached results prevent re-computation
- **3D Structures**: Existing PDB files are not re-downloaded

---

## 🔒 Security & Credentials

### **Protected Files** (gitignored)
```
config/login/
├── 📄 bacdive_info.txt     # BacDive username/password
└── 📄 ncbi_info.txt        # NCBI email/API key
```

### **Template Files** (version controlled)
```
config/login/
├── 📄 bacdive_info.example.txt
└── 📄 ncbi_info.example.txt
```

**Note**: Never commit actual credentials to version control. Always use the example files as templates.

---

## 📊 Performance Considerations

### **Storage Requirements**
- **Minimal Run**: ~1-5 GB (selected genes only)
- **Full Analysis**: ~10-50 GB (complete dataset with structures)
- **Archive Storage**: ~100+ GB (with all intermediate files and plots)

### **Network Usage**
- **API Calls**: Thousands of requests to BacDive, QuickGO, UniProt, NCBI, PDB
- **Rate Limiting**: Built-in delays respect API limits
- **Caching**: Reduces redundant network requests

### **Computational Intensity**
- **CPU**: MSA generation and epitope prediction are CPU-intensive
- **Memory**: Large datasets may require 16+ GB RAM
- **I/O**: Frequent file operations during download and processing phases

---

## 🎯 Git Tracking Policy

### ✅ **Tracked Files** (Essential for reproducibility)
```
📄 Core Workflow
├── Snakefile                    # Main pipeline definition
├── env.yml                      # Environment specification
├── config/config.yaml           # Central configuration
└── config/login/*.example.txt   # Credential templates

📄 Active Scripts
├── scripts/classify_gram.py
├── scripts/fetch_quickgo_data.py
├── scripts/predict_epitopes.py
└── [all core pipeline scripts]

📄 Documentation
├── README.md
├── CLAUDE.md
└── FOLDER_ORGANIZATION.md
```

### ❌ **Ignored Files** (Generated/sensitive content)
```
📁 Data & Results
├── data/                        # All runtime data
├── results/                     # All pipeline outputs
├── logs/                        # Execution logs
└── backup_before_clean/         # Backup files

📁 Credentials & Cache
├── config/login/*.txt           # Actual API credentials
├── .snakemake/                  # Snakemake metadata
└── cache/                       # API response cache

📁 Deprecated & Archive
├── scripts/download_proteins_*.py  # Superseded scripts
├── scripts/analyze_conservation.py # Old versions
└── "Snakefile copy"             # Manual backups
```

---

## 🔧 Development Workflow

### **Adding New Features**
1. Develop new scripts in `scripts/`
2. Test thoroughly with small datasets
3. Update `Snakefile` with new rules
4. Update documentation and configuration
5. Commit with descriptive messages

### **Managing Deprecated Code**
1. Scripts are gitignored when superseded
2. Working versions remain in the repository
3. Major changes documented in commit messages
4. Archive folders for manual cleanup

### **Configuration Management**
1. All parameters centralized in `config/config.yaml`
2. API credentials use template system
3. Environment dependencies managed via `env.yml`
4. Multiple analysis configurations supported

---

**📁 This organization supports scalable, reproducible, and maintainable bioinformatics workflows with comprehensive tracking and quality control.**