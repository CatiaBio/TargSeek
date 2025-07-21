# 🔬 QuickGO Configuration for TargSeek Pipeline

This folder contains **test configuration files** for the QuickGO API integration within the TargSeek protein discovery pipeline. The pipeline uses Gene Ontology (GO) terms to identify surface-accessible proteins suitable for vaccine development and diagnostic applications.

## 📁 Files in This Directory

### 🧪 **Test Configuration** (Tracked in Git)
- `test_params.json` - Sample parameter set for pipeline testing with 5 well-studied bacteria
- `README.txt` - This documentation file

### 🔧 **Production Configuration** (User-Provided, Gitignored)
Users must provide their own production configuration files:
- `go_ids.tsv` - Gene Ontology terms of interest
- `taxon_ids.tsv` - Taxonomic restrictions (NCBI taxonomy IDs)
- `surface_accessible.txt` - Surface accessibility criteria
- `params_1.json`, `params_2.json` - Parameter sets for different analyses

## 🎯 Pipeline Integration

The TargSeek pipeline uses QuickGO to:
1. **Retrieve gene symbols** based on specified GO terms and taxonomic constraints
2. **Filter proteins** for surface accessibility and cellular localization
3. **Assess gene coverage** across target bacterial species
4. **Select candidate proteins** for downstream epitope prediction

## 🧬 Target GO Terms for Surface Proteins

The pipeline focuses on proteins with cellular localizations relevant for **vaccine targets** and **diagnostic markers**:

### **Cell Surface & Extracellular**
- `GO:0005576` - Extracellular region
- `GO:0005615` - Extracellular space  
- `GO:0019861` - Outer membrane
- `GO:0030312` - External encapsulating structure
- `GO:0030313` - Cell envelope
- `GO:0016020` - Membrane (generic)

### **Cell Wall Components**
- `GO:0005618` - Cell wall
- `GO:0009274` - Peptidoglycan-based cell wall
- `GO:0009273` - Gram-positive cell wall
- `GO:0009275` - Gram-negative cell wall

### **Motility & Adhesion**
- `GO:0009278` - Bacterial-type flagellum
- `GO:0009288` - Fimbriae
- `GO:0009289` - Pilus
- `GO:0007155` - Cell adhesion

### **Transport & Virulence**
- `GO:0015288` - Porin activity
- `GO:0090729` - Toxin activity

## 📊 File Format Requirements

### **go_ids.tsv**
```tsv
go_id    go_name
GO:0005576    extracellular region
GO:0019861    flagellum
...
```

### **taxon_ids.tsv**
```tsv
taxon_id    taxon_name
2    Bacteria
1224    Proteobacteria
...
```

### **params_N.json**
```json
{
  "selectedIds": ["GO:0005576", "GO:0019861"],
  "taxonId": "2",
  "taxonUsage": "descendants",
  "geneProductType": "protein"
}
```

## 🚀 Getting Started

### **For Testing:**
1. Use the provided `test_params.json` with the test species list
2. Run: `snakemake fetch_quickgo_data --cores 4`

### **For Production:**
1. Create your GO term selection in `go_ids.tsv`
2. Define taxonomic scope in `taxon_ids.tsv`  
3. Configure parameters in `params_1.json`
4. Update `config/config.yaml` to reference your parameter files

## 🔗 Data Source & API

**QuickGO REST API**: https://www.ebi.ac.uk/QuickGO/api/  
**Interactive Interface**: https://www.ebi.ac.uk/QuickGO/annotations

### **API Rate Limits:**
- The pipeline includes built-in delays to respect EBI rate limits
- Large queries are automatically batched
- Failed requests are retried with exponential backoff

## 📖 Citation

If you use the TargSeek pipeline with QuickGO data, please cite:

> **TargSeek**: Baptista, C. (2024). TargSeek: Protein Discovery and Epitope Prediction Pipeline. GitHub. https://github.com/CatiaBio/PureMilk

> **QuickGO**: Binns D, et al. QuickGO: a web-based tool for Gene Ontology searching. Bioinformatics. 2009;25(22):3045-6.

## 🛠️ Troubleshooting

### **Common Issues:**
- **Empty gene lists**: Check GO term relevance for your target organisms
- **API timeouts**: Reduce batch sizes or add delays
- **Missing annotations**: Verify taxonomic IDs are correct

### **Parameter Optimization:**
- Start with broad GO terms, then refine based on results
- Use `taxonUsage: "descendants"` for taxonomic groups
- Filter by evidence codes if needed (e.g., experimental only)

---

**📁 This directory structure supports flexible, reproducible protein discovery workflows with comprehensive GO-based filtering.**