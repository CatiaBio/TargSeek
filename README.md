# 🧬 Protein Discovery Pipeline

This repository implements a pipeline that, given a list of **taxa** and **Gene Ontology (GO) terms**, identifies and extracts **candidate proteins** for downstream analysis.

The goal is to support the exploration of **conserved, functionally relevant proteins** within specific microbial groups — useful for tasks such as diagnostic marker discovery, target validation, or phylogenomic studies.

---

## 📋 Features

- ✅ Accepts custom lists of GO terms (from the [EBI QuickGO Annotations API](https://www.ebi.ac.uk/QuickGO/annotations)) and taxonomic identifiers
- ✅ Uses [BacDive](https://bacdive.dsmz.de/) to classify the provided taxa as **Gram-positive** or **Gram-negative**
- ✅ Extracts gene symbols associated with the input GO terms and filters them based on naming heuristics (e.g., excludes symbols starting with uppercase letters)
- ✅ Checks how many of the provided taxa contain proteins corresponding to each gene symbol using the **NCBI Protein** database



