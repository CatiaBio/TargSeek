# Cow Milk Microbial Genera and Species Dataset

This dataset (genera_species_cow_milk.xlsx) includes microbial genera and species identified in cow milk, based on metagenomic analyses using MG-RAST and PathoScope tools. The original data was retrieved from the following scientific article:

> Hoque, M.N., Istiaq, A., Clement, R.A. et al. Metagenomic deep sequencing reveals association of microbiome signature with functional biases in bovine mastitis. Sci Rep 9, 13536 (2019). https://doi.org/10.1038/s41598-019-49468-4

---

## Dataset Contents

### Genera

genera.txt: Contains the combined genera identified by both MG-RAST and PathoScope.
genera_mg-rast.txt: Genera identified exclusively by MG-RAST.
genera_pathoscope.txt: Genera identified exclusively by PathoScope.
genera_taxids.txt: NCBI Taxonomy IDs associated with each genus.

### Species

species_CM.txt: Species identified specifically in cow milk samples from cows with clinical mastitis (CM). 
species_H.txt: Species identified specifically in cow milk samples from healthy (H) cows. 

---

## Methodology

The microbial community in cow milk was analyzed using two complementary bioinformatic tools:

- **MG-RAST**: Metagenomics Rapid Annotation using Subsystem Technology, which annotates metagenomic sequences against public databases.
- **PathoScope**: Software for identifying pathogens in metagenomic sequencing data.

Results from these analyses were combined to provide a comprehensive understanding of the microbial diversity within cow milk samples.

---

## Data Format

All files provided are in plain text (`.txt`) format:

- Genera and species lists: One entry per line.
- Genus taxids: Tab-delimited format with columns (`genus` and `taxid`).

---

## Data Retrieval Date
2025-04-04


