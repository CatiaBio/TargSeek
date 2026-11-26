# TargSeek

A Snakemake pipeline that finds conserved, surface-exposed bacterial proteins and
predicts B-cell epitopes on them — for diagnostic markers, vaccine candidates and
target validation.

It works in two steps: a **download** pipeline that collects the data, and an
**analysis** pipeline that studies it.

## What it does

**Download pipeline** (`Snakefile_download`)

1. Classifies your species as Gram-positive or Gram-negative (BacDive)
2. Finds candidate genes for surface/extracellular GO terms (QuickGO)
3. Counts how many of your species carry each gene (NCBI), keeps those above the
   coverage threshold
4. Downloads protein sequences (NCBI) and 3D structures (PDB / AlphaFold)

**Analysis pipeline** (`Snakefile_analysis`)

1. Picks sequences for alignment and the best 3D structure per protein
2. Predicts membrane topology (DeepTMHMM) and keeps surface proteins
3. Aligns sequences (MAFFT) and checks alignment quality
4. Scores conservation (ConSurf)
5. Predicts epitopes (BepiPred 3.0), filters them to extracellular regions, and
   checks cross-reactivity against reference proteomes (BLAST)

## Install

```bash
conda env create -f env.yml
conda activate targseek
```

Some tools install separately: **ConSurf**, **BepiPred 3.0** and **DeepTMHMM** are
required; **DeepLocPro** (localization) and **DiscoTope-3.0** (structure-based
epitopes) are optional. See [docs/external_tools_setup.md](docs/external_tools_setup.md)
and set their paths in `config/config_analysis.yaml`.

## Set up your input

These files are not in the repository — create them before running:

| File | What it is |
|---|---|
| `config/login/bacdive_info.txt` | BacDive email (line 1), password (line 2) |
| `config/login/ncbi_info.txt` | NCBI email, API key recommended |
| `config/species/<batch>.txt` | your species, one per line |
| `config/quickgo/<paramset>.json` | the GO terms to search |

Example credential files are provided as `*.example.txt`.

Then choose which dataset to run, in **both** config files:

```yaml
species_batches:
  - analysis1          # which species list
quickgo_paramsets:
  - params1            # which GO term set
```

## Run

```bash
# always check first - shows what would run, changes nothing
snakemake -s Snakefile_download all_download_data -n

# download, then analyse
snakemake -s Snakefile_download all_download_data --cores 4
snakemake -s Snakefile_analysis all_analysis --cores 4
```

Interrupted? Run the same command again — finished work is skipped and downloads
resume. Add `--keep-going` so one failing step doesn't stop the others.

Run a single stage if you prefer:

```bash
snakemake -s Snakefile_download stage1_species_classification --cores 4
snakemake -s Snakefile_download stage2_gene_discovery --cores 4
snakemake -s Snakefile_download stage3_protein_selection --cores 4
snakemake -s Snakefile_download stage4_downloads --cores 4

snakemake -s Snakefile_analysis stage1_sequence_preparation --cores 4
snakemake -s Snakefile_analysis stage2_alignments_topology --cores 4
snakemake -s Snakefile_analysis stage3_conservation_analysis --cores 4
snakemake -s Snakefile_analysis stage4_epitope_predictions --cores 4
```

## Results

Everything lands under `results/{batch}_{paramset}/`:

```
gene_selection/              coverage table, selected genes, HTML report
protein_analysis/
  stage1_..._3d_selection/   sequences for alignment, chosen structures
  stage2_..._topology/       topology, surface proteins, alignments, QC report
  stage3_..._conservation/   ConSurf scores
  stage4_..._reports/        epitopes, cross-reactivity, final reports
```

Downloaded data is shared between runs: `data/protein_sequences/`,
`data/protein_structures/`, `data/blast_databases/`.

`cache/` holds API results and can take hours to rebuild — don't delete it.

## Settings

Everything is in `config/config_download.yaml` and `config/config_analysis.yaml`;
you should not need to edit the Snakefiles. The ones you are most likely to change:

| Setting | Meaning |
|---|---|
| `gram_thresholds` | minimum % of species that must carry a gene (default 50) |
| `topology.min_extracellular_aa` | residues needed to count as surface-exposed |
| `quality_control.min_species_per_alignment` | smallest alignment used for conservation |
| `bepipred.threshold` | epitope score cut-off |
| `localization.mode` | how topology and localization evidence combine |

## Documentation

- [USAGE_DOWNLOAD.txt](USAGE_DOWNLOAD.txt) — download pipeline guide
- [USAGE_ANALYSIS.txt](USAGE_ANALYSIS.txt) — analysis pipeline guide
- [BLAST_DATABASES.md](BLAST_DATABASES.md) — reference proteomes
- [docs/external_tools_setup.md](docs/external_tools_setup.md) — external tool setup

## External tools and databases

This pipeline calls third-party tools that are **not included** here. Install them
yourself and follow their licences — several (ConSurf, BepiPred 3.0, DeepLocPro,
DeepTMHMM) are academic / non-commercial only.

Tools: MAFFT, ConSurf (+rate4site), DeepTMHMM, BepiPred 3.0, BLAST+, PyMOL,
optionally DeepLocPro and DiscoTope-3.0.

Databases: NCBI, UniProt, PDB, AlphaFold DB, BacDive, QuickGO (EBI/GOA). Please cite
them when publishing.

## Licence and citation

Original code: [MIT](LICENSE). Third-party tools and databases keep their own terms.

```bibtex
@software{baptista2025targseek,
  title={TargSeek: Surface Protein Discovery and Epitope Prediction Pipeline},
  author={Cátia Baptista},
  year={2025},
  url={https://github.com/CatiaBio/TargSeek}
}
```

See also [CITATION.cff](CITATION.cff).
