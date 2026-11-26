# External tool setup

Tools the pipeline calls but does not install. Each gets its own conda
environment, and the pipeline is pointed at it through
`config/config_analysis.yaml` (no hardcoded paths in the scripts).

| Tool | Used for | Config key | Default |
|---|---|---|---|
| ConSurf standalone + rate4site | conservation (stage 3) | `consurf.script_path` | installed |
| BepiPred 3.0 | sequence-based epitopes (stage 4) | `bepipred.path` | installed |
| DeepTMHMM (biolib) | membrane topology (stage 2) | `topology.method: deeptmhmm` | installed |
| DeepLocPro | subcellular localization (stage 2) | `localization.*` | **not installed** |
| DiscoTope-3.0 | structure-based epitopes (stage 4) | `discotope.*` | **not installed** |
| BLAST+ | epitope cross-reactivity | – | in `env.yml` |

Both new tools are **off by default** (`enabled: false`). With them off the
pipeline behaves exactly as before: the localization step writes an empty table
and the surface filter uses topology alone; the DiscoTope step writes an empty
table and the consensus step reports `no_discotope_data`.

---

## DeepLocPro (subcellular localization)

Adds a second line of evidence for "is this protein on the surface?", next to
DeepTMHMM topology. It is trained on prokaryotes and treats Gram-positive and
Gram-negative envelopes separately. License: CC BY-NC-SA 4.0 (non-commercial).

Already done on this machine: the repo is cloned at
`/home/catia/tools/DeepLocPro/repo` and a conda env `deeplocpro` exists with
CPU PyTorch installed.

Remaining step (installs the package and downloads the ESM-2 650M model, ~2.4 GB):

```bash
conda activate deeplocpro
cd /home/catia/tools/DeepLocPro/repo
pip install .
```

Check it works, then enable it:

```bash
conda run -n deeplocpro deeplocpro -f /home/catia/tools/DeepLocPro/repo/test.fasta -o /tmp/dlp_test
```

In `config/config_analysis.yaml`:

```yaml
localization:
  enabled: true
  mode: "intersect"     # topology_only | intersect | union | localization_only
  deeplocpro_bin: "/home/catia/tools/miniconda3/envs/deeplocpro/bin/deeplocpro"
```

`min_confidence` (default 0.7) stops weak predictions overriding topology:
DeepLocPro always returns a class, so a low probability means it had to pick
something. Genes with a low-confidence non-surface call are kept and flagged in
the log for manual review.

`mode` decides how the two lines of evidence combine:

- `topology_only` – ignore localization (previous behaviour)
- `intersect` – surface only if both agree (strictest, fewest targets). A gene
  with no localization prediction keeps its topology call, so nothing is
  dropped silently.
- `union` – surface if either agrees (most targets)
- `localization_only` – localization alone

Output: `stage2_.../localization_predictions/localization_deeplocpro_gram_{group}.tsv`
(with a `confidence` column), and `protein_filter_summary.json` records which
genes each method added or removed.

### Measured effect on the 2025 analysis1 data (2026-11-26)

Topology alone kept **all 40** Gram-negative genes. DeepLocPro called 18 of them
non-surface. With `mode: intersect` and `min_confidence: 0.7`, **10 are dropped**
and 8 are kept with a warning:

| Dropped (confident) | Call | p |
|---|---|---|
| cls | Cytoplasmic Membrane | 1.00 |
| fliQ | Cytoplasmic Membrane | 1.00 |
| fliO | Cytoplasmic Membrane | 0.995 |
| cheZ | Cytoplasmic | 0.986 |
| eno | Cytoplasmic | 0.980 |
| pepN | Cytoplasmic | 0.951 |
| fliP | Cytoplasmic Membrane | 0.924 |
| ycgR | Cytoplasmic | 0.869 |
| fliJ | Cytoplasmic | 0.807 |
| fliH | Cytoplasmic | 0.703 |

These are all textbook cytoplasmic or inner-membrane proteins, i.e. targets an
antibody cannot reach, so topology alone was keeping them wrongly.

Kept but flagged (p < 0.7): fliN 0.69, lolB 0.66, rlpA 0.52, flgK 0.45,
flgF/flgG/fliR 0.42, bamC 0.42. bamC is a real outer-membrane lipoprotein, which
is why the confidence floor exists.

**One judgement call for you:** `eno` (enolase) is confidently cytoplasmic, and
that is its primary location — but enolase is a well-known "moonlighting"
surface protein in several bacteria and a recognised vaccine/diagnostic target.
If you want to keep it, either use `mode: union` or add it back manually.

---

## DiscoTope-3.0 (structure-based epitopes)

Scores residues from the 3D structure (ESM-IF1 inverse folding) rather than
from the sequence, so it uses the structures the pipeline already downloads.
Its performance is reported to be similar on AlphaFold models and solved
structures, and the pipeline passes `--struc_type alphafold` automatically for
`AF-*` structures. Run as a cross-check on BepiPred, not a replacement.

```bash
git clone https://github.com/Magnushhoie/DiscoTope-3.0/ /home/catia/tools/DiscoTope-3.0
cd /home/catia/tools/DiscoTope-3.0
unzip models.zip
conda create -y -n discotope python=3.10
conda activate discotope
pip install -r requirements.txt
pip install .
```

If `torch-scatter` fails to build:
`conda install -c conda-forge pybind11 gcc cxx-compiler` then `pip install pybind11`.

Then in `config/config_analysis.yaml`:

```yaml
discotope:
  enabled: true
  repo_dir: "/home/catia/tools/DiscoTope-3.0"
  python_bin: "/home/catia/tools/miniconda3/envs/discotope/bin/python"
  score_threshold: 0.9          # residue score counted as epitope-like
  consensus_min_fraction: 0.5   # share of an epitope's residues that must clear it
```

Outputs:

- `stage4_.../epitope_prediction_results/gram_{group}/discotope_scores_gram_{group}.tsv`
  – per-residue scores
- `stage4_.../final_reports/gram_{group}/epitope_consensus_gram_{group}.tsv`
  – one row per BepiPred epitope with its DiscoTope score and a verdict:
  `both_agree`, `discotope_partial`, `bepipred_only`, `unmatched_numbering`,
  `no_discotope_data`
- `epitope_consensus_summary_gram_{group}.json` – counts plus the `both_agree`
  epitopes ranked by DiscoTope score

**Numbering:** BepiPred positions and DiscoTope `res_id` need not share an
origin, so the consensus matches each epitope by its peptide sequence and
records the offset it used in `numbering_offset`. Epitopes it cannot place are
marked `unmatched_numbering` instead of being scored against the wrong residues.

---

## Alignment QC (built in, no install)

`check_alignment_quality.py` runs in stage 2 and flags alignments that are
empty, too short, too gappy or ragged. Genes that fail are skipped by ConSurf
instead of producing meaningless scores. Thresholds live under
`quality_control` in `config/config_analysis.yaml`:

```yaml
quality_control:
  min_species_per_alignment: 3    # consider raising: conservation from <10 sequences is noisy
  min_alignment_length: 50
  max_gap_percentage: 0.8
```

Report: `stage2_.../alignment_qc/alignment_qc_gram_{group}.tsv`.

ConSurf now fails the stage unless `consurf.min_success_fraction` of the genes
produce real grades files; it used to mark itself complete even with zero
successes.
