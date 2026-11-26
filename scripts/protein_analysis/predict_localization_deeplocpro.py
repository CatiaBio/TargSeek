#!/usr/bin/env python3
"""
Predict subcellular localization of candidate proteins with DeepLocPro.

Why: the surface-protein filter previously relied on GO terms plus DeepTMHMM
topology. GO annotations are coarse and inconsistent between species, and
topology alone does not distinguish an outer-membrane protein from an inner-
membrane one. DeepLocPro is trained on prokaryotes specifically and knows
Gram-positive and Gram-negative cell envelopes apart.

Input: the per-gene representative 3D structure sequences (same sequences the
topology prediction uses), from selected_3d_paths_gram_{group}.txt.

Output: localization_deeplocpro_gram_{group}.tsv with one row per gene:
    gene, seq_id, localization, surface_exposed, and one column per class.

DeepLocPro is a separate installation (it needs PyTorch and ESM-2). Set its
executable in config/config_analysis.yaml under localization.deeplocpro_bin.
With localization.enabled = false this writes an empty, clearly-marked table so
the rest of the pipeline still runs.
"""

import os
import subprocess
import sys
from pathlib import Path

CLASSES = [
    "Cell wall & surface",
    "Extracellular",
    "Cytoplasmic",
    "Cytoplasmic Membrane",
    "Outer Membrane",
    "Periplasmic",
]

# Classes that put at least part of the protein on the cell surface, where an
# antibody could reach it. Overridable via localization.surface_classes.
DEFAULT_SURFACE_CLASSES = ["Cell wall & surface", "Extracellular", "Outer Membrane"]

OUTPUT_COLUMNS = ["gene", "seq_id", "localization", "surface_exposed"] + CLASSES


def parse_selected_3d_paths(paths_file):
    """Map gene -> structure FASTA path (same convention as predict_topology.py)."""
    gene_fasta_map = {}
    paths_file = Path(paths_file)
    if not paths_file.exists():
        print(f"Warning: selected 3D paths file not found: {paths_file}")
        return gene_fasta_map

    for line in paths_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("/")
        if len(parts) >= 3:
            gene_fasta_map[parts[-2]] = line
    print(f"Loaded {len(gene_fasta_map)} gene-structure mappings from {paths_file}")
    return gene_fasta_map


def read_first_sequence(fasta_path):
    """Return the first sequence in a FASTA file, gaps removed."""
    seq_chunks, started = [], False
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif started:
                seq_chunks.append(line.strip())
    return "".join(seq_chunks).replace("-", "").upper()


def build_input_fasta(gene_fasta_map, fasta_out):
    """Write one representative sequence per gene. Returns {seq_id: gene}."""
    fasta_out = Path(fasta_out)
    fasta_out.parent.mkdir(parents=True, exist_ok=True)

    seq_id_to_gene = {}
    with open(fasta_out, "w") as handle:
        for gene, path in sorted(gene_fasta_map.items()):
            if not os.path.exists(path):
                print(f"Warning: missing structure FASTA for {gene}: {path}")
                continue
            sequence = read_first_sequence(path)
            if not sequence:
                print(f"Warning: no sequence found for {gene} in {path}")
                continue
            pdb_id = Path(path).name.replace(".fasta", "")
            # DeepLocPro slugifies IDs, so keep them simple.
            seq_id = f"{gene}_{pdb_id}".replace("|", "_").replace(" ", "_")
            seq_id_to_gene[seq_id] = gene
            handle.write(f">{seq_id}\n{sequence}\n")

    print(f"Wrote {len(seq_id_to_gene)} sequences to {fasta_out}")
    return seq_id_to_gene


def run_deeplocpro(deeplocpro_bin, fasta_in, out_dir, group, device="cpu"):
    """Run DeepLocPro and return the path of the CSV it produced."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # DeepLocPro names its output results_<timestamp>.csv, so note what exists
    # beforehand to identify the new file.
    before = set(out_dir.glob("results_*.csv"))

    cmd = [
        str(deeplocpro_bin),
        "-f", str(Path(fasta_in).resolve()),
        "-o", str(out_dir.resolve()),
        "-g", group,          # "positive" or "negative"
        "-d", device,
    ]
    print("Running: " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout.strip():
        print(result.stdout.strip()[-4000:])

    if result.returncode != 0:
        raise RuntimeError(
            f"DeepLocPro failed (exit {result.returncode}) for gram_{group}.\n"
            f"{(result.stderr or result.stdout).strip()[-4000:]}"
        )

    produced = sorted(set(out_dir.glob("results_*.csv")) - before,
                      key=lambda p: p.stat().st_mtime)
    if not produced:
        # Fall back to the newest file in case the run reused a name.
        existing = sorted(out_dir.glob("results_*.csv"), key=lambda p: p.stat().st_mtime)
        if not existing:
            raise RuntimeError(
                f"DeepLocPro exited 0 but wrote no results_*.csv in {out_dir}"
            )
        produced = existing
    return produced[-1]


def write_table(rows, out_path, note=None):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        if note:
            handle.write(f"# {note}\n")
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(
                f"{row.get(c, ''):.4f}" if isinstance(row.get(c), float) else str(row.get(c, ""))
                for c in OUTPUT_COLUMNS
            ) + "\n")
    print(f"Wrote {len(rows)} rows to {out_path}")


def collect_rows(csv_path, seq_id_to_gene, surface_classes):
    """Turn the DeepLocPro CSV into our per-gene rows."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    if "ACC" not in df.columns or "Localization" not in df.columns:
        raise RuntimeError(
            f"Unexpected DeepLocPro output columns in {csv_path}: {list(df.columns)}"
        )

    rows = []
    for _, record in df.iterrows():
        seq_id = str(record["ACC"])
        localization = str(record["Localization"])
        row = {
            "gene": seq_id_to_gene.get(seq_id, seq_id.split("_")[0]),
            "seq_id": seq_id,
            "localization": localization,
            "surface_exposed": localization in surface_classes,
        }
        for cls in CLASSES:
            if cls in df.columns:
                try:
                    row[cls] = float(record[cls])
                except (TypeError, ValueError):
                    row[cls] = ""
        rows.append(row)
    return rows


def main():
    if "snakemake" in globals():
        selected_3d_paths = snakemake.input.selected_3d_paths
        out_table = snakemake.output.localization_table
        work_dir = Path(snakemake.output.localization_table).parent / f"deeplocpro_gram_{snakemake.params.group}"
        group = snakemake.params.group
        cfg = dict(snakemake.params.localization or {})
    else:
        if len(sys.argv) < 4:
            print("Usage: predict_localization_deeplocpro.py "
                  "<selected_3d_paths.txt> <out_table.tsv> <positive|negative> "
                  "[deeplocpro_bin]")
            sys.exit(1)
        selected_3d_paths, out_table, group = sys.argv[1:4]
        cfg = {"enabled": True}
        if len(sys.argv) > 4:
            cfg["deeplocpro_bin"] = sys.argv[4]
        work_dir = Path(out_table).parent / f"deeplocpro_gram_{group}"

    enabled = bool(cfg.get("enabled", False))
    deeplocpro_bin = cfg.get("deeplocpro_bin", "deeplocpro")
    device = cfg.get("device", "cpu")
    surface_classes = list(cfg.get("surface_classes", DEFAULT_SURFACE_CLASSES))

    print(f"=== DeepLocPro localization (gram_{group}) ===")

    if not enabled:
        write_table([], out_table,
                    note="DeepLocPro disabled (localization.enabled = false). "
                         "No localization evidence was used.")
        print("DeepLocPro is disabled in the config; wrote an empty table.")
        return

    gene_fasta_map = parse_selected_3d_paths(selected_3d_paths)
    if not gene_fasta_map:
        write_table([], out_table,
                    note=f"No genes with 3D structures in {selected_3d_paths}")
        print("No input genes; wrote an empty table.")
        return

    work_dir.mkdir(parents=True, exist_ok=True)
    fasta_in = work_dir / f"deeplocpro_input_gram_{group}.fasta"
    seq_id_to_gene = build_input_fasta(gene_fasta_map, fasta_in)
    if not seq_id_to_gene:
        raise RuntimeError(f"Could not build any input sequences from {selected_3d_paths}")

    csv_path = run_deeplocpro(deeplocpro_bin, fasta_in, work_dir, group, device)
    print(f"DeepLocPro results: {csv_path}")

    rows = collect_rows(csv_path, seq_id_to_gene, surface_classes)
    write_table(rows, out_table,
                note=f"DeepLocPro predictions; surface classes = {surface_classes}")

    n_surface = sum(1 for r in rows if r["surface_exposed"])
    print(f"Surface-exposed by localization: {n_surface}/{len(rows)}")
    counts = {}
    for row in rows:
        counts[row["localization"]] = counts.get(row["localization"], 0) + 1
    for cls, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {cls}: {n}")


if __name__ == "__main__":
    main()
