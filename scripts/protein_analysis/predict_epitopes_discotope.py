#!/usr/bin/env python3
"""
Structure-based B-cell epitope prediction with DiscoTope-3.0.

Why: the pipeline already selects a 3D structure per protein but predicts
epitopes from sequence only (BepiPred 3.0). DiscoTope-3.0 scores residues from
the structure itself (ESM-IF1 inverse-folding representations) and reports
similar performance on AlphaFold models as on solved structures, so the
structures we already download can be used directly. It is run alongside
BepiPred, not instead of it: see epitope_consensus_bepipred_discotope.py.

Input: the final FASTA-structure mapping table (gene -> structure + chain) and
the surface-protein filter summary, so only surface proteins are scored.

Output: discotope_scores_gram_{group}.tsv, one row per residue:
    gene, structure_id, chain, res_id, residue, discotope_score, rsa,
    predicted_epitope, struc_type

DiscoTope-3.0 is a separate install (PyTorch + PyTorch Geometric). Point
config/config_analysis.yaml at it under discotope.*. With discotope.enabled =
false this writes an empty, clearly-marked table.
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

OUTPUT_COLUMNS = ["gene", "structure_id", "chain", "res_id", "residue",
                  "discotope_score", "rsa", "predicted_epitope", "struc_type"]


def load_surface_genes(summary_path, group):
    """Genes that passed the surface filter for this group."""
    if not summary_path or not os.path.exists(summary_path):
        return None
    with open(summary_path) as handle:
        summary = json.load(handle)
    genes = summary.get("surface_genes_by_group", {}).get(group)
    return set(genes) if genes is not None else None


def is_alphafold(structure_id):
    """AlphaFold model IDs look like AF-Q8EM16; PDB IDs are 4 characters."""
    return str(structure_id).upper().startswith("AF-")


def load_selected_structures(mapping_path, selected_paths_file, surface_genes):
    """Pick one structure per gene: the one the pipeline selected for the MSA.

    The selected_3d_paths file lists the chosen FASTA per gene; the mapping
    table carries the structure path and chain for it.
    """
    import pandas as pd

    df = pd.read_csv(mapping_path, sep="\t")
    required = {"gene_name", "structure_id", "chain", "structure_path", "fasta_path"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"{mapping_path} is missing columns: {sorted(missing)}")

    # The FASTA paths the pipeline actually selected
    selected_fastas = set()
    selected_paths_file = Path(selected_paths_file)
    if selected_paths_file.exists():
        for line in selected_paths_file.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                selected_fastas.add(line)
                selected_fastas.add(os.path.normpath(line))

    chosen = {}
    for _, row in df.iterrows():
        gene = str(row["gene_name"])
        if surface_genes is not None and gene not in surface_genes:
            continue
        fasta_path = str(row["fasta_path"])
        if selected_fastas and not (
            fasta_path in selected_fastas or os.path.normpath(fasta_path) in selected_fastas
        ):
            continue
        if gene in chosen:
            continue
        chain = row["chain"]
        chain = "A" if (pd.isna(chain) or str(chain).strip() == "") else str(chain).strip()
        chosen[gene] = {
            "structure_id": str(row["structure_id"]),
            "chain": chain,
            "structure_path": str(row["structure_path"]),
        }

    return chosen


def stage_structure(structure_path, structure_id, dest_dir):
    """Copy (and gunzip) one structure into dest_dir as <structure_id>.pdb."""
    src = Path(structure_path)
    if not src.exists():
        alt = Path(str(src) + ".gz")
        if alt.exists():
            src = alt
        else:
            print(f"Warning: structure file not found: {structure_path}")
            return None

    dest = Path(dest_dir) / f"{structure_id}.pdb"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rb") as f_in, open(dest, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy(src, dest)
    return dest


def run_discotope(repo_dir, pdb_dir, out_dir, struc_type, cpu_only=True,
                  python_bin=None, extra_args=None):
    """Run DiscoTope-3.0 over a directory of PDB files."""
    repo_dir = Path(repo_dir)
    main_script = repo_dir / "discotope3" / "main.py"
    if not main_script.exists():
        raise FileNotFoundError(
            f"DiscoTope main.py not found at {main_script}. "
            f"Set discotope.repo_dir in config/config_analysis.yaml."
        )

    cmd = [
        python_bin or sys.executable, str(main_script),
        "--pdb_dir", str(Path(pdb_dir).resolve()),
        "--out_dir", str(Path(out_dir).resolve()),
        "--struc_type", struc_type,
    ]
    if cpu_only:
        cmd.append("--cpu_only")
    if extra_args:
        cmd.extend(extra_args)

    print("Running: " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(repo_dir))
    if result.stdout.strip():
        print(result.stdout.strip()[-4000:])
    if result.returncode != 0:
        raise RuntimeError(
            f"DiscoTope-3.0 failed (exit {result.returncode}) for struc_type={struc_type}.\n"
            f"{(result.stderr or result.stdout).strip()[-4000:]}"
        )
    return Path(out_dir)


def collect_scores(out_dir, structure_to_gene, struc_type, threshold):
    """Read every *_discotope3.csv under out_dir into our row format."""
    import pandas as pd

    rows = []
    csv_files = sorted(Path(out_dir).rglob("*_discotope3.csv"))
    print(f"Found {len(csv_files)} DiscoTope CSV files in {out_dir}")

    for csv_path in csv_files:
        df = pd.read_csv(csv_path)
        score_col = next((c for c in df.columns if c.startswith("DiscoTope-3.0_score")), None)
        if score_col is None or "res_id" not in df.columns:
            print(f"Warning: unexpected columns in {csv_path}: {list(df.columns)}")
            continue

        # The "pdb" column is <structure>_<chain>, e.g. 7c4s_A
        for chain_label, chain_df in df.groupby(df.get("pdb", csv_path.stem)):
            label = str(chain_label)
            structure_id, _, chain = label.rpartition("_")
            if not structure_id:
                structure_id, chain = label, ""
            gene = structure_to_gene.get(structure_id.upper(), "")
            if not gene:
                gene = structure_to_gene.get(structure_id, "")

            for _, record in chain_df.iterrows():
                score = record[score_col]
                try:
                    score = float(score)
                except (TypeError, ValueError):
                    continue
                rows.append({
                    "gene": gene,
                    "structure_id": structure_id,
                    "chain": chain,
                    "res_id": record["res_id"],
                    "residue": record.get("residue", ""),
                    "discotope_score": score,
                    "rsa": record.get("rsa", ""),
                    "predicted_epitope": score >= threshold,
                    "struc_type": struc_type,
                })
    return rows


def write_table(rows, out_path, note=None):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        if note:
            handle.write(f"# {note}\n")
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(c, "")) for c in OUTPUT_COLUMNS) + "\n")
    print(f"Wrote {len(rows)} residue rows to {out_path}")


def main():
    if "snakemake" in globals():
        mapping_path = snakemake.input.structure_mapping
        selected_paths_file = snakemake.input.selected_3d_paths
        summary_path = getattr(snakemake.input, "protein_filter_summary", None)
        out_table = snakemake.output.discotope_table
        group = snakemake.params.group
        cfg = dict(snakemake.params.discotope or {})
    else:
        if len(sys.argv) < 5:
            print("Usage: predict_epitopes_discotope.py <mapping.tsv> "
                  "<selected_3d_paths.txt> <out_table.tsv> <positive|negative> "
                  "[repo_dir]")
            sys.exit(1)
        mapping_path, selected_paths_file, out_table, group = sys.argv[1:5]
        summary_path = None
        cfg = {"enabled": True}
        if len(sys.argv) > 5:
            cfg["repo_dir"] = sys.argv[5]

    print(f"=== DiscoTope-3.0 epitope prediction (gram_{group}) ===")

    if not bool(cfg.get("enabled", False)):
        write_table([], out_table,
                    note="DiscoTope-3.0 disabled (discotope.enabled = false).")
        print("DiscoTope-3.0 is disabled in the config; wrote an empty table.")
        return

    repo_dir = cfg.get("repo_dir", "")
    python_bin = cfg.get("python_bin") or None
    cpu_only = bool(cfg.get("cpu_only", True))
    threshold = float(cfg.get("score_threshold", 0.9))
    keep_work = bool(cfg.get("keep_intermediate", False))

    surface_genes = load_surface_genes(summary_path, group)
    if surface_genes is not None:
        print(f"Surface genes for gram_{group}: {len(surface_genes)}")

    chosen = load_selected_structures(mapping_path, selected_paths_file, surface_genes)
    if not chosen:
        write_table([], out_table,
                    note=f"No selected surface-protein structures for gram_{group}")
        print("No structures to score; wrote an empty table.")
        return
    print(f"Structures to score: {len(chosen)}")

    work_root = Path(out_table).parent / f"discotope_work_gram_{group}"
    if work_root.exists():
        shutil.rmtree(work_root)
    work_root.mkdir(parents=True, exist_ok=True)

    # DiscoTope needs one run per structure type, so stage them separately.
    structure_to_gene, staged = {}, {"solved": 0, "alphafold": 0}
    for gene, info in sorted(chosen.items()):
        struc_type = "alphafold" if is_alphafold(info["structure_id"]) else "solved"
        pdb_dir = work_root / f"input_{struc_type}"
        dest = stage_structure(info["structure_path"], info["structure_id"], pdb_dir)
        if dest is None:
            continue
        structure_to_gene[info["structure_id"].upper()] = gene
        structure_to_gene[info["structure_id"]] = gene
        staged[struc_type] += 1
    print(f"Staged structures: {staged['solved']} solved, {staged['alphafold']} AlphaFold")

    rows = []
    for struc_type in ("solved", "alphafold"):
        pdb_dir = work_root / f"input_{struc_type}"
        if not pdb_dir.exists() or not any(pdb_dir.glob("*.pdb")):
            continue
        out_dir = work_root / f"output_{struc_type}"
        run_discotope(repo_dir, pdb_dir, out_dir, struc_type,
                      cpu_only=cpu_only, python_bin=python_bin,
                      extra_args=cfg.get("extra_args"))
        rows.extend(collect_scores(out_dir, structure_to_gene, struc_type, threshold))

    if not rows:
        raise RuntimeError(
            f"DiscoTope-3.0 produced no scores for gram_{group}. "
            f"Check the logs above and {work_root}."
        )

    write_table(rows, out_table,
                note=f"DiscoTope-3.0 scores; predicted_epitope at score >= {threshold}")

    genes_scored = {r["gene"] for r in rows if r["gene"]}
    n_epitope_residues = sum(1 for r in rows if r["predicted_epitope"])
    print(f"Genes scored: {len(genes_scored)}")
    print(f"Residues above threshold: {n_epitope_residues}/{len(rows)}")

    if not keep_work:
        for struc_type in ("solved", "alphafold"):
            shutil.rmtree(work_root / f"input_{struc_type}", ignore_errors=True)


if __name__ == "__main__":
    main()
