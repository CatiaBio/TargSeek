#!/usr/bin/env python3
"""
Quality gate for MAFFT alignments.

Checks every alignment in a directory and writes:
  - alignment_qc_gram_{group}.tsv   one row per gene with the metrics
  - alignments_passed_gram_{group}.txt   genes good enough for downstream use

Why: empty or very small alignments used to pass through the pipeline
unnoticed (e.g. a 0-byte alignment still got a ConSurf run). Conservation
scores from a handful of sequences are not meaningful, so such genes are
flagged here instead of producing quiet nonsense later.

Thresholds come from config/config_analysis.yaml:
  quality_control.min_species_per_alignment  minimum number of sequences
  quality_control.min_alignment_length       minimum alignment columns
  quality_control.max_gap_percentage         maximum mean gap fraction
"""

import sys
from pathlib import Path


def read_fasta(path):
    """Yield (header, sequence) pairs. Avoids a Biopython dependency."""
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def assess_alignment(path):
    """Return metrics for one alignment file."""
    records = list(read_fasta(path))
    n_seqs = len(records)
    if n_seqs == 0:
        return {
            "n_sequences": 0,
            "n_columns": 0,
            "mean_gap_fraction": 1.0,
            "ragged": False,
        }

    lengths = {len(seq) for _, seq in records}
    total_residues = sum(len(seq) for _, seq in records)
    total_gaps = sum(seq.count("-") for _, seq in records)

    return {
        "n_sequences": n_seqs,
        "n_columns": max(lengths),
        "mean_gap_fraction": (total_gaps / total_residues) if total_residues else 1.0,
        # An alignment must be rectangular; unequal lengths mean it is broken.
        "ragged": len(lengths) > 1,
    }


def gene_name_from_path(path):
    name = Path(path).stem
    return name[: -len("_aligned")] if name.endswith("_aligned") else name


def check_directory(alignments_dir, min_sequences, min_columns, max_gap_fraction):
    """Assess every .fasta alignment in a directory."""
    alignments_dir = Path(alignments_dir)
    rows = []

    for path in sorted(alignments_dir.glob("*.fasta")):
        metrics = assess_alignment(path)
        reasons = []
        if metrics["n_sequences"] == 0:
            reasons.append("empty_alignment")
        elif metrics["n_sequences"] < min_sequences:
            reasons.append(f"too_few_sequences(<{min_sequences})")
        if metrics["n_columns"] < min_columns:
            reasons.append(f"too_short(<{min_columns})")
        if metrics["mean_gap_fraction"] > max_gap_fraction:
            reasons.append(f"too_gappy(>{max_gap_fraction})")
        if metrics["ragged"]:
            reasons.append("unequal_sequence_lengths")

        rows.append({
            "gene": gene_name_from_path(path),
            "file": path.name,
            "passed": not reasons,
            "reasons": ";".join(reasons) if reasons else "",
            **metrics,
        })

    return rows


def write_report(rows, report_path, passed_path):
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    columns = ["gene", "file", "passed", "n_sequences", "n_columns",
               "mean_gap_fraction", "ragged", "reasons"]
    with open(report_path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(
                f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c])
                for c in columns
            ) + "\n")

    passed_path = Path(passed_path)
    passed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(passed_path, "w") as handle:
        for row in rows:
            if row["passed"]:
                handle.write(f"{row['gene']}\n")


def main():
    if "snakemake" in globals():
        alignments_dir = snakemake.input.alignments_dir
        report_path = snakemake.output.qc_report
        passed_path = snakemake.output.passed_list
        qc = dict(snakemake.params.quality_control or {})
        group = snakemake.params.group
    else:
        if len(sys.argv) < 5:
            print("Usage: check_alignment_quality.py <alignments_dir> "
                  "<qc_report.tsv> <passed_list.txt> <group>")
            sys.exit(1)
        alignments_dir, report_path, passed_path, group = sys.argv[1:5]
        qc = {}

    min_sequences = int(qc.get("min_species_per_alignment", 3))
    min_columns = int(qc.get("min_alignment_length", 50))
    max_gap_fraction = float(qc.get("max_gap_percentage", 0.8))

    print(f"=== Alignment QC (gram_{group}) ===")
    print(f"Directory: {alignments_dir}")
    print(f"Thresholds: >={min_sequences} sequences, >={min_columns} columns, "
          f"<={max_gap_fraction} mean gap fraction")

    rows = check_directory(alignments_dir, min_sequences, min_columns, max_gap_fraction)
    if not rows:
        raise RuntimeError(f"No alignments found in {alignments_dir}")

    write_report(rows, report_path, passed_path)

    failed = [r for r in rows if not r["passed"]]
    print(f"Alignments checked: {len(rows)}")
    print(f"Passed: {len(rows) - len(failed)}")
    print(f"Flagged: {len(failed)}")
    for row in failed:
        print(f"  ! {row['gene']}: {row['reasons']} "
              f"({row['n_sequences']} sequences, {row['n_columns']} columns)")
    print(f"Report: {report_path}")

    if len(failed) == len(rows):
        raise RuntimeError(
            f"Every alignment in {alignments_dir} failed QC. "
            f"See {report_path}."
        )


if __name__ == "__main__":
    main()
