#!/usr/bin/env python3
"""
Cross-check BepiPred 3.0 epitopes against DiscoTope-3.0 structure scores.

Why: BepiPred predicts from sequence, DiscoTope from structure. Published
benchmarks disagree about which is better (the DiscoTope authors report a clear
advantage, an independent benchmark finds them nearly tied), so the defensible
use is agreement: an epitope both methods like is a better candidate than one
only a single method likes. This script does not replace either tool's output.

Input:
  - BepiPred linear epitope tables: {gene}/{gene}_{structure}_linear_epitopes.tsv
  - DiscoTope residue scores: discotope_scores_gram_{group}.tsv

Output: epitope_consensus_gram_{group}.tsv, one row per BepiPred epitope with
the DiscoTope score over the same residues and an agreement verdict.

Residue numbering: BepiPred's Start/End and DiscoTope's res_id need not use the
same origin, so the epitope's own peptide sequence is used to find the matching
stretch of DiscoTope residues. The offset actually used is reported in the
numbering_offset column; rows where no match was found are marked
unmatched_numbering rather than silently scored.
"""

import json
import os
import sys
from pathlib import Path

OUTPUT_COLUMNS = [
    "gene", "structure_id", "chain", "epitope_no", "start", "end", "peptide",
    "n_residues", "bepipred_score", "discotope_mean", "discotope_max",
    "discotope_residues_above_threshold", "fraction_above_threshold",
    "numbering_offset", "agreement",
]


def read_bepipred_tables(epitope_dir):
    """Read all BepiPred linear epitope tables under a directory."""
    rows = []
    for path in sorted(Path(epitope_dir).rglob("*_linear_epitopes.tsv")):
        # File name: {gene}_{structure}_linear_epitopes.tsv
        stem = path.name.replace("_linear_epitopes.tsv", "")
        gene = path.parent.name
        structure_id = stem[len(gene) + 1:] if stem.startswith(gene + "_") else stem

        header = None
        for line in path.read_text().splitlines():
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            record = dict(zip(header, fields))
            try:
                rows.append({
                    "gene": gene,
                    "structure_id": structure_id,
                    "chain": (record.get("Chain") or "").strip(),
                    "epitope_no": record.get("No.", ""),
                    "start": int(record["Start"]),
                    "end": int(record["End"]),
                    "peptide": (record.get("Peptide") or "").strip().upper(),
                    "n_residues": record.get("Number_of_Residues", ""),
                    "bepipred_score": record.get("Score", ""),
                })
            except (KeyError, ValueError):
                print(f"Warning: could not parse a row in {path}: {record}")
    return rows


def read_discotope(table_path):
    """Return {(structure_id_upper, chain): {res_id: (residue, score)}}."""
    scores, header = {}, None
    path = Path(table_path)
    if not path.exists():
        return scores

    for line in path.read_text().splitlines():
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        record = dict(zip(header, fields))
        try:
            res_id = int(float(record["res_id"]))
            score = float(record["discotope_score"])
        except (KeyError, ValueError):
            continue
        key = (str(record.get("structure_id", "")).upper(),
               str(record.get("chain", "")).strip())
        scores.setdefault(key, {})[res_id] = (
            str(record.get("residue", "")).strip().upper(), score
        )
    return scores


def find_offset(chain_scores, start, peptide):
    """Find the numbering offset where the peptide matches DiscoTope residues.

    Returns (offset, matched_residue_ids) or (None, []).
    """
    if not peptide or not chain_scores:
        return None, []

    res_ids = sorted(chain_scores)

    def matches(offset):
        ids = [start + offset + i for i in range(len(peptide))]
        if not all(i in chain_scores for i in ids):
            return None
        observed = "".join(chain_scores[i][0] for i in ids)
        return ids if observed == peptide else None

    # Try the numbering as given first, then a search over plausible offsets.
    ids = matches(0)
    if ids:
        return 0, ids

    candidates = sorted(
        {res_id - start for res_id in res_ids},
        key=lambda off: abs(off),
    )
    for offset in candidates:
        ids = matches(offset)
        if ids:
            return offset, ids
    return None, []


def summarise(chain_scores, matched_ids, threshold):
    values = [chain_scores[i][1] for i in matched_ids]
    if not values:
        return {}
    above = [v for v in values if v >= threshold]
    return {
        "discotope_mean": round(sum(values) / len(values), 4),
        "discotope_max": round(max(values), 4),
        "discotope_residues_above_threshold": len(above),
        "fraction_above_threshold": round(len(above) / len(values), 4),
    }


def classify(row, threshold, min_fraction):
    """Agreement verdict for one epitope."""
    if row.get("agreement") == "no_discotope_data":
        return "no_discotope_data"
    if row.get("numbering_offset") == "":
        return "unmatched_numbering"
    fraction = row.get("fraction_above_threshold", 0.0)
    mean = row.get("discotope_mean", 0.0)
    if fraction >= min_fraction:
        return "both_agree"
    if mean >= threshold * 0.75:
        return "discotope_partial"
    return "bepipred_only"


def write_table(rows, out_path, note=None):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        if note:
            handle.write(f"# {note}\n")
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(c, "")) for c in OUTPUT_COLUMNS) + "\n")
    print(f"Wrote {len(rows)} epitope rows to {out_path}")


def main():
    if "snakemake" in globals():
        epitope_dir = Path(snakemake.input.epitope_tables_sentinel).parent
        discotope_table = snakemake.input.discotope_table
        out_table = snakemake.output.consensus_table
        out_summary = snakemake.output.consensus_summary
        group = snakemake.params.group
        cfg = dict(snakemake.params.discotope or {})
    else:
        if len(sys.argv) < 6:
            print("Usage: epitope_consensus_bepipred_discotope.py <epitope_dir> "
                  "<discotope_table.tsv> <out_table.tsv> <out_summary.json> "
                  "<positive|negative>")
            sys.exit(1)
        epitope_dir, discotope_table, out_table, out_summary, group = sys.argv[1:6]
        cfg = {}

    threshold = float(cfg.get("score_threshold", 0.9))
    min_fraction = float(cfg.get("consensus_min_fraction", 0.5))

    print(f"=== BepiPred / DiscoTope consensus (gram_{group}) ===")
    print(f"Epitope tables: {epitope_dir}")
    print(f"DiscoTope scores: {discotope_table}")
    print(f"Thresholds: score >= {threshold}, >= {min_fraction:.0%} of epitope residues")

    epitopes = read_bepipred_tables(epitope_dir)
    print(f"BepiPred epitopes found: {len(epitopes)}")

    discotope = read_discotope(discotope_table)
    print(f"DiscoTope chains with scores: {len(discotope)}")

    rows = []
    for epitope in epitopes:
        row = dict(epitope)
        row["numbering_offset"] = ""

        chain_scores = discotope.get((row["structure_id"].upper(), row["chain"]))
        if chain_scores is None:
            # Fall back to the structure's only chain if the label differs
            same_structure = [
                v for (sid, _), v in discotope.items()
                if sid == row["structure_id"].upper()
            ]
            chain_scores = same_structure[0] if len(same_structure) == 1 else None

        if not chain_scores:
            row["agreement"] = "no_discotope_data"
            rows.append(row)
            continue

        offset, matched_ids = find_offset(chain_scores, row["start"], row["peptide"])
        if offset is None:
            row["agreement"] = "unmatched_numbering"
            rows.append(row)
            continue

        row["numbering_offset"] = offset
        row.update(summarise(chain_scores, matched_ids, threshold))
        row["agreement"] = classify(row, threshold, min_fraction)
        rows.append(row)

    write_table(rows, out_table,
                note=(f"BepiPred epitopes cross-checked against DiscoTope-3.0 "
                      f"(score >= {threshold}, >= {min_fraction:.0%} of residues)"))

    counts = {}
    for row in rows:
        counts[row["agreement"]] = counts.get(row["agreement"], 0) + 1

    agreed = [r for r in rows if r["agreement"] == "both_agree"]
    summary = {
        "group": group,
        "thresholds": {"discotope_score": threshold, "min_fraction": min_fraction},
        "total_bepipred_epitopes": len(rows),
        "verdict_counts": counts,
        "consensus_epitopes": sorted(
            ({
                "gene": r["gene"],
                "structure_id": r["structure_id"],
                "peptide": r["peptide"],
                "start": r["start"],
                "end": r["end"],
                "bepipred_score": r["bepipred_score"],
                "discotope_mean": r.get("discotope_mean"),
            } for r in agreed),
            key=lambda d: -(d["discotope_mean"] or 0),
        ),
    }
    Path(out_summary).parent.mkdir(parents=True, exist_ok=True)
    with open(out_summary, "w") as handle:
        json.dump(summary, handle, indent=2)

    print("Verdicts:")
    for verdict, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {verdict}: {n}")
    print(f"Summary: {out_summary}")

    if counts.get("no_discotope_data", 0) == len(rows) and rows:
        print("NOTE: no DiscoTope scores were available, so no epitope could be "
              "cross-checked. Is discotope.enabled set to true?")


if __name__ == "__main__":
    main()
