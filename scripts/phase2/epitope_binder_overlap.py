#!/usr/bin/env python3
"""
Does a predicted epitope sit where a known binder actually grips the protein?

Phase 1 predicts epitopes; phase 2 measures where real nanobodies, Fabs and
peptides bind. This joins the two: for each epitope, how much of it falls
inside a binder's footprint.

An epitope that overlaps a binder footprint is supported by experimental
evidence that the site is reachable and bindable - the strongest result the
pipeline can produce. An epitope with no overlap is not disproved; it simply
has no structure showing anything bound there.

Both sides use the same PDB residue numbering, which is what makes the join
safe: epitope Start/End come from the structure sequence, and contacts are read
from that structure's coordinates.

Output: epitope_binder_overlap.tsv, one row per epitope/binder pair that
overlaps, plus a row for epitopes with no overlap so they are not lost.
"""

import csv
import glob
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

COLUMNS = ["gene", "pdb_id", "epitope_no", "epitope_start", "epitope_end", "peptide",
           "epitope_length", "bepipred_score", "binder_molecule", "binder_type",
           "overlap_residues", "overlap_pct", "overlapped_positions",
           "mean_consurf_grade_of_overlap", "mean_species_identity_of_overlap", "note"]


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_epitopes(root):
    """Read every *_filtered_epitopes.csv under root -> list of epitope dicts."""
    epitopes = []
    for path in sorted(glob.glob(os.path.join(root, "**", "*_filtered_epitopes.csv"),
                                 recursive=True)):
        stem = Path(path).name.replace("_filtered_epitopes.csv", "")
        # names look like gene_gene_PDB_1 -> take gene and the PDB code
        parts = stem.split("_")
        gene = parts[0]
        pdb = next((p for p in parts if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", p)), "")
        with open(path, errors="replace") as handle:
            for row in csv.DictReader(handle):
                try:
                    start, end = int(row["Start"]), int(row["End"])
                except (KeyError, ValueError):
                    continue
                epitopes.append({
                    "gene": gene,
                    "pdb_id": pdb.upper(),
                    "epitope_no": row.get("No.", ""),
                    "start": start,
                    "end": end,
                    "peptide": row.get("Peptide", ""),
                    "score": row.get("Score", ""),
                })
    return epitopes


def main():
    if "snakemake" in globals():
        conservation_path = snakemake.input.conservation_table
        epitope_root = snakemake.params.epitope_root
        out_path = snakemake.output.overlap_table
    else:
        args = sys.argv[1:]
        conservation_path = args[0] if args else "results/phase2/binder_conservation.tsv"
        epitope_root = args[1] if len(args) > 1 else (
            "results/analysis2_params2/protein_analysis/"
            "stage4_epitope_predictions_and_reports/epitope_prediction_results")
        out_path = args[2] if len(args) > 2 else "results/phase2/epitope_binder_overlap.tsv"

    contacts = read_tsv(conservation_path)
    epitopes = load_epitopes(epitope_root)
    print(f"Epitopes loaded        : {len(epitopes)}")
    print(f"Scored contact residues: {len(contacts)}")

    # footprints keyed by gene: {(gene, pdb, binder, type): {res_seq: row}}
    footprints = defaultdict(dict)
    for row in contacts:
        key = (row["gene"], row["pdb_id"], row["binder_molecule"], row["binder_type"])
        # Epitope Start/End and ConSurf positions are both 1-based indices into
        # the structure's SEQRES sequence; res_seq is PDB author numbering
        # (bamA: position 27 is residue 426). Join on the sequence position.
        try:
            footprints[key][int(row["consurf_position"])] = row
        except (ValueError, KeyError):
            continue

    rows = []
    genes_with_binders = {k[0] for k in footprints}

    for epitope in epitopes:
        span = set(range(epitope["start"], epitope["end"] + 1))
        length = len(span)
        matched = False

        for (gene, pdb, binder, btype), residues in footprints.items():
            if gene != epitope["gene"]:
                continue
            hit = sorted(span & set(residues))
            if not hit:
                continue
            matched = True
            grades = [int(residues[p]["consurf_grade"]) for p in hit
                      if residues[p]["consurf_grade"]]
            idents = [float(residues[p]["species_identity_pct"]) for p in hit
                      if residues[p]["species_identity_pct"]]
            rows.append({
                "gene": gene,
                "pdb_id": epitope["pdb_id"] or pdb,
                "epitope_no": epitope["epitope_no"],
                "epitope_start": epitope["start"],
                "epitope_end": epitope["end"],
                "peptide": epitope["peptide"],
                "epitope_length": length,
                "bepipred_score": epitope["score"],
                "binder_molecule": binder,
                "binder_type": btype,
                "overlap_residues": len(hit),
                "overlap_pct": f"{len(hit) / length * 100:.0f}",
                "overlapped_positions": ",".join(str(p) for p in hit),
                "mean_consurf_grade_of_overlap":
                    f"{sum(grades) / len(grades):.1f}" if grades else "",
                "mean_species_identity_of_overlap":
                    f"{sum(idents) / len(idents):.1f}" if idents else "",
                "note": "",
            })

        if not matched:
            note = ("no binder structure for this protein"
                    if epitope["gene"] not in genes_with_binders
                    else "no overlap with any known binder footprint")
            rows.append({
                "gene": epitope["gene"],
                "pdb_id": epitope["pdb_id"],
                "epitope_no": epitope["epitope_no"],
                "epitope_start": epitope["start"],
                "epitope_end": epitope["end"],
                "peptide": epitope["peptide"],
                "epitope_length": length,
                "bepipred_score": epitope["score"],
                "binder_molecule": "", "binder_type": "",
                "overlap_residues": 0, "overlap_pct": "0",
                "overlapped_positions": "",
                "mean_consurf_grade_of_overlap": "",
                "mean_species_identity_of_overlap": "",
                "note": note,
            })

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")

    overlapping = [r for r in rows if int(r["overlap_residues"]) > 0]
    epitopes_hit = {(r["gene"], r["epitope_start"], r["epitope_end"]) for r in overlapping}
    print(f"Epitope/binder overlaps : {len(overlapping)}")
    print(f"Distinct epitopes hit   : {len(epitopes_hit)}")
    print(f"Genes with binder data  : {sorted(genes_with_binders)}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
