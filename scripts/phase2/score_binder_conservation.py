#!/usr/bin/env python3
"""
Score how conserved each binder's footprint is, and in which species it holds.

Takes the contact residues from extract_binder_contacts.py and answers the
question phase 2 exists for: a nanobody was raised against one strain - would
it still bind the other species we care about?

Two outputs:

  binder_conservation.tsv  one row per contact residue, with its ConSurf grade,
                           score and residue variety
  binder_summary.tsv       one row per binder: how conserved its footprint is
                           overall, and in how many species it is intact

Mapping chain
-------------
ConSurf's grades file already ties PDB residues to alignment positions:

    POS   SEQ   ATOM          SCORE    COLOR   MSA DATA   RESIDUE VARIETY
     27    F    PHE:426:A    -0.320      7     138/139    I,L,V,F,T,S

so a contact at PDB residue 426 maps to position 27, grade 7. The alignment
column is then found by walking the structure's own sequence inside the MSA
(the pipeline adds it there), which lets every other species be read off at the
same column.

Grades run 1-9, 9 being most conserved. Scores are the opposite sign: lower
means more conserved.
"""

import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

CONTACT_COLUMNS = ["gene", "pdb_id", "target_chain", "binder_molecule", "binder_type",
                   "res_seq", "res_name", "min_distance", "consurf_position", "consurf_residue",
                   "consurf_grade", "consurf_score", "msa_coverage", "residue_variety",
                   "species_identical", "species_total", "species_identity_pct"]

SUMMARY_COLUMNS = ["gene", "pdb_id", "binder_molecule", "binder_type", "target_chain",
                   "contact_residues", "mapped_to_conservation", "mean_grade",
                   "highly_conserved", "pct_highly_conserved", "mean_species_identity_pct",
                   "species_fully_conserved", "species_total", "verdict"]

# Conservative substitutions - a different residue that usually preserves binding
SIMILAR = [set("AVLIM"), set("FYW"), set("ST"), set("DE"), set("KR"), set("NQ"), set("GP")]


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_consurf_grades(path):
    """Return {pdb_res_seq: {position, grade, score, coverage, variety}}."""
    mapping = {}
    atom_re = re.compile(r"^([A-Z]{3}):(-?\d+):(\S+)$")
    with open(path, errors="replace") as handle:
        for line in handle:
            fields = [f.strip() for f in line.split("\t") if f.strip()]
            if len(fields) < 5 or not fields[0].isdigit():
                continue
            m = atom_re.match(fields[2])
            if not m:
                continue  # position not resolved in the structure
            try:
                mapping[int(m.group(2))] = {
                    "position": int(fields[0]),
                    "seq": fields[1],
                    "chain": m.group(3),
                    "score": float(fields[3]),
                    "grade": int(re.sub(r"\D", "", fields[4]) or 0),
                    "coverage": fields[5] if len(fields) > 5 else "",
                    "variety": fields[6] if len(fields) > 6 else "",
                }
            except (ValueError, IndexError):
                continue
    return mapping


def read_fasta(path):
    records, name, chunks = [], None, []
    with open(path, errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def column_for_positions(alignment, query_name):
    """Map 1-based ungapped position in the query sequence -> alignment column.

    Only the structure the pipeline selected for the MSA appears in the
    alignment. Other structures of the same protein still have ConSurf
    positions (the grades are in the selected structure's coordinates), so fall
    back to that structure's sequence, which is the one with a bare PDB-style
    header rather than a species header.
    """
    query = next((s for n, s in alignment if n.split()[0].upper() == query_name.upper()), None)
    if query is None:
        query = next((s for n, s in alignment
                      if "|" not in n and len(n.split()[0]) <= 12), None)
    if query is None:
        return None, None
    columns, idx = {}, 0
    for col, residue in enumerate(query):
        if residue != "-":
            idx += 1
            columns[idx] = col
    return columns, query


def species_support(alignment, query_name, column, query_residue):
    """How many species carry the same (or a similar) residue at this column."""
    identical = similar = total = 0
    for name, seq in alignment:
        if name.split()[0].upper() == query_name.upper():
            continue
        if column >= len(seq):
            continue
        residue = seq[column].upper()
        if residue == "-":
            continue
        total += 1
        if residue == query_residue.upper():
            identical += 1
        elif any(residue in group and query_residue.upper() in group for group in SIMILAR):
            similar += 1
    return identical, similar, total


def verdict(pct_conserved, pct_species):
    if pct_conserved >= 70 and pct_species >= 70:
        return "broad - footprint conserved across species"
    if pct_conserved >= 40 or pct_species >= 50:
        return "partial - works in a subset"
    return "narrow - footprint is variable"


def main():
    if "snakemake" in globals():
        contacts_path = snakemake.input.contact_table
        mapping_path = snakemake.input.structure_mapping
        consurf_root = snakemake.params.consurf_root
        msa_root = snakemake.params.msa_root
        out_contacts = snakemake.output.conservation_table
        out_summary = snakemake.output.summary_table
        groups = list(snakemake.params.get("groups", ["positive", "negative"]))
    else:
        args = sys.argv[1:]
        contacts_path = args[0] if args else "results/phase2/binder_contacts.tsv"
        mapping_path = args[1] if len(args) > 1 else None
        consurf_root = args[2] if len(args) > 2 else None
        msa_root = args[3] if len(args) > 3 else None
        out_contacts = args[4] if len(args) > 4 else "results/phase2/binder_conservation.tsv"
        out_summary = args[5] if len(args) > 5 else "results/phase2/binder_summary.tsv"
        groups = ["positive", "negative"]

    contacts = read_tsv(contacts_path)
    print(f"Contact residues to score: {len(contacts)}")

    # structure -> gene, from the final structure mapping
    gene_for_pdb = {}
    if mapping_path and os.path.exists(mapping_path):
        for row in read_tsv(mapping_path):
            sid = str(row.get("structure_id", "")).split("_")[0].upper()
            gene = str(row.get("gene_name", "")).strip()
            if sid and gene:
                gene_for_pdb.setdefault(sid, gene)
    print(f"Structure->gene entries: {len(gene_for_pdb)}")

    grades_cache, alignment_cache = {}, {}

    def load_grades(gene):
        if gene in grades_cache:
            return grades_cache[gene]
        found = None
        for group in groups:
            path = Path(consurf_root) / f"gram_{group}" / gene / "consurf_grades.txt"
            if path.exists():
                found = parse_consurf_grades(path)
                break
        grades_cache[gene] = found
        return found

    def load_alignment(gene):
        if gene in alignment_cache:
            return alignment_cache[gene]
        found = None
        for group in groups:
            path = Path(msa_root) / f"gram_{group}" / f"{gene}_aligned.fasta"
            if path.exists():
                found = read_fasta(path)
                break
        alignment_cache[gene] = found
        return found

    rows, per_binder = [], defaultdict(list)
    unmapped_genes, unmapped_residues = set(), 0

    for contact in contacts:
        pdb = contact["pdb_id"].upper()
        gene = gene_for_pdb.get(pdb, "")
        if not gene:
            unmapped_genes.add(pdb)
            continue

        grades = load_grades(gene)
        if not grades:
            unmapped_genes.add(f"{pdb}({gene}: no ConSurf)")
            continue

        info = grades.get(int(contact["res_seq"]))
        if info is None:
            unmapped_residues += 1
            continue

        identical = total = 0
        alignment = load_alignment(gene)
        if alignment:
            columns, _ = column_for_positions(alignment, pdb)
            if columns and info["position"] in columns:
                identical, _similar, total = species_support(
                    alignment, pdb, columns[info["position"]], info["seq"])

        row = {
            "gene": gene,
            "pdb_id": contact["pdb_id"],
            "target_chain": contact["target_chain"],
            "binder_molecule": contact["binder_molecule"],
            "binder_type": contact["binder_type"],
            "res_seq": contact["res_seq"],
            "res_name": contact["res_name"],
            "min_distance": contact["min_distance"],
            "consurf_position": info["position"],
            "consurf_residue": info["seq"],
            "consurf_grade": info["grade"],
            "consurf_score": f"{info['score']:.3f}",
            "msa_coverage": info["coverage"],
            "residue_variety": info["variety"],
            "species_identical": identical,
            "species_total": total,
            "species_identity_pct": f"{identical / total * 100:.1f}" if total else "",
        }
        rows.append(row)
        per_binder[(gene, contact["pdb_id"], contact["binder_molecule"],
                    contact["binder_type"], contact["target_chain"])].append(row)

    # ---- summaries -------------------------------------------------------
    summaries = []
    for (gene, pdb, molecule, btype, chain), group_rows in sorted(per_binder.items()):
        grades = [r["consurf_grade"] for r in group_rows if r["consurf_grade"]]
        high = [g for g in grades if g >= 7]          # ConSurf 7-9 = conserved
        identities = [float(r["species_identity_pct"]) for r in group_rows
                      if r["species_identity_pct"]]
        full = [i for i in identities if i >= 90]
        pct_high = len(high) / len(grades) * 100 if grades else 0
        mean_identity = sum(identities) / len(identities) if identities else 0
        summaries.append({
            "gene": gene,
            "pdb_id": pdb,
            "binder_molecule": molecule,
            "binder_type": btype,
            "target_chain": chain,
            "contact_residues": len(group_rows),
            "mapped_to_conservation": len(grades),
            "mean_grade": f"{sum(grades) / len(grades):.1f}" if grades else "",
            "highly_conserved": len(high),
            "pct_highly_conserved": f"{pct_high:.0f}",
            "mean_species_identity_pct": f"{mean_identity:.1f}" if identities else "",
            "species_fully_conserved": len(full),
            "species_total": group_rows[0]["species_total"] if group_rows else 0,
            "verdict": verdict(pct_high, mean_identity),
        })

    for path, cols, data in ((out_contacts, CONTACT_COLUMNS, rows),
                             (out_summary, SUMMARY_COLUMNS, summaries)):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as handle:
            handle.write("\t".join(cols) + "\n")
            for item in data:
                handle.write("\t".join(str(item[c]) for c in cols) + "\n")

    print(f"Scored contact residues: {len(rows)}")
    print(f"Binders summarised     : {len(summaries)}")
    if unmapped_residues:
        print(f"Contacts with no ConSurf position (not resolved in structure): {unmapped_residues}")
    if unmapped_genes:
        print(f"Structures without a gene/ConSurf match: {len(unmapped_genes)}")
    print(f"Wrote {out_contacts} and {out_summary}")


if __name__ == "__main__":
    main()
