#!/usr/bin/env python3
"""
Find downloaded structures that were solved with a binder bound.

A "binder" here is any chain that is not from the target bacterium: a nanobody
(Lama glama, Vicugna pacos), a Fab (Homo sapiens, Mus musculus), or a synthetic
peptide or small molecule. These are the structures where someone has already
shown a molecule gripping the protein, and the residues it grips can be scored
for conservation across our species.

Detection uses the source organism, joined to chains through MOL_ID:

    COMPND  MOL_ID: 2;  MOLECULE: NANOF7;  CHAIN: B
    SOURCE  MOL_ID: 2;  ORGANISM_SCIENTIFIC: LAMA GLAMA

Name matching does not work (a nanobody chain is often called something like
"NANOF7", with no keyword), and neither does chain length (bacterial partners
such as BamE are the same size as a VHH domain).

Output: binder_structures.tsv, one row per target/binder pair.
"""

import glob
import gzip
import os
import re
import sys
from pathlib import Path

# Organisms a binder comes from. A bacterial chain beside one of these is the
# thing being bound.
BINDER_ORGANISM = re.compile(
    r"lama\b|vicugna|camelus|camelid|alpaca|"
    r"homo sapiens|mus musculus|oryctolagus|gallus|"
    r"synthetic construct|unidentified",
    re.I,
)

# Classify what kind of binder it is: they mean different things biologically.
NANOBODY = re.compile(r"nano|sdab|vhh|single.domain", re.I)
ANTIBODY = re.compile(r"\bfab\b|antibody|immunoglob|scfv|heavy chain|light chain", re.I)

COLUMNS = ["pdb_id", "target_chain", "target_molecule", "target_organism",
           "binder_chain", "binder_molecule", "binder_organism", "binder_type"]


def header_lines(path):
    """Yield COMPND/SOURCE lines, stopping at the coordinates."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith(("COMPND", "SOURCE")):
                yield line
            elif line.startswith(("ATOM", "SEQRES", "HETATM")):
                return


def parse_entities(path):
    """Return {mol_id: {'name','chains','organism'}} for one structure."""
    entities = {}
    current = None
    for line in header_lines(path):
        record, body = line[:6], line[10:].strip()

        mol = re.match(r"MOL_ID:\s*(\d+);", body)
        if mol:
            current = mol.group(1)
            entities.setdefault(current, {"name": "", "chains": [], "organism": ""})
            continue
        if current is None:
            continue

        if record == "COMPND":
            name = re.match(r"MOLECULE:\s*(.+?);?$", body)
            if name:
                entities[current]["name"] = name.group(1).rstrip(";").strip()
            chains = re.match(r"CHAIN:\s*(.+?);?$", body)
            if chains:
                entities[current]["chains"] = [
                    c.strip() for c in chains.group(1).rstrip(";").split(",") if c.strip()
                ]
        else:
            org = re.match(r"ORGANISM_SCIENTIFIC:\s*(.+?);?$", body)
            if org:
                entities[current]["organism"] = org.group(1).rstrip(";").strip()
    return entities


def classify(name, organism):
    """nanobody | antibody | peptide_or_small_molecule."""
    if NANOBODY.search(name) or re.search(r"lama\b|vicugna|camel|alpaca", organism, re.I):
        return "nanobody"
    if ANTIBODY.search(name):
        return "antibody"
    return "peptide_or_small_molecule"


def find_binder_pairs(path):
    """Return rows pairing each bacterial chain with each binder chain."""
    try:
        entities = parse_entities(path)
    except (OSError, EOFError) as e:
        print(f"  skipping {path}: {e}")
        return []

    binders, targets = [], []
    for info in entities.values():
        if not info["chains"]:
            continue
        if info["organism"] and BINDER_ORGANISM.search(info["organism"]):
            binders.append(info)
        else:
            targets.append(info)

    if not binders or not targets:
        return []

    pdb_id = Path(str(path).replace(".gz", "")).stem
    rows = []
    for target in targets:
        for binder in binders:
            rows.append({
                "pdb_id": pdb_id,
                "target_chain": ",".join(target["chains"]),
                "target_molecule": target["name"],
                "target_organism": target["organism"],
                "binder_chain": ",".join(binder["chains"]),
                "binder_molecule": binder["name"],
                "binder_organism": binder["organism"],
                "binder_type": classify(binder["name"], binder["organism"]),
            })
    return rows


def scan(structure_dir):
    files = sorted(glob.glob(os.path.join(structure_dir, "*.pdb")) +
                   glob.glob(os.path.join(structure_dir, "*.pdb.gz")))
    seen, rows = set(), []
    for path in files:
        pdb_id = Path(path.replace(".gz", "")).stem
        if pdb_id in seen:
            continue
        seen.add(pdb_id)
        rows.extend(find_binder_pairs(path))
    return rows, len(seen)


def write_table(rows, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")


def main():
    if "snakemake" in globals():
        structure_dir = snakemake.params.structure_dir
        out_path = snakemake.output.binder_table
    else:
        structure_dir = sys.argv[1] if len(sys.argv) > 1 else "data/protein_structures/pdb_files"
        out_path = sys.argv[2] if len(sys.argv) > 2 else "results/phase2/binder_structures.tsv"

    print(f"Scanning {structure_dir} for structures solved with a binder")
    rows, n_scanned = scan(structure_dir)
    write_table(rows, out_path)

    structures = {r["pdb_id"] for r in rows}
    by_type = {}
    for row in rows:
        by_type[row["binder_type"]] = by_type.get(row["binder_type"], 0) + 1

    print(f"Structures scanned      : {n_scanned}")
    print(f"Structures with a binder: {len(structures)}")
    for kind, count in sorted(by_type.items(), key=lambda kv: -kv[1]):
        print(f"  {kind}: {count} pairs")
    print(f"Wrote {len(rows)} target/binder pairs to {out_path}")


if __name__ == "__main__":
    main()
