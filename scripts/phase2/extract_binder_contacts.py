#!/usr/bin/env python3
"""
Find the residues on a target protein that a bound binder actually touches.

For every target/binder pair from scan_structures_for_binders.py, measure
heavy-atom distances between the two chains and record every target residue
within the cutoff. That set of residues is the binder's footprint: the patch of
protein a nanobody, Fab or peptide grips.

The footprint is what later steps score for conservation - if the residues a
binder depends on vary across species, the binder will only work in some of
them.

Hydrogens are ignored (often absent from crystal structures), and waters and
ions are skipped. The default 4.5 A cutoff is the usual range for a contact
including hydrogen bonds and salt bridges.

Output: binder_contacts.tsv, one row per contacted target residue.
"""

import gzip
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

SKIP_RESIDUES = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD",
                 "NA", "CL", "MG", "CA", "ZN", "K", "MN", "FE", "CD", "NI"}

COLUMNS = ["pdb_id", "target_chain", "binder_chain", "binder_molecule",
           "binder_type", "res_seq", "res_name", "min_distance", "n_atom_contacts"]


def open_structure(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", errors="replace")


def load_atoms(path, chains):
    """Return {chain: [(res_seq, res_name, x, y, z), ...]} for heavy atoms."""
    wanted = set(chains)
    atoms = defaultdict(list)
    with open_structure(path) as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break  # first model only
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21].strip()
            if chain not in wanted:
                continue
            res_name = line[17:20].strip()
            if res_name in SKIP_RESIDUES:
                continue
            element = line[76:78].strip() or line[12:16].strip()[:1]
            if element.upper() == "H":
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                res_seq = int(line[22:26])
            except ValueError:
                continue
            atoms[chain].append((res_seq, res_name, x, y, z))
    return atoms


def contacts(target_atoms, binder_atoms, cutoff=4.5):
    """Target residues within cutoff of any binder atom."""
    cutoff_sq = cutoff * cutoff
    # Bin binder atoms into a grid so this stays fast on large complexes
    cell = cutoff
    grid = defaultdict(list)
    for _, _, x, y, z in binder_atoms:
        grid[(int(x // cell), int(y // cell), int(z // cell))].append((x, y, z))

    found = {}
    for res_seq, res_name, x, y, z in target_atoms:
        gx, gy, gz = int(x // cell), int(y // cell), int(z // cell)
        best_sq, hits = None, 0
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for bx, by, bz in grid.get((gx + dx, gy + dy, gz + dz), ()):
                        d_sq = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
                        if d_sq <= cutoff_sq:
                            hits += 1
                            if best_sq is None or d_sq < best_sq:
                                best_sq = d_sq
        if hits:
            prev = found.get(res_seq)
            dist = math.sqrt(best_sq)
            if prev is None:
                found[res_seq] = [res_name, dist, hits]
            else:
                prev[1] = min(prev[1], dist)
                prev[2] += hits
    return found


def find_structure_file(structure_dir, pdb_id):
    for name in (f"{pdb_id}.pdb", f"{pdb_id}.pdb.gz"):
        path = Path(structure_dir) / name
        if path.exists():
            return path
    return None


def read_pairs(table_path):
    rows, header = [], None
    with open(table_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    return rows


def main():
    if "snakemake" in globals():
        binder_table = snakemake.input.binder_table
        structure_dir = snakemake.params.structure_dir
        out_path = snakemake.output.contact_table
        cutoff = float(snakemake.params.get("contact_cutoff", 4.5))
    else:
        binder_table = sys.argv[1] if len(sys.argv) > 1 else "results/phase2/binder_structures.tsv"
        structure_dir = sys.argv[2] if len(sys.argv) > 2 else "data/protein_structures/pdb_files"
        out_path = sys.argv[3] if len(sys.argv) > 3 else "results/phase2/binder_contacts.tsv"
        cutoff = float(sys.argv[4]) if len(sys.argv) > 4 else 4.5

    pairs = read_pairs(binder_table)
    print(f"Measuring contacts for {len(pairs)} target/binder pairs at {cutoff} A")

    rows, no_contact, missing = [], 0, 0
    for pair in pairs:
        path = find_structure_file(structure_dir, pair["pdb_id"])
        if path is None:
            missing += 1
            continue

        target_chains = [c for c in pair["target_chain"].split(",") if c]
        binder_chains = [c for c in pair["binder_chain"].split(",") if c]
        atoms = load_atoms(path, target_chains + binder_chains)

        binder_atoms = [a for c in binder_chains for a in atoms.get(c, [])]
        if not binder_atoms:
            continue

        pair_rows = 0
        for tchain in target_chains:
            target_atoms = atoms.get(tchain, [])
            if not target_atoms:
                continue
            for res_seq, (res_name, dist, hits) in sorted(
                    contacts(target_atoms, binder_atoms, cutoff).items()):
                rows.append({
                    "pdb_id": pair["pdb_id"],
                    "target_chain": tchain,
                    "binder_chain": pair["binder_chain"],
                    "binder_molecule": pair["binder_molecule"],
                    "binder_type": pair["binder_type"],
                    "res_seq": res_seq,
                    "res_name": res_name,
                    "min_distance": f"{dist:.2f}",
                    "n_atom_contacts": hits,
                })
                pair_rows += 1
        if pair_rows == 0:
            no_contact += 1

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")

    footprints = {(r["pdb_id"], r["binder_chain"], r["target_chain"]) for r in rows}
    print(f"Contacted residues : {len(rows)}")
    print(f"Distinct footprints: {len(footprints)}")
    if no_contact:
        print(f"Pairs with no contact within {cutoff} A: {no_contact}")
    if missing:
        print(f"Structure file not found for {missing} pairs")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
