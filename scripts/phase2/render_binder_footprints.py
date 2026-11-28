#!/usr/bin/env python3
"""
Render each binder's footprint on the target structure, coloured by conservation.

A table of contact residues is hard to picture. This draws the protein in grey
and colours the residues the binder touches by how conserved they are, so the
difference between a binder gripping a conserved core and one gripping variable
loops is visible at a glance.

Colours follow the conservation, not the binder:
    dark blue   ConSurf grade 9-8   highly conserved
    light blue  grade 7-6
    white       grade 5-4           average
    salmon      grade 3-2
    red         grade 1             highly variable

That is ConSurf's own convention reversed onto a blue-white-red scale, so a
mostly-blue footprint means a binder that could work broadly.

One PNG per binder footprint, plus a .pse session to open in PyMOL.
"""

import gzip
import os
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# ConSurf grade -> colour. 9 is most conserved.
GRADE_COLOUR = {
    9: "density", 8: "slate", 7: "lightblue", 6: "palecyan",
    5: "white", 4: "white",
    3: "wheat", 2: "salmon", 1: "firebrick",
}


def read_tsv(path):
    import csv
    with open(path, errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def decompress(path, workdir):
    """PyMOL will not read .gz directly here; stage a plain copy."""
    path = Path(path)
    if path.suffix != ".gz":
        return path
    out = Path(workdir) / path.stem
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def find_structure(structure_dir, pdb_id):
    for name in (f"{pdb_id}.pdb", f"{pdb_id}.pdb.gz"):
        candidate = Path(structure_dir) / name
        if candidate.exists():
            return candidate
    return None


def render(cmd, structure_path, pdb_id, target_chain, binder_chains, contacts,
           out_png, out_pse, width, height, ray_trace):
    """Draw one footprint. contacts: {res_seq: grade}."""
    cmd.delete("all")
    cmd.reinitialize()
    cmd.load(str(structure_path), "complex")

    target_sel = f"complex and chain {target_chain} and polymer"
    cmd.hide("everything", "all")
    cmd.show("cartoon", target_sel)
    cmd.color("grey80", target_sel)

    # the binder itself, as a thin outline for context
    if binder_chains:
        binder_sel = "complex and chain " + "+".join(binder_chains)
        cmd.show("ribbon", binder_sel)
        cmd.color("grey50", binder_sel)

    # contacts coloured by conservation
    by_colour = defaultdict(list)
    for res_seq, grade in contacts.items():
        by_colour[GRADE_COLOUR.get(grade, "white")].append(str(res_seq))
    for colour, residues in by_colour.items():
        sel = f"{target_sel} and resi {'+'.join(residues)}"
        cmd.show("sticks", sel)
        cmd.color(colour, sel)

    cmd.bg_color("white")
    cmd.orient(target_sel)
    if contacts:
        cmd.zoom(f"{target_sel} and resi "
                 f"{'+'.join(str(r) for r in contacts)}", buffer=8)
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0.25)

    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    cmd.png(str(out_png), width=width, height=height, dpi=150, ray=1 if ray_trace else 0)
    cmd.save(str(out_pse))


def main():
    if "snakemake" in globals():
        conservation_path = snakemake.input.conservation_table
        structure_dir = snakemake.params.structure_dir
        out_dir = snakemake.params.out_dir
        sentinel = snakemake.output.sentinel
        cfg = dict(snakemake.params.get("figure", {}) or {})
    else:
        args = sys.argv[1:]
        conservation_path = args[0] if args else "results/analysis2_params2/phase2/binder_conservation.tsv"
        structure_dir = args[1] if len(args) > 1 else "data/protein_structures/pdb_files"
        out_dir = args[2] if len(args) > 2 else "results/analysis2_params2/phase2/figures"
        sentinel = args[3] if len(args) > 3 else None
        cfg = {}

    width = int(cfg.get("width", 1100))
    height = int(cfg.get("height", 850))
    ray_trace = bool(cfg.get("ray_trace", True))

    rows = read_tsv(conservation_path)
    binder_rows = read_tsv(Path(conservation_path).with_name("binder_structures.tsv"))
    binder_chains_for = {}
    for row in binder_rows:
        binder_chains_for[(row["pdb_id"], row["binder_molecule"])] = [
            c for c in row["binder_chain"].split(",") if c]

    footprints = defaultdict(dict)
    meta = {}
    for row in rows:
        key = (row["gene"], row["pdb_id"], row["binder_molecule"], row["target_chain"])
        try:
            footprints[key][int(row["res_seq"])] = int(row["consurf_grade"] or 0)
        except ValueError:
            continue
        meta[key] = row["binder_type"]

    print(f"Footprints to render: {len(footprints)}")
    if not footprints:
        print("Nothing to draw")
    else:
        import pymol
        pymol.finish_launching(["pymol", "-cq"])
        from pymol import cmd

        drawn, skipped = 0, 0
        with tempfile.TemporaryDirectory() as workdir:
            for (gene, pdb_id, binder, tchain), contacts in sorted(footprints.items()):
                path = find_structure(structure_dir, pdb_id)
                if path is None:
                    skipped += 1
                    continue
                safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in binder)[:28]
                stem = f"{gene}_{pdb_id}_{safe}_{tchain}"
                try:
                    render(cmd, decompress(path, workdir), pdb_id, tchain,
                           binder_chains_for.get((pdb_id, binder), []),
                           contacts,
                           Path(out_dir) / f"{stem}.png",
                           Path(out_dir) / f"{stem}.pse",
                           width, height, ray_trace)
                    drawn += 1
                    print(f"  {stem}: {len(contacts)} contacts")
                except Exception as e:
                    skipped += 1
                    print(f"  FAILED {stem}: {type(e).__name__}: {e}")

        print(f"Rendered {drawn} footprints, skipped {skipped}")

    if sentinel:
        Path(sentinel).parent.mkdir(parents=True, exist_ok=True)
        with open(sentinel, "w") as handle:
            handle.write(f"Rendered footprints: {len(footprints)}\n")


if __name__ == "__main__":
    main()
