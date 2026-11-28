#!/usr/bin/env python3
"""
Phase 2 report: which known binders exist, and would they work broadly?

Answers three questions per protein:
  - was its structure solved with something bound, and what?
  - how conserved are the residues that binder grips, across our species?
  - do any of phase 1's predicted epitopes sit in that same patch?

Uses the same page shell as the phase 1 report so both read alike.
"""

import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

try:
    import report_template as rt
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "protein_analysis"))
    import report_template as rt


BINDER_LABEL = {
    "nanobody": "Nanobody / VHH",
    "antibody": "Antibody fragment",
    "peptide_or_small_molecule": "Peptide / small molecule",
}


def find_figure(figures_dir, gene, pdb_id, binder_molecule, target_chain):
    """Locate the render for one footprint, matching how it was named."""
    if not figures_dir or not os.path.isdir(figures_dir):
        return None
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_"
                   for ch in binder_molecule)[:28]
    exact = Path(figures_dir) / f"{gene}_{pdb_id}_{safe}_{target_chain}.png"
    if exact.exists():
        return exact
    # fall back to any render for this gene/structure
    for candidate in sorted(Path(figures_dir).glob(f"{gene}_{pdb_id}_*.png")):
        return candidate
    return None


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, errors="replace") as handle:
        return [r for r in csv.DictReader(handle, delimiter="\t")]


def num(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def verdict_badge(pct_conserved, broad, partial):
    """Colour the verdict by how broadly the footprint holds."""
    if pct_conserved >= broad:
        return f'<span class="badge badge-good badge-filled">{rt.icon("check-double")}Broad</span>'
    if pct_conserved >= partial:
        return f'<span class="badge badge-warning">{rt.icon("dash")}Partial</span>'
    return f'<span class="badge badge-critical">{rt.icon("cross")}Narrow</span>'


def build(summary_rows, conservation_rows, overlap_rows, binder_rows, cfg, report_cfg,
          analysis, paramset, figures_dir=None):
    broad = float(cfg.get("broad_threshold", 70))
    partial = float(cfg.get("partial_threshold", 40))
    include_sm = bool(report_cfg.get("include_small_molecules", True))
    max_per_gene = int(report_cfg.get("max_binders_per_gene", 12))

    if not include_sm:
        summary_rows = [r for r in summary_rows
                        if r["binder_type"] != "peptide_or_small_molecule"]

    sections, nav_main = [], [("overview", "Overview", None)]

    # ---- headline numbers -------------------------------------------------
    structures = {r["pdb_id"] for r in binder_rows}
    genes = sorted({r["gene"] for r in summary_rows if r["gene"]})
    by_type = defaultdict(list)
    for row in summary_rows:
        if row["pct_highly_conserved"]:
            by_type[row["binder_type"]].append(num(row["pct_highly_conserved"]))

    overlaps = [r for r in overlap_rows if int(r.get("overlap_residues") or 0) > 0]
    epitopes_hit = {(r["gene"], r["epitope_start"], r["epitope_end"]) for r in overlaps}

    tiles = rt.stat_row([
        rt.stat_tile(len(structures), "Structures with a binder"),
        rt.stat_tile(len(summary_rows), "Binder footprints"),
        rt.stat_tile(len(genes), "Proteins covered"),
        rt.stat_tile(len(epitopes_hit), "Epitopes at a binder site"),
    ])

    body = tiles
    if not summary_rows:
        body += rt.note(
            "warning", "No binder structures found",
            "<p>None of the downloaded structures were solved with a nanobody, "
            "antibody or peptide bound, so there is nothing to score. This is a "
            "property of what is in the PDB for these proteins, not a failure.</p>")
    else:
        rows = []
        for btype, values in sorted(by_type.items(), key=lambda kv: -len(kv[1])):
            mean = sum(values) / len(values)
            rows.append([
                rt.esc(BINDER_LABEL.get(btype, btype)),
                rt.esc(len(values)),
                rt.score_num(mean / 100),
                verdict_badge(mean, broad, partial),
            ])
        body += (
            "<h3>What binds, and how conserved its grip is</h3>"
            "<p>Each binder's footprint is the set of residues it touches in the "
            "structure. The share of those residues that are highly conserved "
            "decides whether the binder could work across species, or only in the "
            "strain it came from.</p>"
            + rt.table(["Binder type", "Footprints", "Mean conserved fraction", "Typical verdict"],
                       rows, table_id="by-type", numeric_from=1))

    sections.append(rt.section(
        "overview", "Known binders and their footprints", body,
        subtitle=f"{analysis}_{paramset} · conserved = ConSurf grade "
                 f"{int(cfg.get('conserved_grade', 7))}+"))

    # ---- all footprints ---------------------------------------------------
    if summary_rows:
        rows = []
        for row in sorted(summary_rows, key=lambda r: (-num(r["pct_highly_conserved"]),
                                                       r["gene"])):
            pct = num(row["pct_highly_conserved"])
            rows.append([
                rt.esc(row["gene"]),
                rt.esc(row["pdb_id"]),
                rt.esc(row["binder_molecule"][:34]),
                rt.esc(BINDER_LABEL.get(row["binder_type"], row["binder_type"])),
                rt.esc(row["contact_residues"]),
                rt.score_meter(pct / 100),
                rt.esc(row["mean_species_identity_pct"] or "-"),
                verdict_badge(pct, broad, partial),
            ])
        sections.append(rt.section(
            "footprints", "Every binder footprint",
            '<p class="muted">Sorted by how much of the footprint is conserved. '
            'Species identity is the mean share of aligned species carrying the '
            'same residue at those positions.</p>'
            + rt.table(["Protein", "PDB", "Binder", "Type", "Contacts",
                        "Conserved fraction", "Species ID %", "Verdict"],
                       rows, table_id="footprints", numeric_from=4)))
        nav_main.append(("footprints", "All footprints", len(summary_rows)))

    # ---- epitope overlap --------------------------------------------------
    if overlaps:
        rows = []
        for row in sorted(overlaps, key=lambda r: -int(r["overlap_residues"])):
            rows.append([
                rt.esc(row["gene"]),
                f'<span class="seq">{rt.esc(row["peptide"])}</span>',
                rt.esc(f'{row["epitope_start"]}-{row["epitope_end"]}'),
                rt.esc(row["binder_molecule"][:30]),
                rt.esc(BINDER_LABEL.get(row["binder_type"], row["binder_type"])),
                rt.esc(f'{row["overlap_residues"]}/{row["epitope_length"]}'),
                rt.score_meter(num(row["overlap_pct"]) / 100),
                rt.esc(row["mean_consurf_grade_of_overlap"] or "-"),
                rt.esc(row["mean_species_identity_of_overlap"] or "-"),
            ])
        body = (
            "<p>These predicted epitopes sit where a binder has actually been "
            "observed to bind. Phase 1 found them from sequence and topology "
            "alone, so an overlap here is independent experimental support that "
            "the site is exposed and bindable.</p>"
            + rt.table(["Protein", "Epitope", "Position", "Binder", "Type",
                        "Overlap", "Overlap %", "Mean grade", "Species ID %"],
                       rows, table_id="overlap", numeric_from=5))
        sections.append(rt.section("overlap", "Predicted epitopes at a known binder site", body))
        nav_main.append(("overlap", "Epitope overlap", len(overlaps)))
    else:
        no_structure = sum(1 for r in overlap_rows
                           if "no binder structure" in (r.get("note") or ""))
        sections.append(rt.section(
            "overlap", "Predicted epitopes at a known binder site",
            rt.note("info", "No overlaps found",
                    f"<p>None of the predicted epitopes fall inside a binder "
                    f"footprint. {no_structure} of them belong to proteins with no "
                    f"binder structure at all, so there was nothing to compare "
                    f"against.</p>")))
        nav_main.append(("overlap", "Epitope overlap", None))

    # ---- per protein ------------------------------------------------------
    def figure_html(row):
        path = find_figure(figures_dir, row["gene"], row["pdb_id"],
                           row["binder_molecule"], row["target_chain"])
        if path is None:
            return ""
        return rt.structure_figure(
            path, [],
            caption=(f'{row["binder_molecule"]} on {row["gene"]} ({row["pdb_id"]}). '
                     f'Grey is the protein, the binder is outlined; contact residues '
                     f'are coloured by conservation - blue conserved, red variable.'))

    nav_genes = []
    by_gene = defaultdict(list)
    for row in summary_rows:
        by_gene[row["gene"]].append(row)

    contacts_by_binder = defaultdict(list)
    for row in conservation_rows:
        contacts_by_binder[(row["gene"], row["pdb_id"], row["binder_molecule"])].append(row)

    for gene in sorted(by_gene):
        anchor = "g-" + gene
        gene_rows = sorted(by_gene[gene],
                           key=lambda r: -num(r["pct_highly_conserved"]))[:max_per_gene]
        best = gene_rows[0] if gene_rows else None
        chips = (f'<div class="chiprow">'
                 f'<span class="chip">{len(by_gene[gene])} footprints</span>'
                 f'<span class="chip">best conserved '
                 f'{best["pct_highly_conserved"] if best else "-"}%</span></div>')

        blocks = [chips]
        for row in gene_rows:
            key = (row["gene"], row["pdb_id"], row["binder_molecule"])
            contacts = sorted(contacts_by_binder.get(key, []),
                              key=lambda r: int(r["consurf_position"] or 0))
            if not contacts:
                continue
            positions = ", ".join(
                f'{c["res_name"]}{c["consurf_position"]}'
                f'<span class="muted">(g{c["consurf_grade"]})</span>'
                for c in contacts[:40])
            pct = num(row["pct_highly_conserved"])
            blocks.append(
                f'<h3>{rt.esc(row["binder_molecule"])} '
                f'<span class="muted">· {rt.esc(row["pdb_id"])} · '
                f'{rt.esc(BINDER_LABEL.get(row["binder_type"], row["binder_type"]))}</span></h3>'
                f'<div class="chiprow">'
                f'<span class="chip">{len(contacts)} contacts</span>'
                f'<span class="chip">{row["pct_highly_conserved"]}% conserved</span>'
                f'<span class="chip">species ID {row["mean_species_identity_pct"] or "-"}%</span>'
                f'</div>{verdict_badge(pct, broad, partial)}'
                f'<p class="muted" style="margin-top:8px">Footprint (residue, position, '
                f'ConSurf grade): {positions}'
                f'{" …" if len(contacts) > 40 else ""}</p>'
                + figure_html(row))

        sections.append(rt.section(anchor, gene, "".join(blocks)))
        nav_genes.append((anchor, gene, len(by_gene[gene])))

    # ---- methods ----------------------------------------------------------
    methods = f"""
    <h3>How a binder is found</h3>
    <p>A structure counts as having a binder when one of its chains comes from an
    organism that is not the target bacterium &mdash; a llama or alpaca for a
    nanobody, human or mouse for a Fab, or a synthetic construct for a designed
    peptide. Chains are joined to organisms through <code>MOL_ID</code> in the PDB
    header. Matching on chain names does not work: the BamA nanobody chain is
    simply called <code>NANOF7</code>.</p>

    <h3>Footprint</h3>
    <p>Every target residue with a heavy atom within
    {rt.esc(cfg.get('contact_cutoff', 4.5))}&nbsp;&#8491; of a binder atom. Hydrogens,
    waters and ions are ignored.</p>

    <h3>Conservation</h3>
    <p>Each contact residue is matched to its ConSurf grade through the grades
    file, which ties PDB residues to alignment positions. Grades run 1&ndash;9,
    with {int(cfg.get('conserved_grade', 7))}+ counted as conserved. Species identity
    is the share of aligned species carrying the same residue at that position.</p>

    <h3>Verdicts</h3>
    <ul>
      <li>{verdict_badge(100, broad, partial)} at least {broad:.0f}% of contacts conserved</li>
      <li>{verdict_badge(50, broad, partial)} at least {partial:.0f}%</li>
      <li>{verdict_badge(0, broad, partial)} below {partial:.0f}% &mdash; likely strain-specific</li>
    </ul>
    <p class="muted">These describe how conserved a binding site is, not whether a
    binder works. A conserved footprint is necessary for breadth, not sufficient.</p>
    """
    sections.append(rt.section("methods", "Methods", methods))
    nav_main.append(("methods", "Methods", None))

    nav_groups = [("Report", nav_main)]
    if nav_genes:
        nav_groups.append((f"Proteins ({len(nav_genes)})", nav_genes))

    return rt.render_page(
        title="Known binders | conservation of their footprints",
        subtitle="Phase 2 · binder analysis",
        meta_items=[("Analysis", f"{analysis}_{paramset}"),
                    ("Structures", len(structures)),
                    ("Footprints", len(summary_rows))],
        nav_groups=nav_groups,
        sections_html="\n".join(sections),
    )


def main():
    if "snakemake" in globals():
        summary_path = snakemake.input.summary_table
        conservation_path = snakemake.input.conservation_table
        overlap_path = snakemake.input.overlap_table
        binder_path = snakemake.input.binder_table
        out_path = snakemake.output.report
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        cfg = dict(snakemake.params.conservation or {})
        cfg.setdefault("contact_cutoff", 4.5)
        report_cfg = dict(snakemake.params.report_config or {})
        figures_dir = str(Path(snakemake.output.report).parent / "figures")
    else:
        args = sys.argv[1:]
        base = args[0] if args else "results/phase2"
        summary_path = f"{base}/binder_summary.tsv"
        conservation_path = f"{base}/binder_conservation.tsv"
        overlap_path = f"{base}/epitope_binder_overlap.tsv"
        binder_path = f"{base}/binder_structures.tsv"
        out_path = args[1] if len(args) > 1 else f"{base}/phase2_report.html"
        analysis, paramset = "analysis2", "params2"
        cfg = {"conserved_grade": 7, "broad_threshold": 70,
               "partial_threshold": 40, "contact_cutoff": 4.5}
        report_cfg = {"max_binders_per_gene": 12, "include_small_molecules": True}
        figures_dir = f"{base}/figures"

    html = build(read_tsv(summary_path), read_tsv(conservation_path),
                 read_tsv(overlap_path), read_tsv(binder_path),
                 cfg, report_cfg, analysis, paramset, figures_dir)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)
    print(f"Wrote {out_path} ({out_path.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
