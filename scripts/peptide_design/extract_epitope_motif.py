#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_epitope_motif.py

Given epitope predictions (TSV) that are grouped by protein/PDB headers like "((5OR1))",
load the corresponding PDB from a directory (e.g., data/protein_structures/pdb_files/5OR1.pdb.gz),
analyze the epitope surface, and propose complementary peptide binding motifs.

Outputs:
- CSV summary with one row per epitope region (hotspot)
- JSON report per PDB with detailed feature stats and motif rationale

Dependencies:
- Python 3.8+
- biopython
- pandas
- Optional: freesasa (for ASA)

This version uses a precomputed FASTA→PDB mapping table first, and falls back to
simple heuristics only if needed.
"""

import os
import re
import io
import json
import gzip
import math
from collections import Counter
import csv

import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, is_aa
from Bio.SeqUtils import seq1 as three_to_one

# Optional: freesasa for solvent-accessible surface area
try:
    import freesasa
    HAS_FREESASA = True
except Exception:
    HAS_FREESASA = False

# ---------- Residue property tables ----------
AA_CLASSES = {
    'positive': set('KRH'),
    'negative': set('DE'),
    'aromatic': set('FWY'),
    'hydrophobic': set('AILMVFWY'),
    'polar': set('STNQCY'),
    'special': set('PGC')  # proline/glycine/cysteine
}

COMPLEMENT_RULES = {
    'negative': '[KR]',          # complement negative with positive
    'positive': '[DE]',          # complement positive with negative
    'aromatic': '[FWY]',         # π–π / aromatic contacts
    'hydrophobic': '[LIVAMFWY]',
    'polar': '[STNQHY]',
}

# ---------- Utility ----------
def read_pdb_any(path):
    """Read PDB (possibly .gz or .cif) and return Biopython structure."""
    if not os.path.exists(path):
        if os.path.exists(path + ".gz"):
            path = path + ".gz"
        elif os.path.exists(path.replace(".pdb", ".cif")):
            path = path.replace(".pdb", ".cif")
        elif os.path.exists(path.replace(".pdb", ".cif.gz")):
            path = path.replace(".pdb", ".cif.gz")
        else:
            raise FileNotFoundError(f"PDB not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext == ".gz":
        with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as fh:
            head = fh.read(1000)
            fmt = 'pdb' if ('ATOM' in head or 'HETATM' in head) else ('cif' if 'data_' in head else 'pdb')
        with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as handle:
            data = handle.read()
        bio_handle = io.StringIO(data)
        if fmt == 'pdb':
            parser = PDBParser(QUIET=True)
            return parser.get_structure(os.path.basename(path), bio_handle)
        else:
            parser = MMCIFParser(QUIET=True)
            bio_handle.seek(0)
            return parser.get_structure(os.path.basename(path), bio_handle)
    elif ext == ".pdb":
        parser = PDBParser(QUIET=True)
        return parser.get_structure(os.path.basename(path), path)
    elif ext == ".cif":
        parser = MMCIFParser(QUIET=True)
        return parser.get_structure(os.path.basename(path), path)
    else:
        parser = PDBParser(QUIET=True)
        return parser.get_structure(os.path.basename(path), path)

def parse_epitope_blocks(epi_path):
    """
    Parse a flexible epitope TSV/text file.
    Supports headers like "((5OR1))" followed by lines that may contain columns or simple ranges.
    Returns a dict: pdb_id -> DataFrame with columns [chain, start, end, sequence]
    """
    with open(epi_path, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()

    blocks = re.split(r'\(\(([^)]+)\)\)', text)  # ['', '5OR1', '<content>', ...]
    epi_map = {}
    if len(blocks) <= 1:
        # try as TSV/CSV
        try:
            df = pd.read_csv(io.StringIO(text), sep=None, engine='python')
            if 'pdb' in df.columns or 'PDB' in df.columns:
                key = 'pdb' if 'pdb' in df.columns else 'PDB'
                for pdb_id, sub in df.groupby(key):
                    epi_map[str(pdb_id)] = normalize_epitope_df(sub)
            else:
                pdb_id_guess = guess_pdb_from_filename(epi_path)
                epi_map[pdb_id_guess] = normalize_epitope_df(df)
            return epi_map
        except Exception:
            pass

    i = 1
    while i < len(blocks):
        pdb_id = blocks[i].strip()
        content = blocks[i+1] if i+1 < len(blocks) else ""
        i += 2
        if not pdb_id:
            continue
        df = parse_block_to_df(content)
        epi_map[pdb_id] = normalize_epitope_df(df)
    return epi_map

def guess_pdb_from_filename(path):
    base = os.path.basename(path)
    m = re.search(r'([0-9][A-Za-z0-9]{3})', base)
    return m.group(1).upper() if m else os.path.splitext(base)[0]

def parse_block_to_df(block_text):
    block_text = block_text.strip()
    if not block_text:
        return pd.DataFrame(columns=['chain','start','end','sequence'])
    try:
        df = pd.read_csv(io.StringIO(block_text), sep=None, engine='python')
        return df
    except Exception:
        pass

    rows = []
    for line in block_text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        kv = dict(re.findall(r'(\w+)\s*=\s*([A-Za-z0-9_\-]+)', line))
        chain = kv.get('chain') or kv.get('Chain') or kv.get('CHAIN')
        start = kv.get('start') or kv.get('Start') or kv.get('begin')
        end   = kv.get('end')   or kv.get('End')   or kv.get('stop')
        seq   = kv.get('sequence') or kv.get('SEQ') or kv.get('Sequence')

        if not (start and end):
            parts = [p.strip() for p in re.split(r'[,\t]+', line) if p.strip()]
            if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
                chain, start, end = parts[0], parts[1], parts[2]
            elif len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                start, end = parts[0], parts[1]
            else:
                m = re.match(r'([A-Za-z])?[:\s\-]*?(\d+)\s*[-:]\s*(\d+)', line)
                if m:
                    chain = m.group(1) or None
                    start, end = m.group(2), m.group(3)

        if start and end:
            rows.append({'chain': chain, 'start': int(start), 'end': int(end), 'sequence': seq})
        elif seq:
            rows.append({'chain': chain, 'start': None, 'end': None, 'sequence': seq})

    return pd.DataFrame(rows)

def normalize_epitope_df(df):
    cols = ['chain','start','end','sequence']
    for c in cols:
        if c not in df.columns:
            df[c] = None
    for c in ['start','end']:
        df[c] = pd.to_numeric(df[c], errors='coerce').astype('Int64')
    if 'chain' in df.columns:
        df['chain'] = df['chain'].astype(str).str.upper().replace({'<NA>': None, 'NAN': None})
        df.loc[df['chain'].isin(['NONE','NULL','-','N/A']), 'chain'] = None
    if 'sequence' in df.columns:
        df['sequence'] = df['sequence'].astype(str).str.strip().str.upper().replace({'<NA>': None, 'NAN': None})
        df.loc[df['sequence'].isin(['NONE','NULL','-','N/A','']), 'sequence'] = None
    return df[cols]

def chain_seq_and_map(chain):
    """Return the chain one-letter sequence and a list [(residue_obj, pdb_resnum), ...]."""
    ppb = PPBuilder()
    seq = ""
    idx_map = []
    for pp in ppb.build_peptides(chain):
        for r in pp:
            if not is_aa(r, standard=True):
                continue
            seq += three_to_one(r.get_resname())
            idx_map.append((r, r.id[1]))
    return seq, idx_map

def map_epitope_to_chain(chain_map, start, end, seq):
    """
    Fallback mapper (no alignment): PDB numbering if given, else exact substring.
    Returns a list of residue objects (or []).
    """
    idx_map = chain_map['idx']
    chain_seq = chain_map['seq']

    if start is not None and end is not None:
        residues = [r for (r, resnum) in idx_map if start <= resnum <= end]
        if residues:
            return residues

    if seq:
        pos = chain_seq.find(seq)
        if pos != -1:
            window = idx_map[pos:pos+len(seq)]
            return [r for (r, _) in window]

    return []

def get_ca_coord(residue):
    try:
        return residue['CA'].get_coord()
    except Exception:
        for atom in residue.get_atoms():
            return atom.get_coord()
    return None

def distance(a, b):
    return math.sqrt(((a[0]-b[0])**2)+((a[1]-b[1])**2)+((a[2]-b[2])**2))

def simple_cluster(points, radius=8.0, min_pts=3):
    """Very simple density clustering."""
    if not points:
        return []
    coords = [p[1] for p in points]
    unused = set(range(len(points)))
    clusters = []
    while unused:
        seed = unused.pop()
        cluster = {seed}
        changed = True
        while changed:
            changed = False
            for i in list(unused):
                if any(distance(coords[i], coords[j]) <= radius for j in cluster):
                    cluster.add(i)
                    unused.remove(i)
                    changed = True
        if len(cluster) >= min_pts:
            clusters.append(sorted(list(cluster)))
    return clusters if clusters else [[i] for i in range(len(points))]

def epitope_features(residues, asa_by_res=None):
    comp = Counter()
    for r in residues:
        if not is_aa(r, standard=True):
            continue
        aa = three_to_one(r.get_resname())
        for k, s in AA_CLASSES.items():
            if aa in s:
                comp[k] += 1
    n = sum(comp.values()) or 1
    frac = {k: comp[k]/n for k in AA_CLASSES.keys()}
    return comp, frac

def propose_motif_from_features(frac, length_hint=9):
    ranked = sorted(frac.items(), key=lambda x: x[1], reverse=True)
    top = [k for k, v in ranked if v > 0.15] or [ranked[0][0]]
    L = max(7, min(12, length_hint))
    motif = ['X'] * L

    positions = [0, 2, 4, 6, 8, 10]
    pos_i = 0
    for t in top[:3]:
        comp = COMPLEMENT_RULES.get(t, 'X')
        if pos_i < len(positions):
            motif[positions[pos_i]] = comp
            pos_i += 1

    for i in range(1, L-1, 3):
        if motif[i] == 'X':
            motif[i] = '{P}'  # avoid Pro in helix-like core

    if 'aromatic' in top and any(m == 'X' for m in motif):
        for i in range(L):
            if motif[i] == 'X':
                motif[i] = '[LIVAMFWY]'
                break

    motif = [m if m != 'X' else '[ACSTNQHRK]' for m in motif]
    return '-'.join(motif)

# === Mapping table config ===
MAP_PDB_COL = "structure_id"
MAP_CHAIN_COL = "chain"
MAP_FASTA_START_COL = "chain_start"   # start in FASTA numbering
MAP_FASTA_END_COL = "chain_end"       # end in FASTA numbering

def load_seq_to_pdb_mapping(mapping_tsv):
    """
    Load a mapping from a file with start/end ranges per chain.
    Returns mapping[pdb] -> list of (chain_id, fasta_start, fasta_end).
    """
    mapping = {}
    with open(mapping_tsv, newline='', encoding='utf-8', errors='ignore') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            pdb = str(row[MAP_PDB_COL]).strip().upper()
            chain = str(row[MAP_CHAIN_COL]).strip()
            try:
                start = int(row[MAP_FASTA_START_COL])
                end = int(row[MAP_FASTA_END_COL])
            except Exception:
                continue
            mapping.setdefault(pdb, []).append((chain, start, end))
    return mapping

def residues_from_mapping(pdb_id, chain_maps, mapping, start, end):
    """
    Map a FASTA start–end epitope range to PDB residues using a range mapping table.
    """
    if start is None or end is None:
        return []
    if pdb_id.upper() not in mapping:
        return []

    hits = []
    for chain_id, fasta_start, fasta_end in mapping[pdb_id.upper()]:
        # check overlap
        if start >= fasta_start and end <= fasta_end:
            if chain_id in chain_maps:
                cmap = chain_maps[chain_id]
                residues = [r for (r, pdb_resnum) in cmap['idx']
                            if fasta_start <= start <= fasta_end and
                               fasta_start <= end <= fasta_end and
                               pdb_resnum >= start and pdb_resnum <= end]
                hits.extend((chain_id, r) for r in residues)
    return hits


# ---------- Main per-PDB analysis ----------
def analyze_pdb_epitopes(pdb_id, epidf, pdb_dir, out_dir):
    pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
    if not os.path.exists(pdb_path) and os.path.exists(os.path.join(pdb_dir, f"{pdb_id}.pdb.gz")):
        pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb.gz")
    structure = read_pdb_any(pdb_path)

    # Build chain maps
    chain_maps = {}
    for model in structure:
        for chain in model:
            seq, idx_map = chain_seq_and_map(chain)
            chain_maps[chain.id] = {'chain': chain, 'seq': seq, 'idx': idx_map}

    # Optional ASA init (not used in scoring here; safe to keep)
    if HAS_FREESASA:
        try:
            if pdb_path.endswith('.gz'):
                with gzip.open(pdb_path, 'rt', encoding='utf-8', errors='ignore') as fh:
                    pdb_str = fh.read()
            else:
                with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as fh:
                    pdb_str = fh.read()
            _ = freesasa.calc(freesasa.Structure(pdb_str))
        except Exception:
            pass

    # ---- Gather epitope residues (mapping-first, then fallback) ----
    epitope_residues = []
    mapping_log = []
    entry_idx = 0

    for _, row in epidf.iterrows():
        entry_idx += 1
        chain_hint = row['chain'] if pd.notna(row['chain']) else None
        start = int(row['start']) if pd.notna(row['start']) else None
        end   = int(row['end'])   if pd.notna(row['end'])   else None
        seq   = row['sequence'] if isinstance(row['sequence'], str) and row['sequence'] not in (None, 'None', '') else None

        mapped = residues_from_mapping(pdb_id, chain_maps, SEQ2PDB, start, end)
        status = "mapped_mapping" if mapped else "unmapped_mapping"

        if not mapped:
            chains_to_scan = [chain_hint] if chain_hint and chain_hint in chain_maps else list(chain_maps.keys())
            for cid in chains_to_scan:
                cmap = chain_maps[cid]
                fallback = map_epitope_to_chain(cmap, start, end, seq)
                if fallback:
                    mapped = [(cid, r) for r in fallback]
                    status = "mapped_fallback"
                    break

        if mapped:
            epitope_residues.extend(mapped)
            mapping_log.append({
                'entry': entry_idx, 'pdb': pdb_id,
                'mapped_chains': ''.join(sorted({c for c, _ in mapped})),
                'n_residues': len(mapped),
                'start_in': start, 'end_in': end, 'seq_in': seq,
                'status': status
            })
        else:
            mapping_log.append({
                'entry': entry_idx, 'pdb': pdb_id,
                'mapped_chains': None,
                'n_residues': 0,
                'start_in': start, 'end_in': end, 'seq_in': seq,
                'status': 'unmapped_all'
            })

    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame(mapping_log).to_csv(
        os.path.join(out_dir, f"{pdb_id}_mapping_log.tsv"),
        sep='\t', index=False
    )

    # ---- Cluster epitopes by 3D proximity (hotspots) ----
    points = []
    for i, (cid, res) in enumerate(epitope_residues):
        coord = get_ca_coord(res)
        if coord is not None:
            points.append((i, coord))

    clusters_idx = []
    if points:
        if len(points) < 3:
            clusters_idx = [list(range(len(points)))]
        else:
            clusters_idx = simple_cluster(points, radius=8.0, min_pts=3)
            if not clusters_idx:
                clusters_idx = [list(range(len(points)))]

    # ---- Summaries & motifs ----
    clusters = []
    for cl in clusters_idx:
        members = [epitope_residues[i][1] for i in cl]
        comp, frac = epitope_features(members, asa_by_res=None)
        coords = [get_ca_coord(r) for r in members]
        max_span = 0.0
        for a in coords:
            for b in coords:
                if a is not None and b is not None:
                    max_span = max(max_span, distance(a, b))
        length_hint = max(7, min(12, int(round(max_span / 3.4)) or 9))
        motif = propose_motif_from_features(frac, length_hint=length_hint)
        clusters.append({
            'size': len(members),
            'composition': dict(comp),
            'fractions': frac,
            'span_A': round(max_span, 2) if max_span else None,
            'length_hint': length_hint,
            'motif': motif
        })

    # Fallback: if nothing clustered but we mapped residues, make one patch
    if not clusters and epitope_residues:
        members = [r for (_, r) in epitope_residues]
        comp, frac = epitope_features(members, asa_by_res=None)
        motif = propose_motif_from_features(frac, length_hint=9)
        clusters = [{
            'size': len(members),
            'composition': dict(comp),
            'fractions': frac,
            'span_A': None,
            'length_hint': 9,
            'motif': motif
        }]

    # ---- Write outputs ----
    csv_rows = []
    for i, cl in enumerate(clusters, 1):
        csv_rows.append({
            'pdb': pdb_id,
            'cluster_id': i,
            'n_residues': cl['size'],
            'span_A': cl['span_A'],
            'length_hint': cl['length_hint'],
            'motif': cl['motif'],
            'frac_negative': round(cl['fractions'].get('negative', 0), 3),
            'frac_positive': round(cl['fractions'].get('positive', 0), 3),
            'frac_aromatic': round(cl['fractions'].get('aromatic', 0), 3),
            'frac_hydrophobic': round(cl['fractions'].get('hydrophobic', 0), 3),
            'frac_polar': round(cl['fractions'].get('polar', 0), 3),
        })

    csv_path = os.path.join(out_dir, f"{pdb_id}_epitope_motifs.csv")
    json_path = os.path.join(out_dir, f"{pdb_id}_epitope_motifs.json")

    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump({'pdb': pdb_id, 'n_clusters': len(clusters), 'clusters': clusters}, f, indent=2)

    return csv_path, json_path

# === Load mapping file globally ===
mapping_tsv = "data/protein_structures/analysis1_params1_fasta_structure_mapping_final.tsv"
SEQ2PDB = load_seq_to_pdb_mapping(mapping_tsv)

# ---------- Hardcoded run ----------
if __name__ == "__main__":
    epitope_tsv = "results/analysis1_params1/protein_analysis/sequences_with_structure/epitope_predictions_bepipred/bamA/bamA_5OR1_linear_epitopes.tsv"
    pdb_dir     = "data/protein_structures/pdb_files"
    out_dir     = "results/analysis1_params1/protein_analysis/sequences_with_structure/epitope_motif_suggestions"

    epi_map = parse_epitope_blocks(epitope_tsv)
    if not epi_map:
        raise SystemExit("No epitopes parsed from input. Please check the file format.")

    os.makedirs(out_dir, exist_ok=True)
    manifest = []
    for pdb_id, df in epi_map.items():
        try:
            csv_path, json_path = analyze_pdb_epitopes(pdb_id.upper(), df, pdb_dir, out_dir)
            manifest.append({'pdb': pdb_id.upper(), 'csv': csv_path, 'json': json_path, 'status': 'ok'})
        except Exception as e:
            manifest.append({'pdb': pdb_id.upper(), 'csv': None, 'json': None, 'status': f'error: {e}'})
    man_path = os.path.join(out_dir, "manifest.json")
    with open(man_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, indent=2))
