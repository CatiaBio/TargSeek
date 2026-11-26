#!/usr/bin/env python3
"""
Comprehensive epitope analysis integration script.

This script integrates data from multiple analysis tools to provide a comprehensive
assessment of epitope candidates:
1. BepiPred 3.0 epitope predictions (topology-filtered)
2. ConSurf conservation analysis 
3. Topology predictions (extracellular regions)
4. BLAST cross-reactivity analysis

Creates a comprehensive HTML report with rankings, visualizations, and recommendations
for the best epitope candidates for diagnostic applications.

Input:
- Topology-filtered epitope tables
- ConSurf conservation scores
- Topology mapping data
- BLAST cross-reactivity results

Output:
- Comprehensive HTML report with integrated analysis
- Top epitope candidates ranked by combined score
- Summary statistics and visualizations
"""

import os
import sys
import re
import json
import pandas as pd
import numpy as np
import logging
import glob
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_epitope_data(filtered_epitopes_dir: str) -> Dict[str, pd.DataFrame]:
    """Load topology-filtered epitope tables."""
    epitope_data = {}
    
    pattern = os.path.join(filtered_epitopes_dir, "*", "*_filtered_epitopes.csv")
    epitope_files = glob.glob(pattern)
    
    if not epitope_files:
        pattern = os.path.join(filtered_epitopes_dir, "*_filtered_epitopes.csv")
        epitope_files = glob.glob(pattern)
    
    for file_path in epitope_files:
        try:
            basename = os.path.basename(file_path)
            protein_id = basename.replace('_filtered_epitopes.csv', '')
            
            df = pd.read_csv(file_path)
            if not df.empty:
                epitope_data[protein_id] = df
                logger.info(f"Loaded {len(df)} epitopes for {protein_id}")
                
        except Exception as e:
            logger.error(f"Error loading epitope table {file_path}: {e}")
    
    logger.info(f"Loaded epitope data for {len(epitope_data)} proteins")
    return epitope_data

def load_consurf_data(consurf_dir: str, group: str) -> Dict[str, Dict]:
    """Load ConSurf conservation analysis results."""
    consurf_data = {}
    consurf_group_dir = os.path.join(consurf_dir, f"gram_{group}")
    
    if not os.path.exists(consurf_group_dir):
        logger.warning(f"ConSurf directory not found: {consurf_group_dir}")
        return consurf_data
    
    for gene_dir in glob.glob(os.path.join(consurf_group_dir, "*")):
        if os.path.isdir(gene_dir):
            gene_name = os.path.basename(gene_dir)
            
            # Look for conservation data files
            conservation_file = os.path.join(gene_dir, "msa_aa_variety_percentage.csv")
            grades_file = os.path.join(gene_dir, "consurf_summary.json")
            
            if os.path.exists(conservation_file):
                try:
                    conservation_df = pd.read_csv(conservation_file)
                    
                    gene_conservation = {
                        'conservation_data': conservation_df,
                        'summary': {},
                        'grades': {}
                    }
                    
                    # Load ConSurf grades if available
                    if os.path.exists(grades_file):
                        with open(grades_file, 'r') as f:
                            grades_data = json.load(f)
                            gene_conservation['grades'] = grades_data
                    
                    # Calculate summary statistics
                    if 'Variety_Percentage' in conservation_df.columns:
                        variety_scores = conservation_df['Variety_Percentage'].dropna()
                        gene_conservation['summary'] = {
                            'mean_variety': variety_scores.mean(),
                            'median_variety': variety_scores.median(),
                            'std_variety': variety_scores.std(),
                            'highly_conserved_positions': len(variety_scores[variety_scores < 10]),
                            'variable_positions': len(variety_scores[variety_scores > 50]),
                            'total_positions': len(variety_scores)
                        }
                    
                    consurf_data[gene_name] = gene_conservation
                    logger.info(f"Loaded ConSurf data for {gene_name}")
                    
                except Exception as e:
                    logger.error(f"Error loading ConSurf data for {gene_name}: {e}")
    
    logger.info(f"Loaded ConSurf data for {len(consurf_data)} genes")
    return consurf_data

def load_blast_results(blast_summary_file: str) -> Dict[str, Any]:
    """Load BLAST cross-reactivity analysis results."""
    blast_data = {}
    
    if os.path.exists(blast_summary_file):
        try:
            with open(blast_summary_file, 'r') as f:
                blast_data = json.load(f)
            logger.info("Loaded BLAST cross-reactivity data")
        except Exception as e:
            logger.error(f"Error loading BLAST results: {e}")
    else:
        logger.warning(f"BLAST results file not found: {blast_summary_file}")
    
    return blast_data

def load_topology_mapping(topology_mapping_file: str) -> Dict[str, Any]:
    """Load extracellular topology mapping data."""
    topology_data = {}
    
    if os.path.exists(topology_mapping_file):
        try:
            with open(topology_mapping_file, 'r') as f:
                topology_data = json.load(f)
            logger.info("Loaded topology mapping data")
        except Exception as e:
            logger.error(f"Error loading topology mapping: {e}")
    else:
        logger.warning(f"Topology mapping file not found: {topology_mapping_file}")
    
    return topology_data

def _field(row, *names, default=None):
    """First present, non-empty value among several possible column names.

    The epitope tables written by create_epitope_tables_from_bepipred.py use
    Peptide / Start / End / Score; earlier code looked only for
    Epitope_Sequence / Epitope_Start / Max_Score and so silently fell back to
    defaults for every row.
    """
    for name in names:
        if name in row:
            value = row[name]
            if value is None:
                continue
            try:
                if pd.isna(value):
                    continue
            except (TypeError, ValueError):
                pass
            if value != '':
                return value
    return default


def calculate_epitope_scores(epitope_data: Dict[str, pd.DataFrame],
                           consurf_data: Dict[str, Dict],
                           blast_data: Dict[str, Any],
                           topology_data: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    """Calculate comprehensive scores for each epitope."""
    
    scored_epitopes = {}
    
    for protein_id, epitope_df in epitope_data.items():
        # Create copy of epitope dataframe
        scored_df = epitope_df.copy()
        
        # Extract gene name from protein_id (assuming format like "bamA_5OR1")
        gene_name = protein_id.split('_')[0]
        
        # Initialize score columns
        scored_df['BepiPred_Score'] = 0.0
        scored_df['Conservation_Score'] = 0.0
        scored_df['Topology_Score'] = 0.0
        scored_df['Specificity_Score'] = 0.0
        scored_df['Combined_Score'] = 0.0
        scored_df['Rank'] = 0
        scored_df['Recommendation'] = 'Unknown'
        
        for idx, row in scored_df.iterrows():
            epitope_seq = str(_field(row, 'Peptide', 'Epitope_Sequence', 'sequence',
                                     default='') or '')
            epitope_start = int(_field(row, 'Start', 'Epitope_Start', 'start', default=0) or 0)
            epitope_end = int(_field(row, 'End', 'Epitope_End', 'end', default=0) or 0)
            epitope_length = len(epitope_seq)

            # 1. BepiPred Score (immunogenicity prediction)
            bepipred_score = float(_field(row, 'Score', 'Max_Score', 'max_score',
                                          default=0.5) or 0.5)
            scored_df.at[idx, 'BepiPred_Score'] = bepipred_score
            
            # 2. Conservation Score (from ConSurf)
            conservation_score = 0.5  # Default neutral score
            if gene_name in consurf_data and 'conservation_data' in consurf_data[gene_name]:
                conservation_df = consurf_data[gene_name]['conservation_data']
                
                # Find overlapping positions
                if 'Position' in conservation_df.columns and 'Variety_Percentage' in conservation_df.columns:
                    overlapping_positions = conservation_df[
                        (conservation_df['Position'] >= epitope_start) & 
                        (conservation_df['Position'] <= epitope_end)
                    ]
                    
                    if not overlapping_positions.empty:
                        # Use inverse of variety percentage (higher conservation = lower variety)
                        avg_variety = overlapping_positions['Variety_Percentage'].mean()
                        conservation_score = max(0, (100 - avg_variety) / 100)  # Convert to 0-1 scale
            
            scored_df.at[idx, 'Conservation_Score'] = conservation_score
            
            # 3. Topology Score (extracellular accessibility)
            topology_score = 0.5  # Default neutral score
            if 'mapping' in topology_data and protein_id in topology_data['mapping']:
                protein_topology = topology_data['mapping'][protein_id]
                extracellular_coverage = protein_topology.get('extracellular_coverage', 0)
                
                # Score based on overall protein extracellular coverage
                topology_score = min(1.0, extracellular_coverage * 2)  # Scale to emphasize high coverage
            
            scored_df.at[idx, 'Topology_Score'] = topology_score
            
            # 4. Specificity Score (from BLAST analysis)
            specificity_score = 0.5  # Default neutral score
            if 'protein_summaries' in blast_data:
                for blast_protein_id, protein_summary in blast_data['protein_summaries'].items():
                    if blast_protein_id == protein_id:
                        # Use inverse of average cross-reactivity
                        avg_cross_reactivity = protein_summary.get('average_specificity_score', 0.5)
                        specificity_score = max(0, 1 - avg_cross_reactivity)
                        break
            
            scored_df.at[idx, 'Specificity_Score'] = specificity_score
            
            # 5. Calculate Combined Score (weighted average)
            weights = {
                'bepipred': 0.3,      # 30% - immunogenicity prediction
                'conservation': 0.25,  # 25% - evolutionary conservation
                'topology': 0.25,     # 25% - surface accessibility
                'specificity': 0.2    # 20% - cross-reactivity avoidance
            }
            
            combined_score = (
                weights['bepipred'] * bepipred_score +
                weights['conservation'] * conservation_score +
                weights['topology'] * topology_score +
                weights['specificity'] * specificity_score
            )
            
            scored_df.at[idx, 'Combined_Score'] = combined_score
            
            # 6. Generate Recommendation
            if combined_score >= 0.8:
                recommendation = 'Highly Recommended'
            elif combined_score >= 0.6:
                recommendation = 'Recommended'
            elif combined_score >= 0.4:
                recommendation = 'Moderate Candidate'
            else:
                recommendation = 'Low Priority'
            
            # Additional criteria for recommendations
            if epitope_length < 8:
                recommendation = 'Too Short - ' + recommendation
            elif bepipred_score < 0.3:
                recommendation = 'Low Immunogenicity - ' + recommendation
            elif conservation_score < 0.3:
                recommendation = 'Poorly Conserved - ' + recommendation
            
            scored_df.at[idx, 'Recommendation'] = recommendation
        
        # Rank epitopes by combined score
        scored_df = scored_df.sort_values('Combined_Score', ascending=False).reset_index(drop=True)
        scored_df['Rank'] = range(1, len(scored_df) + 1)
        
        scored_epitopes[protein_id] = scored_df
        logger.info(f"Scored {len(scored_df)} epitopes for {protein_id}")
    
    return scored_epitopes

def create_comprehensive_report(scored_epitopes: Dict[str, pd.DataFrame],
                              consurf_data: Dict[str, Dict],
                              blast_data: Dict[str, Any],
                              topology_data: Dict[str, Any],
                              output_file: str,
                              analysis: str,
                              paramset: str,
                              group: str) -> None:
    """Create the integrated HTML report.

    Layout and styling live in report_template.py; this function only decides
    what goes in the page.
    """
    # Snakemake runs this from a temp copy but adds the real script directory to
    # sys.path; keep a fallback for direct execution.
    try:
        import report_template as rt
    except ImportError:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import report_template as rt

    group_label = f"Gram-{group}"
    generated = datetime.now().strftime('%Y-%m-%d %H:%M')

    # Which score components were actually measured. When a tool did not run,
    # calculate_epitope_scores falls back to a neutral default; printing that as
    # a number would pass an assumption off as a measurement.
    has_conservation = bool(consurf_data)
    has_specificity = bool(
        (blast_data or {}).get('detail') or (blast_data or {}).get('protein_summaries'))
    has_topology = bool((topology_data or {}).get('mapping'))

    def comp(value, measured):
        return rt.score_num(value) if measured else '<span class="muted" title="not measured">n/m</span>'

    # --- gather every epitope, best first ---------------------------------
    all_epitopes = []
    for protein_id, df in scored_epitopes.items():
        if df is None or df.empty:
            continue
        for _, row in df.iterrows():
            rec = row.to_dict()
            rec['Protein_ID'] = protein_id
            all_epitopes.append(rec)
    all_epitopes.sort(key=lambda r: r.get('Combined_Score', 0), reverse=True)

    def seq_of(rec):
        return str(_field(rec, 'Peptide', 'Epitope_Sequence', 'sequence', default='') or '')

    def pos_of(rec):
        start = _field(rec, 'Start', 'Epitope_Start', 'start', default='')
        end = _field(rec, 'End', 'Epitope_End', 'end', default='')
        if start == '' and end == '':
            return '-'
        return f"{start}-{end}"

    total_proteins = len([d for d in scored_epitopes.values() if d is not None and not d.empty])
    total_epitopes = len(all_epitopes)
    highly = sum(1 for r in all_epitopes
                 if 'highly recommended' in str(r.get('Recommendation', '')).lower())
    recommended = sum(1 for r in all_epitopes
                      if str(r.get('Recommendation', '')).strip().lower() == 'recommended')
    avg_score = (sum(r.get('Combined_Score', 0) for r in all_epitopes) / total_epitopes
                 if total_epitopes else 0.0)

    sections = []
    nav_main = [('overview', 'Overview', None)]

    # --- overview ---------------------------------------------------------
    tiles = rt.stat_row([
        rt.stat_tile(total_proteins, 'Proteins'),
        rt.stat_tile(total_epitopes, 'Epitopes'),
        rt.stat_tile(highly, 'Highly recommended'),
        rt.stat_tile(recommended, 'Recommended'),
        rt.stat_tile(f"{avg_score:.2f}", 'Mean combined score', 'of 1.00'),
    ])

    overview_body = tiles
    if total_epitopes == 0:
        overview_body += rt.note(
            'warning', 'No extracellular epitopes to report',
            '<p>Every predicted epitope was removed by the topology filter, so there is '
            'nothing to rank. That is usually worth checking rather than accepting: look at '
            '<code>epitope_filter_report.json</code> for how many epitopes went in, and at '
            'the extracellular mapping for the proteins involved.</p>')
    else:
        overview_body += (
            '<h3>How a candidate is scored</h3>'
            '<p>Each epitope carries four independent scores, combined into one ranking '
            'value. The table columns below show all of them, so a high combined score '
            'built on one strong component is visible rather than hidden.</p>'
            '<div class="chiprow">'
            '<span class="chip">BepiPred 3.0 &mdash; immunogenicity (30%)</span>'
            '<span class="chip">ConSurf &mdash; conservation (25%)</span>'
            '<span class="chip">DeepTMHMM &mdash; surface accessibility (25%)</span>'
            '<span class="chip">BLAST &mdash; specificity vs host proteomes (20%)</span>'
            '</div>')
    missing = [n for n, ok in (('ConSurf conservation', has_conservation),
                              ('DeepTMHMM topology', has_topology),
                              ('BLAST specificity', has_specificity)) if not ok]
    if missing and total_epitopes:
        overview_body += rt.note(
            'warning', 'Some scores were not measured',
            '<p>' + ', '.join(missing) + ' did not run for this analysis, so those columns '
            'show <strong>n/m</strong> (not measured). The combined score still uses the '
            'pipeline\'s neutral defaults for them, so treat the ranking as provisional '
            'until those steps have run.</p>')

    sections.append(rt.section('overview', f'{group_label} epitope candidates',
                               overview_body,
                               subtitle=f'{analysis}_{paramset} | generated {generated}'))

    # --- top candidates ---------------------------------------------------
    if all_epitopes:
        top = all_epitopes[:25]
        rows = []
        for i, rec in enumerate(top, 1):
            seq = seq_of(rec)
            rows.append([
                f'<strong>{i}</strong>',
                rt.esc(rec.get('Protein_ID', '')),
                f'<span class="seq">{rt.esc(seq)}</span>',
                rt.esc(len(seq)),
                rt.esc(pos_of(rec)),
                rt.score_num(rec.get('BepiPred_Score')),
                comp(rec.get('Conservation_Score'), has_conservation),
                comp(rec.get('Topology_Score'), has_topology),
                comp(rec.get('Specificity_Score'), has_specificity),
                rt.score_meter(rec.get('Combined_Score')),
                rt.recommendation_badge(rec.get('Recommendation')),
            ])
        headers = ['#', 'Protein', 'Sequence', 'Len', 'Position', 'BepiPred',
                   'Conservation', 'Topology', 'Specificity', 'Combined', 'Assessment']
        body = (f'<p class="muted">The {len(top)} highest-scoring of {total_epitopes} '
                f'epitopes. Click a column heading to re-sort; use the filter box above '
                f'to narrow by gene or sequence.</p>'
                + rt.table(headers, rows, table_id='top-candidates', numeric_from=3))
        sections.append(rt.section('top-candidates', 'Top candidates', body))
        nav_main.append(('top-candidates', 'Top candidates', len(top)))

    # --- cross-reactivity -------------------------------------------------
    if blast_data:
        summary = blast_data.get('summary', {}) if isinstance(blast_data, dict) else {}
        detail = blast_data.get('detail', {}) if isinstance(blast_data, dict) else {}
        per_protein = detail.get('epitope_analysis', {}) if isinstance(detail, dict) else {}

        # Flatten every epitope's BLAST verdict
        flagged, counts = [], {'highly_specific': 0, 'moderately_specific': 0, 'cross_reactive': 0}
        for protein_id, pdata in per_protein.items():
            for _, epi in (pdata.get('epitopes', {}) or {}).items():
                cls = epi.get('classification', 'unknown')
                if cls in counts:
                    counts[cls] += 1
                if cls == 'highly_specific':
                    continue
                # blast_results is keyed by proteome: {db: {significant_hits, best_hit, hits}}
                per_db = epi.get('blast_results', {}) or {}
                for db, dbres in per_db.items():
                    if not isinstance(dbres, dict):
                        continue
                    hits = dbres.get('hits') or []
                    # a hit only counts if it cleared the identity/coverage filter
                    hit = hits[0] if hits else None
                    if hit is None:
                        continue
                    flagged.append({
                        'protein': protein_id,
                        'sequence': epi.get('sequence', ''),
                        'classification': cls,
                        'db': str(db).replace('_proteome', ''),
                        'identity': hit.get('percent_identity'),
                        'evalue': hit.get('evalue'),
                        'n_hits': dbres.get('significant_hits', len(hits)),
                        'subject': (hit.get('subject_title') or hit.get('subject_id') or '')[:70],
                    })

        body = ('<p class="muted">Every epitope is searched with BLAST against the reference '
                'proteomes in <code>data/blast_databases/</code>. A strong match to a host '
                'protein means an antibody raised against that epitope may bind the host too, '
                'which is what a diagnostic must avoid.</p>')

        if any(counts.values()):
            body += rt.stat_row([
                rt.stat_tile(counts['highly_specific'], 'Highly specific', 'no host match'),
                rt.stat_tile(counts['moderately_specific'], 'Moderately specific'),
                rt.stat_tile(counts['cross_reactive'], 'Cross-reactive', 'review before use'),
            ])

        if flagged:
            flagged.sort(key=lambda r: -(r['identity'] or 0))
            rows = []
            for r in flagged[:40]:
                tag = ('critical' if r['classification'] == 'cross_reactive' else 'warning')
                label = ('Cross-reactive' if r['classification'] == 'cross_reactive'
                         else 'Moderate')
                ident = f"{r['identity']:.0f}%" if isinstance(r['identity'], (int, float)) else '-'
                ev = f"{r['evalue']:.1e}" if isinstance(r['evalue'], (int, float)) else '-'
                rows.append([
                    rt.esc(r['protein']),
                    f'<span class="seq">{rt.esc(r["sequence"])}</span>',
                    rt.esc(r['db']),
                    f'<span class="badge badge-{tag}">{rt.icon("alert")}{label}</span>',
                    rt.esc(ident),
                    rt.esc(ev),
                    rt.esc(r.get('n_hits', '')),
                    f'<span class="muted">{rt.esc(r["subject"])}</span>',
                ])
            body += rt.table(
                ['Protein', 'Epitope', 'Proteome', 'Verdict', 'Identity', 'E-value',
                 'Matches', 'Best match'],
                rows, table_id='xreact', numeric_from=4)
            if len(flagged) > 40:
                body += f'<p class="muted">Showing the 40 strongest of {len(flagged)} host matches.</p>'
        elif per_protein:
            body += rt.note('good', 'No host matches flagged',
                            '<p>No epitope matched a reference proteome strongly enough to be '
                            'classed as cross-reactive.</p>')

        if summary:
            items = ''.join(
                f'<dt>{rt.esc(k.replace("_", " ").capitalize())}</dt><dd>{rt.esc(v)}</dd>'
                for k, v in summary.items() if not isinstance(v, (dict, list)))
            if items:
                body += f'<h3>Run summary</h3><dl class="dl">{items}</dl>'

        sections.append(rt.section('cross-reactivity', 'Cross-reactivity screen', body))
        nav_main.append(('cross-reactivity', 'Cross-reactivity',
                         counts['cross_reactive'] or None))

    # --- per protein ------------------------------------------------------
    nav_proteins = []
    for protein_id, df in sorted(scored_epitopes.items()):
        if df is None or df.empty:
            continue
        anchor = 'p-' + re.sub(r'[^A-Za-z0-9_-]', '-', str(protein_id))
        p_high = len(df[df['Recommendation'].astype(str).str.contains('Highly Recommended', na=False)])
        p_mean = df['Combined_Score'].mean()
        rows = []
        for i, (_, row) in enumerate(df.head(15).iterrows(), 1):
            rec = row.to_dict()
            seq = seq_of(rec)
            rows.append([
                f'<strong>{i}</strong>',
                f'<span class="seq">{rt.esc(seq)}</span>',
                rt.esc(pos_of(rec)),
                rt.esc(len(seq)),
                rt.score_num(rec.get('BepiPred_Score')),
                comp(rec.get('Conservation_Score'), has_conservation),
                comp(rec.get('Topology_Score'), has_topology),
                comp(rec.get('Specificity_Score'), has_specificity),
                rt.score_meter(rec.get('Combined_Score')),
                rt.recommendation_badge(rec.get('Recommendation')),
            ])
        headers = ['#', 'Sequence', 'Position', 'Len', 'BepiPred', 'Conservation',
                   'Topology', 'Specificity', 'Combined', 'Assessment']
        chips = (f'<div class="chiprow"><span class="chip">{len(df)} epitopes</span>'
                 f'<span class="chip">mean {p_mean:.2f}</span>'
                 f'<span class="chip">{p_high} highly recommended</span></div>')
        body = chips + rt.table(headers, rows, table_id=anchor + '-t', numeric_from=3)
        if len(df) > 15:
            body += f'<p class="muted">Showing the top 15 of {len(df)} epitopes.</p>'
        sections.append(rt.section(anchor, str(protein_id), body, level=2))
        nav_proteins.append((anchor, str(protein_id), len(df)))

    # --- methods ----------------------------------------------------------
    proteins_listed = ', '.join(sorted(scored_epitopes.keys())) or 'none'
    methods = f"""
    <h3>Scoring</h3>
    <p>The combined score is a weighted sum of four normalised components:</p>
    <ul>
      <li><strong>BepiPred 3.0 (30%)</strong> &mdash; B-cell epitope immunogenicity.</li>
      <li><strong>ConSurf (25%)</strong> &mdash; evolutionary conservation, taken as the
          inverse of amino-acid variability across the alignment.</li>
      <li><strong>DeepTMHMM (25%)</strong> &mdash; surface accessibility from the predicted
          membrane topology.</li>
      <li><strong>BLAST (20%)</strong> &mdash; specificity: similarity to host and
          reference proteomes reduces this score.</li>
    </ul>
    <h3>Assessment thresholds</h3>
    <ul>
      <li>{rt.recommendation_badge('Highly Recommended')} combined score &ge; 0.80</li>
      <li>{rt.recommendation_badge('Recommended')} &ge; 0.60</li>
      <li>{rt.recommendation_badge('Moderate')} &ge; 0.40</li>
      <li>{rt.recommendation_badge('Low Priority')} below 0.40</li>
    </ul>
    <p class="muted">These are ranking aids, not evidence of immunogenicity. Any candidate
    taken forward needs experimental validation.</p>
    <h3>This run</h3>
    <dl class="dl">
      <dt>Analysis</dt><dd>{rt.esc(analysis)}_{rt.esc(paramset)}</dd>
      <dt>Group</dt><dd>{rt.esc(group_label)}</dd>
      <dt>Proteins</dt><dd>{rt.esc(proteins_listed)}</dd>
      <dt>Generated</dt><dd>{rt.esc(generated)}</dd>
    </dl>
    """
    sections.append(rt.section('methods', 'Methods', methods))
    nav_main.append(('methods', 'Methods', None))

    nav_groups = [('Report', nav_main)]
    if nav_proteins:
        nav_groups.append((f'Proteins ({len(nav_proteins)})', nav_proteins))

    page = rt.render_page(
        title=f'Epitope analysis | {group_label}',
        subtitle=f'{group_label} epitope report',
        meta_items=[('Analysis', f'{analysis}_{paramset}'),
                    ('Group', group_label),
                    ('Epitopes', total_epitopes)],
        nav_groups=nav_groups,
        sections_html='\n'.join(sections),
        footer_note=f'TargSeek | {generated}',
    )

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(page)

    logger.info(f"Created comprehensive HTML report: {output_file}")


def main():
    """Main function to integrate comprehensive epitope analysis."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    group = snakemake.params.group
    output_dir = snakemake.params.output_dir
    
    # Input files
    filtered_epitopes_sentinel = snakemake.input.filtered_epitopes_sentinel
    consurf_sentinel = snakemake.input.consurf_sentinel
    blast_results_sentinel = snakemake.input.blast_results_sentinel
    extracellular_mapping = snakemake.input.extracellular_mapping
    
    # Output files
    integrated_sentinel = snakemake.output.integrated_sentinel
    comprehensive_report = snakemake.output.comprehensive_report
    
    logger.info(f"Starting comprehensive epitope analysis integration for {analysis}_{paramset} gram_{group}")
    
    # Find data directories
    filtered_epitopes_dir = os.path.dirname(filtered_epitopes_sentinel)
    consurf_base_dir = os.path.dirname(os.path.dirname(consurf_sentinel))  # Go up two levels
    blast_results_dir = os.path.dirname(blast_results_sentinel)
    
    # Load all data sources
    logger.info("Loading epitope data...")
    epitope_data = load_epitope_data(filtered_epitopes_dir)
    
    logger.info("Loading ConSurf conservation data...")
    consurf_data = load_consurf_data(consurf_base_dir, group)
    
    logger.info("Loading BLAST cross-reactivity data...")
    blast_summary_file = os.path.join(blast_results_dir, "blast_cross_reactivity_summary.json")
    blast_data = load_blast_results(blast_summary_file)
    # The detailed file carries the per-epitope hits (which proteome, which
    # protein, what identity) - that is the part a reader actually needs.
    blast_detail_file = os.path.join(blast_results_dir, "detailed_blast_analysis.json")
    blast_detail = load_blast_results(blast_detail_file)
    if blast_detail:
        blast_data = dict(blast_data or {})
        blast_data["detail"] = blast_detail
    
    logger.info("Loading topology mapping data...")
    topology_data = load_topology_mapping(extracellular_mapping)
    
    if not epitope_data:
        logger.error("No epitope data available for integration")
        sys.exit(1)
    
    # Calculate comprehensive scores
    logger.info("Calculating comprehensive epitope scores...")
    scored_epitopes = calculate_epitope_scores(epitope_data, consurf_data, blast_data, topology_data)
    
    # Create comprehensive report
    logger.info("Creating comprehensive HTML report...")
    create_comprehensive_report(scored_epitopes, consurf_data, blast_data, topology_data, 
                               comprehensive_report, analysis, paramset, group)
    
    # Save scored epitope data as CSV for further analysis
    csv_output_file = os.path.join(output_dir, "comprehensive_epitope_scores.csv")
    all_scored_epitopes = []
    
    for protein_id, df in scored_epitopes.items():
        df_copy = df.copy()
        df_copy['Protein_ID'] = protein_id
        all_scored_epitopes.append(df_copy)
    
    if all_scored_epitopes:
        combined_df = pd.concat(all_scored_epitopes, ignore_index=True)
        combined_df.to_csv(csv_output_file, index=False)
        logger.info(f"Saved comprehensive epitope scores: {csv_output_file}")
    
    # Create sentinel file
    os.makedirs(os.path.dirname(integrated_sentinel), exist_ok=True)
    with open(integrated_sentinel, 'w') as f:
        f.write(f"Comprehensive epitope analysis completed for {analysis}_{paramset} gram_{group}\\n")
        f.write(f"Proteins analyzed: {len(scored_epitopes)}\\n")
        f.write(f"Total epitopes: {sum(len(df) for df in scored_epitopes.values())}\\n")
        f.write(f"Report generated: {comprehensive_report}\\n")
    
    logger.info("Comprehensive epitope analysis integration completed successfully")

if __name__ == "__main__":
    main()