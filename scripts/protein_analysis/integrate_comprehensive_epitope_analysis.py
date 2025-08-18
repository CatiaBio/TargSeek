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
            epitope_seq = row.get('Epitope_Sequence', row.get('sequence', ''))
            epitope_start = int(row.get('Epitope_Start', row.get('start', 0)))
            epitope_end = int(row.get('Epitope_End', row.get('end', 0)))
            epitope_length = len(epitope_seq)
            
            # 1. BepiPred Score (immunogenicity prediction)
            bepipred_score = float(row.get('Max_Score', row.get('max_score', 0.5)))
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
    """Create comprehensive HTML report with integrated analysis."""
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Comprehensive Epitope Analysis Report - {analysis}_{paramset} gram_{group}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 0;
                padding: 20px;
                background-color: #f8f9fa;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 30px;
                border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            .header {{
                text-align: center;
                border-bottom: 3px solid #007bff;
                padding-bottom: 20px;
                margin-bottom: 30px;
            }}
            .header h1 {{
                color: #007bff;
                margin-bottom: 10px;
            }}
            .summary-stats {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin-bottom: 30px;
            }}
            .stat-card {{
                background: linear-gradient(135deg, #007bff, #0056b3);
                color: white;
                padding: 20px;
                border-radius: 8px;
                text-align: center;
            }}
            .stat-number {{
                font-size: 2em;
                font-weight: bold;
                margin-bottom: 5px;
            }}
            .stat-label {{
                font-size: 0.9em;
                opacity: 0.9;
            }}
            .section {{
                margin-bottom: 40px;
            }}
            .section h2 {{
                color: #333;
                border-left: 4px solid #007bff;
                padding-left: 15px;
                margin-bottom: 20px;
            }}
            .protein-section {{
                background-color: #f8f9fa;
                padding: 20px;
                border-radius: 8px;
                margin-bottom: 20px;
            }}
            .protein-header {{
                background-color: #007bff;
                color: white;
                padding: 10px 15px;
                border-radius: 5px;
                margin-bottom: 15px;
            }}
            .epitope-table {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
                font-size: 0.9em;
            }}
            .epitope-table th {{
                background-color: #343a40;
                color: white;
                padding: 10px;
                text-align: left;
            }}
            .epitope-table td {{
                padding: 8px 10px;
                border-bottom: 1px solid #ddd;
            }}
            .epitope-table tr:nth-child(even) {{
                background-color: #f8f9fa;
            }}
            .score-cell {{
                text-align: center;
                font-weight: bold;
            }}
            .score-high {{ color: #28a745; }}
            .score-medium {{ color: #ffc107; }}
            .score-low {{ color: #dc3545; }}
            .recommendation-high {{ 
                background-color: #d4edda; 
                color: #155724; 
                padding: 4px 8px; 
                border-radius: 4px; 
                font-weight: bold;
            }}
            .recommendation-medium {{ 
                background-color: #fff3cd; 
                color: #856404; 
                padding: 4px 8px; 
                border-radius: 4px; 
            }}
            .recommendation-low {{ 
                background-color: #f8d7da; 
                color: #721c24; 
                padding: 4px 8px; 
                border-radius: 4px; 
            }}
            .analysis-summary {{
                background-color: #e7f3ff;
                padding: 20px;
                border-radius: 8px;
                margin-bottom: 20px;
            }}
            .methodology {{
                background-color: #f0f0f0;
                padding: 20px;
                border-radius: 8px;
                margin-top: 30px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>Comprehensive Epitope Analysis Report</h1>
                <p><strong>Analysis:</strong> {analysis}_{paramset} | <strong>Gram Type:</strong> {group.title()} | <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
    """
    
    # Calculate overall statistics
    total_epitopes = sum(len(df) for df in scored_epitopes.values())
    total_proteins = len(scored_epitopes)
    
    highly_recommended = sum(
        len(df[df['Recommendation'].str.contains('Highly Recommended', na=False)])
        for df in scored_epitopes.values()
    )
    
    recommended = sum(
        len(df[df['Recommendation'].str.contains('^Recommended$', na=False, regex=True)])
        for df in scored_epitopes.values()
    )
    
    avg_combined_score = np.mean([
        df['Combined_Score'].mean() 
        for df in scored_epitopes.values() if not df.empty
    ]) if scored_epitopes else 0
    
    # Summary statistics section
    html_content += f"""
            <div class="summary-stats">
                <div class="stat-card">
                    <div class="stat-number">{total_proteins}</div>
                    <div class="stat-label">Proteins Analyzed</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{total_epitopes}</div>
                    <div class="stat-label">Total Epitopes</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{highly_recommended}</div>
                    <div class="stat-label">Highly Recommended</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{recommended}</div>
                    <div class="stat-label">Recommended</div>
                </div>
                <div class="stat-card">
                    <div class="stat-number">{avg_combined_score:.2f}</div>
                    <div class="stat-label">Avg Combined Score</div>
                </div>
            </div>

            <div class="analysis-summary">
                <h3>Analysis Overview</h3>
                <p>This report integrates results from multiple complementary analysis tools:</p>
                <ul>
                    <li><strong>BepiPred 3.0:</strong> B-cell epitope immunogenicity prediction</li>
                    <li><strong>ConSurf:</strong> Evolutionary conservation analysis</li>
                    <li><strong>DeepTMHMM:</strong> Membrane topology and surface accessibility</li>
                    <li><strong>BLAST:</strong> Cross-reactivity assessment against host proteomes</li>
                </ul>
                <p>Epitopes are ranked using a combined score incorporating all four analyses, with recommendations for diagnostic applications.</p>
            </div>
    """
    
    # Top candidates section
    html_content += """
            <div class="section">
                <h2>🏆 Top Epitope Candidates</h2>
    """
    
    # Collect all epitopes and sort by combined score
    all_epitopes = []
    for protein_id, df in scored_epitopes.items():
        for idx, row in df.iterrows():
            epitope_data = row.to_dict()
            epitope_data['Protein_ID'] = protein_id
            all_epitopes.append(epitope_data)
    
    all_epitopes.sort(key=lambda x: x['Combined_Score'], reverse=True)
    top_epitopes = all_epitopes[:20]  # Top 20
    
    if top_epitopes:
        html_content += """
                <table class="epitope-table">
                    <tr>
                        <th>Rank</th>
                        <th>Protein</th>
                        <th>Sequence</th>
                        <th>Length</th>
                        <th>Position</th>
                        <th>BepiPred</th>
                        <th>Conservation</th>
                        <th>Topology</th>
                        <th>Specificity</th>
                        <th>Combined</th>
                        <th>Recommendation</th>
                    </tr>
        """
        
        for i, epitope in enumerate(top_epitopes):
            def score_class(score):
                if score >= 0.7: return 'score-high'
                elif score >= 0.4: return 'score-medium'
                else: return 'score-low'
            
            def recommendation_class(rec):
                if 'Highly Recommended' in rec: return 'recommendation-high'
                elif 'Recommended' in rec: return 'recommendation-medium'
                else: return 'recommendation-low'
            
            sequence = epitope.get('Epitope_Sequence', epitope.get('sequence', 'Unknown'))[:20]
            start = epitope.get('Epitope_Start', epitope.get('start', 0))
            end = epitope.get('Epitope_End', epitope.get('end', 0))
            
            html_content += f"""
                    <tr>
                        <td><strong>{i+1}</strong></td>
                        <td>{epitope['Protein_ID']}</td>
                        <td style="font-family: monospace;">{sequence}{'...' if len(sequence) == 20 else ''}</td>
                        <td>{len(epitope.get('Epitope_Sequence', epitope.get('sequence', '')))}</td>
                        <td>{start}-{end}</td>
                        <td class="score-cell {score_class(epitope['BepiPred_Score'])}">{epitope['BepiPred_Score']:.2f}</td>
                        <td class="score-cell {score_class(epitope['Conservation_Score'])}">{epitope['Conservation_Score']:.2f}</td>
                        <td class="score-cell {score_class(epitope['Topology_Score'])}">{epitope['Topology_Score']:.2f}</td>
                        <td class="score-cell {score_class(epitope['Specificity_Score'])}">{epitope['Specificity_Score']:.2f}</td>
                        <td class="score-cell {score_class(epitope['Combined_Score'])}">{epitope['Combined_Score']:.2f}</td>
                        <td><span class="{recommendation_class(epitope['Recommendation'])}">{epitope['Recommendation']}</span></td>
                    </tr>
            """
        
        html_content += """
                </table>
            </div>
        """
    
    # Detailed protein analysis
    html_content += """
            <div class="section">
                <h2>📊 Detailed Protein Analysis</h2>
    """
    
    for protein_id, df in scored_epitopes.items():
        if df.empty:
            continue
            
        protein_stats = {
            'total_epitopes': len(df),
            'avg_combined_score': df['Combined_Score'].mean(),
            'highly_recommended': len(df[df['Recommendation'].str.contains('Highly Recommended', na=False)]),
            'recommended': len(df[df['Recommendation'].str.contains('^Recommended$', na=False, regex=True)])
        }
        
        html_content += f"""
                <div class="protein-section">
                    <div class="protein-header">
                        <h3>{protein_id}</h3>
                    </div>
                    <p><strong>Total Epitopes:</strong> {protein_stats['total_epitopes']} | 
                       <strong>Avg Score:</strong> {protein_stats['avg_combined_score']:.2f} | 
                       <strong>Highly Recommended:</strong> {protein_stats['highly_recommended']} | 
                       <strong>Recommended:</strong> {protein_stats['recommended']}</p>
        """
        
        # Show top 10 epitopes for this protein
        top_protein_epitopes = df.head(10)
        
        html_content += """
                    <table class="epitope-table">
                        <tr>
                            <th>Rank</th>
                            <th>Sequence</th>
                            <th>Position</th>
                            <th>Length</th>
                            <th>Combined Score</th>
                            <th>Recommendation</th>
                        </tr>
        """
        
        for idx, row in top_protein_epitopes.iterrows():
            sequence = str(row.get('Epitope_Sequence', row.get('sequence', 'Unknown')))
            start = row.get('Epitope_Start', row.get('start', 0))
            end = row.get('Epitope_End', row.get('end', 0))
            
            def recommendation_class(rec):
                if 'Highly Recommended' in rec: return 'recommendation-high'
                elif 'Recommended' in rec: return 'recommendation-medium'
                else: return 'recommendation-low'
            
            def score_class(score):
                if score >= 0.7: return 'score-high'
                elif score >= 0.4: return 'score-medium'
                else: return 'score-low'
            
            html_content += f"""
                        <tr>
                            <td><strong>{row['Rank']}</strong></td>
                            <td style="font-family: monospace;">{sequence[:30]}{'...' if len(sequence) > 30 else ''}</td>
                            <td>{start}-{end}</td>
                            <td>{len(sequence)}</td>
                            <td class="score-cell {score_class(row['Combined_Score'])}">{row['Combined_Score']:.2f}</td>
                            <td><span class="{recommendation_class(row['Recommendation'])}">{row['Recommendation']}</span></td>
                        </tr>
            """
        
        html_content += """
                    </table>
                </div>
        """
    
    # Methodology section
    html_content += f"""
            <div class="methodology">
                <h2>🔬 Methodology</h2>
                <h3>Scoring System</h3>
                <ul>
                    <li><strong>BepiPred Score (30%):</strong> B-cell epitope immunogenicity prediction from BepiPred 3.0</li>
                    <li><strong>Conservation Score (25%):</strong> Evolutionary conservation from ConSurf analysis (inverse of amino acid variability)</li>
                    <li><strong>Topology Score (25%):</strong> Surface accessibility based on DeepTMHMM membrane topology predictions</li>
                    <li><strong>Specificity Score (20%):</strong> Cross-reactivity avoidance based on BLAST similarity to host proteomes</li>
                </ul>
                
                <h3>Recommendation Categories</h3>
                <ul>
                    <li><strong>Highly Recommended (≥0.8):</strong> Excellent candidates for diagnostic applications</li>
                    <li><strong>Recommended (≥0.6):</strong> Good candidates worth further investigation</li>
                    <li><strong>Moderate Candidate (≥0.4):</strong> May be useful but require careful evaluation</li>
                    <li><strong>Low Priority (<0.4):</strong> Not recommended for diagnostic use</li>
                </ul>
                
                <h3>Data Sources</h3>
                <p><strong>Analysis ID:</strong> {analysis}_{paramset}<br>
                <strong>Gram Classification:</strong> {group.title()}<br>
                <strong>Proteins Analyzed:</strong> {', '.join(scored_epitopes.keys())}<br>
                <strong>Generation Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
        </div>
    </body>
    </html>
    """
    
    # Write HTML report
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
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