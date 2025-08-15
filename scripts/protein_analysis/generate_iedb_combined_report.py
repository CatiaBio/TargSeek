#!/usr/bin/env python3
"""
Generate Combined IEDB Validation Report
=========================================

This script generates a comprehensive HTML report combining IEDB validation
results from both Gram-positive and Gram-negative epitope analyses.

Input: 
    - IEDB validation results JSON files (positive and negative)
Output:
    - Combined HTML report with comparative analysis
"""

import json
import pandas as pd
from pathlib import Path
import sys
from datetime import datetime

def load_iedb_results(json_file):
    """Load IEDB validation results from JSON file"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        print(f"Error loading {json_file}: {e}")
        return None

def create_summary_dataframe(results_data, gram_type):
    """Create summary DataFrame from IEDB results"""
    summary_rows = []
    
    if not results_data or 'results' not in results_data:
        return pd.DataFrame()
    
    for gene, gene_data in results_data['results'].items():
        for epitope in gene_data.get('combined_scores', []):
            row = {
                'Gram_Type': gram_type,
                'Gene': gene,
                'Epitope_ID': epitope['epitope_id'],
                'Sequence': epitope['sequence'],
                'Length': epitope['length'],
                'BepiPred_Score': epitope['bepipred_score'],
                'Jalview_Similarity': epitope['current_conservation'].get('jalview_similarity', 0),
                'BLOSUM62': epitope['current_conservation'].get('blosum62', 0),
                'Identity': epitope['current_conservation'].get('identity', 0),
                'IEDB_100%': epitope['iedb_scores'].get('iedb_100', 0),
                'IEDB_90%': epitope['iedb_scores'].get('iedb_90', 0),
                'IEDB_80%': epitope['iedb_scores'].get('iedb_80', 0),
                'IEDB_70%': epitope['iedb_scores'].get('iedb_70', 0),
                'Combined_Score': epitope['combined_score'],
                'Rank': epitope['rank']
            }
            summary_rows.append(row)
    
    return pd.DataFrame(summary_rows)

def generate_html_report(positive_results, negative_results, output_file, analysis, paramset):
    """Generate comprehensive HTML report"""
    
    # Load data
    pos_data = load_iedb_results(positive_results) if positive_results.exists() else None
    neg_data = load_iedb_results(negative_results) if negative_results.exists() else None
    
    # Create summary dataframes
    pos_df = create_summary_dataframe(pos_data, "Gram-positive")
    neg_df = create_summary_dataframe(neg_data, "Gram-negative")
    
    # Combine dataframes
    combined_df = pd.concat([pos_df, neg_df], ignore_index=True)
    
    # Generate HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>IEDB Validation Report - {analysis}_{paramset}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            .header {{ background-color: #f0f8ff; padding: 20px; border-radius: 5px; }}
            .section {{ margin: 20px 0; }}
            .summary-box {{ background-color: #f9f9f9; padding: 15px; border-left: 4px solid #007acc; }}
            table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #007acc; color: white; }}
            .highlight {{ background-color: #fff3cd; }}
            .metric {{ display: inline-block; margin: 10px 20px 10px 0; }}
            .metric-value {{ font-size: 24px; font-weight: bold; color: #007acc; }}
            .metric-label {{ font-size: 12px; color: #666; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>🧬 IEDB Epitope Validation Report</h1>
            <h2>Analysis: {analysis}_{paramset}</h2>
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    """
    
    # Executive Summary
    html_content += """
        <div class="section">
            <h2>📊 Executive Summary</h2>
            <div class="summary-box">
    """
    
    if not combined_df.empty:
        total_epitopes = len(combined_df)
        total_genes = combined_df['Gene'].nunique()
        top_candidates = combined_df.nlargest(10, 'Combined_Score')
        avg_combined_score = combined_df['Combined_Score'].mean()
        high_iedb_90 = len(combined_df[combined_df['IEDB_90%'] > 0.5])
        high_iedb_80 = len(combined_df[combined_df['IEDB_80%'] > 0.5])
        
        html_content += f"""
                <div class="metric">
                    <div class="metric-value">{total_epitopes}</div>
                    <div class="metric-label">Total Epitopes</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{total_genes}</div>
                    <div class="metric-label">Genes Analyzed</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{avg_combined_score:.3f}</div>
                    <div class="metric-label">Avg Combined Score</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{high_iedb_90}</div>
                    <div class="metric-label">High IEDB 90% (>50%)</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{high_iedb_80}</div>
                    <div class="metric-label">High IEDB 80% (>50%)</div>
                </div>
        """
    else:
        html_content += "<p><strong>No epitope data available for analysis.</strong></p>"
    
    html_content += """
            </div>
        </div>
    """
    
    # Top Candidates Table
    if not combined_df.empty:
        top_10 = combined_df.nlargest(10, 'Combined_Score')
        html_content += """
        <div class="section">
            <h2>🏆 Top 10 Epitope Candidates</h2>
        """
        html_content += top_10[['Gram_Type', 'Gene', 'Sequence', 'Length', 'Combined_Score', 
                               'IEDB_90%', 'IEDB_80%', 'Jalview_Similarity', 'BLOSUM62']].to_html(
                                   index=False, classes='highlight', table_id='top-candidates')
        html_content += "</div>"
    
    # Gram Type Comparison
    if not combined_df.empty:
        html_content += """
        <div class="section">
            <h2>⚖️ Gram-Positive vs Gram-Negative Comparison</h2>
        """
        
        if not pos_df.empty and not neg_df.empty:
            comparison_data = {
                'Metric': ['Total Epitopes', 'Avg Combined Score', 'Avg IEDB 90%', 'Avg IEDB 80%', 
                          'Avg Jalview Similarity', 'Avg BLOSUM62'],
                'Gram-Positive': [
                    len(pos_df),
                    f"{pos_df['Combined_Score'].mean():.3f}",
                    f"{pos_df['IEDB_90%'].mean():.3f}",
                    f"{pos_df['IEDB_80%'].mean():.3f}",
                    f"{pos_df['Jalview_Similarity'].mean():.3f}",
                    f"{pos_df['BLOSUM62'].mean():.3f}"
                ],
                'Gram-Negative': [
                    len(neg_df),
                    f"{neg_df['Combined_Score'].mean():.3f}",
                    f"{neg_df['IEDB_90%'].mean():.3f}",
                    f"{neg_df['IEDB_80%'].mean():.3f}",
                    f"{neg_df['Jalview_Similarity'].mean():.3f}",
                    f"{neg_df['BLOSUM62'].mean():.3f}"
                ]
            }
            comparison_df = pd.DataFrame(comparison_data)
            html_content += comparison_df.to_html(index=False)
        else:
            html_content += "<p>Insufficient data for comparison (need both Gram-positive and Gram-negative results).</p>"
        
        html_content += "</div>"
    
    # Complete Results Table
    if not combined_df.empty:
        html_content += """
        <div class="section">
            <h2>📋 Complete IEDB Validation Results</h2>
            <p>All epitopes ranked by combined score (current conservation methods + IEDB validation)</p>
        """
        html_content += combined_df.to_html(index=False, table_id='complete-results')
        html_content += "</div>"
    
    # Methodology
    html_content += """
        <div class="section">
            <h2>🔬 Methodology</h2>
            <div class="summary-box">
                <h3>Scoring Methods Used:</h3>
                <ul>
                    <li><strong>Current Conservation Methods:</strong>
                        <ul>
                            <li>Jalview % Similarity: Position-wise conservation using BLOSUM62-based similarity</li>
                            <li>BLOSUM62: Biochemical similarity conservation analysis</li>
                            <li>Identity: Exact amino acid identity conservation</li>
                        </ul>
                    </li>
                    <li><strong>IEDB Validation:</strong>
                        <ul>
                            <li>IEDB 100%: Exact sequence match conservancy</li>
                            <li>IEDB 90%: 90% identity threshold conservancy</li>
                            <li>IEDB 80%: 80% identity threshold conservancy</li>
                            <li>IEDB 70%: 70% identity threshold conservancy</li>
                        </ul>
                    </li>
                </ul>
                <h3>Combined Scoring Formula:</h3>
                <p><code>Combined Score = 0.3 × Jalview + 0.25 × BLOSUM62 + 0.2 × IEDB_90% + 0.15 × IEDB_80% + 0.1 × Identity</code></p>
            </div>
        </div>
    """
    
    # Analysis Metadata
    html_content += """
        <div class="section">
            <h2>📝 Analysis Metadata</h2>
            <div class="summary-box">
    """
    
    if pos_data:
        pos_metadata = pos_data.get('metadata', {})
        html_content += f"""
                <h3>Gram-Positive Analysis:</h3>
                <ul>
                    <li>Total genes analyzed: {pos_metadata.get('total_genes', 'N/A')}</li>
                    <li>Analysis timestamp: {pos_metadata.get('timestamp', 'N/A')}</li>
                </ul>
        """
    
    if neg_data:
        neg_metadata = neg_data.get('metadata', {})
        html_content += f"""
                <h3>Gram-Negative Analysis:</h3>
                <ul>
                    <li>Total genes analyzed: {neg_metadata.get('total_genes', 'N/A')}</li>
                    <li>Analysis timestamp: {neg_metadata.get('timestamp', 'N/A')}</li>
                </ul>
        """
    
    html_content += """
            </div>
        </div>
    """
    
    # Footer
    html_content += f"""
        <div class="section">
            <hr>
            <p><em>Report generated by TargSeek IEDB Validation Pipeline - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</em></p>
        </div>
    </body>
    </html>
    """
    
    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"IEDB combined report generated: {output_file}")
    
    # Print summary to console
    if not combined_df.empty:
        print(f"\n📊 Summary:")
        print(f"   Total epitopes: {len(combined_df)}")
        print(f"   Genes analyzed: {combined_df['Gene'].nunique()}")
        print(f"   Average combined score: {combined_df['Combined_Score'].mean():.3f}")
        print(f"   Top candidate: {top_10.iloc[0]['Gene']} - {top_10.iloc[0]['Sequence']} (score: {top_10.iloc[0]['Combined_Score']:.3f})")

def main():
    # Get parameters from Snakemake
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    
    # Input files
    positive_results = Path(snakemake.input.positive_results)
    negative_results = Path(snakemake.input.negative_results)
    
    # Output file
    output_file = Path(snakemake.output.combined_report)
    
    # Create output directory
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate report
    generate_html_report(positive_results, negative_results, output_file, analysis, paramset)

if __name__ == "__main__":
    main()