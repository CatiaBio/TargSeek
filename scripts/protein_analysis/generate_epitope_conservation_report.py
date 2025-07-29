#!/usr/bin/env python3
"""
Generate Epitope Conservation Analysis Report
============================================

This script generates a comprehensive HTML report from epitope conservation 
analysis JSON results with the same styling as the download report.

Usage:
    python scripts/protein_analysis/generate_epitope_conservation_report.py
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime
import logging

# Configure logging for Snakemake compatibility
if 'snakemake' in globals():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
else:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def collect_analysis_results(output_dir):
    """Collect all analysis results from JSON files"""
    results = {}
    
    for json_file in Path(output_dir).glob("*.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                
            # Extract gene name from filename
            gene_structure = json_file.stem.replace('_conservation', '')
            gene = gene_structure.split('_')[0]
            
            results[gene] = data
            
        except Exception as e:
            logging.warning(f"Error loading {json_file}: {e}")
    
    return results

def generate_comprehensive_report(results, output_file, analysis, paramset):
    """Generate comprehensive HTML report with download report styling"""
    
    # Calculate summary statistics
    total_genes = len(results)
    total_epitopes = sum(len(data.get('epitope_analyses', [])) for data in results.values())
    
    if total_genes > 0:
        mean_conservation = sum(
            data.get('summary_statistics', {}).get('conservation_statistics', {}).get('mean_conservation', 0)
            for data in results.values()
        ) / total_genes
        
        mean_variants = sum(
            data.get('summary_statistics', {}).get('diversity_statistics', {}).get('mean_variants', 0)
            for data in results.values()
        ) / total_genes
    else:
        mean_conservation = 0
        mean_variants = 0
    
    # HTML template with download report styling
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Epitope Conservation Analysis Report - {analysis}_{paramset}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            margin-bottom: 30px;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 5px;
            margin-top: 30px;
        }}
        h3 {{
            color: #7f8c8d;
            margin-top: 25px;
        }}
        .info-section {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 15px 0;
        }}
        .timestamp {{
            text-align: center;
            color: #7f8c8d;
            font-size: 0.9em;
            margin-top: 30px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            table-layout: fixed;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
            word-wrap: break-word;
            overflow-wrap: break-word;
            vertical-align: top;
        }}
        th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #2c3e50;
            position: sticky;
            top: 0;
        }}
        tr:hover {{
            background-color: #f8f9fa;
        }}
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        .sequence-cell {{
            max-width: 300px;
            word-break: break-word;
            overflow-wrap: break-word;
            hyphens: none;
            font-size: 0.9em;
            line-height: 1.3;
            padding: 8px 12px;
            vertical-align: top;
        }}
        .position-cell {{
            white-space: nowrap;
            min-width: 80px;
        }}
        .conservation-cell {{
            white-space: nowrap;
            min-width: 80px;
            text-align: left;
        }}
        .numeric-cell {{
            white-space: nowrap;
            text-align: left;
            font-size: 0.9em;
        }}
        .conservation-high {{
            background-color: #d4edda;
            color: #155724;
        }}
        .conservation-medium {{
            background-color: #fff3cd;
            color: #856404;
        }}
        .conservation-low {{
            background-color: #f8d7da;
            color: #721c24;
        }}
        .gram-positive {{
            border-left: 4px solid #27ae60;
        }}
        .gram-negative {{
            border-left: 4px solid #e74c3c;
        }}
        .summary-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
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
        .consensus {{
            font-family: 'Courier New', monospace;
            font-weight: bold;
            background-color: #f8f9fa;
            padding: 2px 4px;
            border-radius: 3px;
            letter-spacing: 0.5px;
        }}
        .sequence-text {{
            font-family: 'Courier New', monospace;
            letter-spacing: 0.5px;
            word-spacing: 2px;
        }}
        /* Landscape mode optimization - preferred for long sequences */
        @media screen and (orientation: landscape) {{
            .container {{
                max-width: 95%;
                padding: 20px;
            }}
            .sequence-cell {{
                max-width: 500px;
                word-break: break-word;
                overflow-wrap: break-word;
                hyphens: none;
                font-size: 0.95em;
                line-height: 1.3;
                white-space: normal;
            }}
            .table-container {{
                overflow-x: auto;
                margin: 15px 0;
            }}
            table {{
                min-width: 1200px;
            }}
            th, td {{
                padding: 10px 8px;
            }}
            .summary-stats {{
                grid-template-columns: repeat(4, 1fr);
            }}
        }}
        
        /* Portrait mode - use word breaking for mobile */
        @media screen and (orientation: portrait) and (max-width: 768px) {{
            .table-container {{
                font-size: 0.8em;
            }}
            .sequence-cell {{
                max-width: 150px;
                font-size: 0.8em;
            }}
            th, td {{
                padding: 8px 4px;
            }}
        }}
        
        /* Large landscape screens - optimal for long sequences */
        @media screen and (orientation: landscape) and (min-width: 1024px) {{
            .sequence-cell {{
                max-width: 800px;
                font-size: 1em;
                word-break: break-word;
                overflow-wrap: break-word;
                line-height: 1.3;
            }}
            table {{
                min-width: 1400px;
            }}
            th, td {{
                padding: 12px 10px;
            }}
        }}
        
        /* Print styles - force landscape orientation */
        @media print {{
            @page {{
                size: landscape;
                margin: 0.5in;
            }}
            body {{
                font-size: 10px;
            }}
            .container {{
                max-width: 100%;
                padding: 10px;
                box-shadow: none;
            }}
            .sequence-cell {{
                max-width: 500px;
                word-break: normal;
                white-space: nowrap;
                font-size: 9px;
            }}
            .stat-card {{
                break-inside: avoid;
            }}
            .orientation-notice {{
                display: none;
            }}
            table {{
                break-inside: avoid;
            }}
        }}
        .epitope-summary {{
            margin: 10px 0;
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 5px;
        }}
        .epitope-header {{
            font-weight: bold;
            margin-bottom: 5px;
        }}
        .epitope-details {{
            font-size: 0.9em;
            color: #666;
        }}
        .orientation-notice {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            margin: 20px 0;
            font-weight: 500;
        }}
        .orientation-notice.landscape {{
            background: linear-gradient(135deg, #27ae60 0%, #2ecc71 100%);
        }}
        .visualization-section {{
            margin: 20px 0;
            padding: 15px;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            background: #fafafa;
        }}
        .visualization-container {{
            margin: 15px 0;
            padding: 15px;
            border: 1px solid #ddd;
            border-radius: 5px;
            background: white;
        }}
        .image-container {{
            text-align: center;
            margin: 10px 0;
        }}
        .legend-info {{
            margin-top: 10px;
            font-size: 0.9em;
            background: #f8f9fa;
            padding: 10px;
            border-radius: 3px;
            border-left: 3px solid #3498db;
        }}
        .legend-info pre {{
            font-family: 'Courier New', monospace;
            font-size: 0.8em;
            margin: 5px 0;
            white-space: pre-wrap;
            background: white;
            padding: 8px;
            border-radius: 3px;
            border: 1px solid #e0e0e0;
        }}
        @media screen and (orientation: landscape) {{
            .portrait-notice {{
                display: none;
            }}
        }}
        @media screen and (orientation: portrait) {{
            .landscape-notice {{
                display: none;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Epitope Conservation Analysis Report</h1>
        
        <div class="orientation-notice portrait-notice">
            📱 <strong>Tip:</strong> For better viewing of long epitope sequences, try rotating your device to landscape mode or viewing on a wider screen.
        </div>
        
        
        <div class="info-section">
            <h2>📊 Analysis Overview</h2>
            <p><strong>Analysis:</strong> {analysis}_{paramset}</p>
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Pipeline:</strong> TargSeek Epitope Conservation Analysis</p>
        </div>
        
        <div class="summary-stats">
            <div class="stat-card">
                <div class="stat-number">{total_genes}</div>
                <div class="stat-label">Genes Analyzed</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{total_epitopes}</div>
                <div class="stat-label">Total Epitopes</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{mean_conservation:.3f}</div>
                <div class="stat-label">Mean Conservation</div>
            </div>
            <div class="stat-card">
                <div class="stat-number">{mean_variants:.1f}</div>
                <div class="stat-label">Mean Variants</div>
            </div>
        </div>
        
        <h2>🎯 Gene-Level Summary</h2>
        
        <div class="table-container">
            <table>
                <thead>
                    <tr>
                        <th style="width: 10%;">Gene</th>
                        <th style="width: 15%;">Structure</th>
                        <th style="width: 10%;">Epitopes</th>
                        <th style="width: 15%;">MSA Sequences</th>
                        <th style="width: 15%;">Mean Conservation</th>
                        <th style="width: 15%;">Conservation Range</th>
                        <th style="width: 10%;">Mean Variants</th>
                        <th style="width: 10%;">Conservation Class</th>
                    </tr>
                </thead>
                <tbody>
"""
    
    # Add gene-level data
    for gene, data in sorted(results.items()):
        summary_stats = data.get('summary_statistics', {})
        conservation_stats = summary_stats.get('conservation_statistics', {})
        diversity_stats = summary_stats.get('diversity_statistics', {})
        
        mean_cons = conservation_stats.get('mean_conservation', 0)
        min_cons = conservation_stats.get('min_conservation', 0)
        max_cons = conservation_stats.get('max_conservation', 0)
        num_epitopes = summary_stats.get('total_epitopes', 0)
        num_sequences = summary_stats.get('total_msa_sequences', 0)
        mean_vars = diversity_stats.get('mean_variants', 0)
        structure_id = summary_stats.get('structure_id', 'Unknown')
        
        # Determine conservation class
        if mean_cons >= 0.7:
            cons_class = "conservation-high"
            cons_label = "High"
        elif mean_cons >= 0.4:
            cons_class = "conservation-medium"
            cons_label = "Medium"
        else:
            cons_class = "conservation-low"
            cons_label = "Low"
        
        html_content += f"""
                <tr class="{cons_class}">
                    <td><strong>{gene}</strong></td>
                    <td>{structure_id}</td>
                    <td>{num_epitopes}</td>
                    <td>{num_sequences}</td>
                    <td>{mean_cons:.3f}</td>
                    <td>{min_cons:.3f} - {max_cons:.3f}</td>
                    <td>{mean_vars:.1f}</td>
                    <td>{cons_label}</td>
                </tr>"""
    
    html_content += """
            </tbody>
        </table>
        </div>
        
        <h2>🔍 Detailed Epitope Analysis</h2>
"""
    
    # Add detailed epitope information
    for gene, data in sorted(results.items()):
        epitope_analyses = data.get('epitope_analyses', [])
        if not epitope_analyses:
            continue
            
        structure_id = data.get('summary_statistics', {}).get('structure_id', 'Unknown')
        
        html_content += f"""
        <h3>{gene} ({structure_id})</h3>
        
        <div class="table-container">
            <table>
                <thead>
                    <tr>
                        <th style="width: 5%;">Epitope</th>
                        <th style="width: 35%;">Sequence</th>
                        <th style="width: 8%;">Position</th>
                        <th style="width: 5%;">Length</th>
                        <th style="width: 7%;">Conservation</th>
                        <th style="width: 5%;">Variants</th>
                        <th style="width: 17%;">Consensus 25%</th>
                        <th style="width: 18%;">Consensus 50%</th>
                    </tr>
                </thead>
                <tbody>
"""
        
        for epitope in epitope_analyses:
            epitope_info = epitope.get('epitope_info', {})
            conservation = epitope.get('conservation_scores', {})
            consensus = epitope.get('consensus_sequences', {})
            
            mean_cons = conservation.get('mean_conservation', 0)
            
            # Conservation class
            if mean_cons >= 0.7:
                cons_class = "conservation-high"
            elif mean_cons >= 0.4:
                cons_class = "conservation-medium"
            else:
                cons_class = "conservation-low"
            
            html_content += f"""
                <tr class="{cons_class}">
                    <td class="numeric-cell">{epitope_info.get('number', '')}</td>
                    <td class="sequence-cell sequence-text">{epitope_info.get('peptide', '')}</td>
                    <td class="position-cell numeric-cell">{epitope_info.get('start', '')}-{epitope_info.get('end', '')}</td>
                    <td class="numeric-cell">{epitope_info.get('length', '')}</td>
                    <td class="conservation-cell numeric-cell">{mean_cons:.3f}</td>
                    <td class="numeric-cell">{epitope.get('variant_count', 0)}</td>
                    <td class="consensus sequence-cell sequence-text">{consensus.get('consensus_25', 'N/A')}</td>
                    <td class="consensus sequence-cell sequence-text">{consensus.get('consensus_50', 'N/A')}</td>
                </tr>"""
        
        html_content += """
            </tbody>
        </table>
        </div>
"""
        
        # Add 3D visualization section for this gene
        html_content += f"""
        <h4>🎨 3D Epitope Visualizations for {gene}</h4>
        <div class="visualization-section">
"""
        
        # Look for visualization images
        visualization_dir = Path(f"results/{analysis}_{paramset}/protein_analysis/sequences_with_structure/epitope_predictions_bepipred/3d_visualizations")
        gene_images = list(visualization_dir.glob(f"{gene}_*_epitopes.png"))
        
        
        if gene_images:
            for img_path in sorted(gene_images):
                # Extract structure ID from filename
                img_name = img_path.name
                structure_from_img = img_name.replace(f"{gene}_", "").replace("_epitopes.png", "")
                
                # Get relative path for HTML (from report location to images)
                relative_img_path = f"protein_analysis/sequences_with_structure/epitope_predictions_bepipred/3d_visualizations/{img_name}"
                
                # Check for corresponding legend file
                legend_file = img_path.with_name(img_name.replace("_epitopes.png", "_legend.txt"))
                
                html_content += f"""
            <div class="visualization-container">
                <h5>Structure: {structure_from_img}</h5>
                <div class="image-container">
                    <img src="{relative_img_path}" alt="3D Epitope Visualization for {gene} ({structure_from_img})" 
                         style="max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
                </div>"""
                
                # Add legend information if available
                if legend_file.exists():
                    try:
                        with open(legend_file, 'r') as f:
                            legend_content = f.read()
                        
                        # Extract epitope details from legend
                        if "Epitope Details:" in legend_content:
                            details_section = legend_content.split("Epitope Details:")[1].split("Visualization Notes:")[0]
                            html_content += f"""
                <div class="legend-info" style="margin-top: 10px; font-size: 0.9em; background: #f8f9fa; padding: 10px; border-radius: 3px;">
                    <strong>Epitope Color Legend:</strong>
                    <pre style="font-family: 'Courier New', monospace; font-size: 0.8em; margin: 5px 0; white-space: pre-wrap;">{details_section.strip()}</pre>
                </div>"""
                    except Exception as e:
                        pass
                
                html_content += """
            </div>
            <br>"""
        else:
            html_content += f"""
            <p style="color: #666; font-style: italic;">No 3D visualizations available for {gene}</p>"""
        
        html_content += """
        </div>
"""
    
    html_content += f"""
        <div class="info-section">
            <h2>📝 Analysis Notes</h2>
            <ul>
                <li><strong>Conservation Score:</strong> Based on Shannon entropy (0-1 scale, 1 = perfect conservation)</li>
                <li><strong>High Conservation:</strong> ≥0.7 (suitable for broad-spectrum targeting)</li>
                <li><strong>Medium Conservation:</strong> 0.4-0.7 (genus or family-specific)</li>
                <li><strong>Low Conservation:</strong> &lt;0.4 (species-specific or highly variable)</li>
                <li><strong>Consensus Sequences:</strong> Amino acids appearing above threshold frequency</li>
                <li><strong>X in Consensus:</strong> No amino acid reaches the threshold (ambiguous position)</li>
                <li><strong>3D Visualizations:</strong> PyMOL-generated images showing epitopes as colored cartoon regions on protein structures</li>
                <li><strong>Epitope Colors:</strong> Each epitope is assigned a unique color (red, blue, green, yellow, etc.) for easy identification</li>
            </ul>
        </div>
        
        <div class="timestamp">
            🚀 Generated by TargSeek Protein Analysis Pipeline<br>
            Report created on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')}
        </div>
    </div>
</body>
</html>
"""
    
    # Write HTML report
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    logging.info(f"Comprehensive report written to: {output_file}")

def main():
    """Main function for both Snakemake and command line usage"""
    
    if 'snakemake' in globals():
        # Running from Snakemake
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        output_dir = Path(snakemake.params.conservation_analysis_dir)
        output_file = snakemake.output.conservation_report
        
        logging.info(f"Generating epitope conservation report for {analysis}_{paramset}")
        
    else:
        # Command line usage (fallback)
        analysis = "analysis1"
        paramset = "params1"
        output_dir = Path("epitope_conservation_analysis")
        output_file = f"{analysis}_{paramset}_epitope_conservation_summary.html"
        
        logging.info(f"Generating epitope conservation report (command line mode)")
    
    # Collect analysis results from JSON files
    results = collect_analysis_results(output_dir)
    
    if not results:
        logging.error(f"No conservation analysis results found in {output_dir}")
        return
    
    logging.info(f"Found results for {len(results)} genes")
    
    # Generate comprehensive HTML report
    generate_comprehensive_report(results, output_file, analysis, paramset)
    
    logging.info(f"Report generation complete!")
    logging.info(f"HTML report: {output_file}")

if __name__ == "__main__":
    main()