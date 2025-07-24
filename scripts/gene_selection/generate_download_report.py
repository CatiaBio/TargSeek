#!/usr/bin/env python3
"""
Generate Download Report
Creates comprehensive HTML report summarizing the data download pipeline results

This script now prioritizes protein names from UniProt data (stored in the summary TSV files)
over QuickGO annotations for more accurate and comprehensive protein descriptions.
"""

import pandas as pd
import json
import os
from pathlib import Path
import glob
from datetime import datetime
import yaml

def count_lines_in_file(filepath):
    """Count number of lines in a file, return 0 if file doesn't exist"""
    try:
        with open(filepath, 'r') as f:
            return sum(1 for line in f if line.strip())
    except FileNotFoundError:
        return 0

def load_annotations(annotations_file):
    """Load QuickGO annotations and create gene to protein name mapping"""
    try:
        with open(annotations_file, 'r') as f:
            annotations = json.load(f)
        
        gene_info = {}
        # Handle both list format and dict format with 'results' key
        entries = annotations if isinstance(annotations, list) else annotations.get('results', [])
        
        for entry in entries:
            gene_symbol = entry.get('geneProductId', 'Unknown')
            protein_name = entry.get('name', 'Unknown protein')
            go_id = entry.get('goId', '')
            
            # Extract gene name from UniProt ID if needed
            if gene_symbol.startswith('UniProtKB:'):
                # Try to get gene name from symbol field if available
                symbol = entry.get('symbol', '')
                if symbol:
                    gene_name = symbol
                else:
                    gene_name = gene_symbol.split(':')[1] if ':' in gene_symbol else gene_symbol
            else:
                gene_name = gene_symbol
            
            # Localization will be loaded from summary.tsv instead
            localization = 'Unknown'
            
            if gene_name not in gene_info:
                gene_info[gene_name] = {
                    'protein_name': protein_name,
                    'localization': localization,
                    'go_terms': [go_id] if go_id else []
                }
            else:
                # Add additional GO terms if gene already exists
                if go_id and go_id not in gene_info[gene_name]['go_terms']:
                    gene_info[gene_name]['go_terms'].append(go_id)
        
        return gene_info
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error loading annotations: {e}")
        return {}


def get_dataset_description(analysis):
    """Get dataset description from the first comment line of species file"""
    try:
        # Load config to get species file path
        with open('config/config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        
        species_file_path = config.get('species_files', {}).get(analysis)
        if not species_file_path:
            return f"Dataset: {analysis}"
        
        # Read first comment line from species file
        with open(species_file_path, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                if stripped_line.startswith('#'):
                    # Remove # and any extra spaces
                    description = stripped_line[1:].strip()
                    return description if description else f"Dataset: {analysis}"
                elif stripped_line:  # If we hit a non-comment line, stop
                    break
        
        return f"Dataset: {analysis}"
    except Exception as e:
        print(f"Error reading dataset description: {e}")
        return f"Dataset: {analysis}"

def load_go_descriptions():
    """Load GO term descriptions from the TSV file"""
    go_descriptions = {}
    try:
        go_file = Path("config/quickgo/go_id_descriptionl.tsv")
        if go_file.exists():
            df = pd.read_csv(go_file, sep='\t', comment='#')
            for _, row in df.iterrows():
                go_id = row['GeneOntology'].strip()
                description = row['Description'].strip()
                go_descriptions[go_id] = description
        return go_descriptions
    except Exception as e:
        print(f"Error loading GO descriptions: {e}")
        return {}

def get_go_parameters_info(paramset):
    """Get GO terms and their descriptions from the parameter JSON file"""
    try:
        params_file = f"config/quickgo/{paramset}.json"
        with open(params_file, 'r') as f:
            params = json.load(f)
        
        go_ids = params.get('goId', '').split(',')
        go_descriptions = load_go_descriptions()
        
        go_info = []
        for go_id in go_ids:
            go_id = go_id.strip()
            if go_id:
                description = go_descriptions.get(go_id, 'Description not available')
                go_info.append({'id': go_id, 'description': description})
        
        return go_info
    except Exception as e:
        print(f"Error loading GO parameters: {e}")
        return []

def generate_go_parameters_list(go_parameters):
    """Generate HTML list items for GO parameters"""
    if not go_parameters:
        return "<li>No GO parameters available</li>"
    
    items = []
    for go_param in go_parameters:
        items.append(f"<li><strong>{go_param['id']}</strong> - {go_param['description']}</li>")
    
    return "\n                    ".join(items)

def check_genes_have_3d_structures(metadata_file_path):
    """Check which genes have 3D structures in the unified metadata file"""
    metadata_file = Path(metadata_file_path)
    
    if not metadata_file.exists():
        print(f"3D structure metadata file not found: {metadata_file}")
        return set()
    
    try:
        # Read metadata TSV
        df = pd.read_csv(metadata_file, sep='\t')
        
        # Handle empty file
        if df.empty:
            print(f"No 3D structure data found in {metadata_file}")
            return set()
        
        # Get unique gene names that have structures
        genes_with_structures = set(df['gene_name'].unique())
        
        print(f"Found 3D structures for {len(genes_with_structures)} genes in {metadata_file}")
        return genes_with_structures
        
    except FileNotFoundError:
        print(f"3D structure metadata file not found: {metadata_file}")
        return set()
    except pd.errors.EmptyDataError:
        print(f"3D structure metadata file is empty: {metadata_file}")
        return set()
    except Exception as e:
        print(f"Error reading 3D structure metadata from {metadata_file}: {e}")
        return set()

def generate_html_report(data, output_file):
    """Generate HTML report"""
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>TargSeek Download Report - {data['analysis']}_{data['paramset']}</title>
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
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
            }}
            h2 {{
                color: #34495e;
                border-left: 4px solid #3498db;
                padding-left: 15px;
                margin-top: 30px;
            }}
            h3 {{
                color: #7f8c8d;
                margin-top: 25px;
            }}
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin: 20px 0;
            }}
            .stat-card {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 20px;
                border-radius: 8px;
                text-align: center;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            }}
            .stat-number {{
                font-size: 2em;
                font-weight: bold;
                display: block;
            }}
            .stat-label {{
                font-size: 0.9em;
                opacity: 0.9;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                background: white;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            th, td {{
                padding: 12px;
                text-align: left;
                border-bottom: 1px solid #ddd;
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
            .gram-positive {{
                border-left: 4px solid #27ae60;
            }}
            .gram-negative {{
                border-left: 4px solid #e74c3c;
            }}
            .coverage-high {{
                background-color: #d4edda;
                color: #155724;
            }}
            .coverage-medium {{
                background-color: #fff3cd;
                color: #856404;
            }}
            .coverage-low {{
                background-color: #f8d7da;
                color: #721c24;
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
                font-style: italic;
                margin-top: 30px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>TargSeek Download Report</h1>
            <div class="info-section">
                <strong>Data:</strong> {data['dataset_description']}<br>
                <strong>Gene Ontology parameters:</strong><br>
                <ul style="margin-left: 20px; margin-top: 5px;">
                    {generate_go_parameters_list(data['go_parameters'])}
                </ul>
                <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </div>

            <h2>📊 Species Classification Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <span class="stat-number">{data['species_stats']['total_species']}</span>
                    <span class="stat-label">Total Species Identified</span>
                </div>
                <div class="stat-card">
                    <span class="stat-number">{data['species_stats']['gram_positive']}</span>
                    <span class="stat-label">Gram-Positive Species</span>
                </div>
                <div class="stat-card">
                    <span class="stat-number">{data['species_stats']['gram_negative']}</span>
                    <span class="stat-label">Gram-Negative Species</span>
                </div>
                <div class="stat-card">
                    <span class="stat-number">{data['species_stats']['classification_rate']:.1f}%</span>
                    <span class="stat-label">Classification Success Rate</span>
                </div>
            </div>

            <h2>🧬 Gene Filtering Pipeline</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <span class="stat-number">{data['gene_stats']['initial_genes']}</span>
                    <span class="stat-label">Initial Genes from QuickGO</span>
                </div>
                <div class="stat-card">
                    <span class="stat-number">{data['gene_stats']['surface_accessible']}</span>
                    <span class="stat-label">Surface-Accessible Genes</span>
                </div>
                <div class="stat-card">
                    <span class="stat-number">{data['gene_stats']['final_filtered']}</span>
                    <span class="stat-label">Final Filtered Genes</span>
                </div>
                <div class="stat-card">
                    <span class="stat-number">{data['gene_stats']['filter_rate']:.1f}%</span>
                    <span class="stat-label">Gene Retention Rate</span>
                </div>
            </div>

            <h2>🦠 Gram-Positive Proteins</h2>
            <div class="gram-positive">
                <div class="info-section">
                    <strong>Total genes found:</strong> {data['gram_positive_stats']['total_genes']}<br>
                    <strong>Genes with ≥50% coverage:</strong> {data['gram_positive_stats']['above_threshold']}<br>
                </div>
                <h3>Top 5 Selected Proteins</h3>
                {generate_protein_table(data['gram_positive_proteins'])}
            </div>

            <h2>🦠 Gram-Negative Proteins</h2>
            <div class="gram-negative">
                <div class="info-section">
                    <strong>Total genes found:</strong> {data['gram_negative_stats']['total_genes']}<br>
                    <strong>Genes with ≥50% coverage:</strong> {data['gram_negative_stats']['above_threshold']}<br>
                </div>
                <h3>Top 5 Selected Proteins</h3>
                {generate_protein_table(data['gram_negative_proteins'])}
            </div>

            <h2>🎯 Final Selected Proteins Based on Coverage Thresholds</h2>
            {generate_selected_proteins_section(data['coverage_thresholds'], data['selected_proteins_summary'])}

            <div class="timestamp">
                Report generated on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')}
            </div>
        </div>
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)

def generate_protein_table(proteins_data):
    """Generate HTML table for proteins"""
    if not proteins_data:
        return "<p>No proteins found for this Gram classification.</p>"
    
    table_html = """
    <table>
        <thead>
            <tr>
                <th>Gene Name</th>
                <th>Protein Name</th>
                <th>Location</th>
                <th>Species Coverage (%)</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for protein in proteins_data:
        coverage_class = get_coverage_class(protein['coverage'])
        table_html += f"""
            <tr>
                <td><strong>{protein['gene_name']}</strong></td>
                <td>{protein['protein_name']}</td>
                <td>{protein['localization']}</td>
                <td class="{coverage_class}">{protein['coverage']:.1f}%</td>
            </tr>
        """
    
    table_html += """
        </tbody>
    </table>
    """
    return table_html

def get_coverage_class(coverage):
    """Determine coverage class based on percentage"""
    if coverage >= 70:
        return "coverage-high"
    elif coverage >= 50:
        return "coverage-medium"
    else:
        return "coverage-low"


def generate_selected_proteins_section(coverage_thresholds, selected_summary):
    """Generate section showing final selected proteins based on thresholds"""
    html = f"""
    
    <div class="stats-grid">
        <div class="stat-card">
            <span class="stat-number">{selected_summary['total_selected']}</span>
            <span class="stat-label">Total Proteins Selected</span>
        </div>
        <div class="stat-card">
            <span class="stat-number">{selected_summary['selected_positive']}</span>
            <span class="stat-label">Gram-Positive Selected</span>
        </div>
        <div class="stat-card">
            <span class="stat-number">{selected_summary['selected_negative']}</span>
            <span class="stat-label">Gram-Negative Selected</span>
        </div>
        <div class="stat-card">
            <span class="stat-number">{selected_summary['avg_coverage']:.1f}%</span>
            <span class="stat-label">Average Coverage</span>
        </div>
    </div>
    
    <h3>Final Selected Proteins Summary</h3>
    {generate_final_selection_table(selected_summary['selected_proteins'])}
    """
    return html

def generate_final_selection_table(selected_proteins):
    """Generate table showing final selected proteins across both Gram groups"""
    if not selected_proteins:
        return "<p>No proteins were selected based on the coverage thresholds.</p>"
    
    table_html = """
    <table>
        <thead>
            <tr>
                <th>Gram Type</th>
                <th>Gene Name</th>
                <th>Protein Name</th>
                <th>Coverage (%)</th>
                <th>3D Structures</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for protein in selected_proteins:
        gram_class = "gram-positive" if protein['gram_type'] == 'positive' else "gram-negative"
        coverage_class = get_coverage_class(protein['coverage'])
        structures_info = "✓" if protein['structures'] > 0 else "✗"
        structures_color = "color: #27ae60;" if protein['structures'] > 0 else "color: #e74c3c;"
        
        table_html += f"""
            <tr class="{gram_class}">
                <td><strong>Gram-{protein['gram_type'].title()}</strong></td>
                <td><strong>{protein['gene_name']}</strong></td>
                <td>{protein['protein_name']}</td>
                <td class="{coverage_class}">{protein['coverage']:.1f}%</td>
                <td style="{structures_color} font-weight: bold; text-align: center;">{structures_info}</td>
            </tr>
        """
    
    table_html += """
        </tbody>
    </table>
    """
    return table_html

def main():
    # Get parameters from snakemake
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    
    # Get dataset description and GO parameters info
    dataset_description = get_dataset_description(analysis)
    go_parameters = get_go_parameters_info(paramset)
    
    # Initialize data dictionary
    data = {
        'analysis': analysis,
        'paramset': paramset,
        'dataset_description': dataset_description,
        'go_parameters': go_parameters,
        'species_stats': {},
        'gene_stats': {},
        'gram_positive_proteins': [],
        'gram_negative_proteins': [],
        'genes_with_structures': set(),
        'coverage_thresholds': {},
        'selected_proteins_summary': {}
    }
    
    # Load species classification stats
    try:
        gram_classification = pd.read_csv(snakemake.input.gram_classification, sep='\t')
        total_species = len(gram_classification)
        gram_positive_count = count_lines_in_file(snakemake.input.gram_positive)
        gram_negative_count = count_lines_in_file(snakemake.input.gram_negative)
        
        data['species_stats'] = {
            'total_species': total_species,
            'gram_positive': gram_positive_count,
            'gram_negative': gram_negative_count,
            'classification_rate': ((gram_positive_count + gram_negative_count) / total_species * 100) if total_species > 0 else 0
        }
    except Exception as e:
        print(f"Error loading species stats: {e}")
        data['species_stats'] = {'total_species': 0, 'gram_positive': 0, 'gram_negative': 0, 'classification_rate': 0}
    
    # Load gene filtering stats
    initial_genes = count_lines_in_file(snakemake.input.initial_genes)
    surface_genes = count_lines_in_file(snakemake.input.surface_genes)
    filtered_genes = count_lines_in_file(snakemake.input.filtered_genes)
    
    data['gene_stats'] = {
        'initial_genes': initial_genes,
        'surface_accessible': surface_genes,
        'final_filtered': filtered_genes,
        'filter_rate': (filtered_genes / initial_genes * 100) if initial_genes > 0 else 0
    }
    
    # Load gene annotations
    gene_info = load_annotations(snakemake.input.annotations)
    
    # Load coverage data
    try:
        coverage_df = pd.read_csv(snakemake.input.coverage, sep='\t')
    except Exception as e:
        print(f"Error loading coverage data: {e}")
        coverage_df = pd.DataFrame()
    
    # Load summary data to get all proteins information
    try:
        summary_df = pd.read_csv(snakemake.input.summary, sep='\t')
        print(f"Loaded summary data with {len(summary_df)} proteins")
    except Exception as e:
        print(f"Error loading summary data: {e}")
        summary_df = pd.DataFrame()
    
    # Process Gram-positive proteins with statistics
    gram_positive_stats = {'total_genes': 0, 'above_threshold': 0, 'selected_count': 0}
    try:
        # Get total genes for Gram-positive from coverage data
        if not coverage_df.empty:
            gram_pos_coverage = coverage_df[coverage_df['gram'] == 'positive']
            gram_positive_stats['total_genes'] = len(gram_pos_coverage)
            
            # Count those above 50% threshold
            above_threshold = gram_pos_coverage[gram_pos_coverage['coverage_percentage'] >= 50]
            gram_positive_stats['above_threshold'] = len(above_threshold)
        
        # Get Gram-positive proteins from summary
        if not summary_df.empty:
            gram_pos_proteins = summary_df[summary_df['gram'] == 'positive']
            gram_positive_stats['selected_count'] = len(gram_pos_proteins)
            
            # Process top 5 proteins for detailed display
            top_proteins = gram_pos_proteins.head(5)
            for _, row in top_proteins.iterrows():
                data['gram_positive_proteins'].append({
                    'gene_name': row['gene_name'],
                    'protein_name': row['protein_name'],
                    'localization': row['location'],
                    'coverage': row['coverage_percentage']
                })
        
        data['gram_positive_stats'] = gram_positive_stats
    except Exception as e:
        print(f"Error processing Gram-positive proteins: {e}")
        data['gram_positive_stats'] = {'total_genes': 0, 'above_threshold': 0, 'selected_count': 0}
    
    # Process Gram-negative proteins with statistics
    gram_negative_stats = {'total_genes': 0, 'above_threshold': 0, 'selected_count': 0}
    try:
        # Get total genes for Gram-negative from coverage data
        if not coverage_df.empty:
            gram_neg_coverage = coverage_df[coverage_df['gram'] == 'negative']
            gram_negative_stats['total_genes'] = len(gram_neg_coverage)
            
            # Count those above 50% threshold
            above_threshold = gram_neg_coverage[gram_neg_coverage['coverage_percentage'] >= 50]
            gram_negative_stats['above_threshold'] = len(above_threshold)
        
        # Get Gram-negative proteins from summary
        if not summary_df.empty:
            gram_neg_proteins = summary_df[summary_df['gram'] == 'negative']
            gram_negative_stats['selected_count'] = len(gram_neg_proteins)
            
            # Process top 5 proteins for detailed display
            top_proteins = gram_neg_proteins.head(5)
            for _, row in top_proteins.iterrows():
                data['gram_negative_proteins'].append({
                    'gene_name': row['gene_name'],
                    'protein_name': row['protein_name'],
                    'localization': row['location'],
                    'coverage': row['coverage_percentage']
                })
        
        data['gram_negative_stats'] = gram_negative_stats
    except Exception as e:
        print(f"Error processing Gram-negative proteins: {e}")
        data['gram_negative_stats'] = {'total_genes': 0, 'above_threshold': 0, 'selected_count': 0}
    
    # Check which genes have 3D structures from unified metadata file
    genes_with_structures = check_genes_have_3d_structures(snakemake.input.metadata)
    data['genes_with_structures'] = genes_with_structures
    
    # Load coverage thresholds from config (need to read config file)
    try:
        import yaml
        with open('config/config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        
        data['coverage_thresholds'] = {
            'positive': config.get('gram_thresholds', {}).get('positive', 50),
            'negative': config.get('gram_thresholds', {}).get('negative', 50)
        }
        
    except Exception as e:
        print(f"Could not load config for thresholds: {e}")
        data['coverage_thresholds'] = {'positive': 50, 'negative': 50}
    
    # Process selected proteins summary
    selected_proteins = []
    total_coverage = 0
    
    # Process selected proteins from summary data
    try:
        if not summary_df.empty:
            for idx, (_, row) in enumerate(summary_df.iterrows(), 1):
                gene_name = row['gene_name']
                gram_type = row['gram']
                structures_count = 0
                
                # Check if gene has 3D structures
                has_structures = gene_name in genes_with_structures
                structures_count = 1 if has_structures else 0
                
                selected_proteins.append({
                    'gram_type': gram_type,
                    'gene_name': gene_name,
                    'protein_name': row['protein_name'],
                    'coverage': row['coverage_percentage'],
                    'structures': structures_count
                })
                total_coverage += row['coverage_percentage']
    except Exception as e:
        print(f"Error processing selected proteins from summary: {e}")
    
    # Calculate summary statistics
    selected_positive = len([p for p in selected_proteins if p['gram_type'] == 'positive'])
    selected_negative = len([p for p in selected_proteins if p['gram_type'] == 'negative'])
    avg_coverage = total_coverage / len(selected_proteins) if selected_proteins else 0
    
    data['selected_proteins_summary'] = {
        'selected_positive': selected_positive,
        'selected_negative': selected_negative,
        'total_selected': len(selected_proteins),
        'avg_coverage': avg_coverage,
        'selected_proteins': selected_proteins
    }
    
    # Generate HTML report
    os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)
    generate_html_report(data, snakemake.output[0])
    
    print(f"Download report generated: {snakemake.output[0]}")

if __name__ == "__main__":
    main()