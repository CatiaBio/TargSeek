#!/usr/bin/env python3
"""
3D Structure Metadata Extractor
================================

This script extracts metadata from downloaded 3D structures and creates a comprehensive
TSV file with structure information per gene including:
- Gene name and structure IDs
- Structure type (experimental vs computed)
- Organism information
- Quality scores (resolution for PDB, confidence for AlphaFold)
- Additional metadata

Usage:
- Run via Snakemake: Automatically processes all structure directories
- Standalone: python extract_3d_structure_metadata.py
"""

import pandas as pd
import json
from pathlib import Path
import logging
import glob
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extract_metadata_for_analysis(structures_base_dir, analysis, paramset):
    """
    Extract metadata for all genes and structures in the analysis
    
    Args:
        structures_base_dir (Path): Base directory containing 3D structures
        analysis (str): Analysis identifier (e.g., analysis_1)
        paramset (str): Parameter set identifier (e.g., params_1)
    
    Returns:
        list: List of dictionaries with metadata for each structure
    """
    
    metadata_records = []
    
    # Find all gene directories in the structures directory
    gene_dirs = [d for d in structures_base_dir.iterdir() if d.is_dir()]
    
    logging.info(f"Processing {len(gene_dirs)} gene directories for {analysis}_{paramset}")
    
    for gene_dir in gene_dirs:
        gene_name = gene_dir.name
        metadata_file = gene_dir / "metadata.json"
        
        # Skip if no metadata file exists
        if not metadata_file.exists():
            logging.debug(f"No metadata file found for gene {gene_name}")
            continue
        
        try:
            # Load metadata
            with open(metadata_file, 'r') as f:
                metadata_list = json.load(f)
            
            # Process each structure's metadata
            for metadata in metadata_list:
                record = {
                    'analysis': analysis,
                    'paramset': paramset,
                    'gene_name': gene_name,
                    'structure_id': metadata.get('id', 'Unknown'),
                    'structure_type': metadata.get('type', 'Unknown'),
                    'database': metadata.get('database', 'Unknown'),
                    'organism': metadata.get('organism', 'Unknown'),
                    'protein_name': metadata.get('protein_name', metadata.get('title', 'Unknown')),
                    'experimental_method': metadata.get('method', 'Unknown'),
                    'resolution': metadata.get('resolution', 'N/A'),
                    'model_confidence': metadata.get('model_confidence', 'N/A'),
                    'deposit_date': metadata.get('deposit_date', metadata.get('created_date', 'Unknown')),
                    'journal': metadata.get('journal', 'N/A'),
                    'authors': ', '.join(metadata.get('authors', [])) if metadata.get('authors') else 'Unknown',
                    'uniprot_accession': metadata.get('uniprot_accession', 'N/A'),
                    'model_version': metadata.get('model_version', 'N/A'),
                    'sequence_length': metadata.get('sequence_length', 'N/A')
                }
                
                # Add quality score based on structure type
                if metadata.get('type') == 'experimental_structure':
                    # For experimental structures, use resolution as quality score (lower is better)
                    resolution = metadata.get('resolution')
                    if resolution and resolution != 'N/A' and resolution != '':
                        try:
                            res_float = float(resolution)
                            record['quality_score'] = f"{res_float:.2f} Å"
                            record['quality_type'] = 'Resolution'
                        except (ValueError, TypeError):
                            record['quality_score'] = 'N/A'
                            record['quality_type'] = 'N/A'
                    else:
                        record['quality_score'] = 'N/A'
                        record['quality_type'] = 'N/A'
                elif metadata.get('type') == 'computed_model':
                    # For computed models, use confidence score (higher is better)
                    confidence = metadata.get('model_confidence')
                    if confidence and confidence != 'Unknown' and confidence != 'N/A':
                        try:
                            conf_float = float(confidence)
                            record['quality_score'] = f"{conf_float:.2f}"
                            record['quality_type'] = 'Confidence'
                        except (ValueError, TypeError):
                            record['quality_score'] = str(confidence)
                            record['quality_type'] = 'Confidence'
                    else:
                        record['quality_score'] = 'N/A'
                        record['quality_type'] = 'N/A'
                else:
                    record['quality_score'] = 'N/A'
                    record['quality_type'] = 'N/A'
                
                metadata_records.append(record)
                
        except Exception as e:
            logging.warning(f"Error processing metadata for gene {gene_name}: {e}")
            continue
    
    logging.info(f"Extracted metadata for {len(metadata_records)} structures across {len(gene_dirs)} genes")
    return metadata_records

def create_summary_statistics(metadata_records):
    """
    Create summary statistics from metadata records
    
    Args:
        metadata_records (list): List of metadata dictionaries
    
    Returns:
        dict: Summary statistics
    """
    
    if not metadata_records:
        return {}
    
    df = pd.DataFrame(metadata_records)
    
    stats = {
        'total_structures': len(df),
        'total_genes': df['gene_name'].nunique(),
        'experimental_structures': len(df[df['structure_type'] == 'experimental_structure']),
        'computed_models': len(df[df['structure_type'] == 'computed_model']),
        'genes_with_experimental': df[df['structure_type'] == 'experimental_structure']['gene_name'].nunique(),
        'genes_with_computed': df[df['structure_type'] == 'computed_model']['gene_name'].nunique(),
        'unique_organisms': df['organism'].nunique(),
        'avg_structures_per_gene': len(df) / df['gene_name'].nunique() if df['gene_name'].nunique() > 0 else 0
    }
    
    # Database breakdown
    db_counts = df['database'].value_counts().to_dict()
    stats['database_breakdown'] = db_counts
    
    # Resolution statistics for experimental structures
    exp_df = df[df['structure_type'] == 'experimental_structure']
    if not exp_df.empty:
        # Extract numeric resolution values
        resolutions = []
        for res in exp_df['resolution']:
            if res != 'N/A' and res is not None:
                try:
                    resolutions.append(float(res))
                except (ValueError, TypeError):
                    continue
        
        if resolutions:
            stats['resolution_stats'] = {
                'mean_resolution': sum(resolutions) / len(resolutions),
                'min_resolution': min(resolutions),
                'max_resolution': max(resolutions),
                'count_with_resolution': len(resolutions)
            }
    
    return stats

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get parameters from Snakemake
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        output_file = Path(snakemake.output.metadata_tsv)
        
        # Derive structures directory from output file path
        structures_base_dir = output_file.parent
        
    except NameError:
        # Running standalone - use default values
        analysis = "analysis1"
        paramset = "params1"
        structures_base_dir = Path("data/proteins_3d_structure")
        output_file = Path(f"data/proteins_3d_structure/{analysis}_{paramset}_3d_structures_metadata.tsv")
    
    logging.info(f"Extracting 3D structure metadata for {analysis}_{paramset}")
    logging.info(f"Structures directory: {structures_base_dir}")
    logging.info(f"Output file: {output_file}")
    
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Extract metadata
    metadata_records = extract_metadata_for_analysis(structures_base_dir, analysis, paramset)
    
    if not metadata_records:
        logging.warning("No metadata records found!")
        # Create empty TSV with headers
        empty_df = pd.DataFrame(columns=[
            'analysis', 'paramset', 'gene_name', 'structure_id', 'structure_type',
            'database', 'organism', 'protein_name', 'experimental_method', 'resolution',
            'model_confidence', 'quality_score', 'quality_type', 'deposit_date', 'journal',
            'authors', 'uniprot_accession', 'model_version', 'sequence_length'
        ])
        empty_df.to_csv(output_file, sep='\t', index=False)
        return
    
    # Create DataFrame and sort by gene name, then by quality score
    df = pd.DataFrame(metadata_records)
    
    # Sort by gene name, then by structure type (experimental first), then by quality
    def sort_key(row):
        gene_name = row['gene_name']
        struct_type = 0 if row['structure_type'] == 'experimental_structure' else 1
        
        # For quality, lower resolution is better for experimental, higher confidence is better for computed
        quality = 999  # Default for unknown quality
        if row['structure_type'] == 'experimental_structure' and row['resolution'] != 'N/A':
            try:
                quality = float(row['resolution'])  # Lower is better
            except (ValueError, TypeError):
                quality = 999
        elif row['structure_type'] == 'computed_model' and row['model_confidence'] != 'N/A':
            try:
                quality = -float(row['model_confidence'])  # Higher is better, so negate for sorting
            except (ValueError, TypeError):
                quality = 999
        
        return (gene_name, struct_type, quality)
    
    df['sort_key'] = df.apply(sort_key, axis=1)
    df = df.sort_values('sort_key').drop('sort_key', axis=1)
    
    # Reorder columns for better readability
    column_order = [
        'analysis', 'paramset', 'gene_name', 'structure_id', 'structure_type',
        'database', 'organism', 'protein_name', 'quality_score', 'quality_type',
        'experimental_method', 'resolution', 'model_confidence', 'deposit_date',
        'journal', 'authors', 'uniprot_accession', 'model_version', 'sequence_length'
    ]
    
    df = df[column_order]
    
    # Save to TSV
    df.to_csv(output_file, sep='\t', index=False)
    
    # Generate summary statistics
    stats = create_summary_statistics(metadata_records)
    
    # Save summary statistics
    stats_file = output_file.parent / f"3d_structures_summary.json"
    stats['generation_date'] = datetime.now().isoformat()
    stats['analysis'] = analysis
    stats['paramset'] = paramset
    
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Log summary
    logging.info(f"=== 3D Structure Metadata Summary ===")
    logging.info(f"Total structures: {stats.get('total_structures', 0)}")
    logging.info(f"Total genes: {stats.get('total_genes', 0)}")
    logging.info(f"  - Experimental structures: {stats.get('experimental_structures', 0)}")
    logging.info(f"  - Computed models: {stats.get('computed_models', 0)}")
    logging.info(f"Genes with structures:")
    logging.info(f"  - With experimental: {stats.get('genes_with_experimental', 0)}")
    logging.info(f"  - With computed: {stats.get('genes_with_computed', 0)}")
    logging.info(f"Average structures per gene: {stats.get('avg_structures_per_gene', 0):.2f}")
    logging.info(f"Unique organisms: {stats.get('unique_organisms', 0)}")
    
    if 'resolution_stats' in stats:
        res_stats = stats['resolution_stats']
        logging.info(f"Resolution statistics:")
        logging.info(f"  - Mean: {res_stats['mean_resolution']:.2f} Å")
        logging.info(f"  - Range: {res_stats['min_resolution']:.2f} - {res_stats['max_resolution']:.2f} Å")
        logging.info(f"  - Structures with resolution: {res_stats['count_with_resolution']}")
    
    logging.info(f"Metadata TSV saved to: {output_file}")
    logging.info(f"Summary statistics saved to: {stats_file}")

if __name__ == "__main__":
    main()