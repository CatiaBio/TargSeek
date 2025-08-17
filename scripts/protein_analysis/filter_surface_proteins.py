#!/usr/bin/env python3
"""
Filter proteins to only include those with sufficient extracellular (surface) topology.

This script filters the selected 3D protein structures based on their topology predictions,
keeping only proteins that have adequate extracellular regions for epitope prediction.
Proteins without sufficient surface exposure are excluded from downstream analysis.

Input:
- Topology prediction results (JSON format from predict_topology.py)
- Selected 3D paths files for positive and negative groups
- Configuration parameters for minimum extracellular coverage

Output:
- Filtered 3D paths files for positive and negative groups (excluding non-surface proteins)
- Summary report of filtering results
"""

import os
import sys
import json
import logging
from pathlib import Path
from typing import Dict, List, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_topology_results(topology_file: str) -> Dict:
    """Load topology prediction results from JSON file."""
    if not os.path.exists(topology_file):
        logger.error(f"Topology results file not found: {topology_file}")
        return {}
    
    with open(topology_file, 'r') as f:
        topology_data = json.load(f)
    
    logger.info(f"Loaded topology data for {len(topology_data)} groups")
    return topology_data

def calculate_extracellular_coverage(regions: List[Dict]) -> float:
    """Calculate the percentage of protein sequence that is extracellular."""
    total_length = 0
    extracellular_length = 0
    
    for region in regions:
        region_length = region.get('length', 0)
        total_length += region_length
        
        # Consider 'outside' and 'signal' regions as extracellular
        if region.get('type') in ['outside', 'signal']:
            extracellular_length += region_length
    
    if total_length == 0:
        return 0.0
    
    return extracellular_length / total_length

def identify_surface_proteins(topology_data: Dict, min_extracellular_coverage: float) -> Dict[str, Set[str]]:
    """
    Identify proteins with sufficient extracellular coverage.
    
    Returns:
        Dict mapping group to set of surface protein gene names
    """
    surface_proteins = {}
    filter_stats = {
        'total_proteins': 0,
        'surface_proteins': 0,
        'filtered_out': 0,
        'by_group': {}
    }
    
    for group, genes in topology_data.items():
        surface_genes = set()
        group_stats = {
            'total_proteins': 0,
            'surface_proteins': 0,
            'filtered_out': 0,
            'gene_details': {}
        }
        
        for gene, sequences in genes.items():
            gene_surface_count = 0
            gene_total_count = 0
            gene_details = {
                'sequences': {},
                'has_surface_structures': False,
                'max_extracellular_coverage': 0.0
            }
            
            for seq_id, regions in sequences.items():
                gene_total_count += 1
                filter_stats['total_proteins'] += 1
                group_stats['total_proteins'] += 1
                
                # Calculate extracellular coverage for this sequence
                extracellular_coverage = calculate_extracellular_coverage(regions)
                gene_details['sequences'][seq_id] = {
                    'extracellular_coverage': extracellular_coverage,
                    'is_surface': extracellular_coverage >= min_extracellular_coverage,
                    'regions': regions
                }
                
                # Track maximum coverage for the gene
                gene_details['max_extracellular_coverage'] = max(
                    gene_details['max_extracellular_coverage'], 
                    extracellular_coverage
                )
                
                # Check if this sequence meets surface criteria
                if extracellular_coverage >= min_extracellular_coverage:
                    gene_surface_count += 1
                    gene_details['has_surface_structures'] = True
                
                logger.debug(f"  {seq_id}: {extracellular_coverage:.2%} extracellular coverage")
            
            # Include gene if it has at least one surface structure
            if gene_details['has_surface_structures']:
                surface_genes.add(gene)
                filter_stats['surface_proteins'] += gene_surface_count
                group_stats['surface_proteins'] += gene_surface_count
                logger.info(f"Gene {gene} (gram_{group}): {gene_surface_count}/{gene_total_count} surface structures, max coverage: {gene_details['max_extracellular_coverage']:.2%}")
            else:
                filter_stats['filtered_out'] += gene_total_count
                group_stats['filtered_out'] += gene_total_count
                logger.warning(f"Gene {gene} (gram_{group}): NO surface structures (max coverage: {gene_details['max_extracellular_coverage']:.2%}) - EXCLUDED")
            
            group_stats['gene_details'][gene] = gene_details
        
        surface_proteins[group] = surface_genes
        filter_stats['by_group'][group] = group_stats
        
        logger.info(f"Gram {group}: {len(surface_genes)} genes with surface structures out of {len(genes)} total genes")
    
    # Log overall statistics
    logger.info(f"Overall filtering results:")
    logger.info(f"  Total proteins analyzed: {filter_stats['total_proteins']}")
    logger.info(f"  Surface proteins kept: {filter_stats['surface_proteins']}")
    logger.info(f"  Proteins filtered out: {filter_stats['filtered_out']}")
    logger.info(f"  Surface coverage threshold: {min_extracellular_coverage:.1%}")
    
    return surface_proteins, filter_stats

def filter_3d_paths_file(input_paths_file: str, output_paths_file: str, surface_genes: Set[str], group: str) -> int:
    """
    Filter a 3D paths file to only include genes with surface topology.
    
    Returns:
        Number of genes kept after filtering
    """
    if not os.path.exists(input_paths_file):
        logger.warning(f"Input paths file not found: {input_paths_file}")
        # Create empty output file
        Path(output_paths_file).touch()
        return 0
    
    kept_count = 0
    filtered_count = 0
    
    os.makedirs(os.path.dirname(output_paths_file), exist_ok=True)
    
    with open(input_paths_file, 'r') as infile, open(output_paths_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line and not line.startswith('#'):
                # Extract gene name from path: data/protein_structures/bamA/5OR1_1.fasta -> bamA
                path_parts = line.split('/')
                if len(path_parts) >= 3:
                    gene = path_parts[-2]  # Gene name is second to last part
                    
                    if gene in surface_genes:
                        outfile.write(line + '\n')
                        kept_count += 1
                    else:
                        filtered_count += 1
                        logger.debug(f"Filtered out {gene} (gram_{group}): no surface topology")
                else:
                    # Keep lines that don't match expected format (comments, etc.)
                    outfile.write(line + '\n')
    
    logger.info(f"Filtered {group} 3D paths: kept {kept_count}, removed {filtered_count}")
    return kept_count

def save_filter_summary(filter_stats: Dict, surface_proteins: Dict, output_file: str):
    """Save filtering summary to JSON file."""
    summary = {
        'filter_parameters': {
            'min_extracellular_coverage': snakemake.params.min_extracellular_coverage,
            'timestamp': str(Path().cwd())
        },
        'filtering_stats': filter_stats,
        'surface_proteins_by_group': {
            group: list(genes) for group, genes in surface_proteins.items()
        },
        'summary': {
            'total_genes_analyzed': sum(len(genes) for genes in surface_proteins.values()),
            'surface_genes_kept': sum(len(genes) for genes in surface_proteins.values()),
            'genes_filtered_out': filter_stats['by_group'].get('positive', {}).get('gene_details', {}) and 
                                 filter_stats['by_group'].get('negative', {}).get('gene_details', {})
        }
    }
    
    # Calculate genes filtered out
    all_genes = set()
    for group_data in filter_stats['by_group'].values():
        all_genes.update(group_data.get('gene_details', {}).keys())
    
    surface_genes_all = set()
    for genes in surface_proteins.values():
        surface_genes_all.update(genes)
    
    summary['summary']['genes_filtered_out'] = len(all_genes) - len(surface_genes_all)
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Saved filtering summary: {output_file}")

def main():
    """Main function to filter proteins by surface topology."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    min_extracellular_coverage = snakemake.params.min_extracellular_coverage
    
    # Input files
    topology_results = snakemake.input.topology_results
    selected_3d_paths_positive = snakemake.input.selected_3d_paths_positive
    selected_3d_paths_negative = snakemake.input.selected_3d_paths_negative
    
    # Output files
    filtered_3d_paths_positive = snakemake.output.filtered_3d_paths_positive
    filtered_3d_paths_negative = snakemake.output.filtered_3d_paths_negative
    topology_filter_summary = snakemake.output.topology_filter_summary
    
    logger.info(f"Starting surface topology filtering for {analysis}_{paramset}")
    logger.info(f"Minimum extracellular coverage threshold: {min_extracellular_coverage:.1%}")
    
    # Load topology prediction results
    topology_data = load_topology_results(topology_results)
    
    if not topology_data:
        logger.error("No topology data available for filtering")
        sys.exit(1)
    
    # Identify proteins with sufficient surface topology
    surface_proteins, filter_stats = identify_surface_proteins(topology_data, min_extracellular_coverage)
    
    # Filter 3D paths files for both groups
    positive_kept = filter_3d_paths_file(
        selected_3d_paths_positive, 
        filtered_3d_paths_positive, 
        surface_proteins.get('positive', set()),
        'positive'
    )
    
    negative_kept = filter_3d_paths_file(
        selected_3d_paths_negative, 
        filtered_3d_paths_negative, 
        surface_proteins.get('negative', set()),
        'negative'
    )
    
    total_kept = positive_kept + negative_kept
    
    if total_kept == 0:
        logger.error("No proteins with sufficient surface topology found!")
        logger.error("Consider lowering the min_extracellular_coverage threshold in config/config_analysis.yaml")
        sys.exit(1)
    
    # Save filtering summary
    save_filter_summary(filter_stats, surface_proteins, topology_filter_summary)
    
    logger.info("Surface topology filtering completed successfully")
    logger.info(f"Total proteins kept: {total_kept} (positive: {positive_kept}, negative: {negative_kept})")
    
    # Warn if many proteins were filtered out
    total_analyzed = filter_stats['total_proteins']
    if total_kept < total_analyzed * 0.5:  # If less than 50% kept
        logger.warning(f"Only {total_kept}/{total_analyzed} ({100*total_kept/total_analyzed:.1f}%) proteins kept after surface filtering")
        logger.warning("Consider reviewing the min_extracellular_coverage threshold in the configuration")

if __name__ == "__main__":
    main()