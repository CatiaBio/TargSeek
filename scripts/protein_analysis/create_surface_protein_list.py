#!/usr/bin/env python3
"""
Create filtered gene list containing only genes with surface proteins.

This script analyzes topology prediction results and creates a simple list
of genes that have sufficient extracellular topology coverage. This list
is then used by subsequent workflow steps to process only surface proteins.

Input:
- Topology prediction results (JSON format from predict_topology.py)
- Configuration parameters for minimum extracellular coverage

Output:
- Filtered gene list files for positive and negative groups
- Summary of filtering statistics
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

def identify_surface_genes(topology_data: Dict, min_extracellular_coverage: float) -> Dict[str, Set[str]]:
    """
    Identify genes with at least one protein structure having sufficient extracellular coverage.
    
    Returns:
        Dict mapping group ('positive'/'negative') to set of surface gene names
    """
    surface_genes = {}
    filter_stats = {
        'total_genes': 0,
        'surface_genes': 0,
        'filtered_out_genes': 0,
        'by_group': {}
    }
    
    for group, genes in topology_data.items():
        surface_gene_set = set()
        group_stats = {
            'total_genes': 0,
            'surface_genes': 0,
            'filtered_out_genes': 0,
            'gene_details': {}
        }
        
        for gene, sequences in genes.items():
            group_stats['total_genes'] += 1
            filter_stats['total_genes'] += 1
            
            has_surface_structure = False
            max_coverage = 0.0
            structure_count = 0
            
            # Check all structures for this gene
            for seq_id, regions in sequences.items():
                structure_count += 1
                extracellular_coverage = calculate_extracellular_coverage(regions)
                max_coverage = max(max_coverage, extracellular_coverage)
                
                if extracellular_coverage >= min_extracellular_coverage:
                    has_surface_structure = True
                
                logger.debug(f"  {seq_id}: {extracellular_coverage:.2%} extracellular coverage")
            
            gene_details = {
                'max_extracellular_coverage': max_coverage,
                'has_surface_structures': has_surface_structure,
                'structure_count': structure_count
            }
            
            if has_surface_structure:
                surface_gene_set.add(gene)
                group_stats['surface_genes'] += 1
                filter_stats['surface_genes'] += 1
                logger.info(f"Gene {gene} (gram_{group}): INCLUDED - max coverage {max_coverage:.2%}")
            else:
                group_stats['filtered_out_genes'] += 1
                filter_stats['filtered_out_genes'] += 1
                logger.warning(f"Gene {gene} (gram_{group}): EXCLUDED - max coverage {max_coverage:.2%}")
            
            group_stats['gene_details'][gene] = gene_details
        
        surface_genes[group] = surface_gene_set
        filter_stats['by_group'][group] = group_stats
        
        logger.info(f"Gram {group}: {len(surface_gene_set)}/{len(genes)} genes have surface structures")
    
    # Log overall statistics
    logger.info(f"Overall filtering results:")
    logger.info(f"  Total genes analyzed: {filter_stats['total_genes']}")
    logger.info(f"  Surface genes kept: {filter_stats['surface_genes']}")
    logger.info(f"  Genes filtered out: {filter_stats['filtered_out_genes']}")
    logger.info(f"  Surface coverage threshold: {min_extracellular_coverage:.1%}")
    
    return surface_genes, filter_stats

def save_gene_lists(surface_genes: Dict[str, Set[str]], output_positive: str, output_negative: str):
    """Save gene lists to text files."""
    
    # Save positive genes
    os.makedirs(os.path.dirname(output_positive), exist_ok=True)
    with open(output_positive, 'w') as f:
        for gene in sorted(surface_genes.get('positive', set())):
            f.write(f"{gene}\n")
    
    logger.info(f"Saved {len(surface_genes.get('positive', set()))} positive surface genes to: {output_positive}")
    
    # Save negative genes
    os.makedirs(os.path.dirname(output_negative), exist_ok=True)
    with open(output_negative, 'w') as f:
        for gene in sorted(surface_genes.get('negative', set())):
            f.write(f"{gene}\n")
    
    logger.info(f"Saved {len(surface_genes.get('negative', set()))} negative surface genes to: {output_negative}")

def save_filter_summary(filter_stats: Dict, surface_genes: Dict, output_file: str):
    """Save filtering summary to JSON file."""
    summary = {
        'filter_parameters': {
            'min_extracellular_coverage': snakemake.params.min_extracellular_coverage,
            'analysis': snakemake.params.analysis,
            'paramset': snakemake.params.paramset
        },
        'filtering_stats': filter_stats,
        'surface_genes_by_group': {
            group: sorted(list(genes)) for group, genes in surface_genes.items()
        },
        'summary': {
            'total_genes_analyzed': filter_stats['total_genes'],
            'surface_genes_kept': filter_stats['surface_genes'],
            'genes_filtered_out': filter_stats['filtered_out_genes'],
            'retention_rate': filter_stats['surface_genes'] / max(filter_stats['total_genes'], 1)
        }
    }
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Saved filtering summary: {output_file}")

def main():
    """Main function to create surface gene list."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    min_extracellular_coverage = snakemake.params.min_extracellular_coverage
    
    # Input files
    topology_results = snakemake.input.topology_results
    
    # Output files
    surface_proteins_positive = snakemake.output.surface_proteins_positive
    surface_proteins_negative = snakemake.output.surface_proteins_negative
    protein_filter_summary = snakemake.output.protein_filter_summary
    
    logger.info(f"Creating surface gene list for {analysis}_{paramset}")
    logger.info(f"Minimum extracellular coverage threshold: {min_extracellular_coverage:.1%}")
    
    # Load topology prediction results
    topology_data = load_topology_results(topology_results)
    
    if not topology_data:
        logger.error("No topology data available for filtering")
        sys.exit(1)
    
    # Identify genes with sufficient surface topology
    surface_genes, filter_stats = identify_surface_genes(topology_data, min_extracellular_coverage)
    
    total_surface_genes = sum(len(genes) for genes in surface_genes.values())
    
    if total_surface_genes == 0:
        logger.error("No genes with sufficient surface topology found!")
        logger.error("Consider lowering the min_extracellular_coverage threshold in config/config_analysis.yaml")
        sys.exit(1)
    
    # Save gene lists
    save_gene_lists(surface_genes, surface_proteins_positive, surface_proteins_negative)
    
    # Save filtering summary
    save_filter_summary(filter_stats, surface_genes, protein_filter_summary)
    
    logger.info("Surface gene list creation completed successfully")
    logger.info(f"Total surface genes: {total_surface_genes}")
    
    # Warn if many genes were filtered out
    retention_rate = filter_stats['surface_genes'] / filter_stats['total_genes']
    if retention_rate < 0.5:  # If less than 50% retained
        logger.warning(f"Only {retention_rate:.1%} of genes retained after surface filtering")
        logger.warning("Consider reviewing the min_extracellular_coverage threshold in the configuration")

if __name__ == "__main__":
    main()