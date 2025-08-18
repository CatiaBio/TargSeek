#!/usr/bin/env python3
"""
Create extracellular region mapping from topology predictions.

This script converts the detailed topology predictions into a simplified mapping
format used by the epitope filtering process. It extracts extracellular regions
for each protein structure, which will be used to filter BepiPred epitope 
predictions to keep only those in accessible (extracellular) regions.

Input:
- Topology prediction results (JSON format from predict_topology.py)
- Configuration parameters for minimum extracellular coverage

Output:
- Extracellular region mapping (JSON format for epitope filtering)
"""

import os
import sys
import json
import logging
from pathlib import Path
from typing import Dict, List, Any

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

def convert_to_extracellular_mapping(topology_data: Dict) -> Dict[str, Dict]:
    """
    Convert topology predictions to extracellular region mapping.
    
    This creates a mapping of protein_structure_id -> extracellular_regions
    that can be used by the epitope filtering script.
    
    Args:
        topology_data: Topology predictions organized by group -> gene -> sequence -> regions
        
    Returns:
        Dict mapping protein_structure_id to extracellular region information
    """
    extracellular_mapping = {}
    
    for group, genes in topology_data.items():
        for gene, sequences in genes.items():
            for seq_id, regions in sequences.items():
                # Create protein structure ID
                protein_structure_id = f"{gene}_{seq_id}" if seq_id != gene else gene
                
                # Extract extracellular regions
                extracellular_regions = []
                total_length = 0
                extracellular_length = 0
                
                for region in regions:
                    region_length = region.get('length', 0)
                    total_length += region_length
                    
                    # Consider 'outside' and 'signal' regions as extracellular
                    if region.get('type') in ['outside', 'signal']:
                        extracellular_regions.append({
                            'start': region.get('start', 1),
                            'end': region.get('end', region_length),
                            'length': region_length,
                            'type': region.get('type'),
                            'description': region.get('description', '')
                        })
                        extracellular_length += region_length
                
                # Calculate coverage
                extracellular_coverage = extracellular_length / max(total_length, 1)
                
                # Store mapping information
                extracellular_mapping[protein_structure_id] = {
                    'gene': gene,
                    'structure_id': seq_id,
                    'group': group,
                    'total_length': total_length,
                    'extracellular_length': extracellular_length,
                    'extracellular_coverage': extracellular_coverage,
                    'extracellular_regions': extracellular_regions,
                    'region_count': len(extracellular_regions),
                    'topology_regions': regions  # Keep original topology data for reference
                }
                
                logger.debug(f"Mapped {protein_structure_id}: {len(extracellular_regions)} extracellular regions, {extracellular_coverage:.2%} coverage")
    
    logger.info(f"Created extracellular mapping for {len(extracellular_mapping)} protein structures")
    return extracellular_mapping

def create_mapping_summary(extracellular_mapping: Dict, min_coverage: float) -> Dict:
    """Create summary statistics for the extracellular mapping."""
    
    summary = {
        'total_proteins': len(extracellular_mapping),
        'proteins_with_extracellular': 0,
        'proteins_above_threshold': 0,
        'average_coverage': 0.0,
        'coverage_distribution': {
            '0-10%': 0,
            '10-25%': 0,
            '25-50%': 0,
            '50-75%': 0,
            '75-90%': 0,
            '90-100%': 0
        },
        'by_group': {},
        'min_coverage_threshold': min_coverage
    }
    
    total_coverage = 0.0
    group_stats = {}
    
    for protein_id, protein_data in extracellular_mapping.items():
        coverage = protein_data['extracellular_coverage']
        group = protein_data['group']
        
        # Initialize group stats if needed
        if group not in group_stats:
            group_stats[group] = {
                'total': 0,
                'with_extracellular': 0,
                'above_threshold': 0,
                'average_coverage': 0.0,
                'total_coverage': 0.0
            }
        
        # Update overall stats
        total_coverage += coverage
        if protein_data['region_count'] > 0:
            summary['proteins_with_extracellular'] += 1
            group_stats[group]['with_extracellular'] += 1
        
        if coverage >= min_coverage:
            summary['proteins_above_threshold'] += 1
            group_stats[group]['above_threshold'] += 1
        
        # Update group stats
        group_stats[group]['total'] += 1
        group_stats[group]['total_coverage'] += coverage
        
        # Coverage distribution
        if coverage == 0:
            summary['coverage_distribution']['0-10%'] += 1
        elif coverage < 0.25:
            summary['coverage_distribution']['10-25%'] += 1
        elif coverage < 0.5:
            summary['coverage_distribution']['25-50%'] += 1
        elif coverage < 0.75:
            summary['coverage_distribution']['50-75%'] += 1
        elif coverage < 0.9:
            summary['coverage_distribution']['75-90%'] += 1
        else:
            summary['coverage_distribution']['90-100%'] += 1
    
    # Calculate averages
    if summary['total_proteins'] > 0:
        summary['average_coverage'] = total_coverage / summary['total_proteins']
    
    # Calculate group averages
    for group, stats in group_stats.items():
        if stats['total'] > 0:
            stats['average_coverage'] = stats['total_coverage'] / stats['total']
    
    summary['by_group'] = group_stats
    
    return summary

def main():
    """Main function to create extracellular mapping."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    min_extracellular_coverage = snakemake.params.min_extracellular_coverage
    
    # Input files
    topology_results = snakemake.input.topology_results
    
    # Output files
    extracellular_mapping_file = snakemake.output.extracellular_mapping
    
    logger.info(f"Creating extracellular mapping for {analysis}_{paramset}")
    logger.info(f"Minimum extracellular coverage threshold: {min_extracellular_coverage:.1%}")
    
    # Load topology prediction results
    topology_data = load_topology_results(topology_results)
    
    if not topology_data:
        logger.error("No topology data available for creating extracellular mapping")
        sys.exit(1)
    
    # Convert topology data to extracellular mapping
    extracellular_mapping = convert_to_extracellular_mapping(topology_data)
    
    if not extracellular_mapping:
        logger.error("No extracellular regions found in topology data")
        sys.exit(1)
    
    # Create summary statistics
    mapping_summary = create_mapping_summary(extracellular_mapping, min_extracellular_coverage)
    
    # Prepare output data
    output_data = {
        'metadata': {
            'analysis': analysis,
            'paramset': paramset,
            'min_extracellular_coverage': min_extracellular_coverage,
            'total_proteins': mapping_summary['total_proteins'],
            'proteins_with_extracellular': mapping_summary['proteins_with_extracellular'],
            'proteins_above_threshold': mapping_summary['proteins_above_threshold'],
            'created_from': topology_results
        },
        'summary': mapping_summary,
        'mapping': extracellular_mapping
    }
    
    # Save extracellular mapping
    os.makedirs(os.path.dirname(extracellular_mapping_file), exist_ok=True)
    with open(extracellular_mapping_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    logger.info(f"Saved extracellular mapping: {extracellular_mapping_file}")
    logger.info(f"Total proteins: {mapping_summary['total_proteins']}")
    logger.info(f"Proteins with extracellular regions: {mapping_summary['proteins_with_extracellular']}")
    logger.info(f"Proteins above coverage threshold ({min_extracellular_coverage:.1%}): {mapping_summary['proteins_above_threshold']}")
    
    # Log group statistics
    for group, stats in mapping_summary['by_group'].items():
        logger.info(f"Gram {group}: {stats['above_threshold']}/{stats['total']} proteins above threshold (avg coverage: {stats['average_coverage']:.1%})")
    
    logger.info("Extracellular mapping creation completed successfully")

if __name__ == "__main__":
    main()