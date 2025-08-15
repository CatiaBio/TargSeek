#!/usr/bin/env python3
"""
Map extracellular regions using topology predictions and PDB coordinates.

This script combines topology predictions with PDB structure mappings to identify
which residue positions in the 3D structures correspond to extracellular regions.

Input:
- Topology prediction results (JSON format)
- PDB numbering mapping from download pipeline 
- Selected 3D structure paths
- Protein sequences

Output:
- Extracellular region mapping with PDB coordinates
- Summary of extracellular coverage per structure
- Residue-level annotations for each protein
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import logging
from typing import Dict, List, Tuple, Any, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_topology_predictions(topology_file: str) -> Dict[str, List[Dict]]:
    """Load topology predictions from JSON file."""
    try:
        with open(topology_file, 'r') as f:
            topology_data = json.load(f)
        logger.info(f"Loaded topology predictions for {len(topology_data)} proteins")
        return topology_data
    except Exception as e:
        logger.error(f"Error loading topology predictions: {e}")
        return {}

def load_pdb_mapping(mapping_file: str) -> pd.DataFrame:
    """Load PDB structure mapping from TSV file."""
    try:
        df = pd.read_csv(mapping_file, sep='\t')
        logger.info(f"Loaded PDB mapping with {len(df)} entries")
        return df
    except Exception as e:
        logger.error(f"Error loading PDB mapping: {e}")
        return pd.DataFrame()

def parse_selected_structures(file_path: str) -> Dict[str, str]:
    """Parse selected 3D paths to get gene -> structure mapping."""
    gene_structure_map = {}
    
    if not os.path.exists(file_path):
        logger.warning(f"Selected structures file not found: {file_path}")
        return gene_structure_map
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                if ':' in line:
                    gene, structure = line.split(':', 1)
                    gene_structure_map[gene.strip()] = structure.strip()
    
    return gene_structure_map

def extract_extracellular_regions(topology_data: Dict[str, List[Dict]]) -> Dict[str, List[Tuple[int, int]]]:
    """Extract extracellular regions from topology predictions."""
    extracellular_regions = {}
    
    for protein_id, regions in topology_data.items():
        extracellular_spans = []
        
        for region in regions:
            if region['type'] == 'outside':  # Extracellular
                start = region['start']
                end = region['end']
                extracellular_spans.append((start, end))
        
        if extracellular_spans:
            extracellular_regions[protein_id] = extracellular_spans
            logger.debug(f"{protein_id}: {len(extracellular_spans)} extracellular regions")
    
    logger.info(f"Found extracellular regions in {len(extracellular_regions)} proteins")
    return extracellular_regions

def map_sequence_to_pdb_positions(
    protein_id: str, 
    extracellular_spans: List[Tuple[int, int]], 
    pdb_mapping_df: pd.DataFrame
) -> Dict[str, List[Tuple[int, int]]]:
    """Map sequence positions to PDB coordinates using the mapping file."""
    
    # Extract gene name and structure ID from protein_id (format: gene_structureID)
    if '_' in protein_id:
        gene_name = protein_id.split('_')[0]
        structure_id = '_'.join(protein_id.split('_')[1:])
    else:
        logger.warning(f"Cannot parse gene and structure from protein_id: {protein_id}")
        return {}
    
    # Find matching entries in PDB mapping
    matching_entries = pdb_mapping_df[
        (pdb_mapping_df['gene_name'] == gene_name) & 
        (pdb_mapping_df['structure_id'].str.contains(structure_id, na=False))
    ]
    
    if matching_entries.empty:
        logger.warning(f"No PDB mapping found for {protein_id}")
        return {}
    
    pdb_extracellular_regions = {}
    
    for _, entry in matching_entries.iterrows():
        pdb_id = entry.get('pdb_id', structure_id)
        chain_id = entry.get('chain_id', 'A')
        
        # Get the sequence alignment/mapping if available
        seq_start = entry.get('seq_start', 1)
        seq_end = entry.get('seq_end', 1000)  # Default large number
        pdb_start = entry.get('pdb_start', 1)
        pdb_end = entry.get('pdb_end', 1000)
        
        # Map extracellular spans to PDB coordinates
        pdb_spans = []
        for span_start, span_end in extracellular_spans:
            # Check if span overlaps with the mapped region
            if span_start <= seq_end and span_end >= seq_start:
                # Calculate PDB coordinates
                overlap_start = max(span_start, seq_start)
                overlap_end = min(span_end, seq_end)
                
                # Map to PDB coordinates (simple linear mapping)
                pdb_span_start = pdb_start + (overlap_start - seq_start)
                pdb_span_end = pdb_start + (overlap_end - seq_start)
                
                pdb_spans.append((pdb_span_start, pdb_span_end))
        
        if pdb_spans:
            pdb_key = f"{pdb_id}_{chain_id}"
            pdb_extracellular_regions[pdb_key] = pdb_spans
    
    return pdb_extracellular_regions

def create_residue_annotations(
    extracellular_regions: Dict[str, List[Tuple[int, int]]],
    pdb_mapping_df: pd.DataFrame
) -> Dict[str, Dict]:
    """Create detailed residue-level annotations for each protein structure."""
    
    annotations = {}
    
    for protein_id, seq_spans in extracellular_regions.items():
        # Map to PDB coordinates
        pdb_regions = map_sequence_to_pdb_positions(protein_id, seq_spans, pdb_mapping_df)
        
        protein_annotation = {
            'protein_id': protein_id,
            'sequence_extracellular_spans': seq_spans,
            'pdb_extracellular_regions': pdb_regions,
            'total_extracellular_residues': 0,
            'extracellular_coverage': 0.0,
            'structures': {}
        }
        
        # Calculate total extracellular residues
        total_extracellular = sum(end - start + 1 for start, end in seq_spans)
        protein_annotation['total_extracellular_residues'] = total_extracellular
        
        # Create structure-specific annotations
        for pdb_key, pdb_spans in pdb_regions.items():
            structure_annotation = {
                'pdb_extracellular_spans': pdb_spans,
                'pdb_extracellular_residues': sum(end - start + 1 for start, end in pdb_spans),
                'extracellular_residue_list': []
            }
            
            # Create list of individual extracellular residues
            for start, end in pdb_spans:
                structure_annotation['extracellular_residue_list'].extend(range(start, end + 1))
            
            protein_annotation['structures'][pdb_key] = structure_annotation
        
        annotations[protein_id] = protein_annotation
    
    return annotations

def calculate_coverage_statistics(annotations: Dict[str, Dict]) -> Dict[str, Any]:
    """Calculate summary statistics for extracellular coverage."""
    
    stats = {
        'total_proteins': len(annotations),
        'proteins_with_extracellular': 0,
        'average_extracellular_residues': 0,
        'coverage_distribution': {},
        'structure_coverage': {}
    }
    
    total_extracellular_residues = 0
    coverage_bins = [0, 0.1, 0.25, 0.5, 0.75, 1.0]
    coverage_counts = {f"{coverage_bins[i]:.1f}-{coverage_bins[i+1]:.1f}": 0 
                      for i in range(len(coverage_bins)-1)}
    
    for protein_id, annotation in annotations.items():
        extracellular_residues = annotation['total_extracellular_residues']
        
        if extracellular_residues > 0:
            stats['proteins_with_extracellular'] += 1
            total_extracellular_residues += extracellular_residues
        
        # Calculate coverage (this would need total protein length for accurate calculation)
        # For now, we'll estimate based on extracellular regions only
        coverage = annotation.get('extracellular_coverage', 0.0)
        
        # Bin the coverage
        for i in range(len(coverage_bins)-1):
            if coverage_bins[i] <= coverage <= coverage_bins[i+1]:
                bin_key = f"{coverage_bins[i]:.1f}-{coverage_bins[i+1]:.1f}"
                coverage_counts[bin_key] += 1
                break
        
        # Structure-specific statistics
        for structure_id, struct_data in annotation['structures'].items():
            if structure_id not in stats['structure_coverage']:
                stats['structure_coverage'][structure_id] = {
                    'extracellular_residues': 0,
                    'total_spans': 0
                }
            
            stats['structure_coverage'][structure_id]['extracellular_residues'] += struct_data['pdb_extracellular_residues']
            stats['structure_coverage'][structure_id]['total_spans'] += len(struct_data['pdb_extracellular_spans'])
    
    if stats['proteins_with_extracellular'] > 0:
        stats['average_extracellular_residues'] = total_extracellular_residues / stats['proteins_with_extracellular']
    
    stats['coverage_distribution'] = coverage_counts
    
    return stats

def create_extracellular_bed_format(annotations: Dict[str, Dict], output_dir: str):
    """Create BED format files for visualization of extracellular regions."""
    
    bed_dir = os.path.join(output_dir, "bed_files")
    os.makedirs(bed_dir, exist_ok=True)
    
    for protein_id, annotation in annotations.items():
        for structure_id, struct_data in annotation['structures'].items():
            bed_file = os.path.join(bed_dir, f"{structure_id}_extracellular.bed")
            
            with open(bed_file, 'w') as f:
                f.write(f"track name=\"{structure_id}_extracellular\" description=\"Extracellular regions for {structure_id}\"\n")
                
                for i, (start, end) in enumerate(struct_data['pdb_extracellular_spans']):
                    f.write(f"{structure_id}\t{start-1}\t{end}\textracellular_region_{i+1}\t1000\t+\n")
    
    logger.info(f"Created BED files in {bed_dir}")

def main():
    """Main function to map extracellular regions."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    topology_dir = snakemake.params.topology_dir
    
    # Input files
    topology_results = snakemake.input.topology_results
    pdb_mapping_file = snakemake.input.pdb_numbering_mapping
    selected_3d_positive = snakemake.input.selected_3d_paths_positive
    selected_3d_negative = snakemake.input.selected_3d_paths_negative
    
    # Output files
    extracellular_mapping = snakemake.output.extracellular_mapping
    extracellular_summary = snakemake.output.extracellular_summary
    
    logger.info(f"Starting extracellular region mapping for {analysis}_{paramset}")
    
    # Create output directory
    os.makedirs(os.path.dirname(extracellular_mapping), exist_ok=True)
    
    # Load topology predictions
    topology_data = load_topology_predictions(topology_results)
    if not topology_data:
        logger.error("No topology data loaded")
        sys.exit(1)
    
    # Load PDB mapping
    pdb_mapping_df = load_pdb_mapping(pdb_mapping_file)
    if pdb_mapping_df.empty:
        logger.error("No PDB mapping data loaded")
        sys.exit(1)
    
    # Extract extracellular regions from topology predictions
    extracellular_regions = extract_extracellular_regions(topology_data)
    
    if not extracellular_regions:
        logger.warning("No extracellular regions found in topology predictions")
        # Create empty output files
        with open(extracellular_mapping, 'w') as f:
            json.dump({}, f)
        with open(extracellular_summary, 'w') as f:
            json.dump({'total_proteins': 0, 'proteins_with_extracellular': 0}, f)
        return
    
    # Create detailed residue annotations
    annotations = create_residue_annotations(extracellular_regions, pdb_mapping_df)
    
    # Save detailed mapping
    with open(extracellular_mapping, 'w') as f:
        json.dump(annotations, f, indent=2)
    
    # Calculate coverage statistics
    coverage_stats = calculate_coverage_statistics(annotations)
    
    # Save summary statistics
    with open(extracellular_summary, 'w') as f:
        json.dump(coverage_stats, f, indent=2)
    
    # Create BED files for visualization
    create_extracellular_bed_format(annotations, topology_dir)
    
    logger.info("Extracellular region mapping completed successfully")
    logger.info(f"Mapped {len(annotations)} proteins with extracellular regions")
    logger.info(f"Average extracellular residues per protein: {coverage_stats['average_extracellular_residues']:.1f}")

if __name__ == "__main__":
    main()