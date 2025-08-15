#!/usr/bin/env python3
"""
Predict membrane protein topology using DeepTMHMM.

This script runs topology prediction on protein sequences from selected 3D structures
to identify transmembrane regions, signal peptides, and extracellular/intracellular domains.

Fixed version that:
1. Extracts the correct 3D structure sequences from MSA files
2. Uses real DeepTMHMM via biolib
3. Properly handles DeepTMHMM output formats

Input:
- Selected 3D structure paths for positive and negative groups
- Protein sequences from MSA files
- Configuration parameters for topology prediction

Output:
- Topology predictions for each protein structure
- Summary statistics of topology predictions
- Confidence scores and region annotations
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import subprocess
import tempfile
import logging
from typing import Dict, List, Tuple, Any
import re

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_selected_3d_paths(paths_file: str) -> Dict[str, str]:
    """Parse the selected 3D paths file to get gene -> FASTA file path mapping."""
    gene_fasta_map = {}
    
    if not os.path.exists(paths_file):
        logger.warning(f"Selected 3D paths file not found: {paths_file}")
        return gene_fasta_map
    
    with open(paths_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Extract gene name from path: data/protein_structures/bamA/5OR1_1.fasta -> bamA
                path_parts = line.split('/')
                if len(path_parts) >= 3:
                    gene = path_parts[-2]  # Gene name is second to last part
                    gene_fasta_map[gene] = line
    
    logger.info(f"Loaded {len(gene_fasta_map)} gene-structure path mappings from {paths_file}")
    return gene_fasta_map

def extract_structure_sequences_from_paths(gene_fasta_map: Dict[str, str]) -> Dict[str, str]:
    """Extract sequences directly from the individual 3D structure FASTA files."""
    sequences = {}
    
    for gene, fasta_path in gene_fasta_map.items():
        if os.path.exists(fasta_path):
            try:
                # Read the sequence from the individual structure FASTA file
                for record in SeqIO.parse(fasta_path, "fasta"):
                    # Extract PDB ID from the file name
                    pdb_id = os.path.basename(fasta_path).replace('.fasta', '')
                    
                    # Get the sequence without gaps
                    structure_sequence = str(record.seq).replace('-', '')
                    
                    seq_id = f"{gene}_{pdb_id}"
                    sequences[seq_id] = structure_sequence
                    
                    logger.info(f"Loaded 3D structure sequence for {gene} ({pdb_id}): {len(structure_sequence)} residues")
                    break  # Only take the first (and usually only) sequence
                    
            except Exception as e:
                logger.error(f"Error reading {fasta_path}: {e}")
        else:
            logger.warning(f"3D structure FASTA file not found: {fasta_path}")
    
    logger.info(f"Extracted {len(sequences)} 3D structure sequences for topology prediction")
    return sequences

def run_deeptmhmm_biolib(sequences: Dict[str, str], output_dir: str) -> bool:
    """Run DeepTMHMM using biolib with proper input/output handling."""
    
    # Create temporary FASTA file
    temp_fasta = os.path.join(output_dir, "topology_input.fasta")
    
    with open(temp_fasta, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n{sequence}\n")
    
    logger.info(f"Created input FASTA with {len(sequences)} sequences: {temp_fasta}")
    
    try:
        # Run DeepTMHMM via biolib
        cmd = [
            "biolib", "run", "DTU/DeepTMHMM",
            "--fasta", temp_fasta,
            "--local"  # Run locally for faster processing
        ]
        
        logger.info(f"Running DeepTMHMM via biolib: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir, timeout=600)  # 10 min timeout
        
        if result.returncode == 0:
            logger.info("DeepTMHMM completed successfully via biolib")
            logger.info(f"DeepTMHMM stdout: {result.stdout[:500]}...")  # First 500 chars
            
            # Check if output files were created
            output_files = [f for f in os.listdir(output_dir) if f.endswith(('.gff3', '.3line', '.png'))]
            if output_files:
                logger.info(f"DeepTMHMM output files: {output_files}")
                return True
            else:
                logger.warning("DeepTMHMM completed but no output files found")
                return False
        else:
            logger.error(f"DeepTMHMM failed with return code {result.returncode}")
            logger.error(f"stderr: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error("DeepTMHMM via biolib timed out")
        return False
    except Exception as e:
        logger.error(f"DeepTMHMM via biolib failed: {e}")
        return False

def parse_deeptmhmm_output(output_dir: str) -> Dict[str, List[Dict]]:
    """Parse DeepTMHMM output files to extract topology predictions."""
    topology_data = {}
    
    # Look for GFF3 output file (standard DeepTMHMM output)
    gff_files = [f for f in os.listdir(output_dir) if f.endswith('.gff3')]
    
    if not gff_files:
        logger.warning("No GFF3 files found in DeepTMHMM output")
        # Try to find 3line format
        line_files = [f for f in os.listdir(output_dir) if f.endswith('.3line')]
        if line_files:
            return parse_deeptmhmm_3line(os.path.join(output_dir, line_files[0]))
        else:
            logger.error("No DeepTMHMM output files found")
            return topology_data
    
    gff_file = os.path.join(output_dir, gff_files[0])
    logger.info(f"Parsing DeepTMHMM GFF3 output: {gff_file}")
    
    with open(gff_file, 'r') as f:
        current_seq = None
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 9:
                seq_id = parts[0]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                attributes = parts[8]
                
                if seq_id not in topology_data:
                    topology_data[seq_id] = []
                
                # Parse region type from attributes
                region_type = 'unknown'
                if 'TOPOLOGY=' in attributes:
                    topology_match = re.search(r'TOPOLOGY=([^;]+)', attributes)
                    if topology_match:
                        region_type = topology_match.group(1).lower()
                
                topology_data[seq_id].append({
                    'start': start,
                    'end': end,
                    'type': region_type,
                    'feature': feature_type,
                    'attributes': attributes,
                    'length': end - start + 1
                })
    
    logger.info(f"Parsed topology data for {len(topology_data)} sequences")
    return topology_data

def parse_deeptmhmm_3line(file_path: str) -> Dict[str, List[Dict]]:
    """Parse DeepTMHMM 3line format output."""
    topology_data = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('>'):
            seq_id = line[1:]  # Remove '>'
            
            if i + 2 < len(lines):
                topology_line = lines[i + 2].strip()  # Third line contains topology
                topology_data[seq_id] = parse_topology_string(topology_line)
            i += 3
        else:
            i += 1
    
    return topology_data

def parse_topology_string(topology_str: str) -> List[Dict]:
    """Parse topology string into regions."""
    regions = []
    current_type = None
    start_pos = 1
    
    for i, char in enumerate(topology_str):
        if char != current_type:
            if current_type is not None:
                # End previous region
                regions.append({
                    'start': start_pos,
                    'end': i,
                    'type': get_region_type(current_type),
                    'length': i - start_pos + 1
                })
            # Start new region
            current_type = char
            start_pos = i + 1
    
    # Add final region
    if current_type is not None:
        regions.append({
            'start': start_pos,
            'end': len(topology_str),
            'type': get_region_type(current_type),
            'length': len(topology_str) - start_pos + 1
        })
    
    return regions

def get_region_type(char: str) -> str:
    """Convert DeepTMHMM character to region type."""
    mapping = {
        'O': 'outside',      # Extracellular
        'M': 'membrane',     # Transmembrane
        'I': 'inside',       # Intracellular
        'S': 'signal',       # Signal peptide
        'P': 'periplasm'     # Periplasmic (for bacteria)
    }
    return mapping.get(char.upper(), 'unknown')

def create_topology_summary(topology_data: Dict[str, List[Dict]], output_file: str):
    """Create summary statistics of topology predictions."""
    summary_stats = {
        'total_proteins': len(topology_data),
        'proteins_with_tm_regions': 0,
        'proteins_extracellular_only': 0,
        'average_tm_regions_per_protein': 0,
        'extracellular_coverage': {},
        'detailed_predictions': {}
    }
    
    total_tm_regions = 0
    
    for seq_id, regions in topology_data.items():
        tm_count = 0
        total_length = 0
        extracellular_length = 0
        
        protein_info = {
            'regions': regions,
            'transmembrane_count': 0,
            'extracellular_coverage': 0.0,
            'has_transmembrane': False,
            'total_length': 0
        }
        
        for region in regions:
            region_length = region['length']
            total_length += region_length
            
            if region['type'] == 'membrane':
                tm_count += 1
                protein_info['has_transmembrane'] = True
            elif region['type'] in ['outside', 'signal']:  # Include signal peptides as extracellular
                extracellular_length += region_length
        
        protein_info['transmembrane_count'] = tm_count
        protein_info['total_length'] = total_length
        
        if total_length > 0:
            protein_info['extracellular_coverage'] = extracellular_length / total_length
        
        summary_stats['detailed_predictions'][seq_id] = protein_info
        summary_stats['extracellular_coverage'][seq_id] = protein_info['extracellular_coverage']
        
        total_tm_regions += tm_count
        
        if tm_count > 0:
            summary_stats['proteins_with_tm_regions'] += 1
        
        if tm_count == 0 and extracellular_length == total_length:
            summary_stats['proteins_extracellular_only'] += 1
    
    if summary_stats['total_proteins'] > 0:
        summary_stats['average_tm_regions_per_protein'] = total_tm_regions / summary_stats['total_proteins']
    
    # Save summary to JSON
    with open(output_file, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    
    logger.info(f"Created topology summary: {output_file}")
    logger.info(f"Total proteins: {summary_stats['total_proteins']}")
    logger.info(f"Proteins with TM regions: {summary_stats['proteins_with_tm_regions']}")
    logger.info(f"Proteins extracellular only: {summary_stats['proteins_extracellular_only']}")
    
    return summary_stats

def main():
    """Main function to run topology prediction."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    topology_method = snakemake.params.topology_method
    topology_dir = snakemake.params.topology_dir
    
    # Input files - use the selected 3D paths files directly from Snakemake input
    selected_3d_paths_positive = snakemake.input.selected_3d_paths_positive
    selected_3d_paths_negative = snakemake.input.selected_3d_paths_negative
    
    # Output files
    topology_results = snakemake.output.topology_results
    topology_summary = snakemake.output.topology_summary
    
    logger.info(f"Starting topology prediction for {analysis}_{paramset}")
    logger.info(f"Using method: {topology_method}")
    
    # Create output directory
    os.makedirs(topology_dir, exist_ok=True)
    os.makedirs(os.path.dirname(topology_results), exist_ok=True)
    
    # Collect all 3D structure sequences for prediction
    all_sequences = {}
    
    # Process positive group
    logger.info("Processing positive group...")
    if os.path.exists(selected_3d_paths_positive):
        pos_gene_map = parse_selected_3d_paths(selected_3d_paths_positive)
        pos_sequences = extract_structure_sequences_from_paths(pos_gene_map)
        all_sequences.update(pos_sequences)
    else:
        logger.warning(f"Positive 3D paths file not found: {selected_3d_paths_positive}")
    
    # Process negative group  
    logger.info("Processing negative group...")
    if os.path.exists(selected_3d_paths_negative):
        neg_gene_map = parse_selected_3d_paths(selected_3d_paths_negative)
        neg_sequences = extract_structure_sequences_from_paths(neg_gene_map)
        all_sequences.update(neg_sequences)
    else:
        logger.warning(f"Negative 3D paths file not found: {selected_3d_paths_negative}")
    
    if not all_sequences:
        logger.error("No 3D structure sequences found for topology prediction")
        sys.exit(1)
    
    logger.info(f"Total 3D structure sequences for topology prediction: {len(all_sequences)}")
    for seq_id, seq in all_sequences.items():
        logger.info(f"  {seq_id}: {len(seq)} residues")
    
    # Run topology prediction using real DeepTMHMM
    success = run_deeptmhmm_biolib(all_sequences, topology_dir)
    
    if not success:
        logger.error("DeepTMHMM topology prediction failed")
        sys.exit(1)
    
    # Parse results
    topology_data = parse_deeptmhmm_output(topology_dir)
    
    if not topology_data:
        logger.error("No topology data could be parsed from DeepTMHMM output")
        sys.exit(1)
    
    # Save detailed results
    with open(topology_results, 'w') as f:
        json.dump(topology_data, f, indent=2)
    
    logger.info(f"Saved detailed topology results: {topology_results}")
    
    # Create summary
    summary_stats = create_topology_summary(topology_data, topology_summary)
    
    logger.info("Topology prediction completed successfully")
    logger.info(f"Summary: {summary_stats['proteins_with_tm_regions']}/{summary_stats['total_proteins']} proteins have transmembrane regions")

if __name__ == "__main__":
    main()