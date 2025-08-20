#!/usr/bin/env python3
"""
Filter epitopes to keep only those in extracellular regions.

This script filters BepiPred epitope predictions to retain only epitopes that are
predicted to be in extracellular regions based on topology analysis. This ensures
that only accessible epitopes are considered for downstream analysis.

Input:
- Epitope tables from BepiPred predictions
- Extracellular region mapping from topology analysis
- PDB numbering mapping
- Configuration parameters for filtering

Output:
- Filtered epitope tables (only extracellular epitopes)
- Topology filter report with statistics
- Updated epitope conservation analysis inputs
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path
import glob
import logging
from typing import Dict, List, Tuple, Any, Set

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_extracellular_mapping(mapping_file: str) -> Dict[str, Dict]:
    """Load extracellular region mapping from JSON file."""
    try:
        with open(mapping_file, 'r') as f:
            mapping_data = json.load(f)
        
        # Extract the actual mapping from the JSON structure
        if 'mapping' in mapping_data:
            protein_mappings = mapping_data['mapping']
        else:
            protein_mappings = mapping_data
        
        logger.info(f"Loaded extracellular mapping for {len(protein_mappings)} proteins")
        return protein_mappings
    except Exception as e:
        logger.error(f"Error loading extracellular mapping: {e}")
        return {}

def find_epitope_files(epitope_dir: str) -> List[str]:
    """Find all epitope table files in the BepiPred predictions directory."""
    epitope_files = []
    
    # Look for epitope table files in various locations
    search_patterns = [
        os.path.join(epitope_dir, "**", "*epitope_table*.csv"),
        os.path.join(epitope_dir, "**", "*epitope_table*.tsv"),
        os.path.join(epitope_dir, "**", "*epitopes*.csv"),
        os.path.join(epitope_dir, "**", "*epitopes*.tsv"),
        os.path.join(epitope_dir, "**", "*linear_epitopes*.tsv"),
        os.path.join(epitope_dir, "**", "*linear_epitopes*.csv")
    ]
    
    for pattern in search_patterns:
        epitope_files.extend(glob.glob(pattern, recursive=True))
    
    # Remove duplicates
    epitope_files = list(set(epitope_files))
    
    logger.info(f"Found {len(epitope_files)} epitope table files")
    return epitope_files

def load_epitope_table(file_path: str) -> pd.DataFrame:
    """Load epitope table from CSV/TSV file."""
    try:
        # Try different separators and handle comment lines
        for sep in ['\t', ',']:
            try:
                df = pd.read_csv(file_path, sep=sep, comment='#')
                if len(df.columns) > 1 and len(df) > 0:  # Valid dataframe with data
                    logger.debug(f"Loaded epitope table with {len(df)} rows from {file_path}")
                    return df
            except Exception as e:
                logger.debug(f"Failed to load with separator '{sep}': {e}")
                continue
        
        logger.warning(f"Could not load epitope table: {file_path}")
        return pd.DataFrame()
    
    except Exception as e:
        logger.error(f"Error loading epitope table {file_path}: {e}")
        return pd.DataFrame()

def extract_protein_structure_info(file_path: str) -> Tuple[str, str]:
    """Extract protein and structure information from epitope file path."""
    
    # Extract from file name (expected format: protein_structure_epitope_table.csv)
    file_name = os.path.basename(file_path)
    
    # Remove file extension
    base_name = os.path.splitext(file_name)[0]
    
    # Try to parse protein and structure
    if '_epitope_table' in base_name:
        protein_structure = base_name.replace('_epitope_table', '')
    elif '_epitopes' in base_name:
        protein_structure = base_name.replace('_epitopes', '')
    else:
        protein_structure = base_name
    
    # Try to split into protein and structure
    if '_' in protein_structure:
        parts = protein_structure.split('_')
        protein = parts[0]
        structure = '_'.join(parts[1:])
    else:
        protein = protein_structure
        structure = ""
    
    return protein, structure

def check_epitope_extracellular_overlap(
    epitope_start: int, 
    epitope_end: int, 
    extracellular_regions: List[Tuple[int, int]],
    min_coverage: float = 0.5
) -> Tuple[bool, float]:
    """Check if epitope overlaps with extracellular regions and calculate coverage."""
    
    epitope_length = epitope_end - epitope_start + 1
    overlap_length = 0
    
    for ext_start, ext_end in extracellular_regions:
        # Calculate overlap
        overlap_start = max(epitope_start, ext_start)
        overlap_end = min(epitope_end, ext_end)
        
        if overlap_start <= overlap_end:
            overlap_length += overlap_end - overlap_start + 1
    
    coverage = overlap_length / epitope_length if epitope_length > 0 else 0.0
    is_extracellular = coverage >= min_coverage
    
    return is_extracellular, coverage

def filter_epitopes_by_topology(
    epitope_df: pd.DataFrame,
    protein_id: str,
    extracellular_mapping: Dict[str, Dict],
    min_extracellular_coverage: float = 0.5
) -> Tuple[pd.DataFrame, Dict]:
    """Filter epitopes based on extracellular region overlap."""
    
    if protein_id not in extracellular_mapping:
        logger.warning(f"No extracellular mapping found for {protein_id}")
        return pd.DataFrame(), {
            'total_epitopes': len(epitope_df),
            'extracellular_epitopes': 0,
            'filtered_epitopes': 0,
            'coverage_stats': {}
        }
    
    protein_mapping = extracellular_mapping[protein_id]
    extracellular_spans = protein_mapping.get('sequence_extracellular_spans', [])
    
    if not extracellular_spans:
        logger.warning(f"No extracellular spans found for {protein_id}")
        return pd.DataFrame(), {
            'total_epitopes': len(epitope_df),
            'extracellular_epitopes': 0,
            'filtered_epitopes': 0,
            'coverage_stats': {}
        }
    
    # Create filtered dataframe
    filtered_epitopes = []
    coverage_stats = {
        'coverage_distribution': [],
        'average_coverage': 0.0,
        'epitopes_by_coverage': {
            '0.0-0.25': 0,
            '0.25-0.5': 0,
            '0.5-0.75': 0,
            '0.75-1.0': 0
        }
    }
    
    total_coverage = 0.0
    
    for idx, row in epitope_df.iterrows():
        # Get epitope coordinates (adapt to your epitope table format)
        epitope_start = row.get('start', row.get('Start', row.get('position', 1)))
        epitope_end = row.get('end', row.get('End', epitope_start + row.get('length', 1) - 1))
        
        # Check overlap with extracellular regions
        is_extracellular, coverage = check_epitope_extracellular_overlap(
            epitope_start, epitope_end, extracellular_spans, min_extracellular_coverage
        )
        
        # Add coverage information to row
        row_dict = row.to_dict()
        row_dict['extracellular_coverage'] = coverage
        row_dict['is_extracellular'] = is_extracellular
        
        if is_extracellular:
            filtered_epitopes.append(row_dict)
        
        # Update statistics
        total_coverage += coverage
        coverage_stats['coverage_distribution'].append(coverage)
        
        # Bin coverage
        if coverage < 0.25:
            coverage_stats['epitopes_by_coverage']['0.0-0.25'] += 1
        elif coverage < 0.5:
            coverage_stats['epitopes_by_coverage']['0.25-0.5'] += 1
        elif coverage < 0.75:
            coverage_stats['epitopes_by_coverage']['0.5-0.75'] += 1
        else:
            coverage_stats['epitopes_by_coverage']['0.75-1.0'] += 1
    
    # Calculate average coverage
    if len(epitope_df) > 0:
        coverage_stats['average_coverage'] = total_coverage / len(epitope_df)
    
    filtered_df = pd.DataFrame(filtered_epitopes) if filtered_epitopes else pd.DataFrame()
    
    filter_stats = {
        'total_epitopes': len(epitope_df),
        'extracellular_epitopes': len(filtered_df),
        'filtered_epitopes': len(epitope_df) - len(filtered_df),
        'coverage_stats': coverage_stats
    }
    
    return filtered_df, filter_stats

def process_all_epitope_files(
    epitope_files: List[str],
    extracellular_mapping: Dict[str, Dict],
    output_dir: str,
    min_extracellular_coverage: float = 0.5
) -> Dict[str, Any]:
    """Process all epitope files and create filtered versions."""
    
    # Create output directory for filtered epitopes
    filtered_dir = os.path.join(output_dir, "filtered_epitope_tables")
    os.makedirs(filtered_dir, exist_ok=True)
    
    processing_stats = {
        'total_files_processed': 0,
        'files_with_epitopes': 0,
        'total_epitopes': 0,
        'total_extracellular_epitopes': 0,
        'protein_stats': {},
        'overall_coverage': 0.0
    }
    
    total_coverage = 0.0
    
    for epitope_file in epitope_files:
        try:
            # Load epitope table
            epitope_df = load_epitope_table(epitope_file)
            
            if epitope_df.empty:
                continue
            
            # Extract protein and structure info
            protein, structure = extract_protein_structure_info(epitope_file)
            protein_id = f"{protein}_{structure}" if structure else protein
            
            # Filter epitopes
            filtered_df, filter_stats = filter_epitopes_by_topology(
                epitope_df, protein_id, extracellular_mapping, min_extracellular_coverage
            )
            
            # Save filtered epitopes
            if not filtered_df.empty:
                output_file = os.path.join(filtered_dir, f"{protein_id}_filtered_epitopes.csv")
                filtered_df.to_csv(output_file, index=False)
                logger.info(f"Saved {len(filtered_df)} filtered epitopes for {protein_id}")
            
            # Update statistics
            processing_stats['total_files_processed'] += 1
            if filter_stats['total_epitopes'] > 0:
                processing_stats['files_with_epitopes'] += 1
            
            processing_stats['total_epitopes'] += filter_stats['total_epitopes']
            processing_stats['total_extracellular_epitopes'] += filter_stats['extracellular_epitopes']
            processing_stats['protein_stats'][protein_id] = filter_stats
            
            total_coverage += filter_stats['coverage_stats']['average_coverage']
            
        except Exception as e:
            logger.error(f"Error processing epitope file {epitope_file}: {e}")
    
    # Calculate overall coverage
    if processing_stats['files_with_epitopes'] > 0:
        processing_stats['overall_coverage'] = total_coverage / processing_stats['files_with_epitopes']
    
    return processing_stats

def create_topology_filter_report(
    processing_stats: Dict[str, Any],
    extracellular_mapping: Dict[str, Dict],
    output_file: str
):
    """Create comprehensive topology filter report."""
    
    report = {
        'summary': {
            'total_proteins_analyzed': len(processing_stats['protein_stats']),
            'total_epitopes_input': processing_stats['total_epitopes'],
            'total_extracellular_epitopes': processing_stats['total_extracellular_epitopes'],
            'filtration_rate': 0.0,
            'average_extracellular_coverage': processing_stats['overall_coverage']
        },
        'protein_details': processing_stats['protein_stats'],
        'topology_summary': {
            'proteins_with_extracellular_regions': len(extracellular_mapping),
            'extracellular_region_stats': {}
        },
        'filtering_parameters': {
            'min_extracellular_coverage': snakemake.params.min_extracellular_coverage
        }
    }
    
    # Calculate filtration rate
    if processing_stats['total_epitopes'] > 0:
        report['summary']['filtration_rate'] = (
            processing_stats['total_extracellular_epitopes'] / processing_stats['total_epitopes']
        )
    
    # Add extracellular region statistics
    total_extracellular_residues = 0
    proteins_with_regions = 0
    
    for protein_id, mapping in extracellular_mapping.items():
        # Use 'extracellular_length' as it exists in the JSON structure
        extracellular_length = mapping.get('extracellular_length', 0)
        if extracellular_length > 0:
            proteins_with_regions += 1
            total_extracellular_residues += extracellular_length
    
    if proteins_with_regions > 0:
        report['topology_summary']['extracellular_region_stats'] = {
            'average_extracellular_residues_per_protein': total_extracellular_residues / proteins_with_regions,
            'proteins_with_extracellular_regions': proteins_with_regions
        }
    
    # Save report
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Created topology filter report: {output_file}")
    logger.info(f"Filtration rate: {report['summary']['filtration_rate']:.2%}")
    logger.info(f"Extracellular epitopes: {report['summary']['total_extracellular_epitopes']}/{report['summary']['total_epitopes_input']}")

def main():
    """Main function to filter epitopes by topology."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    topology_dir = snakemake.params.topology_dir
    min_extracellular_coverage = snakemake.params.min_extracellular_coverage
    
    # Input files
    epitope_tables_sentinel = snakemake.input.epitope_tables_sentinel
    extracellular_mapping_file = snakemake.input.extracellular_mapping
    pdb_mapping_file = snakemake.input.pdb_numbering_mapping
    
    # Output files
    filtered_epitopes_sentinel = snakemake.output.filtered_epitopes_sentinel
    topology_filter_report = snakemake.output.topology_filter_report
    
    logger.info(f"Starting epitope topology filtering for {analysis}_{paramset}")
    logger.info(f"Minimum extracellular coverage: {min_extracellular_coverage}")
    
    # Create output directory
    os.makedirs(os.path.dirname(filtered_epitopes_sentinel), exist_ok=True)
    
    # Load extracellular mapping
    extracellular_mapping = load_extracellular_mapping(extracellular_mapping_file)
    
    if not extracellular_mapping:
        logger.error("No extracellular mapping data loaded")
        sys.exit(1)
    
    # Find epitope table files
    epitope_dir = os.path.dirname(epitope_tables_sentinel)
    epitope_files = find_epitope_files(epitope_dir)
    
    if not epitope_files:
        logger.warning("No epitope table files found")
        # Create empty output files
        with open(filtered_epitopes_sentinel, 'w') as f:
            f.write("No epitope tables found for filtering\n")
        
        empty_report = {
            'summary': {'total_epitopes_input': 0, 'total_extracellular_epitopes': 0},
            'protein_details': {},
            'topology_summary': {}
        }
        with open(topology_filter_report, 'w') as f:
            json.dump(empty_report, f, indent=2)
        return
    
    # Process all epitope files
    processing_stats = process_all_epitope_files(
        epitope_files, extracellular_mapping, topology_dir, min_extracellular_coverage
    )
    
    # Create filter report
    create_topology_filter_report(
        processing_stats, extracellular_mapping, topology_filter_report
    )
    
    # Create sentinel file
    with open(filtered_epitopes_sentinel, 'w') as f:
        f.write(f"Epitope topology filtering completed successfully\n")
        f.write(f"Analysis: {analysis}_{paramset}\n")
        f.write(f"Total epitopes processed: {processing_stats['total_epitopes']}\n")
        f.write(f"Extracellular epitopes retained: {processing_stats['total_extracellular_epitopes']}\n")
        f.write(f"Filtration rate: {processing_stats['total_extracellular_epitopes']/max(processing_stats['total_epitopes'],1):.2%}\n")
    
    logger.info("Epitope topology filtering completed successfully")

if __name__ == "__main__":
    main()