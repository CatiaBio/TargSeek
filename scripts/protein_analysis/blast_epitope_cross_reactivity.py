#!/usr/bin/env python3
"""
BLAST cross-reactivity analysis for filtered epitopes.

This script performs BLAST searches of topology-filtered epitopes against
multiple proteome databases to assess potential cross-reactivity with
host organisms and other microorganisms. This helps identify epitopes
that are specific to target pathogens.

Input:
- Topology-filtered epitope tables (CSV format)
- Extracellular mapping data
- BLAST database directory with formatted proteomes

Output:
- BLAST results for each epitope
- Cross-reactivity summary report
- Analysis of epitope specificity scores
"""

import os
import sys
import json
import pandas as pd
import subprocess
import logging
import glob
from pathlib import Path
from typing import Dict, List, Tuple, Any
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def find_epitope_tables(filtered_dir: str) -> List[str]:
    """Find all topology-filtered epitope table files."""
    pattern = os.path.join(filtered_dir, "*", "*_filtered_epitopes.csv")
    epitope_files = glob.glob(pattern)
    
    if not epitope_files:
        # Try alternative patterns
        pattern = os.path.join(filtered_dir, "*_filtered_epitopes.csv")
        epitope_files = glob.glob(pattern)
    
    logger.info(f"Found {len(epitope_files)} epitope table files")
    return epitope_files

def load_epitope_tables(epitope_files: List[str]) -> Dict[str, pd.DataFrame]:
    """Load all epitope tables into a dictionary."""
    epitope_data = {}
    
    for file_path in epitope_files:
        try:
            # Extract protein/structure info from filename
            basename = os.path.basename(file_path)
            protein_structure_id = basename.replace('_filtered_epitopes.csv', '')
            
            df = pd.read_csv(file_path)
            if not df.empty:
                epitope_data[protein_structure_id] = df
                logger.info(f"Loaded {len(df)} epitopes for {protein_structure_id}")
            
        except Exception as e:
            logger.error(f"Error loading epitope table {file_path}: {e}")
    
    return epitope_data

def get_available_blast_databases(blast_db_dir: str) -> List[Dict[str, str]]:
    """Get list of available BLAST databases."""
    databases = []
    
    # Look for formatted BLAST databases (.phr, .pin, .psq files)
    phr_files = glob.glob(os.path.join(blast_db_dir, "*.phr"))
    
    for phr_file in phr_files:
        db_name = os.path.basename(phr_file).replace('.phr', '')
        fasta_file = os.path.join(blast_db_dir, f"{db_name}.fasta")
        
        if os.path.exists(fasta_file):
            databases.append({
                'name': db_name,
                'path': os.path.join(blast_db_dir, db_name),
                'fasta_path': fasta_file
            })
    
    logger.info(f"Found {len(databases)} formatted BLAST databases")
    for db in databases[:10]:  # Log first 10
        logger.info(f"  - {db['name']}")
    
    return databases

def run_blast_search(epitope_seq: str, database_path: str, temp_dir: str) -> List[Dict]:
    """Run BLAST search for a single epitope sequence."""
    results = []
    
    try:
        # Create temporary files
        query_file = os.path.join(temp_dir, "epitope_query.fasta")
        output_file = os.path.join(temp_dir, "blast_results.xml")
        
        # Write epitope sequence to temporary FASTA file
        with open(query_file, 'w') as f:
            f.write(f">epitope\\n{epitope_seq}\\n")
        
        # Run BLAST search
        blast_cmd = [
            "blastp",
            "-query", query_file,
            "-db", database_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
            "-evalue", "1000",  # High evalue to catch weak matches
            "-word_size", "2",   # Lower word size for short sequences
            "-matrix", "PAM30",  # Better for short sequences
            "-comp_based_stats", "0",  # Turn off composition-based statistics
            "-max_target_seqs", "50"   # Limit number of hits
        ]
        
        result = subprocess.run(
            blast_cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode == 0 and result.stdout.strip():
            # Parse BLAST results
            lines = result.stdout.strip().split('\\n')
            for line in lines:
                if line.strip():
                    parts = line.split('\\t')
                    if len(parts) >= 12:
                        results.append({
                            'query_id': parts[0],
                            'subject_id': parts[1],
                            'percent_identity': float(parts[2]),
                            'alignment_length': int(parts[3]),
                            'mismatches': int(parts[4]),
                            'gap_opens': int(parts[5]),
                            'query_start': int(parts[6]),
                            'query_end': int(parts[7]),
                            'subject_start': int(parts[8]),
                            'subject_end': int(parts[9]),
                            'evalue': float(parts[10]),
                            'bit_score': float(parts[11]),
                            'subject_title': parts[12] if len(parts) > 12 else 'Unknown'
                        })
        
        elif result.returncode != 0:
            logger.warning(f"BLAST search failed: {result.stderr}")
    
    except subprocess.TimeoutExpired:
        logger.warning(f"BLAST search timed out for epitope: {epitope_seq[:20]}...")
    except Exception as e:
        logger.error(f"Error running BLAST search: {e}")
    
    # Clean up temporary files
    for temp_file in [query_file, output_file]:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return results

def analyze_epitope_specificity(epitope_data: Dict[str, pd.DataFrame], 
                               blast_databases: List[Dict[str, str]], 
                               output_dir: str) -> Dict[str, Any]:
    """Analyze epitope specificity using BLAST searches."""
    
    specificity_results = {
        'epitope_analysis': {},
        'database_summary': {},
        'specificity_scores': {},
        'cross_reactive_epitopes': [],
        'highly_specific_epitopes': [],
        'analysis_parameters': {
            'min_identity_threshold': 70.0,
            'min_coverage_threshold': 0.8,
            'databases_searched': [db['name'] for db in blast_databases]
        }
    }
    
    # Create temporary directory for BLAST operations
    with tempfile.TemporaryDirectory() as temp_dir:
        logger.info(f"Starting BLAST analysis with temporary directory: {temp_dir}")
        
        epitope_count = 0
        total_epitopes = sum(len(df) for df in epitope_data.values())
        
        for protein_id, epitope_df in epitope_data.items():
            logger.info(f"Processing epitopes for {protein_id} ({len(epitope_df)} epitopes)")
            
            protein_results = {
                'epitopes': {},
                'summary': {
                    'total_epitopes': len(epitope_df),
                    'cross_reactive_count': 0,
                    'highly_specific_count': 0,
                    'average_specificity_score': 0.0
                }
            }
            
            for idx, epitope_row in epitope_df.iterrows():
                epitope_count += 1
                epitope_seq = str(epitope_row.get('Epitope_Sequence', epitope_row.get('sequence', 'Unknown')))
                epitope_start = epitope_row.get('Epitope_Start', epitope_row.get('start', 0))
                epitope_end = epitope_row.get('Epitope_End', epitope_row.get('end', 0))
                
                if len(epitope_seq) < 6:  # Skip very short sequences
                    continue
                
                logger.info(f"BLAST searching epitope {epitope_count}/{total_epitopes}: {epitope_seq[:20]}...")
                
                epitope_blast_results = {}
                cross_reactive_hits = 0
                total_hits = 0
                
                # Search against each database
                for db_info in blast_databases:
                    db_name = db_info['name']
                    db_path = db_info['path']
                    
                    # Skip microbial databases for now (focus on host organisms)
                    host_dbs = ['human_proteome', 'cow_proteome', 'mouse_proteome', 'pig_proteome', 
                               'chicken_proteome', 'sheep_proteome', 'goat_proteome', 'horse_proteome']
                    
                    if db_name not in host_dbs:
                        continue
                    
                    blast_hits = run_blast_search(epitope_seq, db_path, temp_dir)
                    
                    if blast_hits:
                        # Filter significant hits
                        significant_hits = [
                            hit for hit in blast_hits
                            if hit['percent_identity'] >= 70.0 and 
                               hit['alignment_length'] / len(epitope_seq) >= 0.8
                        ]
                        
                        epitope_blast_results[db_name] = {
                            'total_hits': len(blast_hits),
                            'significant_hits': len(significant_hits),
                            'best_hit': blast_hits[0] if blast_hits else None,
                            'hits': significant_hits[:10]  # Keep top 10 significant hits
                        }
                        
                        cross_reactive_hits += len(significant_hits)
                        total_hits += len(blast_hits)
                
                # Calculate specificity score (lower = more specific)
                specificity_score = cross_reactive_hits / max(len(host_dbs), 1)
                
                epitope_analysis = {
                    'sequence': epitope_seq,
                    'start': epitope_start,
                    'end': epitope_end,
                    'length': len(epitope_seq),
                    'blast_results': epitope_blast_results,
                    'cross_reactive_hits': cross_reactive_hits,
                    'total_blast_hits': total_hits,
                    'specificity_score': specificity_score,
                    'classification': 'highly_specific' if specificity_score < 0.1 else 
                                   'moderately_specific' if specificity_score < 0.3 else 'cross_reactive'
                }
                
                protein_results['epitopes'][f"epitope_{idx}"] = epitope_analysis
                
                # Update summary statistics
                if specificity_score < 0.1:
                    protein_results['summary']['highly_specific_count'] += 1
                    specificity_results['highly_specific_epitopes'].append({
                        'protein': protein_id,
                        'epitope': epitope_seq,
                        'score': specificity_score
                    })
                elif specificity_score > 0.3:
                    protein_results['summary']['cross_reactive_count'] += 1
                    specificity_results['cross_reactive_epitopes'].append({
                        'protein': protein_id,
                        'epitope': epitope_seq,
                        'score': specificity_score
                    })
                
                time.sleep(0.1)  # Small delay between searches
            
            # Calculate average specificity score for protein
            epitope_scores = [epi['specificity_score'] for epi in protein_results['epitopes'].values()]
            if epitope_scores:
                protein_results['summary']['average_specificity_score'] = sum(epitope_scores) / len(epitope_scores)
            
            specificity_results['epitope_analysis'][protein_id] = protein_results
    
    # Generate database summary
    for db_info in blast_databases:
        specificity_results['database_summary'][db_info['name']] = {
            'database_path': db_info['path'],
            'fasta_path': db_info['fasta_path'],
            'searched': db_info['name'] in ['human_proteome', 'cow_proteome', 'mouse_proteome']
        }
    
    logger.info(f"BLAST analysis completed for {epitope_count} epitopes")
    logger.info(f"Highly specific epitopes: {len(specificity_results['highly_specific_epitopes'])}")
    logger.info(f"Cross-reactive epitopes: {len(specificity_results['cross_reactive_epitopes'])}")
    
    return specificity_results

def save_blast_results(specificity_results: Dict[str, Any], output_dir: str, 
                      blast_summary_file: str) -> None:
    """Save BLAST analysis results to files."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save comprehensive results
    detailed_results_file = os.path.join(output_dir, "detailed_blast_analysis.json")
    with open(detailed_results_file, 'w') as f:
        json.dump(specificity_results, f, indent=2)
    
    # Save summary report
    with open(blast_summary_file, 'w') as f:
        summary = {
            'analysis_summary': {
                'total_proteins_analyzed': len(specificity_results['epitope_analysis']),
                'total_epitopes_analyzed': sum(
                    protein_data['summary']['total_epitopes'] 
                    for protein_data in specificity_results['epitope_analysis'].values()
                ),
                'highly_specific_epitopes': len(specificity_results['highly_specific_epitopes']),
                'cross_reactive_epitopes': len(specificity_results['cross_reactive_epitopes']),
                'databases_searched': specificity_results['analysis_parameters']['databases_searched']
            },
            'protein_summaries': {
                protein_id: protein_data['summary']
                for protein_id, protein_data in specificity_results['epitope_analysis'].items()
            },
            'top_specific_epitopes': specificity_results['highly_specific_epitopes'][:20],
            'most_cross_reactive_epitopes': specificity_results['cross_reactive_epitopes'][:20]
        }
        json.dump(summary, f, indent=2)
    
    # Create CSV report for easy analysis
    csv_report_file = os.path.join(output_dir, "epitope_specificity_report.csv")
    epitope_rows = []
    
    for protein_id, protein_data in specificity_results['epitope_analysis'].items():
        for epitope_id, epitope_data in protein_data['epitopes'].items():
            epitope_rows.append({
                'Protein_ID': protein_id,
                'Epitope_ID': epitope_id,
                'Epitope_Sequence': epitope_data['sequence'],
                'Epitope_Length': epitope_data['length'],
                'Epitope_Start': epitope_data['start'],
                'Epitope_End': epitope_data['end'],
                'Cross_Reactive_Hits': epitope_data['cross_reactive_hits'],
                'Total_BLAST_Hits': epitope_data['total_blast_hits'],
                'Specificity_Score': epitope_data['specificity_score'],
                'Classification': epitope_data['classification']
            })
    
    if epitope_rows:
        epitope_df = pd.DataFrame(epitope_rows)
        epitope_df.to_csv(csv_report_file, index=False)
        logger.info(f"Saved epitope specificity report: {csv_report_file}")
    
    logger.info(f"Saved BLAST analysis results to {output_dir}")

def main():
    """Main function to run BLAST cross-reactivity analysis."""
    
    # Get Snakemake parameters
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    group = snakemake.params.group
    blast_databases_dir = snakemake.params.blast_databases_dir
    output_dir = snakemake.params.output_dir
    
    # Input files
    filtered_epitopes_sentinel = snakemake.input.filtered_epitopes_sentinel
    
    # Output files
    blast_results_sentinel = snakemake.output.blast_results_sentinel
    blast_summary = snakemake.output.blast_summary
    
    logger.info(f"Starting BLAST cross-reactivity analysis for {analysis}_{paramset} gram_{group}")
    
    # Find filtered epitope directory
    filtered_dir = os.path.dirname(filtered_epitopes_sentinel)
    if not os.path.exists(filtered_dir):
        logger.error(f"Filtered epitopes directory not found: {filtered_dir}")
        sys.exit(1)
    
    # Find epitope table files
    epitope_files = find_epitope_tables(filtered_dir)
    if not epitope_files:
        logger.error(f"No filtered epitope files found in {filtered_dir}")
        sys.exit(1)
    
    # Load epitope data
    epitope_data = load_epitope_tables(epitope_files)
    if not epitope_data:
        logger.error("No valid epitope data loaded")
        sys.exit(1)
    
    # Get available BLAST databases
    blast_databases = get_available_blast_databases(blast_databases_dir)
    if not blast_databases:
        logger.error(f"No formatted BLAST databases found in {blast_databases_dir}")
        sys.exit(1)
    
    # Run BLAST analysis
    logger.info("Starting BLAST specificity analysis...")
    specificity_results = analyze_epitope_specificity(epitope_data, blast_databases, output_dir)
    
    # Save results
    save_blast_results(specificity_results, output_dir, blast_summary)
    
    # Create sentinel file
    os.makedirs(os.path.dirname(blast_results_sentinel), exist_ok=True)
    with open(blast_results_sentinel, 'w') as f:
        f.write(f"BLAST cross-reactivity analysis completed for {analysis}_{paramset} gram_{group}\\n")
        f.write(f"Analyzed {len(epitope_data)} proteins\\n")
        f.write(f"Total epitopes: {sum(len(df) for df in epitope_data.values())}\\n")
        f.write(f"Databases searched: {len(blast_databases)}\\n")
    
    logger.info("BLAST cross-reactivity analysis completed successfully")

if __name__ == "__main__":
    main()