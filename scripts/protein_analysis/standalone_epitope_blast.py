#!/usr/bin/env python3
"""
Standalone BLAST cross-reactivity analysis for epitopes.

This script performs BLAST searches of provided epitope sequences against
multiple proteome databases to assess potential cross-reactivity with
host organisms and other microorganisms. This helps identify epitopes
that are specific to target pathogens.

Usage:
    python standalone_epitope_blast.py "EPITOPE1,EPITOPE2,EPITOPE3"

Input:
- Comma-separated epitope sequences as string argument
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
import argparse
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

def parse_epitope_string(epitope_string: str) -> Dict[str, pd.DataFrame]:
    """Parse comma-separated epitope string into format expected by the analysis functions."""
    epitopes = [epitope.strip().upper() for epitope in epitope_string.split(',')]
    epitopes = [epitope for epitope in epitopes if epitope]  # Remove empty strings

    # Validate epitopes (basic check for amino acid sequences)
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    valid_epitopes = []

    for epitope in epitopes:
        if all(aa in valid_aa for aa in epitope):
            valid_epitopes.append(epitope)
            logger.info(f"Valid epitope: {epitope} ({len(epitope)} aa)")
        else:
            logger.warning(f"Invalid epitope (contains non-amino acid characters): {epitope}")

    if not valid_epitopes:
        logger.error("No valid epitopes provided")
        return {}

    # Create a DataFrame in the format expected by the analysis functions
    epitope_data = {}
    for i, epitope in enumerate(valid_epitopes):
        epitope_id = f"epitope_{i+1}"

        # Create a DataFrame with the epitope data
        df = pd.DataFrame({
            'Epitope_ID': [epitope_id],
            'Epitope_Sequence': [epitope],
            'Epitope_Start': [1],  # Dummy start position
            'Epitope_End': [len(epitope)],  # Dummy end position
            'Epitope_Length': [len(epitope)],
            'BepiPred_Score': [0.5],  # Dummy score
            'Extracellular': [True]  # Assume all are extracellular
        })

        epitope_data[epitope_id] = df

    logger.info(f"Parsed {len(valid_epitopes)} valid epitopes")
    return epitope_data

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

    # Look for formatted BLAST databases (.fasta.phr files)
    phr_files = glob.glob(os.path.join(blast_db_dir, "*.fasta.phr"))

    for phr_file in phr_files:
        # Extract database name (remove .fasta.phr extension)
        db_name = os.path.basename(phr_file).replace('.fasta.phr', '')
        fasta_file = os.path.join(blast_db_dir, f"{db_name}.fasta")
        db_path = os.path.join(blast_db_dir, f"{db_name}.fasta")

        if os.path.exists(fasta_file):
            databases.append({
                'name': db_name,
                'path': db_path,  # Use .fasta as the database path for BLAST
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
            f.write(f">epitope\n{epitope_seq}\n")
        
        # Run BLAST search
        blast_cmd = [
            "blastp",
            "-query", query_file,
            "-db", database_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
            "-evalue", "1000",  # High evalue to catch weak matches
            "-word_size", "2",   # Lower word size for short sequences
            "-matrix", "PAM30",  # Better for short sequences
            "-comp_based_stats", "0"  # Turn off composition-based statistics
            # Removed max_target_seqs to get all hits
        ]
        
        result = subprocess.run(
            blast_cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode == 0 and result.stdout.strip():
            # Parse BLAST results
            lines = result.stdout.strip().split('\n')
            for line in lines:
                if line.strip():
                    parts = line.split('\t')
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

                    # Search all databases (removed filtering to include all species)
                    # host_dbs = ['human_proteome', 'cow_proteome', 'mouse_proteome', 'pig_proteome',
                    #            'chicken_proteome', 'sheep_proteome', 'goat_proteome', 'horse_proteome']

                    # if db_name not in host_dbs:
                    #     continue
                    
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
                            'significant_hits': significant_hits,  # Keep all significant hits
                            'significant_count': len(significant_hits),
                            'best_hit': blast_hits[0] if blast_hits else None,
                            'all_hits': blast_hits  # Keep all hits
                        }
                        
                        cross_reactive_hits += len(significant_hits)
                        total_hits += len(blast_hits)
                
                # Calculate specificity score (lower = more specific)
                specificity_score = cross_reactive_hits / max(len(blast_databases), 1)
                
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
    detailed_csv_file = os.path.join(output_dir, "epitope_detailed_cross_reactivity.csv")

    epitope_rows = []
    detailed_rows = []

    for protein_id, protein_data in specificity_results['epitope_analysis'].items():
        for epitope_id, epitope_data in protein_data['epitopes'].items():
            # Basic summary row
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

            # Detailed species-specific rows
            blast_results = epitope_data.get('blast_results', {})
            for db_name, db_results in blast_results.items():
                if db_results and db_results.get('significant_hits'):
                    species_name = db_name.replace('_proteome', '').replace('_proteins', '').title()

                    # Get best hit for this species
                    best_hit = db_results.get('best_hit', {})
                    significant_hits = db_results.get('significant_hits', [])

                    # Calculate statistics for this species
                    if significant_hits:
                        identities = [hit.get('percent_identity', 0) for hit in significant_hits]
                        evalues = [hit.get('evalue', 1.0) for hit in significant_hits]
                        bitscores = [hit.get('bitscore', 0) for hit in significant_hits]

                        detailed_rows.append({
                            'Epitope_Sequence': epitope_data['sequence'],
                            'Species': species_name,
                            'Hit_Count': len(significant_hits),
                            'Best_Identity_%': max(identities) if identities else 0,
                            'Best_E_value': min(evalues) if evalues else 1.0,
                            'Best_Bitscore': max(bitscores) if bitscores else 0,
                            'Avg_Identity_%': sum(identities) / len(identities) if identities else 0,
                            'Min_E_value': min(evalues) if evalues else 1.0,
                            'Max_Bitscore': max(bitscores) if bitscores else 0,
                            'Best_Hit_ID': best_hit.get('subject_id', 'N/A'),
                            'Classification': epitope_data['classification']
                        })

    # Save basic summary
    if epitope_rows:
        epitope_df = pd.DataFrame(epitope_rows)
        epitope_df.to_csv(csv_report_file, index=False)
        logger.info(f"Saved epitope specificity report: {csv_report_file}")

    # Save detailed species breakdown
    if detailed_rows:
        detailed_df = pd.DataFrame(detailed_rows)
        # Sort by epitope and hit count
        detailed_df = detailed_df.sort_values(['Epitope_Sequence', 'Hit_Count'], ascending=[True, False])
        detailed_df.to_csv(detailed_csv_file, index=False)
        logger.info(f"Saved detailed cross-reactivity report: {detailed_csv_file}")
    
    logger.info(f"Saved BLAST analysis results to {output_dir}")

def main():
    """Main function to run standalone BLAST cross-reactivity analysis."""

    parser = argparse.ArgumentParser(
        description="Standalone BLAST cross-reactivity analysis for epitopes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s "DFQADDADLSDYTN,LSNMQPQVAMWRYLYSMGEHPSTSDQDNSFK,GKEMPFY"
  %(prog)s "EPITOPE1,EPITOPE2" --output-dir results
  %(prog)s "PHQASNYDPDYDYESATQDGAKDLSKSD,NWDSSQYSGYPDYSDP" --blast-db data/blast_databases
        """
    )

    parser.add_argument('epitopes',
                       help='Epitope sequences (comma-separated)')

    parser.add_argument('--blast-db', '-b',
                       default='data/blast_databases',
                       help='BLAST databases directory (default: data/blast_databases)')

    parser.add_argument('--output-dir', '-o',
                       default='epitope_blast_results',
                       help='Output directory for results (default: epitope_blast_results)')

    args = parser.parse_args()

    # Parameters
    blast_databases_dir = args.blast_db
    output_dir = args.output_dir

    logger.info(f"Starting standalone BLAST cross-reactivity analysis")
    logger.info(f"BLAST databases directory: {blast_databases_dir}")
    logger.info(f"Output directory: {output_dir}")

    # Parse epitope string
    epitope_data = parse_epitope_string(args.epitopes)
    if not epitope_data:
        logger.error("No valid epitopes provided")
        sys.exit(1)

    # Get available BLAST databases
    blast_databases = get_available_blast_databases(blast_databases_dir)
    if not blast_databases:
        logger.error(f"No formatted BLAST databases found in {blast_databases_dir}")
        sys.exit(1)

    # Run BLAST analysis
    logger.info("Starting BLAST specificity analysis...")
    specificity_results = analyze_epitope_specificity(epitope_data, blast_databases, output_dir)

    # Generate output file names
    blast_summary = os.path.join(output_dir, "blast_summary.json")

    # Save results
    save_blast_results(specificity_results, output_dir, blast_summary)

    logger.info("Standalone BLAST cross-reactivity analysis completed successfully")
    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"Summary report: {blast_summary}")

if __name__ == "__main__":
    main()