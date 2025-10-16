#!/usr/bin/env python3
"""
Epitope Cross-Reactivity Analysis using BLAST
============================================

Standalone script to check cross-reactivity of epitopes against various databases
using BLAST to identify potential off-target matches.

Usage:
    python epitope_cross_reactivity_blast.py "PEPTIDE1,PEPTIDE2,PEPTIDE3"
    python epitope_cross_reactivity_blast.py "ALKYLPHENYLALANINE" --databases human,viral
    python epitope_cross_reactivity_blast.py "PEPTIDE" --evalue 1.0 --word-size 3

Features:
- Multiple epitope analysis (comma-separated)
- Multiple BLAST database support
- Configurable BLAST parameters
- Detailed hit analysis with organism classification
- Comprehensive summary table output
- Risk classification system
"""

import argparse
import subprocess
import tempfile
import os
import pandas as pd
from pathlib import Path
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class EpitopeCrossReactivityAnalyzer:
    """BLAST-based epitope cross-reactivity analyzer"""

    def __init__(self):
        # Database descriptions based on available local databases
        self.available_databases = {
            'human': 'Human proteins (Homo sapiens)',
            'viral': 'Viral proteins',
            'mouse': 'Mouse proteins (Mus musculus)',
            'cow': 'Cow proteins (Bos taurus)',
            'pig': 'Pig proteins (Sus scrofa)',
            'chicken': 'Chicken proteins (Gallus gallus)',
            'ecoli': 'E. coli proteins (Escherichia coli)',
            'salmonella': 'Salmonella proteins',
            'staph': 'Staphylococcus proteins',
            'fungal': 'Fungal proteins',
            'rat': 'Rat proteins (Rattus norvegicus)',
            'horse': 'Horse proteins (Equus caballus)',
            'zebrafish': 'Zebrafish proteins (Danio rerio)',
            'arabidopsis': 'Arabidopsis proteins (Arabidopsis thaliana)',
            'buffalo': 'Buffalo proteins (Bubalus bubalis)',
            'goat': 'Goat proteins (Capra hircus)',
            'rice': 'Rice proteins (Oryza sativa)',
            'salmon': 'Salmon proteins (Salmo salar)',
            'sheep': 'Sheep proteins (Ovis aries)',
            'soybean': 'Soybean proteins (Glycine max)'
        }

        # Local database paths in data folder
        self.database_paths = {
            'human': 'data/blast_databases/human_proteome.fasta',
            'viral': 'data/blast_databases/viral_proteins.fasta',
            'mouse': 'data/blast_databases/mouse_proteome.fasta',
            'cow': 'data/blast_databases/cow_proteome.fasta',
            'pig': 'data/blast_databases/pig_proteome.fasta',
            'chicken': 'data/blast_databases/chicken_proteome.fasta',
            'ecoli': 'data/blast_databases/ecoli_proteome.fasta',
            'salmonella': 'data/blast_databases/salmonella_proteome.fasta',
            'staph': 'data/blast_databases/staph_proteome.fasta',
            'fungal': 'data/blast_databases/fungal_proteins.fasta',
            'rat': 'data/blast_databases/rat_proteome.fasta',
            'horse': 'data/blast_databases/horse_proteome.fasta',
            'zebrafish': 'data/blast_databases/zebrafish_proteome.fasta',
            'arabidopsis': 'data/blast_databases/arabidopsis_proteome.fasta',
            'buffalo': 'data/blast_databases/buffalo_proteome.fasta',
            'goat': 'data/blast_databases/goat_proteome.fasta',
            'rice': 'data/blast_databases/rice_proteome.fasta',
            'salmon': 'data/blast_databases/salmon_proteome.fasta',
            'sheep': 'data/blast_databases/sheep_proteome.fasta',
            'soybean': 'data/blast_databases/soybean_proteome.fasta'
        }

    def parse_epitopes(self, epitope_string):
        """Parse comma-separated epitope string"""
        epitopes = [epitope.strip().upper() for epitope in epitope_string.split(',')]
        epitopes = [epitope for epitope in epitopes if epitope]  # Remove empty strings

        # Validate epitopes (basic check for amino acid sequences)
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        valid_epitopes = []

        for epitope in epitopes:
            if all(aa in valid_aa for aa in epitope):
                valid_epitopes.append(epitope)
                logging.info(f"Valid epitope: {epitope} ({len(epitope)} aa)")
            else:
                logging.warning(f"Invalid epitope (contains non-amino acid characters): {epitope}")

        return valid_epitopes

    def create_fasta_file(self, epitopes):
        """Create temporary FASTA file with epitopes"""
        temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)

        for i, epitope in enumerate(epitopes, 1):
            temp_fasta.write(f">epitope_{i}|{epitope}\n{epitope}\n")

        temp_fasta.close()
        return temp_fasta.name

    def run_blast(self, fasta_file, database, evalue=10.0, word_size=3, max_targets=100):
        """Run BLAST search against specified database"""

        # Check if database exists
        db_path = self.database_paths.get(database)
        if not db_path:
            logging.error(f"Unknown database: {database}")
            return None

        # Create output file
        output_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        output_file.close()

        # BLAST command
        blast_cmd = [
            'blastp',
            '-query', fasta_file,
            '-db', db_path,
            '-evalue', str(evalue),
            '-word_size', str(word_size),
            '-max_target_seqs', str(max_targets),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
            '-out', output_file.name
        ]

        try:
            logging.info(f"Running BLAST against {database} database...")
            result = subprocess.run(blast_cmd, capture_output=True, text=True, timeout=300)

            if result.returncode != 0:
                logging.error(f"BLAST failed: {result.stderr}")
                # Check if it's a database not found error
                if "No alias or index file found" in result.stderr:
                    logging.warning(f"Database '{database}' not found. This script requires pre-installed BLAST databases.")
                    logging.info("To set up BLAST databases, run: update_blastdb.pl --decompress nr")
                return None

            return output_file.name

        except subprocess.TimeoutExpired:
            logging.error("BLAST search timed out")
            return None
        except FileNotFoundError:
            logging.error("BLAST not found. Please ensure BLAST+ is installed and in PATH")
            return None

    def parse_blast_results(self, blast_output_file, epitopes):
        """Parse BLAST output and create results dataframe"""

        if not os.path.exists(blast_output_file) or os.path.getsize(blast_output_file) == 0:
            logging.warning("No BLAST results found")
            return pd.DataFrame()

        # Column names for BLAST output format 6
        columns = [
            'query_id', 'subject_id', 'identity_pct', 'alignment_length',
            'mismatches', 'gap_opens', 'query_start', 'query_end',
            'subject_start', 'subject_end', 'evalue', 'bitscore', 'subject_title'
        ]

        try:
            df = pd.read_csv(blast_output_file, sep='\t', names=columns)

            # Add epitope sequence
            epitope_map = {f"epitope_{i+1}": epitope for i, epitope in enumerate(epitopes)}
            df['epitope_sequence'] = df['query_id'].map(epitope_map)

            # Extract organism from subject title
            df['organism'] = df['subject_title'].str.extract(r'\[([^\]]+)\]')

            # Calculate match quality score
            df['match_score'] = (df['identity_pct'] * df['alignment_length']) / 100

            return df

        except Exception as e:
            logging.error(f"Error parsing BLAST results: {e}")
            return pd.DataFrame()

    def analyze_cross_reactivity(self, results_df, identity_threshold=70, length_threshold=0.8):
        """Analyze cross-reactivity potential"""

        if results_df.empty:
            return {
                'high_risk': pd.DataFrame(),
                'medium_risk': pd.DataFrame(),
                'low_risk': pd.DataFrame(),
                'summary': {
                    'total_hits': 0,
                    'high_risk_hits': 0,
                    'medium_risk_hits': 0,
                    'low_risk_hits': 0,
                    'unique_organisms': 0,
                    'top_organisms': {}
                }
            }

        # Calculate epitope coverage
        results_df['epitope_length'] = results_df['epitope_sequence'].str.len()
        results_df['coverage'] = results_df['alignment_length'] / results_df['epitope_length']

        # Risk classification
        high_risk = results_df[
            (results_df['identity_pct'] >= identity_threshold) &
            (results_df['coverage'] >= length_threshold)
        ]

        medium_risk = results_df[
            ((results_df['identity_pct'] >= identity_threshold - 20) &
             (results_df['coverage'] >= length_threshold)) |
            ((results_df['identity_pct'] >= identity_threshold) &
             (results_df['coverage'] >= length_threshold - 0.3))
        ]
        medium_risk = medium_risk[~medium_risk.index.isin(high_risk.index)]

        low_risk = results_df[
            ~results_df.index.isin(high_risk.index) &
            ~results_df.index.isin(medium_risk.index)
        ]

        # Summary statistics
        summary = {
            'total_hits': len(results_df),
            'high_risk_hits': len(high_risk),
            'medium_risk_hits': len(medium_risk),
            'low_risk_hits': len(low_risk),
            'unique_organisms': results_df['organism'].nunique(),
            'top_organisms': results_df['organism'].value_counts().head(10).to_dict()
        }

        return {
            'high_risk': high_risk,
            'medium_risk': medium_risk,
            'low_risk': low_risk,
            'summary': summary
        }

    def generate_summary_table(self, epitopes, analysis_results, output_dir='.'):
        """Generate comprehensive summary table of all epitope matches across databases"""

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        table_file = Path(output_dir) / f"epitope_summary_table_{timestamp}.tsv"

        # Collect all matches across databases
        all_matches = []

        for db_name, db_analysis in analysis_results.items():
            if db_analysis is None:
                # Add row showing no matches for this database
                for epitope in epitopes:
                    all_matches.append({
                        'Epitope': epitope,
                        'Database': db_name,
                        'Matches_Found': 0,
                        'Best_Identity_%': 'N/A',
                        'Best_E_value': 'N/A',
                        'Best_Match_ID': 'N/A',
                        'Best_Match_Organism': 'N/A',
                        'High_Risk_Matches': 0,
                        'Medium_Risk_Matches': 0,
                        'Total_Matches': 0
                    })
                continue

            # Get results for this database
            summary = db_analysis['summary']
            high_risk_df = db_analysis['high_risk']
            medium_risk_df = db_analysis['medium_risk']

            # Group by epitope to get per-epitope statistics
            if not high_risk_df.empty or not medium_risk_df.empty:
                # Combine all results for analysis
                all_results = []
                if not high_risk_df.empty:
                    all_results.append(high_risk_df)
                if not medium_risk_df.empty:
                    all_results.append(medium_risk_df)
                if not db_analysis['low_risk'].empty:
                    all_results.append(db_analysis['low_risk'])

                if all_results:
                    combined_df = pd.concat(all_results, ignore_index=True)

                    # Group by epitope
                    for epitope in epitopes:
                        epitope_matches = combined_df[combined_df['epitope_sequence'] == epitope]

                        if epitope_matches.empty:
                            all_matches.append({
                                'Epitope': epitope,
                                'Database': db_name,
                                'Matches_Found': 0,
                                'Best_Identity_%': 'N/A',
                                'Best_E_value': 'N/A',
                                'Best_Match_ID': 'N/A',
                                'Best_Match_Organism': 'N/A',
                                'High_Risk_Matches': 0,
                                'Medium_Risk_Matches': 0,
                                'Total_Matches': 0
                            })
                        else:
                            # Find best match (highest identity)
                            best_match = epitope_matches.loc[epitope_matches['identity_pct'].idxmax()]

                            # Count risk levels for this epitope
                            epitope_high = len(high_risk_df[high_risk_df['epitope_sequence'] == epitope])
                            epitope_medium = len(medium_risk_df[medium_risk_df['epitope_sequence'] == epitope])

                            all_matches.append({
                                'Epitope': epitope,
                                'Database': db_name,
                                'Matches_Found': 1,
                                'Best_Identity_%': f"{best_match['identity_pct']:.1f}",
                                'Best_E_value': f"{best_match['evalue']:.2e}",
                                'Best_Match_ID': best_match['subject_id'],
                                'Best_Match_Organism': best_match.get('organism', 'Unknown'),
                                'High_Risk_Matches': epitope_high,
                                'Medium_Risk_Matches': epitope_medium,
                                'Total_Matches': len(epitope_matches)
                            })
                else:
                    # No matches found
                    for epitope in epitopes:
                        all_matches.append({
                            'Epitope': epitope,
                            'Database': db_name,
                            'Matches_Found': 0,
                            'Best_Identity_%': 'N/A',
                            'Best_E_value': 'N/A',
                            'Best_Match_ID': 'N/A',
                            'Best_Match_Organism': 'N/A',
                            'High_Risk_Matches': 0,
                            'Medium_Risk_Matches': 0,
                            'Total_Matches': 0
                        })
            else:
                # No matches found
                for epitope in epitopes:
                    all_matches.append({
                        'Epitope': epitope,
                        'Database': db_name,
                        'Matches_Found': 0,
                        'Best_Identity_%': 'N/A',
                        'Best_E_value': 'N/A',
                        'Best_Match_ID': 'N/A',
                        'Best_Match_Organism': 'N/A',
                        'High_Risk_Matches': 0,
                        'Medium_Risk_Matches': 0,
                        'Total_Matches': 0
                    })

        # Create DataFrame and save
        summary_df = pd.DataFrame(all_matches)
        summary_df.to_csv(table_file, sep='\t', index=False)

        logging.info(f"Summary table saved to: {table_file}")
        return table_file, summary_df

    def generate_report(self, epitopes, databases, analysis_results, output_dir='.'):
        """Generate comprehensive cross-reactivity report"""

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = Path(output_dir) / f"epitope_cross_reactivity_report_{timestamp}.txt"

        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("EPITOPE CROSS-REACTIVITY ANALYSIS REPORT\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Epitopes analyzed: {len(epitopes)}\n")
            f.write(f"Databases searched: {', '.join(databases)}\n\n")

            # Epitope summary
            f.write("EPITOPES ANALYZED:\n")
            f.write("-" * 20 + "\n")
            for i, epitope in enumerate(epitopes, 1):
                f.write(f"{i}. {epitope} ({len(epitope)} amino acids)\n")
            f.write("\n")

            # Cross-reactivity summary for each database
            for db_name, db_analysis in analysis_results.items():
                f.write(f"DATABASE: {db_name.upper()}\n")
                f.write("-" * 40 + "\n")

                if db_analysis is None:
                    f.write("No analysis results (BLAST search failed)\n\n")
                    continue

                summary = db_analysis['summary']
                f.write(f"Total hits: {summary['total_hits']}\n")
                f.write(f"High risk hits: {summary['high_risk_hits']}\n")
                f.write(f"Medium risk hits: {summary['medium_risk_hits']}\n")
                f.write(f"Low risk hits: {summary['low_risk_hits']}\n")
                f.write(f"Unique organisms: {summary['unique_organisms']}\n\n")

                # Top organisms
                if summary['top_organisms']:
                    f.write("Top organisms with hits:\n")
                    for org, count in list(summary['top_organisms'].items())[:5]:
                        f.write(f"  {org}: {count} hits\n")
                f.write("\n")

                # High risk details
                if not db_analysis['high_risk'].empty:
                    f.write("HIGH RISK MATCHES:\n")
                    for _, hit in db_analysis['high_risk'].head(10).iterrows():
                        f.write(f"  Epitope: {hit['epitope_sequence']}\n")
                        f.write(f"  Match: {hit['subject_id']} | {hit['subject_title']}\n")
                        f.write(f"  Identity: {hit['identity_pct']:.1f}% | Coverage: {hit['coverage']:.1f}\n")
                        f.write(f"  Organism: {hit.get('organism', 'Unknown')}\n")
                        f.write("\n")

            f.write("="*80 + "\n")
            f.write("Analysis completed successfully.\n")

        logging.info(f"Report saved to: {report_file}")
        return report_file

def main():
    parser = argparse.ArgumentParser(
        description="Analyze epitope cross-reactivity using BLAST",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s "ALKYLPHENYLALANINE"
  %(prog)s "PEPTIDE1,PEPTIDE2,PEPTIDE3" --databases human,viral
  %(prog)s "EPITOPE" --evalue 1.0 --word-size 3 --identity-threshold 80
        """
    )

    parser.add_argument('epitopes',
                       help='Epitope sequences (comma-separated)')

    parser.add_argument('--databases', '-d',
                       default='all',
                       help='BLAST databases to search (comma-separated). Use "all" to search all available databases. Options: human,viral,mouse,cow,pig,chicken,ecoli,salmonella,staph,fungal,rat,horse,zebrafish,arabidopsis,buffalo,goat,rice,salmon,sheep,soybean (default: all)')

    parser.add_argument('--evalue', '-e',
                       type=float, default=10.0,
                       help='BLAST E-value threshold (default: 10.0)')

    parser.add_argument('--word-size', '-w',
                       type=int, default=3,
                       help='BLAST word size (default: 3)')

    parser.add_argument('--max-targets', '-m',
                       type=int, default=100,
                       help='Maximum target sequences per epitope (default: 100)')

    parser.add_argument('--identity-threshold', '-i',
                       type=float, default=70.0,
                       help='Identity threshold for high-risk classification (default: 70.0)')

    parser.add_argument('--coverage-threshold', '-c',
                       type=float, default=0.8,
                       help='Coverage threshold for high-risk classification (default: 0.8)')

    parser.add_argument('--output-dir', '-o',
                       default='.',
                       help='Output directory for reports (default: current directory)')

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = EpitopeCrossReactivityAnalyzer()

    # Parse epitopes
    epitopes = analyzer.parse_epitopes(args.epitopes)
    if not epitopes:
        logging.error("No valid epitopes provided")
        return

    # Parse databases
    if args.databases.lower() == 'all':
        databases = list(analyzer.available_databases.keys())
        logging.info(f"Searching all available databases: {', '.join(databases)}")
    else:
        databases = [db.strip() for db in args.databases.split(',')]
        invalid_dbs = [db for db in databases if db not in analyzer.available_databases]
        if invalid_dbs:
            logging.error(f"Invalid databases: {', '.join(invalid_dbs)}")
            logging.info(f"Available databases: {', '.join(analyzer.available_databases.keys())}")
            return

    # Create FASTA file
    fasta_file = analyzer.create_fasta_file(epitopes)

    try:
        # Run BLAST against each database
        analysis_results = {}

        for database in databases:
            logging.info(f"Analyzing against {database} database...")

            # Run BLAST
            blast_output = analyzer.run_blast(
                fasta_file, database,
                evalue=args.evalue,
                word_size=args.word_size,
                max_targets=args.max_targets
            )

            if blast_output:
                # Parse results
                results_df = analyzer.parse_blast_results(blast_output, epitopes)

                # Analyze cross-reactivity
                analysis = analyzer.analyze_cross_reactivity(
                    results_df,
                    identity_threshold=args.identity_threshold,
                    length_threshold=args.coverage_threshold
                )

                analysis_results[database] = analysis

                # Cleanup
                os.unlink(blast_output)
            else:
                logging.warning(f"BLAST search failed for database: {database}")
                analysis_results[database] = None

        # Generate comprehensive summary table
        if analysis_results:
            table_file, summary_df = analyzer.generate_summary_table(
                epitopes, analysis_results, args.output_dir
            )

            # Generate detailed report if there are actual results
            report_file = None
            if any(analysis_results.values()):
                report_file = analyzer.generate_report(
                    epitopes, databases, analysis_results, args.output_dir
                )

            # Print summary table to console
            print("\n" + "="*80)
            print("EPITOPE CROSS-REACTIVITY ANALYSIS SUMMARY TABLE")
            print("="*80)
            print(summary_df.to_string(index=False))

            print(f"\n\nFiles generated:")
            print(f"Summary table: {table_file}")
            if report_file:
                print(f"Detailed report: {report_file}")

            # Print quick statistics
            print(f"\nQuick Statistics:")
            print(f"Epitopes analyzed: {len(epitopes)}")
            print(f"Databases searched: {len(databases)}")

            matches_found = summary_df[summary_df['Matches_Found'] == 1]
            if not matches_found.empty:
                print(f"Epitope-Database combinations with matches: {len(matches_found)}")
                high_risk_total = summary_df['High_Risk_Matches'].sum()
                medium_risk_total = summary_df['Medium_Risk_Matches'].sum()
                print(f"High-risk matches across all databases: {high_risk_total}")
                print(f"Medium-risk matches across all databases: {medium_risk_total}")
            else:
                print("No matches found across any database")
        else:
            logging.error("No analysis results generated")

    finally:
        # Cleanup
        if os.path.exists(fasta_file):
            os.unlink(fasta_file)

if __name__ == "__main__":
    main()