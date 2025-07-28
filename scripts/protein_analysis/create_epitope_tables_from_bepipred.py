#!/usr/bin/env python3
"""
Create ElliPro-style epitope tables from BepiPred raw CSV output files
=====================================================================

This script processes BepiPred 3.0 raw CSV output files and creates:
1. Linear epitopes table (similar to ElliPro format)
2. Raw scores table with all residue information

Usage:
    python create_epitope_tables_from_bepipred.py <bepipred_output_dir>
    
Example:
    python scripts/protein_analysis/create_epitope_tables_from_bepipred.py results/analysis1_params1/protein_analysis/epitope_predictions_bepipred
"""

import pandas as pd
import json
from pathlib import Path
import logging
from Bio import SeqIO
import argparse
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class BepiPredTableGenerator:
    """Generate epitope tables from BepiPred CSV output"""
    
    def __init__(self, min_epitope_length=6):
        """
        Initialize table generator
        
        Args:
            min_epitope_length: Minimum length for linear epitopes (default: 6)
        """
        self.min_epitope_length = min_epitope_length
    
    def parse_bepipred_csv(self, csv_file, sequence_id, sequence):
        """Parse BepiPred CSV output to extract epitope information"""
        try:
            # Read CSV file
            df = pd.read_csv(csv_file)
            
            # Handle different CSV formats from BepiPred
            if 'Pos' in df.columns and 'AA' in df.columns:
                # Expected format: Pos, AA, Score, Prediction
                pos_col, aa_col, score_col = 'Pos', 'AA', 'Score'
                pred_col = 'Prediction' if 'Prediction' in df.columns else None
            elif 'Accession' in df.columns and 'Residue' in df.columns:
                # Actual BepiPred format: Accession, Residue, BepiPred-3.0 score, BepiPred-3.0 linear epitope score
                # Extract position from Accession (format: protein_pos)
                pos_col, aa_col = 'Accession', 'Residue'
                score_col = 'BepiPred-3.0 score'
                pred_col = 'BepiPred-3.0 linear epitope score'
            else:
                logging.warning(f"Unexpected CSV format in {csv_file}")
                logging.warning(f"Found columns: {list(df.columns)}")
                return None
            
            # Get all residues with scores
            all_residues = []
            for i, row in df.iterrows():
                if pos_col == 'Accession':
                    # Position is row index + 1 (BepiPred uses sequential numbering)
                    position = i + 1
                else:
                    position = int(row[pos_col])
                
                # Determine if this is an epitope residue based on score
                score = float(row[score_col])
                if pred_col and pred_col in df.columns:
                    # Use linear epitope score if available (BepiPred 3.0 range seems to be 0-1 with lower values)
                    epitope_score = float(row[pred_col])
                    prediction = 'E' if epitope_score > 0.15 else 'N'  # Lower threshold for BepiPred 3.0
                else:
                    # Use main score with threshold (BepiPred 3.0 scores typically < 0.5)
                    prediction = 'E' if score > 0.25 else 'N'
                
                all_residues.append({
                    'position': position,
                    'amino_acid': str(row[aa_col]),
                    'score': score,
                    'prediction': prediction
                })
            
            # Find epitope regions (consecutive residues with positive predictions)
            linear_epitopes = self._find_linear_epitopes(all_residues, sequence)
            
            return {
                'all_residues': all_residues,
                'linear_epitopes': linear_epitopes
            }
            
        except Exception as e:
            logging.error(f"Error parsing BepiPred CSV {csv_file}: {e}")
            return None
    
    def _find_linear_epitopes(self, all_residues, sequence):
        """Find linear epitope regions from BepiPred predictions"""
        epitopes = []
        current_epitope = []
        
        for residue in all_residues:
            if residue['prediction'] == 'E':  # Epitope residue
                current_epitope.append(residue)
            else:
                # End of epitope region
                if len(current_epitope) >= self.min_epitope_length:
                    epitopes.append(self._create_epitope_from_residues(current_epitope, sequence))
                current_epitope = []
        
        # Check final epitope
        if len(current_epitope) >= self.min_epitope_length:
            epitopes.append(self._create_epitope_from_residues(current_epitope, sequence))
        
        return epitopes
    
    def _create_epitope_from_residues(self, residues, sequence):
        """Create epitope entry from list of residues"""
        start_pos = residues[0]['position']
        end_pos = residues[-1]['position']
        
        # Extract peptide sequence
        peptide = sequence[start_pos-1:end_pos]  # Convert to 0-based indexing
        
        # Calculate average score
        avg_score = sum(r['score'] for r in residues) / len(residues)
        
        return {
            'start': start_pos,
            'end': end_pos,
            'peptide': peptide,
            'length': len(residues),
            'avg_score': avg_score,
            'residues': residues
        }
    
    def write_linear_epitopes_table(self, epitopes, output_file, sequence_id):
        """Write linear epitopes table similar to ElliPro format"""
        try:
            with open(output_file, 'w') as f:
                f.write("# Predicted Linear Epitope(s):\n")
                f.write("No.\tChain\tStart\tEnd\tPeptide\tNumber_of_Residues\tScore\t3D_Structure\n")
                
                for i, epitope in enumerate(epitopes, 1):
                    f.write(f"{i}\tA\t{epitope['start']}\t{epitope['end']}\t{epitope['peptide']}\t{epitope['length']}\t{epitope['avg_score']:.3f}\t{sequence_id}\n")
            
            logging.debug(f"Wrote {len(epitopes)} linear epitopes to {output_file}")
            
        except Exception as e:
            logging.error(f"Error writing linear epitopes table: {e}")
    
    def write_raw_scores_table(self, all_residues, output_file, sequence_id):
        """Write raw scores table with all residue information"""
        try:
            with open(output_file, 'w') as f:
                f.write("# BepiPred 3.0 Raw Scores and Predictions\n")
                f.write("No.\tChain\tResidue_Number\tResidue_Name\tScore\tPrediction\t3D_Structure\n")
                
                for i, residue in enumerate(all_residues, 1):
                    f.write(f"{i}\tA\t{residue['position']}\t{residue['amino_acid']}\t{residue['score']:.3f}\t{residue['prediction']}\t{sequence_id}\n")
            
            logging.debug(f"Wrote {len(all_residues)} residue scores to {output_file}")
            
        except Exception as e:
            logging.error(f"Error writing raw scores table: {e}")

def find_bepipred_csv_files(output_dir):
    """Find all BepiPred raw CSV output files"""
    output_path = Path(output_dir)
    csv_files = []
    
    # Look for CSV files in gene subdirectories
    for gene_dir in output_path.iterdir():
        if gene_dir.is_dir():
            for csv_file in gene_dir.glob("*_raw_output.csv"):
                csv_files.append(csv_file)
    
    # Also look for CSV files directly in output directory
    for csv_file in output_path.glob("*_raw_output.csv"):
        csv_files.append(csv_file)
    
    return sorted(csv_files)

def find_corresponding_fasta(csv_file):
    """Find the corresponding FASTA file for a CSV file"""
    # Extract info from CSV filename: gene_pdb_raw_output.csv
    csv_name = csv_file.stem  # Remove .csv
    if csv_name.endswith('_raw_output'):
        sequence_id = csv_name[:-11]  # Remove '_raw_output'
        
        # Look for FASTA files in common locations
        possible_paths = [
            # In the same directory
            csv_file.parent / f"{sequence_id}.fasta",
            # In data/protein_structures/gene/
            Path("data") / "protein_structures" / csv_file.parent.name / f"{sequence_id.split('_')[-1]}.fasta",
            # Try variations
            csv_file.parent / f"{sequence_id.split('_')[-1]}.fasta"
        ]
        
        for fasta_path in possible_paths:
            if fasta_path.exists():
                return fasta_path
    
    return None

def extract_sequence_from_fasta(fasta_file, sequence_id):
    """Extract the first sequence from a FASTA file"""
    try:
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        if sequences:
            return str(sequences[0].seq)
        else:
            logging.warning(f"No sequences found in {fasta_file}")
            return None
    except Exception as e:
        logging.error(f"Error reading FASTA file {fasta_file}: {e}")
        return None

def process_bepipred_outputs(output_dir, min_epitope_length=6):
    """Process all BepiPred outputs in a directory"""
    
    logging.info(f"Processing BepiPred outputs in: {output_dir}")
    
    # Find all CSV files
    csv_files = find_bepipred_csv_files(output_dir)
    
    if not csv_files:
        logging.error(f"No BepiPred CSV files found in {output_dir}")
        return
    
    logging.info(f"Found {len(csv_files)} BepiPred CSV files")
    
    # Initialize table generator
    generator = BepiPredTableGenerator(min_epitope_length=min_epitope_length)
    
    processed = 0
    failed = 0
    
    for csv_file in csv_files:
        logging.info(f"Processing: {csv_file}")
        
        # Extract sequence ID from filename
        csv_name = csv_file.stem
        if csv_name.endswith('_raw_output'):
            sequence_id = csv_name[:-11]
        else:
            sequence_id = csv_name
        
        # Find corresponding FASTA file
        fasta_file = find_corresponding_fasta(csv_file)
        if not fasta_file:
            logging.warning(f"Could not find FASTA file for {csv_file}")
            # Create a dummy sequence based on CSV length
            try:
                df = pd.read_csv(csv_file)
                sequence = 'A' * len(df)  # Dummy sequence
                logging.info(f"Using dummy sequence of length {len(sequence)}")
            except:
                failed += 1
                continue
        else:
            sequence = extract_sequence_from_fasta(fasta_file, sequence_id)
            if not sequence:
                failed += 1
                continue
        
        # Parse CSV and create tables
        epitope_data = generator.parse_bepipred_csv(csv_file, sequence_id, sequence)
        
        if epitope_data:
            # Create output files in the same directory as CSV
            output_dir_path = csv_file.parent
            
            # Linear epitopes table
            linear_file = output_dir_path / f"{sequence_id}_linear_epitopes.tsv"
            generator.write_linear_epitopes_table(epitope_data['linear_epitopes'], linear_file, sequence_id)
            
            # Raw scores table  
            raw_file = output_dir_path / f"{sequence_id}_raw_scores.tsv"
            generator.write_raw_scores_table(epitope_data['all_residues'], raw_file, sequence_id)
            
            logging.info(f"✓ Created tables for {sequence_id}: {len(epitope_data['linear_epitopes'])} linear epitopes")
            processed += 1
        else:
            logging.error(f"✗ Failed to process {sequence_id}")
            failed += 1
    
    logging.info(f"\n=== Summary ===")
    logging.info(f"Processed: {processed} structures")
    logging.info(f"Failed: {failed} structures")
    logging.info(f"Success rate: {processed/(processed+failed)*100:.1f}%" if (processed+failed) > 0 else "No files processed")

def main():
    """Main function for both command line and Snakemake usage"""
    
    # Check if running from Snakemake
    if 'snakemake' in globals():
        # Running from Snakemake
        bepipred_sentinel = snakemake.input.bepipred_sentinel
        sentinel_file = snakemake.output.epitope_tables_sentinel
        min_length = snakemake.params.get('min_epitope_length', 6)
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        
        # Derive BepiPred directory from sentinel file path
        from pathlib import Path
        bepipred_dir = Path(bepipred_sentinel).parent
        
        logging.info(f"Running epitope table creation for {analysis}_{paramset}")
        logging.info(f"BepiPred directory: {bepipred_dir}")
        logging.info(f"Sentinel file: {sentinel_file}")
        
        # Process BepiPred outputs
        process_bepipred_outputs(bepipred_dir, min_length)
        
        # Create sentinel file
        with open(sentinel_file, 'w') as f:
            f.write(f"Epitope tables creation completed for {analysis}_{paramset}\n")
            f.write(f"Minimum epitope length: {min_length}\n")
            f.write(f"Processed directory: {bepipred_dir}\n")
        
        logging.info(f"Sentinel file created: {sentinel_file}")
        
    else:
        # Running from command line
        parser = argparse.ArgumentParser(description='Create epitope tables from BepiPred CSV output')
        parser.add_argument('output_dir', help='Directory containing BepiPred output files')
        parser.add_argument('--min-length', type=int, default=6, help='Minimum epitope length (default: 6)')
        
        args = parser.parse_args()
        
        if not Path(args.output_dir).exists():
            logging.error(f"Output directory does not exist: {args.output_dir}")
            sys.exit(1)
        
        process_bepipred_outputs(args.output_dir, args.min_length)

if __name__ == "__main__":
    main()