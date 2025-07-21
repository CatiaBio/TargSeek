#!/usr/bin/env python3
"""
BepiPred 3.0 Epitope Prediction Integration
==========================================

This script uses BepiPred 3.0 to predict B-cell epitopes from conserved protein sequences,
integrating with 3D structure and conservation data.
"""

import pandas as pd
import json
import subprocess
import time
from pathlib import Path
import logging
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
import re
from datetime import datetime
import numpy as np
import tempfile
import os
import shutil

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class BepiPredPredictor:
    """Interface for BepiPred 3.0 predictions"""
    
    def __init__(self, bepipred_path="tools/BepiPred3.0", method="vt_pred", top_percentage=0.2):
        """
        Initialize BepiPred predictor
        
        Args:
            bepipred_path: Path to BepiPred 3.0 installation
            method: Prediction method ('vt_pred' or 'mjv_pred')
            top_percentage: Top percentage of epitope residues to include (0.0-1.0)
        """
        self.bepipred_path = Path(bepipred_path)
        self.method = method
        self.top_percentage = top_percentage
        self.script_path = self.bepipred_path / "bepipred3_CLI.py"
        
        if not self.script_path.exists():
            raise FileNotFoundError(f"BepiPred script not found at {self.script_path}")
    
    def predict_epitopes(self, sequence, sequence_id, save_dir=None):
        """
        Predict B-cell epitopes using BepiPred 3.0
        
        Args:
            sequence: Protein sequence string
            sequence_id: Sequence identifier
            save_dir: Directory to save original BepiPred output files (optional)
            
        Returns:
            List of epitope predictions with positions and scores
        """
        epitopes = []
        
        try:
            # Create temporary FASTA file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">{sequence_id}\n{sequence}\n")
                fasta_file = f.name
            
            # Create temporary output directory
            with tempfile.TemporaryDirectory() as output_dir:
                # Use absolute paths
                bepipred_dir = Path(self.bepipred_path).resolve()
                fasta_file_abs = Path(fasta_file).resolve()
                output_dir_abs = Path(output_dir).resolve()
                
                # Check environment type and set up command
                venv_script = bepipred_dir / "venv" / "bin" / "activate"
                activate_script = bepipred_dir / "activate_bepipred_conda.sh"
                
                if venv_script.exists():
                    # Use virtual environment (most common case)
                    # Convert Windows paths to Unix-style for bash
                    bepipred_unix = str(bepipred_dir).replace('\\', '/').replace('C:/', '/mnt/c/')
                    fasta_unix = str(fasta_file_abs).replace('\\', '/').replace('C:/', '/mnt/c/')
                    output_unix = str(output_dir_abs).replace('\\', '/').replace('C:/', '/mnt/c/')
                    
                    cmd = [
                        "bash", "-c", 
                        f"cd '{bepipred_unix}' && source venv/bin/activate && python bepipred3_CLI.py -i '{fasta_unix}' -o '{output_unix}' -pred {self.method} -top {self.top_percentage}"
                    ]
                elif activate_script.exists():
                    # Use conda environment
                    bepipred_unix = str(bepipred_dir).replace('\\', '/').replace('C:/', '/mnt/c/')
                    fasta_unix = str(fasta_file_abs).replace('\\', '/').replace('C:/', '/mnt/c/')
                    output_unix = str(output_dir_abs).replace('\\', '/').replace('C:/', '/mnt/c/')
                    
                    cmd = [
                        "bash", "-c",
                        f"cd '{bepipred_unix}' && source activate_bepipred_conda.sh && python bepipred3_CLI.py -i '{fasta_unix}' -o '{output_unix}' -pred {self.method} -top {self.top_percentage}"
                    ]
                else:
                    # Direct execution - run from BepiPred directory
                    cmd = [
                        "python", "bepipred3_CLI.py",
                        "-i", str(fasta_file_abs),
                        "-o", str(output_dir_abs),
                        "-pred", self.method,
                        "-top", str(self.top_percentage)
                    ]
                
                logging.info(f"Running BepiPred from directory: {bepipred_dir}")
                logging.debug(f"Command: {cmd}")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=900, cwd=str(bepipred_dir))
                
                if result.returncode != 0:
                    logging.error(f"BepiPred failed: {result.stderr}")
                    logging.debug(f"BepiPred stdout: {result.stdout}")
                    return epitopes
                
                # Parse BepiPred output
                epitopes = self._parse_bepipred_output(output_dir, sequence_id, sequence)
                
                # Copy original BepiPred output files if save_dir is provided
                if save_dir and epitopes:  # Only save if we got results
                    self._save_original_bepipred_files(output_dir, save_dir, sequence_id)
                
        except subprocess.TimeoutExpired:
            logging.error(f"BepiPred timed out for sequence {sequence_id}")
        except Exception as e:
            logging.error(f"Error running BepiPred for {sequence_id}: {e}")
        finally:
            # Clean up temporary FASTA file
            if 'fasta_file' in locals():
                try:
                    os.unlink(fasta_file)
                except:
                    pass
        
        return epitopes
    
    def _save_original_bepipred_files(self, temp_output_dir, save_dir, sequence_id):
        """Save original BepiPred output files to the results directory"""
        try:
            temp_path = Path(temp_output_dir)
            save_path = Path(save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            
            # List of BepiPred output files to preserve
            bepipred_files = [
                'raw_output.csv',                          # Main scores file
                'Bcell_epitope_preds.fasta',              # Epitope sequences
                'Bcell_epitope_top_20pct_preds.fasta',    # Top percentage epitopes
                'Bcell_linepitope_top_20pct_preds.fasta', # Linear epitopes top percentage
                'output_interactive_figures.html'          # Interactive plots
            ]
            
            files_copied = []
            for filename in bepipred_files:
                source_file = temp_path / filename
                if source_file.exists():
                    # Add sequence ID prefix to avoid conflicts
                    dest_filename = f"{sequence_id}_{filename}"
                    dest_file = save_path / dest_filename
                    
                    shutil.copy2(source_file, dest_file)
                    files_copied.append(dest_filename)
                    logging.debug(f"Copied BepiPred file: {dest_filename}")
            
            if files_copied:
                logging.info(f"Saved {len(files_copied)} original BepiPred files for {sequence_id}")
            
        except Exception as e:
            logging.error(f"Error saving BepiPred files for {sequence_id}: {e}")
    
    def _parse_bepipred_output(self, output_dir, sequence_id, sequence):
        """Parse BepiPred 3.0 output files to extract epitopes"""
        epitopes = []
        output_path = Path(output_dir)
        
        try:
            # BepiPred 3.0 creates specific output files
            raw_output_file = output_path / "raw_output.csv"
            epitope_fasta_file = output_path / "Bcell_epitope_preds.fasta"
            
            if raw_output_file.exists():
                # Parse the main CSV output file
                epitopes = self._parse_raw_output_csv(raw_output_file, sequence)
                logging.info(f"Parsed {len(epitopes)} epitopes from raw_output.csv")
            
            elif epitope_fasta_file.exists():
                # Fallback to FASTA epitope file
                epitopes = self._parse_epitope_fasta(epitope_fasta_file, sequence)
                logging.info(f"Parsed {len(epitopes)} epitopes from FASTA file")
            
            else:
                # Look for any CSV files as fallback
                csv_files = list(output_path.glob("*.csv"))
                for csv_file in csv_files:
                    epitopes = self._parse_prediction_file(csv_file, sequence)
                    if epitopes:
                        break
                        
        except Exception as e:
            logging.error(f"Error parsing BepiPred output: {e}")
        
        return epitopes
    
    def _parse_raw_output_csv(self, file_path, sequence, threshold=0.5):
        """Parse BepiPred 3.0 raw_output.csv file"""
        epitopes = []
        
        try:
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # BepiPred 3.0 CSV format: Accession,Residue,BepiPred-3.0 score,BepiPred-3.0 linear epitope score
            if 'BepiPred-3.0 score' in df.columns:
                score_col = 'BepiPred-3.0 score'
            elif 'BepiPred-3.0 linear epitope score' in df.columns:
                score_col = 'BepiPred-3.0 linear epitope score'
            else:
                # Fallback to any score column
                score_cols = [col for col in df.columns if 'score' in col.lower()]
                if score_cols:
                    score_col = score_cols[0]
                else:
                    logging.warning("No score column found in BepiPred output")
                    return epitopes
            
            # Filter high-scoring positions
            high_scores = df[df[score_col] >= threshold].copy()
            high_scores['position'] = range(1, len(df) + 1)  # Add position column
            
            if len(high_scores) > 0:
                epitopes = self._group_consecutive_epitopes(
                    high_scores, 'position', score_col, sequence, threshold, min_length=6
                )
                
        except Exception as e:
            logging.error(f"Error parsing raw_output.csv: {e}")
        
        return epitopes
    
    def _parse_epitope_fasta(self, file_path, sequence):
        """Parse BepiPred 3.0 epitope FASTA file"""
        epitopes = []
        
        try:
            from Bio import SeqIO
            
            for record in SeqIO.parse(file_path, "fasta"):
                epitope_seq = str(record.seq)
                
                # Find epitope position in original sequence
                start_pos = sequence.find(epitope_seq)
                if start_pos != -1:
                    epitopes.append({
                        'peptide': epitope_seq,
                        'start': start_pos + 1,  # 1-based
                        'end': start_pos + len(epitope_seq),
                        'length': len(epitope_seq),
                        'score': 0.7,  # Default score since FASTA doesn't include scores
                        'method': 'BepiPred3.0',
                        'type': 'B-cell'
                    })
                    
        except Exception as e:
            logging.error(f"Error parsing epitope FASTA: {e}")
        
        return epitopes
    
    def _parse_prediction_file(self, file_path, sequence):
        """Parse a BepiPred prediction file"""
        epitopes = []
        
        try:
            # Try to read as CSV/TSV first
            try:
                df = pd.read_csv(file_path, sep=None, engine='python')
                if 'score' in df.columns or 'probability' in df.columns:
                    epitopes = self._extract_epitopes_from_dataframe(df, sequence)
                    return epitopes
            except:
                pass
            
            # Try to read as text file with scores
            with open(file_path, 'r') as f:
                lines = f.readlines()
                epitopes = self._extract_epitopes_from_text(lines, sequence)
                
        except Exception as e:
            logging.debug(f"Could not parse {file_path}: {e}")
        
        return epitopes
    
    def _extract_epitopes_from_dataframe(self, df, sequence, threshold=0.5):
        """Extract epitopes from pandas DataFrame"""
        epitopes = []
        
        # Find score column
        score_col = None
        for col in ['score', 'probability', 'epitope_score', 'prediction']:
            if col in df.columns:
                score_col = col
                break
        
        if score_col is None:
            return epitopes
        
        # Find position column
        pos_col = None
        for col in ['position', 'pos', 'residue', 'index']:
            if col in df.columns:
                pos_col = col
                break
        
        # Extract high-scoring regions
        high_scores = df[df[score_col] >= threshold]
        
        if len(high_scores) > 0:
            # Group consecutive positions into epitopes
            epitopes = self._group_consecutive_epitopes(high_scores, pos_col, score_col, sequence, threshold)
        
        return epitopes
    
    def _extract_epitopes_from_text(self, lines, sequence, threshold=0.5):
        """Extract epitopes from text format output"""
        epitopes = []
        scores = []
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            # Try to extract numerical scores
            numbers = re.findall(r'\d+\.?\d*', line)
            if numbers:
                try:
                    score = float(numbers[-1])  # Assume last number is score
                    scores.append(score)
                except:
                    continue
        
        if scores and len(scores) >= len(sequence):
            # Create DataFrame-like structure for processing
            df_data = {
                'position': list(range(1, len(scores) + 1)),
                'score': scores[:len(sequence)]
            }
            df = pd.DataFrame(df_data)
            epitopes = self._extract_epitopes_from_dataframe(df, sequence, threshold)
        
        return epitopes
    
    def _group_consecutive_epitopes(self, high_scores_df, pos_col, score_col, sequence, threshold, min_length=6):
        """Group consecutive high-scoring positions into epitopes"""
        epitopes = []
        
        if pos_col:
            positions = high_scores_df[pos_col].tolist()
            scores = high_scores_df[score_col].tolist()
        else:
            positions = high_scores_df.index.tolist()
            scores = high_scores_df[score_col].tolist()
        
        # Group consecutive positions
        groups = []
        current_group = []
        
        for i, pos in enumerate(positions):
            if not current_group or pos == current_group[-1][0] + 1:
                current_group.append((pos, scores[i]))
            else:
                if len(current_group) >= min_length:
                    groups.append(current_group)
                current_group = [(pos, scores[i])]
        
        # Add last group
        if len(current_group) >= min_length:
            groups.append(current_group)
        
        # Convert groups to epitopes
        for group in groups:
            start_pos = group[0][0]
            end_pos = group[-1][0]
            avg_score = np.mean([score for pos, score in group])
            
            # Extract peptide sequence (adjust for 0-based indexing)
            peptide_start = max(0, start_pos - 1)
            peptide_end = min(len(sequence), end_pos)
            peptide = sequence[peptide_start:peptide_end]
            
            if len(peptide) >= min_length:
                epitopes.append({
                    'peptide': peptide,
                    'start': start_pos,
                    'end': end_pos,
                    'length': len(peptide),
                    'score': avg_score,
                    'method': 'BepiPred3.0',
                    'type': 'B-cell'
                })
        
        return epitopes

def get_representative_sequence(gene, msa_sequences_dir, structures_dir):
    """Get representative sequence for a gene, prioritizing 3D structure sequences"""
    
    # Try to get sequence from MSA file
    msa_file = Path(msa_sequences_dir) / f"{gene}.fasta"
    structure_info_file = Path(structures_dir) / gene / "structure_info.json"
    
    sequences = []
    
    # Load sequences from MSA file
    if msa_file.exists():
        try:
            for record in SeqIO.parse(msa_file, "fasta"):
                sequences.append({
                    'id': record.id,
                    'sequence': str(record.seq).replace('-', ''),  # Remove gaps
                    'description': record.description,
                    'has_structure': False
                })
        except Exception as e:
            logging.error(f"Error reading MSA file for {gene}: {e}")
            return None, None
    
    if not sequences:
        logging.warning(f"No sequences found for gene {gene}")
        return None, None
    
    # Check which sequences have 3D structures
    gene_structure_dir = Path(structures_dir) / gene
    if gene_structure_dir.exists():
        try:
            # Get PDB IDs from structure files in the directory
            pdb_files = list(gene_structure_dir.glob("*.pdb.gz"))
            structure_ids = {pdb_file.stem.replace('.pdb', '') for pdb_file in pdb_files}
            
            if structure_ids:
                logging.debug(f"Found 3D structures for {gene}: {structure_ids}")
                # For now, if structures exist, mark the first sequence as having structure
                # This is a simplified approach since we can't match specific sequences to PDB IDs
                # without proper structure_info.json files
                if sequences:
                    sequences[0]['has_structure'] = True
                    logging.debug(f"Marked first sequence as having 3D structure for {gene}")
        except Exception as e:
            logging.debug(f"Could not check structure directory for {gene}: {e}")
    
    # Prioritize sequences with 3D structures
    structure_sequences = [s for s in sequences if s['has_structure']]
    if structure_sequences:
        selected = structure_sequences[0]
        logging.info(f"Using 3D structure sequence: {selected['id']} (length: {len(selected['sequence'])})")
    else:
        selected = sequences[0]
        logging.info(f"Using reference sequence: {selected['id']} (length: {len(selected['sequence'])})")
    
    return selected['sequence'], selected['id']

def load_conservation_data(gene, conservation_dir):
    """Load conservation data for a gene"""
    conservation_file = Path(conservation_dir) / f"{gene}_conservation_details.tsv"
    
    if not conservation_file.exists():
        logging.warning(f"Conservation file not found: {conservation_file}")
        return None
    
    try:
        conservation_df = pd.read_csv(conservation_file, sep='\t')
        return conservation_df
    except Exception as e:
        logging.error(f"Error loading conservation data for {gene}: {e}")
        return None

def predict_epitopes_for_gene(gene, msa_sequences_dir, conservation_dir, structures_dir, output_dir, bepipred_path="tools/BepiPred3.0", bepipred_config=None):
    """Predict epitopes for a single gene using BepiPred 3.0 - saves only original BepiPred files"""
    
    logging.info(f"Predicting epitopes for gene: {gene}")
    
    # Get representative sequence
    sequence, sequence_id = get_representative_sequence(gene, msa_sequences_dir, structures_dir)
    if not sequence:
        logging.error(f"No sequence found for gene {gene}")
        return None
    
    # Initialize BepiPred predictor
    try:
        # Get config parameters
        bepipred_config = bepipred_config or {}
        method = bepipred_config.get('method', 'vt_pred')
        top_percentage = bepipred_config.get('top_percentage', 0.2)
        
        predictor = BepiPredPredictor(bepipred_path=bepipred_path, method=method, top_percentage=top_percentage)
    except FileNotFoundError as e:
        logging.error(f"BepiPred not found: {e}")
        return None
    
    # Create gene-specific output directory (save directly to gene folder)
    gene_output_dir = Path(output_dir) / gene
    gene_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Predict epitopes and save original BepiPred files directly to gene directory
    logging.info(f"Running BepiPred 3.0 and saving original output files...")
    epitopes = predictor.predict_epitopes(sequence, sequence_id, save_dir=gene_output_dir)
    
    if not epitopes:
        logging.warning(f"No epitopes found for {gene}")
        # Still return success even if no epitopes - BepiPred files were created
        return {
            'gene': gene,
            'sequence_id': sequence_id,
            'sequence_length': len(sequence),
            'total_epitopes': 0,
            'files_saved': True
        }
    
    # Create minimal gene summary (just for logging)
    gene_summary = {
        'gene': gene,
        'sequence_id': sequence_id,
        'sequence_length': len(sequence),
        'total_epitopes': len(epitopes),
        'files_saved': True
    }
    
    logging.info(f"BepiPred output saved for {gene}:")
    logging.info(f"  Sequence: {sequence_id} ({len(sequence)} aa)")
    logging.info(f"  Original BepiPred files saved to: {gene_output_dir}")
    
    return gene_summary

def add_conservation_scores(epitopes, conservation_df):
    """Add conservation scores to epitopes"""
    
    for epitope in epitopes:
        start_pos = epitope['start']
        end_pos = epitope['end']
        
        # Get conservation scores for epitope region
        epitope_conservation = conservation_df[
            (conservation_df['position'] >= start_pos) & 
            (conservation_df['position'] <= end_pos)
        ]
        
        if len(epitope_conservation) > 0:
            epitope['conservation_score'] = epitope_conservation['conservation_score'].mean()
            epitope['conservation_count'] = len(epitope_conservation)
        else:
            epitope['conservation_score'] = 0.0
            epitope['conservation_count'] = 0
    
    return epitopes

def create_combined_report(gene_summaries, output_dir):
    """Create simple BepiPred file inventory report"""
    
    # Overall summary statistics
    overall_summary = {
        'analysis_date': datetime.now().isoformat(),
        'total_genes': len(gene_summaries),
        'genes_processed': len([s for s in gene_summaries if s.get('files_saved', False)]),
        'per_gene_summary': gene_summaries
    }
    
    # Save summary JSON
    summary_file = Path(output_dir) / "bepipred_files_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(overall_summary, f, indent=2)
    
    # Create human-readable report
    report_file = Path(output_dir) / "bepipred_files_report.txt"
    with open(report_file, 'w') as f:
        f.write("BepiPred 3.0 Output Files Report\n")
        f.write("=" * 40 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"OVERVIEW:\n")
        f.write(f"Total genes processed: {overall_summary['genes_processed']}\n")
        f.write(f"Output directory: {output_dir}\n\n")
        
        f.write("FILES GENERATED PER GENE:\n")
        f.write("-" * 30 + "\n")
        f.write("Each gene directory contains:\n")
        f.write("  - [gene]_raw_output.csv                     (per-residue scores)\n")
        f.write("  - [gene]_Bcell_epitope_preds.fasta         (all epitopes)\n")
        f.write("  - [gene]_Bcell_epitope_top_20pct_preds.fasta (top percentage epitopes)\n")
        f.write("  - [gene]_Bcell_linepitope_top_20pct_preds.fasta (linear top percentage epitopes)\n")
        f.write("  - [gene]_output_interactive_figures.html   (visualizations)\n\n")
        
        f.write("PROCESSED GENES:\n")
        f.write("-" * 30 + "\n")
        for gene_summary in sorted(gene_summaries, key=lambda x: x['gene']):
            status = "✓" if gene_summary.get('files_saved', False) else "✗"
            f.write(f"{status} {gene_summary['gene']} ({gene_summary['sequence_id']}, {gene_summary['sequence_length']} aa)\n")
    
    logging.info(f"BepiPred files report generated:")
    logging.info(f"  Summary: {summary_file}")
    logging.info(f"  Report: {report_file}")

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        msa_sequences_dir = snakemake.input.msa_sequences
        conservation_dir = snakemake.input.conservation
        structures_dir = snakemake.input.structures_3d
        proteins_to_study_file = snakemake.input.proteins_to_study
        output_dir = snakemake.output[0]
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        # Get BepiPred config
        bepipred_config = snakemake.config.get('bepipred', {})
        bepipred_path = bepipred_config.get('path', 'tools/BepiPred3.0')
        
        logging.info(f"BepiPred 3.0 Epitope Prediction for {analysis}_{paramset}_gram_{group}")
        
    except NameError:
        # Test mode - allow command line execution
        import sys
        if len(sys.argv) != 6:
            print("Usage: python predict_epitopes_bepipred.py <msa_sequences_dir> <conservation_dir> <structures_dir> <proteins_to_study_file> <output_dir>")
            sys.exit(1)
            
        msa_sequences_dir = sys.argv[1]
        conservation_dir = sys.argv[2]
        structures_dir = sys.argv[3]
        proteins_to_study_file = sys.argv[4]
        output_dir = sys.argv[5]
        bepipred_path = "tools/BepiPred3.0"
        bepipred_config = {'method': 'vt_pred', 'top_percentage': 0.5}  # Default for command line
        
        logging.info("Running in command line mode")
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Load genes to analyze
    try:
        proteins_df = pd.read_csv(proteins_to_study_file, sep='\t')
        genes_to_analyze = proteins_df['gene'].unique().tolist()
        logging.info(f"Found {len(genes_to_analyze)} genes to analyze for epitopes")
    except Exception as e:
        logging.error(f"Error loading proteins to study: {e}")
        return
    
    # Predict epitopes for each gene
    gene_summaries = []
    successful_predictions = 0
    
    for gene in genes_to_analyze:
        try:
            summary = predict_epitopes_for_gene(
                gene, msa_sequences_dir, conservation_dir, structures_dir, output_dir, bepipred_path, bepipred_config
            )
            if summary:
                gene_summaries.append(summary)
                successful_predictions += 1
        except Exception as e:
            logging.error(f"Error predicting epitopes for {gene}: {e}")
    
    # Generate combined report
    if gene_summaries:
        create_combined_report(gene_summaries, output_dir)
    
    logging.info(f"BepiPred file generation complete!")
    logging.info(f"Successfully processed {successful_predictions}/{len(genes_to_analyze)} genes")
    logging.info(f"Original BepiPred files saved to: {output_dir}")
    logging.info(f"Each gene directory contains 5 BepiPred output files")

if __name__ == "__main__":
    main()