#!/usr/bin/env python3
"""
BepiPred 3.0 Epitope Prediction for Selected 3D Structures
=========================================================

This simplified script runs BepiPred 3.0 only on the selected 3D structure sequences
from the MSA creation step. It reads the paths from selected_3d_paths.txt and
processes each sequence individually.
"""

import pandas as pd
import json
import subprocess
import time
from pathlib import Path
import logging
from Bio import SeqIO
import tempfile
import os
import shutil
from datetime import datetime

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
            save_dir: Directory to save original BepiPred output files
            
        Returns:
            bool: True if prediction succeeded, False otherwise
        """
        
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
                    # Use virtual environment
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
                    # Direct execution
                    cmd = [
                        "python", "bepipred3_CLI.py",
                        "-i", str(fasta_file_abs),
                        "-o", str(output_dir_abs),
                        "-pred", self.method,
                        "-top", str(self.top_percentage)
                    ]
                
                logging.info(f"Running BepiPred from directory: {bepipred_dir}")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=900, cwd=str(bepipred_dir))
                
                if result.returncode != 0:
                    logging.error(f"BepiPred failed: {result.stderr}")
                    return False
                
                # Copy original BepiPred output files if save_dir is provided
                if save_dir:
                    self._save_bepipred_files(output_dir, save_dir, sequence_id)
                
                return True
                
        except subprocess.TimeoutExpired:
            logging.error(f"BepiPred timed out for sequence {sequence_id}")
            return False
        except Exception as e:
            logging.error(f"Error running BepiPred for {sequence_id}: {e}")
            return False
        finally:
            # Clean up temporary FASTA file
            if 'fasta_file' in locals():
                try:
                    os.unlink(fasta_file)
                except:
                    pass
    
    def _save_bepipred_files(self, temp_output_dir, save_dir, sequence_id):
        """Save original BepiPred output files to the results directory"""
        try:
            temp_path = Path(temp_output_dir)
            save_path = Path(save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            
            # List of BepiPred output files to preserve
            bepipred_files = [
                'raw_output.csv',
                'Bcell_epitope_preds.fasta',
                'Bcell_epitope_top_20pct_preds.fasta',
                'Bcell_linepitope_top_20pct_preds.fasta',
                'output_interactive_figures.html'
            ]
            
            for filename in bepipred_files:
                source_file = temp_path / filename
                if source_file.exists():
                    # Add sequence ID prefix to avoid conflicts
                    dest_filename = f"{sequence_id}_{filename}"
                    dest_file = save_path / dest_filename
                    
                    shutil.copy2(source_file, dest_file)
                    logging.debug(f"Copied BepiPred file: {dest_filename}")
            
            logging.info(f"Saved BepiPred output files for {sequence_id}")
            
        except Exception as e:
            logging.error(f"Error saving BepiPred files for {sequence_id}: {e}")

def predict_epitopes_for_3d_structure(fasta_path, predictor, output_dir):
    """
    Predict epitopes for a single 3D structure FASTA file
    
    Args:
        fasta_path: Path to FASTA file
        predictor: BepiPredPredictor instance
        output_dir: Output directory for results
        
    Returns:
        dict: Summary information about the prediction
    """
    fasta_path = Path(fasta_path)
    
    if not fasta_path.exists():
        logging.warning(f"FASTA file not found: {fasta_path}")
        return None
    
    # Parse filename to extract information
    filename = fasta_path.stem
    parts = filename.split('_')
    pdb_id = parts[0] if parts else filename
    
    # Extract gene name from path (parent directory)
    gene_name = fasta_path.parent.name
    
    # Read sequence
    try:
        sequences = list(SeqIO.parse(fasta_path, "fasta"))
        if not sequences:
            logging.warning(f"No sequences found in {fasta_path}")
            return None
        
        seq_record = sequences[0]  # Take first sequence
        sequence = str(seq_record.seq)
        
        # Create unique identifier
        sequence_id = f"{gene_name}_{pdb_id}"
        
        logging.info(f"Processing {sequence_id}: {len(sequence)} aa")
        
        # Create gene-specific output directory
        gene_output_dir = output_dir / gene_name
        
        # Run BepiPred prediction
        success = predictor.predict_epitopes(sequence, sequence_id, save_dir=gene_output_dir)
        
        if success:
            return {
                'gene': gene_name,
                'pdb_id': pdb_id,
                'sequence_id': sequence_id,
                'sequence_length': len(sequence),
                'fasta_path': str(fasta_path),
                'success': True
            }
        else:
            return {
                'gene': gene_name,
                'pdb_id': pdb_id,
                'sequence_id': sequence_id,
                'sequence_length': len(sequence),
                'fasta_path': str(fasta_path),
                'success': False
            }
            
    except Exception as e:
        logging.error(f"Error processing {fasta_path}: {e}")
        return None

def load_selected_3d_paths(paths_file):
    """Load the list of selected 3D structure FASTA paths"""
    paths = []
    
    try:
        with open(paths_file, 'r') as f:
            for line in f:
                path = line.strip()
                if path:
                    paths.append(path)
        
        logging.info(f"Loaded {len(paths)} 3D structure paths")
        return paths
        
    except Exception as e:
        logging.error(f"Error loading paths file: {e}")
        return []

def create_summary_report(results, output_dir):
    """Create a summary report of all predictions"""
    
    # Create summary statistics
    total_structures = len(results)
    successful_predictions = len([r for r in results if r and r.get('success', False)])
    failed_predictions = total_structures - successful_predictions
    
    # Group by gene
    by_gene = {}
    for result in results:
        if result:
            gene = result['gene']
            if gene not in by_gene:
                by_gene[gene] = []
            by_gene[gene].append(result)
    
    # Create summary
    summary = {
        'analysis_date': datetime.now().isoformat(),
        'total_3d_structures': total_structures,
        'successful_predictions': successful_predictions,
        'failed_predictions': failed_predictions,
        'genes_analyzed': len(by_gene),
        'per_gene_summary': {}
    }
    
    for gene, gene_results in by_gene.items():
        summary['per_gene_summary'][gene] = {
            'total_structures': len(gene_results),
            'successful': len([r for r in gene_results if r.get('success', False)]),
            'structures': [r['pdb_id'] for r in gene_results]
        }
    
    # Save summary JSON
    summary_file = output_dir / "bepipred_3d_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Create human-readable report
    report_file = output_dir / "bepipred_3d_report.txt"
    with open(report_file, 'w') as f:
        f.write("BepiPred 3.0 Analysis of Selected 3D Structures\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"SUMMARY:\n")
        f.write(f"Total 3D structures analyzed: {total_structures}\n")
        f.write(f"Successful predictions: {successful_predictions}\n")
        f.write(f"Failed predictions: {failed_predictions}\n")
        f.write(f"Success rate: {successful_predictions/max(1,total_structures)*100:.1f}%\n\n")
        
        f.write(f"BY GENE:\n")
        f.write("-" * 30 + "\n")
        for gene in sorted(by_gene.keys()):
            gene_data = summary['per_gene_summary'][gene]
            f.write(f"\n{gene}:\n")
            f.write(f"  Structures analyzed: {gene_data['total_structures']}\n")
            f.write(f"  Successful: {gene_data['successful']}\n")
            f.write(f"  PDB IDs: {', '.join(gene_data['structures'])}\n")
    
    logging.info(f"Created summary report: {report_file}")
    logging.info(f"Created summary JSON: {summary_file}")

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        selected_3d_paths_file = snakemake.input.selected_3d_paths
        output_dir = snakemake.output[0]
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        # Get BepiPred config
        bepipred_config = snakemake.config.get('bepipred', {})
        bepipred_path = bepipred_config.get('path', 'tools/BepiPred3.0')
        
    except NameError:
        # Test mode
        import sys
        if len(sys.argv) != 3:
            print("Usage: python predict_epitopes_bepipred_3d_only.py <selected_3d_paths.txt> <output_dir>")
            sys.exit(1)
            
        selected_3d_paths_file = sys.argv[1]
        output_dir = sys.argv[2]
        bepipred_path = "tools/BepiPred3.0"
        bepipred_config = {'method': 'vt_pred', 'top_percentage': 0.2}
    
    logging.info("=== BepiPred 3.0 Epitope Prediction for 3D Structures ===")
    logging.info(f"Selected 3D paths file: {selected_3d_paths_file}")
    logging.info(f"Output directory: {output_dir}")
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load selected 3D structure paths
    paths = load_selected_3d_paths(selected_3d_paths_file)
    
    if not paths:
        logging.error("No 3D structure paths found")
        return
    
    # Initialize BepiPred predictor
    try:
        method = bepipred_config.get('method', 'vt_pred')
        top_percentage = bepipred_config.get('top_percentage', 0.2)
        
        predictor = BepiPredPredictor(
            bepipred_path=bepipred_path, 
            method=method, 
            top_percentage=top_percentage
        )
    except FileNotFoundError as e:
        logging.error(f"BepiPred not found: {e}")
        return
    
    # Process each 3D structure
    results = []
    
    for i, fasta_path in enumerate(paths, 1):
        logging.info(f"\n[{i}/{len(paths)}] Processing: {fasta_path}")
        
        result = predict_epitopes_for_3d_structure(fasta_path, predictor, output_dir)
        if result:
            results.append(result)
    
    # Create summary report
    create_summary_report(results, output_dir)
    
    logging.info(f"\n=== BepiPred 3D Analysis Complete ===")
    logging.info(f"Processed {len(results)} structures")
    logging.info(f"Output saved to: {output_dir}")

if __name__ == "__main__":
    main()