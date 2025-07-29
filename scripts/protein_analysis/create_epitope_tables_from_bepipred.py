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
    
    def __init__(self, min_epitope_length=6, pdb_mapping=None, selected_3d_paths=None):
        """
        Initialize table generator
        
        Args:
            min_epitope_length: Minimum length for linear epitopes (default: 6)
            pdb_mapping: PDB numbering mapping dictionary
            selected_3d_paths: Dictionary mapping gene names to selected 3D structure paths
        """
        self.min_epitope_length = min_epitope_length
        self.pdb_mapping = pdb_mapping or {}
        self.selected_3d_paths = selected_3d_paths or {}
    
    def get_structure_coverage(self, sequence_id):
        """
        Get structure coverage range for a given sequence ID
        
        Args:
            sequence_id: Sequence identifier (e.g., 'bamA_5D0Q')
        
        Returns:
            tuple: (start, end) positions or (None, None) if no structure coverage
        """
        # Extract gene from sequence ID
        if '_' in sequence_id:
            gene = sequence_id.split('_')[0]
        else:
            gene = sequence_id
        
        # Get selected 3D path for this gene
        if gene not in self.selected_3d_paths:
            logging.warning(f"No selected 3D structure found for gene: {gene}")
            return None, None
        
        selected_path = self.selected_3d_paths[gene]
        
        # Get structure coverage from PDB mapping
        if selected_path in self.pdb_mapping:
            pdb_info = self.pdb_mapping[selected_path]
            return pdb_info['start'], pdb_info['end']
        else:
            logging.warning(f"No PDB mapping found for path: {selected_path}")
            return None, None
    
    def filter_epitopes_by_coverage(self, epitopes, start_pos, end_pos):
        """Filter epitopes to only include those within structure coverage"""
        if start_pos is None or end_pos is None:
            return epitopes  # No filtering if no coverage info
        
        filtered_epitopes = []
        for epitope in epitopes:
            # Check if epitope overlaps with structure coverage
            epitope_start = epitope['start']
            epitope_end = epitope['end']
            
            # Keep epitope if it has any overlap with structure coverage
            if epitope_start <= end_pos and epitope_end >= start_pos:
                filtered_epitopes.append(epitope)
            else:
                logging.debug(f"Filtered out epitope {epitope_start}-{epitope_end} (outside structure coverage {start_pos}-{end_pos})")
        
        return filtered_epitopes
    
    def filter_residues_by_coverage(self, residues, start_pos, end_pos):
        """Filter residues to only include those within structure coverage"""
        if start_pos is None or end_pos is None:
            return residues  # No filtering if no coverage info
        
        filtered_residues = []
        for residue in residues:
            if start_pos <= residue['position'] <= end_pos:
                filtered_residues.append(residue)
        
        return filtered_residues
    
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
            
            # Get structure coverage for filtering
            coverage_start, coverage_end = self.get_structure_coverage(sequence_id)
            
            # Filter residues to structure coverage
            filtered_residues = self.filter_residues_by_coverage(all_residues, coverage_start, coverage_end)
            
            # Find epitope regions (consecutive residues with positive predictions)
            linear_epitopes = self._find_linear_epitopes(filtered_residues, sequence)
            
            # Filter epitopes to structure coverage (additional safety check)
            filtered_epitopes = self.filter_epitopes_by_coverage(linear_epitopes, coverage_start, coverage_end)
            
            if coverage_start is not None and coverage_end is not None:
                logging.info(f"Applied structure coverage filtering ({coverage_start}-{coverage_end}): "
                           f"{len(all_residues)} -> {len(filtered_residues)} residues, "
                           f"{len(linear_epitopes)} -> {len(filtered_epitopes)} epitopes")
            
            return {
                'all_residues': filtered_residues,
                'linear_epitopes': filtered_epitopes
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
            # Get PDB information for this sequence
            pdb_info = self.pdb_mapping.get(sequence_id, {})
            chain = pdb_info.get('chain', 'A')
            structure_id = pdb_info.get('structure_id', sequence_id)
            pdb_start = pdb_info.get('start', 1)
            
            # Extract gene name and structure from sequence_id (e.g., bamA_5D0Q -> bamA, 5D0Q)
            if '_' in sequence_id:
                gene_name = sequence_id.split('_')[0]
                structure_name = sequence_id.split('_')[1]
            else:
                gene_name = sequence_id
                structure_name = structure_id
            
            with open(output_file, 'w') as f:
                f.write(f"# Predicted Linear Epitope(s) for {gene_name} ({structure_name}):\n")
                f.write("No.\tChain\tStart\tEnd\tPeptide\tNumber_of_Residues\tScore\n")
                
                for i, epitope in enumerate(epitopes, 1):
                    f.write(f"{i}\t{chain}\t{epitope['start']}\t{epitope['end']}\t{epitope['peptide']}\t{epitope['length']}\t{epitope['avg_score']:.3f}\n")
            
            logging.debug(f"Wrote {len(epitopes)} linear epitopes to {output_file}")
            
        except Exception as e:
            logging.error(f"Error writing linear epitopes table: {e}")
    
    def write_raw_scores_table(self, all_residues, output_file, sequence_id):
        """Write raw scores table with all residue information"""
        try:
            # Get PDB information for this sequence
            pdb_info = self.pdb_mapping.get(sequence_id, {})
            chain = pdb_info.get('chain', 'A')
            structure_id = pdb_info.get('structure_id', sequence_id)
            pdb_start = pdb_info.get('start', 1)
            
            # Extract gene name and structure from sequence_id (e.g., bamA_5D0Q -> bamA, 5D0Q)
            if '_' in sequence_id:
                gene_name = sequence_id.split('_')[0]
                structure_name = sequence_id.split('_')[1]
            else:
                gene_name = sequence_id
                structure_name = structure_id
            
            with open(output_file, 'w') as f:
                f.write(f"# BepiPred 3.0 Raw Scores and Predictions for {gene_name} ({structure_name})\n")
                f.write("No.\tChain\tResidue_Number\tResidue_Name\tScore\tPrediction\n")
                
                for i, residue in enumerate(all_residues, 1):
                    f.write(f"{i}\t{chain}\t{residue['position']}\t{residue['amino_acid']}\t{residue['score']:.3f}\t{residue['prediction']}\n")
            
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

def load_selected_3d_paths(gram_positive_file, gram_negative_file):
    """Load selected 3D structure paths for each gram type"""
    selected_paths = {}
    
    for gram_type, file_path in [('gram_positive', gram_positive_file), ('gram_negative', gram_negative_file)]:
        if file_path and Path(file_path).exists():
            try:
                with open(file_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            # Extract gene name from path: data/protein_structures/bamA/5D0Q_chain_1.fasta -> bamA
                            parts = line.split('/')
                            if len(parts) >= 3 and parts[1] == 'protein_structures':
                                gene = parts[2]
                                selected_paths[gene] = line
            except Exception as e:
                logging.warning(f"Error reading {file_path}: {e}")
    
    logging.info(f"Loaded selected 3D paths for {len(selected_paths)} genes")
    return selected_paths

def load_pdb_numbering_mapping(mapping_file):
    """Load PDB numbering mapping from TSV file"""
    if not mapping_file or not Path(mapping_file).exists():
        logging.warning(f"PDB numbering mapping file not found: {mapping_file}")
        return {}
    
    try:
        df = pd.read_csv(mapping_file, sep='\t')
        mapping = {}
        
        for _, row in df.iterrows():
            gene = row['Gene']
            structure_id = row['Structure_ID']
            chain = row['Chain']
            start = row['Start']
            end = row['End']
            fasta_path = row['FASTA_Path']
            
            # Create key from FASTA path
            key = fasta_path  # Use full path as key for exact matching
            
            mapping[key] = {
                'gene': gene,
                'structure_id': structure_id,
                'chain': chain,
                'start': start,
                'end': end,
                'fasta_path': fasta_path
            }
        
        logging.info(f"Loaded PDB numbering mapping for {len(mapping)} sequences")
        return mapping
    
    except Exception as e:
        logging.error(f"Error loading PDB numbering mapping: {e}")
        return {}

def process_bepipred_outputs(output_dir, min_epitope_length=6, pdb_numbering_mapping=None, selected_3d_paths_positive=None, selected_3d_paths_negative=None):
    """Process all BepiPred outputs in a directory"""
    
    logging.info(f"Processing BepiPred outputs in: {output_dir}")
    
    # Load PDB numbering mapping if provided
    pdb_mapping = load_pdb_numbering_mapping(pdb_numbering_mapping) if pdb_numbering_mapping else {}
    
    # Load selected 3D paths if provided
    selected_3d_paths = load_selected_3d_paths(selected_3d_paths_positive, selected_3d_paths_negative) if (selected_3d_paths_positive or selected_3d_paths_negative) else {}
    
    # Find all CSV files
    csv_files = find_bepipred_csv_files(output_dir)
    
    if not csv_files:
        logging.error(f"No BepiPred CSV files found in {output_dir}")
        return
    
    logging.info(f"Found {len(csv_files)} BepiPred CSV files")
    
    # Initialize table generator
    generator = BepiPredTableGenerator(min_epitope_length=min_epitope_length, pdb_mapping=pdb_mapping, selected_3d_paths=selected_3d_paths)
    
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
        
        # Find corresponding FASTA file or extract sequence from CSV
        fasta_file = find_corresponding_fasta(csv_file)
        if not fasta_file:
            logging.warning(f"Could not find FASTA file for {csv_file}")
            # Extract sequence from CSV data
            try:
                df = pd.read_csv(csv_file)
                if 'Residue' in df.columns:
                    sequence = ''.join(df['Residue'].astype(str))
                    logging.info(f"Extracted sequence from CSV: length {len(sequence)}")
                elif 'AA' in df.columns:
                    sequence = ''.join(df['AA'].astype(str))
                    logging.info(f"Extracted sequence from CSV: length {len(sequence)}")
                else:
                    logging.error(f"No residue column found in CSV {csv_file}")
                    failed += 1
                    continue
            except Exception as e:
                logging.error(f"Error reading CSV {csv_file}: {e}")
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

def load_structure_mapping(mapping_file):
    """Load structure mapping with chain information from final mapping file"""
    try:
        import pandas as pd
        df = pd.read_csv(mapping_file, sep='\t')
        
        # Create mapping dictionaries
        structure_info = {}
        selected_paths = {}
        
        for _, row in df.iterrows():
            gene_name = row['gene_name']
            structure_path = row['structure_path']
            chain_start = row.get('chain_start')
            chain_end = row.get('chain_end')
            
            # Store structure coverage information
            if chain_start is not None and chain_end is not None:
                structure_info[structure_path] = {
                    'start': int(chain_start),
                    'end': int(chain_end)
                }
            
            # Store selected path for each gene (use the first/primary structure)
            if gene_name not in selected_paths:
                selected_paths[gene_name] = structure_path
        
        logging.info(f"Loaded structure mapping for {len(selected_paths)} genes")
        logging.info(f"Structure coverage info for {len(structure_info)} structures")
        
        return structure_info, selected_paths
        
    except Exception as e:
        logging.error(f"Error loading structure mapping: {e}")
        return {}, {}

def process_bepipred_outputs_with_mapping(output_dir, min_epitope_length=6, structure_mapping_file=None):
    """Process BepiPred outputs using the final structure mapping file"""
    
    logging.info(f"Processing BepiPred outputs in: {output_dir}")
    
    # Load structure mapping
    if structure_mapping_file:
        structure_info, selected_paths = load_structure_mapping(structure_mapping_file)
    else:
        logging.warning("No structure mapping file provided")
        structure_info, selected_paths = {}, {}
    
    # Find all CSV files
    csv_files = find_bepipred_csv_files(output_dir)
    
    if not csv_files:
        logging.error(f"No BepiPred CSV files found in {output_dir}")
        return
    
    logging.info(f"Found {len(csv_files)} BepiPred CSV files")
    
    # Initialize table generator
    generator = BepiPredTableGenerator(
        min_epitope_length=min_epitope_length, 
        pdb_mapping=structure_info, 
        selected_3d_paths=selected_paths
    )
    
    processed = 0
    failed = 0
    
    for csv_file in csv_files:
        logging.info(f"Processing: {csv_file}")
        
        # Extract sequence ID from filename
        csv_name = csv_file.stem
        if csv_name.endswith('_raw_output'):
            sequence_id = csv_name.replace('_raw_output', '')
        else:
            sequence_id = csv_name
        
        try:
            # Extract sequence from CSV data
            df = pd.read_csv(csv_file)
            if 'Residue' in df.columns:
                sequence = ''.join(df['Residue'].astype(str))
            elif 'AA' in df.columns:
                sequence = ''.join(df['AA'].astype(str))
            else:
                logging.error(f"No residue column found in CSV {csv_file}")
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
            
        except Exception as e:
            logging.error(f"Failed to process {csv_file}: {e}")
            failed += 1
    
    logging.info(f"Processed: {processed} structures")
    logging.info(f"Failed: {failed} structures")
    logging.info(f"Success rate: {processed/(processed+failed)*100:.1f}%" if (processed+failed) > 0 else "No files processed")

def main():
    """Main function for both command line and Snakemake usage"""
    
    # Check if running from Snakemake
    if 'snakemake' in globals():
        # Running from Snakemake
        bepipred_sentinel = snakemake.input.bepipred_sentinel
        structure_mapping = snakemake.input.structure_mapping
        sentinel_file = snakemake.output.epitope_tables_sentinel
        min_length = snakemake.params.get('min_epitope_length', 6)
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        
        # Derive BepiPred directory from sentinel file path
        bepipred_dir = Path(bepipred_sentinel).parent
        
        logging.info(f"Running epitope table creation for {analysis}_{paramset}")
        logging.info(f"BepiPred directory: {bepipred_dir}")
        logging.info(f"Structure mapping file: {structure_mapping}")
        logging.info(f"Sentinel file: {sentinel_file}")
        
        # Process BepiPred outputs with new structure mapping
        process_bepipred_outputs_with_mapping(bepipred_dir, min_length, structure_mapping)
        
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
        parser.add_argument('--pdb-mapping', help='PDB numbering mapping TSV file (optional)')
        parser.add_argument('--selected-3d-positive', help='Selected 3D paths file for gram positive (optional)')
        parser.add_argument('--selected-3d-negative', help='Selected 3D paths file for gram negative (optional)')
        
        args = parser.parse_args()
        
        if not Path(args.output_dir).exists():
            logging.error(f"Output directory does not exist: {args.output_dir}")
            sys.exit(1)
        
        process_bepipred_outputs(args.output_dir, args.min_length, args.pdb_mapping, args.selected_3d_positive, args.selected_3d_negative)

if __name__ == "__main__":
    main()