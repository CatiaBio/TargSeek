#!/usr/bin/env python3
"""
Simplified MSA FASTA Creation Script
====================================

This script creates MSA FASTA files using pre-selected 3D structures.
The 3D structure selection is now handled by a separate rule.

Creates two types of MSA files:
1. No-3D: Only protein sequences (when enabled)
2. With-3D: Protein sequences + pre-selected 3D structures (when enabled)
"""

import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Maximum sequences to include in MSA
MAX_SEQUENCES_PER_GENE = 100  # Reasonable limit for MAFFT

def extract_genus(species_name):
    """Extract genus from species name"""
    parts = species_name.replace('_', ' ').split()
    return parts[0] if parts else species_name

def select_representative_sequence(fasta_path: Path) -> SeqRecord:
    """Select the sequence closest to median length"""
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    if not sequences:
        raise ValueError(f"No sequences in {fasta_path}")
    
    # If only one sequence, return it
    if len(sequences) == 1:
        return sequences[0]
    
    # Calculate median length
    lengths = sorted([len(seq.seq) for seq in sequences])
    median_len = lengths[len(lengths) // 2]
    
    # Return sequence closest to median length
    return min(sequences, key=lambda rec: abs(len(rec.seq) - median_len))

def select_diverse_sequences(fasta_paths, max_sequences=MAX_SEQUENCES_PER_GENE):
    """
    Select diverse sequences ensuring taxonomic diversity
    
    Returns: List of (fasta_path, SeqRecord) tuples
    """
    # Group by genus
    genus_groups = defaultdict(list)
    
    for fasta_path in fasta_paths:
        species_name = fasta_path.stem
        genus = extract_genus(species_name)
        genus_groups[genus].append(fasta_path)
    
    selected = []
    
    # First round: one sequence per genus
    for genus, paths in genus_groups.items():
        if paths:
            # Randomly select one species from this genus
            selected_path = random.choice(paths)
            try:
                seq = select_representative_sequence(selected_path)
                selected.append((selected_path, seq))
            except Exception as e:
                logging.debug(f"Could not select from {selected_path}: {e}")
    
    # If we need more sequences and have room, add more diversity
    remaining_slots = max_sequences - len(selected)
    if remaining_slots > 0:
        # Get all paths not yet selected
        selected_paths = {path for path, _ in selected}
        remaining_paths = [p for paths in genus_groups.values() 
                          for p in paths if p not in selected_paths]
        
        # Shuffle and select additional sequences
        random.shuffle(remaining_paths)
        for path in remaining_paths[:remaining_slots]:
            try:
                seq = select_representative_sequence(path)
                selected.append((path, seq))
            except Exception as e:
                logging.debug(f"Could not select from {path}: {e}")
    
    return selected[:max_sequences]

def load_filelist(filelist_path: Path) -> list:
    if not filelist_path.exists():
        logging.warning(f"Missing: {filelist_path}")
        return []
    return [Path(line.strip()) for line in filelist_path.read_text().splitlines() if line.strip()]

def load_selected_3d_structures(selected_3d_paths_file: Path):
    """Load pre-selected 3D structure paths"""
    if not selected_3d_paths_file.exists():
        logging.debug(f"No pre-selected 3D structures file: {selected_3d_paths_file}")
        return {}
    
    selected_3d_by_gene = defaultdict(list)
    
    try:
        with open(selected_3d_paths_file, 'r') as f:
            for line in f:
                structure_path = Path(line.strip())
                if structure_path.exists():
                    # Extract gene name from path structure
                    # Expected format: data/protein_structures/{gene_name}/{structure_file}
                    if len(structure_path.parts) >= 3:
                        gene_name = structure_path.parts[-2]  # Parent directory name
                        selected_3d_by_gene[gene_name].append(structure_path)
    except Exception as e:
        logging.warning(f"Error loading selected 3D structures: {e}")
    
    return selected_3d_by_gene

def create_gene_msa(gene_dir: Path, output_dir_msa_no_3d: Path, output_dir_msa_with_3d: Path, 
                   selected_3d_by_gene: dict, enable_no_3d: bool, enable_with_3d: bool):
    """Create MSA FASTA files for a gene using pre-selected 3D structures"""
    gene_name = gene_dir.name
    filelist_path = gene_dir / f"{gene_name}_filelist.txt"
    
    # Load regular sequence paths
    fasta_paths = load_filelist(filelist_path)
    if not fasta_paths:
        logging.warning(f"No sequences found for {gene_name}")
        return None
    
    # Select diverse sequences
    selected_sequences = select_diverse_sequences(fasta_paths)
    
    # Prepare records for no-3D version
    records_no_3d = []
    for fasta_path, seq_record in selected_sequences:
        # Update record metadata
        species_name = fasta_path.stem
        seq_record.id = f"{species_name}|{seq_record.id}"
        seq_record.description = f"{species_name} {seq_record.description}"
        records_no_3d.append(seq_record)
    
    # Write no-3D MSA file only if enabled
    if enable_no_3d and records_no_3d:
        msa_no_3d_out = output_dir_msa_no_3d / f"{gene_name}.fasta"
        output_dir_msa_no_3d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records_no_3d, msa_no_3d_out, "fasta")
        logging.info(f"  No-3D MSA: {len(records_no_3d)} sequences → {msa_no_3d_out.name}")
    elif not enable_no_3d:
        logging.info(f"  Skipping no-3D MSA creation for {gene_name} (disabled)")
    
    # Create with-3D version only if enabled
    if not enable_with_3d:
        logging.info(f"  Skipping with-3D MSA creation for {gene_name} (disabled)")
        return None
    
    records_with_3d = records_no_3d.copy()
    
    # Add pre-selected 3D structures
    selected_3d_paths = selected_3d_by_gene.get(gene_name, [])
    
    if selected_3d_paths:
        logging.info(f"  Adding {len(selected_3d_paths)} pre-selected 3D structures for {gene_name}")
        
        for structure_path in selected_3d_paths:
            try:
                sequences = list(SeqIO.parse(structure_path, "fasta"))
                if sequences:
                    seq_record = sequences[0]
                    # Update record metadata for 3D
                    pdb_info = structure_path.stem
                    seq_record.id = f"3D|{pdb_info}|{seq_record.id}"
                    seq_record.description = f"[3D-STRUCTURE] {pdb_info} {seq_record.description}"
                    records_with_3d.append(seq_record)
                    logging.info(f"    Added 3D structure: {structure_path.name}")
            except Exception as e:
                logging.warning(f"Error loading 3D structure {structure_path}: {e}")
        
        # Write with-3D MSA file
        msa_with_3d_out = output_dir_msa_with_3d / f"{gene_name}.fasta"
        output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records_with_3d, msa_with_3d_out, "fasta")
        logging.info(f"  With-3D MSA: {len(records_with_3d)} sequences ({len(selected_3d_paths)} 3D) → {msa_with_3d_out.name}")
    else:
        # No 3D structures available, copy the no-3D version if it exists
        if records_no_3d:
            logging.info(f"  No 3D structures found for {gene_name}, copying no-3D version")
            msa_with_3d_out = output_dir_msa_with_3d / f"{gene_name}.fasta"
            output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)
            SeqIO.write(records_no_3d, msa_with_3d_out, "fasta")
        else:
            logging.warning(f"  No sequences available for {gene_name}")

def main():
    """Main function with improved logging"""
    try:
        input_dir = Path(snakemake.input.reference_dir)
        selected_3d_paths_file = Path(snakemake.input.selected_3d_paths) if hasattr(snakemake.input, 'selected_3d_paths') else None
        output_dir_msa_no_3d = Path(snakemake.output.msa_no_3d)
        output_dir_msa_with_3d = Path(snakemake.output.msa_with_3d)
        
        # Get processing options
        enable_no_3d = snakemake.params.enable_no_3d
        enable_with_3d = snakemake.params.enable_with_3d
        
        # Get analysis parameters for output file naming
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        # Get config parameters if available
        try:
            global MAX_SEQUENCES_PER_GENE
            MAX_SEQUENCES_PER_GENE = snakemake.config.get('msa', {}).get('max_sequences_per_gene', 100)
        except:
            pass
            
    except NameError:
        # For standalone testing
        input_dir = Path(sys.argv[1])
        selected_3d_paths_file = Path(sys.argv[2]) if len(sys.argv) > 2 else None
        output_dir_msa_no_3d = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("output_no_3d")
        output_dir_msa_with_3d = Path(sys.argv[4]) if len(sys.argv) > 4 else Path("output_with_3d")
        enable_no_3d = True
        enable_with_3d = True
        analysis = "analysis_1"
        paramset = "params_1"
        group = "unknown"

    logging.info("=== MSA FASTA Creation with Pre-selected 3D Structures ===")
    logging.info(f"Reference directory: {input_dir}")
    logging.info(f"Pre-selected 3D paths file: {selected_3d_paths_file}")
    logging.info(f"Output no-3D: {output_dir_msa_no_3d}")
    logging.info(f"Output with-3D: {output_dir_msa_with_3d}")
    logging.info(f"Enable no-3D: {enable_no_3d}")
    logging.info(f"Enable with-3D: {enable_with_3d}")
    logging.info(f"Max sequences per gene: {MAX_SEQUENCES_PER_GENE}")
    
    output_dir_msa_no_3d.mkdir(parents=True, exist_ok=True)
    output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)

    # Load pre-selected 3D structures
    selected_3d_by_gene = {}
    if selected_3d_paths_file and enable_with_3d:
        selected_3d_by_gene = load_selected_3d_structures(selected_3d_paths_file)
        logging.info(f"Loaded pre-selected 3D structures for {len(selected_3d_by_gene)} genes")

    # Process each gene
    gene_dirs = [d for d in sorted(input_dir.iterdir()) if d.is_dir() and not d.name.startswith('_')]
    logging.info(f"\nProcessing {len(gene_dirs)} genes...")
    
    for i, gene_dir in enumerate(gene_dirs, 1):
        logging.info(f"\n[{i}/{len(gene_dirs)}] Processing {gene_dir.name}")
        create_gene_msa(gene_dir, output_dir_msa_no_3d, output_dir_msa_with_3d, 
                       selected_3d_by_gene, enable_no_3d, enable_with_3d)
    
    logging.info(f"\n=== MSA FASTA creation completed successfully! ===")

if __name__ == "__main__":
    main()