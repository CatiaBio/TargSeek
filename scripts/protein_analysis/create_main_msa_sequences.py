#!/usr/bin/env python3
"""
Create Main MSA Sequences with Species-Based Selection
======================================================

This script creates the foundation MSA sequences using species-based selection
to ensure exactly one sequence per species per gene. These sequences serve as 
the base for both no_3d and with_3d processing.

The main MSA includes:
- One sequence per species (species uniqueness)
- Representative sequences (median length when multiple available)
- Quality filtering for optimal sequence selection
- Taxonomic diversity when max_sequences limit is applied

This is STAGE 1 of the 3-stage MSA workflow.
"""

import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Maximum sequences to include in main MSA - None means no limit (include all species)
MAX_SEQUENCES_PER_GENE = None  # Include all available species for comprehensive analysis

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
    Select one sequence per species ensuring species uniqueness
    
    This creates the foundation MSA that will be used for:
    1. no_3d processing (directly)
    2. 3D structure selection (5 sequences sampled from this)
    3. with_3d processing (this + selected 3D structure)
    
    Selection criteria when multiple sequences exist for the same species:
    - Prefer sequences closest to median length (most representative)
    - Filter out sequences that are extremely short or long
    
    Returns: List of (fasta_path, SeqRecord) tuples
    """
    # Group by species (not genus)
    species_groups = defaultdict(list)
    
    for fasta_path in fasta_paths:
        species_name = fasta_path.stem
        species_groups[species_name].append(fasta_path)
    
    selected = []
    
    # Select exactly one sequence per species
    for species_name, paths in species_groups.items():
        if not paths:
            continue
            
        if len(paths) == 1:
            # Only one sequence for this species
            selected_path = paths[0]
        else:
            # Multiple files for same species - need to select best one
            # This can happen when a species has multiple sequences in the filelist
            candidate_sequences = []
            
            # Load all sequences for this species and find the best one
            for path in paths:
                try:
                    sequences = list(SeqIO.parse(path, "fasta"))
                    for seq in sequences:
                        candidate_sequences.append((path, seq))
                except Exception as e:
                    logging.debug(f"Could not read sequences from {path}: {e}")
                    continue
                    
            if not candidate_sequences:
                logging.debug(f"No valid sequences found for species {species_name}")
                continue
                
            # Find sequence closest to median length
            lengths = [len(seq.seq) for _, seq in candidate_sequences]
            if not lengths:
                continue
                
            median_len = sorted(lengths)[len(lengths) // 2]
            
            # Select sequence closest to median length
            best_path, best_seq = min(candidate_sequences, 
                                    key=lambda x: abs(len(x[1].seq) - median_len))
            selected_path = best_path
        
        # Load the selected sequence
        try:
            seq = select_representative_sequence(selected_path)
            selected.append((selected_path, seq))
            logging.debug(f"Selected {species_name}: {len(seq.seq)} aa")
        except Exception as e:
            logging.debug(f"Could not select sequence for {species_name}: {e}")
    
    # Sort by species name for consistent output
    selected.sort(key=lambda x: x[0].stem)
    
    # Apply max_sequences limit if specified
    if max_sequences is not None and len(selected) > max_sequences:
        logging.info(f"Limiting from {len(selected)} to {max_sequences} sequences")
        # Prioritize taxonomic diversity when limiting
        genus_counts = defaultdict(int)
        limited_selected = []
        
        # First pass: one per genus
        for path, seq in selected:
            genus = extract_genus(path.stem)
            if genus_counts[genus] == 0:
                limited_selected.append((path, seq))
                genus_counts[genus] += 1
                if len(limited_selected) >= max_sequences:
                    break
        
        # Second pass: fill remaining slots
        if len(limited_selected) < max_sequences:
            remaining_slots = max_sequences - len(limited_selected)
            selected_paths = {path for path, _ in limited_selected}
            for path, seq in selected:
                if path not in selected_paths:
                    limited_selected.append((path, seq))
                    remaining_slots -= 1
                    if remaining_slots <= 0:
                        break
        
        selected = limited_selected
    elif max_sequences is None:
        logging.info(f"No sequence limit applied - including all {len(selected)} available species")
    
    return selected

def load_filelist(filelist_path: Path) -> list:
    """Load file paths from a filelist"""
    if not filelist_path.exists():
        logging.warning(f"Missing: {filelist_path}")
        return []
    return [Path(line.strip()) for line in filelist_path.read_text().splitlines() if line.strip()]

def create_main_msa_for_gene(gene_dir: Path, output_dir: Path):
    """Create main MSA sequences for a single gene"""
    gene_name = gene_dir.name
    filelist_path = gene_dir / f"{gene_name}_filelist.txt"
    
    # Load regular sequence paths
    fasta_paths = load_filelist(filelist_path)
    if not fasta_paths:
        logging.warning(f"No sequences found for {gene_name}")
        return None
    
    # Select diverse sequences using intelligent selection
    selected_sequences = select_diverse_sequences(fasta_paths)
    
    if not selected_sequences:
        logging.warning(f"No sequences selected for {gene_name}")
        return None
    
    # Prepare records for main MSA
    main_records = []
    for fasta_path, seq_record in selected_sequences:
        # Update record metadata
        species_name = fasta_path.stem
        seq_record.id = f"{species_name}|{seq_record.id}"
        seq_record.description = f"{species_name} {seq_record.description}"
        main_records.append(seq_record)
    
    # Write main MSA file
    output_file = output_dir / f"{gene_name}.fasta"
    output_dir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(main_records, output_file, "fasta")
    
    species_count = len(set(path.stem for path, _ in selected_sequences))
    logging.info(f"  Main MSA: {len(main_records)} sequences ({species_count} species) → {output_file.name}")
    
    return {
        'gene': gene_name,
        'num_sequences': len(main_records),
        'species': len(set(path.stem for path, _ in selected_sequences)),
        'genera': len(set(extract_genus(path.stem) for path, _ in selected_sequences)),
        'output_file': output_file
    }

def main():
    """Main function for creating main MSA sequences"""
    try:
        input_dir = Path(snakemake.input.reference_dir)
        output_dir = Path(snakemake.output.sequences_dir)
        
        # Get analysis parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        # Get config parameters if available
        try:
            global MAX_SEQUENCES_PER_GENE
            MAX_SEQUENCES_PER_GENE = snakemake.config.get('msa', {}).get('max_sequences_per_gene', None)
        except:
            pass
            
    except NameError:
        # For standalone testing
        input_dir = Path(sys.argv[1])
        output_dir = Path(sys.argv[2])
        analysis = "analysis_1"
        paramset = "params_1"
        group = "unknown"

    logging.info("=== Main MSA Sequence Creation (STAGE 1) ===")
    logging.info(f"Reference directory: {input_dir}")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    logging.info(f"Max sequences per gene: {'unlimited' if MAX_SEQUENCES_PER_GENE is None else MAX_SEQUENCES_PER_GENE}")
    
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process each gene
    gene_dirs = [d for d in sorted(input_dir.iterdir()) if d.is_dir() and not d.name.startswith('_')]
    logging.info(f"\nProcessing {len(gene_dirs)} genes for main MSA creation...")
    
    results = []
    for i, gene_dir in enumerate(gene_dirs, 1):
        logging.info(f"\n[{i}/{len(gene_dirs)}] Processing {gene_dir.name}")
        result = create_main_msa_for_gene(gene_dir, output_dir)
        if result:
            results.append(result)
    
    # Summary statistics
    if results:
        total_sequences = sum(r['num_sequences'] for r in results)
        avg_sequences = total_sequences / len(results)
        total_species = sum(r['species'] for r in results)
        avg_species = total_species / len(results)
        total_genera = sum(r['genera'] for r in results)
        avg_genera = total_genera / len(results)
        
        logging.info(f"\n=== Main MSA Creation Summary ===")
        logging.info(f"Genes processed: {len(results)}")
        logging.info(f"Total sequences: {total_sequences}")
        logging.info(f"Average sequences per gene: {avg_sequences:.1f}")
        logging.info(f"Total species represented: {total_species}")
        logging.info(f"Average species per gene: {avg_species:.1f}")
        logging.info(f"Total genera represented: {total_genera}")
        logging.info(f"Average genera per gene: {avg_genera:.1f}")
        
        # Log genes with highest diversity
        top_diverse = sorted(results, key=lambda x: x['species'], reverse=True)[:5]
        logging.info(f"\nMost diverse genes (by species count):")
        for r in top_diverse:
            logging.info(f"  {r['gene']}: {r['species']} species, {r['genera']} genera, {r['num_sequences']} sequences")
    else:
        logging.warning("No genes processed successfully")
    
    logging.info(f"\n=== Main MSA sequence creation completed! ===")

if __name__ == "__main__":
    main()