#!/usr/bin/env python3
"""
Create Structure MSA Files
==========================

This script creates structure MSA files for MAFFT processing by combining:
- Main MSA sequences (already available as main_msa_sequences)
- Selected 3D structures with length-based optimization

This is STAGE 3 of the 3-stage MSA workflow:
1. Main MSA created ✅ (main_msa_sequences)
2. 3D structures selected from main MSA ✅
3. Create structure MSA files (this script) → structure_msa_sequences

The structure MSA files contain main sequences + 3D structures for alignment.
"""

import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import shutil

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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

def create_structure_msa_for_gene(gene_name: str, main_msa_file: Path, 
                                 output_structure_dir: Path, selected_3d_by_gene: dict):
    """Create structure MSA file for a single gene (main MSA + 3D structures)"""
    
    if not main_msa_file.exists():
        logging.warning(f"Main MSA file not found for {gene_name}: {main_msa_file}")
        return None
    
    # Load main MSA sequences
    try:
        main_sequences = list(SeqIO.parse(main_msa_file, "fasta"))
    except Exception as e:
        logging.error(f"Error loading main MSA for {gene_name}: {e}")
        return None
    
    if not main_sequences:
        logging.warning(f"No sequences in main MSA for {gene_name}")
        return None
    
    logging.info(f"  Loaded {len(main_sequences)} sequences from main MSA")
    
    # Start with main MSA sequences
    final_sequences = main_sequences.copy()
    
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
                    final_sequences.append(seq_record)
                    logging.info(f"    Added 3D structure: {structure_path.name}")
            except Exception as e:
                logging.warning(f"Error loading 3D structure {structure_path}: {e}")
        
        # Write structure MSA file
        output_file = output_structure_dir / f"{gene_name}.fasta"
        output_structure_dir.mkdir(parents=True, exist_ok=True)
        SeqIO.write(final_sequences, output_file, "fasta")
        logging.info(f"  Structure MSA: {len(final_sequences)} sequences ({len(selected_3d_paths)} 3D) → {output_file.name}")
        
        return {
            'gene': gene_name,
            'main_sequences': len(main_sequences),
            '3d_structures_added': len(selected_3d_paths),
            'final_sequences': len(final_sequences)
        }
    else:
        # No 3D structures available, copy the main MSA to structure output
        if main_sequences:
            logging.info(f"  No 3D structures found for {gene_name}, copying main MSA to structure output")
            output_file = output_structure_dir / f"{gene_name}.fasta"
            output_structure_dir.mkdir(parents=True, exist_ok=True)
            SeqIO.write(main_sequences, output_file, "fasta")
            logging.info(f"  Structure MSA: {len(main_sequences)} sequences (0 3D) → {output_file.name}")
            
            return {
                'gene': gene_name,
                'main_sequences': len(main_sequences),
                '3d_structures_added': 0,
                'final_sequences': len(main_sequences)
            }
        else:
            logging.warning(f"  No sequences available for {gene_name}")
            return None

def main():
    """Main function for creating structure MSA files"""
    try:
        main_msa_dir = Path(snakemake.input.sequences_dir)
        selected_3d_paths_file = Path(snakemake.input.selected_3d_paths) if hasattr(snakemake.input, 'selected_3d_paths') else None
        output_structure_dir = Path(snakemake.output.sequences_with_structure_dir)
        
        # Get analysis parameters for output file naming
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
    except NameError:
        # For standalone testing
        main_msa_dir = Path(sys.argv[1])
        selected_3d_paths_file = Path(sys.argv[2]) if len(sys.argv) > 2 else None
        output_structure_dir = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("output_structure")
        analysis = "analysis_1"
        paramset = "params_1"
        group = "unknown"

    logging.info("=== Structure MSA Files Creation (STAGE 3) ===")
    logging.info(f"Main MSA directory: {main_msa_dir}")
    logging.info(f"Pre-selected 3D paths file: {selected_3d_paths_file}")
    logging.info(f"Output structure MSA: {output_structure_dir}")
    
    output_structure_dir.mkdir(parents=True, exist_ok=True)

    # Load pre-selected 3D structures
    selected_3d_by_gene = {}
    if selected_3d_paths_file:
        selected_3d_by_gene = load_selected_3d_structures(selected_3d_paths_file)
        logging.info(f"Loaded pre-selected 3D structures for {len(selected_3d_by_gene)} genes")
    
    # Process each gene based on main MSA files
    main_msa_files = list(main_msa_dir.glob("*.fasta"))
    if not main_msa_files:
        logging.error(f"No main MSA files found in {main_msa_dir}")
        return
    
    logging.info(f"\nProcessing {len(main_msa_files)} genes for structure MSA creation...")
    
    results = []
    for i, main_msa_file in enumerate(sorted(main_msa_files), 1):
        gene_name = main_msa_file.stem
        logging.info(f"\n[{i}/{len(main_msa_files)}] Processing {gene_name}")
        
        result = create_structure_msa_for_gene(
            gene_name, main_msa_file, output_structure_dir, selected_3d_by_gene
        )
        
        if result:
            results.append(result)
    
    # Summary statistics
    if results:
        total_main_sequences = sum(r['main_sequences'] for r in results)
        total_3d_added = sum(r['3d_structures_added'] for r in results)
        genes_with_3d = sum(1 for r in results if r['3d_structures_added'] > 0)
        
        logging.info(f"\n=== Structure MSA Creation Summary ===")
        logging.info(f"Genes processed: {len(results)}")
        logging.info(f"Total main sequences: {total_main_sequences}")
        logging.info(f"Total 3D structures added: {total_3d_added}")
        logging.info(f"Genes with 3D structures: {genes_with_3d}/{len(results)} ({genes_with_3d/len(results)*100:.1f}%)")
        logging.info(f"Structure MSA files created: {len(results)}")
        
        # Show genes with most 3D structures added
        genes_with_most_3d = sorted(results, key=lambda x: x['3d_structures_added'], reverse=True)[:5]
        if genes_with_most_3d[0]['3d_structures_added'] > 0:
            logging.info(f"\nGenes with most 3D structures:")
            for r in genes_with_most_3d:
                if r['3d_structures_added'] > 0:
                    logging.info(f"  {r['gene']}: +{r['3d_structures_added']} 3D structures ({r['main_sequences']} → {r['final_sequences']} total)")
    else:
        logging.warning("No genes processed successfully")
    
    logging.info(f"\n=== Structure MSA files creation completed! ===")

if __name__ == "__main__":
    main()