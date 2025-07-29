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
import pandas as pd
import re

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

def load_structure_organism_mapping(mapping_file: Path):
    """Load FASTA structure mapping to get organism information for each structure"""
    if not mapping_file.exists():
        logging.warning(f"FASTA structure mapping file not found: {mapping_file}")
        return {}
    
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        
        # Create a mapping from FASTA path to organism information
        structure_organism_map = {}
        for _, row in mapping_df.iterrows():
            fasta_path = row['fasta_path']
            species = row['species']
            structure_organism_map[fasta_path] = species
            
        logging.info(f"Loaded organism mapping for {len(structure_organism_map)} structures")
        return structure_organism_map
        
    except Exception as e:
        logging.error(f"Error loading structure organism mapping: {e}")
        return {}

def extract_organism_from_species(species_text):
    """Extract organism name from species text like 'ESCHERICHIA COLI (562)'"""
    if not species_text:
        return None
    
    # Extract organism name before parentheses and normalize
    organism = re.sub(r'\s*\([^)]*\)', '', species_text).strip()
    return organism.lower() if organism else None

def should_exclude_sequence(seq_record, excluded_organisms):
    """Check if a sequence should be excluded based on organism information"""
    if not excluded_organisms:
        return False
    
    # Extract organism information from sequence description/header
    # Typical format: ">GeneSymbol|GeneName|Species_name|Taxonomy_ID|..."
    description = seq_record.description.lower()
    seq_id = seq_record.id.lower()
    
    # Check if any excluded organism appears in the sequence description or ID
    for organism in excluded_organisms:
        if organism in description or organism in seq_id:
            return True
    
    return False

def create_structure_msa_for_gene(gene_name: str, main_msa_file: Path, 
                                 output_structure_dir: Path, selected_3d_by_gene: dict,
                                 structure_organism_map: dict = None):
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
    
    # Identify organisms to exclude based on selected 3D structures
    excluded_organisms = set()
    selected_3d_paths = selected_3d_by_gene.get(gene_name, [])
    
    if selected_3d_paths and structure_organism_map:
        for structure_path in selected_3d_paths:
            structure_path_str = str(structure_path)
            if structure_path_str in structure_organism_map:
                species_info = structure_organism_map[structure_path_str]
                organism = extract_organism_from_species(species_info)
                if organism:
                    excluded_organisms.add(organism)
                    logging.info(f"  Will exclude sequences from organism: {organism} (structure: {structure_path.name})")
    
    # Filter main MSA sequences to exclude organisms with selected 3D structures
    if excluded_organisms:
        original_count = len(main_sequences)
        filtered_sequences = []
        excluded_count = 0
        
        for seq_record in main_sequences:
            if should_exclude_sequence(seq_record, excluded_organisms):
                excluded_count += 1
                logging.debug(f"    Excluding sequence: {seq_record.id} (organism match)")
            else:
                filtered_sequences.append(seq_record)
        
        final_sequences = filtered_sequences
        logging.info(f"  Filtered sequences: {original_count} → {len(final_sequences)} (excluded {excluded_count} duplicate organisms)")
    else:
        # No organisms to exclude
        final_sequences = main_sequences.copy()
    
    # Add pre-selected 3D structures
    
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
            'final_sequences': len(final_sequences),
            'excluded_organisms': excluded_count if excluded_organisms else 0,
            'excluded_organism_names': list(excluded_organisms) if excluded_organisms else []
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
                'final_sequences': len(main_sequences),
                'excluded_organisms': 0,
                'excluded_organism_names': []
            }
        else:
            logging.warning(f"  No sequences available for {gene_name}")
            return None

def main():
    """Main function for creating structure MSA files"""
    try:
        main_msa_dir = Path(snakemake.input.sequences_dir)
        selected_3d_paths_file = Path(snakemake.input.selected_3d_paths) if hasattr(snakemake.input, 'selected_3d_paths') else None
        fasta_structure_mapping_file = Path(snakemake.input.fasta_structure_mapping) if hasattr(snakemake.input, 'fasta_structure_mapping') else None
        output_structure_dir = Path(snakemake.output.sequences_with_structure_dir)
        
        # Get analysis parameters for output file naming
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
    except NameError:
        # For standalone testing
        main_msa_dir = Path(sys.argv[1])
        selected_3d_paths_file = Path(sys.argv[2]) if len(sys.argv) > 2 else None
        fasta_structure_mapping_file = Path(sys.argv[3]) if len(sys.argv) > 3 else None
        output_structure_dir = Path(sys.argv[4]) if len(sys.argv) > 4 else Path("output_structure")
        analysis = "analysis_1"
        paramset = "params_1"
        group = "unknown"

    logging.info("=== Structure MSA Files Creation (STAGE 3) ===")
    logging.info(f"Main MSA directory: {main_msa_dir}")
    logging.info(f"Pre-selected 3D paths file: {selected_3d_paths_file}")
    logging.info(f"FASTA structure mapping file: {fasta_structure_mapping_file}")
    logging.info(f"Output structure MSA: {output_structure_dir}")
    
    output_structure_dir.mkdir(parents=True, exist_ok=True)

    # Load pre-selected 3D structures
    selected_3d_by_gene = {}
    if selected_3d_paths_file:
        selected_3d_by_gene = load_selected_3d_structures(selected_3d_paths_file)
        logging.info(f"Loaded pre-selected 3D structures for {len(selected_3d_by_gene)} genes")
    
    # Load structure organism mapping
    structure_organism_map = {}
    if fasta_structure_mapping_file:
        structure_organism_map = load_structure_organism_mapping(fasta_structure_mapping_file)
        logging.info(f"Loaded organism mapping for {len(structure_organism_map)} structures")
    
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
            gene_name, main_msa_file, output_structure_dir, selected_3d_by_gene, structure_organism_map
        )
        
        if result:
            results.append(result)
    
    # Summary statistics
    if results:
        total_main_sequences = sum(r['main_sequences'] for r in results)
        total_3d_added = sum(r['3d_structures_added'] for r in results)
        total_excluded_sequences = sum(r['excluded_organisms'] for r in results)
        genes_with_3d = sum(1 for r in results if r['3d_structures_added'] > 0)
        genes_with_exclusions = sum(1 for r in results if r['excluded_organisms'] > 0)
        
        logging.info(f"\n=== Structure MSA Creation Summary ===")
        logging.info(f"Genes processed: {len(results)}")
        logging.info(f"Total main sequences: {total_main_sequences}")
        logging.info(f"Total 3D structures added: {total_3d_added}")
        logging.info(f"Total sequences excluded (duplicate organisms): {total_excluded_sequences}")
        logging.info(f"Genes with 3D structures: {genes_with_3d}/{len(results)} ({genes_with_3d/len(results)*100:.1f}%)")
        logging.info(f"Genes with organism exclusions: {genes_with_exclusions}/{len(results)} ({genes_with_exclusions/len(results)*100:.1f}%)")
        logging.info(f"Structure MSA files created: {len(results)}")
        
        # Show genes with most 3D structures added
        genes_with_most_3d = sorted(results, key=lambda x: x['3d_structures_added'], reverse=True)[:5]
        if genes_with_most_3d[0]['3d_structures_added'] > 0:
            logging.info(f"\nGenes with most 3D structures:")
            for r in genes_with_most_3d:
                if r['3d_structures_added'] > 0:
                    excluded_info = f" (excluded {r['excluded_organisms']} {', '.join(r['excluded_organism_names'][:2])})" if r['excluded_organisms'] > 0 else ""
                    logging.info(f"  {r['gene']}: +{r['3d_structures_added']} 3D structures ({r['main_sequences']} → {r['final_sequences']} total){excluded_info}")
        
        # Show genes with most organism exclusions
        genes_with_most_exclusions = sorted(results, key=lambda x: x['excluded_organisms'], reverse=True)[:5]
        if genes_with_most_exclusions[0]['excluded_organisms'] > 0:
            logging.info(f"\nGenes with most organism exclusions:")
            for r in genes_with_most_exclusions:
                if r['excluded_organisms'] > 0:
                    organisms_str = ', '.join(r['excluded_organism_names'][:3])
                    if len(r['excluded_organism_names']) > 3:
                        organisms_str += f" (+{len(r['excluded_organism_names'])-3} more)"
                    logging.info(f"  {r['gene']}: excluded {r['excluded_organisms']} sequences from [{organisms_str}]")
    else:
        logging.warning("No genes processed successfully")
    
    logging.info(f"\n=== Structure MSA files creation completed! ===")

if __name__ == "__main__":
    main()