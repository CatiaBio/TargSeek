#!/usr/bin/env python3
"""
Select proteins for MSA from shared data directories
===================================================

This script selects representative proteins for MSA using:
- Analysis-specific gene lists from data/{analysis}_{paramset}/genes_species/gram_{group}/
- Shared protein FASTA files from data/proteins_fasta/{gene}/
- Shared 3D structures from data/proteins_3d_structure/{gene}/

Outputs analysis-specific MSA sequences to results/msa_sequences/{analysis}_{paramset}_gram_{group}/
"""

import sys
import json
import logging
import pandas as pd
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_gene_species_lists(gene_lists_dir):
    """
    Load gene-specific species lists for this analysis
    
    Returns:
        dict: {gene_name: [species_list]}
    """
    gene_lists_path = Path(gene_lists_dir)
    gene_species_mapping = {}
    
    if not gene_lists_path.exists():
        logging.error(f"Gene lists directory not found: {gene_lists_dir}")
        return {}
    
    # Find all .txt files (gene species lists)
    gene_files = list(gene_lists_path.glob("*.txt"))
    logging.info(f"Found {len(gene_files)} gene files to process")
    
    for gene_file in gene_files:
        if gene_file.name.startswith("_"):
            continue  # Skip summary files
            
        gene_name = gene_file.stem
        
        # Load species list for this gene
        with open(gene_file, 'r') as f:
            species_list = [line.strip() for line in f if line.strip()]
        
        gene_species_mapping[gene_name] = species_list
        logging.info(f"Loaded {len(species_list)} target species for gene '{gene_name}'")
    
    return gene_species_mapping

def get_target_genes(protein_list_file):
    """Get list of target genes from protein study file"""
    try:
        proteins_df = pd.read_csv(protein_list_file, sep='\t')
        target_genes = proteins_df['gene'].unique().tolist()
        logging.info(f"Target genes from protein list: {target_genes}")
        return target_genes
    except Exception as e:
        logging.error(f"Error reading protein list: {e}")
        return []

def find_available_sequences(gene_name, target_species):
    """
    Find available FASTA sequences for a gene from shared data directory
    
    Returns:
        dict: {species: fasta_file_path}
    """
    gene_fasta_dir = Path("data/proteins_fasta") / gene_name
    available_sequences = {}
    
    if not gene_fasta_dir.exists():
        logging.warning(f"No FASTA directory found for gene {gene_name}")
        return {}
    
    # Look for FASTA files for target species
    for species in target_species:
        safe_species = species.replace(' ', '_').replace('/', '_')
        fasta_file = gene_fasta_dir / f"{safe_species}.fasta"
        
        if fasta_file.exists() and fasta_file.stat().st_size > 0:
            available_sequences[species] = fasta_file
    
    logging.info(f"Found {len(available_sequences)}/{len(target_species)} FASTA sequences for {gene_name}")
    return available_sequences

def find_available_structures(gene_name, target_species):
    """
    Find available 3D structure files for a gene from shared data directory
    
    Returns:
        dict: {species: pdb_file_path}
    """
    gene_structures_dir = Path("data/proteins_3d_structure") / gene_name
    available_structures = {}
    
    if not gene_structures_dir.exists():
        logging.debug(f"No structures directory found for gene {gene_name}")
        return {}
    
    # Look for PDB files for target species
    for species in target_species:
        safe_species = species.replace(' ', '_').replace('/', '_')
        pdb_file = gene_structures_dir / f"{safe_species}.pdb"
        
        if pdb_file.exists() and pdb_file.stat().st_size > 0:
            available_structures[species] = pdb_file
    
    logging.info(f"Found {len(available_structures)}/{len(target_species)} 3D structures for {gene_name}")
    return available_structures

def select_representative_sequence(gene_name, available_sequences, available_structures, priority_species=None):
    """
    Select one representative sequence per species for MSA, preferring 3D structure species
    
    Args:
        gene_name: Name of the gene
        available_sequences: Dict of {species: fasta_file}
        available_structures: Dict of {species: pdb_file}
        priority_species: List of species to prioritize (optional)
        
    Returns:
        dict: {species: selected_fasta_file}
    """
    selected_sequences = {}
    
    # For each species with FASTA sequences
    for species, fasta_file in available_sequences.items():
        # Check if this species also has a 3D structure
        has_structure = species in available_structures
        
        # For now, select all available sequences
        # In the future, this could implement more sophisticated selection logic
        selected_sequences[species] = {
            'fasta_file': fasta_file,
            'has_structure': has_structure,
            'structure_file': available_structures.get(species) if has_structure else None
        }
    
    # Log structure priority
    structure_species = [s for s in selected_sequences.keys() if selected_sequences[s]['has_structure']]
    if structure_species:
        logging.info(f"  Species with 3D structures: {len(structure_species)} ({', '.join(structure_species[:5])}{'...' if len(structure_species) > 5 else ''})")
    
    return selected_sequences

def copy_sequences_for_msa(gene_name, selected_sequences, output_dir):
    """
    Copy selected sequences to MSA directory with proper naming
    
    Args:
        gene_name: Name of the gene
        selected_sequences: Dict from select_representative_sequence
        output_dir: Output directory for MSA sequences
        
    Returns:
        dict: Summary statistics
    """
    
    # Create gene-specific MSA directory
    gene_msa_dir = output_dir / gene_name
    gene_msa_dir.mkdir(parents=True, exist_ok=True)
    
    sequences_copied = 0
    sequences_with_structures = 0
    
    # Copy each selected sequence
    for species, seq_info in selected_sequences.items():
        fasta_file = seq_info['fasta_file']
        has_structure = seq_info['has_structure']
        
        # Read the original FASTA file
        try:
            with open(fasta_file, 'r') as f:
                sequences = list(SeqIO.parse(f, "fasta"))
            
            if not sequences:
                logging.warning(f"  Empty FASTA file for {species}: {fasta_file}")
                continue
            
            # Take the first sequence (should only be one per file)
            seq_record = sequences[0]
            
            # Create new filename for MSA
            safe_species = species.replace(' ', '_').replace('/', '_')
            msa_filename = f"{safe_species}.fasta"
            msa_file = gene_msa_dir / msa_filename
            
            # Modify sequence header to include structure info
            original_id = seq_record.id
            structure_tag = "[3D]" if has_structure else "[SEQ]"
            new_description = f"{seq_record.description} {structure_tag}"
            
            # Create new sequence record
            msa_record = SeqRecord(
                seq_record.seq,
                id=original_id,
                description=new_description
            )
            
            # Write to MSA directory
            with open(msa_file, 'w') as f:
                SeqIO.write(msa_record, f, "fasta")
            
            sequences_copied += 1
            if has_structure:
                sequences_with_structures += 1
                
        except Exception as e:
            logging.error(f"  Error copying sequence for {species}: {e}")
            continue
    
    # Create gene summary file
    summary_file = gene_msa_dir / f"{gene_name}_msa_summary.json"
    summary = {
        "gene": gene_name,
        "total_sequences": sequences_copied,
        "sequences_with_structures": sequences_with_structures,
        "sequences_seq_only": sequences_copied - sequences_with_structures,
        "structure_percentage": (sequences_with_structures / sequences_copied * 100) if sequences_copied > 0 else 0,
        "created_at": datetime.now().isoformat()
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logging.info(f"  MSA sequences: {sequences_copied} total, {sequences_with_structures} with 3D structures ({summary['structure_percentage']:.1f}%)")
    
    return summary

def select_proteins_for_msa_shared(gene_lists_dir, protein_list_file, output_dir, analysis, paramset, group):
    """
    Main function to select proteins for MSA from shared data directories
    """
    
    logging.info(f"=== Selecting Proteins for MSA from Shared Data ===")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    logging.info(f"Gene lists directory: {gene_lists_dir}")
    logging.info(f"Output directory: {output_dir}")
    
    # Load gene-species mappings for this analysis
    gene_species_mapping = load_gene_species_lists(gene_lists_dir)
    
    # Get target genes from protein list
    target_genes = get_target_genes(protein_list_file)
    
    # Filter to only target genes
    filtered_mapping = {gene: species_list for gene, species_list in gene_species_mapping.items() 
                       if gene in target_genes}
    
    if not filtered_mapping:
        logging.error("No overlapping genes found between gene lists and target proteins")
        return
    
    logging.info(f"Processing {len(filtered_mapping)} genes for MSA sequence selection")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Process each gene
    all_summaries = []
    total_sequences = 0
    total_sequences_with_structures = 0
    
    for gene_idx, (gene_name, target_species) in enumerate(filtered_mapping.items(), 1):
        logging.info(f"\n=== Processing Gene {gene_idx}/{len(filtered_mapping)}: {gene_name} ===")
        logging.info(f"Target species: {len(target_species)}")
        
        # Find available sequences and structures
        available_sequences = find_available_sequences(gene_name, target_species)
        available_structures = find_available_structures(gene_name, target_species)
        
        if not available_sequences:
            logging.warning(f"No FASTA sequences available for gene {gene_name}, skipping")
            continue
        
        # Select representative sequences
        selected_sequences = select_representative_sequence(
            gene_name, available_sequences, available_structures, target_species
        )
        
        # Copy sequences to MSA directory
        gene_summary = copy_sequences_for_msa(gene_name, selected_sequences, output_path)
        all_summaries.append(gene_summary)
        
        total_sequences += gene_summary['total_sequences']
        total_sequences_with_structures += gene_summary['sequences_with_structures']
    
    # Create overall summary
    overall_summary = {
        "analysis": analysis,
        "paramset": paramset,
        "group": group,
        "total_genes_processed": len(all_summaries),
        "total_sequences": total_sequences,
        "total_sequences_with_structures": total_sequences_with_structures,
        "total_sequences_seq_only": total_sequences - total_sequences_with_structures,
        "overall_structure_percentage": (total_sequences_with_structures / total_sequences * 100) if total_sequences > 0 else 0,
        "gene_summaries": all_summaries,
        "created_at": datetime.now().isoformat()
    }
    
    summary_file = output_path / "_msa_selection_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(overall_summary, f, indent=2)
    
    # Print final summary
    logging.info(f"\n=== MSA Selection Summary ===")
    logging.info(f"Genes processed: {len(all_summaries)}")
    logging.info(f"Total sequences selected: {total_sequences}")
    logging.info(f"Sequences with 3D structures: {total_sequences_with_structures}")
    logging.info(f"Overall structure coverage: {overall_summary['overall_structure_percentage']:.1f}%")
    logging.info(f"Output directory: {output_path}")
    logging.info(f"Summary file: {summary_file}")

def main():
    """Main function for Snakemake integration"""
    logging.info("=== MSA Sequence Selector (Shared Data) ===")
    
    try:
        # Get inputs from Snakemake
        gene_lists_dir = snakemake.input.gene_lists
        protein_list_file = snakemake.input.protein_list
        output_dir = snakemake.output[0]
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        
        # Select proteins for MSA from shared directories
        select_proteins_for_msa_shared(gene_lists_dir, protein_list_file, output_dir, analysis, paramset, group)
        
        logging.info("MSA sequence selection completed successfully")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        print("Usage: This script should be run from Snakemake")
    except Exception as e:
        logging.error(f"Error during MSA sequence selection: {e}")
        raise

if __name__ == "__main__":
    main()