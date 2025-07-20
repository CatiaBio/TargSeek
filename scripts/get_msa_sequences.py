#!/usr/bin/env python3
"""
Select Representative Proteins for MSA
======================================

This script processes downloaded protein sequences and selects one representative
protein per species for each gene, preparing them for multiple sequence alignment.
"""

from Bio import SeqIO
import os
import re
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import json

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extract_species_from_description(description):
    """Extract species name from FASTA description with multiple strategies"""
    
    # Strategy 1: Look for species name in brackets [Genus species]
    bracket_match = re.search(r"\[(.*?)\]", description)
    if bracket_match:
        species = bracket_match.group(1)
        # Check if it's a proper binomial name (Genus species)
        if len(species.split()) >= 2:
            # Return first two words (Genus species)
            words = species.split()
            return f"{words[0]} {words[1]}"
    
    # Strategy 2: Look for "ORGANISM=" pattern (sometimes used)
    organism_match = re.search(r"ORGANISM=([^,\]]+)", description)
    if organism_match:
        species = organism_match.group(1).strip()
        words = species.split()
        if len(words) >= 2:
            return f"{words[0]} {words[1]}"
    
    # Strategy 3: Extract from beginning of description (after protein name)
    # Pattern: "protein_name [Genus species]" or "protein_name, Genus species"
    patterns = [
        r"\b([A-Z][a-z]+ [a-z]+)\b",  # Binomial nomenclature pattern
        r"([A-Z][a-z]+\s+[a-z]+)",    # Alternative pattern
    ]
    
    for pattern in patterns:
        matches = re.findall(pattern, description)
        for match in matches:
            words = match.split()
            if len(words) == 2 and words[0][0].isupper() and words[1][0].islower():
                return match
    
    return None

def is_partial_sequence(description):
    """Check if sequence is marked as partial"""
    return "partial" in description.lower()

def calculate_sequence_score(record):
    """Calculate quality score for sequence selection (higher = better)"""
    score = 0
    desc = record.description.lower()
    
    # Base score: sequence length (normalized to 0-100 range)
    score += min(len(record.seq) / 10, 100)  # Longer sequences get higher scores
    
    # Source quality bonuses
    if 'reviewed' in desc or 'swiss-prot' in desc:
        score += 50  # Curated/reviewed sequences
    elif 'uniprot' in desc:
        score += 20  # UniProt but not reviewed
    elif 'refseq' in desc or 'reference' in desc:
        score += 30  # RefSeq/reference sequences
    
    # Annotation quality bonuses
    if 'hypothetical' in desc or 'unknown' in desc or 'putative' in desc:
        score -= 20  # Penalize poorly annotated sequences
    
    # Completeness bonus
    if not is_partial_sequence(record.description):
        score += 30  # Complete sequences preferred
    
    # Strain preference (reference strains often have better annotation)
    strain_indicators = ['type strain', 'atcc', 'dsm ', 'reference', 'model organism']
    if any(indicator in desc for indicator in strain_indicators):
        score += 15
    
    return score

def select_best_sequence(records):
    """Select the best sequence from a list of records for the same species"""
    if not records:
        return None
    
    if len(records) == 1:
        return records[0]
    
    # Calculate scores for all records
    scored_records = [(record, calculate_sequence_score(record)) for record in records]
    
    # Sort by score (highest first) and return the best
    scored_records.sort(key=lambda x: x[1], reverse=True)
    
    best_record, best_score = scored_records[0]
    logging.debug(f"Selected sequence with score {best_score:.1f}: {best_record.description[:100]}...")
    
    return best_record

def filter_sequences_by_length(records, gene_name, max_length_ratio=2.0):
    """Filter sequences to remove outliers that would create poor alignments"""
    if len(records) <= 2:
        return records  # Don't filter if we have very few sequences
    
    lengths = [len(record.seq) for record in records]
    median_length = np.median(lengths)
    
    # Remove sequences that are too short or too long compared to median
    min_length = median_length / max_length_ratio
    max_length = median_length * max_length_ratio
    
    filtered_records = []
    removed_count = 0
    
    for record in records:
        seq_length = len(record.seq)
        if min_length <= seq_length <= max_length:
            filtered_records.append(record)
        else:
            removed_count += 1
            logging.debug(f"  Removed {record.id}: length {seq_length} outside range [{min_length:.0f}-{max_length:.0f}]")
    
    if removed_count > 0:
        logging.info(f"  {gene_name}: Filtered out {removed_count} length outliers, kept {len(filtered_records)} sequences")
    
    return filtered_records if filtered_records else records  # Return original if all filtered out

def process_gene_sequences(gene_name, gene_dir, structures_dir, output_dir, expected_species=None):
    """Process all sequences for a single gene, including 3D structures first"""
    
    gene_path = Path(gene_dir) / gene_name
    structures_gene_path = Path(structures_dir) / gene_name
    
    if not gene_path.exists():
        logging.warning(f"Gene directory not found: {gene_path}")
        return 0, None  # Return tracking info as None
    
    species_to_records = {}
    total_sequences = 0
    unmatched_descriptions = []
    
    # Initialize 3D structure tracking
    structure_tracking = {
        "gene": gene_name,
        "has_3d_structures": False,
        "available_3d_files": [],
        "selected_3d_file": None,
        "selected_3d_id": None,
        "selected_3d_species": None,
        "reason_no_3d": None
    }
    
    logging.info(f"Processing gene {gene_name}")
    
    # FIRST: Process 3D structure sequences if available
    if structures_gene_path.exists():
        structure_files = list(structures_gene_path.glob("*.fasta"))
        structure_tracking["available_3d_files"] = [f.name for f in structure_files]
        
        if structure_files:
            structure_tracking["has_3d_structures"] = True
            logging.info(f"  Found {len(structure_files)} 3D structure sequences")
            for fasta_file in structure_files:
                try:
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        total_sequences += 1
                        # For 3D structures, use "3D_Structure" as the "species" to ensure they're included first
                        species = "3D_Structure"
                        if species not in species_to_records:
                            species_to_records[species] = []
                        species_to_records[species].append(record)
                        logging.info(f"    Added 3D structure: {record.id}")
                except Exception as e:
                    logging.warning(f"Error processing 3D structure {fasta_file}: {e}")
        else:
            # Check for no_structures_found.txt file
            no_structures_file = structures_gene_path / "no_structures_found.txt"
            if no_structures_file.exists():
                structure_tracking["reason_no_3d"] = "No 3D structures found in PDB (confirmed by marker file)"
                logging.info(f"  No 3D structures available for {gene_name} (confirmed by marker file)")
            else:
                structure_tracking["reason_no_3d"] = "No 3D structure FASTA files in directory"
                logging.info(f"  No 3D structure FASTA files found for {gene_name}")
    else:
        structure_tracking["reason_no_3d"] = "3D structures directory not found"
        logging.info(f"  3D structures directory not found for {gene_name}")
    
    # SECOND: Process regular protein sequences from species
    logging.info(f"  Processing species sequences from {gene_path}")
    
    # Process all FASTA files for this gene
    for fasta_file in gene_path.glob("*.fasta"):
        try:
            logging.info(f"  Processing file: {fasta_file.name}")
            file_sequences = 0
            
            for record in SeqIO.parse(fasta_file, "fasta"):
                total_sequences += 1
                file_sequences += 1
                
                # Debug: show first few descriptions
                if file_sequences <= 3:
                    logging.info(f"    Example description: {record.description}")
                
                # Skip partial sequences initially
                if is_partial_sequence(record.description):
                    continue
                
                species = extract_species_from_description(record.description)
                if species:
                    if species not in species_to_records:
                        species_to_records[species] = []
                    species_to_records[species].append(record)
                else:
                    # Store unmatched for debugging
                    if len(unmatched_descriptions) < 5:  # Limit output
                        unmatched_descriptions.append(record.description)
                
            logging.info(f"    Found {file_sequences} sequences in {fasta_file.name}")
                
        except Exception as e:
            logging.warning(f"Error processing {fasta_file}: {e}")
    
    # Debug output
    if species_to_records:
        logging.info(f"  Found species: {list(species_to_records.keys())}")
    else:
        logging.warning(f"  No species matches found for {gene_name}")
        if unmatched_descriptions:
            logging.warning("  Examples of unmatched descriptions:")
            for desc in unmatched_descriptions:
                logging.warning(f"    - {desc}")
    
    # If no species found, try more aggressive matching including partial sequences
    if not species_to_records:
        logging.warning(f"  No complete sequences matched. Trying with partial sequences...")
        
        for fasta_file in gene_path.glob("*.fasta"):
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    # Try species extraction without filtering partial sequences
                    species = extract_species_from_description(record.description)
                    if species:
                        if species not in species_to_records:
                            species_to_records[species] = []
                        species_to_records[species].append(record)
                    else:
                        # Last resort: try to get species from filename
                        filename = fasta_file.stem  # Remove .fasta extension
                        # Clean filename (remove underscores, etc.)
                        potential_species = filename.replace('_', ' ')
                        # Check if it looks like a species name
                        words = potential_species.split()
                        if len(words) >= 2 and words[0][0].isupper() and words[1][0].islower():
                            species = f"{words[0]} {words[1]}"
                            if species not in species_to_records:
                                species_to_records[species] = []
                            species_to_records[species].append(record)
                            logging.info(f"    Using filename-based species: {species}")
                            
            except Exception as e:
                logging.warning(f"Error in fallback processing {fasta_file}: {e}")
        
        if species_to_records:
            logging.info(f"  Fallback found species: {list(species_to_records.keys())}")
    
    if not species_to_records:
        logging.error(f"  Still no species found for {gene_name} after all strategies")
        return 0
    
    # Select best sequence per species, ensuring only one 3D structure closest to average length
    selected_records = []
    
    # First, add best sequence per regular species (to calculate average length)
    species_records = []
    for species, records in species_to_records.items():
        if species == "3D_Structure":
            continue  # Process later
            
        best_record = select_best_sequence(records)
        if best_record:
            # Clean up the header to include species info
            best_record.id = f"{gene_name}_{species.replace(' ', '_')}"
            best_record.description = f"{gene_name} [{species}] | {best_record.description}"
            species_records.append(best_record)
    
    # Calculate average length of species sequences
    if species_records:
        avg_length = sum(len(record.seq) for record in species_records) / len(species_records)
        logging.info(f"  Average length of species sequences: {avg_length:.1f}")
        
        # Select the 3D structure sequence closest to average length
        selected_3d_species = None
        if "3D_Structure" in species_to_records:
            structure_records = species_to_records["3D_Structure"]
            best_structure = None
            min_diff = float('inf')
            
            for record in structure_records:
                length_diff = abs(len(record.seq) - avg_length)
                if length_diff < min_diff:
                    min_diff = length_diff
                    best_structure = record
            
            if best_structure:
                # Extract species from 3D structure description to check for duplicates
                selected_3d_species = extract_species_from_description(best_structure.description)
                
                # Update tracking info
                structure_tracking["selected_3d_id"] = best_structure.id
                structure_tracking["selected_3d_species"] = selected_3d_species
                # Find which file this structure came from
                for file in structure_tracking["available_3d_files"]:
                    structure_tracking["selected_3d_file"] = file  # For now, just use the last one
                
                # Keep original ID for 3D structure but add gene prefix
                best_structure.id = f"{gene_name}_{best_structure.id}"
                best_structure.description = f"{gene_name} [3D_Structure] | {best_structure.description}"
                selected_records.append(best_structure)
                logging.info(f"  Selected 3D structure {best_structure.id} (length: {len(best_structure.seq)}, diff from avg: {min_diff:.1f})")
                
                if selected_3d_species:
                    logging.info(f"  3D structure is from species: {selected_3d_species}")
        
        # Add species sequences, but skip if already represented by 3D structure
        excluded_count = 0
        for record in species_records:
            # Extract species from the record's description or ID
            record_species = extract_species_from_description(record.description)
            
            # Check if this species is already represented by the selected 3D structure
            if selected_3d_species and record_species == selected_3d_species:
                excluded_count += 1
                logging.info(f"  Skipping {record_species} sequence - already represented by 3D structure")
                continue
            
            selected_records.append(record)
        
        if excluded_count > 0:
            logging.info(f"  Excluded {excluded_count} species sequences to avoid duplication with 3D structure")
    
    elif "3D_Structure" in species_to_records:
        # No species sequences, just add the first 3D structure
        structure_record = species_to_records["3D_Structure"][0]
        
        # Update tracking info
        structure_tracking["selected_3d_id"] = structure_record.id
        structure_tracking["selected_3d_species"] = extract_species_from_description(structure_record.description)
        if structure_tracking["available_3d_files"]:
            structure_tracking["selected_3d_file"] = structure_tracking["available_3d_files"][0]
        
        structure_record.id = f"{gene_name}_{structure_record.id}"
        structure_record.description = f"{gene_name} [3D_Structure] | {structure_record.description}"
        selected_records.append(structure_record)
        logging.info(f"  Added 3D structure sequence (no species sequences for comparison)")
    
    else:
        selected_records.extend(species_records)
    
    # Filter sequences by length similarity for better alignment
    if len(selected_records) > 1:
        selected_records = filter_sequences_by_length(selected_records, gene_name)
    
    # Write merged FASTA file
    if selected_records:
        output_file = Path(output_dir) / f"{gene_name}.fasta"
        with open(output_file, 'w') as f:
            SeqIO.write(selected_records, f, "fasta")
        
        logging.info(f"Gene {gene_name}: {len(selected_records)} sequences from {len(species_to_records)} species (processed {total_sequences} total)")
        return len(selected_records), structure_tracking
    
    return 0, structure_tracking

def main():
    """Main function for Snakemake integration"""
    logging.info("MSA Sequence Selection")
    logging.info("="*50)
    
    try:
        # Get inputs from Snakemake
        protein_dir = snakemake.input.protein_dir
        protein_list_file = snakemake.input.protein_list
        structures_dir = snakemake.input.structures_dir
        output_dir = snakemake.output[0]
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        logging.info(f"Protein directory: {protein_dir}")
        logging.info(f"3D structures directory: {structures_dir}")
        logging.info(f"Output directory: {output_dir}")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        protein_dir = "results/protein_fasta/analysis_1_params_1_gram_positive"
        protein_list_file = "results/proteins_to_study/analysis_1_params_1_gram_positive.tsv"
        structures_dir = "results/3d_structures/analysis_1_params_1_gram_positive"
        output_dir = "results/msa_sequences/test"
        
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the list of genes to process
    try:
        proteins_df = pd.read_csv(protein_list_file, sep='\t')
        genes_to_process = proteins_df['gene'].unique().tolist()
        logging.info(f"Found {len(genes_to_process)} genes to process")
    except Exception as e:
        logging.error(f"Error reading protein list: {e}")
        return
    
    # Check if protein directory exists
    if not Path(protein_dir).exists():
        logging.error(f"Protein directory not found: {protein_dir}")
        return
    
    # Process each gene
    total_sequences_selected = 0
    successful_genes = 0
    all_structure_tracking = []
    
    for gene in genes_to_process:
        logging.info(f"Processing gene: {gene}")
        sequences_selected, structure_tracking = process_gene_sequences(gene, protein_dir, structures_dir, output_dir)
        
        if structure_tracking:
            all_structure_tracking.append(structure_tracking)
        
        if sequences_selected > 0:
            successful_genes += 1
            total_sequences_selected += sequences_selected
        else:
            logging.warning(f"No sequences selected for gene {gene}")
    
    # Generate 3D structure tracking report
    report_file = Path(output_dir) / "3d_structure_selection_report.json"
    with open(report_file, 'w') as f:
        json.dump(all_structure_tracking, f, indent=2)
    
    # Generate human-readable summary
    summary_file = Path(output_dir) / "3d_structure_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("3D Structure Selection Summary\n")
        f.write("=" * 50 + "\n\n")
        
        genes_with_3d = [t for t in all_structure_tracking if t["has_3d_structures"]]
        genes_without_3d = [t for t in all_structure_tracking if not t["has_3d_structures"]]
        
        f.write(f"Total genes processed: {len(all_structure_tracking)}\n")
        f.write(f"Genes with 3D structures: {len(genes_with_3d)}\n")
        f.write(f"Genes without 3D structures: {len(genes_without_3d)}\n\n")
        
        if genes_with_3d:
            f.write("GENES WITH 3D STRUCTURES SELECTED FOR MSA:\n")
            f.write("-" * 50 + "\n")
            for tracking in genes_with_3d:
                f.write(f"Gene: {tracking['gene']}\n")
                f.write(f"  Selected 3D file: {tracking['selected_3d_file']}\n")
                f.write(f"  Selected 3D ID: {tracking['selected_3d_id']}\n")
                f.write(f"  Species: {tracking['selected_3d_species']}\n")
                f.write(f"  Available files: {', '.join(tracking['available_3d_files'])}\n")
                f.write("\n")
        
        if genes_without_3d:
            f.write("GENES WITHOUT 3D STRUCTURES:\n")
            f.write("-" * 50 + "\n")
            for tracking in genes_without_3d:
                f.write(f"Gene: {tracking['gene']}\n")
                f.write(f"  Reason: {tracking['reason_no_3d']}\n")
                f.write("\n")
    
    # Summary
    logging.info(f"\nSummary:")
    logging.info(f"Genes requested: {len(genes_to_process)}")
    logging.info(f"Genes with sequences: {successful_genes}")
    logging.info(f"Total sequences selected for MSA: {total_sequences_selected}")
    logging.info(f"3D structure reports saved to:")
    logging.info(f"  - {report_file}")
    logging.info(f"  - {summary_file}")
    
    if successful_genes == 0:
        logging.error("No sequences were selected for any gene!")
        logging.info("Available genes in protein directory:")
        protein_path = Path(protein_dir)
        if protein_path.exists():
            for gene_dir in protein_path.iterdir():
                if gene_dir.is_dir():
                    logging.info(f"  - {gene_dir.name}")

if __name__ == "__main__":
    main()