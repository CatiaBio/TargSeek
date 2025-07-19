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
from pathlib import Path
import logging

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

def select_best_sequence(records):
    """Select the best sequence from a list of records for the same species"""
    if not records:
        return None
    
    # Filter out partial sequences first
    complete_records = [r for r in records if not is_partial_sequence(r.description)]
    
    # If we have complete sequences, use those; otherwise use all
    candidates = complete_records if complete_records else records
    
    # Select the longest sequence
    best_record = max(candidates, key=lambda x: len(x.seq))
    return best_record

def process_gene_sequences(gene_name, gene_dir, output_dir, expected_species=None):
    """Process all sequences for a single gene"""
    
    gene_path = Path(gene_dir) / gene_name
    if not gene_path.exists():
        logging.warning(f"Gene directory not found: {gene_path}")
        return 0
    
    species_to_records = {}
    total_sequences = 0
    unmatched_descriptions = []
    
    logging.info(f"Processing gene {gene_name} in {gene_path}")
    
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
    
    # Select best sequence per species
    selected_records = []
    for species, records in species_to_records.items():
        best_record = select_best_sequence(records)
        if best_record:
            # Clean up the header to include species info
            best_record.id = f"{gene_name}_{species.replace(' ', '_')}"
            best_record.description = f"{gene_name} [{species}] | {best_record.description}"
            selected_records.append(best_record)
    
    # Write merged FASTA file
    if selected_records:
        output_file = Path(output_dir) / f"{gene_name}.fasta"
        with open(output_file, 'w') as f:
            SeqIO.write(selected_records, f, "fasta")
        
        logging.info(f"Gene {gene_name}: {len(selected_records)} sequences from {len(species_to_records)} species (processed {total_sequences} total)")
        return len(selected_records)
    
    return 0

def main():
    """Main function for Snakemake integration"""
    logging.info("MSA Sequence Selection")
    logging.info("="*50)
    
    try:
        # Get inputs from Snakemake
        protein_dir = snakemake.input.protein_dir
        protein_list_file = snakemake.input.protein_list
        output_dir = snakemake.output[0]
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"Processing analysis={analysis}, paramset={paramset}, group={group}")
        logging.info(f"Protein directory: {protein_dir}")
        logging.info(f"Output directory: {output_dir}")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        protein_dir = "results/proteins/analysis_1_params_1_gram_positive"
        protein_list_file = "results/proteins_to_study/analysis_1_params_1_gram_positive.tsv"
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
    
    for gene in genes_to_process:
        logging.info(f"Processing gene: {gene}")
        sequences_selected = process_gene_sequences(gene, protein_dir, output_dir)
        
        if sequences_selected > 0:
            successful_genes += 1
            total_sequences_selected += sequences_selected
        else:
            logging.warning(f"No sequences selected for gene {gene}")
    
    # Summary
    logging.info(f"\nSummary:")
    logging.info(f"Genes requested: {len(genes_to_process)}")
    logging.info(f"Genes with sequences: {successful_genes}")
    logging.info(f"Total sequences selected for MSA: {total_sequences_selected}")
    
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