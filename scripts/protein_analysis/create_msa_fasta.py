#!/usr/bin/env python3
"""
Improved MSA FASTA Creation Script with Intelligent Selection
============================================================

This script creates MSA FASTA files by intelligently selecting sequences from:
- Regular sequences: Prioritizes taxonomic diversity and sequence quality
- 3D structures: Selects best quality structures from different organisms

Improved Selection Criteria:
1. For regular sequences: 
   - Ensures taxonomic diversity (different genera)
   - Selects sequences closest to median length
   - Avoids redundant sequences from same species
2. For 3D structures:
   - Prioritizes high-quality structures from model organisms
   - Ensures organism diversity
   - Limits to reasonable number for alignment
"""

import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Priority organisms for 3D structure selection
PRIORITY_ORGANISMS = [
    'Escherichia coli',
    'Salmonella enterica',
    'Bacillus subtilis',
    'Staphylococcus aureus',
    'Pseudomonas aeruginosa',
    'Mycobacterium tuberculosis',
    'Streptococcus pneumoniae',
    'Klebsiella pneumoniae',
    'Listeria monocytogenes',
    'Enterococcus faecalis'
]

# Maximum sequences to include in MSA
MAX_SEQUENCES_PER_GENE = 100  # Reasonable limit for MAFFT
MAX_3D_STRUCTURES = 10  # Limit 3D structures to avoid overwhelming the alignment

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

def select_best_3d_structures(filelist_3d_path, max_structures=MAX_3D_STRUCTURES):
    """
    Select best 3D structures prioritizing quality and diversity
    
    Returns: List of (fasta_path, SeqRecord) tuples
    """
    if not filelist_3d_path.exists():
        return []
    
    fasta_paths = load_filelist(filelist_3d_path)
    if not fasta_paths:
        return []
    
    # Score each structure
    scored_structures = []
    
    for fasta_path in fasta_paths:
        try:
            # Read first sequence from structure file
            sequences = list(SeqIO.parse(fasta_path, "fasta"))
            if not sequences:
                continue
                
            seq = sequences[0]
            filename = fasta_path.stem
            
            # Base score
            score = 0
            
            # Check for priority organisms
            desc_lower = seq.description.lower()
            filename_lower = filename.lower()
            
            for i, organism in enumerate(PRIORITY_ORGANISMS):
                organism_key = organism.lower().replace(' ', '_')
                if organism_key in desc_lower or organism_key in filename_lower:
                    # Higher priority organisms get higher scores
                    score += (len(PRIORITY_ORGANISMS) - i) * 10
                    break
            
            # Prefer structures with standard PDB naming
            if len(filename.split('_')[0]) == 4:  # Standard PDB ID
                score += 5
            
            # Add some randomness to ensure diversity
            score += random.random()
            
            scored_structures.append((score, fasta_path, seq))
            
        except Exception as e:
            logging.debug(f"Could not score {fasta_path}: {e}")
    
    # Sort by score and select top structures
    scored_structures.sort(key=lambda x: x[0], reverse=True)
    
    # Select structures ensuring some diversity
    selected = []
    seen_prefixes = set()  # Track PDB prefixes to ensure diversity
    
    for score, path, seq in scored_structures:
        if len(selected) >= max_structures:
            break
            
        # Extract PDB prefix
        prefix = path.stem.split('_')[0][:3]  # First 3 chars of PDB ID
        
        # Take high-scoring structures, but ensure some diversity
        if len(selected) < 3 or prefix not in seen_prefixes:
            selected.append((path, seq))
            seen_prefixes.add(prefix)
    
    return selected

def load_filelist(filelist_path: Path) -> list:
    if not filelist_path.exists():
        logging.warning(f"Missing: {filelist_path}")
        return []
    return [Path(line.strip()) for line in filelist_path.read_text().splitlines() if line.strip()]

def create_gene_msa(gene_dir: Path, output_dir_msa_no_3d: Path, output_dir_msa_with_3d: Path):
    """Create MSA FASTA files for a gene with improved selection"""
    gene_name = gene_dir.name
    filelist_path = gene_dir / f"{gene_name}_filelist.txt"
    filelist_3d_path = gene_dir / f"{gene_name}_3d_filelist.txt"
    
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
    
    # Write no-3D MSA file
    if records_no_3d:
        msa_no_3d_out = output_dir_msa_no_3d / f"{gene_name}.fasta"
        output_dir_msa_no_3d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records_no_3d, msa_no_3d_out, "fasta")
        logging.info(f"  No-3D MSA: {len(records_no_3d)} sequences → {msa_no_3d_out.name}")
    
    # Create with-3D version
    records_with_3d = records_no_3d.copy()
    
    # Select and add 3D structures
    selected_3d = select_best_3d_structures(filelist_3d_path)
    
    # Initialize tracking info
    selected_3d_info = {
        'structures': [],
        'paths': []
    }
    
    if selected_3d:
        logging.info(f"  Found {len(selected_3d)} 3D structures for {gene_name}")
        
        for fasta_path, seq_record in selected_3d:
            # Update record metadata for 3D
            pdb_info = fasta_path.stem
            seq_record.id = f"3D|{pdb_info}|{seq_record.id}"
            seq_record.description = f"[3D-STRUCTURE] {pdb_info} {seq_record.description}"
            records_with_3d.append(seq_record)
            
            # Extract PDB ID and organism for tracking
            pdb_id = pdb_info.split('_')[0] if '_' in pdb_info else pdb_info
            organism = extract_organism_from_fasta(fasta_path, seq_record)
            
            # Add to tracking info
            selected_3d_info['structures'].append({
                'gene': gene_name,
                'pdb': pdb_id,
                'bacteria': organism
            })
            selected_3d_info['paths'].append(str(fasta_path))
        
        # Write with-3D MSA file
        msa_with_3d_out = output_dir_msa_with_3d / f"{gene_name}.fasta"
        output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records_with_3d, msa_with_3d_out, "fasta")
        logging.info(f"  With-3D MSA: {len(records_with_3d)} sequences ({len(selected_3d)} 3D) → {msa_with_3d_out.name}")
    else:
        # No 3D structures available, copy the no-3D version
        logging.info(f"  No 3D structures found for {gene_name}")
        msa_with_3d_out = output_dir_msa_with_3d / f"{gene_name}.fasta"
        output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records_no_3d, msa_with_3d_out, "fasta")
        logging.info(f"  Copied no-3D version to with-3D directory")
    
    return selected_3d_info if selected_3d else None


def extract_organism_from_fasta(fasta_path, seq_record):
    """Extract organism name from FASTA header or filename"""
    # Try to extract from description
    desc_lower = seq_record.description.lower()
    
    # Check for common patterns
    if 'organism:' in desc_lower:
        # Pattern: organism:Escherichia coli
        parts = desc_lower.split('organism:')
        if len(parts) > 1:
            organism = parts[1].split('[')[0].split('|')[0].strip()
            return organism.title()
    
    if 'os=' in desc_lower:
        # UniProt pattern: OS=Escherichia coli
        parts = desc_lower.split('os=')
        if len(parts) > 1:
            organism = parts[1].split('ox=')[0].split('gn=')[0].strip()
            return organism.title()
    
    # Check for organism names in description
    for priority_org in PRIORITY_ORGANISMS:
        if priority_org.lower() in desc_lower:
            return priority_org
    
    # Fallback: try to extract from filename
    filename = fasta_path.stem
    # Common pattern: 1ABC_Escherichia_coli
    parts = filename.split('_')
    if len(parts) >= 3:
        # Try to reconstruct organism name
        potential_organism = ' '.join(parts[1:3]).replace('_', ' ')
        return potential_organism.title()
    
    return "Unknown"

def main():
    """Main function with improved logging and 3D tracking"""
    try:
        input_dir = Path(snakemake.input.reference_dir)
        output_dir_msa_no_3d = Path(snakemake.output.msa_no_3d)
        output_dir_msa_with_3d = Path(snakemake.output.msa_with_3d)
        
        # Get analysis parameters for output file naming
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        # Get config parameters if available
        try:
            global MAX_SEQUENCES_PER_GENE, MAX_3D_STRUCTURES
            MAX_SEQUENCES_PER_GENE = snakemake.config.get('msa', {}).get('max_sequences_per_gene', 100)
            MAX_3D_STRUCTURES = snakemake.config.get('msa', {}).get('max_3d_structures_per_gene', 10)
        except:
            pass
            
    except NameError:
        input_dir = Path(sys.argv[1])
        output_dir_msa_no_3d = Path(sys.argv[2])
        output_dir_msa_with_3d = Path(sys.argv[3])
        analysis = "analysis_1"
        paramset = "params_1"
        group = "unknown"

    logging.info("=== MSA FASTA Creation with Intelligent Selection ===")
    logging.info(f"Reference directory: {input_dir}")
    logging.info(f"Output no-3D: {output_dir_msa_no_3d}")
    logging.info(f"Output with-3D: {output_dir_msa_with_3d}")
    logging.info(f"Max sequences per gene: {MAX_SEQUENCES_PER_GENE}")
    logging.info(f"Max 3D structures per gene: {MAX_3D_STRUCTURES}")
    
    output_dir_msa_no_3d.mkdir(parents=True, exist_ok=True)
    output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)

    # Initialize tracking for selected 3D structures
    all_selected_3d = []
    all_3d_paths = []

    # Process each gene
    gene_dirs = [d for d in sorted(input_dir.iterdir()) if d.is_dir() and not d.name.startswith('_')]
    logging.info(f"\nProcessing {len(gene_dirs)} genes...")
    
    for i, gene_dir in enumerate(gene_dirs, 1):
        logging.info(f"\n[{i}/{len(gene_dirs)}] Processing {gene_dir.name}")
        selected_3d_info = create_gene_msa(gene_dir, output_dir_msa_no_3d, output_dir_msa_with_3d)
        
        # Track 3D selections
        if selected_3d_info:
            all_selected_3d.extend(selected_3d_info['structures'])
            all_3d_paths.extend(selected_3d_info['paths'])
    
    # Use Snakemake outputs for tracking files
    tsv_path = Path(snakemake.output.selected_3d_tsv)
    paths_file = Path(snakemake.output.selected_3d_paths)
    
    # Create parent directories for output files
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    paths_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Write selected_3d_fasta.tsv (group-specific)
    if all_selected_3d:
        with open(tsv_path, 'w') as f:
            f.write("gene\tpdb\tbacteria\n")
            for entry in all_selected_3d:
                f.write(f"{entry['gene']}\t{entry['pdb']}\t{entry['bacteria']}\n")
        logging.info(f"\nCreated 3D structure tracking TSV: {tsv_path}")
        logging.info(f"Total 3D structures selected: {len(all_selected_3d)}")
    else:
        # Create empty file if no 3D structures
        tsv_path.touch()
        logging.info(f"\nNo 3D structures found, created empty TSV: {tsv_path}")
    
    # Write selected_3d_paths.txt (group-specific)
    if all_3d_paths:
        with open(paths_file, 'w') as f:
            for path in all_3d_paths:
                f.write(f"{path}\n")
        logging.info(f"Created 3D structure paths file: {paths_file}")
    else:
        # Create empty file if no 3D paths
        paths_file.touch()
        logging.info(f"No 3D paths found, created empty paths file: {paths_file}")
    
    logging.info(f"\n=== MSA FASTA creation completed successfully! ===")

if __name__ == "__main__":
    main()
