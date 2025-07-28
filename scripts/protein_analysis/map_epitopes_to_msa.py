#!/usr/bin/env python3
"""
Map epitope positions from structure sequences to MSA alignment positions
=========================================================================

This script maps epitope positions from the ungapped 3D structure sequence
to the corresponding positions in the multiple sequence alignment (MSA).

Usage:
    python map_epitopes_to_msa.py <gene> <epitope_start> <epitope_end>
    
Example:
    python scripts/protein_analysis/map_epitopes_to_msa.py bamA 49 54
"""

from Bio import SeqIO
from pathlib import Path
import argparse
import sys

def find_structure_sequence_in_msa(msa_file, structure_id):
    """
    Find the 3D structure sequence in the MSA alignment
    
    Args:
        msa_file: Path to MSA FASTA file
        structure_id: Structure identifier (e.g., '5D0Q')
    
    Returns:
        Bio.SeqRecord or None: The structure sequence record from MSA
    """
    try:
        sequences = list(SeqIO.parse(msa_file, "fasta"))
        
        # Look for sequence with 3D structure identifier
        for seq_record in sequences:
            if structure_id in seq_record.description and "3D" in seq_record.description:
                return seq_record
                
        # Fallback: look for any sequence containing the structure ID
        for seq_record in sequences:
            if structure_id in seq_record.description:
                return seq_record
                
        return None
        
    except Exception as e:
        print(f"Error reading MSA file {msa_file}: {e}")
        return None

def map_ungapped_to_gapped_positions(gapped_sequence, start_pos, end_pos):
    """
    Map positions from ungapped sequence to gapped (MSA) sequence
    
    Args:
        gapped_sequence: The aligned sequence with gaps (-)
        start_pos: Start position in ungapped sequence (1-based)
        end_pos: End position in ungapped sequence (1-based)
    
    Returns:
        tuple: (msa_start, msa_end) positions in the gapped alignment (1-based)
    """
    ungapped_pos = 0
    msa_start = None
    msa_end = None
    
    for i, char in enumerate(gapped_sequence):
        if char != '-':  # Not a gap
            ungapped_pos += 1
            
            # Mark start position
            if ungapped_pos == start_pos:
                msa_start = i + 1  # Convert to 1-based
                
            # Mark end position
            if ungapped_pos == end_pos:
                msa_end = i + 1  # Convert to 1-based
                break
    
    return msa_start, msa_end

def extract_epitope_from_msa(gapped_sequence, msa_start, msa_end):
    """
    Extract epitope sequence from MSA alignment (may contain gaps)
    
    Args:
        gapped_sequence: The aligned sequence with gaps
        msa_start: Start position in MSA (1-based)
        msa_end: End position in MSA (1-based)
    
    Returns:
        str: Epitope sequence from MSA (may contain gaps)
    """
    if msa_start is None or msa_end is None:
        return None
        
    # Convert to 0-based indexing for slicing
    return str(gapped_sequence[msa_start-1:msa_end])

def get_msa_context(msa_file, msa_start, msa_end, context_size=10):
    """
    Get alignment context around epitope region for all sequences
    
    Args:
        msa_file: Path to MSA FASTA file
        msa_start: Start position in MSA (1-based)
        msa_end: End position in MSA (1-based)
        context_size: Number of positions to show before/after epitope
    
    Returns:
        list: List of (seq_id, context_sequence) tuples
    """
    if msa_start is None or msa_end is None:
        return []
        
    try:
        sequences = list(SeqIO.parse(msa_file, "fasta"))
        context_sequences = []
        
        # Calculate context boundaries
        context_start = max(1, msa_start - context_size)
        context_end = min(len(sequences[0].seq) if sequences else 0, msa_end + context_size)
        
        for seq_record in sequences:
            seq_id = seq_record.description.split('|')[0].replace('>', '')
            context_seq = str(seq_record.seq[context_start-1:context_end])
            
            # Mark epitope region with brackets
            epitope_start_in_context = msa_start - context_start
            epitope_end_in_context = msa_end - context_start + 1
            
            marked_seq = (context_seq[:epitope_start_in_context] + 
                         '[' + context_seq[epitope_start_in_context:epitope_end_in_context] + ']' +
                         context_seq[epitope_end_in_context:])
            
            context_sequences.append((seq_id, marked_seq))
            
        return context_sequences
        
    except Exception as e:
        print(f"Error getting MSA context: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description='Map epitope positions to MSA alignment')
    parser.add_argument('gene', help='Gene name (e.g., bamA)')
    parser.add_argument('epitope_start', type=int, help='Epitope start position (1-based)')
    parser.add_argument('epitope_end', type=int, help='Epitope end position (1-based)')
    parser.add_argument('--gram-type', choices=['positive', 'negative'], default='negative',
                       help='Gram stain type (default: negative)')
    parser.add_argument('--context', type=int, default=10,
                       help='Context size around epitope (default: 10)')
    parser.add_argument('--analysis', default='analysis1', help='Analysis name')
    parser.add_argument('--paramset', default='params1', help='Parameter set')
    
    args = parser.parse_args()
    
    # Construct paths
    msa_file = Path(f"results/{args.analysis}_{args.paramset}/protein_analysis/sequences_with_structure/msa_sequences/gram_{args.gram_type}/{args.gene}.fasta")
    selected_3d_file = Path(f"results/{args.analysis}_{args.paramset}/protein_analysis/sequences_with_structure/selected_3d_paths_gram_{args.gram_type}.txt")
    
    if not msa_file.exists():
        print(f"Error: MSA file not found: {msa_file}")
        sys.exit(1)
    
    # Get selected 3D structure for this gene
    structure_id = None
    if selected_3d_file.exists():
        with open(selected_3d_file, 'r') as f:
            for line in f:
                if f"protein_structures/{args.gene}/" in line:
                    # Extract structure ID from path like: data/protein_structures/bamA/5D0Q_chain_1.fasta
                    structure_id = line.strip().split('/')[-1].split('_')[0].replace('.fasta', '')
                    break
    
    if not structure_id:
        print(f"Error: Could not find 3D structure for gene {args.gene}")
        sys.exit(1)
    
    print(f"Mapping epitope for {args.gene} (structure: {structure_id})")
    print(f"Epitope position: {args.epitope_start}-{args.epitope_end}")
    print(f"MSA file: {msa_file}")
    print()
    
    # Find structure sequence in MSA
    structure_seq_record = find_structure_sequence_in_msa(msa_file, structure_id)
    if not structure_seq_record:
        print(f"Error: Could not find structure sequence {structure_id} in MSA")
        sys.exit(1)
    
    print(f"Found structure sequence: {structure_seq_record.description}")
    print()
    
    # Map positions
    msa_start, msa_end = map_ungapped_to_gapped_positions(
        structure_seq_record.seq, args.epitope_start, args.epitope_end
    )
    
    if msa_start is None or msa_end is None:
        print(f"Error: Could not map epitope positions {args.epitope_start}-{args.epitope_end}")
        sys.exit(1)
    
    print(f"Mapped MSA positions: {msa_start}-{msa_end}")
    
    # Extract epitope from MSA
    epitope_msa = extract_epitope_from_msa(structure_seq_record.seq, msa_start, msa_end)
    print(f"Epitope in MSA: {epitope_msa}")
    print()
    
    # Show context in alignment
    print(f"Alignment context (±{args.context} positions, epitope in brackets):")
    print("-" * 80)
    
    context_sequences = get_msa_context(msa_file, msa_start, msa_end, args.context)
    
    for seq_id, context_seq in context_sequences[:10]:  # Show first 10 sequences
        print(f"{seq_id[:30]:<30} {context_seq}")
    
    if len(context_sequences) > 10:
        print(f"... and {len(context_sequences) - 10} more sequences")
    
    print("-" * 80)
    print(f"Total sequences in alignment: {len(context_sequences)}")

if __name__ == "__main__":
    main()