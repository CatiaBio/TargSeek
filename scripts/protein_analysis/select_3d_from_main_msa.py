#!/usr/bin/env python3
"""
Select 3D Structures from Main MSA
===================================

This script selects the best 3D structures for each gene by:
1. Randomly selecting 5 sequences from the main MSA
2. Creating MSA alignments with available 3D structure sequences using MAFFT
3. Calculating similarity scores using the specified similarity method
4. Selecting the best 3D structure based on min_similarity_score threshold

This is STAGE 2 of the 3-stage MSA workflow.
"""

import sys
import logging
import random
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import pandas as pd
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_main_msa_sequences(main_msa_file: Path, num_sequences: int = 5) -> list:
    """Load and randomly select sequences from main MSA"""
    if not main_msa_file.exists():
        logging.warning(f"Main MSA file not found: {main_msa_file}")
        return []

    sequences = list(SeqIO.parse(main_msa_file, "fasta"))
    if not sequences:
        return []

    # Randomly select sequences (or all if fewer than requested)
    selected_count = min(num_sequences, len(sequences))
    selected_sequences = random.sample(sequences, selected_count)

    logging.debug(f"Selected {selected_count} sequences from {len(sequences)} total")
    return selected_sequences

def load_3d_structure_sequences(structure_dir: Path) -> list:
    """Load all 3D structure FASTA sequences for a gene"""
    if not structure_dir.exists():
        logging.debug(f"3D structure directory not found: {structure_dir}")
        return []

    structure_sequences = []
    for fasta_file in structure_dir.glob("*.fasta"):
        try:
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            for seq in sequences:
                # Add metadata about the structure file
                seq.structure_file = fasta_file.name
                seq.structure_path = fasta_file
                structure_sequences.append(seq)
        except Exception as e:
            logging.debug(f"Could not read {fasta_file}: {e}")

    logging.debug(f"Found {len(structure_sequences)} 3D structure sequences")
    return structure_sequences

def run_mafft_alignment(sequences: list) -> dict:
    """Run MAFFT alignment on a set of sequences"""
    if len(sequences) < 2:
        return {}

    # Create temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
        SeqIO.write(sequences, temp_input, "fasta")
        temp_input_path = temp_input.name

    try:
        # Run MAFFT
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_output:
            temp_output_path = temp_output.name

        cmd = ['mafft', '--auto', '--quiet', temp_input_path]
        with open(temp_output_path, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            logging.warning(f"MAFFT failed: {result.stderr}")
            return {}

        # Load aligned sequences
        aligned_sequences = list(SeqIO.parse(temp_output_path, "fasta"))
        alignment_dict = {seq.id: seq for seq in aligned_sequences}

        # Clean up temp files
        Path(temp_input_path).unlink()
        Path(temp_output_path).unlink()

        return alignment_dict

    except Exception as e:
        logging.warning(f"MAFFT alignment failed: {e}")
        # Clean up temp files
        try:
            Path(temp_input_path).unlink()
            if 'temp_output_path' in locals():
                Path(temp_output_path).unlink()
        except:
            pass
        return {}


def calculate_length_based_score(structure_seq, main_sequences: list) -> float:
    """
    Calculate length-based similarity score comparing structure length to MSA average
    Returns a score between 0.0 and 1.0 where 1.0 is perfect length match
    """
    if not main_sequences:
        return 0.0

    # Calculate average length of main MSA sequences (ungapped)
    main_lengths = []
    for seq in main_sequences:
        # Remove gaps to get actual sequence length
        ungapped_length = len(str(seq.seq).replace('-', ''))
        main_lengths.append(ungapped_length)

    if not main_lengths:
        return 0.0

    average_msa_length = sum(main_lengths) / len(main_lengths)

    # Get ungapped structure sequence length
    structure_length = len(str(structure_seq.seq).replace('-', ''))

    # Calculate length similarity score
    # Score = 1.0 - (absolute_difference / max_length)
    # This gives higher scores to structures with lengths closer to the MSA average
    max_length = max(average_msa_length, structure_length)
    if max_length == 0:
        return 0.0

    length_difference = abs(structure_length - average_msa_length)
    length_score = max(0.0, 1.0 - (length_difference / max_length))

    return length_score

def calculate_msa_based_similarity_score(aligned_sequences: dict, structure_id: str) -> float:
    """
    Calculate similarity score based on MSA alignment positions
    Uses local alignment scoring which is most suitable for protein homology detection
    """
    if structure_id not in aligned_sequences:
        return 0.0

    structure_seq = aligned_sequences[structure_id]
    other_sequences = [seq for seq_id, seq in aligned_sequences.items() if seq_id != structure_id]

    if not other_sequences:
        return 0.0

    # Calculate alignment-based similarity score
    total_positions = len(str(structure_seq.seq))
    similarity_scores = []

    for other_seq in other_sequences:
        matches = 0
        valid_positions = 0

        # Compare aligned positions
        for i, (aa1, aa2) in enumerate(zip(str(structure_seq.seq), str(other_seq.seq))):
            if aa1 != '-' and aa2 != '-':  # Both have amino acids at this position
                valid_positions += 1
                if aa1 == aa2:
                    matches += 1

        # Calculate similarity for this pair
        if valid_positions > 0:
            pair_similarity = matches / valid_positions
            similarity_scores.append(pair_similarity)

    # Return average similarity across all main MSA sequences
    return sum(similarity_scores) / len(similarity_scores) if similarity_scores else 0.0

def calculate_coverage_score(aligned_sequences: dict, structure_id: str) -> float:
    """
    Calculate coverage score based on how much of the structure sequence aligns
    with non-gap positions in the MSA sequences
    """
    if structure_id not in aligned_sequences:
        return 0.0

    structure_seq = aligned_sequences[structure_id]
    other_sequences = [seq for seq_id, seq in aligned_sequences.items() if seq_id != structure_id]

    if not other_sequences:
        return 0.0

    structure_str = str(structure_seq.seq)
    total_structure_positions = len(structure_str.replace('-', ''))

    if total_structure_positions == 0:
        return 0.0

    # Count positions where structure has amino acid and at least one MSA sequence also has amino acid
    covered_positions = 0

    for i, structure_aa in enumerate(structure_str):
        if structure_aa != '-':  # Structure has amino acid at this position
            # Check if any MSA sequence has amino acid at this position
            has_coverage = any(str(seq.seq)[i] != '-' for seq in other_sequences if i < len(str(seq.seq)))
            if has_coverage:
                covered_positions += 1

    coverage_score = covered_positions / total_structure_positions
    return coverage_score

def select_best_3d_structure(gene_name: str, main_msa_dir: Path, structure_dir: Path,
                           similarity_method: str, min_similarity_score: float,
                           max_structures: int = 1, scoring_weights: dict = None,
                           enable_prefilter: bool = False, max_length_deviation: float = 0.5) -> list:
    """Select the best 3D structure for a gene based on similarity to main MSA"""

    # Default scoring weights if not provided
    if scoring_weights is None:
        scoring_weights = {'length': 0.5, 'similarity': 0.3, 'coverage': 0.2}

    # Step 1: Load main MSA sequences (randomly select 5) - SAME SET FOR ALL STRUCTURE EVALUATIONS
    main_msa_file = main_msa_dir / f"{gene_name}.fasta"
    main_sequences = load_main_msa_sequences(main_msa_file, num_sequences=5)

    if not main_sequences:
        logging.warning(f"No main MSA sequences found for {gene_name}")
        return []

    logging.info(f"  Using {len(main_sequences)} sequences from main MSA for all structure evaluations")

    # Step 2: Load 3D structure sequences
    structure_sequences = load_3d_structure_sequences(structure_dir)

    if not structure_sequences:
        logging.debug(f"No 3D structure sequences found for {gene_name}")
        return []

    # Optional: Length-based pre-filtering
    if enable_prefilter:
        # Calculate MSA average length
        msa_lengths = [len(str(seq.seq).replace('-', '')) for seq in main_sequences]
        avg_msa_length = sum(msa_lengths) / len(msa_lengths) if msa_lengths else 0

        # Filter structures by length deviation
        filtered_structures = []
        for structure_seq in structure_sequences:
            struct_length = len(str(structure_seq.seq).replace('-', ''))
            length_deviation = abs(struct_length - avg_msa_length) / avg_msa_length if avg_msa_length > 0 else 1.0

            if length_deviation <= max_length_deviation:
                filtered_structures.append(structure_seq)
            else:
                logging.debug(f"    Pre-filtered out {structure_seq.structure_file}: "
                            f"length {struct_length} vs MSA avg {avg_msa_length:.1f} "
                            f"(deviation: {length_deviation:.2f})")

        logging.debug(f"  Pre-filtering: {len(filtered_structures)}/{len(structure_sequences)} structures passed")
        structure_sequences = filtered_structures

        if not structure_sequences:
            logging.debug(f"No structures passed length pre-filtering for {gene_name}")
            return []

    # Step 3: Evaluate each 3D structure using the SAME main MSA sequences
    structure_scores = []

    for i, structure_seq in enumerate(structure_sequences, 1):
        logging.debug(f"    Evaluating structure {i}/{len(structure_sequences)}: {structure_seq.structure_file}")

        # Create MSA: 5 selected sequences + 1 structure sequence
        msa_sequences = main_sequences + [structure_seq]

        # Run MAFFT alignment on this combination (5 + 1 = 6 sequences total)
        alignment = run_mafft_alignment(msa_sequences)

        if len(alignment) != len(msa_sequences):
            logging.debug(f"      Alignment failed for structure {structure_seq.structure_file}")
            continue

        # Find structure sequence ID in alignment
        structure_id = None
        for seq_id in alignment.keys():
            # Try to match by original sequence ID or file name
            if (hasattr(structure_seq, 'id') and structure_seq.id in seq_id) or \
               structure_seq.structure_file.replace('.fasta', '') in seq_id:
                structure_id = seq_id
                break

        if not structure_id:
            # Fallback: find the sequence that doesn't match main sequences
            main_seq_ids = [seq.id for seq in main_sequences]
            for seq_id in alignment.keys():
                if not any(main_id in seq_id for main_id in main_seq_ids):
                    structure_id = seq_id
                    break

        if not structure_id:
            logging.debug(f"      Could not identify structure sequence in alignment")
            continue

        # Step 4: Calculate multiple scores for comprehensive evaluation
        # Length-based score (compare structure length to MSA average)
        length_score = calculate_length_based_score(structure_seq, main_sequences)

        # MSA-based similarity score (sequence identity)
        similarity_score = calculate_msa_based_similarity_score(alignment, structure_id)

        # Coverage score (how well structure aligns with MSA)
        coverage_score = calculate_coverage_score(alignment, structure_id)

        # Combined score: weighted average of all metrics using configurable weights
        combined_score = (scoring_weights['length'] * length_score +
                         scoring_weights['similarity'] * similarity_score +
                         scoring_weights['coverage'] * coverage_score)

        logging.debug(f"      Length score: {length_score:.3f}")
        logging.debug(f"      Similarity score: {similarity_score:.3f}")
        logging.debug(f"      Coverage score: {coverage_score:.3f}")
        logging.debug(f"      Combined score: {combined_score:.3f}")

        # Step 5: Apply minimum threshold to combined score
        if combined_score >= min_similarity_score:
            # Extract metadata from structure sequence
            pdb_id = structure_seq.structure_file.replace('.fasta', '').split('_')[0]
            organism = getattr(structure_seq, 'description', 'Unknown').split(' ')[0:2]
            organism_name = ' '.join(organism) if len(organism) >= 2 else 'Unknown'

            structure_info = {
                'gene': gene_name,
                'pdb': pdb_id,
                'chain': '',  # Could be extracted from description if available
                'bacteria': organism_name,
                'length_score': length_score,
                'similarity_score': similarity_score,
                'coverage_score': coverage_score,
                'combined_score': combined_score,
                'consensus_similarity': combined_score,  # For backward compatibility
                'total_score': combined_score,
                'structure_length': len(str(structure_seq.seq).replace('-', '')),
                'msa_avg_length': sum(len(str(seq.seq).replace('-', '')) for seq in main_sequences) / len(main_sequences),
                'file_name': structure_seq.structure_file,
                'file_path': str(structure_seq.structure_path),
                'selection_method': 'length_similarity_coverage'
            }

            structure_scores.append(structure_info)
        else:
            logging.debug(f"      Below threshold ({min_similarity_score}), skipping")

    # Step 6: Sort by combined score (descending) and select the best structure
    structure_scores.sort(key=lambda x: x['combined_score'], reverse=True)
    selected_structures = structure_scores[:max_structures]

    logging.info(f"  Selected {len(selected_structures)}/{len(structure_sequences)} 3D structures for {gene_name}")
    if selected_structures:
        best_struct = selected_structures[0]
        logging.info(f"  Best structure: {best_struct['pdb']} (combined={best_struct['combined_score']:.3f}, "
                    f"length={best_struct['length_score']:.3f}, similarity={best_struct['similarity_score']:.3f}, "
                    f"coverage={best_struct['coverage_score']:.3f})")
        logging.info(f"  Structure length: {int(best_struct['structure_length'])} aa, "
                    f"MSA average: {int(best_struct['msa_avg_length'])} aa")

    return selected_structures

def main():
    """Main function for 3D structure selection"""
    try:
        reference_dir = Path(snakemake.input.reference_dir)
        main_msa_dir =Path(snakemake.input.reference_dir)
        selected_3d_paths = Path(snakemake.output.selected_3d_paths)
        selected_3d_tsv = Path(snakemake.output.selected_3d_tsv)

        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        # Using combined length-similarity-coverage scoring
        min_similarity_score = snakemake.params.min_similarity_score
        max_structures = snakemake.params.max_structures_per_gene

        # Get scoring weights from config
        scoring_weights = snakemake.config.get('structure_selection', {}).get('scoring_weights',
                                             {'length': 0.5, 'similarity': 0.3, 'coverage': 0.2})

        # Get length filtering parameters
        length_config = snakemake.config.get('structure_selection', {}).get('length_filtering', {})
        enable_prefilter = length_config.get('enable_prefilter', False)
        max_length_deviation = length_config.get('max_length_deviation', 0.5)

        # Set random seed for deterministic structure selection
        random_seed = snakemake.config.get('structure_selection', {}).get('random_seed', 42)
        random.seed(random_seed)
        logging.info(f"Random seed set to: {random_seed}")

        # Structure directory path
        structure_base_dir = Path("data/protein_structures")

    except NameError:
        # For standalone testing
        reference_dir = Path(sys.argv[1])
        main_msa_dir = Path(sys.argv[2])
        selected_3d_paths = Path(sys.argv[3])
        selected_3d_tsv = Path(sys.argv[4])

        analysis = "test_analysis"
        paramset = "test_params"
        group = "test_group"
        similarity_method = "local_alignment"
        min_similarity_score = 0.1
        max_structures = 1
        structure_base_dir = Path("data/protein_structures")

    logging.info("=== 3D Structure Selection from Main MSA (STAGE 2) ===")
    logging.info(f"Reference directory: {reference_dir}")
    logging.info(f"Main MSA directory: {main_msa_dir}")
    logging.info(f"Analysis: {analysis}, Paramset: {paramset}, Group: {group}")
    logging.info(f"Selection method: Length-based + Similarity + Coverage (combined scoring)")
    logging.info(f"Scoring weights: Length={scoring_weights['length']:.1%}, "
                f"Similarity={scoring_weights['similarity']:.1%}, "
                f"Coverage={scoring_weights['coverage']:.1%}")
    logging.info(f"Length pre-filtering: {'Enabled' if enable_prefilter else 'Disabled'}")
    if enable_prefilter:
        logging.info(f"Max length deviation: {max_length_deviation:.1%}")
    logging.info(f"Min combined score threshold: {min_similarity_score}")
    logging.info(f"Max structures per gene: {max_structures}")

    # Create output directories
    selected_3d_paths.parent.mkdir(parents=True, exist_ok=True)
    selected_3d_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Find all genes with main MSA files
    main_msa_files = list(main_msa_dir.glob("*.fasta"))
    gene_names = [f.stem for f in main_msa_files]

    logging.info(f"\nProcessing {len(gene_names)} genes for 3D structure selection...")

    all_selected_structures = []
    selected_paths = []

    for i, gene_name in enumerate(gene_names, 1):
        logging.info(f"\n[{i}/{len(gene_names)}] Processing {gene_name}")

        structure_dir = structure_base_dir / gene_name

        selected_structures = select_best_3d_structure(
            gene_name, main_msa_dir, structure_dir,
            "length_similarity_coverage", min_similarity_score, max_structures,
            scoring_weights, enable_prefilter, max_length_deviation
        )

        if selected_structures:
            all_selected_structures.extend(selected_structures)
            for struct in selected_structures:
                selected_paths.append(struct['file_path'])
        else:
            logging.info(f"  No suitable 3D structures found for {gene_name}")

    # Write output files
    with open(selected_3d_paths, 'w') as f:
        for path in selected_paths:
            f.write(f"{path}\n")

    # Write TSV file
    if all_selected_structures:
        df = pd.DataFrame(all_selected_structures)
        # Reorder columns to include new scoring metrics
        column_order = ['gene', 'pdb', 'chain', 'bacteria', 'length_score', 'similarity_score',
                       'coverage_score', 'combined_score', 'consensus_similarity', 'total_score',
                       'structure_length', 'msa_avg_length', 'file_name', 'selection_method']
        # Only include columns that exist in the dataframe
        available_columns = [col for col in column_order if col in df.columns]
        df = df[available_columns]
        df.to_csv(selected_3d_tsv, sep='\t', index=False)
    else:
        # Create empty TSV with headers
        empty_df = pd.DataFrame(columns=['gene', 'pdb', 'chain', 'bacteria', 'length_score',
                                        'similarity_score', 'coverage_score', 'combined_score',
                                        'consensus_similarity', 'total_score', 'structure_length',
                                        'msa_avg_length', 'file_name', 'selection_method'])
        empty_df.to_csv(selected_3d_tsv, sep='\t', index=False)

    # Summary
    logging.info(f"\n=== 3D Structure Selection Summary ===")
    logging.info(f"Genes processed: {len(gene_names)}")
    logging.info(f"Structures selected: {len(all_selected_structures)}")
    logging.info(f"Output paths file: {selected_3d_paths}")
    logging.info(f"Output TSV file: {selected_3d_tsv}")

    if all_selected_structures:
        avg_combined = sum(s['combined_score'] for s in all_selected_structures) / len(all_selected_structures)
        avg_length = sum(s['length_score'] for s in all_selected_structures) / len(all_selected_structures)
        avg_similarity = sum(s['similarity_score'] for s in all_selected_structures) / len(all_selected_structures)
        avg_coverage = sum(s['coverage_score'] for s in all_selected_structures) / len(all_selected_structures)

        logging.info(f"Average combined score: {avg_combined:.3f}")
        logging.info(f"Average length score: {avg_length:.3f}")
        logging.info(f"Average similarity score: {avg_similarity:.3f}")
        logging.info(f"Average coverage score: {avg_coverage:.3f}")

    logging.info(f"\n=== 3D structure selection completed! ===")

if __name__ == "__main__":
    main()
