#!/usr/bin/env python3
"""
New MSA FASTA creation script
- Loads gene FASTA path lists from: results/{analysis}_{paramset}/msa_sequence_refs/gram_{group}/{gene}/{gene}_filelist.txt
- For each species FASTA, selects 1 sequence (closest to average length)
- Outputs:
  - msa_alignments/gram_{group}/{gene}.fasta (no 3D structure)
  - msa_alignments_msa_with_3d_fasta/gram_{group}/{gene}.fasta (with 3D structure)
"""

import sys
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

COMMON_BACTERIA = [
    'Escherichia coli',
    'Salmonella enterica',
    'Bacillus subtilis',
    'Staphylococcus aureus',
    'Pseudomonas aeruginosa',
    'Mycobacterium tuberculosis'
]

def select_representative_sequence(fasta_path: Path) -> SeqRecord:
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    if not sequences:
        raise ValueError(f"No sequences in {fasta_path}")
    avg_len = sum(len(seq.seq) for seq in sequences) / len(sequences)
    return min(sequences, key=lambda rec: abs(len(rec.seq) - avg_len))

def select_3d_structure_sequence(gene_name: str) -> SeqRecord:
    structure_dir = Path("data/proteins_3d_structure") / gene_name
    fasta_files = sorted(structure_dir.glob("*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(f"No structure FASTA files for {gene_name}")

    def priority(f):
        header = next(SeqIO.parse(f, "fasta")).description.lower()
        for i, species in enumerate(COMMON_BACTERIA):
            if species.lower() in header:
                return i
        return len(COMMON_BACTERIA)

    sorted_files = sorted(fasta_files, key=priority)
    return next(SeqIO.parse(sorted_files[0], "fasta"))

def load_filelist(filelist_path: Path) -> list:
    if not filelist_path.exists():
        logging.warning(f"Missing: {filelist_path}")
        return []
    return [Path(line.strip()) for line in filelist_path.read_text().splitlines() if line.strip()]

def create_gene_msa(gene_dir: Path, output_dir_msa_no_3d: Path, output_dir_msa_with_3d: Path):
    gene_name = gene_dir.name
    filelist_path = gene_dir / f"{gene_name}_filelist.txt"
    fasta_paths = load_filelist(filelist_path)

    selected_records = []
    for fasta_path in fasta_paths:
        try:
            record = select_representative_sequence(fasta_path)
            record.id = fasta_path.stem
            record.description += " [SEQ]"
            selected_records.append(record)
        except Exception as e:
            logging.warning(f"Failed to select from {fasta_path}: {e}")

    # Write no-3D MSA file
    if selected_records:
        msa_no_3d_out = output_dir_msa_no_3d / f"{gene_name}.fasta"
        output_dir_msa_no_3d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(selected_records, msa_no_3d_out, "fasta")
        logging.info(f"Written no-3D MSA: {msa_no_3d_out}")

    # Add 3D structure
    try:
        structure_record = select_3d_structure_sequence(gene_name)
        structure_record.id = f"3D_{structure_record.id}"
        structure_record.description += " [3D-STRUCTURE]"
        msa_with_3d_records = selected_records + [structure_record]
        msa_with_3d_out = output_dir_msa_with_3d / f"{gene_name}.fasta"
        output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(msa_with_3d_records, msa_with_3d_out, "fasta")
        logging.info(f"Written with-3D MSA: {msa_with_3d_out}")
    except Exception as e:
        logging.warning(f"No 3D structure added for {gene_name}: {e}")


def main():
    try:
        input_dir = Path(snakemake.input.reference_dir)
        output_dir_msa_no_3d = Path(snakemake.output.msa_no_3d)
        output_dir_msa_with_3d = Path(snakemake.output.msa_with_3d)
    except NameError:
        input_dir = Path(sys.argv[1])
        output_dir_msa_no_3d = Path(sys.argv[2])
        output_dir_msa_with_3d = Path(sys.argv[3])

    output_dir_msa_no_3d.mkdir(parents=True, exist_ok=True)
    output_dir_msa_with_3d.mkdir(parents=True, exist_ok=True)

    for gene_dir in sorted(input_dir.iterdir()):
        if not gene_dir.is_dir():
            continue
        create_gene_msa(gene_dir, output_dir_msa_no_3d, output_dir_msa_with_3d)

if __name__ == "__main__":
    main()
