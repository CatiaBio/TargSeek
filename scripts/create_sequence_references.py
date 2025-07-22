#!/usr/bin/env python3
"""
Script to generate only the *_filelist.txt for one gene
Given a gene list directory and a specific gene, it finds valid FASTA sequences
and writes the list of paths to {gene}_filelist.txt under the output directory.
"""

import sys
import logging
from pathlib import Path
from select_msa_proteins import (
    load_gene_species_lists,
    find_available_sequences,
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    try:
        gene_lists_dir = Path(snakemake.input.gene_lists)
        output_dir = Path(snakemake.output[0])
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
    except NameError:
        gene_lists_dir = Path(sys.argv[1])
        output_dir = Path(sys.argv[2])
        analysis = "analysis_1"
        paramset = "params_1"
        group = "gram_negative"

    output_dir.mkdir(parents=True, exist_ok=True)

    gene_species_mapping = load_gene_species_lists(gene_lists_dir)

    for gene_name, species_list in gene_species_mapping.items():
        logging.info(f"Generating filelist for gene: {gene_name}")
        sequences = find_available_sequences(gene_name, species_list)

        gene_dir = output_dir / gene_name
        gene_dir.mkdir(parents=True, exist_ok=True)

        filelist_path = gene_dir / f"{gene_name}_filelist.txt"
        with open(filelist_path, 'w') as f:
            for species, fasta_path in sequences.items():
                f.write(f"{fasta_path}\n")

        logging.info(f"Wrote {len(sequences)} entries to {filelist_path}")

if __name__ == "__main__":
    main()
