# Snakefile
# ------------------------------------------------------
# PureMilk Snakemake Workflow
# 
# This pipeline downloads BacDive metadata, classifies species by Gram stain,
# fetches QuickGO annotations, selects important proteins, downloads protein sequences,
# performs MSA (MAFFT), and finds conserved amino acid regions.
#
# Configurations are loaded from config/config.yaml.
# ------------------------------------------------------


# Load configuration
configfile: "config/config.yaml"

# Helper function to dynamically load selected protein names
def get_selected_genes(group):
    file_path = f"results/{group}_proteins.txt"
    with open(file_path) as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes

rule all:
    input:
        expand(
            "results/{group}/fasta_merged/{gene}",
            group=["gram_positive"],
            gene=get_selected_genes("gram_positive")
        )


# ---------------------
# Final Target
# ---------------------
# rule all:
#     input:
#         expand(
#             "results/positive_proteins_aligned/{gene}_aligned.fasta",
#             gene=get_selected_genes("positive")
#         ) +
#         expand(
#             "results/negative_proteins_aligned/{gene}_aligned.fasta",
#             gene=get_selected_genes("negative")
#         )


# ---------------------
# Rules
# ---------------------

# Classify species as Gram-positive or Gram-negative using BacDive API
rule classify_taxa_by_gram:
    input:
        species=config["bacdive"]["species"],
        bacdive_info=config["login"]["bacdive_info"]
    output:
        config["bacdive"]["all_json"],
        config["bacdive"]["not_found"],
        config["bacdive"]["errors"],
        config["bacdive"]["downloaded"],
        config["bacdive"]["gram_classification"]
    script:
        config["scripts"]["classify_gram"]

# Split species list into Gram-positive and Gram-negative text files
rule split_species_by_gram:
    input:
        config["bacdive"]["gram_classification"]
    output:
        gram_positive=config["bacdive"]["gram_positive"],
        gram_negative=config["bacdive"]["gram_negative"]
    shell:
        """
        mkdir -p data/bacdive
        awk -F'\t' 'NR > 1 && $2 == "positive" {{ print $1 }}' {input} > {output.gram_positive}
        awk -F'\t' 'NR > 1 && $2 == "negative" {{ print $1 }}' {input} > {output.gram_negative}
        """

# Fetch GO annotations and gene symbols from QuickGO
# There is something wrong with the QuickGO API, so we need to use a local file for now
# rule fetch_quickgo_annotations:
#     input:
#         go_ids = config["quickgo"]["go_ids"],
#         taxon_ids = config["quickgo"]["taxon_ids"]
#     output:
#         annotations = config["quickgo"]["annotations"],
#         genes = config["quickgo"]["genes"]
#     script:
#         config["scripts"]["fetch_quickgo"]

# Assess how many species have each gene annotated
rule assess_gene_taxa_coverage:
    input:
        species=lambda wildcards: config["bacdive"]["gram_species_template"].format(group=wildcards.group),
        genes=config["quickgo"]["genes"],
        ncbi_info=config["login"]["ncbi_info"]
    output:
        coverage="results/gene_coverage_gram_{group}.tsv"
    script:
        config["scripts"]["gene_coverage"]

# Filter genes by minimum species coverage and sort them
rule filter_and_sort_coverage:
    input:
        "results/gene_coverage_gram_{group}.tsv"
    output:
        "results/gene_coverage_gram_{group}_filtered.tsv"
    params:
        threshold=lambda wildcards: config["gram_thresholds"][wildcards.group]
    shell:
        """
        (head -n 1 {input} && awk -F'\t' 'NR > 1 && $2 > {params.threshold}' {input} | sort -k2,2nr) > {output}
        """

# Select top proteins to analyse further (after coverage filtering)
rule select_proteins_to_analyse:
    input:
        "results/gene_coverage_gram_{group}_filtered.tsv"
    output:
        "results/proteins_unique_in_{group}.txt"
    params:
        num_proteins_gram_positive=config["protein_selection"]["positive"],
        num_proteins_gram_negative=config["protein_selection"]["negative"]
    script:
        config["scripts"]["select_proteins"]

# Download protein sequences from NCBI for selected genes and species
rule download_proteins_to_analyse:
    input:
        proteins="results/{group}_proteins.txt",
        species=lambda wildcards: config["bacdive"]["gram_species_template"].format(group=wildcards.group),
        ncbi_info=config["login"]["ncbi_info"]
    output:
        output_folder=directory("results/{group}")
    params:
        group="{group}"
    script:
        config["scripts"]["download_proteins"]

# Select one representative protein per species for MSA
rule select_proteins_for_msa:
    input:
        input_folder="results/{group}/{gene}"
    output:
        output_folder="results/{group}/fasta_merged/{gene}.fasta",
    params:
        group="{group}"
    script:
        config["scripts"]["get_msa_sequences"]

# Perform multiple sequence alignment (MSA) using MAFFT with threading
rule run_mafft:
    input:
        fasta="results/{group}/fasta_merged/{gene}.fasta"
    output:
        msa="results/{group}/fasta_msa/{gene}.fasta"
    threads: config["mafft"]["threads"]
    shell:
        """
        mkdir -p $(dirname {output.msa})
        mafft --thread {threads} --maxiterate 10000 --localpair {input.fasta} > {output.msa}
        """

# MSA quality check using AliStat
rule check_alignment_quality:
    input:
        msa="results/{group}/fasta_msa/{gene}.fasta"
    output:
        quality="results/{group}/fasta_msa/alistat/{gene}_aligned.stat"
    log:
        "logs/{group}/fasta_msa_alistat/{gene}_alistat.log"
    shell:
        """
        mkdir -p $(dirname {output.quality})
        mkdir -p $(dirname {log})
        alistat {input.msa} > {output.quality} 2> {log}
        """

# Trim bad regions with trimAl
rule trim_alignment:
    input:
        msa="results/{group}/fasta_msa/{gene}.fasta"
    output:
        trimmed="results/{group}/fasta_msa_trimmed/{gene}.fasta"
    shell:
        """
        mkdir -p $(dirname {output.trimmed})
        trimal -in {input.msa} -out {output.trimmed} -automated1
        """


# Find conserved amino acid positions in the MSA
# rule find_conserved_amino_acids:
#     input:
#         expand(
#             "results/positive_proteins_trimmed/{gene}_trimmed.fasta",
#             gene=get_selected_genes("positive")
#         ) +
#         expand(
#             "results/negative_proteins_trimmed/{gene}_trimmed.fasta",
#             gene=get_selected_genes("negative")
#         )
#     script:
#         config["scripts"]["conserved_aminoacids"]


