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

# Define wildcard constraints to prevent ambiguous matches
wildcard_constraints:
    analysis="analysis_[0-9]+",
    paramset="params_[0-9]+"

rule all_coverage_filtered:
    input:
        expand(
        "results/coverage/{analysis}_{paramset}_gram_{group}_coverage_count.tsv",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_uniprot_info:
    input:
        expand(
        "results/uniprot_info/{analysis}_{paramset}_gram_{group}_uniprot_info",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_bacteria_protein_info:
    input:
        expand(
        "results/protein_info/{analysis}_{paramset}_gram_{group}_protein_info_bacteria.tsv",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

# Removed: all_selected_proteins rule - protein selection step has been removed

rule all_surface_accessible_proteins:
    input:
        expand(
        "results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_downloaded_proteins:
    input:
        expand(
        "results/proteins/{analysis}_{paramset}_gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_sequences:
    input:
        expand(
        "results/msa_sequences/{analysis}_{paramset}_gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_alignments:
    input:
        expand(
        "results/msa_alignments/{analysis}_{paramset}_gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_quality:
    input:
        expand(
        "results/msa_quality/{analysis}_{paramset}_gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

# ---------------------
# Rules
# ---------------------

# Classify species as Gram-positive or Gram-negative using BacDive API
rule classify_taxa_by_gram:
    input:
        species=lambda wildcards: config["species_files"][wildcards.analysis],
        bacdive_info=config["login"]["bacdive_info"]
    output:
        all_json="data/bacdive/{analysis}/gram_raw.json",
        not_found="data/bacdive/{analysis}/not_found.txt",
        errors="data/bacdive/{analysis}/errors.txt",
        downloaded="data/bacdive/{analysis}/downloaded.txt",
        gram_classification="data/bacdive/{analysis}/gram.tsv"
    script:
        "scripts_test/classify_gram.py"

# Supplement BacDive Gram classification using genus-based inference
rule supplement_bacdive_gram_classification:
    input:
        bacdive_classification = "data/bacdive/{analysis}/gram.tsv",
        not_found = "data/bacdive/{analysis}/not_found.txt",
        bacdive_info = config["login"]["bacdive_info"]
    output:
        updated_classification = "data/bacdive/{analysis}/updated_gram.tsv",
        updated_not_found = "data/bacdive/{analysis}/updated_not_found.txt"
    script:
        "scripts_test/supplement_gram_classification.py"


# Split species list into Gram-positive and Gram-negative text files 
rule split_species_by_gram:
    input:
        gram_classification = "data/bacdive/{analysis}/updated_gram.tsv"
    output:
        gram_positive = "data/bacdive/{analysis}/gram_positive.txt",
        gram_negative = "data/bacdive/{analysis}/gram_negative.txt"
    shell:
        """
        awk -F'\t' 'NR > 1 && $2 == "positive" {{ print $1 }}' {input.gram_classification} > {output.gram_positive}
        awk -F'\t' 'NR > 1 && $2 == "negative" {{ print $1 }}' {input.gram_classification} > {output.gram_negative}
        """
     
# Fetch GO annotations and gene symbols from QuickGO
rule fetch_quickgo_annotations:
    input:
        params_file = "config/quickgo/{paramset}.json"
    output:
        annotations = "data/quickgo/{paramset}/annotations.tsv",
        genes = "data/quickgo/{paramset}/gene_symbols.txt"
    script:
        "scripts_test/fetch_quickgo_data.py"


# Match genes to taxon names from QuickGO annotations
rule filter_quickgo_genes:
    input:
        genes = "data/quickgo/{paramset}/gene_symbols.txt",
        annotations = "data/quickgo/{paramset}/annotations.tsv"
    output:
        filtered_genes = "data/quickgo/{paramset}/gene_symbols_filtered.txt"
    script:
        "scripts_test/filter_quickgo_data.py"

# Note: Cache initialization is now handled within assess_gene_taxa_coverage rule

# Assess how many species have each gene annotated
rule assess_gene_taxa_coverage_positive:
    input:
        species_list="data/bacdive/{analysis}/gram_positive.txt",
        gene_list="data/quickgo/{paramset}/gene_symbols_filtered.txt",
        ncbi_info=config["login"]["ncbi_info"]
    output:
        coverage="results/coverage/{analysis}_{paramset}_gram_positive_coverage.tsv"
    priority: 10
    log:
        "logs/coverage/{analysis}_{paramset}_gram_positive_coverage.log"
    script:
        "scripts/gene_taxa_coverage_cached.py"

rule assess_gene_taxa_coverage_negative:
    input:
        species_list="data/bacdive/{analysis}/gram_negative.txt",
        gene_list="data/quickgo/{paramset}/gene_symbols_filtered.txt",
        ncbi_info=config["login"]["ncbi_info"]
    output:
        coverage="results/coverage/{analysis}_{paramset}_gram_negative_coverage.tsv"
    priority: 5
    log:
        "logs/coverage/{analysis}_{paramset}_gram_negative_coverage.log"
    script:
        "scripts/gene_taxa_coverage_cached.py" 

# Filter genes by minimum species coverage and create count summary
rule filter_and_sort_coverage:
    input:
        "results/coverage/{analysis}_{paramset}_gram_{group}_coverage.tsv"
    output:
        "results/coverage/{analysis}_{paramset}_gram_{group}_coverage_count.tsv"
    params:
        threshold=lambda wildcards: config["gram_thresholds"][wildcards.group]
    script:
        config["scripts"]["filter_coverage"]

# Fetch comprehensive protein information from UniProt (bacteria-specific)
rule fetch_uniprot_info:
    input:
        "results/coverage/{analysis}_{paramset}_gram_{group}_coverage_count.tsv"
    output:
        directory("results/uniprot_info/{analysis}_{paramset}_gram_{group}_uniprot_info")
    params:
        max_species_per_gene=3  # Limit to 3 species per gene for faster processing
    script:
        "scripts/fetch_uniprot_info.py"

# Filter proteins by surface accessibility based on GO cellular component
rule filter_surface_accessible_proteins:
    input:
        coverage="results/uniprot_info/{analysis}_{paramset}_gram_{group}_uniprot_info/{analysis}_{paramset}_gram_{group}_coverage_count_location.tsv",
        surface_terms="config/quickgo/surface_accessible.txt"
    output:
        "results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv"
    script:
        "scripts/filter_surface_accessible_proteins.py"

# Download protein sequences from NCBI for surface-accessible proteins
rule download_proteins_to_analyse:
    input:
        proteins="results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv",
        species="data/bacdive/{analysis}/gram_{group}.txt",
        ncbi_info=config["login"]["ncbi_info"]
    output:
        output_folder=directory("results/proteins/{analysis}_{paramset}_gram_{group}"),
        complete_flag="results/proteins/{analysis}_{paramset}_gram_{group}/.download_complete"
    params:
        group="{group}",
        analysis="{analysis}",
        paramset="{paramset}"
    script:
        config["scripts"]["download_proteins"]

# Select one representative protein per species for MSA
rule select_proteins_for_msa:
    input:
        protein_dir="results/proteins/{analysis}_{paramset}_gram_{group}",
        protein_list="results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv",
        download_complete="results/proteins/{analysis}_{paramset}_gram_{group}/.download_complete"
    output:
        directory("results/msa_sequences/{analysis}_{paramset}_gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}", 
        group="{group}"
    script:
        config["scripts"]["get_msa_sequences"]

# Perform multiple sequence alignment (MSA) using MAFFT with threading
rule run_mafft:
    input:
        fasta="results/msa_sequences/{analysis}_{paramset}_gram_{group}/{gene}.fasta"
    output:
        msa="results/msa_alignments/{analysis}_{paramset}_gram_{group}/{gene}_aligned.fasta"
    threads: config["mafft"]["threads"]
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        gene="{gene}"
    shell:
        """
        mkdir -p $(dirname {output.msa})
        mafft --thread {threads} --auto {input.fasta} > {output.msa}
        """

# Run MAFFT alignments for all genes in a group
rule run_all_mafft_for_group:
    input:
        msa_dir="results/msa_sequences/{analysis}_{paramset}_gram_{group}",
        protein_list="results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv"
    output:
        directory("results/msa_alignments/{analysis}_{paramset}_gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        threads=config["mafft"]["threads"]
    script:
        "scripts/run_mafft_alignments.py"

# MSA quality assessment
rule assess_alignment_quality:
    input:
        "results/msa_alignments/{analysis}_{paramset}_gram_{group}"
    output:
        directory("results/msa_quality/{analysis}_{paramset}_gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/assess_alignment_quality.py"

# # MSA quality check using AliStat
# rule check_alignment_quality:
#     input:
#         msa="results/{group}/fasta_msa/{gene}.fasta"
#     output:
#         quality="results/{group}/fasta_msa/alistat/{gene}_aligned.stat"
#     log:
#         "logs/{group}/fasta_msa_alistat/{gene}_alistat.log"
#     shell:
#         """
#         mkdir -p $(dirname {output.quality})
#         mkdir -p $(dirname {log})
#         alistat {input.msa} > {output.quality} 2> {log}
#         """

# # Trim bad regions with trimAl
# rule trim_alignment:
#     input:
#         msa="results/{group}/fasta_msa/{gene}.fasta"
#     output:
#         trimmed="results/{group}/fasta_msa_trimmed/{gene}.fasta"
#     shell:
#         """
#         mkdir -p $(dirname {output.trimmed})
#         trimal -in {input.msa} -out {output.trimmed} -automated1
#         """


# # Find conserved amino acid positions in the MSA
# # rule find_conserved_amino_acids:
# #     input:
# #         expand(
# #             "results/positive_proteins_trimmed/{gene}_trimmed.fasta",
# #             gene=get_selected_genes("positive")
# #         ) +
# #         expand(
# #             "results/negative_proteins_trimmed/{gene}_trimmed.fasta",
# #             gene=get_selected_genes("negative")
# #         )
# #     script:
# #         config["scripts"]["conserved_aminoacids"]


