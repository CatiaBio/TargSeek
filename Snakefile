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

rule all_coverage_unified:
    input:
        expand(
        "results/{analysis}_{paramset}/coverage",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"]
        )

rule all_epitope_predictions_bepipred:
    input:
        expand(
        "results/{analysis}_{paramset}/epitope_predictions_bepipred/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_uniprot_info:
    input:
        expand(
        "results/{analysis}_{paramset}/uniprot_info/gram_{group}_uniprot_info",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_bacteria_protein_info:
    input:
        expand(
        "results/{analysis}_{paramset}/protein_info/gram_{group}_protein_info_bacteria.tsv",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

# Removed: all_selected_proteins rule - protein selection step has been removed

rule all_surface_accessible_proteins:
    input:
        expand(
        "results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_gene_species_lists:
    input:
        expand(
        "data/{analysis}_{paramset}/genes_species/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_downloaded_proteins:
    input:
        expand(
        "data/proteins_fasta/.{analysis}_{paramset}_{group}_download_complete",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_sequences:
    input:
        expand(
        "results/{analysis}_{paramset}/msa_sequences/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_alignments:
    input:
        expand(
        "results/{analysis}_{paramset}/msa_alignments/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_trimmed:
    input:
        expand(
        "results/{analysis}_{paramset}/msa_trimmed/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_quality:
    input:
        expand(
        "results/{analysis}_{paramset}/msa_quality/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_conservation:
    input:
        expand(
        "results/{analysis}_{paramset}/conservation/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_download_summaries:
    input:
        expand(
        "results/{analysis}_{paramset}/download_summary/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_epitope_predictions:
    input:
        expand(
        "results/{analysis}_{paramset}/epitope_predictions/gram_{group}",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_reports:
    input:
        expand(
        "results/{analysis}_{paramset}/reports/final_report.html",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"]
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
        gram_classification="data/bacdive/{analysis}/gram.tsv",
        all_identified="data/bacdive/{analysis}/all_identified.txt"
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
        annotations = "data/quickgo/{paramset}/annotations.json",
        genes = "data/quickgo/{paramset}/gene_symbols.txt"
    script:
        "scripts_test/fetch_quickgo_data.py"


# Fetch gene aliases/synonyms from NCBI before filtering
rule fetch_gene_aliases:
    input:
        genes = "data/quickgo/{paramset}/gene_symbols.txt"
    output:
        aliases_file = "data/quickgo/{paramset}/gene_aliases.txt"
    script:
        "scripts/fetch_gene_aliases.py"

# Validate GO term assignments for proteins with alias support (moved from later in pipeline)
rule validate_go_assignments:
    input:
        aliases = "data/quickgo/{paramset}/gene_aliases.txt"
    output:
        validation_report = "data/quickgo/{paramset}/protein_go_validation_report.tsv" 
    script:
        "scripts/validate_protein_go_assignments.py"

# Filter genes by surface accessibility using GO validation
rule filter_surface_accessible_genes:
    input:
        go_validation = "data/quickgo/{paramset}/protein_go_validation_report.tsv",
        surface_accessible = "config/quickgo/surface_accessible.txt"
    output:
        surface_genes = "data/quickgo/{paramset}/surface_accessible_proteins.txt"
    script:
        "scripts_test/filter_surface_accessible.py"

# Filter genes by taxa coverage using surface accessible genes
rule filter_quickgo_genes:
    input:
        genes = "data/quickgo/{paramset}/surface_accessible_proteins.txt",
        annotations = "data/quickgo/{paramset}/annotations.json",
        aliases = "data/quickgo/{paramset}/gene_aliases.txt"
    output:
        proteins_to_test = "data/quickgo/{paramset}/proteins_to_be_tested.txt"
    script:
        "scripts_test/filter_quickgo_data_simplified.py"

# Note: Cache initialization is now handled within assess_gene_taxa_coverage rule

# Assess coverage for all species (both Gram-positive and Gram-negative) in unified analysis
rule assess_gene_taxa_coverage_unified:
    input:
        all_species="data/bacdive/{analysis}/all_identified.txt",
        gram_positive="data/bacdive/{analysis}/gram_positive.txt",
        gram_negative="data/bacdive/{analysis}/gram_negative.txt",
        gene_list="data/quickgo/{paramset}/proteins_to_be_tested.txt",
        aliases="data/quickgo/{paramset}/gene_aliases.txt",
        ncbi_info=config["login"]["ncbi_info"]
    output:
        directory("results/{analysis}_{paramset}/coverage")
    priority: 10
    log:
        "logs/{analysis}_{paramset}_unified_coverage.log"
    params:
        analysis="{analysis}",
        paramset="{paramset}"
    script:
        "scripts/gene_taxa_coverage_unified.py" 

# Note: proteins_to_study file is now created by select_proteins_to_study rule
# using coverage data and the pre-filtered proteins_to_be_tested.txt list

# Note: filter_and_sort_coverage rule is now handled by assess_gene_taxa_coverage_unified
# which directly creates the unified coverage_count.tsv file with filtering and sorting

# Select proteins to study from unified coverage data 
rule select_proteins_to_study:
    input:
        coverage="results/{analysis}_{paramset}/coverage/coverage_count.tsv",
        tested_proteins="data/quickgo/{paramset}/proteins_to_be_tested.txt"
    output:
        "results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv"
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        max_proteins=lambda wildcards: config["protein_selection"][wildcards.group]
    conda:
        "env.yml"
    script:
        "scripts/select_proteins_to_study.py"

# Create gene-specific species lists from unified coverage data (species that actually have each gene)
# Now saves to analysis-specific data directory: data/{analysis}_{paramset}/genes_species/{group}/
rule create_gene_species_lists:
    input:
        coverage="results/{analysis}_{paramset}/coverage/coverage_count.tsv",
        proteins="results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv"
    output:
        directory("data/{analysis}_{paramset}/genes_species/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/create_gene_species_lists_from_coverage.py"

# Download proteins with multi-stage approach: UniProt batch -> UniProt individual -> NCBI
# Now downloads to shared data/proteins_fasta/ directory with caching and alias fallback
rule download_proteins_to_analyse:
    input:
        protein_lists="data/{analysis}_{paramset}/genes_species/gram_{group}",
        aliases="data/quickgo/{paramset}/gene_aliases.txt",
        ncbi_info=config["login"]["ncbi_info"]
    output:
        sentinel=touch("data/proteins_fasta/.{analysis}_{paramset}_{group}_download_complete")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/download_proteins_cached_shared.py"

# Download 3D structures and integrate into shared protein directories
# Now downloads to shared data/proteins_3d_structure/ directory with caching
rule download_3d_structures:
    input:
        protein_list="results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv"
    output:
        sentinel=touch("data/proteins_3d_structure/.{analysis}_{paramset}_{group}_structures_complete")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/download_3d_structures_cached_shared.py"

# Select one representative protein per species for MSA
# Now uses shared protein data directories and analysis-specific gene lists
rule select_proteins_for_msa:
    input:
        gene_lists="data/{analysis}_{paramset}/genes_species/gram_{group}",
        protein_list="results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv",
        protein_download_sentinel="data/proteins_fasta/.{analysis}_{paramset}_{group}_download_complete",
        structures_download_sentinel="data/proteins_3d_structure/.{analysis}_{paramset}_{group}_structures_complete"
    output:
        directory("results/{analysis}_{paramset}/msa_sequences/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}", 
        group="{group}"
    script:
        "scripts/select_proteins_for_msa_shared.py"

# Run MAFFT alignments for all genes in a group
rule run_all_mafft_for_group:
    input:
        msa_dir="results/{analysis}_{paramset}/msa_sequences/gram_{group}",
        protein_list="results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv"
    output:
        directory("results/{analysis}_{paramset}/msa_alignments/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        threads=config["mafft"]["threads"]
    script:
        "scripts/run_mafft_alignments.py"

# Trim poorly aligned regions with trimAl
rule trim_alignments:
    input:
        "results/{analysis}_{paramset}/msa_alignments/gram_{group}"
    output:
        directory("results/{analysis}_{paramset}/msa_trimmed/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/trim_alignments.py"

# MSA quality assessment (compare before/after trimming)
rule assess_alignment_quality:
    input:
        raw_alignments="results/{analysis}_{paramset}/msa_alignments/gram_{group}",
        trimmed_alignments="results/{analysis}_{paramset}/msa_trimmed/gram_{group}"
    output:
        directory("results/{analysis}_{paramset}/msa_quality/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/assess_alignment_quality_comparison.py"

# Analyze amino acid conservation using best alignments (raw vs trimmed based on quality)
rule analyze_conservation:
    input:
        raw_alignments="results/msa_alignments/{analysis}_{paramset}_gram_{group}",
        trimmed_alignments="results/msa_trimmed/{analysis}_{paramset}_gram_{group}",
        quality_assessment="results/msa_quality/{analysis}_{paramset}_gram_{group}"
    output:
        directory("results/conservation/{analysis}_{paramset}_gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        create_logos=False  # Set to True to enable logo plots
    script:
        "scripts/analyze_conservation_adaptive.py"

# Predict epitopes using IEDB API for conserved sequences with 3D structure data
rule predict_epitopes:
    input:
        msa_sequences="results/msa_sequences/{analysis}_{paramset}_gram_{group}",
        conservation="results/conservation/{analysis}_{paramset}_gram_{group}",
        structures_3d="results/3d_structures/{analysis}_{paramset}_gram_{group}",
        proteins_to_study="results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv"
    output:
        directory("results/epitope_predictions/{analysis}_{paramset}_gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/predict_epitopes.py"

# Predict epitopes using BepiPred 3.0 for B-cell epitope prediction
rule predict_epitopes_bepipred:
    input:
        msa_sequences="results/msa_sequences/{analysis}_{paramset}_gram_{group}",
        conservation="results/conservation/{analysis}_{paramset}_gram_{group}",
        structures_3d="results/3d_structures/{analysis}_{paramset}_gram_{group}",
        proteins_to_study="results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv"
    output:
        directory("results/epitope_predictions_bepipred/{analysis}_{paramset}_gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/predict_epitopes_bepipred.py"

# Generate download summary with actual vs expected species counts
rule generate_download_summary:
    input:
        coverage="results/uniprot_info/{analysis}_{paramset}_gram_{group}_uniprot_info/{analysis}_{paramset}_gram_{group}_coverage_count_location.tsv",
        protein_fasta="results/protein_fasta/{analysis}_{paramset}_gram_{group}",
        proteins_to_study="results/proteins_to_study/{analysis}_{paramset}_gram_{group}.tsv"
    output:
        directory("results/download_summary/{analysis}_{paramset}_gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/generate_download_summary.py"

# Generate comprehensive final report
rule generate_final_report:
    input:
        quality_positive="results/msa_quality/{analysis}_{paramset}_gram_positive",
        quality_negative="results/msa_quality/{analysis}_{paramset}_gram_negative",
        conservation_positive="results/conservation/{analysis}_{paramset}_gram_positive",
        conservation_negative="results/conservation/{analysis}_{paramset}_gram_negative",
        download_summary_positive="results/download_summary/{analysis}_{paramset}_gram_positive",
        download_summary_negative="results/download_summary/{analysis}_{paramset}_gram_negative"
    output:
        "results/reports/{analysis}_{paramset}_final_report.html"
    params:
        analysis="{analysis}",
        paramset="{paramset}"
    script:
        "scripts/generate_final_report.py"