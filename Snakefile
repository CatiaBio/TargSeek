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

# Helper function to get alignment directory based on config
def get_alignment_dir(wildcards):
    """Get the appropriate alignment directory based on config setting"""
    use_3d = config.get("mafft", {}).get("use_3d_alignments", "no_3d")
    if use_3d == "with_3d":
        return f"results/{wildcards.analysis}_{wildcards.paramset}/msa_alignments_with_3d_fasta/gram_{wildcards.group}"
    else:
        return f"results/{wildcards.analysis}_{wildcards.paramset}/msa_alignments/gram_{wildcards.group}"

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

# Download all bacterial 3D structures for all analyses and groups
# Convenience rule to download 3D structures for both gram-positive and gram-negative
# Structures are stored in shared data/proteins_3d_structure/ directory for reuse
rule all_3d_structures:
    input:
        expand(
        "data/proteins_3d_structure/.{analysis}_{paramset}_{group}_structures_complete",
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
        ["results/{analysis}_{paramset}/no_3d/msa_sequences/gram_{group}",
         "results/{analysis}_{paramset}/with_3d/msa_sequences/gram_{group}"],
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_alignments:
    input:
        expand(
        ["results/{analysis}_{paramset}/no_3d/msa_alignments/gram_{group}",
         "results/{analysis}_{paramset}/with_3d/msa_alignments/gram_{group}"],
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_trimmed:
    input:
        expand(
        ["results/{analysis}_{paramset}/no_3d/msa_trimmed/gram_{group}",
         "results/{analysis}_{paramset}/with_3d/msa_trimmed/gram_{group}"],
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_msa_quality:
    input:
        expand(
        ["results/{analysis}_{paramset}/no_3d/msa_quality/gram_{group}",
         "results/{analysis}_{paramset}/with_3d/msa_quality/gram_{group}"],
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          group=["positive", "negative"]
        )

rule all_conservation:
    input:
        expand(
        ["results/{analysis}_{paramset}/no_3d/conservation/gram_{group}",
         "results/{analysis}_{paramset}/with_3d/conservation/gram_{group}"],
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

# REMOVED: all_epitope_predictions rule - use all_epitope_predictions_bepipred instead

rule all_reports:
    input:
        expand(
        "results/{analysis}_{paramset}/reports/final_report_{use_3d_dir}.html",
            analysis=config["species_batches"],
          paramset=config["quickgo_paramsets"],
          use_3d_dir=["no_3d", "with_3d"]
        )

# ---------------------
# Rules
# ---------------------

# Classify bacterial species by Gram staining properties
# Queries BacDive API for comprehensive Gram stain classification data
# Handles species identification and taxonomic validation with error tracking
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
        "scripts/classify_gram.py"

# Enhance Gram classification using genus-based inference
# Supplements BacDive results by inferring Gram status from well-characterized genus patterns
# Reduces unclassified species by leveraging taxonomic relationships
rule supplement_bacdive_gram_classification:
    input:
        bacdive_classification = "data/bacdive/{analysis}/gram.tsv",
        not_found = "data/bacdive/{analysis}/not_found.txt",
        bacdive_info = config["login"]["bacdive_info"]
    output:
        updated_classification = "data/bacdive/{analysis}/updated_gram.tsv",
        updated_not_found = "data/bacdive/{analysis}/updated_not_found.txt"
    script:
        "scripts/supplement_gram_classification.py"


# Separate species by Gram classification for parallel processing
# Creates Gram-specific species lists from unified classification data
# Enables independent analysis of Gram-positive and Gram-negative bacterial groups
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
     
# Retrieve Gene Ontology annotations and symbols from QuickGO database
# Downloads comprehensive GO term assignments and associated gene nomenclature
# Filters by specified taxonomic groups and GO categories for targeted protein discovery
rule fetch_quickgo_annotations:
    input:
        params_file = "config/quickgo/{paramset}.json"
    output:
        annotations = "data/quickgo/{paramset}/annotations.json",
        genes = "data/quickgo/{paramset}/gene_symbols.txt"
    script:
        "scripts/fetch_quickgo_data.py"


# Collect gene aliases and synonyms from NCBI Gene database
# Builds comprehensive alias mapping to improve protein search success rates
# Essential for handling gene nomenclature variations across different databases
rule fetch_gene_aliases:
    input:
        genes = "data/quickgo/{paramset}/gene_symbols.txt"
    output:
        aliases_file = "data/quickgo/{paramset}/gene_aliases.txt"
    script:
        "scripts/fetch_gene_aliases.py"

# Validate Gene Ontology term assignments using multiple gene identifiers
# Cross-references GO annotations with gene aliases to ensure assignment accuracy
# Generates validation reports for quality control of GO-based protein selection
rule validate_go_assignments:
    input:
        aliases = "data/quickgo/{paramset}/gene_aliases.txt"
    output:
        validation_report = "data/quickgo/{paramset}/protein_go_validation_report.tsv" 
    script:
        "scripts/validate_protein_go_assignments.py"

# Select surface-accessible proteins using GO term validation
# Filters candidates by cellular localization GO terms (membrane, extracellular, etc.)
# Prioritizes proteins likely to be accessible for diagnostic or therapeutic targeting
rule filter_surface_accessible_genes:
    input:
        go_validation = "data/quickgo/{paramset}/protein_go_validation_report.tsv",
        surface_accessible = "config/quickgo/surface_accessible.txt"
    output:
        surface_genes = "data/quickgo/{paramset}/surface_accessible_proteins.txt"
    script:
        "scripts/filter_surface_accessible.py"

# Apply taxonomic coverage filters to surface-accessible proteins
# Removes genes with insufficient representation across target bacterial species
# Ensures selected proteins have broad applicability within the study group
rule filter_quickgo_genes:
    input:
        genes = "data/quickgo/{paramset}/surface_accessible_proteins.txt",
        annotations = "data/quickgo/{paramset}/annotations.json",
        aliases = "data/quickgo/{paramset}/gene_aliases.txt"
    output:
        proteins_to_test = "data/quickgo/{paramset}/proteins_to_be_tested.txt"
    script:
        "scripts/filter_quickgo_data_simplified.py"

# Assess gene coverage across all bacterial species in unified analysis
# Systematically queries NCBI Protein database to determine gene presence/absence
# Uses intelligent caching and alias fallback to maximize coverage detection accuracy
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

# Select proteins to study based on coverage thresholds
# Filters candidate proteins by Gram-specific coverage requirements (50% for both groups)
# Prioritizes proteins with highest taxonomic coverage within each Gram classification
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
    group: "gram_group_{analysis}_{paramset}"
    conda:
        "env.yml"
    script:
        "scripts/select_proteins_to_study.py"

# Create gene-specific species lists for protein downloads
# Extracts species that actually contain each selected gene from coverage data
# Organizes species lists by gene in analysis-specific directories for downstream processing
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
    group: "gram_group_{analysis}_{paramset}"
    script:
        "scripts/create_gene_species_lists_from_coverage.py"

# Download protein sequences with intelligent search strategy
# Multi-stage approach: UniProt primary search → gene alias fallback → NCBI bulk validation
# Employs persistent caching and smart completion detection to avoid redundant downloads
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
    group: "gram_group_{analysis}_{paramset}"
    script:
        "scripts/download_proteins.py"

# Download bacterial 3D structures and sequences from PDB
# Searches UniProt for bacterial proteins with available 3D structures, prioritizing by resolution
# Downloads top 3 structures per gene with concurrent processing and intelligent completion checking
rule download_3d_structures:
    input:
        protein_list="results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv"
    output:
        sentinel="data/proteins_3d_structure/.{analysis}_{paramset}_{group}_structures_complete"
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        max_structures=config.get("3d_structures", {}).get("max_structures", 3)
    group: "gram_group_{analysis}_{paramset}"
    script:
        "scripts/download_3d_structures.py"
    

# Select representative protein sequences for multiple sequence alignment
# Chooses one protein per species per gene

## Step 1
rule create_msa_sequence_references:
    input:
        gene_lists="data/{analysis}_{paramset}/genes_species/gram_{group}",
        protein_download_sentinel="data/proteins_fasta/.{analysis}_{paramset}_{group}_download_complete",
        structures_download_sentinel="data/proteins_3d_structure/.{analysis}_{paramset}_{group}_structures_complete"
    output:
        directory("results/{analysis}_{paramset}/msa_sequence_refs/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/create_sequence_references.py"

## Step 2
rule create_msa_fasta_files:
    input:
        reference_dir="results/{analysis}_{paramset}/msa_sequence_refs/gram_{group}"
    output:
        msa_no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_sequences/gram_{group}"),
        msa_with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_sequences/gram_{group}"),
        selected_3d_paths="results/{analysis}_{paramset}/selected_3d_paths_gram_{group}.txt",
        selected_3d_tsv="results/{analysis}_{paramset}/selected_3d_fasta_gram_{group}.tsv"
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/create_msa_fasta.py"

# Perform multiple sequence alignments using MAFFT algorithm
# Generates high-quality amino acid sequence alignments for each gene within Gram groups
# Employs multi-threading for efficient processing of large protein datasets
rule run_all_mafft_for_group:
    input:
        msa_no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_sequences/gram_{group}"),
        msa_with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_sequences/gram_{group}")
    output:
        no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_alignments/gram_{group}"),
        with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_alignments/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        threads=config["mafft"]["threads"]
    script:
        "scripts/run_mafft_alignments.py"


# Remove poorly aligned regions using trimAl automated trimming
# Eliminates gaps and ambiguously aligned positions to improve alignment quality
# Applies statistical algorithms to retain only reliably aligned sequence regions
rule trim_alignments:
    input:
        no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_alignments/gram_{group}"),
        with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_alignments/gram_{group}")
    output:
        trim_no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_trimmed/gram_{group}"),
        trim_with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_trimmed/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        use_3d=config.get("mafft", {}).get("use_3d_alignments", "no_3d")
    script:
        "scripts/trim_alignments.py"

# Evaluate multiple sequence alignment quality before and after trimming
# Compares alignment statistics (gaps, conservation, length) between raw and trimmed versions
# Provides quality metrics to guide selection of optimal alignments for downstream analysis
rule assess_alignment_quality:
    input:
        no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_alignments/gram_{group}"),
        with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_alignments/gram_{group}"),
        trim_no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_trimmed/gram_{group}"),
        trim_with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_trimmed/gram_{group}")
    output:
        no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_quality/gram_{group}"),
        with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_quality/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        use_3d=config.get("mafft", {}).get("use_3d_alignments", "no_3d")
    script:
        "scripts/assess_alignment_quality_comparison.py"

# Calculate amino acid conservation patterns from optimal alignments
# Identifies highly conserved residues and regions across aligned protein sequences
# Adaptively selects best alignment (raw or trimmed) based on quality assessment results
rule analyze_conservation:
    input:
        no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_alignments/gram_{group}"),
        with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_alignments/gram_{group}"),
        trim_no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_trimmed/gram_{group}"),
        trim_with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_trimmed/gram_{group}"),
        quality_assessment_no_3d=directory("results/{analysis}_{paramset}/no_3d/msa_quality/gram_{group}"),
        quality_assessment_with_3d=directory("results/{analysis}_{paramset}/with_3d/msa_quality/gram_{group}")
    output:
        no_3d=directory("results/{analysis}_{paramset}/no_3d/conservation/gram_{group}"),
        with_3d=directory("results/{analysis}_{paramset}/with_3d/conservation/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}",
        use_3d=config.get("mafft", {}).get("use_3d_alignments", "no_3d"),
        create_logos=False
    script:
        "scripts/analyze_conservation_adaptive.py"


# Predict B-cell epitopes using BepiPred 3.0 machine learning algorithm
# Employs state-of-the-art deep learning models for linear B-cell epitope identification
# Focuses on surface-accessible regions with high immunogenic potential

# Predict B-cell epitopes using BepiPred 3.0 on selected 3D structure sequences
rule predict_epitopes_bepipred:
    input:
        selected_3d_paths="results/{analysis}_{paramset}/selected_3d_paths_gram_{group}.txt"
    output:
        directory("results/{analysis}_{paramset}/epitope_predictions_bepipred/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/predict_epitopes_bepipred_3d_only.py"

# Create comprehensive protein download summary and statistics
# Compares actual protein recovery rates against expected coverage from database searches
# Identifies download gaps and success patterns across genes and species
rule generate_download_summary:
    input:
        coverage="results/{analysis}_{paramset}/coverage/coverage_count.tsv",
        protein_fasta="data/proteins_fasta",
        proteins_to_study="results/{analysis}_{paramset}/proteins_to_study/gram_{group}.tsv"
    output:
        directory("results/{analysis}_{paramset}/download_summary/gram_{group}")
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        group="{group}"
    script:
        "scripts/generate_download_summary.py"

# Compile integrated analysis report combining all pipeline results
# Synthesizes alignment quality, conservation patterns, and epitope predictions
# Produces HTML report with visualizations and summary statistics for both Gram groups
rule generate_final_report:
    input:
        quality_positive="results/{analysis}_{paramset}/{use_3d_dir}/msa_quality/gram_positive",
        quality_negative="results/{analysis}_{paramset}/{use_3d_dir}/msa_quality/gram_negative",
        conservation_positive="results/{analysis}_{paramset}/{use_3d_dir}/conservation/gram_positive",
        conservation_negative="results/{analysis}_{paramset}/{use_3d_dir}/conservation/gram_negative",
        download_summary_positive="results/{analysis}_{paramset}/download_summary/gram_positive",
        download_summary_negative="results/{analysis}_{paramset}/download_summary/gram_negative"
    output:
        "results/{analysis}_{paramset}/reports/final_report_{use_3d_dir}.html"
    params:
        analysis="{analysis}",
        paramset="{paramset}",
        use_3d_dir="{use_3d_dir}"
    script:
        "scripts/generate_final_report.py"