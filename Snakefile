# Load configuration
configfile: "config/config.yaml"

# Helper function to load selected genes dynamically
def get_selected_genes(group):
    file_path = f"results/proteins_unique_in_{group}.txt"
    with open(file_path) as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes

# ---------------------
# Final Target
# ---------------------
rule all:
    input:
        expand(
            "results/positive_proteins_aligned/{gene}_aligned.fasta",
            gene=get_selected_genes("positive")
        ) +
        expand(
            "results/negative_proteins_aligned/{gene}_aligned.fasta",
            gene=get_selected_genes("negative")
        )

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
rule fetch_quickgo_annotations:
    input:
        config["quickgo"]["go_ids"],
        config["quickgo"]["taxon_ids"]
    output:
        config["quickgo"]["annotations"],
        config["quickgo"]["genes"]
    script:
        config["scripts"]["fetch_quickgo"]

# Assess how many species have each gene annotated
rule assess_gene_taxa_coverage:
    input:
        species="data/bacdive/gram_{group}.txt",
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

# Download protein sequences from NCBI for selected genes and species
rule download_proteins_to_analyse:
    input:
        proteins="results/proteins_unique_in_{group}.txt",
        species=lambda wildcards: config["bacdive"]["gram_species_template"].format(group=wildcards.group),
        ncbi_info=config["login"]["ncbi_info"]
    output:
        output_folder=directory("results/{group}_proteins_downloaded")
    params:
        group="{group}"
    script:
        config["scripts"]["download_proteins"]


# Download protein sequences from NCBI for selected genes and species
rule download_proteins_to_analyse:
    input:
        protein_selection="results/proteins_unique_in_{group}.txt",
        species="data/bacdive/gram_{group}.txt",
        ncbi_info=config["login"]["ncbi_info"]
    output:
        complete_flag="results/proteins_downloaded_{group}.txt"
    script:
        config["scripts"]["download_proteins"]

# Select one representative protein per species for MSA
rule select_proteins_for_msa:
    input:
        input_folder="results/{group}_proteins_downloaded",
    output:
        output_folder=directory("results/{group}_proteins_merged"),
    params:
        group="{group}"
    script:
        config["scripts"]["get_msa_sequences"]

# Perform multiple sequence alignment (MSA) using MAFFT with threading
rule run_mafft:
    input:
        fasta="results/{group}_proteins_merged/{gene}.fasta"
    output:
        msa="results/{group}_proteins_aligned/{gene}_aligned.fasta"
    threads: config["mafft"]["threads"]
    shell:
        """
        mkdir -p results/{wildcards.group}_proteins_aligned
        mafft --thread {threads} --maxiterate 10000 --localpair {input.fasta} > {output.msa}
        """

# Find conserved amino acid positions in the MSA
rule find_conserved_amino_acids:
    input:
        expand(
            "results/positive_proteins_aligned/{gene}_aligned.fasta",
            gene=get_selected_genes("positive")
        ) +
        expand(
            "results/negative_proteins_aligned/{gene}_aligned.fasta",
            gene=get_selected_genes("negative")
        )
    script:
        config["scripts"]["conserved_aminoacids"]


