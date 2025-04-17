# ---------------------
# Final targets
# ---------------------
rule all:
    input:
        "results/proteins_unique_in_positive.txt",
        "results/proteins_unique_in_negative.txt"

# ---------------------
# Download BacDive metadata and classify taxa
# ---------------------
rule classify_taxa_by_gram:
    """Use BacDive API to classify taxa into Gram-positive/negative."""
    input:
        species="config/microbiome/cow_milk/unique_species.txt",
        bacdive_info="config/login/bacdive_info.txt"
    output:
        all_json="data/bacdive/all_species.json",
        not_found="data/bacdive/not_found.txt",
        errors="data/bacdive/errors.txt",
        downloaded="data/bacdive/downloaded.txt",
        gram_classification="data/bacdive/gram_classification.tsv"
    script:
        "scripts/bacdive_classification.py"

# ---------------------
# Split species into Gram-positive and Gram-negative lists
# ---------------------
rule split_species_by_gram:
    """Split species by Gram classification into two files without headers."""
    input:
        "data/bacdive/gram_classification.tsv"
    output:
        gram_positive="data/bacdive/gram_positive.txt",
        gram_negative="data/bacdive/gram_negative.txt"
    shell:
        """
        mkdir -p data/bacdive

        awk -F'\t' 'NR > 1 && $2 == "positive" {{ print $1 }}' {input} > {output.gram_positive}
        awk -F'\t' 'NR > 1 && $2 == "negative" {{ print $1 }}' {input} > {output.gram_negative}
        """

# ---------------------
# Download GO annotations & extract gene symbols
# ---------------------
rule fetch_quickgo_annotations:
    """Fetch GO annotations and extract gene symbols using QuickGO API."""
    input:
        go_ids="config/quickgo/go_ids.tsv",
        taxon_ids="config/quickgo/taxon_ids.tsv"
    output:
        annotations="data/quickgo/annotations_all.json",
        symbols="data/quickgo/gene_symbols.txt"
    script:
        "scripts/fetch_quickgo_data.py"

# ---------------------
# Download protein sequences by gene and species
# ---------------------
rule fetch_ncbi_proteins:
    """Download proteins for each gene/taxon pair using NCBI Entrez."""
    input:
        proteins="data/gene_symbols.txt",
        species=expand("data/bacdive/gram_{group}.txt", group=["positive", "negative"]),
        ncbi_info="config/login/ncbi_info.txt"
    output:
        complete_flag="data/proteins/.download_complete"
    script:
        "scripts/fetch_ncbi_proteins.py"

# ---------------------
# Check coverage: how many taxa have each gene?
# ---------------------
rule assess_gene_taxa_coverage:
    """Check how many taxa have a protein hit for each gene."""
    input:
        species="data/bacdive/gram_{group}.txt",
        genes="data/quickgo/gene_symbols.txt",
        ncbi_info="config/login/ncbi_info.txt"
    output:
        coverage="results/gene_coverage_gram_{group}.tsv"
    script:
        "scripts/gene_taxa_coverage.py"

# ---------------------
# Filter and sort gene coverage by threshold
# ---------------------
GRAM_THRESHOLDS = {
    "positive": 5,
    "negative": 7,
}

rule filter_and_sort_coverage:
    """Filter genes by abundance threshold and sort by number of taxa, keeping the header."""
    input:
        "results/gene_coverage_gram_{group}.tsv"
    output:
        "results/gene_coverage_gram_{group}_filtered.tsv"
    params:
        threshold=lambda wildcards: GRAM_THRESHOLDS[wildcards.group]
    shell:
        """
        (head -n 1 {input} && awk -F'\t' 'NR > 1 && $2 > {params.threshold}' {input} | sort -k2,2nr) > {output}
        """

# ---------------------
# Compare Gram groups to find unique gene symbols
# ---------------------
rule find_unique_genes:
    """Find genes unique to each Gram group by comparing filtered gene lists directly."""
    input:
        positive="results/gene_coverage_gram_positive_filtered.tsv",
        negative="results/gene_coverage_gram_negative_filtered.tsv"
    output:
        unique_in_positive="results/unique_in_positive.txt",
        unique_in_negative="results/unique_in_negative.txt"
    shell:
        """
        comm -13 <(cut -f1 {input.negative} | tail -n +2 | sort) <(cut -f1 {input.positive} | tail -n +2 | sort) > {output.unique_in_positive}
        comm -23 <(cut -f1 {input.negative} | tail -n +2 | sort) <(cut -f1 {input.positive} | tail -n +2 | sort) > {output.unique_in_negative}
        """
    
# ---------------------
# Fetch protein names for unique gene sets
# ---------------------
rule get_protein_names_for_unique_genes:
    """Search NCBI for protein names using gene names only."""
    input:
        gene_file="results/unique_in_{group}.txt",
        ncbi_info="config/login/ncbi_info.txt"
    output:
        protein_names="results/proteins_unique_in_{group}.txt"
    script:
        "scripts/get_protein_names_from_ncbi.py"
