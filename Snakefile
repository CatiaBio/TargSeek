# ---------------------
# Final targets
# ---------------------


# ---------------------
# Download taxonomy databases from NCBI
# ---------------------
rule download_taxdump_accession2taxid:
    output: 
        taxdump="other/taxdump.tar.gz",
        accession2taxid="other/nucl_gb.accession2taxid.gz"
    shell: 
        """
        mkdir -p other

        wget -q -O {output.taxdump} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        wget -q -O {output.accession2taxid} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        """

# ---------------------
# Unpack the taxdump files
# ---------------------
rule unpack_taxdump:
    input:
        "other/taxdump.tar.gz"
    output:
        "other/nodes.dmp",
        "other/names.dmp"
    shell:
        """
        mkdir -p other
        tar -zxvf {input} -C other nodes.dmp names.dmp
        """

# ---------------------
# Generate taxonomy and lineage files
# ---------------------
rule generate_lineage_taxonomy:
    input:
        "other/nodes.dmp",
        "other/names.dmp"
    output:
        "other/taxonomy.tsv",
        "other/lineage.tsv"
    shell:
        """
        python scripts/get_taxonomy_lineage.py --taxid 2
        """
        
rule download_bacdive_data:
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
        "scripts/bacdive_data.py"
        
# ---------------------
# Download GO annotations and extract gene symbols
# ---------------------
rule download_quickgo_data:
    input:
        go_ids="config/quickgo/go_ids_test.txt",
        taxon_ids="config/quickgo/taxon_ids_test.txt"
    output:
        annotations="data/quickgo_test/annotations_all.json",
        symbols="data/quickgo_test/gene_symbols.txt"
    script:
        "scripts/quickgo_data.py"


# ---------------------
# Get protein names from gene symbols
# ---------------------
rule get_protein_names:
    input:
        ncbi_info="config/login/ncbi_info.txt",
        gene_file="data/quickgo_test/gene_symbols.txt"
    output:
        protein_names="data/quickgo_test/protein_name.tsv"
    params:
        species="bacteria"
    script:
        "scripts/get_protein_name_from_gene.py"


# ---------------------
# Use gram information and protein gene information to get the protein sequences
# ---------------------
rule download_proteins_by_gene:
    input:
        proteins="data/quickgo/gene_symbols.txt",
        species="data/bacdive/gram_negative.txt",
        ncbi_info="config/login/ncbi_info.txt"
    output:
        complete_flag="data/proteins/gram_negative/.download_complete"
    script:
        "scripts/download_protein_per_species.py"

# # ---------------------
# # Run MAFFT on all downloaded FASTAs (after ALL downloads)
# # ---------------------
# rule msa_protein_alignment:
#     input:
#         done_flag="logs/_all_downloads_complete.txt"
#     output:
#         "results/{query}_alignment.fasta"
#     params:
#         data_dir=lambda wildcards: f"data/{wildcards.query}"
#     shell:
#         """
#         mkdir -p results
#         gunzip -c {params.data_dir}/*.fasta.gz > results/{wildcards.query}_combined.fasta
#         mafft --auto results/{wildcards.query}_combined.fasta > {output}
#         """
