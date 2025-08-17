#!/bin/bash

# Script to download and format BLAST databases for epitope specificity analysis
# Run this script to set up the databases needed for BLAST analysis

DB_DIR="data/blast_databases"
mkdir -p "$DB_DIR"

echo "Setting up BLAST databases for epitope specificity analysis..."

# 1. Human proteome (from UniProt)
echo "Downloading human proteome..."
if [ ! -f "$DB_DIR/human_proteome.fasta.gz" ]; then
    wget -O "$DB_DIR/human_proteome.fasta.gz" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
fi

# Decompress and format
if [ ! -f "$DB_DIR/human_proteome.fasta" ]; then
    gunzip -c "$DB_DIR/human_proteome.fasta.gz" > "$DB_DIR/human_proteome.fasta"
fi

if [ ! -f "$DB_DIR/human_proteome.fasta.phr" ]; then
    echo "Formatting human proteome database..."
    makeblastdb -in "$DB_DIR/human_proteome.fasta" -dbtype prot -title "Human Proteome"
fi

# 2. Viral proteins (from NCBI)
echo "Downloading viral proteins..."
if [ ! -f "$DB_DIR/viral_proteins.fasta.gz" ]; then
    wget -O "$DB_DIR/viral_proteins.fasta.gz" \
    "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
fi

# Decompress and format
if [ ! -f "$DB_DIR/viral_proteins.fasta" ]; then
    gunzip -c "$DB_DIR/viral_proteins.fasta.gz" > "$DB_DIR/viral_proteins.fasta"
fi

if [ ! -f "$DB_DIR/viral_proteins.fasta.phr" ]; then
    echo "Formatting viral proteins database..."
    makeblastdb -in "$DB_DIR/viral_proteins.fasta" -dbtype prot -title "Viral Proteins"
fi

# 3. Fungal proteins (from NCBI RefSeq)
echo "Downloading fungal proteins..."
if [ ! -f "$DB_DIR/fungal_proteins.fasta.gz" ]; then
    wget -O "$DB_DIR/fungal_proteins.fasta.gz" \
    "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.1.protein.faa.gz"
fi

# Decompress and format
if [ ! -f "$DB_DIR/fungal_proteins.fasta" ]; then
    gunzip -c "$DB_DIR/fungal_proteins.fasta.gz" > "$DB_DIR/fungal_proteins.fasta"
fi

if [ ! -f "$DB_DIR/fungal_proteins.fasta.phr" ]; then
    echo "Formatting fungal proteins database..."
    makeblastdb -in "$DB_DIR/fungal_proteins.fasta" -dbtype prot -title "Fungal Proteins"
fi

# 4. Farm animal proteomes (crucial for milk project)
echo "Downloading farm animal proteomes..."

# Cow (Bos taurus) - most important for milk project
if [ ! -f "$DB_DIR/cow_proteome.fasta.gz" ]; then
    echo "  Downloading cow proteome..."
    wget -O "$DB_DIR/cow_proteome.fasta.gz" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000009136/UP000009136_9913.fasta.gz"
fi

if [ ! -f "$DB_DIR/cow_proteome.fasta" ]; then
    gunzip -c "$DB_DIR/cow_proteome.fasta.gz" > "$DB_DIR/cow_proteome.fasta"
fi

if [ ! -f "$DB_DIR/cow_proteome.fasta.phr" ]; then
    echo "  Formatting cow proteome database..."
    makeblastdb -in "$DB_DIR/cow_proteome.fasta" -dbtype prot -title "Cow Proteome"
fi

# Sheep (Ovis aries)
if [ ! -f "$DB_DIR/sheep_proteome.fasta.gz" ]; then
    echo "  Downloading sheep proteome..."
    wget -O "$DB_DIR/sheep_proteome.fasta.gz" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000772810/UP000772810_9940.fasta.gz"
fi

if [ ! -f "$DB_DIR/sheep_proteome.fasta" ]; then
    gunzip -c "$DB_DIR/sheep_proteome.fasta.gz" > "$DB_DIR/sheep_proteome.fasta"
fi

if [ ! -f "$DB_DIR/sheep_proteome.fasta.phr" ]; then
    echo "  Formatting sheep proteome database..."
    makeblastdb -in "$DB_DIR/sheep_proteome.fasta" -dbtype prot -title "Sheep Proteome"
fi

# Goat (Capra hircus)
if [ ! -f "$DB_DIR/goat_proteome.fasta.gz" ]; then
    echo "  Downloading goat proteome..."
    wget -O "$DB_DIR/goat_proteome.fasta.gz" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000291000/UP000291000_9925.fasta.gz"
fi

if [ ! -f "$DB_DIR/goat_proteome.fasta" ]; then
    gunzip -c "$DB_DIR/goat_proteome.fasta.gz" > "$DB_DIR/goat_proteome.fasta"
fi

if [ ! -f "$DB_DIR/goat_proteome.fasta.phr" ]; then
    echo "  Formatting goat proteome database..."
    makeblastdb -in "$DB_DIR/goat_proteome.fasta" -dbtype prot -title "Goat Proteome"
fi

echo ""
echo "BLAST databases setup complete!"
echo ""
echo "Database locations:"
echo "  Human: $PWD/$DB_DIR/human_proteome.fasta"
echo "  Viral: $PWD/$DB_DIR/viral_proteins.fasta"
echo "  Fungal: $PWD/$DB_DIR/fungal_proteins.fasta"
echo "  Cow: $PWD/$DB_DIR/cow_proteome.fasta"
echo "  Sheep: $PWD/$DB_DIR/sheep_proteome.fasta"
echo "  Goat: $PWD/$DB_DIR/goat_proteome.fasta"
echo ""
echo "Usage examples:"
echo ""
echo "# Standard analysis (human, viral, fungal):"
echo "python scripts/protein_analysis/create_comprehensive_epitope_visualizations.py \\"
echo "  --gene bamA --structure-id 5OR1 --gram-type negative \\"
echo "  --blast-db-human $PWD/$DB_DIR/human_proteome.fasta \\"
echo "  --blast-db-viral $PWD/$DB_DIR/viral_proteins.fasta \\"
echo "  --blast-db-fungal $PWD/$DB_DIR/fungal_proteins.fasta"
echo ""
echo "# Milk project analysis (include farm animals):"
echo "python scripts/protein_analysis/create_comprehensive_epitope_visualizations.py \\"
echo "  --gene bamA --structure-id 5OR1 --gram-type negative \\"
echo "  --blast-db-human $PWD/$DB_DIR/human_proteome.fasta \\"
echo "  --blast-db-viral $PWD/$DB_DIR/viral_proteins.fasta \\"
echo "  --blast-db-fungal $PWD/$DB_DIR/fungal_proteins.fasta \\"
echo "  --blast-db-cow $PWD/$DB_DIR/cow_proteome.fasta \\"
echo "  --blast-db-sheep $PWD/$DB_DIR/sheep_proteome.fasta \\"
echo "  --blast-db-goat $PWD/$DB_DIR/goat_proteome.fasta"