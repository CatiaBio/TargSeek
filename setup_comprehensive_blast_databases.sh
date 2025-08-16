#!/bin/bash

# Comprehensive BLAST databases setup for epitope cross-reactivity analysis
# Downloads proteomes from multiple categories for thorough specificity testing

DB_DIR="data/blast_databases"
mkdir -p "$DB_DIR"

echo "Setting up comprehensive BLAST databases for epitope cross-reactivity analysis..."
echo "This will download ~15-20GB of data. Ensure you have sufficient disk space."
echo ""

# Function to download and format a proteome
download_and_format() {
    local name="$1"
    local url="$2"
    local description="$3"
    
    echo "🔽 Downloading $description ($name)..."
    
    if [ ! -f "$DB_DIR/${name}_proteome.fasta.gz" ]; then
        wget -O "$DB_DIR/${name}_proteome.fasta.gz" "$url"
        if [ $? -ne 0 ]; then
            echo "  ❌ Failed to download $name proteome"
            return 1
        fi
    else
        echo "  ✅ $name proteome already downloaded"
    fi
    
    # Decompress
    if [ ! -f "$DB_DIR/${name}_proteome.fasta" ]; then
        gunzip -c "$DB_DIR/${name}_proteome.fasta.gz" > "$DB_DIR/${name}_proteome.fasta"
        if [ $? -ne 0 ]; then
            echo "  ❌ Failed to decompress $name proteome"
            return 1
        fi
    fi
    
    # Format for BLAST
    if [ ! -f "$DB_DIR/${name}_proteome.fasta.phr" ]; then
        echo "  📚 Formatting $name proteome database..."
        makeblastdb -in "$DB_DIR/${name}_proteome.fasta" -dbtype prot -title "$description"
        if [ $? -ne 0 ]; then
            echo "  ❌ Failed to format $name proteome database"
            return 1
        fi
    else
        echo "  ✅ $name proteome already formatted"
    fi
    
    echo "  ✅ $description complete"
    echo ""
}

# ================================
# 1. STANDARD DATABASES (already set up)
# ================================
echo "=== STANDARD DATABASES ==="

# Human (Homo sapiens)
download_and_format "human" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz" \
    "Human Proteome"

# Viral proteins (NCBI RefSeq)
echo "🔽 Downloading Viral proteins..."
if [ ! -f "$DB_DIR/viral_proteins.fasta.gz" ]; then
    wget -O "$DB_DIR/viral_proteins.fasta.gz" \
    "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
fi

if [ ! -f "$DB_DIR/viral_proteins.fasta" ]; then
    gunzip -c "$DB_DIR/viral_proteins.fasta.gz" > "$DB_DIR/viral_proteins.fasta"
fi

if [ ! -f "$DB_DIR/viral_proteins.fasta.phr" ]; then
    echo "  📚 Formatting viral proteins database..."
    makeblastdb -in "$DB_DIR/viral_proteins.fasta" -dbtype prot -title "Viral Proteins"
fi

# Fungal proteins (NCBI RefSeq)
echo "🔽 Downloading Fungal proteins..."
if [ ! -f "$DB_DIR/fungal_proteins.fasta.gz" ]; then
    wget -O "$DB_DIR/fungal_proteins.fasta.gz" \
    "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.1.protein.faa.gz"
fi

if [ ! -f "$DB_DIR/fungal_proteins.fasta" ]; then
    gunzip -c "$DB_DIR/fungal_proteins.fasta.gz" > "$DB_DIR/fungal_proteins.fasta"
fi

if [ ! -f "$DB_DIR/fungal_proteins.fasta.phr" ]; then
    echo "  📚 Formatting fungal proteins database..."
    makeblastdb -in "$DB_DIR/fungal_proteins.fasta" -dbtype prot -title "Fungal Proteins"
fi

# ================================
# 2. FARM ANIMALS & LIVESTOCK
# ================================
echo "=== FARM ANIMALS & LIVESTOCK ==="

# Cow (Bos taurus) - most important for milk
download_and_format "cow" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000009136/UP000009136_9913.fasta.gz" \
    "Cow Proteome"

# Goat (Capra hircus)
download_and_format "goat" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000291000/UP000291000_9925.fasta.gz" \
    "Goat Proteome"

# Horse (Equus caballus)
download_and_format "horse" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002281/UP000002281_9796.fasta.gz" \
    "Horse Proteome"

# Pig (Sus scrofa)
download_and_format "pig" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000008227/UP000008227_9823.fasta.gz" \
    "Pig Proteome"

# Chicken (Gallus gallus)
download_and_format "chicken" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000539/UP000000539_9031.fasta.gz" \
    "Chicken Proteome"

# Water buffalo (Bubalus bubalis)
download_and_format "buffalo" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000006522/UP000006522_89462.fasta.gz" \
    "Water Buffalo Proteome"

# ================================
# 3. RESEARCH MODEL ORGANISMS
# ================================
echo "=== RESEARCH MODEL ORGANISMS ==="

# Mouse (Mus musculus)
download_and_format "mouse" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.fasta.gz" \
    "Mouse Proteome"

# Rat (Rattus norvegicus)
download_and_format "rat" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002494/UP000002494_10116.fasta.gz" \
    "Rat Proteome"

# ================================
# 4. AQUATIC FOOD ORGANISMS
# ================================
echo "=== AQUATIC FOOD ORGANISMS ==="

# Atlantic salmon (Salmo salar)
download_and_format "salmon" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000087266/UP000087266_8030.fasta.gz" \
    "Atlantic Salmon Proteome"

# Zebrafish (Danio rerio) - model organism
download_and_format "zebrafish" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000437/UP000000437_7955.fasta.gz" \
    "Zebrafish Proteome"

# ================================
# 5. PLANT CROPS (Food allergens/contamination)
# ================================
echo "=== PLANT CROPS ==="

# Arabidopsis thaliana (model plant)
download_and_format "arabidopsis" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000006548/UP000006548_3702.fasta.gz" \
    "Arabidopsis Proteome"

# Rice (Oryza sativa)
download_and_format "rice" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000059680/UP000059680_39947.fasta.gz" \
    "Rice Proteome"

# Soybean (Glycine max)
download_and_format "soybean" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000008827/UP000008827_3847.fasta.gz" \
    "Soybean Proteome"

# ================================
# 6. BACTERIAL PATHOGENS
# ================================
echo "=== BACTERIAL PATHOGENS ==="

# E. coli (Escherichia coli K-12)
download_and_format "ecoli" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz" \
    "E. coli K-12 Proteome"

# Salmonella (Salmonella enterica)
download_and_format "salmonella" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000541/UP000000541_28901.fasta.gz" \
    "Salmonella enterica Proteome"

# Staphylococcus aureus
download_and_format "staph" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000008816/UP000008816_93061.fasta.gz" \
    "Staphylococcus aureus Proteome"

echo ""
echo "==============================================="
echo "📊 BLAST DATABASES SETUP COMPLETE!"
echo "==============================================="
echo ""
echo "Database categories installed:"
echo "  🧑 Human & Research Models: human, mouse, rat"
echo "  🐄 Farm Animals: cow, goat, horse, pig, chicken, buffalo"
echo "  🐟 Aquatic: salmon, zebrafish"
echo "  🌱 Plants: arabidopsis, rice, soybean"
echo "  🦠 Pathogens: ecoli, salmonella, staph"
echo "  🦠 Viral & Fungal: viral_proteins, fungal_proteins"
echo ""
echo "Database locations in: $PWD/$DB_DIR/"
echo ""
echo "💡 Usage example (milk project - comprehensive analysis):"
echo "python scripts/protein_analysis/create_comprehensive_epitope_visualizations.py \\"
echo "  --gene bamA --structure-id 5OR1 --gram-type negative \\"
echo "  --blast-db-human $PWD/$DB_DIR/human_proteome.fasta \\"
echo "  --blast-db-viral $PWD/$DB_DIR/viral_proteins.fasta \\"
echo "  --blast-db-fungal $PWD/$DB_DIR/fungal_proteins.fasta \\"
echo "  --blast-db-cow $PWD/$DB_DIR/cow_proteome.fasta \\"
echo "  --blast-db-goat $PWD/$DB_DIR/goat_proteome.fasta"
echo ""
echo "💡 For full cross-reactivity screening, add additional databases as needed."
echo "Total databases: ~18 proteomes covering major organism groups"