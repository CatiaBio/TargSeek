import csv
from collections import defaultdict

# ========== CONFIGURATION ========== #
INPUT_FILE = "data/data.tsv"
OUTPUT_FILE = "gene_summary.tsv"
EXPECTED_COLUMNS = [
    "Entry", "Entry Name", "Protein names", "Gene Names", "Organism",
    "Gene Ontology (cellular component)", "Organism (ID)"
]

# ========== DATA STRUCTURE ========== #
# gene → set of unique organisms
gene_data = defaultdict(set)

# ========== MAIN PROCESSING ========== #

with open(INPUT_FILE, newline='', encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter='\t')

    if not all(col in reader.fieldnames for col in EXPECTED_COLUMNS):
        print("❌ File header does not match expected UniProt TSV format.")
        exit(1)

    for row in reader:
        gene_names = row["Gene Names"].strip()
        organism = row["Organism"].strip()

        if not gene_names or not organism:
            continue

        gene = gene_names.split()[0].strip()
        organism = organism.replace("\n", " ").strip()

        gene_data[gene].add(organism)

print(f"✅ Processed {len(gene_data)} unique genes.")

# ========== SAVE OUTPUT (Sorted by Count Descending) ========== #

with open(OUTPUT_FILE, "w", encoding="utf-8", newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["Gene", "Count", "Organisms"])

    for gene in sorted(gene_data, key=lambda g: len(gene_data[g]), reverse=True):
        organisms = sorted(gene_data[gene])
        organism_string = ", ".join(organisms)
        writer.writerow([gene, len(organisms), organism_string])

print(f"✅ Summary saved to {OUTPUT_FILE}")
