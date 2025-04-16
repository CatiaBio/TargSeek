import os

species_file = "data/bacdive/gram_stain/gram_negative.txt"
title_dir = "results/raw_titles"
output_file = "results/gene_species_filtered.tsv"

with open(species_file) as sf:
    species_list = [line.strip().lower() for line in sf if line.strip()]

results = []

for filename in os.listdir(title_dir):
    if not filename.endswith("_titles.txt"):
        continue

    gene = filename.replace("_titles.txt", "")
    path = os.path.join(title_dir, filename)

    with open(path) as f:
        titles = f.readlines()

    matched_species = set()
    for title in titles:
        for species in species_list:
            if species in title.lower():
                matched_species.add(species)

    results.append((gene, len(matched_species)))

# Save summary
with open(output_file, "w") as out:
    out.write("gene\tspecies_count\n")
    for gene, count in results:
        out.write(f"{gene}\t{count}\n")

print(f"✅ Filtered species counts saved to: {output_file}")
