import pandas as pd

# ----------------------
# Snakemake I/O
# ----------------------
lineage_file = snakemake.input.lineage
taxid_file = snakemake.input.taxids
output_species = snakemake.output.species

# ----------------------
# Load taxid list
# ----------------------
with open(taxid_file, "r") as f:
    bacvac_taxids = set(line.strip() for line in f if line.strip().isdigit())

# ----------------------
# Load lineage table
# ----------------------
lineage_df = pd.read_csv(lineage_file, sep="\t", dtype=str)

# ----------------------
# Filter and extract species names
# ----------------------
filtered_df = lineage_df[lineage_df["taxid"].isin(bacvac_taxids)]

species_names = filtered_df["species"].dropna().unique()
species_names_sorted = sorted(species_names)

# ----------------------
# Save output
# ----------------------
with open(output_species, "w") as f:
    for name in species_names_sorted:
        f.write(f"{name}\n")

print(f"✅ Saved {len(species_names_sorted)} species names to {output_species}")