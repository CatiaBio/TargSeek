from Bio import SeqIO
import os
import re

# ---------------------- CONFIGURATION ----------------------
#input_folder = snakemake.input.input_folder
input_folder = "results/gram_positive/yfkN"
#output_folder = snakemake.output.output_folder
output_folder = "results/gram_positive/fasta_merged"
#group = snakemake.params.group
group = "gram_positive"

os.makedirs(output_folder, exist_ok=True)

# Automatically detect gene name from folder name
gene_name = os.path.basename(input_folder)

# ---------------------- FUNCTIONS ----------------------
def extract_species(description):
    match = re.search(r"\[(.*?)\]", description)
    if match:
        species = match.group(1)
        if len(species.split()) == 2:
            return species
    return None

def is_partial(description):
    return "partial" in description.lower()

# ---------------------- PROCESS ----------------------
species_to_record = {}

for subfile in os.listdir(input_folder):
    fasta_path = os.path.join(input_folder, subfile)
    if os.path.isdir(fasta_path):
        continue
    if not subfile.endswith(".fasta"):
        continue

    for record in SeqIO.parse(fasta_path, "fasta"):
        if is_partial(record.description):
            continue

        species = extract_species(record.description)
        if not species:
            continue

        if species not in species_to_record:
            species_to_record[species] = record
        else:
            if len(record.seq) > len(species_to_record[species].seq):
                species_to_record[species] = record

# Write merged fasta for this protein
if species_to_record:
    output_path = os.path.join(output_folder, f"{gene_name}.fasta")
    with open(output_path, "w") as f:
        SeqIO.write(species_to_record.values(), f, "fasta")
    print(f"Wrote {len(species_to_record)} sequences to {output_path} ({group} Clean)")
else:
    print(f"No sequences selected for {group} ({input_folder})")
