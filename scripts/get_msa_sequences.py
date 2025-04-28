from Bio import SeqIO
import os
import re

# ---------------------- CONFIGURATION ----------------------
input_folder = snakemake.input.input_folder
output_folder = snakemake.output.output_folder
group = snakemake.params.group

os.makedirs(output_folder, exist_ok=True)

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

# ---------------------- PROCESS PER PROTEIN FOLDER ----------------------
for protein_folder in os.listdir(input_folder):
    full_folder = os.path.join(input_folder, protein_folder)
    if not os.path.isdir(full_folder):
        continue

    species_to_record = {}

    for subfile in os.listdir(full_folder):
        fasta_path = os.path.join(full_folder, subfile)
        if os.path.isdir(fasta_path):
            continue  # skip unexpected folders
        if not subfile.endswith(".fasta"):
            continue

        for record in SeqIO.parse(fasta_path, "fasta"):
            if is_partial(record.description):
                continue

            species = extract_species(record.description)
            if not species:
                continue

            # Select the sequence with the most amino acids for each species
            if species not in species_to_record:
                species_to_record[species] = record
            else:
                if len(record.seq) > len(species_to_record[species].seq):
                    species_to_record[species] = record

    # Write merged fasta for this protein
    if species_to_record:
        output_path = os.path.join(output_folder, f"{protein_folder}.fasta")
        with open(output_path, "w") as f:
            SeqIO.write(species_to_record.values(), f, "fasta")
        print(f"Wrote {len(species_to_record)} sequences for {protein_folder} ({group} Clean)")
    else:
        print(f"No sequences selected for {protein_folder} ({group} Clean)")
