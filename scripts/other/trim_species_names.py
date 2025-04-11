from collections import Counter

input_file = "data/microbiome_taxids/merged_species_taxids.tsv"
trimmed_output_file  = "data/microbiome_taxids/merged_species_taxids_trimmed.tsv"
species_list_file = "data/microbiome_taxids/all_species.txt"
taxid_list_file = "data/microbiome_taxids/all_taxids.txt"
unique_species_file = "data/microbiome_taxids/unique_species.txt"

species_list = []
taxid_list = []

with open(input_file, "r", encoding="utf-8") as infile, \
     open(trimmed_output_file, "w", encoding="utf-8") as trimmed_out:

    for line in infile:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            full_name = parts[0].strip()
            taxid = parts[1].strip()

            # Trim to first two words
            trimmed_name = " ".join(full_name.split()[:2])

            # Write trimmed line to output
            trimmed_out.write(f"{trimmed_name}\t{taxid}\n")

            # Collect species and taxid
            species_list.append(trimmed_name)
            taxid_list.append(taxid)

# Write full species list (may include duplicates)
with open(species_list_file, "w", encoding="utf-8") as species_out:
    for name in species_list:
        species_out.write(name + "\n")

# Write full taxid list
with open(taxid_list_file, "w", encoding="utf-8") as taxid_out:
    for tid in taxid_list:
        taxid_out.write(tid + "\n")

# Count species occurrences
species_counter = Counter(species_list)
duplicates = {name: count for name, count in species_counter.items() if count > 1}

# Get all species (including duplicates), but remove repeated ones
unique_species = sorted(set(species_list))

# Save unique species names (with one entry per name)
with open(unique_species_file, "w", encoding="utf-8") as unique_out:
    for name in unique_species:
        unique_out.write(name + "\n")

# Report
if duplicates:
    print("⚠️ Duplicated species names found in the trimmed list:")
    for name, count in duplicates.items():
        print(f"- {name}: {count} times")
else:
    print("✅ No duplicates found in trimmed species names.")

print(f"🧾 Unique species (1 per name) saved to: {unique_species_file}")