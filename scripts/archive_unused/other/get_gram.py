import bacdive
import time

# ----------------------
# Snakemake I/O
# ----------------------
species_file = snakemake.input.species
output_positive = snakemake.output.positive
output_negative = snakemake.output.negative
output_unknown = snakemake.output.unknown
credentials_file = snakemake.input.bacdive_info  
# ----------------------
# BacDive Auth
# ----------------------
with open(credentials_file) as bacdive_info:
    lines = bacdive_info.readlines()
    email = lines[0].strip()
    password = lines[1].strip() if len(lines) > 1 else None
   
# Initialize the client
client = bacdive.BacdiveClient(email, password)

# ----------------------
# Load species list
# ----------------------
with open(species_file, "r") as f:
    species_list = [line.strip() for line in f if line.strip()]

# ----------------------
# Classify by Gram stain
# ----------------------
gram_positive = []
gram_negative = []
gram_unknown = []

for species in species_list:
    print(f"🔍 Querying: {species}")
    try:
        client.search(taxonomy=species)

        for strain in client.retrieve():
            general_info = strain.get("General", {})
            morphology_info = strain.get("Morphology", {})
            cell_morph = morphology_info.get("cell morphology", {})
            gram_stain = cell_morph.get("gram stain")

            if gram_stain == "positive":
                gram_positive.append(species)
                break
            elif gram_stain == "negative":
                gram_negative.append(species)
                break
        else:
            gram_unknown.append(species)

        time.sleep(1.5)

    except Exception as e:
        print(f"❌ Error with {species}: {e}")
        gram_unknown.append(species)

# ----------------------
# Save results
# ----------------------
def save_list(filename, data):
    with open(filename, "w") as f:
        for item in sorted(set(data)):
            f.write(f"{item}\n")

save_list(output_positive, gram_positive)
save_list(output_negative, gram_negative)
save_list(output_unknown, gram_unknown)

print("✅ Done! Results saved to:")
print("-", output_positive)
print("-", output_negative)
print("-", output_unknown)
