import os
import json

# ----------------------
# Snakemake I/O
# ----------------------
json_dir = snakemake.input.json_dir
not_found_path = snakemake.input.not_found
output_positive = snakemake.output.positive
output_negative = snakemake.output.negative
output_unknown = snakemake.output.unknown

# ----------------------
# Read not found species list
# ----------------------
with open(not_found_path, "r") as nf:
    not_found = {line.strip() for line in nf if line.strip()}

# Format for matching JSON filenames
not_found_safe = {name.replace(" ", "_") + ".json" for name in not_found}

# ----------------------
# Gram stain extraction helpers
# ----------------------
def extract_gram_from_morphology(entry):
    try:
        morph = entry.get("Morphology", {})
        cell_morph = morph.get("cell morphology", {})

        if isinstance(cell_morph, list):
            return [m.get("gram stain", "").lower() for m in cell_morph if isinstance(m, dict)]
        elif isinstance(cell_morph, dict):
            gram = cell_morph.get("gram stain", "").lower()
            return [gram] if gram else []
    except Exception:
        pass
    return []

def search_gram_stain_recursively(entry):
    grams = set()
    def recurse(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k.lower() == "gram stain":
                    if isinstance(v, str):
                        grams.add(v.lower())
                    elif isinstance(v, list):
                        grams.update(x.lower() for x in v if isinstance(x, str))
                else:
                    recurse(v)
        elif isinstance(obj, list):
            for item in obj:
                recurse(item)
    recurse(entry)
    return list(grams)

# ----------------------
# Classify all species with valid .json
# ----------------------
positive, negative, unknown = [], [], []

for file in sorted(os.listdir(json_dir)):
    if not file.endswith(".json"):
        continue
    if file in not_found_safe:
        continue

    species = file.replace("_", " ").replace(".json", "")
    json_path = os.path.join(json_dir, file)

    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, list):
            data = [data]

        gram_result = "unknown"

        for entry in data:
            grams = extract_gram_from_morphology(entry)
            if not grams:
                grams = search_gram_stain_recursively(entry)

            if "positive" in grams:
                gram_result = "positive"
                break
            elif "negative" in grams:
                gram_result = "negative"
                break

        if gram_result == "positive":
            positive.append(species)
        elif gram_result == "negative":
            negative.append(species)
        else:
            unknown.append(species)

    except Exception as e:
        print(f"❌ Error reading {file}: {e}")
        unknown.append(species)

# ----------------------
# Save classified lists
# ----------------------
def write_list(path, values):
    with open(path, "w") as f:
        for item in sorted(values):
            f.write(f"{item}\n")

write_list(output_positive, positive)
write_list(output_negative, negative)
write_list(output_unknown, unknown)

print("✅ Gram classification complete.")
print(f"  Positive: {len(positive)}")
print(f"  Negative: {len(negative)}")
print(f"  Unknown : {len(unknown)}")