import sys
from pathlib import Path
import bacdive
import json
import os
from datetime import datetime

# -------------------------
# USE FOR TESTING OUTSIDE SNAKEMAKE
# -------------------------
# Uncomment the following lines to run standalone:
# bacdive_gram_file = Path("data/bacdive/gram_classification.tsv")
# not_found_file = Path("data/bacdive/not_found.txt")
# fallback_file = Path("config/microbiome/fallback_gram_lookup.tsv")
# updated_gram_file = Path("data/bacdive/updated_gram_classification.tsv")
# updated_not_found_file = Path("data/bacdive/updated_not_found.txt")

# -------------------------
# USE INSIDE SNAKEMAKE
# -------------------------
bacdive_gram_file = Path(snakemake.input.bacdive_classification)
not_found_file = Path(snakemake.input.not_found)
updated_gram_file = Path(snakemake.output.updated_classification)
updated_not_found_file = Path(snakemake.output.updated_not_found)

# -----------------------------
# BacDive authentication
# -----------------------------
try:
    bacdive_info_path = snakemake.input.bacdive_info
except NameError:
    bacdive_info_path = "config/login/bacdive_info.txt"  # fallback for standalone

with open(bacdive_info_path, "r", encoding="utf-8") as f:
    lines = f.readlines()
    email = lines[0].strip()
    password = lines[1].strip() if len(lines) > 1 else None

client = bacdive.BacdiveClient(email, password)

# --- Load not-found species ---
with not_found_file.open("r", encoding="utf-8") as f:
    missing_species = [line.strip() for line in f if line.strip()]

# --- Load existing BacDive Gram classifications ---
existing_gram_data = {}
with bacdive_gram_file.open("r", encoding="utf-8") as f:
    next(f)  # skip header
    for line in f:
        species, gram = line.strip().split("\t")
        existing_gram_data[species] = gram.lower()

print(f"Found {len(missing_species)} species in not_found.txt")
print(f"Existing gram data has {len(existing_gram_data)} species")

# --- Check the BacDive cache to see if missing species were actually resolved ---
cache_file = "cache/bacdive/bacdive_cache.json"
cache_resolved = []

if os.path.exists(cache_file):
    try:
        with open(cache_file, 'r', encoding='utf-8') as f:
            cache = json.load(f)
        
        # Check which "missing" species are actually in the cache with valid data
        for species in missing_species:
            if species in cache:
                cache_entry = cache[species]
                if (cache_entry.get("status") == "found" and 
                    "gram_classification" in cache_entry and 
                    cache_entry["gram_classification"] in ["positive", "negative"]):
                    
                    # This species was resolved via cache - add to existing data
                    existing_gram_data[species] = cache_entry["gram_classification"]
                    cache_resolved.append(species)
                    print(f"✓ Resolved {species} from cache: {cache_entry['gram_classification']}")
        
        print(f"Found {len(cache_resolved)} species resolved via cache")
        
    except Exception as e:
        print(f"Warning: Could not check cache: {e}")
        cache_resolved = []

# --- Remove cache-resolved species from missing list ---
truly_missing = [s for s in missing_species if s not in cache_resolved]
print(f"After cache check: {len(truly_missing)} species still need genus-based lookup")

# Load the full BacDive cache for updating
def load_bacdive_cache():
    """Load BacDive cache"""
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: Could not load BacDive cache: {e}")
            return {}
    return {}

def save_bacdive_cache(cache_data):
    """Save BacDive cache"""
    try:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(cache_data, f, indent=2, ensure_ascii=False)
        print(f"✓ Updated BacDive cache saved to: {cache_file}")
    except Exception as e:
        print(f"Warning: Could not save BacDive cache: {e}")

# Load the main BacDive cache for updating
cache = load_bacdive_cache()
print(f"Loaded BacDive cache with {len(cache)} entries for updating")

# --- Known gram-negative genera (expanded from analysis1 data) ---
KNOWN_GRAM_NEGATIVE_GENERA = {
    'achromobacter', 'acidovorax', 'acinetobacter', 'aeromonas', 'afipia', 'agrobacterium',
    'alcanivorax', 'aliiglaciecola', 'alishewanella', 'alteromonas', 'aquabacterium',
    'atlantibacter', 'aureimonas', 'azorhizobium', 'azospira', 'azospirillum',
    'bacteroides', 'bibersteinia', 'bosea', 'brachymonas', 'bradyrhizobium',
    'brevundimonas', 'burkholderia', 'campylobacter', 'cardiobacterium', 'cedecea',
    'chryseobacterium', 'citrobacter', 'comamonas', 'cronobacter', 'cupriavidus', 'curvibacter',
    'delftia', 'dickeya', 'elizabethkingia', 'empedobacter', 'enterobacter',
    'erwinia', 'erythrobacter', 'escherichia', 'franconibacter', 'fusobacterium',
    'gallaecimonas', 'haematobacter', 'haemophilus', 'hydrogenophaga', 'hyphomicrobium',
    'klebsiella', 'kluyvera', 'kosakonia', 'leclercia', 'luteimonas',
    'mangrovibacter', 'mannheimia', 'massilia', 'mesorhizobium', 'methylobacterium',
    'methylophaga', 'methylorubrum', 'microvirga', 'moraxella', 'morganella',
    'nitratireductor', 'obesumbacterium', 'ornithobacterium', 'pannonibacter', 'pantoea',
    'paraburkholderia', 'paracoccus', 'paraglaciecola', 'pelomonas', 'pluralibacter', 'porphyromonas',
    'prevotella', 'proteus', 'pseudescherichia', 'pseudochrobactrum', 'pseudomonas',
    'pseudopedobacter', 'psychrobacter', 'rahnella', 'ralstonia', 'rheinheimera',
    'rouxiella', 'salmonella', 'selenomonas', 'serratia', 'shewanella',
    'shigella', 'shinella', 'siccibacter', 'skermanella', 'sphingobacterium', 'sphingomonas',
    'spongiibacter', 'stenotrophomonas', 'tepidimonas', 'thauera', 'thermus',
    'veillonella', 'vibrio', 'weeksella', 'xanthomonas', 'yokenella'
}

# --- Known gram-positive genera (expanded from analysis1 data) ---
KNOWN_GRAM_POSITIVE_GENERA = {
    'actinomyces', 'actinopolyspora', 'aerococcus', 'arthrobacter', 'bacillus',
    'bavariicoccus', 'bifidobacterium', 'brachybacterium', 'catellicoccus', 'clostridium',
    'corynebacterium', 'curtobacterium', 'cutibacterium', 'deinococcus', 'dermacoccus',
    'dolosigranulum', 'enterococcus', 'frigoribacterium', 'gemella', 'geobacillus',
    'janibacter', 'jeotgalicoccus', 'knoellia', 'kocuria', 'kytococcus',
    'lactobacillus', 'lactococcus', 'leifsonia', 'leucobacter', 'leuconostoc', 'listeria',
    'marinococcus', 'microbacterium', 'micrococcus', 'micromonospora', 'mycobacterium',
    'neisseria', 'paenarthrobacter', 'paenibacillus', 'peptoniphilus', 'psychromicrobium',
    'rothia', 'staphylococcus', 'streptococcus', 'streptomyces', 'terrabacter',
    'tessaracoccus', 'trueperella', 'weissella'
}

# --- All conflicting genera have been resolved and moved to their respective known lists above ---

# --- Only do genus-based search for truly missing species ---
genus_cache = {}
updated_gram_data = existing_gram_data.copy()

for species in truly_missing:
    genus = species.split()[0].lower()

    # Check known genera first (prevents BacDive misclassifications)
    if genus in KNOWN_GRAM_NEGATIVE_GENERA:
        gram_result = "negative"
        updated_gram_data[species] = gram_result
        print(f"✓ {species} classified as negative (known gram-negative genus: {genus})")
    elif genus in KNOWN_GRAM_POSITIVE_GENERA:
        gram_result = "positive"
        updated_gram_data[species] = gram_result
        print(f"✓ {species} classified as positive (known gram-positive genus: {genus})")
    else:
        # Fall back to BacDive genus search for unknown genera
        try:
            if genus in genus_cache:
                records = genus_cache[genus]
            else:
                client.search(taxonomy=genus)
                records = list(client.retrieve())
                genus_cache[genus] = records

            strain_grams = set()

            def recurse(obj):
                if isinstance(obj, dict):
                    for k, v in obj.items():
                        if k.lower() == "gram stain":
                            if isinstance(v, str):
                                strain_grams.add(v.lower())
                            elif isinstance(v, list):
                                strain_grams.update(x.lower() for x in v if isinstance(x, str))
                        else:
                            recurse(v)
                elif isinstance(obj, list):
                    for item in obj:
                        recurse(item)

            for record in records:
                recurse(record)

            # Improved logic: prefer definitive results, avoid ambiguous data
            positive_count = sum(1 for g in strain_grams if "positive" in g)
            negative_count = sum(1 for g in strain_grams if "negative" in g)

            if negative_count > positive_count:
                gram_result = "negative"
                updated_gram_data[species] = gram_result
                print(f"✓ {species} classified as negative via BacDive (neg:{negative_count}, pos:{positive_count})")
            elif positive_count > negative_count:
                gram_result = "positive"
                updated_gram_data[species] = gram_result
                print(f"✓ {species} classified as positive via BacDive (neg:{negative_count}, pos:{positive_count})")
            else:
                gram_result = "unknown"
                updated_gram_data[species] = gram_result
                print(f"? {species} classification unclear via BacDive (neg:{negative_count}, pos:{positive_count})")

        except Exception as e:
            gram_result = "unknown"
            updated_gram_data[species] = gram_result
            print(f"✗ BacDive lookup failed for {species}: {e}")

    # Update cache with classification result
    try:

        # Update BacDive cache with genus-based result to prevent future lookups
        if species in cache:
            cache[species]["gram_classification"] = gram_result
            cache[species]["timestamp"] = datetime.now().isoformat()
            cache[species]["genus_lookup"] = True  # Mark as genus-based lookup
            print(f"✓ Updated BacDive cache for {species}: {gram_result} (genus-based)")
        else:
            # Create cache entry for species found via genus lookup
            cache[species] = {
                "status": "found_via_genus",
                "gram_classification": gram_result,
                "timestamp": datetime.now().isoformat(),
                "genus_lookup": True,
                "genus": genus
            }
            print(f"✓ Added to BacDive cache: {species}: {gram_result} (genus-based)")

    except Exception as e:
        updated_gram_data[species] = "unknown"
        
        # Cache the failed genus lookup to prevent retries
        if species in cache:
            cache[species]["gram_classification"] = "unknown"
            cache[species]["timestamp"] = datetime.now().isoformat()
            cache[species]["genus_lookup_failed"] = True
        else:
            cache[species] = {
                "status": "genus_lookup_failed",
                "gram_classification": "unknown", 
                "timestamp": datetime.now().isoformat(),
                "genus_lookup_failed": True,
                "genus": genus,
                "error": str(e)
            }
        print(f"✗ Genus lookup failed for {species}, cached as unknown")

# --- Save updated Gram classification with Unix line endings ---
with updated_gram_file.open("w", encoding="utf-8", newline="\n") as out:
    out.write("species\tgram_classification\n")
    for species in sorted(updated_gram_data):
        out.write(f"{species}\t{updated_gram_data[species]}\n")

print(f"✅ Updated Gram classification saved to: {updated_gram_file}")

# --- Update not_found file with still unknowns from original not_found list ---
still_unknown = [species for species in missing_species if updated_gram_data.get(species, "") == "unknown"]
with updated_not_found_file.open("w", encoding="utf-8") as nf:
    for species in sorted(still_unknown):
        nf.write(f"{species}\n")

print(f"📄 Updated list of unknowns saved to: {updated_not_found_file}")

# Save the updated BacDive cache with genus-based results
save_bacdive_cache(cache)

# Print summary
print(f"\n=== Supplementation Summary ===")
print(f"Original missing species: {len(missing_species)}")
print(f"Resolved via cache: {len(cache_resolved)}")
print(f"Required genus lookup: {len(truly_missing)}")
print(f"Still unknown: {len(still_unknown)}")
print(f"Total classified species: {len([s for s in updated_gram_data.values() if s != 'unknown'])}")
print(f"BacDive cache entries after update: {len(cache)}")