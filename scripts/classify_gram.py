#!/usr/bin/env python3

import os
import json
import time
import bacdive
import shutil
from datetime import datetime

# ----------------------
# Snakemake I/O
# ----------------------
species_file = snakemake.input.species
credentials_file = snakemake.input.bacdive_info
output_json_path = snakemake.output[0]        # all_json
not_found_path = snakemake.output[1]          # not_found
output_errors_path = snakemake.output[2]      # errors
output_downloaded_path = snakemake.output[3]  # downloaded
output_gram_path = snakemake.output[4]        # gram_classification

# ----------------------
# BacDive Auth
# ----------------------
with open(credentials_file) as f:
    lines = f.readlines()
    email = lines[0].strip()
    password = lines[1].strip() if len(lines) > 1 else None

client = bacdive.BacdiveClient(email, password)

# ----------------------
# Load species list
# ----------------------
with open(species_file, "r") as f:
    species_list = [line.strip() for line in f if line.strip()]

# ----------------------
# Cache management
# ----------------------
cache_file = "cache/bacdive/bacdive_cache.json"
os.makedirs(os.path.dirname(cache_file), exist_ok=True)

def load_cache():
    """Load existing cache or return empty dict"""
    if os.path.exists(cache_file):
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, FileNotFoundError):
            return {}
    return {}

def save_cache(cache_data):
    """Save cache to file"""
    with open(cache_file, "w", encoding="utf-8") as f:
        json.dump(cache_data, f, indent=2, ensure_ascii=False)

def is_cache_valid(cache_entry, max_age_days=180):
    """Check if cache entry is still valid (not older than max_age_days)"""
    if "timestamp" not in cache_entry:
        return False
    
    try:
        cache_time = datetime.fromisoformat(cache_entry["timestamp"])
        age_days = (datetime.now() - cache_time).days
        return age_days <= max_age_days
    except (ValueError, TypeError):
        return False

# Load existing cache
cache = load_cache()

# ----------------------
# Helper to flatten retrieve()
# ----------------------
def flatten_strain_entries(strains):
    for s in strains:
        if isinstance(s, dict):
            yield s
        elif isinstance(s, list):
            for item in flatten_strain_entries(s):
                yield item

# ----------------------
# Data structures
# ----------------------
not_found = []
errors = []
downloaded = []
combined_data = {}

# Process species sequentially to avoid authentication issues
for species in species_list:
    # Check cache first
    if species in cache and is_cache_valid(cache[species]):
        cache_entry = cache[species]
        if cache_entry["status"] == "found":
            combined_data[species] = cache_entry["data"]
            downloaded.append(species)
            print(f"Using cached data for {species}")
            continue
        elif cache_entry["status"] == "not_found":
            not_found.append(species)
            print(f"Using cached 'not found' for {species}")
            continue
        elif cache_entry["status"] == "error":
            not_found.append(species)
            errors.append(f"{species}\t{cache_entry.get('error', 'Cached error')}")
            print(f"Using cached error for {species}")
            continue
    
    # Not in cache or cache expired - fetch from API
    try:
        print(f"Fetching {species} from BacDive API...")
        client.search(taxonomy=species)
        raw_data = list(client.retrieve())
        all_strains = list(flatten_strain_entries(raw_data))

        if not all_strains:
            not_found.append(species)
            # Cache the not found result
            cache[species] = {
                "status": "not_found",
                "timestamp": datetime.now().isoformat()
            }
            continue

        combined_data[species] = all_strains
        downloaded.append(species)
        
        # Cache the successful result
        cache[species] = {
            "status": "found",
            "data": all_strains,
            "timestamp": datetime.now().isoformat()
        }
        
        time.sleep(0.5)  # Reduced from 1.5s to 0.5s

    except Exception as e:
        not_found.append(species)
        errors.append(f"{species}\t{str(e)}")
        
        # Cache the error result
        cache[species] = {
            "status": "error",
            "error": str(e),
            "timestamp": datetime.now().isoformat()
        }

# Save updated cache
save_cache(cache)


with open(output_json_path, "w", encoding="utf-8") as out_json:
    json.dump(combined_data, out_json, indent=2, ensure_ascii=False)

# ----------------------
# Save species that raised errors
# ----------------------
with open(output_errors_path, "w", encoding="utf-8") as ef:
    for err in errors:
        ef.write(f"{err}\n")

# ----------------------
# Save successfully downloaded species
# ----------------------
with open(output_downloaded_path, "w", encoding="utf-8") as df:
    for species in sorted(downloaded):
        df.write(f"{species}\n")

# ----------------------
# Gram stain extraction helpers
# ----------------------
def extract_gram_from_morphology(entry):
    try:
        morph = entry.get("Morphology", {})
        cell_morph = morph.get("cell morphology", {})

        if isinstance(cell_morph, list):
            grams = []
            for m in cell_morph:
                if isinstance(m, dict):
                    gram = m.get("gram stain", "").lower()
                    if gram:
                        grams.append(gram)
                        # Early exit if we found a definitive result
                        if "positive" in gram or "negative" in gram:
                            return grams
            return grams
        elif isinstance(cell_morph, dict):
            gram = cell_morph.get("gram stain", "").lower()
            return [gram] if gram else []
    except Exception:
        pass
    return []

def search_gram_stain_recursively(entry):
    grams = set()
    found_definitive = False
    
    def recurse(obj):
        nonlocal found_definitive
        if found_definitive:
            return
            
        if isinstance(obj, dict):
            for k, v in obj.items():
                if found_definitive:
                    return
                if k.lower() == "gram stain":
                    if isinstance(v, str):
                        gram = v.lower()
                        grams.add(gram)
                        if "positive" in gram or "negative" in gram:
                            found_definitive = True
                            return
                    elif isinstance(v, list):
                        for x in v:
                            if isinstance(x, str):
                                gram = x.lower()
                                grams.add(gram)
                                if "positive" in gram or "negative" in gram:
                                    found_definitive = True
                                    return
                else:
                    recurse(v)
        elif isinstance(obj, list):
            for item in obj:
                if found_definitive:
                    return
                recurse(item)
    
    recurse(entry)
    return list(grams)

# ----------------------
# Classify Gram stain
# ----------------------
gram_classification = {}

for species, entries in combined_data.items():
    # Check if gram classification is already cached
    if (species in cache and 
        is_cache_valid(cache[species]) and 
        cache[species]["status"] == "found" and 
        "gram_classification" in cache[species]):
        
        gram_result = cache[species]["gram_classification"]
        print(f"Using cached gram classification for {species}: {gram_result}")
    else:
        # Classify gram stain from data
        gram_result = "unknown"
        for entry in entries:
            grams = extract_gram_from_morphology(entry)
            if not grams:
                grams = search_gram_stain_recursively(entry)
            if "positive" in grams:
                gram_result = "positive"
                break
            elif "negative" in grams:
                gram_result = "negative"
                break
        
        # Update cache with gram classification
        if species in cache and cache[species]["status"] == "found":
            cache[species]["gram_classification"] = gram_result
            cache[species]["timestamp"] = datetime.now().isoformat()

    if gram_result == "unknown":
        not_found.append(species)  # ✅ move to not_found
    else:
        gram_classification[species] = gram_result

# Save updated cache with gram classifications
save_cache(cache)

# ----------------------
# Save Gram classification as TSV (skip unknowns)
# ----------------------
with open(output_gram_path, "w", encoding="utf-8") as gf:
    gf.write("species\tgram_classification\n")
    for species, gram in sorted(gram_classification.items()):
        gf.write(f"{species}\t{gram}\n")

# ----------------------
# Save not found or unknown Gram species
# ----------------------
with open(not_found_path, "w", encoding="utf-8") as nf:
    for name in sorted(set(not_found)):  # de-duplicate
        nf.write(f"{name}\n")

# ----------------------
# Save all identified species (both Gram-positive and Gram-negative)
# ----------------------
all_identified_path = output_gram_path.replace("gram.tsv", "all_identified.txt")
with open(all_identified_path, "w", encoding="utf-8") as af:
    for species in sorted(gram_classification.keys()):
        af.write(f"{species}\n")

print(f"✓ Saved {len(gram_classification)} species to all_identified.txt")

# ----------------------
# Optimization complete - no temporary files to cleanup
# ----------------------
