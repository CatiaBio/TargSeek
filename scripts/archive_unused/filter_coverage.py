#!/usr/bin/env python3

import pandas as pd
import sys
import os
from pathlib import Path

# Get inputs from Snakemake
input_file = snakemake.input[0]
output_file = snakemake.output[0]
threshold = snakemake.params.threshold

# Use pathlib for better cross-platform compatibility
output_path = Path(output_file)
print(f"DEBUG: Output file path: {output_path}")
print(f"DEBUG: Output directory: {output_path.parent}")

# Ensure output directory exists using pathlib
try:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"DEBUG: Successfully created/confirmed directory using pathlib")
except Exception as e:
    print(f"DEBUG: Error creating directory with pathlib: {e}")
    
print(f"DEBUG: Output directory exists: {output_path.parent.exists()}")

# Read the coverage file
df = pd.read_csv(input_file, sep='\t')

# Ensure the count column is numeric
df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')

# Filter by threshold (count column is the 3rd column, index 2)
filtered_df = df[df.iloc[:, 2] > threshold]

# Check if any rows pass the filter
if len(filtered_df) > 0:
    # Sort by coverage (column 3, index 2) in descending order
    filtered_df = filtered_df.sort_values(by=filtered_df.columns[2], ascending=False)
    
    # Write to output file using pathlib
    try:
        with output_path.open('w', newline='') as f:
            filtered_df.to_csv(f, sep='\t', index=False)
        print(f"DEBUG: Successfully wrote filtered file using pathlib: {output_path}")
    except Exception as e:
        print(f"DEBUG: Error writing filtered file with pathlib: {e}")
        # Fallback to pandas with string path
        try:
            filtered_df.to_csv(str(output_path), sep='\t', index=False)
            print(f"DEBUG: Successfully wrote with string path")
        except Exception as e2:
            print(f"DEBUG: String path also failed: {e2}")
            # Last resort: write to temp file and copy
            import tempfile
            import shutil
            with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_file:
                filtered_df.to_csv(tmp_file.name, sep='\t', index=False)
                shutil.copy2(tmp_file.name, str(output_path))
                os.unlink(tmp_file.name)
            print(f"DEBUG: Successfully wrote via temp file")
    print(f"Filtered {len(filtered_df)} genes with coverage > {threshold}")
else:
    # Create an empty file with just the header
    header_df = df.iloc[0:0]  # Empty dataframe with headers only
    try:
        with output_path.open('w', newline='') as f:
            header_df.to_csv(f, sep='\t', index=False)
        print(f"DEBUG: Successfully wrote empty file using pathlib: {output_path}")
    except Exception as e:
        print(f"DEBUG: Error writing empty file with pathlib: {e}")
        # Fallback to pandas with string path
        try:
            header_df.to_csv(str(output_path), sep='\t', index=False)
            print(f"DEBUG: Successfully wrote empty file with string path")
        except Exception as e2:
            print(f"DEBUG: String path also failed for empty file: {e2}")
            # Last resort: write to temp file and copy
            import tempfile
            import shutil
            with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_file:
                header_df.to_csv(tmp_file.name, sep='\t', index=False)
                shutil.copy2(tmp_file.name, str(output_path))
                os.unlink(tmp_file.name)
            print(f"DEBUG: Successfully wrote empty file via temp file")
    print(f"No genes passed the threshold of {threshold} for {snakemake.wildcards.analysis}_{snakemake.wildcards.paramset}")
    print("Created empty file with header only")