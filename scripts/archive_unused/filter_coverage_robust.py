#!/usr/bin/env python3
"""
Robust Coverage Filtering Script
===============================

This script filters gene coverage data by minimum species coverage threshold
and sorts results by coverage count. Includes comprehensive error handling
and cross-platform compatibility.
"""

import pandas as pd
import sys
from pathlib import Path
from typing import Optional

def ensure_output_directory(output_file: str):
    """Ensure output directory exists"""
    output_path = Path(output_file)
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Output directory ready: {output_path.parent}")
    except Exception as e:
        raise ValueError(f"Could not create output directory: {e}")

def load_coverage_data(input_file: str) -> pd.DataFrame:
    """Load coverage data from TSV file"""
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Loaded {len(df)} records from {input_file}")
        
        # Ensure the count column is numeric
        if 'count' in df.columns:
            df['count'] = pd.to_numeric(df['count'], errors='coerce')
        else:
            # Assume third column is count if 'count' column doesn't exist
            df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')
        
        return df
    except Exception as e:
        raise ValueError(f"Could not load coverage data from {input_file}: {e}")

def filter_and_sort_coverage(df: pd.DataFrame, threshold: int) -> pd.DataFrame:
    """Filter coverage data by threshold and sort by count"""
    try:
        # Filter by threshold
        if 'count' in df.columns:
            filtered_df = df[df['count'] > threshold]
            # Sort by count in descending order
            filtered_df = filtered_df.sort_values('count', ascending=False)
        else:
            # Use third column (index 2) as count
            filtered_df = df[df.iloc[:, 2] > threshold]
            # Sort by third column in descending order
            filtered_df = filtered_df.sort_values(df.columns[2], ascending=False)
        
        print(f"Filtered {len(filtered_df)} records above threshold {threshold}")
        return filtered_df
    except Exception as e:
        raise ValueError(f"Could not filter coverage data: {e}")

def save_filtered_data(df: pd.DataFrame, output_file: str):
    """Save filtered data to TSV file"""
    ensure_output_directory(output_file)
    
    try:
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Filtered data saved to {output_file}")
    except Exception as e:
        raise ValueError(f"Could not save filtered data to {output_file}: {e}")

def create_empty_output(input_file: str, output_file: str):
    """Create empty output file with header only"""
    try:
        # Load input to get header
        df = pd.read_csv(input_file, sep='\t')
        header_df = df.iloc[0:0]  # Empty dataframe with headers only
        
        ensure_output_directory(output_file)
        header_df.to_csv(output_file, sep='\t', index=False)
        print(f"Created empty output file with header: {output_file}")
    except Exception as e:
        raise ValueError(f"Could not create empty output file: {e}")

def main():
    """Main function"""
    print("Starting Coverage Filtering...")
    
    # Get inputs from Snakemake or fallback to test values
    try:
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        threshold = snakemake.params.threshold
        wildcards = snakemake.wildcards
        print(f"Processing {input_file} with threshold {threshold}")
    except NameError:
        # Fallback for testing
        print("Snakemake object not available, using test values")
        input_file = "results/coverage/test_coverage.tsv"
        output_file = "results/coverage/test_coverage_filtered.tsv"
        threshold = 25
        wildcards = None
    
    # Load coverage data
    df = load_coverage_data(input_file)
    
    # Filter and sort
    filtered_df = filter_and_sort_coverage(df, threshold)
    
    # Save results
    if len(filtered_df) > 0:
        save_filtered_data(filtered_df, output_file)
        print(f"Filtered {len(filtered_df)} genes with coverage > {threshold}")
    else:
        create_empty_output(input_file, output_file)
        if wildcards:
            print(f"No genes passed the threshold of {threshold} for {wildcards.analysis}_{wildcards.paramset}")
        else:
            print(f"No genes passed the threshold of {threshold}")
        print("Created empty file with header only")

if __name__ == "__main__":
    main()