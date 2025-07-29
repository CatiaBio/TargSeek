#!/usr/bin/env python3
"""
Update the structure mapping file with chain start and end information.
Adds 'chain_start' and 'chain_end' columns after the 'chain' column.
"""

import pandas as pd
import os
from pathlib import Path
import sys


def load_chain_info(structure_id, chains_dir):
    """Load chain information for a specific structure."""
    chain_file = chains_dir / f"{structure_id}_chains.tsv"
    
    if not chain_file.exists():
        print(f"Warning: Chain file not found for {structure_id}: {chain_file}", file=sys.stderr)
        return {}
    
    try:
        df = pd.read_csv(chain_file, sep='\t')
        # Create a dictionary mapping chain -> (start, end)
        chain_dict = {}
        for _, row in df.iterrows():
            # Convert chain to string to handle both string and numeric chains
            chain_key = str(row['chain'])
            chain_dict[chain_key] = {
                'start': row['start'],
                'end': row['end']
            }
        return chain_dict
    except Exception as e:
        print(f"Error reading chain file {chain_file}: {e}", file=sys.stderr)
        return {}


def update_mapping_file(mapping_file, chains_dir, output_file=None):
    """Update the mapping file with chain start and end information."""
    
    if output_file is None:
        output_file = mapping_file
    
    # Load the mapping file
    print(f"Loading mapping file: {mapping_file}")
    df = pd.read_csv(mapping_file, sep='\t')
    
    print(f"Original file has {len(df)} rows and {len(df.columns)} columns")
    print(f"Columns: {list(df.columns)}")
    
    # Add new columns for chain start and end
    df['chain_start'] = None
    df['chain_end'] = None
    
    # Process each row
    processed_structures = set()
    missing_chains = []
    
    for idx, row in df.iterrows():
        structure_id = row['structure_id'].split('_')[0]  # Remove the _1 suffix if present
        chain_id = str(row['chain'])  # Convert to string to handle numeric chains
        
        # Load chain info for this structure (cache to avoid repeated loading)
        if structure_id not in processed_structures:
            print(f"Processing structure {structure_id}...")
            processed_structures.add(structure_id)
        
        chain_info = load_chain_info(structure_id, chains_dir)
        
        if chain_id in chain_info:
            df.at[idx, 'chain_start'] = chain_info[chain_id]['start']
            df.at[idx, 'chain_end'] = chain_info[chain_id]['end']
        else:
            missing_chains.append(f"{structure_id}:{chain_id}")
            print(f"Warning: Chain {chain_id} not found in {structure_id}", file=sys.stderr)
    
    # Report statistics
    filled_rows = df['chain_start'].notna().sum()
    total_rows = len(df)
    
    print(f"\nResults:")
    print(f"- Total rows: {total_rows}")
    print(f"- Rows with chain info added: {filled_rows}")
    print(f"- Missing chain info: {total_rows - filled_rows}")
    
    if missing_chains:
        print(f"\nMissing chains ({len(missing_chains)}):")
        for missing in missing_chains[:10]:  # Show first 10
            print(f"  - {missing}")
        if len(missing_chains) > 10:
            print(f"  ... and {len(missing_chains) - 10} more")
    
    # Reorder columns to place chain_start and chain_end after chain
    columns = list(df.columns)
    chain_idx = columns.index('chain')
    
    # Remove the new columns from their current position
    columns.remove('chain_start')
    columns.remove('chain_end')
    
    # Insert them after the chain column
    columns.insert(chain_idx + 1, 'chain_start')
    columns.insert(chain_idx + 2, 'chain_end')
    
    df = df[columns]
    
    # Save the updated file
    print(f"\nSaving updated file to: {output_file}")
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Updated file has {len(df)} rows and {len(df.columns)} columns")
    print(f"New columns: {list(df.columns)}")


def main():
    """Main function."""
    # Get parameters from snakemake if available, otherwise use defaults
    if 'snakemake' in globals():
        mapping_file = Path(snakemake.input.mapping_file)
        output_file = Path(snakemake.output.final_mapping)
        chains_dir = Path("data/protein_structures/pdb_files")
    else:
        # Fallback for direct execution
        base_dir = Path(".")
        mapping_file = base_dir / "data" / "protein_structures" / "analysis1_params1_fasta_structure_mapping.tsv"
        output_file = base_dir / "data" / "protein_structures" / "analysis1_params1_fasta_structure_mapping_final.tsv"
        chains_dir = base_dir / "data" / "protein_structures" / "pdb_files"
    
    if not mapping_file.exists():
        print(f"Error: Mapping file not found: {mapping_file}", file=sys.stderr)
        sys.exit(1)
    
    if not chains_dir.exists():
        print(f"Error: Chains directory not found: {chains_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Update the mapping file
    update_mapping_file(mapping_file, chains_dir, output_file)
    
    print("\nDone!")


if __name__ == "__main__":
    main()