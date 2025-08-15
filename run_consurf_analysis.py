#!/usr/bin/env python3
"""
Automated ConSurf analysis script for protein MSAs and structures.
This script processes all MSA files and their corresponding 3D structures to run ConSurf analysis.
"""

import os
import sys
import subprocess
import pandas as pd
import re
from pathlib import Path

def extract_pdb_id_from_tsv(tsv_file):
    """Extract PDB IDs and associated data from the selected 3D structures TSV file."""
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        pdb_data = {}
        for _, row in df.iterrows():
            gene = row['gene']
            pdb = row['pdb']
            # Handle NaN chains - extract from bacteria column or default to A
            chain = row.get('chain', 'A')
            if pd.isna(chain) or chain == '':
                # Try to extract chain from bacteria column
                bacteria_info = str(row.get('bacteria', ''))
                if 'Chain ' in bacteria_info:
                    chain_match = re.search(r'Chain ([A-Z])', bacteria_info)
                    if chain_match:
                        chain = chain_match.group(1)
                    else:
                        chain = 'A'
                else:
                    chain = 'A'
            pdb_data[gene] = {'pdb': pdb, 'chain': chain}
        return pdb_data
    except Exception as e:
        print(f"Error reading TSV file {tsv_file}: {e}")
        return {}

def extract_query_from_msa(msa_file, pdb_id):
    """Extract the query sequence identifier from MSA file based on PDB ID."""
    try:
        with open(msa_file, 'r') as f:
            for line in f:
                if line.startswith('>') and pdb_id in line:
                    # Extract the full header line without the '>'
                    query = line.strip()[1:]
                    # Try different query formats that ConSurf might accept
                    possible_queries = [
                        query,  # Full header
                        pdb_id + '_1',  # Simple PDB format
                        pdb_id,  # Just PDB ID
                        query.split('|')[0] if '|' in query else query,  # First part before |
                    ]
                    return possible_queries
        return None
    except Exception as e:
        print(f"Error reading MSA file {msa_file}: {e}")
        return None

def find_structure_file(pdb_id, structure_dirs):
    """Find the corresponding structure file for a PDB ID."""
    for structure_dir in structure_dirs:
        # Try different naming patterns including compressed files
        patterns = [
            f"{pdb_id}.pdb",
            f"{pdb_id}.pdb.gz",
            f"{pdb_id}_1.pdb", 
            f"{pdb_id}_1.pdb.gz",
            f"AF-{pdb_id}.pdb",
            f"AF-{pdb_id}.pdb.gz"
        ]
        
        for pattern in patterns:
            structure_path = structure_dir / pattern
            if structure_path.exists():
                return structure_path
    return None

def decompress_if_needed(structure_file):
    """Decompress .gz file if needed and return path to decompressed file."""
    if str(structure_file).endswith('.gz'):
        import gzip
        import shutil
        
        decompressed_path = Path(str(structure_file)[:-3])  # Remove .gz extension
        
        if not decompressed_path.exists():
            print(f"Decompressing {structure_file}...")
            with gzip.open(structure_file, 'rb') as f_in:
                with open(decompressed_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        return decompressed_path
    else:
        return structure_file

def run_consurf_analysis(msa_file, structure_file, chain, possible_queries, output_dir, consurf_script):
    """Run ConSurf analysis for a single protein, trying different query formats."""
    try:
        # Create output directory for this analysis
        gene_name = Path(msa_file).stem
        gene_output_dir = output_dir / gene_name
        gene_output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Running ConSurf for {gene_name}...")
        
        # Try different query formats
        for i, query in enumerate(possible_queries):
            print(f"  Attempt {i+1}: Trying query format '{query}'")
            
            # Construct ConSurf command
            cmd = [
                'python3', str(consurf_script),
                '--dir', str(gene_output_dir),
                '--msa', str(msa_file),
                '--structure', str(structure_file),
                '--chain', chain,
                '--query', query
            ]
            
            # Run the command
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=gene_output_dir)
            
            if result.returncode == 0 and "The query sequence is not in the msa" not in result.stdout:
                print(f"  ✓ ConSurf analysis completed for {gene_name} with query '{query}'")
                return True
            else:
                print(f"  ✗ Query '{query}' failed: {result.stdout.strip()}")
        
        print(f"✗ All query formats failed for {gene_name}")
        return False
            
    except Exception as e:
        print(f"Error running ConSurf for {gene_name}: {e}")
        return False

def main():
    # Configuration
    base_dir = Path('/mnt/c/Users/catia/Projects/PureMilk')
    results_dir = base_dir / 'results' / 'analysis1_params1' / 'protein_analysis'
    consurf_script = base_dir / 'tools' / 'stand_alone_consurf' / 'stand_alone_consurf.py'
    
    # Structure directories to search
    structure_dirs = [
        base_dir / 'data' / 'protein_structures' / 'pdb_files',
        base_dir / 'data' / 'protein_structures' / 'structures'
    ]
    
    # Output directory for ConSurf results
    output_dir = base_dir / 'results' / 'analysis1_params1' / 'consurf_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if ConSurf script exists
    if not consurf_script.exists():
        print(f"Error: ConSurf script not found at {consurf_script}")
        sys.exit(1)
    
    # Process both gram-positive and gram-negative
    for gram_type in ['gram_negative', 'gram_positive']:
        print(f"\n=== Processing {gram_type} proteins ===")
        
        # Paths for this gram type
        tsv_file = results_dir / f'selected_3d_fasta_{gram_type}.tsv'
        msa_dir = results_dir / 'mafft_alignment_with_3d' / gram_type
        
        if not tsv_file.exists():
            print(f"Warning: TSV file not found: {tsv_file}")
            continue
            
        if not msa_dir.exists():
            print(f"Warning: MSA directory not found: {msa_dir}")
            continue
        
        # Extract PDB data from TSV
        pdb_data = extract_pdb_id_from_tsv(tsv_file)
        if not pdb_data:
            print(f"Warning: No PDB data found in {tsv_file}")
            continue
        
        # Process each gene
        success_count = 0
        total_count = 0
        
        for gene, info in pdb_data.items():
            pdb_id = info['pdb']
            chain = info['chain'] if info['chain'] else 'A'
            
            # Find MSA file
            msa_file = msa_dir / f"{gene}.fasta"
            if not msa_file.exists():
                print(f"Warning: MSA file not found: {msa_file}")
                continue
            
            # Extract query from MSA
            possible_queries = extract_query_from_msa(msa_file, pdb_id)
            if not possible_queries:
                print(f"Warning: Could not extract query for {gene} with PDB {pdb_id}")
                continue
            
            # Find structure file
            structure_file = find_structure_file(pdb_id, structure_dirs)
            if not structure_file:
                print(f"Warning: Structure file not found for PDB {pdb_id}")
                continue
            
            # Decompress if needed
            structure_file = decompress_if_needed(structure_file)
            
            # Create specific output directory
            gene_gram_output_dir = output_dir / gram_type / gene
            
            total_count += 1
            if run_consurf_analysis(msa_file, structure_file, chain, possible_queries, gene_gram_output_dir, consurf_script):
                success_count += 1
        
        print(f"\n{gram_type} summary: {success_count}/{total_count} analyses completed successfully")
    
    print(f"\nAll ConSurf analyses completed!")
    print(f"Results are saved in: {output_dir}")

if __name__ == "__main__":
    main()