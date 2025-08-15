#!/usr/bin/env python3
"""
ConSurf conservation analysis script for Snakemake pipeline.
This script processes MSA files with 3D structures to run ConSurf analysis.

Usage: Called by Snakemake rule run_consurf_conservation_analysis
Input: MSA directory, selected 3D structures TSV, structures download sentinel
Output: ConSurf results directory and completion sentinel
"""

import os
import sys
import subprocess
import pandas as pd
import re
from pathlib import Path
import gzip
import shutil

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

def find_structure_file(pdb_name, structure_dirs):
    """Find the corresponding structure file for a PDB name."""
    for structure_dir in structure_dirs:
        # Try different file extensions and formats
        patterns = [
            f"{pdb_name}.pdb.gz",  # Compressed PDB (preferred)
            f"{pdb_name}.pdb",     # Uncompressed PDB
            f"{pdb_name}.cif.gz",  # Compressed CIF
            f"{pdb_name}.cif",     # Uncompressed CIF
        ]
        
        for pattern in patterns:
            structure_path = structure_dir / pattern
            if structure_path.exists():
                return structure_path
    return None

def decompress_if_needed(structure_file):
    """Decompress .gz file if needed and return path to decompressed file."""
    if str(structure_file).endswith('.gz'):
        decompressed_path = Path(str(structure_file)[:-3])  # Remove .gz extension
        
        if not decompressed_path.exists():
            print(f"  Decompressing {structure_file}...")
            with gzip.open(structure_file, 'rb') as f_in:
                with open(decompressed_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"  Decompressed to {decompressed_path}")
        
        return decompressed_path
    else:
        return structure_file

def prepare_msa_for_consurf(msa_file, pdb_name, output_dir):
    """Create a modified MSA with the structure sequence header renamed to match PDB name."""
    try:
        gene_name = Path(msa_file).stem
        modified_msa = output_dir / f"{gene_name}_modified.fasta"
        
        print(f"  Looking for structure sequence with PDB name '{pdb_name}' in MSA headers...")
        
        with open(msa_file, 'r') as infile, open(modified_msa, 'w') as outfile:
            found_structure_sequence = False
            for line in infile:
                if line.startswith('>') and pdb_name in line:
                    # Found the structure sequence - rename it to just the PDB name
                    print(f"  Found structure sequence: {line.strip()}")
                    print(f"  Renaming to: >{pdb_name}")
                    outfile.write(f">{pdb_name}\n")
                    found_structure_sequence = True
                else:
                    outfile.write(line)
        
        if not found_structure_sequence:
            print(f"  Warning: No sequence containing '{pdb_name}' found in MSA")
            print(f"  Available headers (first 5):")
            with open(msa_file, 'r') as f:
                count = 0
                for line in f:
                    if line.startswith('>'):
                        print(f"    {line.strip()}")
                        count += 1
                        if count >= 5:
                            break
            return None
            
        return modified_msa
        
    except Exception as e:
        print(f"Error preparing MSA for {gene_name}: {e}")
        return None

def run_consurf_analysis(msa_file, structure_file, chain, pdb_name, output_dir, consurf_script):
    """Run ConSurf analysis for a single protein after preparing the MSA."""
    try:
        # Create gene-specific directory
        gene_name = Path(msa_file).stem
        gene_output_dir = output_dir / gene_name
        gene_output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Running ConSurf for {gene_name}...")
        print(f"  Gene: {gene_name}")
        print(f"  PDB: {pdb_name}")
        print(f"  Structure: {structure_file}")
        print(f"  Chain: {chain}")
        print(f"  Gene output directory: {gene_output_dir}")
        
        # Prepare MSA with renamed header in the gene directory
        modified_msa = prepare_msa_for_consurf(msa_file, pdb_name, gene_output_dir)
        if not modified_msa:
            print(f"  ✗ Failed to prepare MSA for {gene_name}")
            return False
        
        # Use the PDB name as query (this should now match the renamed header)
        query = pdb_name
        print(f"  Using query: '{query}' with modified MSA")
        
        # Run ConSurf from the project root directory (not gene directory)
        cmd = [
            'python3', str(consurf_script),
            '--dir', str(gene_output_dir),
            '--msa', str(modified_msa),
            '--structure', str(structure_file),
            '--chain', chain,
            '--query', query
        ]
        
        print(f"  ConSurf command: {' '.join(cmd)}")
        print(f"  Working directory: {Path.cwd()}")
        
        # Run the command from the project root directory
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())
        
        if result.returncode == 0 and "The query sequence is not in the msa" not in result.stdout:
            print(f"  ✓ ConSurf analysis completed for {gene_name}")
            
            # Create success marker file in the gene directory
            success_file = gene_output_dir / "consurf_success.txt"
            with open(success_file, 'w') as f:
                f.write(f"ConSurf analysis completed successfully\n")
                f.write(f"Gene: {gene_name}\n")
                f.write(f"PDB structure name: {pdb_name}\n")
                f.write(f"Structure file: {structure_file}\n")
                f.write(f"Chain: {chain}\n")
                f.write(f"Query used: {query}\n")
                f.write(f"Modified MSA: {modified_msa}\n")
                f.write(f"ConSurf results in: {gene_output_dir}\n")
            
            return True
        else:
            print(f"  ✗ ConSurf failed for {gene_name}")
            print(f"  ✗ stdout: {result.stdout.strip()}")
            if result.stderr:
                print(f"  ✗ stderr: {result.stderr.strip()}")
            return False
            
    except Exception as e:
        print(f"Error running ConSurf for {gene_name}: {e}")
        return False

def main():
    # Get parameters from Snakemake
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    group = snakemake.params.group
    consurf_script = Path(snakemake.params.consurf_script)
    
    # Input files
    alignments_dir = Path(snakemake.input.alignments_with_3d_dir)
    selected_3d_tsv = Path(snakemake.input.selected_3d_tsv)
    
    # Output files
    output_dir = Path(snakemake.output.consurf_results_dir)
    sentinel_file = Path(snakemake.output.consurf_sentinel)
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Structure directories to search
    base_dir = Path.cwd()
    structure_dirs = [
        base_dir / 'data' / 'protein_structures' / 'pdb_files',
        base_dir / 'data' / 'protein_structures' / 'structures'
    ]
    
    # Check if ConSurf script exists
    if not consurf_script.exists():
        print(f"Error: ConSurf script not found at {consurf_script}")
        sys.exit(1)
    
    print(f"=== Processing {group} proteins for {analysis}_{paramset} ===")
    
    # Check if TSV file exists
    if not selected_3d_tsv.exists():
        print(f"Warning: TSV file not found: {selected_3d_tsv}")
        # Create empty sentinel and exit
        sentinel_file.parent.mkdir(parents=True, exist_ok=True)
        with open(sentinel_file, 'w') as f:
            f.write("No 3D structures found for analysis\n")
        return
    
    # Check if MSA directory exists
    if not alignments_dir.exists():
        print(f"Warning: MSA directory not found: {alignments_dir}")
        # Create empty sentinel and exit
        sentinel_file.parent.mkdir(parents=True, exist_ok=True)
        with open(sentinel_file, 'w') as f:
            f.write("No MSA directory found\n")
        return
    
    # Extract PDB data from TSV
    pdb_data = extract_pdb_id_from_tsv(selected_3d_tsv)
    if not pdb_data:
        print(f"Warning: No PDB data found in {selected_3d_tsv}")
        # Create empty sentinel and exit
        sentinel_file.parent.mkdir(parents=True, exist_ok=True)
        with open(sentinel_file, 'w') as f:
            f.write("No PDB data found in TSV file\n")
        return
    
    # Process each gene: 
    # 1. Get MSA file from alignments_with_3d_dir
    # 2. Look up PDB structure name in selected_3d_tsv
    # 3. Find PDB sequence in MSA headers and rename to structure name
    # 4. Find structure file in data/protein_structures/pdb_files/
    
    success_count = 0
    total_count = 0
    results_summary = []
    
    print(f"Found {len(pdb_data)} genes with 3D structures in TSV")
    
    for gene, info in pdb_data.items():
        pdb_name = info['pdb']  # This is the structure name from TSV
        chain = info['chain'] if info['chain'] else 'A'
        
        print(f"\n--- Processing gene: {gene} ---")
        print(f"Selected PDB structure: {pdb_name}")
        print(f"Chain: {chain}")
        
        # 1. Find MSA file for this gene
        msa_file = alignments_dir / f"{gene}.fasta"
        if not msa_file.exists():
            print(f"✗ MSA file not found: {msa_file}")
            results_summary.append(f"{gene}: MSA file not found")
            continue
        
        print(f"✓ Found MSA file: {msa_file}")
        
        # 2. Find structure file in data/protein_structures/pdb_files/
        structure_file = find_structure_file(pdb_name, structure_dirs)
        if not structure_file:
            print(f"✗ Structure file not found for PDB {pdb_name}")
            print(f"  Searched in: {structure_dirs}")
            results_summary.append(f"{gene}: Structure file not found for {pdb_name}")
            continue
        
        print(f"✓ Found structure file: {structure_file}")
        
        # 3. Decompress structure if needed (ConSurf prefers uncompressed)
        structure_file = decompress_if_needed(structure_file)
        
        # 4. Run ConSurf analysis
        total_count += 1
        if run_consurf_analysis(msa_file, structure_file, chain, pdb_name, output_dir, consurf_script):
            success_count += 1
            results_summary.append(f"{gene}: ConSurf analysis completed successfully")
        else:
            results_summary.append(f"{gene}: ConSurf analysis failed")
    
    print(f"\n{group} summary: {success_count}/{total_count} analyses completed successfully")
    
    # Create sentinel file with summary
    sentinel_file.parent.mkdir(parents=True, exist_ok=True)
    with open(sentinel_file, 'w') as f:
        f.write(f"ConSurf analysis completed for {group} group\n")
        f.write(f"Analysis: {analysis}_{paramset}\n")
        f.write(f"Success rate: {success_count}/{total_count}\n")
        f.write("\nDetailed results:\n")
        for result in results_summary:
            f.write(f"  {result}\n")
    
    print(f"ConSurf analysis completed for {group}!")
    print(f"Results saved in: {output_dir}")
    print(f"Summary saved in: {sentinel_file}")

if __name__ == "__main__":
    main()