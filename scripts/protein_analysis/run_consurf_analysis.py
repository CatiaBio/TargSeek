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
            
            # Construct ConSurf command with conda environment and predefined model to avoid ProtTest dependency
            # Run from ConSurf directory and use absolute paths
            consurf_dir = Path("/home/catia/tools/stand_alone_consurf/stand_alone_consurf-1.00")
            conda_cmd = f"cd {consurf_dir} && source /home/catia/tools/miniconda3/etc/profile.d/conda.sh && conda activate targseek && python3 stand_alone_consurf.py --dir {gene_output_dir} --msa {msa_file} --pdb {structure_file} --chain {chain} --query '{query}' --model LG"
            cmd = ['bash', '-c', conda_cmd]
            
            # Run the command
            result = subprocess.run(cmd, capture_output=True, text=True)
            
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
    base_dir = Path('/home/catia/Projects/TargSeek')
    results_dir = base_dir / 'results' / 'analysis1_params1' / 'protein_analysis'
    consurf_script = Path('/home/catia/tools/stand_alone_consurf/stand_alone_consurf-1.00/stand_alone_consurf.py')
    
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
            msa_file = msa_dir / f"{gene}_aligned.fasta"
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

def run_single_consurf_analysis(alignment_file, pdb_file, gene_name, output_dir, consurf_script_path):
    """Run ConSurf analysis for a single gene"""
    try:
        # Determine chain from the alignment file if possible
        chain = "A"  # Default chain
        
        # Create ConSurf command
        cmd = [
            "python3", str(consurf_script_path),
            "--msa", str(alignment_file),
            "--pdb", str(pdb_file), 
            "--chain", chain,
            "--dir", str(output_dir),
            "--query", gene_name
        ]
        
        print(f"Running ConSurf command: {' '.join(cmd)}")
        
        # Run ConSurf
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)  # 30 min timeout
        
        if result.returncode == 0:
            print(f"ConSurf analysis successful for {gene_name}")
            
            # Create a summary file
            summary_file = output_dir / f"{gene_name}_consurf_summary.txt"
            with open(summary_file, 'w') as f:
                f.write(f"ConSurf Analysis Summary for {gene_name}\n")
                f.write(f"================================\n")
                f.write(f"PDB: {pdb_file.name}\n")
                f.write(f"Chain: {chain}\n")
                f.write(f"Alignment: {alignment_file.name}\n")
                f.write(f"Output directory: {output_dir}\n")
                f.write(f"Status: Completed successfully\n\n")
                f.write("STDOUT:\n")
                f.write(result.stdout)
                if result.stderr:
                    f.write("\nSTDERR:\n")
                    f.write(result.stderr)
            
            return True
        else:
            print(f"ConSurf analysis failed for {gene_name}")
            print(f"Return code: {result.returncode}")
            print(f"STDERR: {result.stderr}")
            
            # Create error log
            error_file = output_dir / f"{gene_name}_consurf_error.txt"
            with open(error_file, 'w') as f:
                f.write(f"ConSurf Analysis Error for {gene_name}\n")
                f.write(f"Return code: {result.returncode}\n")
                f.write(f"STDERR: {result.stderr}\n")
                f.write(f"STDOUT: {result.stdout}\n")
                
            return False
            
    except subprocess.TimeoutExpired:
        print(f"ConSurf analysis timed out for {gene_name}")
        return False
    except Exception as e:
        print(f"Error running ConSurf for {gene_name}: {e}")
        return False

def snakemake_main():
    """Main function for Snakemake integration"""
    try:
        # Get inputs from Snakemake
        alignments_with_3d_dir = snakemake.input.alignments_with_3d_dir
        selected_3d_tsv = snakemake.input.selected_3d_tsv
        protein_filter_summary = snakemake.input.protein_filter_summary
        
        # Get outputs
        output_dir = Path(snakemake.output.consurf_results_dir)
        sentinel_file = snakemake.output.consurf_sentinel
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        consurf_script_path = snakemake.params.consurf_script
        
        print(f"=== ConSurf Conservation Analysis ===")
        print(f"Group: {group}")
        print(f"Alignments directory: {alignments_with_3d_dir}")
        print(f"Selected 3D TSV: {selected_3d_tsv}")
        print(f"ConSurf script: {consurf_script_path}")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if ConSurf script exists
        if not Path(consurf_script_path).exists():
            raise FileNotFoundError(f"ConSurf script not found at {consurf_script_path}")
            
        # Run actual ConSurf analysis
        success_count = 0
        total_count = 0
        
        # Read the TSV file to get structure mappings
        import pandas as pd
        try:
            df = pd.read_csv(selected_3d_tsv, sep='\t')
            print(f"Found {len(df)} entries in TSV file")
            
            # Process each gene with structure
            for _, row in df.iterrows():
                gene = row['gene']
                pdb = row['pdb']
                file_name = row.get('file_name', f"{pdb}.fasta")
                
                # Check if we have alignment for this gene
                alignment_file = Path(alignments_with_3d_dir) / f"{gene}_aligned.fasta"
                if not alignment_file.exists():
                    print(f"Warning: No alignment file for {gene}, skipping")
                    continue
                
                # Find corresponding PDB structure using original function
                structure_dirs = [
                    Path("/home/catia/Projects/TargSeek/data/protein_structures/pdb_files"),
                ]
                
                pdb_file = find_structure_file(pdb, structure_dirs)
                if not pdb_file:
                    print(f"Warning: No PDB structure found for {gene} ({pdb}), skipping")
                    continue
                
                # Decompress if needed
                pdb_file = decompress_if_needed(pdb_file)
                
                total_count += 1
                print(f"Running ConSurf analysis for {gene} ({pdb})")
                
                # Create gene-specific output directory
                gene_output_dir = output_dir / gene
                gene_output_dir.mkdir(parents=True, exist_ok=True)
                
                # ConSurf needs the directory to exist before running
                (gene_output_dir / "temp").mkdir(exist_ok=True)
                
                # Determine chain - default to A, but try to get from TSV if available
                chain = row.get('chain', 'A')
                if pd.isna(chain) or chain == '':
                    chain = 'A'
                
                # Use cleaned PDB ID as query (matches the simplified header we created)
                # Remove _1, _2 suffixes but keep names with - (like AF-C5D794)
                import re
                query_id = file_name.replace('.fasta', '') if file_name else pdb
                query_id = re.sub(r'_[0-9]+$', '', query_id)  # Remove _1, _2 etc.
                
                possible_queries = [query_id, query_id.upper(), query_id.lower()]
                print(f"  Expected query ID: {query_id} (cleaned from filename: {file_name})")
                
                # Run ConSurf analysis using original function
                if run_consurf_analysis(alignment_file, pdb_file, chain, possible_queries, output_dir, consurf_script_path):
                    success_count += 1
                    print(f"✓ ConSurf analysis completed for {gene}")
                else:
                    print(f"✗ ConSurf analysis failed for {gene}")
            
        except Exception as e:
            print(f"Error processing TSV file: {e}")
            
        print(f"ConSurf analysis completed: {success_count}/{total_count} analyses successful")
        
        # Create sentinel file
        with open(sentinel_file, 'w') as f:
            f.write(f"ConSurf analysis completed for {group} group\n")
            f.write(f"Analysis: {analysis}\n")
            f.write(f"Paramset: {paramset}\n")
            
        return True
        
    except NameError:
        # Not running under Snakemake, run standalone
        main()
    except Exception as e:
        print(f"Error in ConSurf analysis: {e}")
        return False

if __name__ == "__main__":
    snakemake_main()