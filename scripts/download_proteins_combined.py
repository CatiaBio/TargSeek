#!/usr/bin/env python3
"""
Combined protein download approach:
1. UniProt iterative (batch + individual queries)
2. NCBI for remaining missing species
3. Update not_found file with final missing species
"""

import argparse
import subprocess
import sys
from pathlib import Path
import logging
import json

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def run_uniprot_iterative(gene: str, species_file: str, output_dir: str) -> dict:
    """
    Run the UniProt iterative download script.
    
    Returns:
        Dictionary with download statistics
    """
    logging.info("=== STEP 1: UniProt Iterative Download ===")
    
    cmd = [
        sys.executable, 
        "scripts/download_proteins_uniprot_iterative.py",
        gene,
        species_file,
        output_dir
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info("UniProt iterative download completed successfully")
        
        # Try to read the summary file
        summary_file = Path(output_dir) / gene / "download_summary.json"
        if summary_file.exists():
            with open(summary_file, 'r') as f:
                return json.load(f)
        else:
            return {}
            
    except subprocess.CalledProcessError as e:
        logging.error(f"UniProt download failed: {e}")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        return {}


def run_ncbi_missing(gene: str, output_dir: str, ncbi_config: str = "config/login/ncbi_info.txt") -> dict:
    """
    Run NCBI download for missing species.
    
    Returns:
        Dictionary with download statistics
    """
    gene_dir = Path(output_dir) / gene
    not_found_file = gene_dir / "not_found_species.txt"
    
    if not not_found_file.exists() or not not_found_file.stat().st_size:
        logging.info("=== STEP 2: NCBI Download - No missing species ===")
        return {"ncbi_downloaded": 0, "still_missing": 0}
    
    # Count missing species before NCBI
    with open(not_found_file, 'r') as f:
        missing_before = len([line.strip() for line in f if line.strip()])
    
    logging.info(f"=== STEP 2: NCBI Download for {missing_before} missing species ===")
    
    cmd = [
        sys.executable,
        "scripts/download_proteins_ncbi_missing.py", 
        gene,
        str(not_found_file),
        str(gene_dir),
        "--ncbi-config", ncbi_config
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info("NCBI download completed successfully")
        
        # Count missing species after NCBI
        missing_after = 0
        if not_found_file.exists():
            with open(not_found_file, 'r') as f:
                missing_after = len([line.strip() for line in f if line.strip()])
        
        ncbi_found = missing_before - missing_after
        
        return {
            "ncbi_downloaded": ncbi_found,
            "still_missing": missing_after
        }
        
    except subprocess.CalledProcessError as e:
        logging.error(f"NCBI download failed: {e}")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        return {"ncbi_downloaded": 0, "still_missing": missing_before}


def generate_final_summary(gene: str, output_dir: str, uniprot_stats: dict, ncbi_stats: dict) -> None:
    """
    Generate and display final download summary.
    """
    gene_dir = Path(output_dir) / gene
    
    # Count total FASTA files
    fasta_files = list(gene_dir.glob("*.fasta"))
    total_downloaded = len(fasta_files)
    
    # Count final missing
    not_found_file = gene_dir / "not_found_species.txt"
    final_missing = 0
    if not_found_file.exists():
        with open(not_found_file, 'r') as f:
            final_missing = len([line.strip() for line in f if line.strip()])
    
    # Calculate totals
    total_requested = uniprot_stats.get('total_species', 0)
    already_existed = uniprot_stats.get('already_existed', 0)
    uniprot_batch = uniprot_stats.get('batch_found', 0)
    uniprot_individual = uniprot_stats.get('individual_found', 0)
    ncbi_downloaded = ncbi_stats.get('ncbi_downloaded', 0)
    
    coverage = (total_downloaded / total_requested * 100) if total_requested > 0 else 0
    
    print(f"\n{'='*50}")
    print(f"FINAL DOWNLOAD SUMMARY FOR {gene.upper()}")
    print(f"{'='*50}")
    print(f"Total species requested: {total_requested}")
    print(f"Already had files: {already_existed}")
    print(f"")
    print(f"UniProt batch query: {uniprot_batch} species")
    print(f"UniProt individual queries: {uniprot_individual} species")
    print(f"NCBI download: {ncbi_downloaded} species")
    print(f"")
    print(f"Total downloaded: {total_downloaded}")
    print(f"Final missing: {final_missing}")
    print(f"Coverage: {coverage:.1f}%")
    
    if final_missing > 0:
        print(f"\nRemaining missing species:")
        with open(not_found_file, 'r') as f:
            missing_species = [line.strip() for line in f if line.strip()]
            for i, species in enumerate(missing_species[:10]):
                print(f"  - {species}")
            if len(missing_species) > 10:
                print(f"  ... and {len(missing_species) - 10} more")
    else:
        print(f"\n🎉 All species found! Complete coverage achieved.")
    
    # Update summary file
    final_stats = {
        "total_species": total_requested,
        "already_existed": already_existed,
        "uniprot_batch_found": uniprot_batch,
        "uniprot_individual_found": uniprot_individual,
        "ncbi_found": ncbi_downloaded,
        "total_downloaded": total_downloaded,
        "final_missing": final_missing,
        "coverage_percent": coverage
    }
    
    summary_file = gene_dir / "final_download_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(final_stats, f, indent=2)


def main():
    parser = argparse.ArgumentParser(description='Combined UniProt + NCBI protein download')
    parser.add_argument('gene', help='Gene name (e.g., bamA, eno)')
    parser.add_argument('species_file', help='File containing species names (one per line)')
    parser.add_argument('output_dir', help='Output directory for FASTA files')
    parser.add_argument('--ncbi-config', default='config/login/ncbi_info.txt',
                       help='Path to NCBI credentials file')
    
    args = parser.parse_args()
    
    logging.info(f"Starting combined download for gene: {args.gene}")
    
    # Step 1: UniProt iterative download
    uniprot_stats = run_uniprot_iterative(args.gene, args.species_file, args.output_dir)
    
    # Step 2: NCBI for missing species
    ncbi_stats = run_ncbi_missing(args.gene, args.output_dir, args.ncbi_config)
    
    # Step 3: Generate final summary
    generate_final_summary(args.gene, args.output_dir, uniprot_stats, ncbi_stats)


def run_snakemake():
    """Run the download with Snakemake inputs"""
    # Get parameters from Snakemake
    protein_lists_dir = snakemake.input.protein_lists
    ncbi_config = snakemake.input.ncbi_info
    output_dir = snakemake.output.download_dir
    
    analysis = snakemake.params.analysis
    paramset = snakemake.params.paramset
    group = snakemake.params.group
    
    logging.info(f"Starting Snakemake download for analysis={analysis}, paramset={paramset}, group={group}")
    logging.info(f"Protein lists directory: {protein_lists_dir}")
    logging.info(f"Output directory: {output_dir}")
    
    # Process all gene files in the protein lists directory
    protein_lists_path = Path(protein_lists_dir)
    if not protein_lists_path.exists():
        logging.error(f"Protein lists directory not found: {protein_lists_dir}")
        return
    
    # Find all .txt files (gene species lists)
    gene_files = list(protein_lists_path.glob("*.txt"))
    logging.info(f"Found {len(gene_files)} gene files to process")
    
    # Process each gene
    for gene_file in gene_files:
        gene = gene_file.stem  # filename without extension
        logging.info(f"Processing gene: {gene}")
        
        # Step 1: UniProt iterative download
        uniprot_stats = run_uniprot_iterative(gene, str(gene_file), output_dir)
        
        # Step 2: NCBI for missing species
        ncbi_stats = run_ncbi_missing(gene, output_dir, ncbi_config)
        
        # Step 3: Generate final summary
        generate_final_summary(gene, output_dir, uniprot_stats, ncbi_stats)


if __name__ == '__main__':
    # Check if running from Snakemake
    if 'snakemake' in globals():
        run_snakemake()
    else:
        main()