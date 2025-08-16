#!/usr/bin/env python3
"""
Batch script to run integrated epitope analysis for all proteins with complete data from:
1. BepiPred epitope predictions
2. ConSurf conservation analysis
3. DeepTMHMM topology predictions

This script identifies all genes that have data from all three tools and runs the
integrated analysis for each.
"""

import os
import sys
import glob
import subprocess
import json
from pathlib import Path
import argparse

def find_genes_with_complete_data(analysis, paramset, base_dir="results"):
    """
    Find genes that have data from all three analysis tools.
    
    Returns:
        Dictionary with gram_type as keys and list of (gene, structure_id) tuples as values
    """
    results_dir = f"{base_dir}/{analysis}_{paramset}/protein_analysis"
    
    genes_with_complete_data = {"positive": [], "negative": []}
    
    for gram_type in ["positive", "negative"]:
        print(f"\nScanning for complete data in gram_{gram_type}...")
        
        # Find ConSurf results
        consurf_dir = f"{results_dir}/consurf_analysis/gram_{gram_type}"
        if not os.path.exists(consurf_dir):
            print(f"Warning: ConSurf directory not found: {consurf_dir}")
            continue
            
        consurf_genes = set()
        for gene_dir in glob.glob(f"{consurf_dir}/*/"):
            gene = os.path.basename(gene_dir.rstrip('/'))
            conservation_file = f"{gene_dir}/msa_aa_variety_percentage.csv"
            if os.path.exists(conservation_file):
                consurf_genes.add(gene)
        
        print(f"  ConSurf genes: {len(consurf_genes)}")
        
        # Find Topology results
        topology_dir = f"{results_dir}/topology_analysis/gram_{gram_type}"
        if not os.path.exists(topology_dir):
            print(f"Warning: Topology directory not found: {topology_dir}")
            continue
            
        topology_genes = set()
        for gene_dir in glob.glob(f"{topology_dir}/*/"):
            gene = os.path.basename(gene_dir.rstrip('/'))
            topology_file = f"{gene_dir}/biolib_results/deeptmhmm_results.md"
            if os.path.exists(topology_file):
                topology_genes.add(gene)
        
        print(f"  Topology genes: {len(topology_genes)}")
        
        # Find BepiPred results
        bepipred_dir = f"{results_dir}/bepipred_epitope_predictions"
        if not os.path.exists(bepipred_dir):
            print(f"Warning: BepiPred directory not found: {bepipred_dir}")
            continue
            
        bepipred_genes = {}
        for gene_dir in glob.glob(f"{bepipred_dir}/*/"):
            gene = os.path.basename(gene_dir.rstrip('/'))
            # Find epitope files for this gene
            epitope_files = glob.glob(f"{gene_dir}/*_linear_epitopes.tsv")
            if epitope_files:
                # Extract structure ID from filename
                for epitope_file in epitope_files:
                    filename = os.path.basename(epitope_file)
                    # Format: {gene}_{structure_id}_linear_epitopes.tsv
                    parts = filename.replace('_linear_epitopes.tsv', '').split('_')
                    if len(parts) >= 2:
                        structure_id = '_'.join(parts[1:])  # Handle structure IDs with underscores
                        bepipred_genes[gene] = structure_id
        
        print(f"  BepiPred genes: {len(bepipred_genes)}")
        
        # Find intersection of all three
        common_genes = consurf_genes & topology_genes & set(bepipred_genes.keys())
        
        print(f"  Common genes (all 3 tools): {len(common_genes)}")
        
        for gene in common_genes:
            structure_id = bepipred_genes[gene]
            genes_with_complete_data[gram_type].append((gene, structure_id))
            print(f"    {gene} ({structure_id})")
    
    return genes_with_complete_data

def run_integrated_analysis(gene, structure_id, gram_type, analysis, paramset, output_dir):
    """
    Run integrated analysis for a single gene.
    """
    print(f"\nRunning integrated analysis for {gene} ({structure_id}) - gram_{gram_type}")
    
    cmd = [
        "python", "scripts/protein_analysis/integrate_epitope_analysis.py",
        "--gene", gene,
        "--structure-id", structure_id,
        "--gram-type", gram_type,
        "--analysis", analysis,
        "--paramset", paramset,
        "--output-dir", f"{output_dir}/gram_{gram_type}/{gene}"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"  ✓ Success: {gene}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Failed: {gene}")
        print(f"    Error: {e.stderr}")
        return False

def create_summary_report(genes_with_complete_data, analysis, paramset, output_dir):
    """
    Create a summary report of all integrated analyses.
    """
    from datetime import datetime
    summary = {
        "analysis": analysis,
        "paramset": paramset,
        "timestamp": str(datetime.now()),
        "summary": {
            "total_genes_positive": len(genes_with_complete_data["positive"]),
            "total_genes_negative": len(genes_with_complete_data["negative"]),
            "genes_by_gram_type": {}
        }
    }
    
    for gram_type in ["positive", "negative"]:
        summary["summary"]["genes_by_gram_type"][f"gram_{gram_type}"] = [
            {"gene": gene, "structure_id": structure_id}
            for gene, structure_id in genes_with_complete_data[gram_type]
        ]
    
    # Save summary
    summary_file = f"{output_dir}/integrated_analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nSummary report saved to: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Run integrated epitope analysis for all genes with complete data')
    parser.add_argument('--analysis', default='analysis1', help='Analysis name')
    parser.add_argument('--paramset', default='params1', help='Parameter set')
    parser.add_argument('--output-dir', help='Output directory (optional)')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be processed without running')
    
    args = parser.parse_args()
    
    # Set output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = f"results/{args.analysis}_{args.paramset}/protein_analysis/integrated_epitope_analysis"
    
    print(f"TargSeek Integrated Epitope Analysis")
    print(f"Analysis: {args.analysis}")
    print(f"Parameter set: {args.paramset}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find genes with complete data
    genes_with_complete_data = find_genes_with_complete_data(args.analysis, args.paramset)
    
    total_genes = sum(len(genes) for genes in genes_with_complete_data.values())
    print(f"\nTotal genes with complete data: {total_genes}")
    
    if total_genes == 0:
        print("No genes found with complete data from all three tools.")
        sys.exit(1)
    
    if args.dry_run:
        print("\n=== DRY RUN - Would process the following genes ===")
        for gram_type, genes in genes_with_complete_data.items():
            if genes:
                print(f"\nGram {gram_type}:")
                for gene, structure_id in genes:
                    print(f"  {gene} ({structure_id})")
        sys.exit(0)
    
    # Run integrated analysis for each gene
    successful = 0
    failed = 0
    
    for gram_type, genes in genes_with_complete_data.items():
        for gene, structure_id in genes:
            success = run_integrated_analysis(
                gene, structure_id, gram_type, 
                args.analysis, args.paramset, output_dir
            )
            if success:
                successful += 1
            else:
                failed += 1
    
    print(f"\n=== RESULTS ===")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Total: {successful + failed}")
    
    # Create summary report
    create_summary_report(genes_with_complete_data, args.analysis, args.paramset, output_dir)
    
    # Create sentinel file
    sentinel_file = f"{output_dir}/integrated_analysis_complete.sentinel"
    with open(sentinel_file, 'w') as f:
        f.write(f"Integrated epitope analysis completed\n")
        f.write(f"Analysis: {args.analysis}\n")
        f.write(f"Parameter set: {args.paramset}\n")
        f.write(f"Successful analyses: {successful}\n")
        f.write(f"Failed analyses: {failed}\n")
    
    print(f"Sentinel file created: {sentinel_file}")
    
    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    main()