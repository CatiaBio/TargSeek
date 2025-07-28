#!/usr/bin/env python3

"""
Script to generate a summary report of analyzed vs not analyzed genes
Creates a comprehensive report in the BepiPred folder
"""

import pandas as pd
import os
import sys
from pathlib import Path
from datetime import datetime

def load_used_structures_mapping(mapping_file):
    """Load the used structures mapping file"""
    try:
        df = pd.read_csv(mapping_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading used structures mapping: {e}")
        return pd.DataFrame()

def get_all_bepipred_files(bepipred_dir):
    """Get all available BepiPred raw output files"""
    bepipred_files = []
    bepipred_path = Path(bepipred_dir)
    
    for csv_file in bepipred_path.rglob("*_raw_output.csv"):
        # Extract gene and structure info from filename
        filename = csv_file.stem  # Remove .csv extension
        parts = filename.replace('_raw_output', '').split('_')
        
        if len(parts) >= 2:
            gene = parts[0]
            structure_id = '_'.join(parts[1:])
            bepipred_files.append({
                'gene': gene,
                'structure_id': structure_id,
                'file_path': str(csv_file)
            })
    
    return bepipred_files

def get_generated_epitope_files(bepipred_dir):
    """Get all generated epitope TSV files"""
    epitope_files = []
    bepipred_path = Path(bepipred_dir)
    
    for tsv_file in bepipred_path.rglob("*_epitopes.tsv"):
        # Extract gene and structure info from filename
        filename = tsv_file.stem  # Remove .tsv extension
        parts = filename.replace('_epitopes', '').split('_')
        
        if len(parts) >= 2:
            gene = parts[0]
            structure_id = '_'.join(parts[1:])
            epitope_files.append({
                'gene': gene,
                'structure_id': structure_id,
                'file_path': str(tsv_file),
                'folder': tsv_file.parent.name
            })
    
    return epitope_files

def count_epitopes_in_file(file_path):
    """Count the number of epitopes in a TSV file"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            # Skip header lines and count data lines
            epitope_count = 0
            for line in lines[2:]:  # Skip "Predicted Linear Epitope(s):" and header
                if line.strip() and not line.startswith('#'):
                    epitope_count += 1
            return epitope_count
    except Exception as e:
        print(f"Error counting epitopes in {file_path}: {e}")
        return 0

def generate_analysis_summary(bepipred_dir, output_file):
    """Generate comprehensive analysis summary report"""
    
    # Load data
    used_structures_file = os.path.join(bepipred_dir, "used_structures_mapping.tsv")
    used_structures_df = load_used_structures_mapping(used_structures_file)
    
    all_bepipred_files = get_all_bepipred_files(bepipred_dir)
    generated_epitope_files = get_generated_epitope_files(bepipred_dir)
    
    # Create analysis summary
    with open(output_file, 'w') as f:
        f.write("# BepiPred-3.0 Epitope Analysis Summary Report\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Analysis Directory: {bepipred_dir}\n")
        f.write("#\n")
        f.write("# This report summarizes the epitope prediction analysis results\n")
        f.write("# including successfully analyzed genes and any that were not processed\n")
        f.write("\n")
        
        # Section 1: Overall Statistics
        f.write("=" * 80 + "\n")
        f.write("OVERALL ANALYSIS STATISTICS\n")
        f.write("=" * 80 + "\n")
        
        total_selected_structures = len(used_structures_df)
        total_bepipred_files = len(all_bepipred_files)
        total_generated_reports = len(generated_epitope_files)
        
        f.write(f"Total structures selected for analysis: {total_selected_structures}\n")
        f.write(f"Total BepiPred output files available: {total_bepipred_files}\n")
        f.write(f"Total epitope reports generated: {total_generated_reports}\n")
        f.write(f"Success rate: {total_generated_reports/total_selected_structures*100:.1f}%\n")
        f.write("\n")
        
        # Section 2: Successfully Analyzed Genes
        f.write("=" * 80 + "\n")
        f.write("SUCCESSFULLY ANALYZED GENES\n")
        f.write("=" * 80 + "\n")
        f.write("Gene\tStructure_ID\tChain\tEpitopes_Found\tReport_File\n")
        
        analyzed_genes = set()
        total_epitopes = 0
        
        for epitope_file in sorted(generated_epitope_files, key=lambda x: (x['gene'], x['structure_id'])):
            gene = epitope_file['gene']
            structure_id = epitope_file['structure_id']
            
            # Find corresponding chain from used structures
            chain = "Unknown"
            matching_structure = used_structures_df[
                (used_structures_df['Gene'] == gene) & 
                (used_structures_df['Structure_ID'] == structure_id)
            ]
            if not matching_structure.empty:
                chain = matching_structure.iloc[0]['Chain']
            
            # Count epitopes in the file
            epitope_count = count_epitopes_in_file(epitope_file['file_path'])
            total_epitopes += epitope_count
            
            f.write(f"{gene}\t{structure_id}\t{chain}\t{epitope_count}\t{epitope_file['folder']}/{os.path.basename(epitope_file['file_path'])}\n")
            analyzed_genes.add(gene)
        
        f.write(f"\nTotal epitopes found across all genes: {total_epitopes}\n")
        f.write(f"Number of unique genes successfully analyzed: {len(analyzed_genes)}\n")
        f.write("\n")
        
        # Section 3: Not Analyzed Structures
        f.write("=" * 80 + "\n")
        f.write("STRUCTURES NOT ANALYZED\n")
        f.write("=" * 80 + "\n")
        
        not_analyzed = []
        
        for _, row in used_structures_df.iterrows():
            gene = row['Gene']
            structure_id = row['Structure_ID']
            chain = row['Chain']
            
            # Check if epitope report was generated
            found_report = any(
                ef['gene'] == gene and ef['structure_id'] == structure_id 
                for ef in generated_epitope_files
            )
            
            if not found_report:
                # Determine reason
                bepipred_exists = any(
                    bf['gene'] == gene and bf['structure_id'] == structure_id 
                    for bf in all_bepipred_files
                )
                
                if not bepipred_exists:
                    reason = "BepiPred output file not found"
                else:
                    reason = "Structure range mapping issue or no epitopes found"
                
                not_analyzed.append({
                    'gene': gene,
                    'structure_id': structure_id,
                    'chain': chain,
                    'reason': reason
                })
        
        if not_analyzed:
            f.write("Gene\tStructure_ID\tChain\tReason\n")
            for item in sorted(not_analyzed, key=lambda x: (x['gene'], x['structure_id'])):
                f.write(f"{item['gene']}\t{item['structure_id']}\t{item['chain']}\t{item['reason']}\n")
        else:
            f.write("All selected structures were successfully analyzed!\n")
        
        f.write(f"\nNumber of structures not analyzed: {len(not_analyzed)}\n")
        f.write("\n")
        
        # Section 4: Gene Summary
        f.write("=" * 80 + "\n")
        f.write("GENE-LEVEL SUMMARY\n")
        f.write("=" * 80 + "\n")
        
        # Count epitopes by gene
        gene_epitope_counts = {}
        gene_structure_counts = {}
        
        for epitope_file in generated_epitope_files:
            gene = epitope_file['gene']
            epitope_count = count_epitopes_in_file(epitope_file['file_path'])
            
            if gene not in gene_epitope_counts:
                gene_epitope_counts[gene] = 0
                gene_structure_counts[gene] = 0
            
            gene_epitope_counts[gene] += epitope_count
            gene_structure_counts[gene] += 1
        
        f.write("Gene\tStructures_Analyzed\tTotal_Epitopes\tAvg_Epitopes_per_Structure\n")
        
        for gene in sorted(gene_epitope_counts.keys()):
            structures = gene_structure_counts[gene]
            epitopes = gene_epitope_counts[gene]
            avg_epitopes = epitopes / structures if structures > 0 else 0
            
            f.write(f"{gene}\t{structures}\t{epitopes}\t{avg_epitopes:.1f}\n")
        
        f.write("\n")
        
        # Section 5: Files and Directories
        f.write("=" * 80 + "\n")
        f.write("OUTPUT FILES AND DIRECTORIES\n")
        f.write("=" * 80 + "\n")
        f.write("All epitope reports are organized in gene-specific folders:\n")
        
        gene_folders = set(ef['folder'] for ef in generated_epitope_files)
        for folder in sorted(gene_folders):
            gene_files = [ef for ef in generated_epitope_files if ef['folder'] == folder]
            f.write(f"\n{folder}/ ({len(gene_files)} files):\n")
            for ef in sorted(gene_files, key=lambda x: x['structure_id']):
                filename = os.path.basename(ef['file_path'])
                epitope_count = count_epitopes_in_file(ef['file_path'])
                f.write(f"  - {filename} ({epitope_count} epitopes)\n")
        
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")
    
    print(f"Analysis summary report created: {output_file}")
    print(f"Successfully analyzed: {total_generated_reports}/{total_selected_structures} structures")
    print(f"Total epitopes found: {total_epitopes}")
    print(f"Unique genes analyzed: {len(analyzed_genes)}")

def main():
    bepipred_dir = "results/analysis1_params1/protein_analysis/epitope_predictions_bepipred"
    output_file = os.path.join(bepipred_dir, "analysis_summary_report.txt")
    
    if not os.path.exists(bepipred_dir):
        print(f"Error: BepiPred directory not found: {bepipred_dir}")
        sys.exit(1)
    
    print(f"Generating analysis summary report...")
    print(f"BepiPred directory: {bepipred_dir}")
    print(f"Output file: {output_file}")
    print()
    
    generate_analysis_summary(bepipred_dir, output_file)

if __name__ == "__main__":
    main()