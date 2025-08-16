#!/usr/bin/env python3
"""
Run comprehensive epitope visualizations for all genes in the analysis.
"""

import os
import sys
import argparse
import pandas as pd
import subprocess
from pathlib import Path

def find_integrated_analysis_files(base_dir):
    """Find all integrated epitope analysis files"""

    integrated_files = []

    # Look for integrated analysis files in both gram types
    for gram_type in ['positive', 'negative']:
        gram_dir = os.path.join(base_dir, f"gram_{gram_type}")

        if not os.path.exists(gram_dir):
            continue

        # Find all gene directories
        for gene_dir in os.listdir(gram_dir):
            gene_path = os.path.join(gram_dir, gene_dir)

            if not os.path.isdir(gene_path):
                continue

            # Look for integrated analysis TSV files
            for file in os.listdir(gene_path):
                if file.endswith('_integrated_epitope_analysis.tsv'):
                    # Extract gene and structure_id from filename
                    # Format: {gene}_{structure_id}_integrated_epitope_analysis.tsv
                    base_name = file.replace('_integrated_epitope_analysis.tsv', '')
                    parts = base_name.split('_')

                    if len(parts) >= 2:
                        gene = parts[0]
                        structure_id = '_'.join(parts[1:])  # Handle multi-part structure IDs

                        integrated_files.append({
                            'gene': gene,
                            'structure_id': structure_id,
                            'gram_type': gram_type,
                            'file_path': os.path.join(gene_path, file)
                        })

    return integrated_files

def run_comprehensive_visualizations(analysis, paramset):
    """Run comprehensive visualizations for all genes"""

    # Find all integrated analysis files
    base_dir = f"results/{analysis}_{paramset}/protein_analysis/integrated_epitope_analysis"

    if not os.path.exists(base_dir):
        print(f"Integrated analysis directory not found: {base_dir}")
        return

    integrated_files = find_integrated_analysis_files(base_dir)

    if not integrated_files:
        print("No integrated epitope analysis files found")
        return

    print(f"Found {len(integrated_files)} genes with integrated epitope analysis:")

    success_count = 0
    error_count = 0

    for file_info in integrated_files:
        gene = file_info['gene']
        structure_id = file_info['structure_id']
        gram_type = file_info['gram_type']

        print(f"\nProcessing {gene} ({structure_id}) - gram {gram_type}...")

        try:
            # Run comprehensive visualization script
            cmd = [
                'python',
                'scripts/protein_analysis/create_comprehensive_epitope_visualizations.py',
                '--gene', gene,
                '--structure-id', structure_id,
                '--gram-type', gram_type,
                '--analysis', analysis,
                '--paramset', paramset
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            if result.stdout:
                print(f"  Output: {result.stdout.strip()}")

            success_count += 1
            print(f"  ✅ Completed successfully")

        except subprocess.CalledProcessError as e:
            error_count += 1
            print(f"  ❌ Failed: {e}")
            if e.stdout:
                print(f"    stdout: {e.stdout}")
            if e.stderr:
                print(f"    stderr: {e.stderr}")

        except Exception as e:
            error_count += 1
            print(f"  ❌ Unexpected error: {e}")

    print(f"\n" + "="*60)
    print("COMPREHENSIVE VISUALIZATION SUMMARY")
    print("="*60)
    print(f"Total genes processed: {len(integrated_files)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {error_count}")
    print("="*60)

    if error_count > 0:
        print(f"⚠️  {error_count} genes failed visualization generation")
    else:
        print("🎉 All genes completed successfully!")

    # Create summary file
    summary_file = f"results/{analysis}_{paramset}/protein_analysis/integrated_epitope_analysis/visualization_summary.txt"

    with open(summary_file, 'w') as f:
        f.write(f"Comprehensive Epitope Visualization Summary\n")
        f.write(f"Analysis: {analysis}_{paramset}\n")
        f.write(f"Generated: {pd.Timestamp.now()}\n\n")

        f.write(f"Total genes: {len(integrated_files)}\n")
        f.write(f"Successful: {success_count}\n")
        f.write(f"Failed: {error_count}\n\n")

        f.write("Processed genes:\n")
        for file_info in integrated_files:
            status = "✅" if file_info in integrated_files[:success_count] else "❌"
            f.write(f"{status} {file_info['gene']} ({file_info['structure_id']}) - gram {file_info['gram_type']}\n")

        f.write(f"\nVisualization types created for each epitope:\n")
        f.write(f"• conservation_heatmap.png: ConSurf conservation scores\n")
        f.write(f"• aa_composition.png: Amino acid composition percentages\n")
        f.write(f"• species_count_absolute.png: Absolute species counts\n")
        f.write(f"• species_count_percentage.png: Percentage species presence\n")
        f.write(f"• conservation_vs_presence.png: Rate4Site conservation vs species presence\n")

    print(f"\nSummary saved to: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Run comprehensive epitope visualizations for all genes')
    parser.add_argument('--analysis', required=True, help='Analysis name (e.g., analysis1)')
    parser.add_argument('--paramset', required=True, help='Parameter set (e.g., params1)')

    args = parser.parse_args()

    run_comprehensive_visualizations(args.analysis, args.paramset)

if __name__ == "__main__":
    main()
