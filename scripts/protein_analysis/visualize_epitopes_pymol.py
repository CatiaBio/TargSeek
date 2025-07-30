#!/usr/bin/env python3
"""
PyMOL 3D Visualization of Predicted Epitopes
============================================

This script creates 3D visualizations of predicted epitopes on protein structures
using PyMOL. It highlights epitope regions with different colors and creates
publication-ready images with legends.

Usage:
    python visualize_epitopes_pymol.py <epitope_tables_dir> <structures_mapping_file>
    
Example:
    python scripts/protein_analysis/visualize_epitopes_pymol.py results/analysis1_params1/protein_analysis/sequences_with_structure/epitope_predictions_bepipred results/analysis1_params1/protein_analysis/sequences_with_structure/epitope_predictions_bepipred/used_structures_mapping.tsv

Requirements:
    - PyMOL (pymol-open-source or PyMOL)
    - pandas
    - biopython
"""

import pandas as pd
import json
import os
import sys
from pathlib import Path
import logging
import argparse
from typing import Dict, List, Tuple

# PyMOL imports - handle both open-source and commercial versions
try:
    import pymol
    from pymol import cmd, stored
    PYMOL_AVAILABLE = True
except ImportError:
    print("Warning: PyMOL not available. Please install pymol-open-source:")
    print("conda install -c conda-forge pymol-open-source")
    PYMOL_AVAILABLE = False

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class EpitopeVisualizer:
    """Create 3D visualizations of epitopes using PyMOL"""
    
    def __init__(self, output_dir: str, pdb_dir: str = "data/protein_structures"):
        """
        Initialize epitope visualizer
        
        Args:
            output_dir: Directory to save visualization outputs
            pdb_dir: Directory containing PDB structure files
        """
        self.output_dir = Path(output_dir)
        self.pdb_dir = Path(pdb_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Color scheme for epitopes (up to 10 different epitopes)
        self.epitope_colors = [
            "red", "blue", "green", "yellow", "magenta",
            "cyan", "orange", "purple", "brown", "pink"
        ]
        
        # Initialize PyMOL
        if PYMOL_AVAILABLE:
            pymol.finish_launching(['pymol', '-cq'])  # -c for command line, -q for quiet
    
    def load_epitope_tables(self, epitope_dir: Path) -> Dict[str, pd.DataFrame]:
        """Load all linear epitope tables from directory"""
        epitope_data = {}
        
        for gene_dir in epitope_dir.iterdir():
            if not gene_dir.is_dir():
                continue
                
            gene_name = gene_dir.name
            
            # Find linear epitopes table
            for epitope_file in gene_dir.glob("*_linear_epitopes.tsv"):
                try:
                    # Extract structure ID from filename (e.g., bamA_5OR1_linear_epitopes.tsv)
                    raw_structure_id = epitope_file.stem.replace(f"{gene_name}_", "").replace("_linear_epitopes", "")
                    
                    # Convert to match structure mapping format
                    # Experimental structures: 5OR1 -> 5OR1_1
                    # Computed structures: AF-O66043 -> AF-O66043 (no change)
                    if raw_structure_id.startswith('AF-'):
                        # AlphaFold computed structure - use as is
                        structure_id = raw_structure_id
                    else:
                        # Experimental structure - add _1 suffix if not present
                        if not raw_structure_id.endswith('_1'):
                            structure_id = f"{raw_structure_id}_1"
                        else:
                            structure_id = raw_structure_id
                    
                    # Read epitope table
                    df = pd.read_csv(epitope_file, sep='\t', comment='#')
                    
                    # Store with combined key
                    key = f"{gene_name}_{structure_id}"
                    epitope_data[key] = {
                        'gene': gene_name,
                        'structure_id': structure_id,
                        'epitopes': df,
                        'epitope_file': epitope_file
                    }
                    
                    logging.info(f"Loaded {len(df)} epitopes for {key}")
                    
                except Exception as e:
                    logging.error(f"Error loading epitope table {epitope_file}: {e}")
        
        return epitope_data
    
    def load_structure_mapping(self, mapping_file: Path) -> Dict[str, Dict]:
        """Load structure mapping file"""
        try:
            df = pd.read_csv(mapping_file, sep='\t')
            mapping = {}
            
            for _, row in df.iterrows():
                gene = row['gene_name']
                structure_id = row['structure_id']
                structure_type = row.get('structure_type', 'unknown')
                chain = row['chain']
                chain_start = row.get('chain_start')
                chain_end = row.get('chain_end')
                fasta_path = row['fasta_path']
                structure_path = row['structure_path']
                
                key = f"{gene}_{structure_id}"
                mapping[key] = {
                    'gene': gene,
                    'structure_id': structure_id,
                    'structure_type': structure_type,
                    'chain': chain,
                    'chain_start': chain_start,
                    'chain_end': chain_end,
                    'fasta_path': fasta_path,
                    'pdb_file': Path(structure_path)  # Use structure_path directly from mapping
                }
            
            logging.info(f"Loaded structure mapping for {len(mapping)} structures")
            return mapping
            
        except Exception as e:
            logging.error(f"Error loading structure mapping: {e}")
            return {}
    
    def _find_pdb_file(self, gene: str, structure_id: str) -> Path:
        """Find the corresponding PDB file for a gene/structure"""
        gene_dir = self.pdb_dir / gene
        
        # Look for PDB files with structure ID (most common first)
        possible_files = [
            gene_dir / f"{structure_id}.pdb.gz",  # Most common format
            gene_dir / f"{structure_id.lower()}.pdb.gz",
            gene_dir / f"{structure_id.upper()}.pdb.gz",
            gene_dir / f"{structure_id}.pdb",
            gene_dir / f"{structure_id.lower()}.pdb",
            gene_dir / f"{structure_id.upper()}.pdb",
        ]
        
        # Also check for other formats
        for ext in ['.cif', '.cif.gz', '.mmcif', '.mmcif.gz']:
            possible_files.extend([
                gene_dir / f"{structure_id}{ext}",
                gene_dir / f"{structure_id.lower()}{ext}",
                gene_dir / f"{structure_id.upper()}{ext}"
            ])
        
        for pdb_file in possible_files:
            if pdb_file.exists():
                return pdb_file
        
        # If not found, log warning but return expected path
        logging.warning(f"PDB file not found for {gene}/{structure_id}")
        return gene_dir / f"{structure_id}.pdb.gz"
    
    def create_visualization(self, gene_structure_key: str, epitope_data: Dict, structure_mapping: Dict) -> bool:
        """Create PyMOL visualization for a specific gene/structure combination"""
        
        if not PYMOL_AVAILABLE:
            logging.error("PyMOL not available - skipping visualization")
            return False
        
        if gene_structure_key not in epitope_data:
            logging.warning(f"No epitope data for {gene_structure_key}")
            return False
        
        if gene_structure_key not in structure_mapping:
            logging.warning(f"No structure mapping for {gene_structure_key}")
            return False
        
        epitopes = epitope_data[gene_structure_key]['epitopes']
        structure_info = structure_mapping[gene_structure_key]
        
        gene = structure_info['gene']
        structure_id = structure_info['structure_id']
        structure_type = structure_info.get('structure_type', 'unknown')
        chain = structure_info['chain']
        chain_start = structure_info.get('chain_start')
        chain_end = structure_info.get('chain_end')
        pdb_file = structure_info['pdb_file']
        
        # Check if the PDB file exists
        if not pdb_file.exists():
            logging.error(f"PDB file not found: {pdb_file}")
            logging.error(f"  Expected path from structure mapping: {pdb_file}")
            return False
        
        try:
            # Clear PyMOL session
            cmd.delete("all")
            cmd.reinitialize()
            
            # Load structure
            cmd.load(str(pdb_file), structure_id)
            cmd.hide("everything", "all")  # Hide everything initially

            # Extract only the desired chain
            if chain and len(chain) > 0:
                chain_selection = f"{structure_id} and chain {chain}"
                cmd.extract("chain_obj", chain_selection)
                cmd.delete(structure_id)  # Delete original multi-chain structure
                structure_id = "chain_obj"  # Replace reference
            else:
                logging.warning(f"No chain information for {gene_structure_key}, using full structure")
                structure_id = structure_id  # fallback (should rarely happen)

            # Show cartoon for the structure (only the selected chain is now loaded)
            cmd.show("cartoon", structure_id)
            cmd.color("gray80", structure_id)

            # Select range if available
            if chain_start is not None and chain_end is not None:
                cmd.select("chain_range", f"{structure_id} and resi {chain_start}-{chain_end}")
                cmd.color("gray80", "chain_range")
                logging.info(f"Focusing on chain {chain} residues {chain_start}-{chain_end}")

            # Create output file names first
            output_prefix = self.output_dir / f"{gene}_{structure_id}"
            png_file = f"{output_prefix}_epitopes.png"
            pse_file = f"{output_prefix}_epitopes.pse"
            
            # Create epitope info for legend
            epitope_info = []
            
            # Epitope coloring - updated to use only chain_obj
            for i, (_, epitope) in enumerate(epitopes.iterrows()):
                if i >= len(self.epitope_colors):
                    break
                start = int(epitope['Start'])
                end = int(epitope['End'])
                color = self.epitope_colors[i]
                selection_name = f"epitope_{i+1}"

                # Adjust if within range
                if chain_start is not None and chain_end is not None:
                    epitope_start = max(start, chain_start)
                    epitope_end = min(end, chain_end)
                    if epitope_start > epitope_end:
                        continue
                    selection = f"{structure_id} and resi {epitope_start}-{epitope_end}"
                else:
                    selection = f"{structure_id} and resi {start}-{end}"

                cmd.select(selection_name, selection)
                cmd.color(color, selection_name)
                
                # Add to epitope info for legend
                epitope_info.append({
                    'number': i + 1,
                    'start': start,
                    'end': end,
                    'peptide': epitope['Peptide'],
                    'score': float(epitope['Score']),
                    'color': color
                })

            # Orientation and export
            cmd.orient(structure_id)
            cmd.zoom(structure_id)
            cmd.set("ray_opaque_background", "off")
            cmd.set("ray_shadows", "on")
            cmd.set("antialias", 2)
            cmd.set("orthoscopic", "on")
            cmd.bg_color("white")

            # Save output files
            cmd.ray(1200, 900)
            cmd.png(str(png_file), dpi=300)
            cmd.save(str(pse_file))
            
            # Create legend/summary file
            self._create_legend(output_prefix, gene, structure_id, epitope_info, chain, chain_start, chain_end)
            
            logging.info(f"✓ Created visualization for {gene}_{structure_id}")
            logging.info(f"  Files: {png_file}, {pse_file}")
            
            return True
            
        except Exception as e:
            logging.error(f"Error creating visualization for {gene_structure_key}: {e}")
            return False
    
    def _create_legend(self, output_prefix: Path, gene: str, structure_id: str, 
                      epitope_info: List[Dict], chain: str, chain_start=None, chain_end=None):
        """Create a legend file with epitope information"""
        
        legend_file = f"{output_prefix}_legend.txt"
        
        with open(legend_file, 'w') as f:
            f.write(f"3D Epitope Visualization Legend\n")
            f.write(f"===============================\n\n")
            f.write(f"Gene: {gene}\n")
            f.write(f"Structure: {structure_id}\n")
            f.write(f"Chain: {chain}\n")
            if chain_start is not None and chain_end is not None:
                f.write(f"Chain Range: {chain_start}-{chain_end}\n")
            f.write(f"Total Epitopes: {len(epitope_info)}\n\n")
            
            f.write("Epitope Details:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'No.':<4} {'Color':<10} {'Start':<6} {'End':<6} {'Length':<7} {'Score':<8} {'Peptide':<20}\n")
            f.write("-" * 80 + "\n")
            
            for epitope in epitope_info:
                f.write(f"{epitope['number']:<4} {epitope['color']:<10} "
                       f"{epitope['start']:<6} {epitope['end']:<6} "
                       f"{len(epitope['peptide']):<7} {epitope['score']:<8.3f} "
                       f"{epitope['peptide']:<20}\n")
            
            f.write("\nVisualization Notes:\n")
            f.write(f"- Visualization focuses on chain {chain}\n")
            if chain_start is not None and chain_end is not None:
                f.write(f"- Only residues {chain_start}-{chain_end} from the chain are shown\n")
            f.write("- Protein backbone is shown in gray cartoon representation\n")
            f.write("- Epitope regions are colored in cartoon representation\n")
            f.write("- Each epitope has a unique color as listed above\n")
            f.write("- Epitopes are only shown if they fall within the chain range\n")
            f.write("- Cartoon coloring shows epitope location on secondary structure\n")
            f.write("- Scores represent BepiPred 3.0 average epitope scores\n")
        
        # Also create JSON format for programmatic access
        json_file = f"{output_prefix}_epitopes.json"
        epitope_json = {
            'gene': gene,
            'structure_id': structure_id,
            'chain': chain,
            'total_epitopes': len(epitope_info),
            'epitopes': epitope_info,
            'files': {
                'image': f"{output_prefix.name}_epitopes.png",
                'session': f"{output_prefix.name}_epitopes.pse",
                'legend': f"{output_prefix.name}_legend.txt"
            }
        }
        
        with open(json_file, 'w') as f:
            json.dump(epitope_json, f, indent=2)
    
    def create_all_visualizations(self, epitope_dir: Path, mapping_file: Path) -> Dict[str, bool]:
        """Create visualizations for all available structures"""
        
        if not PYMOL_AVAILABLE:
            logging.warning("PyMOL not available - creating placeholder files instead")
            return self._create_placeholder_files(epitope_dir, mapping_file)
        
        # Load data
        epitope_data = self.load_epitope_tables(epitope_dir)
        structure_mapping = self.load_structure_mapping(mapping_file)
        
        if not epitope_data:
            logging.error("No epitope data loaded")
            return {}
        
        if not structure_mapping:
            logging.error("No structure mapping loaded")
            return {}
        
        # Create visualizations
        results = {}
        successful = 0
        failed = 0
        
        for gene_structure_key in epitope_data.keys():
            if gene_structure_key in structure_mapping:
                success = self.create_visualization(gene_structure_key, epitope_data, structure_mapping)
                results[gene_structure_key] = success
                
                if success:
                    successful += 1
                else:
                    failed += 1
            else:
                logging.warning(f"No structure mapping for {gene_structure_key}")
                results[gene_structure_key] = False
                failed += 1
        
        # Create summary report
        self._create_summary_report(results, successful, failed)
        
        logging.info(f"\n=== Visualization Summary ===")
        logging.info(f"Successful: {successful}")
        logging.info(f"Failed: {failed}")
        logging.info(f"Success rate: {successful/(successful+failed)*100:.1f}%" if (successful+failed) > 0 else "No structures processed")
        
        return results
    
    def _create_placeholder_files(self, epitope_dir: Path, mapping_file: Path) -> Dict[str, bool]:
        """Create placeholder files when PyMOL is not available"""
        
        # Load data
        epitope_data = self.load_epitope_tables(epitope_dir)
        structure_mapping = self.load_structure_mapping(mapping_file)
        
        if not epitope_data:
            logging.error("No epitope data loaded")
            return {}
        
        results = {}
        
        for gene_structure_key in epitope_data.keys():
            if gene_structure_key in structure_mapping:
                try:
                    epitopes = epitope_data[gene_structure_key]['epitopes']
                    structure_info = structure_mapping[gene_structure_key]
                    
                    gene = structure_info['gene']
                    structure_id = structure_info['structure_id']
                    chain = structure_info['chain']
                    chain_start = structure_info.get('chain_start')
                    chain_end = structure_info.get('chain_end')
                    
                    # Create epitope info
                    epitope_info = []
                    for i, (_, epitope) in enumerate(epitopes.iterrows()):
                        epitope_info.append({
                            'number': i + 1,
                            'start': int(epitope['Start']),
                            'end': int(epitope['End']),
                            'peptide': epitope['Peptide'],
                            'score': float(epitope['Score']),
                            'color': self.epitope_colors[i % len(self.epitope_colors)]
                        })
                    
                    # Create output files
                    output_prefix = self.output_dir / f"{gene}_{structure_id}"
                    
                    # Create placeholder image file
                    placeholder_png = f"{output_prefix}_epitopes.png"
                    with open(placeholder_png, 'w') as f:
                        f.write("# PyMOL visualization placeholder\n")
                        f.write("# PyMOL not available - install with: conda install -c conda-forge pymol-open-source\n")
                    
                    # Create legend and JSON files
                    self._create_legend(output_prefix, gene, structure_id, epitope_info, chain, chain_start, chain_end)
                    
                    results[gene_structure_key] = True
                    logging.info(f"✓ Created placeholder files for {gene}_{structure_id}")
                    
                except Exception as e:
                    logging.error(f"Error creating placeholder for {gene_structure_key}: {e}")
                    results[gene_structure_key] = False
            else:
                results[gene_structure_key] = False
        
        # Create summary report
        successful = sum(1 for success in results.values() if success)
        failed = len(results) - successful
        self._create_summary_report(results, successful, failed)
        
        return results
    
    def _create_summary_report(self, results: Dict[str, bool], successful: int, failed: int):
        """Create a summary report of all visualizations"""
        
        summary_file = self.output_dir / "visualization_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("3D Epitope Visualization Summary Report\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total structures processed: {len(results)}\n")
            f.write(f"Successful visualizations: {successful}\n")
            f.write(f"Failed visualizations: {failed}\n")
            f.write(f"Success rate: {successful/(successful+failed)*100:.1f}%\n\n" if (successful+failed) > 0 else "Success rate: 0%\n\n")
            
            f.write("Detailed Results:\n")
            f.write("-" * 50 + "\n")
            
            for gene_structure, success in sorted(results.items()):
                status = "✓ SUCCESS" if success else "✗ FAILED"
                f.write(f"{gene_structure:<30} {status}\n")
            
            f.write(f"\nOutput files saved to: {self.output_dir}\n")
            f.write("\nFile naming convention:\n")
            f.write("- {gene}_{structure_id}_epitopes.png - 3D visualization image\n")
            f.write("- {gene}_{structure_id}_epitopes.pse - PyMOL session file\n")
            f.write("- {gene}_{structure_id}_legend.txt - Epitope details and legend\n")
            f.write("- {gene}_{structure_id}_epitopes.json - Machine-readable epitope data\n")

def main():
    """Main function for both command line and Snakemake usage"""
    
    # Check if running from Snakemake
    if 'snakemake' in globals():
        # Running from Snakemake
        epitope_tables_sentinel = snakemake.input.epitope_tables_sentinel
        structures_mapping = snakemake.input.structures_mapping
        visualization_sentinel = snakemake.output.visualization_sentinel
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        
        # Derive epitope directory from sentinel file
        epitope_dir = Path(epitope_tables_sentinel).parent
        
        logging.info(f"Running PyMOL epitope visualization for {analysis}_{paramset}")
        logging.info(f"Epitope directory: {epitope_dir}")
        logging.info(f"Structure mapping: {structures_mapping}")
        
        # Create output directory
        output_dir = epitope_dir / "3d_visualizations"
        
        # Create visualizer and run
        visualizer = EpitopeVisualizer(output_dir)
        results = visualizer.create_all_visualizations(epitope_dir, Path(structures_mapping))
        
        # Create sentinel file
        with open(visualization_sentinel, 'w') as f:
            f.write(f"3D epitope visualizations completed for {analysis}_{paramset}\n")
            f.write(f"Processed directory: {epitope_dir}\n")
            f.write(f"Output directory: {output_dir}\n")
            successful = sum(1 for success in results.values() if success)
            failed = len(results) - successful
            f.write(f"Successful: {successful}\n")
            f.write(f"Failed: {failed}\n")
        
        logging.info(f"Sentinel file created: {visualization_sentinel}")
        
    else:
        # Running from command line
        parser = argparse.ArgumentParser(description='Create 3D PyMOL visualizations of predicted epitopes')
        parser.add_argument('epitope_dir', help='Directory containing epitope tables')
        parser.add_argument('mapping_file', help='Structure mapping TSV file')
        parser.add_argument('--output-dir', help='Output directory for visualizations (default: epitope_dir/3d_visualizations)')
        parser.add_argument('--pdb-dir', default='data/protein_structures', help='Directory containing PDB files (default: data/protein_structures)')
        
        args = parser.parse_args()
        
        if not Path(args.epitope_dir).exists():
            logging.error(f"Epitope directory does not exist: {args.epitope_dir}")
            sys.exit(1)
        
        if not Path(args.mapping_file).exists():
            logging.error(f"Mapping file does not exist: {args.mapping_file}")
            sys.exit(1)
        
        # Set output directory
        if args.output_dir:
            output_dir = Path(args.output_dir)
        else:
            output_dir = Path(args.epitope_dir) / "3d_visualizations"
        
        # Create visualizer and run
        visualizer = EpitopeVisualizer(output_dir, args.pdb_dir)
        results = visualizer.create_all_visualizations(Path(args.epitope_dir), Path(args.mapping_file))
        
        print(f"\nVisualization completed. Results saved to: {output_dir}")

if __name__ == "__main__":
    main()