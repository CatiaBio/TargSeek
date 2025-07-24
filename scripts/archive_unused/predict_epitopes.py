#!/usr/bin/env python3
"""
IEDB Epitope Prediction Integration
===================================

This script uses the IEDB Analysis Resource API to predict B-cell and T-cell epitopes
from conserved protein sequences, integrating with 3D structure and conservation data.
"""

import pandas as pd
import json
import requests
import time
from pathlib import Path
import logging
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
import re
from datetime import datetime
import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Common HLA alleles for population coverage
COMMON_HLA_ALLELES = [
    # HLA Class I (most common worldwide)
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01", "HLA-B*40:01",
    "HLA-C*03:04", "HLA-C*04:01", "HLA-C*06:02", "HLA-C*07:01", "HLA-C*07:02"
]

COMMON_HLA_CLASS_II = [
    # HLA Class II (most common worldwide)
    "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", 
    "HLA-DRB1*11:01", "HLA-DRB1*13:01", "HLA-DRB1*15:01",
    "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*03:01/DQB1*03:02",
    "HLA-DPA1*01:03/DPB1*02:01", "HLA-DPA1*01:03/DPB1*04:01"
]

class IEDBPredictor:
    """Interface for IEDB API predictions"""
    
    def __init__(self, delay=1.0):
        self.delay = delay  # Delay between API calls to be respectful
        self.base_url = "http://tools-cluster-interface.iedb.org/tools_api"
    
    def predict_mhc_class_i(self, sequence, alleles=None, peptide_lengths=[9]):
        """Predict MHC Class I binding epitopes"""
        if alleles is None:
            alleles = COMMON_HLA_ALLELES[:5]  # Use top 5 for speed
        
        results = []
        
        for length in peptide_lengths:
            url = f"{self.base_url}/mhci/"
            data = {
                'method': 'netmhcpan_el',
                'sequence_text': str(sequence),
                'allele': ','.join(alleles),
                'length': str(length)
            }
            
            try:
                response = requests.post(url, data=data, timeout=30)
                response.raise_for_status()
                
                # Parse IEDB response
                lines = response.text.strip().split('\n')
                header_found = False
                
                for line in lines:
                    if line.startswith('allele') or line.startswith('Allele'):
                        header_found = True
                        continue
                    
                    if header_found and line.strip() and not line.startswith('#'):
                        parts = line.split('\t')
                        if len(parts) >= 6:
                            try:
                                result = {
                                    'allele': parts[0],
                                    'peptide': parts[1],
                                    'start': int(parts[2]),
                                    'end': int(parts[3]),
                                    'length': int(parts[4]),
                                    'percentile_rank': float(parts[5]) if parts[5] != 'NA' else None,
                                    'prediction_type': 'MHC_Class_I'
                                }
                                results.append(result)
                            except (ValueError, IndexError):
                                continue
                
                logging.info(f"MHC Class I: Found {len([r for r in results if r['length'] == length])} epitopes (length {length})")
                time.sleep(self.delay)
                
            except Exception as e:
                logging.warning(f"MHC Class I prediction failed for length {length}: {e}")
        
        return results
    
    def predict_mhc_class_ii(self, sequence, alleles=None):
        """Predict MHC Class II binding epitopes"""
        if alleles is None:
            alleles = COMMON_HLA_CLASS_II[:3]  # Use top 3 for speed
        
        url = f"{self.base_url}/mhcii/"
        data = {
            'method': 'NetMHCIIpan',
            'sequence_text': str(sequence),
            'allele': ','.join(alleles)
        }
        
        results = []
        
        try:
            response = requests.post(url, data=data, timeout=60)
            response.raise_for_status()
            
            # Parse IEDB response
            lines = response.text.strip().split('\n')
            header_found = False
            
            for line in lines:
                if line.startswith('allele') or line.startswith('Allele'):
                    header_found = True
                    continue
                
                if header_found and line.strip() and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 6:
                        try:
                            result = {
                                'allele': parts[0],
                                'peptide': parts[1],
                                'start': int(parts[2]),
                                'end': int(parts[3]),
                                'length': int(parts[4]),
                                'percentile_rank': float(parts[5]) if parts[5] != 'NA' else None,
                                'prediction_type': 'MHC_Class_II'
                            }
                            results.append(result)
                        except (ValueError, IndexError):
                            continue
            
            logging.info(f"MHC Class II: Found {len(results)} epitopes")
            time.sleep(self.delay)
            
        except Exception as e:
            logging.warning(f"MHC Class II prediction failed: {e}")
        
        return results
    
    def predict_bcell_epitopes(self, sequence):
        """Predict B-cell epitopes using BepiPred-like approach"""
        # IEDB doesn't have a direct B-cell API, so we'll use a simple approach
        # based on hydrophilicity and surface accessibility predictions
        
        results = []
        min_length = 6
        max_length = 20
        
        # Simple sliding window approach for potential B-cell epitopes
        # This is a basic implementation - in practice, you might want to use
        # more sophisticated methods or local BepiPred installation
        
        sequence_str = str(sequence)
        
        for length in range(min_length, min_length + 5):  # 6-10 amino acids
            for i in range(len(sequence_str) - length + 1):
                peptide = sequence_str[i:i+length]
                
                # Simple scoring based on hydrophilic residues
                hydrophilic_score = sum(1 for aa in peptide if aa in 'DEHKNQRS') / len(peptide)
                
                if hydrophilic_score > 0.4:  # Threshold for potential B-cell epitope
                    result = {
                        'peptide': peptide,
                        'start': i + 1,
                        'end': i + length,
                        'length': length,
                        'hydrophilic_score': hydrophilic_score,
                        'prediction_type': 'B_cell_linear'
                    }
                    results.append(result)
        
        # Remove overlapping epitopes, keep highest scoring
        filtered_results = []
        for result in sorted(results, key=lambda x: x['hydrophilic_score'], reverse=True):
            overlaps = any(
                not (result['end'] <= existing['start'] or result['start'] >= existing['end'])
                for existing in filtered_results
            )
            if not overlaps:
                filtered_results.append(result)
        
        logging.info(f"B-cell epitopes: Found {len(filtered_results)} potential epitopes")
        return filtered_results

def load_conservation_data(conservation_dir, gene):
    """Load conservation scores for a gene"""
    conservation_file = Path(conservation_dir) / f"{gene}_conservation.tsv"
    
    if not conservation_file.exists():
        logging.warning(f"Conservation file not found: {conservation_file}")
        return None
    
    try:
        conservation_df = pd.read_csv(conservation_file, sep='\t')
        return conservation_df
    except Exception as e:
        logging.warning(f"Error loading conservation data for {gene}: {e}")
        return None

def load_3d_structure_info(structures_dir, gene, structures_tsv_data=None):
    """Load 3D structure information, optionally filtered by TSV data"""
    gene_structures_dir = Path(structures_dir) / gene
    
    if not gene_structures_dir.exists():
        return None
    
    # If we have TSV data, only use structures mentioned there
    if structures_tsv_data is not None:
        gene_structures = structures_tsv_data[structures_tsv_data['gene'] == gene]
        if gene_structures.empty:
            return None
        
        # Load only the specific PDB sequences mentioned in TSV
        structure_info = {
            'structures': [],
            'has_3d': True,
            'source': '3D_TSV'
        }
        
        for _, row in gene_structures.iterrows():
            pdb_id = row['pdb_id']
            species = row['species']
            
            # Look for the corresponding FASTA file
            fasta_files = list(gene_structures_dir.glob(f"{pdb_id}*.fasta"))
            
            for fasta_file in fasta_files:
                try:
                    # Read sequence from FASTA
                    with open(fasta_file, 'r') as f:
                        sequences = list(SeqIO.parse(f, "fasta"))
                        if sequences:
                            seq_record = sequences[0]
                            structure_info['structures'].append({
                                'pdb_id': pdb_id,
                                'species': species,
                                'sequence': str(seq_record.seq),
                                'description': seq_record.description,
                                'file': str(fasta_file)
                            })
                except Exception as e:
                    logging.warning(f"Error reading 3D structure file {fasta_file}: {e}")
        
        return structure_info if structure_info['structures'] else None
    
    # Original behavior when no TSV data
    # Look for structure summary or PDB files
    summary_file = gene_structures_dir / "structure_summary.json"
    if summary_file.exists():
        try:
            with open(summary_file, 'r') as f:
                return json.load(f)
        except:
            pass
    
    # Count PDB files as a simple indicator
    pdb_files = list(gene_structures_dir.glob("*.pdb"))
    fasta_files = list(gene_structures_dir.glob("*.fasta"))
    
    if pdb_files or fasta_files:
        return {
            "has_structure": True,
            "pdb_count": len(pdb_files),
            "sequence_count": len(fasta_files)
        }
    
    return None

def score_epitopes_with_conservation(epitopes, conservation_data):
    """Score epitopes based on conservation data"""
    if conservation_data is None or epitopes is None:
        return epitopes
    
    # Create a position-to-conservation mapping
    position_conservation = {}
    for _, row in conservation_data.iterrows():
        position_conservation[row['position']] = row.get('conservation_score', 0)
    
    # Add conservation scores to epitopes
    for epitope in epitopes:
        start_pos = epitope['start']
        end_pos = epitope['end']
        
        # Calculate average conservation for epitope region
        conservation_scores = []
        for pos in range(start_pos, end_pos + 1):
            if pos in position_conservation:
                conservation_scores.append(position_conservation[pos])
        
        if conservation_scores:
            epitope['avg_conservation'] = np.mean(conservation_scores)
            epitope['min_conservation'] = np.min(conservation_scores)
            epitope['conservation_scores'] = conservation_scores
        else:
            epitope['avg_conservation'] = 0
            epitope['min_conservation'] = 0
            epitope['conservation_scores'] = []
    
    return epitopes

def predict_epitopes_for_gene(gene, msa_sequences_dir, conservation_dir, structures_dir, output_dir):
    """Predict epitopes for a single gene"""
    
    logging.info(f"Predicting epitopes for gene: {gene}")
    
    # Load MSA sequences
    gene_fasta = Path(msa_sequences_dir) / f"{gene}.fasta"
    if not gene_fasta.exists():
        logging.warning(f"MSA file not found: {gene_fasta}")
        return None
    
    # Load sequences
    sequences = []
    try:
        for record in SeqIO.parse(gene_fasta, "fasta"):
            sequences.append(record)
    except Exception as e:
        logging.error(f"Error loading sequences for {gene}: {e}")
        return None
    
    if not sequences:
        logging.warning(f"No sequences found for {gene}")
        return None
    
    # Use the first sequence (usually consensus or 3D structure)
    reference_sequence = sequences[0]
    logging.info(f"Using reference sequence: {reference_sequence.id} (length: {len(reference_sequence.seq)})")
    
    # Load conservation data
    conservation_data = load_conservation_data(conservation_dir, gene)
    
    # Load 3D structure info
    structure_info = load_3d_structure_info(structures_dir, gene)
    
    # Initialize IEDB predictor
    predictor = IEDBPredictor()
    
    # Predict epitopes
    all_epitopes = []
    
    # MHC Class I epitopes
    logging.info("Predicting MHC Class I epitopes...")
    mhc_i_epitopes = predictor.predict_mhc_class_i(reference_sequence.seq, peptide_lengths=[9, 10])
    all_epitopes.extend(mhc_i_epitopes)
    
    # MHC Class II epitopes
    logging.info("Predicting MHC Class II epitopes...")
    mhc_ii_epitopes = predictor.predict_mhc_class_ii(reference_sequence.seq)
    all_epitopes.extend(mhc_ii_epitopes)
    
    # B-cell epitopes
    logging.info("Predicting B-cell epitopes...")
    bcell_epitopes = predictor.predict_bcell_epitopes(reference_sequence.seq)
    all_epitopes.extend(bcell_epitopes)
    
    # Add conservation scores
    all_epitopes = score_epitopes_with_conservation(all_epitopes, conservation_data)
    
    # Add gene and sequence information
    for epitope in all_epitopes:
        epitope['gene'] = gene
        epitope['reference_sequence_id'] = reference_sequence.id
        epitope['sequence_length'] = len(reference_sequence.seq)
        epitope['has_3d_structure'] = structure_info is not None
        epitope['molecular_weight'] = molecular_weight(epitope['peptide'], seq_type='protein')
    
    # Create gene output directory
    gene_output_dir = Path(output_dir) / gene
    gene_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save detailed results
    epitopes_df = pd.DataFrame(all_epitopes)
    epitopes_file = gene_output_dir / f"{gene}_epitopes.tsv"
    epitopes_df.to_csv(epitopes_file, sep='\t', index=False)
    
    # Create summary
    summary = {
        'gene': gene,
        'reference_sequence_id': reference_sequence.id,
        'sequence_length': len(reference_sequence.seq),
        'has_conservation_data': conservation_data is not None,
        'has_3d_structure': structure_info is not None,
        'total_epitopes': len(all_epitopes),
        'mhc_class_i_epitopes': len([e for e in all_epitopes if e['prediction_type'] == 'MHC_Class_I']),
        'mhc_class_ii_epitopes': len([e for e in all_epitopes if e['prediction_type'] == 'MHC_Class_II']),
        'bcell_epitopes': len([e for e in all_epitopes if e['prediction_type'] == 'B_cell_linear']),
        'high_conservation_epitopes': len([e for e in all_epitopes if e.get('avg_conservation', 0) > 0.8])
    }
    
    # Save summary
    summary_file = gene_output_dir / f"{gene}_epitope_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logging.info(f"Epitope prediction complete for {gene}:")
    logging.info(f"  Total epitopes: {summary['total_epitopes']}")
    logging.info(f"  MHC Class I: {summary['mhc_class_i_epitopes']}")
    logging.info(f"  MHC Class II: {summary['mhc_class_ii_epitopes']}")
    logging.info(f"  B-cell: {summary['bcell_epitopes']}")
    
    return summary

def create_combined_report(gene_summaries, output_dir):
    """Create a combined epitope prediction report"""
    
    # Create overall summary
    overall_summary = {
        'analysis_date': datetime.now().isoformat(),
        'total_genes': len(gene_summaries),
        'genes_with_epitopes': len([s for s in gene_summaries if s['total_epitopes'] > 0]),
        'total_epitopes_found': sum(s['total_epitopes'] for s in gene_summaries),
        'genes_with_3d_structure': len([s for s in gene_summaries if s['has_3d_structure']]),
        'genes_with_conservation': len([s for s in gene_summaries if s['has_conservation_data']]),
        'per_gene_summary': gene_summaries
    }
    
    # Save JSON summary
    summary_file = Path(output_dir) / "epitope_prediction_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(overall_summary, f, indent=2)
    
    # Create TSV summary table
    summary_df = pd.DataFrame(gene_summaries)
    summary_table = Path(output_dir) / "epitope_prediction_summary.tsv"
    summary_df.to_csv(summary_table, sep='\t', index=False)
    
    # Create human-readable report
    report_file = Path(output_dir) / "epitope_prediction_report.txt"
    with open(report_file, 'w') as f:
        f.write("Epitope Prediction Report\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"OVERVIEW:\n")
        f.write(f"Total genes analyzed: {overall_summary['total_genes']}\n")
        f.write(f"Genes with epitopes: {overall_summary['genes_with_epitopes']}\n")
        f.write(f"Total epitopes found: {overall_summary['total_epitopes_found']}\n")
        f.write(f"Genes with 3D structures: {overall_summary['genes_with_3d_structure']}\n")
        f.write(f"Genes with conservation data: {overall_summary['genes_with_conservation']}\n\n")
        
        f.write("TOP GENES BY EPITOPE COUNT:\n")
        f.write("-" * 30 + "\n")
        top_genes = sorted(gene_summaries, key=lambda x: x['total_epitopes'], reverse=True)[:10]
        for gene_summary in top_genes:
            f.write(f"{gene_summary['gene']}: {gene_summary['total_epitopes']} epitopes ")
            f.write(f"(MHC-I: {gene_summary['mhc_class_i_epitopes']}, ")
            f.write(f"MHC-II: {gene_summary['mhc_class_ii_epitopes']}, ")
            f.write(f"B-cell: {gene_summary['bcell_epitopes']})\n")
        
        f.write(f"\nGENES WITH HIGH CONSERVATION EPITOPES:\n")
        f.write("-" * 30 + "\n")
        high_conservation_genes = [s for s in gene_summaries if s['high_conservation_epitopes'] > 0]
        for gene_summary in sorted(high_conservation_genes, key=lambda x: x['high_conservation_epitopes'], reverse=True):
            f.write(f"{gene_summary['gene']}: {gene_summary['high_conservation_epitopes']} highly conserved epitopes\n")
    
    logging.info(f"Combined epitope report generated:")
    logging.info(f"  Summary: {summary_file}")
    logging.info(f"  Table: {summary_table}")
    logging.info(f"  Report: {report_file}")

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        msa_sequences_dir = snakemake.input.msa_sequences
        conservation_dir = snakemake.input.conservation
        structures_dir = snakemake.input.structures_3d
        proteins_to_study_file = snakemake.input.proteins_to_study
        output_dir = snakemake.output[0]
        
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        group = snakemake.params.group
        
        logging.info(f"IEDB Epitope Prediction for {analysis}_{paramset}_gram_{group}")
        
    except NameError:
        # Test mode
        logging.info("Running in test mode")
        msa_sequences_dir = "results/msa_sequences/analysis_1_params_1_gram_positive"
        conservation_dir = "results/conservation/analysis_1_params_1_gram_positive"
        structures_dir = "results/3d_structures/analysis_1_params_1_gram_positive"
        proteins_to_study_file = "results/analysis_1_params_1/proteins_to_study/gram_positive.tsv"
        output_dir = "results/epitope_predictions/test"
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Load genes to analyze
    try:
        proteins_df = pd.read_csv(proteins_to_study_file, sep='\t')
        genes_to_analyze = proteins_df['gene'].unique().tolist()
        logging.info(f"Found {len(genes_to_analyze)} genes to analyze for epitopes")
    except Exception as e:
        logging.error(f"Error loading proteins to study: {e}")
        return
    
    # Predict epitopes for each gene
    gene_summaries = []
    successful_predictions = 0
    
    for gene in genes_to_analyze:
        try:
            summary = predict_epitopes_for_gene(
                gene, msa_sequences_dir, conservation_dir, structures_dir, output_dir
            )
            if summary:
                gene_summaries.append(summary)
                successful_predictions += 1
        except Exception as e:
            logging.error(f"Error predicting epitopes for {gene}: {e}")
    
    # Generate combined report
    if gene_summaries:
        create_combined_report(gene_summaries, output_dir)
    
    logging.info(f"Epitope prediction complete!")
    logging.info(f"Successfully analyzed {successful_predictions}/{len(genes_to_analyze)} genes")
    logging.info(f"Results saved to: {output_dir}")

if __name__ == "__main__":
    main()