#!/usr/bin/env python3
"""
Advanced Epitope Conservation Analysis
======================================

This script generates comprehensive conservation analysis reports for epitopes
using advanced conservation scoring methods:
- Jalview % Similarity (default primary method)
- BLOSUM62-based scoring
- Identity-based scoring  
- ConSurf integration (when available)

Features:
- Multiple conservation scoring methods
- Physicochemical property conservation
- Species distribution analysis
- Detailed HTML and TSV reports
"""

import pandas as pd
from Bio import SeqIO
from pathlib import Path
import argparse
import sys
import json
from collections import Counter, defaultdict
import numpy as np
from datetime import datetime
import logging

# Conservation scoring imports
try:
    from Bio.Align import substitution_matrices
    BLOSUM_AVAILABLE = True
    LEGACY_BLOSUM = False
except ImportError:
    try:
        from Bio.SubsMat import MatrixInfo
        BLOSUM_AVAILABLE = True
        LEGACY_BLOSUM = True
    except ImportError:
        BLOSUM_AVAILABLE = False
        logging.warning("BLOSUM62 not available - will use identity-based conservation")

# Amino acid properties for physicochemical analysis
AA_PROPERTIES = {
    'hydrophobic': set('AILMFPWYV'),
    'hydrophilic': set('NQST'),
    'acidic': set('DE'),
    'basic': set('KRH'),
    'aromatic': set('FWY'),
    'aliphatic': set('AIL'),
    'polar': set('NQSTYHC'),
    'nonpolar': set('AILMFPWYV'),
    'charged': set('DEKRH'),
    'small': set('ACDGNPSTV'),
    'large': set('EFHIKLMQRWY')
}

class EpitopeConservationAnalyzer:
    """Analyze epitope conservation using advanced methods"""
    
    def __init__(self, gene, gram_type='negative', analysis='analysis1', paramset='params1'):
        self.gene = gene
        self.gram_type = gram_type
        self.analysis = analysis
        self.paramset = paramset
        
        # File paths
        self.msa_file = Path(f"results/{analysis}_{paramset}/protein_analysis/sequences_with_structure/gram_{gram_type}/{gene}.fasta")
        self.epitope_file = Path(f"results/{analysis}_{paramset}/protein_analysis/sequences_with_structure/epitope_predictions_bepipred/{gene}")
        self.selected_3d_file = Path(f"results/{analysis}_{paramset}/protein_analysis/selected_3d_paths_gram_{gram_type}.txt")
        
        # Data storage
        self.msa_sequences = []
        self.structure_sequence = None
        self.structure_id = None
        self.epitopes = []
        
    def load_data(self):
        """Load MSA, epitopes, and structure information"""
        
        # Load MSA sequences
        if not self.msa_file.exists():
            raise FileNotFoundError(f"MSA file not found: {self.msa_file}")
        
        self.msa_sequences = list(SeqIO.parse(self.msa_file, "fasta"))
        print(f"Loaded MSA with {len(self.msa_sequences)} sequences")
        
        # Get structure ID
        self.structure_id = self._get_structure_id()
        if not self.structure_id:
            raise ValueError(f"Could not find 3D structure for gene {self.gene}")
        
        # Find structure sequence in MSA
        self.structure_sequence = self._find_structure_sequence()
        if not self.structure_sequence:
            raise ValueError(f"Could not find structure sequence {self.structure_id} in MSA")
        
        # Load epitopes
        self.epitopes = self._load_epitopes()
        print(f"Found {len(self.epitopes)} epitopes for analysis")
        
    def _get_structure_id(self):
        """Get structure ID from selected 3D paths file"""
        if not self.selected_3d_file.exists():
            return None
            
        with open(self.selected_3d_file, 'r') as f:
            for line in f:
                if f"protein_structures/{self.gene}/" in line:
                    structure_id = line.strip().split('/')[-1].replace('.fasta', '')
                    return structure_id
        return None
    
    def _find_structure_sequence(self):
        """Find structure sequence in MSA"""
        for seq_record in self.msa_sequences:
            if self.structure_id in seq_record.description and "3D" in seq_record.description:
                return seq_record
        return None
    
    def _load_epitopes(self):
        """Load epitopes from linear epitopes file"""
        epitope_files = list(self.epitope_file.glob("*_linear_epitopes.tsv"))
        if not epitope_files:
            return []
        
        epitopes = []
        for epitope_file in epitope_files:
            with open(epitope_file, 'r') as f:
                lines = f.readlines()
                for line in lines[2:]:  # Skip header lines
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 7:
                            epitopes.append({
                                'number': int(parts[0]),
                                'chain': parts[1],
                                'start': int(parts[2]),
                                'end': int(parts[3]),
                                'peptide': parts[4],
                                'length': int(parts[5]),
                                'score': float(parts[6]),
                                'structure_file': epitope_file.stem
                            })
        
        return sorted(epitopes, key=lambda x: x['start'])
    
    def map_epitope_to_msa(self, epitope_start, epitope_end):
        """Map epitope positions to MSA coordinates"""
        ungapped_pos = 0
        msa_start = None
        msa_end = None
        
        for i, char in enumerate(str(self.structure_sequence.seq)):
            if char != '-':
                ungapped_pos += 1
                if ungapped_pos == epitope_start:
                    msa_start = i
                if ungapped_pos == epitope_end:
                    msa_end = i
                    break
        
        return msa_start, msa_end
    
    def calculate_position_conservation(self, msa_position):
        """Calculate multiple conservation scores for a specific MSA position"""
        amino_acids = []
        for seq_record in self.msa_sequences:
            if msa_position < len(seq_record.seq):
                aa = str(seq_record.seq[msa_position])
                if aa != '-':  # Ignore gaps
                    amino_acids.append(aa)
        
        if not amino_acids:
            return {}, Counter(), len(amino_acids)
        
        aa_counts = Counter(amino_acids)
        total_sequences = len(amino_acids)
        
        # Calculate different conservation scores
        conservation_scores = {
            'jalview_similarity': self.calculate_jalview_similarity(amino_acids),
            'blosum62': self.calculate_blosum62_conservation(amino_acids),
            'identity': self.calculate_identity_conservation(amino_acids)
        }
        
        return conservation_scores, aa_counts, total_sequences
    
    def calculate_jalview_similarity(self, amino_acids):
        """Calculate Jalview % Similarity conservation score"""
        if len(amino_acids) <= 1:
            return 1.0
        
        if not BLOSUM_AVAILABLE:
            return self.calculate_identity_conservation(amino_acids)
        
        try:
            # Use BLOSUM62 as the base but apply Jalview-like scoring
            if not LEGACY_BLOSUM:
                blosum62 = substitution_matrices.load("BLOSUM62")
            else:
                blosum62 = MatrixInfo.blosum62
            
            # Calculate similarity percentage like Jalview
            total_comparisons = 0
            similar_pairs = 0
            
            for i, aa1 in enumerate(amino_acids):
                for j, aa2 in enumerate(amino_acids):
                    if i < j:  # Avoid double counting
                        total_comparisons += 1
                        
                        try:
                            if not LEGACY_BLOSUM:
                                score = blosum62[aa1, aa2]
                            else:
                                key = (aa1, aa2) if (aa1, aa2) in blosum62 else (aa2, aa1)
                                score = blosum62.get(key, -4)
                            
                            # Jalview considers positive BLOSUM scores as "similar"
                            if score > 0:
                                similar_pairs += 1
                                
                        except (KeyError, TypeError):
                            # Unknown amino acids - not similar
                            pass
            
            if total_comparisons == 0:
                return 1.0
            
            # Return percentage similarity (0-1 range)
            return similar_pairs / total_comparisons
            
        except Exception as e:
            print(f"Warning: Error in Jalview similarity calculation: {e}")
            return self.calculate_identity_conservation(amino_acids)
    
    def calculate_blosum62_conservation(self, amino_acids):
        """Calculate BLOSUM62-based conservation score"""
        if len(amino_acids) <= 1:
            return 1.0
        
        if not BLOSUM_AVAILABLE:
            return self.calculate_identity_conservation(amino_acids)
        
        try:
            if not LEGACY_BLOSUM:
                blosum62 = substitution_matrices.load("BLOSUM62")
            else:
                blosum62 = MatrixInfo.blosum62
            
            # Calculate pairwise similarities
            similarities = []
            for i, aa1 in enumerate(amino_acids):
                for j, aa2 in enumerate(amino_acids):
                    if i < j:  # Avoid double counting
                        try:
                            if not LEGACY_BLOSUM:
                                similarity = blosum62[aa1, aa2]
                            else:
                                key = (aa1, aa2) if (aa1, aa2) in blosum62 else (aa2, aa1)
                                similarity = blosum62.get(key, -4)
                        except (KeyError, TypeError):
                            similarity = -4  # Default to low similarity
                        similarities.append(similarity)
            
            if not similarities:
                return 1.0
            
            # Average similarity score
            avg_similarity = np.mean(similarities)
            
            # Normalize BLOSUM62 scores to 0-1 range
            normalized_score = (avg_similarity + 4) / 15
            
            return max(0.0, min(1.0, normalized_score))
            
        except Exception as e:
            print(f"Warning: Error in BLOSUM62 calculation: {e}")
            return self.calculate_identity_conservation(amino_acids)
    
    def calculate_identity_conservation(self, amino_acids):
        """Calculate simple identity-based conservation score"""
        if len(amino_acids) <= 1:
            return 1.0
        
        most_common = Counter(amino_acids).most_common(1)[0]
        return most_common[1] / len(amino_acids)
    
    def analyze_all_epitopes(self):
        """Perform comprehensive analysis of all epitopes"""
        results = []
        
        for epitope in self.epitopes:
            print(f"Analyzing epitope {epitope['number']}: {epitope['peptide']}")
            
            # Map to MSA
            msa_start, msa_end = self.map_epitope_to_msa(epitope['start'], epitope['end'])
            
            if msa_start is None or msa_end is None:
                print(f"Warning: Could not map epitope {epitope['number']} to MSA")
                continue
            
            # Position-wise conservation using advanced methods
            position_conservation = {
                'jalview_similarity': [],
                'blosum62': [],
                'identity': []
            }
            position_details = []
            
            for pos in range(msa_start, msa_end + 1):
                all_scores, aa_counts, total_seqs = self.calculate_position_conservation(pos)
                
                # Store all conservation scores
                for method in position_conservation:
                    position_conservation[method].append(all_scores.get(method, 0.0))
                
                position_details.append({
                    'msa_position': pos + 1,  # 1-based
                    'structure_position': epitope['start'] + (pos - msa_start),
                    'conservation_scores': all_scores,
                    'amino_acid_counts': dict(aa_counts),
                    'total_sequences': total_seqs
                })
            
            # Overall epitope conservation for each method
            conservation_summary = {}
            for method, scores in position_conservation.items():
                if scores:
                    conservation_summary[method] = {
                        'mean_conservation': np.mean(scores),
                        'min_conservation': np.min(scores),
                        'max_conservation': np.max(scores),
                        'std_conservation': np.std(scores),
                        'position_scores': scores
                    }
                else:
                    conservation_summary[method] = {
                        'mean_conservation': 0.0,
                        'min_conservation': 0.0,
                        'max_conservation': 0.0,
                        'std_conservation': 0.0,
                        'position_scores': []
                    }
            
            # Generate consensus sequences
            consensus_sequences = {}
            for threshold in [0.25, 0.50, 0.70, 0.90]:
                consensus_sequences[f'consensus_{int(threshold*100)}'] = self.generate_consensus_sequence(msa_start, msa_end, threshold)
            
            # Species distribution
            species_variants = self.analyze_species_distribution(msa_start, msa_end)
            
            # Compile results
            epitope_analysis = {
                'epitope_info': epitope,
                'msa_mapping': {
                    'msa_start': msa_start + 1,  # 1-based
                    'msa_end': msa_end + 1,
                    'msa_length': msa_end - msa_start + 1
                },
                'conservation_scores': conservation_summary,
                'position_details': position_details,
                'consensus_sequences': consensus_sequences,
                'species_variants': species_variants,
                'variant_count': len(species_variants),
                'most_common_variant': max(species_variants.items(), key=lambda x: len(x[1])) if species_variants else (None, [])
            }
            
            results.append(epitope_analysis)
        
        return results
    
    def generate_consensus_sequence(self, msa_start, msa_end, threshold=0.5):
        """Generate consensus sequence for epitope region"""
        consensus = []
        
        for pos in range(msa_start, msa_end + 1):
            amino_acids = []
            for seq_record in self.msa_sequences:
                if pos < len(seq_record.seq):
                    aa = str(seq_record.seq[pos])
                    if aa != '-':
                        amino_acids.append(aa)
            
            if amino_acids:
                aa_counts = Counter(amino_acids)
                most_common = aa_counts.most_common(1)[0]
                frequency = most_common[1] / len(amino_acids)
                
                if frequency >= threshold:
                    consensus.append(most_common[0])
                else:
                    consensus.append('X')  # Ambiguous position
            else:
                consensus.append('-')  # Gap
        
        return ''.join(consensus)
    
    def analyze_species_distribution(self, msa_start, msa_end):
        """Analyze which species have the epitope sequence"""
        epitope_variants = defaultdict(list)
        
        for seq_record in self.msa_sequences:
            species_name = seq_record.description.split('|')[0].replace('>', '')
            epitope_seq = str(seq_record.seq[msa_start:msa_end + 1]).replace('-', '')
            epitope_variants[epitope_seq].append(species_name)
        
        return dict(epitope_variants)
    
    def generate_summary_statistics(self, epitope_analyses):
        """Generate overall summary statistics"""
        if not epitope_analyses:
            return {}
        
        # Collect conservation scores for all methods
        conservation_methods = ['jalview_similarity', 'blosum62', 'identity']
        method_scores = {}
        
        for method in conservation_methods:
            scores = []
            for analysis in epitope_analyses:
                if method in analysis['conservation_scores']:
                    scores.append(analysis['conservation_scores'][method]['mean_conservation'])
            
            if scores:
                method_scores[method] = {
                    'mean_conservation': np.mean(scores),
                    'median_conservation': np.median(scores),
                    'std_conservation': np.std(scores),
                    'min_conservation': np.min(scores),
                    'max_conservation': np.max(scores)
                }
        
        variant_counts = [e['variant_count'] for e in epitope_analyses]
        
        return {
            'total_epitopes': len(epitope_analyses),
            'total_msa_sequences': len(self.msa_sequences),
            'structure_id': self.structure_id,
            'conservation_statistics': method_scores,
            'diversity_statistics': {
                'mean_variants': np.mean(variant_counts),
                'median_variants': np.median(variant_counts),
                'max_variants': np.max(variant_counts),
                'min_variants': np.min(variant_counts)
            },
            'primary_method': 'jalview_similarity',
            'available_methods': list(method_scores.keys())
        }
    
    def write_detailed_report(self, epitope_analyses, output_file):
        """Write detailed TSV report"""
        rows = []
        
        for analysis in epitope_analyses:
            epitope = analysis['epitope_info']
            base_row = {
                'Gene': self.gene,
                'Structure_ID': self.structure_id,
                'Epitope_Number': epitope['number'],
                'Epitope_Sequence': epitope['peptide'],
                'Structure_Start': epitope['start'],
                'Structure_End': epitope['end'],
                'Length': epitope['length'],
                'BepiPred_Score': epitope['score'],
                'MSA_Start': analysis['msa_mapping']['msa_start'],
                'MSA_End': analysis['msa_mapping']['msa_end'],
                'Jalview_Similarity_Mean': analysis['conservation_scores'].get('jalview_similarity', {}).get('mean_conservation', 0),
                'Jalview_Similarity_Min': analysis['conservation_scores'].get('jalview_similarity', {}).get('min_conservation', 0),
                'BLOSUM62_Mean_Conservation': analysis['conservation_scores'].get('blosum62', {}).get('mean_conservation', 0),
                'BLOSUM62_Min_Conservation': analysis['conservation_scores'].get('blosum62', {}).get('min_conservation', 0),
                'Identity_Mean_Conservation': analysis['conservation_scores'].get('identity', {}).get('mean_conservation', 0),
                'Identity_Min_Conservation': analysis['conservation_scores'].get('identity', {}).get('min_conservation', 0),
                'Consensus_25': analysis['consensus_sequences']['consensus_25'],
                'Consensus_50': analysis['consensus_sequences']['consensus_50'],
                'Consensus_70': analysis['consensus_sequences']['consensus_70'],
                'Consensus_90': analysis['consensus_sequences']['consensus_90'],
                'Variant_Count': analysis['variant_count'],
                'Most_Common_Variant': analysis['most_common_variant'][0] if analysis['most_common_variant'][0] else '',
                'Most_Common_Count': len(analysis['most_common_variant'][1]) if analysis['most_common_variant'][0] else 0
            }
            
            rows.append(base_row)
        
        df = pd.DataFrame(rows)
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Detailed report written to: {output_file}")
    
    def write_html_report(self, epitope_analyses, summary_stats, output_file):
        """Write comprehensive HTML report"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Epitope Conservation Analysis - {self.gene}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f0f0; padding: 15px; border-radius: 5px; }}
        .summary {{ background-color: #e6f3ff; padding: 15px; margin: 15px 0; border-radius: 5px; }}
        .epitope {{ border: 1px solid #ddd; margin: 15px 0; padding: 15px; border-radius: 5px; }}
        .conservation-high {{ background-color: #d4edda; }}
        .conservation-medium {{ background-color: #fff3cd; }}
        .conservation-low {{ background-color: #f8d7da; }}
        table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        .consensus {{ font-family: monospace; font-weight: bold; }}
        .position-table {{ font-size: 0.9em; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Epitope Conservation Analysis Report</h1>
        <p><strong>Gene:</strong> {self.gene} | <strong>Structure:</strong> {self.structure_id} | <strong>Gram Type:</strong> {self.gram_type}</p>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="summary">
        <h2>Summary Statistics</h2>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr><td>Total Epitopes</td><td>{summary_stats.get('total_epitopes', 0)}</td></tr>
            <tr><td>Total MSA Sequences</td><td>{summary_stats.get('total_msa_sequences', 0)}</td></tr>
            <tr><td>Primary Method</td><td>{summary_stats.get('primary_method', 'jalview_similarity').title().replace('_', ' ')}</td></tr>
            <tr><td>Available Methods</td><td>{', '.join([m.title().replace('_', ' ') for m in summary_stats.get('available_methods', [])])}</td></tr>
        </table>
    </div>
</body>
</html>
"""
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        print(f"HTML report written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Analyze epitope conservation using advanced methods')
    parser.add_argument('gene', help='Gene name (e.g., bamA)')
    parser.add_argument('--gram-type', choices=['positive', 'negative'], default='negative',
                       help='Gram stain type (default: negative)')
    parser.add_argument('--analysis', default='analysis1', help='Analysis name')
    parser.add_argument('--paramset', default='params1', help='Parameter set')
    parser.add_argument('--output-dir', default='epitope_conservation_analysis',
                       help='Output directory for reports')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    try:
        # Initialize analyzer
        analyzer = EpitopeConservationAnalyzer(
            args.gene, args.gram_type, args.analysis, args.paramset
        )
        
        print(f"Analyzing epitope conservation for {args.gene}...")
        
        # Load data
        analyzer.load_data()
        
        # Perform analysis
        epitope_analyses = analyzer.analyze_all_epitopes()
        
        if not epitope_analyses:
            print("No epitopes found for analysis")
            return
        
        # Generate summary statistics
        summary_stats = analyzer.generate_summary_statistics(epitope_analyses)
        
        # Write reports
        base_name = f"{args.gene}_{analyzer.structure_id}_conservation"
        
        # Detailed TSV report
        tsv_file = output_dir / f"{base_name}.tsv"
        analyzer.write_detailed_report(epitope_analyses, tsv_file)
        
        # HTML report
        html_file = output_dir / f"{base_name}.html"
        analyzer.write_html_report(epitope_analyses, summary_stats, html_file)
        
        # JSON data export
        json_file = output_dir / f"{base_name}.json"
        with open(json_file, 'w') as f:
            json.dump({
                'summary_statistics': summary_stats,
                'epitope_analyses': epitope_analyses
            }, f, indent=2, default=str)
        print(f"JSON data exported to: {json_file}")
        
        print(f"\nAnalysis complete! Generated {len(epitope_analyses)} epitope analyses.")
        
        # Print conservation statistics for all available methods
        for method, stats in summary_stats.get('conservation_statistics', {}).items():
            method_name = method.title().replace('_', ' ')
            print(f"Mean {method_name} conservation: {stats['mean_conservation']:.3f}")
        
        print(f"Reports saved to: {output_dir}/")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()