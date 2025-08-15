#!/usr/bin/env python3
"""
Epitope Conservation Analysis Report Generator
==============================================

This script generates comprehensive conservation analysis reports for epitopes
by analyzing their positions within multiple sequence alignments (MSAs).

Features:
- Conservation scores for each position
- Consensus sequence analysis
- Physicochemical property conservation
- Species distribution analysis
- Epitope ranking by conservation
- Detailed HTML and TSV reports

Usage:
    python analyze_epitope_conservation.py <gene> [options]
    
Example:
    python scripts/protein_analysis/analyze_epitope_conservation.py bamA --output-dir epitope_analysis
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

# Try to import AAindex if available
try:
    import requests
    AAINDEX_AVAILABLE = True
except ImportError:
    AAINDEX_AVAILABLE = False
    logging.warning("requests not available - AAindex scoring disabled")

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
    """Analyze epitope conservation within MSAs"""
    
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
        self.conservation_analysis = {}
        
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
    
    def calculate_position_conservation(self, msa_position, method='jalview_similarity'):
        """Calculate conservation score for a specific MSA position using advanced methods"""
        amino_acids = []
        for seq_record in self.msa_sequences:
            if msa_position < len(seq_record.seq):
                aa = str(seq_record.seq[msa_position])
                if aa != '-':  # Ignore gaps
                    amino_acids.append(aa)
        
        if not amino_acids:
            return 0.0, {}, len(amino_acids), {}
        
        # Calculate frequency distribution
        aa_counts = Counter(amino_acids)
        total_sequences = len(amino_acids)
        
        # Calculate different conservation scores
        conservation_scores = {
            'blosum62': self.calculate_blosum62_conservation(amino_acids),
            'jalview_similarity': self.calculate_jalview_similarity(amino_acids),
            'identity': self.calculate_identity_conservation(amino_acids)
        }
        
        # Primary score based on method
        primary_score = conservation_scores.get(method, conservation_scores['identity'])
        
        return primary_score, aa_counts, total_sequences, conservation_scores
    
    def analyze_physicochemical_conservation(self, msa_start, msa_end):
        """Analyze conservation of physicochemical properties"""
        property_conservation = {}
        
        for prop_name, prop_aas in AA_PROPERTIES.items():
            position_scores = []
            
            for pos in range(msa_start, msa_end + 1):
                amino_acids = []
                for seq_record in self.msa_sequences:
                    if pos < len(seq_record.seq):
                        aa = str(seq_record.seq[pos])
                        if aa != '-':
                            amino_acids.append(aa)
                
                if amino_acids:
                    # Calculate percentage with this property
                    with_property = sum(1 for aa in amino_acids if aa in prop_aas)
                    conservation = with_property / len(amino_acids)
                    position_scores.append(conservation)
            
            if position_scores:
                property_conservation[prop_name] = {
                    'mean': np.mean(position_scores),
                    'min': np.min(position_scores),
                    'max': np.max(position_scores),
                    'positions': position_scores
                }
        
        return property_conservation
    
    def calculate_blosum62_conservation(self, amino_acids):
        """Calculate BLOSUM62-based conservation score considering amino acid similarity"""
        if len(amino_acids) <= 1:
            return 1.0  # Fully conserved or no data
        
        if not BLOSUM_AVAILABLE:
            # Fallback to simple identity scoring
            most_common = Counter(amino_acids).most_common(1)[0]
            return most_common[1] / len(amino_acids)
        
        try:
            # Use modern Biopython substitution matrices
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
                                # Legacy format uses tuple keys
                                key = (aa1, aa2) if (aa1, aa2) in blosum62 else (aa2, aa1)
                                similarity = blosum62.get(key, -4)
                        except (KeyError, TypeError):
                            similarity = -4  # Default to low similarity for unknown AAs
                        similarities.append(similarity)
            
            if not similarities:
                return 1.0
            
            # Average similarity score
            avg_similarity = np.mean(similarities)
            
            # Normalize BLOSUM62 scores to 0-1 range
            # BLOSUM62 scores typically range from -4 to 11
            normalized_score = (avg_similarity + 4) / 15
            
            # Ensure score is in [0, 1] range
            return max(0.0, min(1.0, normalized_score))
            
        except Exception as e:
            print(f"Warning: Error in BLOSUM62 calculation: {e}, falling back to identity")
            # Fallback to simple identity scoring
            most_common = Counter(amino_acids).most_common(1)[0]
            return most_common[1] / len(amino_acids)
    
    def calculate_jalview_similarity(self, amino_acids):
        """Calculate Jalview % Similarity conservation score"""
        if len(amino_acids) <= 1:
            return 1.0
        
        # Jalview similarity considers both identity and physicochemical similarity
        # Uses a similarity matrix approach similar to BLOSUM but with different scoring
        
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
            similarity_percentage = similar_pairs / total_comparisons
            return similarity_percentage
            
        except Exception as e:
            print(f"Warning: Error in Jalview similarity calculation: {e}")
            return self.calculate_identity_conservation(amino_acids)
    
    def calculate_identity_conservation(self, amino_acids):
        """Calculate simple identity-based conservation score"""
        if len(amino_acids) <= 1:
            return 1.0
        
        # Most frequent amino acid percentage
        most_common = Counter(amino_acids).most_common(1)[0]
        identity_score = most_common[1] / len(amino_acids)
        return identity_score
    
    def calculate_aaindex_conservation(self, amino_acids, property_index='HOPT810101'):
        """Calculate AAindex-based conservation score (placeholder for future implementation)"""
        # This would implement AAindex-based scoring using physicochemical properties
        # For now, fallback to identity conservation
        return self.calculate_identity_conservation(amino_acids)
    
    def get_consurf_scores(self, gene_name, epitope_start, epitope_end):
        """Get ConSurf conservation scores for epitope region (placeholder)"""
        # This would integrate with ConSurf results if available
        # Check for ConSurf results file
        consurf_file = Path(f"results/{self.analysis}_{self.paramset}/protein_analysis/conservation/consurf/{gene_name}_conservation.json")
        
        if consurf_file.exists():
            try:
                with open(consurf_file, 'r') as f:
                    consurf_data = json.load(f)
                
                # Extract scores for epitope region
                epitope_scores = []
                for pos in range(epitope_start, epitope_end + 1):
                    for position_data in consurf_data.get('positions', []):
                        if position_data.get('position') == pos:
                            # ConSurf scores are 1-9, normalize to 0-1
                            normalized_score = (position_data.get('consurf_score', 1) - 1) / 8
                            epitope_scores.append(normalized_score)
                            break
                    else:
                        epitope_scores.append(0.5)  # Default moderate conservation
                
                return np.mean(epitope_scores) if epitope_scores else 0.5
                
            except Exception as e:
                print(f"Warning: Error reading ConSurf data for {gene_name}: {e}")
        
        return None  # ConSurf data not available
    
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
                primary_score, aa_counts, total_seqs, all_scores = self.calculate_position_conservation(pos, 'jalview_similarity')
                
                # Store all conservation scores
                for method in position_conservation:
                    position_conservation[method].append(all_scores.get(method, 0.0))
                
                position_details.append({
                    'msa_position': pos + 1,  # 1-based
                    'structure_position': epitope['start'] + (pos - msa_start),
                    'conservation_scores': all_scores,
                    'primary_conservation': primary_score,
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
            
            # Check for ConSurf scores if available
            consurf_score = self.get_consurf_scores(self.gene, epitope['start'], epitope['end'])
            if consurf_score is not None:
                conservation_summary['consurf'] = {
                    'mean_conservation': consurf_score,
                    'source': 'phylogenetic_analysis'
                }
            
            # Physicochemical conservation
            physchem_conservation = self.analyze_physicochemical_conservation(msa_start, msa_end)
            
            # Consensus sequence
            consensus_25 = self.generate_consensus_sequence(msa_start, msa_end, 0.25)
            consensus_50 = self.generate_consensus_sequence(msa_start, msa_end, 0.5)
            consensus_70 = self.generate_consensus_sequence(msa_start, msa_end, 0.7)
            consensus_90 = self.generate_consensus_sequence(msa_start, msa_end, 0.9)
            
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
                'consensus_sequences': {
                    'consensus_25': consensus_25,
                    'consensus_50': consensus_50,
                    'consensus_70': consensus_70,
                    'consensus_90': consensus_90
                },
                'physicochemical_conservation': physchem_conservation,
                'species_variants': species_variants,
                'variant_count': len(species_variants),
                'most_common_variant': max(species_variants.items(), key=lambda x: len(x[1])) if species_variants else (None, [])
            }
            
            results.append(epitope_analysis)
        
        return results
    
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
        
        # Check for ConSurf scores
        consurf_scores = []
        for analysis in epitope_analyses:
            if 'consurf' in analysis['conservation_scores']:
                consurf_scores.append(analysis['conservation_scores']['consurf']['mean_conservation'])
        
        if consurf_scores:
            method_scores['consurf'] = {
                'mean_conservation': np.mean(consurf_scores),
                'median_conservation': np.median(consurf_scores),
                'std_conservation': np.std(consurf_scores),
                'min_conservation': np.min(consurf_scores),
                'max_conservation': np.max(consurf_scores),
                'source': 'phylogenetic_analysis'
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
            
            # Add physicochemical properties
            for prop_name, prop_data in analysis['physicochemical_conservation'].items():
                base_row[f'{prop_name.title()}_Conservation'] = prop_data['mean']
            
            # Add ConSurf score if available
            if 'consurf' in analysis['conservation_scores']:
                base_row['ConSurf_Score'] = analysis['conservation_scores']['consurf']['mean_conservation']
            
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
        
        </table>
    </div>
"""
            <tr><td>Length</td><td>{epitope['length']} residues</td></tr>
            <tr><td>BepiPred Score</td><td>{epitope['score']:.3f}</td></tr>
            <tr><td>Shannon Conservation</td><td>{shannon_score:.3f}</td></tr>
            <tr><td>BLOSUM62 Conservation</td><td>{blosum_score:.3f}</td></tr>
            <tr><td>Average Conservation</td><td>{avg_score:.3f}</td></tr>
            <tr><td>Sequence Variants</td><td>{analysis['variant_count']}</td></tr>
        </table>
        
        <h4>Consensus Sequences</h4>
        <table>
            <tr><th>Threshold</th><th>Consensus</th></tr>
            <tr><td>25%</td><td class="consensus">{analysis['consensus_sequences']['consensus_25']}</td></tr>
            <tr><td>50%</td><td class="consensus">{analysis['consensus_sequences']['consensus_50']}</td></tr>
            <tr><td>70%</td><td class="consensus">{analysis['consensus_sequences']['consensus_70']}</td></tr>
            <tr><td>90%</td><td class="consensus">{analysis['consensus_sequences']['consensus_90']}</td></tr>
        </table>
        
        <h4>Position-wise Conservation</h4>
        <table class="position-table">
            <tr><th>Structure Pos</th><th>MSA Pos</th><th>Shannon</th><th>BLOSUM62</th><th>Most Common AA</th><th>Frequency</th></tr>
"""
            
            for pos_detail in analysis['position_details']:
                if pos_detail['amino_acid_counts']:
                    most_common = max(pos_detail['amino_acid_counts'].items(), key=lambda x: x[1])
                    frequency = most_common[1] / pos_detail['total_sequences']
                    html_content += f"""
            <tr>
                <td>{pos_detail['structure_position']}</td>
                <td>{pos_detail['msa_position']}</td>
                <td>{pos_detail['shannon_conservation']:.3f}</td>
                <td>{pos_detail['blosum62_conservation']:.3f}</td>
                <td>{most_common[0]}</td>
                <td>{frequency:.2%}</td>
            </tr>"""
            
            html_content += """
        </table>
        
        <h4>Top Sequence Variants</h4>
        <table>
            <tr><th>Variant Sequence</th><th>Species Count</th><th>Representative Species</th></tr>
"""
            
            # Show top 5 variants
            sorted_variants = sorted(analysis['species_variants'].items(), key=lambda x: len(x[1]), reverse=True)[:5]
            for variant_seq, species_list in sorted_variants:
                html_content += f"""
            <tr>
                <td class="consensus">{variant_seq}</td>
                <td>{len(species_list)}</td>
                <td>{', '.join(species_list[:3])}{'...' if len(species_list) > 3 else ''}</td>
            </tr>"""
            
            html_content += """
        </table>
    </div>
"""
        
        html_content += """
</body>
</html>
"""
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        print(f"HTML report written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Analyze epitope conservation in MSA')
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
        print(f"Mean Shannon conservation: {summary_stats['conservation_statistics']['shannon']['mean_conservation']:.3f}")
        print(f"Mean BLOSUM62 conservation: {summary_stats['conservation_statistics']['blosum62']['mean_conservation']:.3f}")
        print(f"Reports saved to: {output_dir}/")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()