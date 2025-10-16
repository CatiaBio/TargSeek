#!/usr/bin/env python3
"""
Epitope Conservation Analysis using ConSurf Results
Evaluates epitope sequences against conservation scores from ConSurf analysis
"""

import pandas as pd
import numpy as np
import os
import sys
import re
from pathlib import Path
import argparse
from typing import Dict, List, Tuple, Optional

def parse_consurf_grades(consurf_file: str) -> Dict[int, Dict]:
    """
    Parse ConSurf grades file to extract conservation scores
    Returns: {position: {'score': float, 'grade': int, 'aa': str}}
    """
    conservation_data = {}

    if not os.path.exists(consurf_file):
        print(f"ConSurf file not found: {consurf_file}")
        return conservation_data

    with open(consurf_file, 'r') as f:
        lines = f.readlines()

    # Find the start of the data section
    data_start = False
    for i, line in enumerate(lines):
        if "POS	AA	SCORE	COLOR" in line:
            data_start = i + 1
            break

    if not data_start:
        print("Could not find conservation data section in ConSurf file")
        return conservation_data

    # Parse conservation data
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('*'):
            continue

        parts = line.split('\t')
        if len(parts) >= 4:
            try:
                pos = int(parts[0])
                aa = parts[1]
                score = float(parts[2])
                color = int(parts[3]) if parts[3].isdigit() else 5

                conservation_data[pos] = {
                    'aa': aa,
                    'score': score,
                    'grade': color,
                    'conservation_level': get_conservation_level(color)
                }
            except (ValueError, IndexError):
                continue

    return conservation_data

def get_conservation_level(grade: int) -> str:
    """Convert ConSurf grade to conservation level"""
    if grade >= 8:
        return "Highly Conserved"
    elif grade >= 6:
        return "Moderately Conserved"
    elif grade >= 4:
        return "Variable"
    else:
        return "Highly Variable"

def analyze_epitope_conservation(epitope_seq: str, start_pos: int, conservation_data: Dict) -> Dict:
    """
    Analyze conservation of an epitope sequence
    """
    epitope_length = len(epitope_seq)
    epitope_positions = list(range(start_pos, start_pos + epitope_length))

    conservation_scores = []
    conservation_grades = []
    missing_positions = []

    for i, pos in enumerate(epitope_positions):
        if pos in conservation_data:
            conservation_scores.append(conservation_data[pos]['score'])
            conservation_grades.append(conservation_data[pos]['grade'])
        else:
            missing_positions.append(pos)

    if not conservation_scores:
        return {
            'epitope_sequence': epitope_seq,
            'start_position': start_pos,
            'end_position': start_pos + epitope_length - 1,
            'length': epitope_length,
            'avg_conservation_score': None,
            'avg_conservation_grade': None,
            'conservation_level': "Unknown",
            'highly_conserved_residues': 0,
            'variable_residues': 0,
            'missing_positions': missing_positions,
            'conservation_profile': []
        }

    avg_score = np.mean(conservation_scores)
    avg_grade = np.mean(conservation_grades)

    # Count conservation levels
    highly_conserved = sum(1 for grade in conservation_grades if grade >= 8)
    moderately_conserved = sum(1 for grade in conservation_grades if 6 <= grade < 8)
    variable = sum(1 for grade in conservation_grades if 4 <= grade < 6)
    highly_variable = sum(1 for grade in conservation_grades if grade < 4)

    # Overall conservation assessment
    if avg_grade >= 7:
        overall_conservation = "Highly Conserved"
    elif avg_grade >= 5.5:
        overall_conservation = "Moderately Conserved"
    elif avg_grade >= 4:
        overall_conservation = "Variable"
    else:
        overall_conservation = "Highly Variable"

    # Detailed conservation profile
    conservation_profile = []
    for i, pos in enumerate(epitope_positions):
        if pos in conservation_data:
            conservation_profile.append({
                'position': pos,
                'residue': epitope_seq[i] if i < len(epitope_seq) else 'X',
                'score': conservation_data[pos]['score'],
                'grade': conservation_data[pos]['grade'],
                'level': conservation_data[pos]['conservation_level']
            })

    return {
        'epitope_sequence': epitope_seq,
        'start_position': start_pos,
        'end_position': start_pos + epitope_length - 1,
        'length': epitope_length,
        'avg_conservation_score': round(avg_score, 3),
        'avg_conservation_grade': round(avg_grade, 1),
        'conservation_level': overall_conservation,
        'highly_conserved_residues': highly_conserved,
        'moderately_conserved_residues': moderately_conserved,
        'variable_residues': variable,
        'highly_variable_residues': highly_variable,
        'conservation_percentage': round((highly_conserved + moderately_conserved) / len(conservation_scores) * 100, 1),
        'missing_positions': missing_positions,
        'conservation_profile': conservation_profile
    }

def main():
    parser = argparse.ArgumentParser(description='Analyze epitope conservation using ConSurf results')
    parser.add_argument('--consurf-file', required=True, help='Path to ConSurf grades file')
    parser.add_argument('--epitopes', required=True, help='Comma-separated epitope sequences')
    parser.add_argument('--positions', help='Comma-separated start positions (if known)')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--protein-name', default='protein', help='Protein name for output files')

    args = parser.parse_args()

    # Parse ConSurf results
    print(f"Loading ConSurf conservation data from: {args.consurf_file}")
    conservation_data = parse_consurf_grades(args.consurf_file)

    if not conservation_data:
        print("No conservation data found. Exiting.")
        return

    print(f"Loaded conservation data for {len(conservation_data)} positions")

    # Parse epitopes
    epitope_sequences = [seq.strip() for seq in args.epitopes.split(',')]

    # Parse positions if provided
    positions = []
    if args.positions:
        positions = [int(pos.strip()) for pos in args.positions.split(',')]

    # If positions not provided, search for epitopes in sequence
    if not positions:
        positions = [1] * len(epitope_sequences)  # Default to position 1
        print("Warning: Epitope positions not provided, using default position 1")

    # Ensure we have same number of epitopes and positions
    while len(positions) < len(epitope_sequences):
        positions.append(1)

    # Analyze each epitope
    results = []
    for i, epitope_seq in enumerate(epitope_sequences):
        if i < len(positions):
            start_pos = positions[i]
        else:
            start_pos = 1

        print(f"\\nAnalyzing epitope {i+1}: {epitope_seq} (starting at position {start_pos})")

        result = analyze_epitope_conservation(epitope_seq, start_pos, conservation_data)
        results.append(result)

        # Print summary
        print(f"  Conservation level: {result['conservation_level']}")
        print(f"  Average conservation grade: {result['avg_conservation_grade']}")
        print(f"  Conservation percentage: {result['conservation_percentage']}%")
        print(f"  Highly conserved residues: {result['highly_conserved_residues']}/{result['length']}")

    # Create output DataFrame
    df_results = pd.DataFrame(results)

    # Save results
    os.makedirs(args.output_dir, exist_ok=True)

    output_file = os.path.join(args.output_dir, f"{args.protein_name}_epitope_conservation_analysis.csv")
    df_results.to_csv(output_file, index=False)

    # Create detailed report
    report_file = os.path.join(args.output_dir, f"{args.protein_name}_epitope_conservation_report.txt")
    with open(report_file, 'w') as f:
        f.write(f"Epitope Conservation Analysis Report - {args.protein_name}\\n")
        f.write("=" * 60 + "\\n\\n")

        for i, result in enumerate(results):
            f.write(f"Epitope {i+1}: {result['epitope_sequence']}\\n")
            f.write(f"Position: {result['start_position']}-{result['end_position']}\\n")
            f.write(f"Length: {result['length']} residues\\n")
            f.write(f"Conservation Level: {result['conservation_level']}\\n")
            f.write(f"Average Conservation Grade: {result['avg_conservation_grade']}\\n")
            f.write(f"Average Conservation Score: {result['avg_conservation_score']}\\n")
            f.write(f"Conservation Percentage: {result['conservation_percentage']}%\\n")
            f.write(f"Highly Conserved Residues: {result['highly_conserved_residues']}\\n")
            f.write(f"Moderately Conserved Residues: {result['moderately_conserved_residues']}\\n")
            f.write(f"Variable Residues: {result['variable_residues']}\\n")
            f.write(f"Highly Variable Residues: {result['highly_variable_residues']}\\n")

            if result['conservation_profile']:
                f.write("\\nDetailed Conservation Profile:\\n")
                for residue_data in result['conservation_profile']:
                    f.write(f"  {residue_data['position']:3d} {residue_data['residue']} "
                           f"{residue_data['score']:6.3f} {residue_data['grade']:2d} "
                           f"{residue_data['level']}\\n")

            f.write("\\n" + "-" * 40 + "\\n\\n")

    print(f"\\nResults saved to:")
    print(f"  {output_file}")
    print(f"  {report_file}")

if __name__ == "__main__":
    main()