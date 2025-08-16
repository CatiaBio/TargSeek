#!/usr/bin/env python3
"""
Create comprehensive epitope visualizations for any gene and gram type.
Generates all the plot types we developed for bamA analysis.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re
import argparse
import sys
import subprocess
import tempfile
import json

def parse_consurf_detailed(consurf_file):
    """Parse ConSurf data with all amino acid percentages"""
    
    with open(consurf_file, 'r') as f:
        lines = f.readlines()
    
    # Find the header line
    header_line = None
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith('pos,'):
            header_line = line.strip()
            data_start = i + 1
            break
    
    if header_line is None:
        raise ValueError("Could not find header line in ConSurf file")
    
    # Parse header to get column names
    columns = [col.strip('"') for col in header_line.split(',')]
    
    # Read data lines
    data_rows = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        values = [val.strip('"') for val in line.split(',')]
        if len(values) == len(columns):
            data_rows.append(values)
    
    # Create DataFrame
    df = pd.DataFrame(data_rows, columns=columns)
    
    # Convert position and conservation score to numeric
    df['pos'] = df['pos'].astype(int)
    df['ConSurf grade'] = df['ConSurf grade'].astype(float)
    
    return df

def parse_consurf_aa_data(consurf_file):
    """Parse ConSurf amino acid percentage data"""
    
    with open(consurf_file, 'r') as f:
        lines = f.readlines()
    
    # Find the header line
    header_line = None
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith('pos,'):
            header_line = line.strip()
            data_start = i + 1
            break
    
    if header_line is None:
        raise ValueError("Could not find header line in ConSurf file")
    
    # Parse header to get column names
    columns = [col.strip('"') for col in header_line.split(',')]
    
    # Read data lines
    data_rows = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        values = [val.strip('"') for val in line.split(',')]
        if len(values) == len(columns):
            data_rows.append(values)
    
    # Create DataFrame
    df = pd.DataFrame(data_rows, columns=columns)
    
    # Convert position to numeric
    df['pos'] = df['pos'].astype(int)
    
    # Convert amino acid percentages to numeric
    aa_columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    for col in aa_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    return df, aa_columns

def parse_r4s_file(r4s_file):
    """Parse ConSurf r4s.res file to extract amino acid and species count data"""
    
    data = []
    with open(r4s_file, 'r') as f:
        lines = f.readlines()
    
    # Find data lines (skip comments and headers)
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#') and not line.startswith('The alpha'):
            # Parse lines like: "    1     M  -1.337   [-1.361,-1.323] 0.02233  153/157"
            parts = line.split()
            if len(parts) >= 6:
                try:
                    pos = int(parts[0])
                    aa = parts[1]
                    score = float(parts[2])
                    # Extract species count from last part (e.g., "153/157")
                    species_info = parts[-1]
                    if '/' in species_info:
                        present, total = species_info.split('/')
                        present = int(present)
                        total = int(total)
                        percentage = (present / total) * 100
                        
                        data.append({
                            'position': pos,
                            'amino_acid': aa,
                            'conservation_score': score,
                            'species_present': present,
                            'total_species': total,
                            'percentage_present': percentage
                        })
                except (ValueError, IndexError):
                    continue
    
    return pd.DataFrame(data)

def get_epitope_regions(integrated_file):
    """Get epitope regions from integrated analysis file"""
    
    if not os.path.exists(integrated_file):
        raise FileNotFoundError(f"Integrated analysis file not found: {integrated_file}")
    
    epitopes = []
    with open(integrated_file, 'r') as f:
        lines = f.readlines()
        
    for line in lines[1:]:  # Skip header
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                epitope_no = int(parts[3])
                start = int(parts[5])
                end = int(parts[6])
                peptide = parts[7]
                epitopes.append({
                    "no": epitope_no, 
                    "start": start, 
                    "end": end, 
                    "peptide": peptide
                })
    
    # Sort by epitope number
    epitopes.sort(key=lambda x: x['no'])
    return epitopes

def run_blast_specificity_analysis(epitopes, gene, output_dir, blast_db_human=None, blast_db_viral=None, blast_db_fungal=None, blast_db_cow=None, blast_db_sheep=None, blast_db_goat=None, blast_db_horse=None, blast_db_pig=None, blast_db_chicken=None, blast_db_buffalo=None, blast_db_mouse=None, blast_db_custom=None, skip_blast=False):
    """Run BLAST analysis to check epitope specificity against non-bacterial organisms"""
    
    print(f"  🔍 Running BLAST specificity analysis for {len(epitopes)} epitopes...")
    
    # Check if BLAST analysis should be skipped
    if skip_blast:
        print(f"  ⏭️  BLAST specificity analysis skipped (--skip-blast flag set)")
        return []
    
    # Define databases to search against for specificity
    databases = {}
    if blast_db_human:
        databases['human'] = blast_db_human
    if blast_db_viral:
        databases['viral'] = blast_db_viral
    if blast_db_fungal:
        databases['fungal'] = blast_db_fungal
    if blast_db_cow:
        databases['cow'] = blast_db_cow
    if blast_db_sheep:
        databases['sheep'] = blast_db_sheep
    if blast_db_goat:
        databases['goat'] = blast_db_goat
    if blast_db_horse:
        databases['horse'] = blast_db_horse
    if blast_db_pig:
        databases['pig'] = blast_db_pig
    if blast_db_chicken:
        databases['chicken'] = blast_db_chicken
    if blast_db_buffalo:
        databases['buffalo'] = blast_db_buffalo
    if blast_db_mouse:
        databases['mouse'] = blast_db_mouse
    
    # Support for custom databases (comma-separated list of db_name:path pairs)
    if blast_db_custom:
        for custom_db in blast_db_custom.split(','):
            if ':' in custom_db:
                db_name, db_path = custom_db.split(':', 1)
                databases[db_name.strip()] = db_path.strip()
    
    if not databases:
        print(f"  ⚠️  No BLAST databases provided. Use --blast-db-* arguments or --blast-db-custom")
        return []
    
    blast_results = []
    
    for epitope in epitopes:
        epitope_seq = epitope['peptide']
        epitope_id = f"{gene}_epitope_{epitope['no']}"
        
        print(f"    Analyzing epitope {epitope['no']}: {epitope_seq}")
        
        # Create temporary FASTA file for this epitope
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            temp_fasta.write(f">{epitope_id}\n{epitope_seq}\n")
            temp_fasta_path = temp_fasta.name
        
        epitope_results = {
            'epitope_no': epitope['no'],
            'epitope_seq': epitope_seq,
            'start': epitope['start'],
            'end': epitope['end'],
            'hits': {}
        }
        
        # Run BLAST against each database
        for db_name, db_path in databases.items():
            if not os.path.exists(db_path):
                print(f"    Warning: Database {db_name} not found at {db_path}")
                epitope_results['hits'][db_name] = {'status': 'database_not_found', 'hits': []}
                continue
            
            try:
                # Run BLAST search
                blast_cmd = [
                    'blastp',
                    '-query', temp_fasta_path,
                    '-db', db_path,
                    '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle',
                    '-evalue', '1000',  # Relaxed e-value to catch potential cross-reactions
                    '-max_target_seqs', '10',
                    '-word_size', '2'   # Lower word size for short sequences
                ]
                
                result = subprocess.run(blast_cmd, capture_output=True, text=True, timeout=60)
                
                if result.returncode == 0:
                    hits = []
                    for line in result.stdout.strip().split('\n'):
                        if line:
                            fields = line.split('\t')
                            if len(fields) >= 13:
                                hit = {
                                    'subject_id': fields[1],
                                    'identity_pct': float(fields[2]),
                                    'alignment_length': int(fields[3]),
                                    'mismatches': int(fields[4]),
                                    'gaps': int(fields[5]),
                                    'evalue': float(fields[10]),
                                    'bitscore': float(fields[11]),
                                    'description': fields[12] if len(fields) > 12 else 'No description'
                                }
                                
                                # Filter for significant hits (high identity, good coverage)
                                coverage = hit['alignment_length'] / len(epitope_seq) * 100
                                if hit['identity_pct'] >= 70 and coverage >= 80:
                                    hit['coverage_pct'] = coverage
                                    hits.append(hit)
                    
                    epitope_results['hits'][db_name] = {'status': 'success', 'hits': hits}
                    print(f"      {db_name}: {len(hits)} significant hits found")
                else:
                    epitope_results['hits'][db_name] = {'status': 'blast_error', 'hits': []}
                    print(f"      {db_name}: BLAST error")
            
            except subprocess.TimeoutExpired:
                epitope_results['hits'][db_name] = {'status': 'timeout', 'hits': []}
                print(f"      {db_name}: BLAST timeout")
            except Exception as e:
                epitope_results['hits'][db_name] = {'status': 'error', 'hits': []}
                print(f"      {db_name}: Error - {str(e)}")
        
        # Clean up temporary file
        os.unlink(temp_fasta_path)
        
        # Calculate specificity score
        total_hits = sum(len(result['hits']) for result in epitope_results['hits'].values())
        epitope_results['specificity_score'] = max(0, 100 - total_hits * 10)  # Simple scoring
        
        blast_results.append(epitope_results)
    
    # Save BLAST results to JSON file
    blast_output_file = f"{output_dir}/{gene}_blast_specificity.json"
    with open(blast_output_file, 'w') as f:
        json.dump(blast_results, f, indent=2)
    
    # Create summary table
    create_blast_summary_table(blast_results, gene, output_dir)
    
    print(f"  ✅ BLAST specificity analysis completed")
    return blast_results

def create_blast_summary_table(blast_results, gene, output_dir):
    """Create a summary table of BLAST specificity results"""
    
    summary_data = []
    for result in blast_results:
        row = {
            'Epitope': f"Epitope {result['epitope_no']}",
            'Sequence': result['epitope_seq'],
            'Position': f"{result['start']}-{result['end']}",
            'Specificity_Score': result['specificity_score']
        }
        
        # Add hit counts for each database
        for db_name, db_result in result['hits'].items():
            row[f'{db_name.capitalize()}_Hits'] = len(db_result['hits']) if db_result['status'] == 'success' else 'N/A'
        
        # Add risk assessment
        total_hits = sum(len(db_result['hits']) for db_result in result['hits'].values() 
                        if db_result['status'] == 'success')
        if total_hits == 0:
            risk = 'Low'
        elif total_hits <= 2:
            risk = 'Medium'
        else:
            risk = 'High'
        row['Cross_Reactivity_Risk'] = risk
        
        summary_data.append(row)
    
    # Create DataFrame and save
    summary_df = pd.DataFrame(summary_data)
    summary_file = f"{output_dir}/{gene}_blast_specificity_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    
    print(f"    📊 BLAST summary table saved: {summary_file}")

def create_conservation_heatmap(epitope, df, output_dir):
    """Create conservation heatmap for a single epitope"""
    
    # Get epitope data
    epitope_data = df[(df['pos'] >= epitope['start']) & (df['pos'] <= epitope['end'])].copy()
    
    if epitope_data.empty:
        print(f"No conservation data found for epitope {epitope['no']}")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(8, len(epitope_data) * 0.5), 4))
    
    # Prepare conservation data
    conservation_scores = epitope_data['ConSurf grade'].values
    positions = epitope_data['pos'].values
    
    # Create heatmap
    heatmap_data = conservation_scores.reshape(1, -1)
    im = ax.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto', vmin=1, vmax=9)
    
    # Set x-axis with amino acid positions
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(positions, fontsize=10)
    ax.set_xlabel('Amino Acid Position', fontsize=12, fontweight='bold')
    
    # Set y-axis with vertical label
    ax.set_yticks([0])
    ax.set_yticklabels(['Conservation'], fontsize=14, fontweight='bold', rotation=90, va='center')
    
    # Add title
    ax.set_title(f'Epitope {epitope["no"]} Conservation Analysis\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.05)
    cbar.set_label('ConSurf Conservation Score\n(1 = Variable, 9 = Highly Conserved)', 
                   fontsize=10, fontweight='bold')
    
    # Clean up spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    output_file = f"{output_dir}/conservation_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def create_aa_composition_heatmap(epitope, df, aa_columns, output_dir):
    """Create amino acid composition heatmap for a single epitope"""
    
    # Get epitope data
    epitope_data = df[(df['pos'] >= epitope['start']) & (df['pos'] <= epitope['end'])].copy()
    
    if epitope_data.empty:
        print(f"No amino acid data found for epitope {epitope['no']}")
        return
    
    # Prepare amino acid percentage matrix
    aa_matrix = []
    positions = epitope_data['pos'].values
    
    for _, row in epitope_data.iterrows():
        aa_percentages = [float(row[aa]) if aa in row else 0 for aa in aa_columns]
        aa_matrix.append(aa_percentages)
    
    aa_matrix = np.array(aa_matrix).T  # Transpose so amino acids are rows, positions are columns
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(10, len(positions) * 0.8), 8))
    
    # Create heatmap
    im = ax.imshow(aa_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=100)
    
    # Set x-axis (positions)
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels(positions, fontsize=10)
    ax.set_xlabel('Amino Acid Position', fontsize=12, fontweight='bold')
    
    # Set y-axis (amino acids)
    ax.set_yticks(range(len(aa_columns)))
    ax.set_yticklabels(aa_columns, fontsize=10)
    ax.set_ylabel('Amino Acid', fontsize=12, fontweight='bold')
    
    # Add title
    ax.set_title(f'Epitope {epitope["no"]} Amino Acid Composition\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add percentage values on heatmap
    for i in range(len(aa_columns)):
        for j in range(len(positions)):
            value = aa_matrix[i, j]
            if value > 0:  # Only show non-zero values
                text_color = 'white' if value > 50 else 'black'
                ax.text(j, i, f'{value:.0f}%', ha='center', va='center', 
                       fontsize=8, color=text_color, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.05)
    cbar.set_label('Amino Acid Percentage (%)', fontsize=12, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    output_file = f"{output_dir}/aa_composition.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def create_species_count_plots(epitope, df, output_dir):
    """Create species count plots for a single epitope"""
    
    # Get epitope data
    epitope_data = df[(df['position'] >= epitope['start']) & (df['position'] <= epitope['end'])].copy()
    
    if epitope_data.empty:
        print(f"No species count data found for epitope {epitope['no']}")
        return
    
    positions = epitope_data['position'].values
    amino_acids = epitope_data['amino_acid'].values
    species_present = epitope_data['species_present'].values
    total_species = epitope_data['total_species'].values
    percentages = epitope_data['percentage_present'].values
    conservation_scores = epitope_data['conservation_score'].values
    
    # Plot 1: Species count (absolute numbers)
    fig1, ax1 = plt.subplots(figsize=(max(10, len(epitope_data) * 0.8), 6))
    bars1 = ax1.bar(range(len(positions)), species_present, 
                    color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axhline(y=total_species[0], color='red', linestyle='--', alpha=0.7, linewidth=2)
    ax1.set_ylabel('Number of Species', fontsize=12, fontweight='bold')
    ax1.set_title(f'Epitope {epitope["no"]} - Species Count per Position\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Calculate appropriate y-axis limit with padding
    max_species = max(species_present)
    total_species_count = total_species[0]
    max_value = max(max_species, total_species_count)
    
    # Add 10 and round up to next 10
    y_limit = ((max_value + 10) // 10 + 1) * 10
    ax1.set_ylim(0, y_limit)
    
    # Add total species caption in bottom right corner (within plot area)
    ax1.text(0.98, 0.05, f'Total species: {total_species_count}', 
             transform=ax1.transAxes, fontsize=11, fontweight='bold',
             ha='right', va='bottom', 
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='red'))
    
    # Add amino acid labels and values inside bars
    for i, (bar, aa, count) in enumerate(zip(bars1, amino_acids, species_present)):
        height = bar.get_height()
        # Position text in the middle of the bar
        ax1.text(bar.get_x() + bar.get_width()/2., height/2,
                f'{aa}\n{count}', ha='center', va='center', fontsize=10, fontweight='bold', color='white')
    
    # Set x-axis
    ax1.set_xticks(range(len(positions)))
    ax1.set_xticklabels(positions, fontsize=10)
    ax1.set_xlabel('Position', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    output1 = f"{output_dir}/species_count_absolute.png"
    plt.savefig(output1, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    # Plot 2: Percentage present
    fig2, ax2 = plt.subplots(figsize=(max(10, len(epitope_data) * 0.8), 6))
    bars2 = ax2.bar(range(len(positions)), percentages, 
                    color='forestgreen', alpha=0.7, edgecolor='black')
    ax2.set_ylabel('Percentage of Species (%)', fontsize=12, fontweight='bold')
    ax2.set_title(f'Epitope {epitope["no"]} - Percentage of Species with Each Amino Acid\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                  fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 100)
    
    # Add percentage values inside bars
    for i, (bar, aa, pct) in enumerate(zip(bars2, amino_acids, percentages)):
        height = bar.get_height()
        # Position text in the middle of the bar
        ax2.text(bar.get_x() + bar.get_width()/2., height/2,
                f'{aa}\n{pct:.1f}%', ha='center', va='center', fontsize=10, fontweight='bold', color='white')
    
    # Set x-axis
    ax2.set_xticks(range(len(positions)))
    ax2.set_xticklabels(positions, fontsize=10)
    ax2.set_xlabel('Position', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    output2 = f"{output_dir}/species_count_percentage.png"
    plt.savefig(output2, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    # Plot 3: Conservation score vs species presence
    fig3, ax3 = plt.subplots(figsize=(10, 8))
    scatter = ax3.scatter(percentages, conservation_scores, 
                         c=range(len(positions)), cmap='viridis', 
                         s=100, alpha=0.7, edgecolors='black')
    ax3.set_xlabel('Percentage of Species with Amino Acid (%)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Normalized Rate4Site Conservation Score', fontsize=12, fontweight='bold')
    ax3.set_title(f'Epitope {epitope["no"]} - Normalized Rate4Site Conservation Score vs Species Presence\n{epitope["peptide"]} (positions {epitope["start"]}-{epitope["end"]})', 
                  fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Add amino acid labels with smart angle adjustment to avoid overlaps
    import math
    
    # Calculate all dot positions in figure coordinates
    fig_width, fig_height = fig3.get_size_inches()
    x_range = ax3.get_xlim()[1] - ax3.get_xlim()[0]
    y_range = ax3.get_ylim()[1] - ax3.get_ylim()[0]
    
    # Convert data coordinates to figure coordinates for distance calculations
    dot_positions = []
    for pct, score in zip(percentages, conservation_scores):
        x_fig = (pct - ax3.get_xlim()[0]) / x_range * fig_width * 72  # Convert to points
        y_fig = (score - ax3.get_ylim()[0]) / y_range * fig_height * 72
        dot_positions.append((x_fig, y_fig))
    
    used_angles = []  # Track used angles for each dot
    offset_distance = 20  # Distance in points (reduced from 35)
    
    for i, (pct, score, aa, pos) in enumerate(zip(percentages, conservation_scores, amino_acids, positions)):
        # Find best angle for this label to avoid overlaps
        best_angle = None
        min_conflicts = float('inf')
        
        # Try different angles in 15-degree increments
        for test_angle in range(0, 360, 15):
            test_angle_rad = math.radians(test_angle)
            
            # Calculate where this label would be positioned
            label_x = dot_positions[i][0] + offset_distance * math.cos(test_angle_rad)
            label_y = dot_positions[i][1] + offset_distance * math.sin(test_angle_rad)
            
            # Count conflicts with existing labels
            conflicts = 0
            for j, other_angles in enumerate(used_angles):
                if j != i:  # Don't compare with self
                    for other_angle in other_angles:
                        other_angle_rad = math.radians(other_angle)
                        other_label_x = dot_positions[j][0] + offset_distance * math.cos(other_angle_rad)
                        other_label_y = dot_positions[j][1] + offset_distance * math.sin(other_angle_rad)
                        
                        # Calculate distance between labels
                        distance = ((label_x - other_label_x)**2 + (label_y - other_label_y)**2)**0.5
                        if distance < 25:  # Minimum separation in points
                            conflicts += 1
            
            # Choose angle with fewest conflicts
            if conflicts < min_conflicts:
                min_conflicts = conflicts
                best_angle = test_angle
        
        # If no good angle found, use default
        if best_angle is None:
            best_angle = (i * 45) % 360  # Fallback
        
        # Store this angle as used
        used_angles.append([best_angle])
        
        # Convert angle to radians
        angle_rad = math.radians(best_angle)
        
        # Calculate initial label position
        label_offset_x = offset_distance * math.cos(angle_rad)
        label_offset_y = offset_distance * math.sin(angle_rad)
        
        # Check if label would go outside plot boundaries and adjust
        # Get plot boundaries in points
        x_min, x_max = ax3.get_xlim()
        y_min, y_max = ax3.get_ylim()
        
        # Convert current dot position to normalized coordinates (0-1)
        x_norm = (pct - x_min) / (x_max - x_min)
        y_norm = (score - y_min) / (y_max - y_min)
        
        # Adjust label position if too close to boundaries
        if x_norm > 0.85 and label_offset_x > 0:  # Too close to right edge, move left
            label_offset_x = -abs(label_offset_x)
        elif x_norm < 0.15 and label_offset_x < 0:  # Too close to left edge, move right
            label_offset_x = abs(label_offset_x)
            
        if y_norm > 0.85 and label_offset_y > 0:  # Too close to top edge, move down
            label_offset_y = -abs(label_offset_y)
        elif y_norm < 0.15 and label_offset_y < 0:  # Too close to bottom edge, move up
            label_offset_y = abs(label_offset_y)
        
        ax3.annotate(f'{aa}{pos}', 
                    xy=(pct, score),  # Point to annotate (the dot)
                    xytext=(label_offset_x, label_offset_y),  # Adjusted label offset
                    textcoords='offset points',  # Use points for offset
                    fontsize=8, fontweight='bold',
                    ha='center', va='center',
                    arrowprops=dict(arrowstyle='-', 
                                  color='black', alpha=0.0, lw=1.5))  # Made invisible (alpha=0.0)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Position Index', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    output3 = f"{output_dir}/conservation_vs_presence.png"
    plt.savefig(output3, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def create_general_conservation_plot(gene, conservation_df, epitopes, output_dir):
    """Create general conservation scores plot for protein region around epitopes"""
    
    # Calculate plot range: 10 aa before first epitope to 10 aa after last epitope
    min_epitope_pos = min(epitope['start'] for epitope in epitopes)
    max_epitope_pos = max(epitope['end'] for epitope in epitopes)
    
    plot_start = max(1, min_epitope_pos - 10)  # Don't go below position 1
    plot_end = max_epitope_pos + 10
    
    # Filter data to plot range
    plot_data = conservation_df[(conservation_df['pos'] >= plot_start) & (conservation_df['pos'] <= plot_end)].copy()
    
    if plot_data.empty:
        print(f"No conservation data in plot range {plot_start}-{plot_end}")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 6))
    
    # Prepare data for plot region
    plot_positions = plot_data['pos'].values
    # Use original ConSurf scores (1-9 scale)
    plot_conservation = plot_data['ConSurf grade'].values
    
    # Create line plot without dots
    ax.plot(plot_positions, plot_conservation, '-', linewidth=2, color='blue')
    
    # Add conservation threshold lines
    ax.axhline(y=7, color='red', linestyle='--', linewidth=1.5, alpha=0.7, 
              label='Highly Conserved (>7)')
    ax.axhline(y=4, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, 
              label='Moderately Conserved (4-7)')
    
    # Add epitope markers and highlights with original style
    epitope_colors = ['red', 'green', 'orange', 'purple', 'brown', 'pink']
    
    for epitope in epitopes:
        color = epitope_colors[(epitope['no'] - 1) % len(epitope_colors)]
        
        # Add background shading for epitope region
        ax.axvspan(epitope['start'], epitope['end'], alpha=0.2, color=color)
        
        # Add epitope label with simpler style
        center_pos = (epitope['start'] + epitope['end']) / 2
        max_conservation = max(plot_conservation)
        ax.text(center_pos, max_conservation * 0.95, f"E{epitope['no']}", 
               ha='center', va='center', fontsize=10, fontweight='bold', color='black')
    
    ax.set_xlabel(f'Amino Acid Position in {gene} Protein', fontsize=12, fontweight='bold')
    ax.set_ylabel('ConSurf Conservation Score', fontsize=12, fontweight='bold')
    ax.set_title(f'{gene} - Conservation Scores (Positions {plot_start}-{plot_end})', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 10)
    
    # Add legend for conservation thresholds
    ax.legend(loc='upper right', fontsize=10)
    
    # Set x-axis to show epitope region with buffer
    ax.set_xlim(plot_start, plot_end)
    
    plt.tight_layout()
    output_file = f"{output_dir}/{gene}_1_conservation_scores.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def create_general_diversity_plot(gene, aa_df, epitopes, output_dir):
    """Create general diversity index plot for protein region around epitopes"""
    
    # Calculate plot range: 10 aa before first epitope to 10 aa after last epitope
    min_epitope_pos = min(epitope['start'] for epitope in epitopes)
    max_epitope_pos = max(epitope['end'] for epitope in epitopes)
    
    plot_start = max(1, min_epitope_pos - 10)  # Don't go below position 1
    plot_end = max_epitope_pos + 10
    
    # Filter data to plot range
    plot_data = aa_df[(aa_df['pos'] >= plot_start) & (aa_df['pos'] <= plot_end)].copy()
    
    if plot_data.empty:
        print(f"No diversity data in plot range {plot_start}-{plot_end}")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 6))
    
    # Calculate diversity index from amino acid variety percentages
    aa_columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    
    plot_positions = plot_data['pos'].values
    plot_diversity = []
    
    for _, row in plot_data.iterrows():
        # Get percentages for all amino acids at this position
        aa_percentages = [float(row[aa]) if aa in row else 0 for aa in aa_columns]
        
        # Calculate diversity using Shannon entropy or variety measure
        # Simple diversity: 100 - max percentage (higher = more diverse)
        max_percentage = max(aa_percentages)
        diversity_index = (100 - max_percentage) / 100  # Normalize to 0-1
        plot_diversity.append(diversity_index)
    
    plot_diversity = np.array(plot_diversity)
    
    # Create line plot without dots
    ax.plot(plot_positions, plot_diversity, '-', linewidth=2, color='red')
    
    # Add diversity threshold lines
    ax.axhline(y=0.7, color='green', linestyle='--', linewidth=1.5, alpha=0.7, 
              label='Good Conservation (<0.7)')
    ax.axhline(y=0.8, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, 
              label='Poor Conservation (>0.8)')
    
    # Add epitope markers and highlights with original style
    epitope_colors = ['red', 'green', 'orange', 'purple', 'brown', 'pink']
    
    for epitope in epitopes:
        color = epitope_colors[(epitope['no'] - 1) % len(epitope_colors)]
        
        # Add background shading for epitope region
        ax.axvspan(epitope['start'], epitope['end'], alpha=0.2, color=color)
        
        # Add epitope label with simpler style
        center_pos = (epitope['start'] + epitope['end']) / 2
        max_diversity = max(plot_diversity)
        ax.text(center_pos, max_diversity * 0.95, f"E{epitope['no']}", 
               ha='center', va='center', fontsize=10, fontweight='bold', color='black')
    
    ax.set_xlabel(f'Amino Acid Position in {gene} Protein', fontsize=12, fontweight='bold')
    ax.set_ylabel('Amino acid diversity index', fontsize=12, fontweight='bold')
    ax.set_title(f'{gene} - Amino Acid Diversity (Positions {plot_start}-{plot_end})', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.0)
    
    # Add legend for diversity thresholds
    ax.legend(loc='upper right', fontsize=10)
    
    # Set x-axis to show epitope region with buffer
    ax.set_xlim(plot_start, plot_end)
    
    plt.tight_layout()
    output_file = f"{output_dir}/{gene}_2_diversity_index.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()


def create_combined_conservation_heatmap(gene, conservation_df, epitopes, output_dir):
    """Create combined plot with full sequence and epitope zoom conservation heatmaps"""
    
    # Sort conservation data by position
    conservation_df = conservation_df.sort_values('pos')
    
    # Create figure with two subplots: full sequence on top, epitopes below
    # Reduced spacing between plots
    fig = plt.figure(figsize=(20, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.15)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    
    # ==== TOP PANEL: Full sequence conservation heatmap ====
    positions = conservation_df['pos'].values
    conservation_scores = conservation_df['ConSurf grade'].values
    
    # Create full sequence heatmap
    heatmap_data = conservation_scores.reshape(1, -1)
    im1 = ax1.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto', vmin=1, vmax=9)
    
    # Set x-axis for full sequence
    tick_interval = max(1, len(positions) // 20)
    tick_positions = range(0, len(positions), tick_interval)
    tick_labels = [str(positions[i]) for i in tick_positions]
    
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels, fontsize=10)
    ax1.set_xlabel(f'Amino Acid Position in {gene} Protein', fontsize=12, fontweight='bold')
    
    # Remove y-axis ticks and labels
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    
    # Add title for full sequence almost attached to plot
    ax1.set_title(f'{gene} - Conservation Analysis: Full Sequence & Epitope Regions', 
                 fontsize=16, fontweight='bold', pad=1)
    
    # Remove epitope colored boxes from full sequence plot
    epitope_colors = ['red', 'green', 'orange', 'purple', 'brown', 'pink']
    
    ax1.set_xlim(-0.5, len(positions) - 0.5)
    ax1.set_ylim(-0.1, 1.1)
    
    # ==== BOTTOM PANEL: Epitope regions heatmap (zoomed) ====
    # Prepare data for all epitopes with spacing between them
    all_data = []
    all_positions = []
    epitope_indices = []
    epitope_boundaries = []  # Track where each epitope starts and ends
    current_pos = 0
    
    for i, epitope in enumerate(epitopes):
        epitope_data = conservation_df[(conservation_df['pos'] >= epitope['start']) & (conservation_df['pos'] <= epitope['end'])].copy()
        
        epitope_start_pos = current_pos
        
        for _, row in epitope_data.iterrows():
            pos = int(row['pos'])
            conservation = row['ConSurf grade']
            all_data.append(conservation)
            all_positions.append(pos)
            epitope_indices.append(epitope['no'])
            current_pos += 1
        
        epitope_end_pos = current_pos - 1
        epitope_boundaries.append((epitope_start_pos, epitope_end_pos, epitope['no']))
        
        # Add gap between epitopes (except after the last one)
        if i < len(epitopes) - 1:
            # Add 2 empty positions as spacer
            all_data.extend([0, 0])  # Use 0 as placeholder for gap
            all_positions.extend([-1, -1])  # Use -1 as placeholder position
            epitope_indices.extend([0, 0])  # Use 0 as placeholder epitope
            current_pos += 2
    
    # Create epitope heatmap
    if all_data:  # Only create if we have epitope data
        # Replace gap values (0) with NaN for proper visualization
        epitope_data_with_gaps = np.array(all_data, dtype=float)
        epitope_data_with_gaps[epitope_data_with_gaps == 0] = np.nan
        
        epitope_heatmap_data = epitope_data_with_gaps.reshape(1, -1)
        im2 = ax2.imshow(epitope_heatmap_data, cmap='RdYlBu_r', aspect='auto', vmin=1, vmax=9)
        
        # Remove x-axis numbering for epitopes
        ax2.set_xticks([])
        ax2.set_xlabel('')
        
        # Remove y-axis ticks and labels
        ax2.set_yticks([])
        ax2.set_yticklabels([])
        
        # Add epitope names with start/end positions below heatmap
        for start_idx, end_idx, epitope_no in epitope_boundaries:
            center_idx = (start_idx + end_idx) / 2
            
            # Find the corresponding epitope to get start/end positions
            epitope_info = epitopes[epitope_no - 1]
            epitope_start = epitope_info['start']
            epitope_end = epitope_info['end']
            
            # Add epitope name with AA positions with spacing between plot and epitope
            ax2.text(center_idx, -0.6, f"Epitope {epitope_no}\n({epitope_start}-{epitope_end})", 
                   ha='center', va='top', fontsize=10, fontweight='bold', color='black')
        
        ax2.set_xlim(-1, len(all_data))
        ax2.set_ylim(-0.9, 1.1)  # Extended bottom margin for epitope labels with proper spacing
    
    # Add smaller colorbar positioned properly to the right
    # Adjust subplot positions to make room for colorbar
    plt.subplots_adjust(bottom=0.12, right=0.85)
    
    # Get updated positions after adjustment
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    
    # Position much smaller colorbar to the right of the plots, reduced height
    # Calculate colorbar height to be 60% of the total plot height
    colorbar_height = (pos1.y1 - pos2.y0) * 0.6
    colorbar_y_center = pos2.y0 + (pos1.y1 - pos2.y0) / 2
    colorbar_y_start = colorbar_y_center - colorbar_height / 2
    
    cbar_ax = fig.add_axes([0.87, colorbar_y_start, 0.01, colorbar_height])
    cbar = plt.colorbar(im1, cax=cbar_ax)
    cbar.set_label('ConSurf Conservation Score\n(1 = Variable, 9 = Highly Conserved)', 
                   fontsize=9, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)
    
    # Remove epitope legend since we removed the colored boxes
    
    # Clean up spines for both axes
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    
    # Save plot
    output_file = f"{output_dir}/{gene}_conservation_heatmap_combined.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.4)
    plt.close()

def create_comprehensive_epitope_visualizations(gene, structure_id, gram_type, analysis, paramset, blast_args=None):
    """Create comprehensive epitope visualizations for a gene"""
    
    # Construct file paths
    base_dir = f"results/{analysis}_{paramset}/protein_analysis"
    
    integrated_file = f"{base_dir}/integrated_epitope_analysis/gram_{gram_type}/{gene}/{gene}_{structure_id}_integrated_epitope_analysis.tsv"
    consurf_file = f"{base_dir}/consurf_analysis/gram_{gram_type}/{gene}/msa_aa_variety_percentage.csv"
    r4s_file = f"{base_dir}/consurf_analysis/gram_{gram_type}/{gene}/r4s.res"
    
    # Check if input files exist
    for file_path, file_type in [(integrated_file, "Integrated analysis"), 
                                (consurf_file, "ConSurf"), 
                                (r4s_file, "Rate4Site")]:
        if not os.path.exists(file_path):
            print(f"Warning: {file_type} file not found: {file_path}")
            if file_type == "Integrated analysis":
                print(f"Skipping {gene} - no integrated analysis found")
                return
    
    # Load data
    epitopes = get_epitope_regions(integrated_file)
    
    if not epitopes:
        print(f"No epitopes found for {gene}")
        return
    
    # Gene-level output directory
    gene_output_dir = f"{base_dir}/integrated_epitope_analysis/gram_{gram_type}/{gene}"
    os.makedirs(gene_output_dir, exist_ok=True)
    
    # Load ConSurf data if available
    conservation_df = None
    aa_df = None
    aa_columns = None
    if os.path.exists(consurf_file):
        try:
            conservation_df = parse_consurf_detailed(consurf_file)
            aa_df, aa_columns = parse_consurf_aa_data(consurf_file)
        except Exception as e:
            print(f"Warning: Could not parse ConSurf file for {gene}: {e}")
    
    # Load Rate4Site data if available
    r4s_df = None
    if os.path.exists(r4s_file):
        try:
            r4s_df = parse_r4s_file(r4s_file)
        except Exception as e:
            print(f"Warning: Could not parse Rate4Site file for {gene}: {e}")
    
    # Create general (gene-level) plots
    print(f"Creating general plots for {gene}...")
    
    if conservation_df is not None:
        try:
            create_general_conservation_plot(gene, conservation_df, epitopes, gene_output_dir)
            print(f"  ✅ General conservation plot created")
        except Exception as e:
            print(f"  ❌ General conservation plot failed: {e}")
        
        try:
            create_combined_conservation_heatmap(gene, conservation_df, epitopes, gene_output_dir)
            print(f"  ✅ Combined conservation heatmap created")
        except Exception as e:
            print(f"  ❌ Combined conservation heatmap failed: {e}")
    
    if aa_df is not None:
        try:
            create_general_diversity_plot(gene, aa_df, epitopes, gene_output_dir)
            print(f"  ✅ General diversity plot created")
        except Exception as e:
            print(f"  ❌ General diversity plot failed: {e}")
    
    # Create epitope-specific visualizations
    for epitope in epitopes:
        # Create epitope-specific directory
        epitope_output_dir = f"{base_dir}/integrated_epitope_analysis/gram_{gram_type}/{gene}/epitope{epitope['no']}_{epitope['peptide']}"
        os.makedirs(epitope_output_dir, exist_ok=True)
        
        print(f"Creating epitope-specific plots for {gene} epitope {epitope['no']}: {epitope['peptide']}")
        
        # Create amino acid composition heatmap (per epitope)
        if aa_df is not None and aa_columns is not None:
            try:
                create_aa_composition_heatmap(epitope, aa_df, aa_columns, epitope_output_dir)
                print(f"  ✅ AA composition heatmap created")
            except Exception as e:
                print(f"  ❌ AA composition heatmap failed: {e}")
        
        # Create species count plots (per epitope)
        if r4s_df is not None:
            try:
                create_species_count_plots(epitope, r4s_df, epitope_output_dir)
                print(f"  ✅ Species count plots created")
            except Exception as e:
                print(f"  ❌ Species count plots failed: {e}")
    
    # Run BLAST specificity analysis
    if blast_args and not blast_args.skip_blast:
        try:
            blast_results = run_blast_specificity_analysis(
                epitopes, gene, gene_output_dir, 
                blast_db_human=blast_args.blast_db_human,
                blast_db_viral=blast_args.blast_db_viral, 
                blast_db_fungal=blast_args.blast_db_fungal,
                blast_db_cow=blast_args.blast_db_cow,
                blast_db_sheep=blast_args.blast_db_sheep,
                blast_db_goat=blast_args.blast_db_goat,
                blast_db_horse=blast_args.blast_db_horse,
                blast_db_pig=blast_args.blast_db_pig,
                blast_db_chicken=blast_args.blast_db_chicken,
                blast_db_buffalo=blast_args.blast_db_buffalo,
                blast_db_mouse=blast_args.blast_db_mouse,
                blast_db_custom=blast_args.blast_db_custom,
                skip_blast=blast_args.skip_blast
            )
            print(f"  ✅ BLAST specificity analysis completed")
        except Exception as e:
            print(f"  ⚠️  BLAST specificity analysis failed: {e}")
    else:
        print(f"  ⏭️  BLAST specificity analysis skipped")
    
    print(f"Comprehensive visualizations completed for {gene} ({len(epitopes)} epitopes)")
    print(f"  📁 General plots saved in: {gene_output_dir}")
    print(f"  📁 Epitope-specific plots saved in individual epitope folders")
    print(f"  🗂️ Combined conservation heatmap created")
    print(f"  🔍 BLAST specificity results saved")

def main():
    parser = argparse.ArgumentParser(description='Create comprehensive epitope visualizations')
    parser.add_argument('--gene', required=True, help='Gene name')
    parser.add_argument('--structure-id', required=True, help='Structure ID (e.g., 5OR1)')
    parser.add_argument('--gram-type', required=True, choices=['positive', 'negative'], 
                       help='Gram type (positive/negative)')
    parser.add_argument('--analysis', default='analysis1', help='Analysis name')
    parser.add_argument('--paramset', default='params1', help='Parameter set')
    parser.add_argument('--blast-db-human', help='Path to human proteome BLAST database')
    parser.add_argument('--blast-db-viral', help='Path to viral proteins BLAST database')  
    parser.add_argument('--blast-db-fungal', help='Path to fungal proteins BLAST database')
    parser.add_argument('--blast-db-cow', help='Path to cow proteome BLAST database')
    parser.add_argument('--blast-db-sheep', help='Path to sheep proteome BLAST database')
    parser.add_argument('--blast-db-goat', help='Path to goat proteome BLAST database')
    parser.add_argument('--blast-db-horse', help='Path to horse proteome BLAST database')
    parser.add_argument('--blast-db-pig', help='Path to pig proteome BLAST database')
    parser.add_argument('--blast-db-chicken', help='Path to chicken proteome BLAST database')
    parser.add_argument('--blast-db-buffalo', help='Path to buffalo proteome BLAST database')
    parser.add_argument('--blast-db-mouse', help='Path to mouse proteome BLAST database')
    parser.add_argument('--blast-db-custom', help='Custom databases as comma-separated db_name:path pairs (e.g., "rat:path/to/rat.fasta,salmon:path/to/salmon.fasta")')
    parser.add_argument('--skip-blast', action='store_true', help='Skip BLAST specificity analysis')
    
    args = parser.parse_args()
    
    create_comprehensive_epitope_visualizations(
        args.gene, args.structure_id, args.gram_type, args.analysis, args.paramset, blast_args=args
    )

if __name__ == "__main__":
    main()