#!/usr/bin/env python3
"""
Protein Structure Validation Script
===================================

Validates that downloaded protein structures match the intended gene functions
by cross-checking gene names with protein names and functions.
"""

import json
import pandas as pd
from pathlib import Path
import requests
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_gene_annotations_from_pipeline_data(paramset="params1"):
    """
    Load gene annotations and GO terms from existing pipeline data
    Returns dict mapping gene_name -> expected keywords/functions
    """
    base_dir = Path(__file__).parent.parent.parent
    gene_annotations = {}
    
    try:
        # Load GO validation report with functional annotations
        go_report_path = base_dir / "data" / "quickgo" / paramset / "go_validation_report.tsv"
        if go_report_path.exists():
            go_df = pd.read_csv(go_report_path, sep='\t')
            for _, row in go_df.iterrows():
                gene = row.get('gene_symbol', '').lower()
                go_terms = row.get('go_terms', '')
                if gene and go_terms:
                    # Extract keywords from GO terms
                    keywords = [gene]  # Always include gene name
                    import re
                    go_keywords = re.findall(r'\b\w+\b', go_terms.lower())
                    keywords.extend([kw for kw in go_keywords if len(kw) > 3])  # Filter short words
                    gene_annotations[gene] = keywords
        
        # Load protein information from UniProt data
        protein_info_path = base_dir / "data" / "uniprot" / "protein_info.json"
        if protein_info_path.exists():
            with open(protein_info_path, 'r') as f:
                protein_data = json.load(f)
            
            for gene, info in protein_data.items():
                gene_lower = gene.lower()
                if gene_lower not in gene_annotations:
                    gene_annotations[gene_lower] = []
                
                # Add gene name and aliases
                gene_annotations[gene_lower].append(gene_lower)
                aliases = info.get('aliases', [])
                gene_annotations[gene_lower].extend([alias.lower() for alias in aliases])
                
                # Add protein names and functions if available
                protein_names = info.get('protein_names', [])
                gene_annotations[gene_lower].extend([name.lower() for name in protein_names])
        
        # Load gene aliases
        gene_aliases_path = base_dir / "data" / "quickgo" / paramset / "gene_aliases.json"
        if gene_aliases_path.exists():
            with open(gene_aliases_path, 'r') as f:
                aliases_data = json.load(f)
            
            for gene, aliases in aliases_data.items():
                gene_lower = gene.lower()
                if gene_lower not in gene_annotations:
                    gene_annotations[gene_lower] = []
                gene_annotations[gene_lower].extend([alias.lower() for alias in aliases])
        
        logging.info(f"Loaded validation data for {len(gene_annotations)} genes from pipeline data")
        
    except Exception as e:
        logging.warning(f"Could not load gene annotations from pipeline data: {e}")
        # Fallback to empty mappings
        return {}
    
    return gene_annotations

def validate_protein_structure(gene_name, structure_id, protein_name, organism, paramset="params1"):
    """
    Validate if a protein structure matches the expected gene function using pipeline data
    
    Args:
        gene_name: Target gene name
        structure_id: PDB/AlphaFold structure ID
        protein_name: Protein name from structure metadata
        organism: Source organism
        paramset: Parameter set to load data from
    
    Returns:
        dict: Validation results with confidence score
    """
    # Load expected functions from pipeline data
    gene_annotations = load_gene_annotations_from_pipeline_data(paramset)
    return validate_protein_structure_with_data(gene_name, structure_id, protein_name, organism, gene_annotations)

def validate_protein_structure_with_data(gene_name, structure_id, protein_name, organism, gene_annotations):
    """
    Validate if a protein structure matches the expected gene function using pre-loaded data
    
    Args:
        gene_name: Target gene name
        structure_id: PDB/AlphaFold structure ID
        protein_name: Protein name from structure metadata
        organism: Source organism
        gene_annotations: Pre-loaded gene annotations dict
    
    Returns:
        dict: Validation results with confidence score
    """
    expected_keywords = gene_annotations.get(gene_name.lower(), [])
    
    if not expected_keywords:
        # Fallback: use gene name itself
        expected_keywords = [gene_name.lower()]
    
    protein_name_lower = protein_name.lower()
    warnings = []
    matches = 0
    
    # Check if protein name contains expected keywords
    for keyword in expected_keywords:
        if keyword in protein_name_lower:
            matches += 1
    
    # Special case validations for known problematic patterns
    if gene_name.lower() == 'cls':
        # CLS should be cardiolipin synthase, not LsrB
        if 'lsrb' in protein_name_lower or 'lsr' in protein_name_lower:
            return {
                'valid': False,
                'confidence': 0.1,
                'reason': f'Gene {gene_name} should be cardiolipin synthase, but found LsrB protein',
                'warnings': ['Wrong protein family - LsrB instead of cardiolipin synthase']
            }
        if 'cardiolipin' not in protein_name_lower and 'synthase' not in protein_name_lower:
            warnings.append('CLS gene but no cardiolipin synthase keywords found')
    
    # Calculate confidence based on matches
    confidence = min(matches / max(len(expected_keywords), 1), 1.0) if expected_keywords else 0.5
    
    # Boost confidence if gene name appears in protein description
    if gene_name.lower() in protein_name_lower:
        confidence = min(confidence + 0.3, 1.0)
    
    # Additional checks for flagellar proteins
    if gene_name.lower().startswith('flg') or gene_name.lower().startswith('fli'):
        if 'flagellar' not in protein_name_lower and 'flg' not in protein_name_lower and 'fli' not in protein_name_lower:
            warnings.append(f'Flagellar gene {gene_name} but protein name lacks flagellar keywords')
            confidence *= 0.7
    
    valid = confidence >= 0.3  # Threshold for validity
    
    return {
        'valid': valid,
        'confidence': confidence,
        'reason': f'Matched {matches}/{len(expected_keywords)} expected keywords from pipeline data',
        'warnings': warnings
    }

def create_validation_report(structures_file, output_file):
    """
    Create validation report for all downloaded structures
    
    Args:
        structures_file: Path to downloaded structures TSV
        output_file: Path for validation report
    """
    df = pd.read_csv(structures_file, sep='\t')
    
    # Extract paramset from filename once
    paramset = "params1"  # Default
    if "params" in structures_file:
        import re
        paramset_match = re.search(r'params\d+', structures_file)
        if paramset_match:
            paramset = paramset_match.group()
    
    # Load gene annotations once for efficiency
    gene_annotations = load_gene_annotations_from_pipeline_data(paramset)
    logging.info(f"Loaded validation data for {len(gene_annotations)} genes - using for all validations")
    
    validation_results = []
    
    for _, row in df.iterrows():
        gene_name = row.get('gene_name', '')
        structure_id = row.get('structure_id', '')
        protein_name = row.get('protein_name', '')
        organism = row.get('organism', '')
        
        validation = validate_protein_structure_with_data(gene_name, structure_id, protein_name, organism, gene_annotations)
        
        validation_results.append({
            'gene_name': gene_name,
            'structure_id': structure_id,
            'protein_name': protein_name,
            'organism': organism,
            'valid': validation['valid'],
            'confidence': validation['confidence'],
            'reason': validation['reason'],
            'warnings': '; '.join(validation['warnings'])
        })
    
    # Create validation report
    validation_df = pd.DataFrame(validation_results)
    validation_df.to_csv(output_file, sep='\t', index=False)
    
    # Print summary
    total = len(validation_df)
    valid = len(validation_df[validation_df['valid']])
    invalid = total - valid
    
    print(f"\nProtein Structure Validation Report")
    print(f"=" * 40)
    print(f"Total structures: {total}")
    print(f"Valid structures: {valid} ({valid/total*100:.1f}%)")
    print(f"Invalid structures: {invalid} ({invalid/total*100:.1f}%)")
    
    if invalid > 0:
        print(f"\nInvalid structures:")
        invalid_df = validation_df[~validation_df['valid']]
        for _, row in invalid_df.iterrows():
            print(f"  {row['gene_name']} -> {row['structure_id']}: {row['reason']}")
    
    return validation_df

def main():
    """Main function"""
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python validate_protein_structures.py <structures_tsv> <output_report>")
        print("Example: python validate_protein_structures.py data/protein_structures/analysis1_params1_downloaded_structures.tsv validation_report.tsv")
        sys.exit(1)
    
    structures_file = sys.argv[1]
    output_file = sys.argv[2]
    
    validation_df = create_validation_report(structures_file, output_file)
    print(f"\nValidation report saved to: {output_file}")

if __name__ == "__main__":
    main()