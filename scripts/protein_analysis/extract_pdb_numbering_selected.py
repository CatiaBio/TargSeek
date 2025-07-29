#!/usr/bin/env python3

"""
Extract PDB numbering for selected 3D structures only
This version reads from selected_3d_paths files instead of the full mapping file
Used in the analysis pipeline to process only the structures selected for epitope analysis
"""

import os
import sys
import gzip
import pandas as pd
import re
from pathlib import Path

def extract_structure_chain_info(structure_file_path):
    """
    Extract chain information (amino acid residue numbering) from PDB or mmCIF structure files
    
    Args:
        structure_file_path: Path to structure file (PDB or mmCIF format, compressed or uncompressed)
    
    Returns:
        dict: Dictionary mapping chain IDs to residue information
              Format: {chain_id: {'min_res': int, 'max_res': int}}
    """
    
    structure_path = Path(structure_file_path)
    
    # Check if file exists
    if not structure_path.exists():
        print(f"Structure file not found: {structure_file_path}")
        return {}
    
    chains = {}
    
    try:
        # Handle compressed files
        if structure_path.name.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
        
        with open_func(structure_path, mode) as f:
            # Determine file format
            first_line = f.readline().strip()
            f.seek(0)  # Reset to beginning
            
            if structure_path.name.endswith(('.cif', '.cif.gz')) or first_line.startswith('data_'):
                # mmCIF format
                chains = _parse_mmcif_chains(f)
            else:
                # PDB format
                chains = _parse_pdb_chains(f)
                
    except Exception as e:
        print(f"Error reading structure file {structure_file_path}: {e}")
        return {}
    
    return chains

def _parse_pdb_chains(file_handle):
    """Parse PDB format file for chain information"""
    chains = {}
    
    for line in file_handle:
        if line.startswith(('ATOM', 'HETATM')):
            # Extract chain ID and residue number
            chain_id = line[21].strip()
            if not chain_id:
                chain_id = 'A'  # Default chain if missing
            
            try:
                res_num = int(line[22:26].strip())
                
                if chain_id not in chains:
                    chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                else:
                    chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                    chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                    
            except ValueError:
                continue  # Skip invalid residue numbers
    
    return chains

def _parse_mmcif_chains(file_handle):
    """Parse mmCIF format file for chain information"""
    chains = {}
    in_atom_site = False
    chain_col = None
    seq_id_col = None
    
    for line in file_handle:
        line = line.strip()
        
        if line.startswith('_atom_site.'):
            in_atom_site = True
            # Find column indices
            if 'auth_asym_id' in line:
                chain_col = len([x for x in line.split('.')[0].split('_') if x]) - 1
            elif 'auth_seq_id' in line:
                seq_id_col = len([x for x in line.split('.')[0].split('_') if x]) - 1
        elif in_atom_site and not line.startswith('_') and not line.startswith('#'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    parts = line.split()
                    if chain_col is not None and seq_id_col is not None:
                        chain_id = parts[chain_col] if chain_col < len(parts) else 'A'
                        res_num = int(parts[seq_id_col]) if seq_id_col < len(parts) else 0
                        
                        if chain_id not in chains:
                            chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                        else:
                            chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                            chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                except (ValueError, IndexError):
                    continue
        elif line.startswith('#') and in_atom_site:
            break  # End of atom_site loop
    
    return chains

def extract_enhanced_chain_info(structure_file_path):
    """
    Enhanced chain information extraction for problematic structures
    Uses alternative parsing strategies for mmCIF files and complex structures
    
    Args:
        structure_file_path: Path to structure file
    
    Returns:
        dict: Enhanced chain information mapping
    """
    structure_path = Path(structure_file_path)
    
    if not structure_path.exists():
        return {}
    
    chains = {}
    
    try:
        # Handle compressed files
        if structure_path.name.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
        
        with open_func(structure_path, mode) as f:
            # Determine file format
            first_line = f.readline().strip()
            f.seek(0)  # Reset to beginning
            
            if structure_path.name.endswith(('.cif', '.cif.gz')) or first_line.startswith('data_'):
                # Enhanced mmCIF parsing
                chains = _parse_mmcif_chains_enhanced(f)
            else:
                # Enhanced PDB parsing
                chains = _parse_pdb_chains_enhanced(f)
                
    except Exception as e:
        print(f"    Enhanced parsing error for {structure_file_path}: {e}")
        return {}
    
    return chains

def _parse_pdb_chains_enhanced(file_handle):
    """Enhanced PDB format parsing with better chain detection"""
    chains = {}
    
    for line in file_handle:
        if line.startswith(('ATOM', 'HETATM')):
            # Extract chain ID and residue number with better error handling
            try:
                chain_id = line[21].strip()
                if not chain_id or chain_id == ' ':
                    chain_id = 'A'  # Default chain if missing
                
                # Try different residue number positions
                res_num_str = line[22:26].strip()
                if not res_num_str:
                    res_num_str = line[23:27].strip()  # Alternative position
                
                if res_num_str:
                    res_num = int(res_num_str)
                    
                    if chain_id not in chains:
                        chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                    else:
                        chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                        chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                        
            except (ValueError, IndexError):
                continue  # Skip invalid lines
    
    return chains

def _parse_mmcif_chains_enhanced(file_handle):
    """Enhanced mmCIF format parsing with multiple strategies"""
    chains = {}
    lines = file_handle.read().split('\n')
    
    # Strategy 1: Look for _atom_site.label_seq_id and _atom_site.label_asym_id
    chain_col = None
    seq_id_col = None
    headers = []
    in_atom_site = False
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Find the atom_site section
        if line.startswith('_atom_site.'):
            in_atom_site = True
            headers.append(line)
            
            if 'label_asym_id' in line or 'auth_asym_id' in line:
                chain_col = len(headers) - 1
            elif 'label_seq_id' in line or 'auth_seq_id' in line:
                seq_id_col = len(headers) - 1
                
        elif in_atom_site and not line.startswith('_') and not line.startswith('#') and line:
            if line.startswith(('ATOM', 'HETATM')) or (chain_col is not None and seq_id_col is not None):
                try:
                    parts = line.split()
                    if len(parts) > max(chain_col or 0, seq_id_col or 0):
                        chain_id = parts[chain_col] if chain_col is not None else 'A'
                        
                        # Try multiple sequence ID positions
                        res_num = None
                        if seq_id_col is not None and seq_id_col < len(parts):
                            try:
                                res_num = int(parts[seq_id_col])
                            except ValueError:
                                pass
                        
                        # Fallback: try to find any numeric field that could be residue number
                        if res_num is None:
                            for part in parts:
                                try:
                                    test_num = int(part)
                                    if 1 <= test_num <= 10000:  # Reasonable residue range
                                        res_num = test_num
                                        break
                                except ValueError:
                                    continue
                        
                        if res_num is not None:
                            if chain_id not in chains:
                                chains[chain_id] = {'min_res': res_num, 'max_res': res_num}
                            else:
                                chains[chain_id]['min_res'] = min(chains[chain_id]['min_res'], res_num)
                                chains[chain_id]['max_res'] = max(chains[chain_id]['max_res'], res_num)
                                
                except (ValueError, IndexError):
                    continue
                    
        elif line.startswith('#') and in_atom_site:
            break  # End of atom_site section
    
    return chains

def extract_chain_info_from_fasta_header(fasta_path):
    """
    Extract chain information from FASTA header to determine which chain this FASTA represents
    
    Args:
        fasta_path: Path to FASTA file
    
    Returns:
        tuple: (structure_id, target_chain_id, all_chains_in_header)
    """
    try:
        with open(fasta_path, 'r') as f:
            header_line = f.readline().strip()
    except Exception as e:
        print(f"Error reading FASTA header from {fasta_path}: {e}")
        return None, None, []
    
    # Extract structure ID from header (e.g., >2HQS_2HQS_1 -> 2HQS)
    header_parts = header_line.split()
    if len(header_parts) < 2:
        return None, None, []
    
    # Get structure ID from PDB: field or first part
    structure_id = None
    for part in header_parts:
        if part.startswith('PDB:'):
            structure_id = part.split(':')[1]
            break
        elif part.startswith('AlphaFold:AF-'):
            structure_id = part.split(':')[1]
            break
    
    if not structure_id:
        # Fallback: extract from first part of header
        first_part = header_parts[0].lstrip('>')
        if '_' in first_part:
            structure_id = first_part.split('_')[0]
        else:
            structure_id = first_part
    
    # Parse chain information from header
    all_chains = []
    target_chain = None
    
    # Look for chain information in the header
    # Handle multi-chain case: "Chains AI[auth B1], AJ[auth B2], AK[auth B7], ..." or "Chains A, E[auth F]" or "Chains A, B"
    chains_match = re.search(r'Chains\s+([^|]+)', header_line)
    if chains_match:
        chains_text = chains_match.group(1).strip()
        # Extract all chain IDs (both with and without [auth ...] notation)
        # This regex captures: A, E[auth F], NA[auth k], B, etc.
        # Updated to handle 1-2 letter chain IDs including numbers
        chain_matches = re.findall(r'([A-Z]{1,2}[A-Z0-9]*)(?:\[auth\s+[^]]+\])?', chains_text)
        if chain_matches:
            all_chains = chain_matches
    else:
        # Handle single chain case: "Chain A[auth Y]" or just "Chain A" or "Chain B"
        single_chain_match = re.search(r'Chain\s+([A-Z]{1,2}[A-Z0-9]*)(?:\[auth\s+[^]]+\])?', header_line)
        if single_chain_match:
            all_chains = [single_chain_match.group(1)]
    
    # Determine target chain based on filename
    fasta_name = Path(fasta_path).stem
    if '_chain_' in fasta_name:
        # Extract chain number from filename (e.g., "2HQS_chain_1" -> chain 1)
        chain_part = fasta_name.split('_chain_')[-1]
        if chain_part.isdigit():
            chain_index = int(chain_part) - 1  # Convert to 0-based index
            if 0 <= chain_index < len(all_chains):
                target_chain = all_chains[chain_index]
            else:
                # If index is out of range, use first chain
                target_chain = all_chains[0] if all_chains else 'A'
        else:
            target_chain = all_chains[0] if all_chains else 'A'
    else:
        # Single chain file
        target_chain = all_chains[0] if all_chains else 'A'
    
    return structure_id, target_chain, all_chains

def extract_structure_info_from_fasta_path(fasta_path):
    """
    Extract structure ID and gene name from FASTA path
    
    Args:
        fasta_path: Path like 'data/protein_structures/pal/2HQS_chain_1.fasta'
    
    Returns:
        tuple: (gene_name, structure_id, target_chain_id)
    """
    path_parts = Path(fasta_path).parts
    gene_name = path_parts[-2]  # Parent directory name
    
    # Use header parsing to get accurate chain information
    structure_id, target_chain, _ = extract_chain_info_from_fasta_header(fasta_path)
    
    if not structure_id:
        # Fallback to filename parsing
        fasta_name = Path(fasta_path).stem
        if '_chain_' in fasta_name:
            structure_id = fasta_name.split('_chain_')[0]
        else:
            structure_id = fasta_name
    
    if not target_chain:
        target_chain = 'A'  # Default fallback
    
    return gene_name, structure_id, target_chain

def find_structure_file(structure_id, structures_dir):
    """
    Find the actual structure file for a given structure ID
    
    Args:
        structure_id: Structure ID (e.g., '2HQS', 'AF-P12345')
        structures_dir: Directory containing structure files
    
    Returns:
        Path or None: Path to structure file if found
    """
    structures_path = Path(structures_dir)
    
    # Try different file extensions and compression
    possible_files = [
        f"{structure_id}.pdb",
        f"{structure_id}.pdb.gz",
        f"{structure_id}.cif",
        f"{structure_id}.cif.gz"
    ]
    
    for filename in possible_files:
        file_path = structures_path / filename
        if file_path.exists():
            return file_path
    
    return None

def create_selected_pdb_numbering_mapping(selected_paths_positive, selected_paths_negative, structures_dir, output_tsv, analysis, paramset):
    """
    Create PDB numbering mapping for selected 3D structures only
    
    Args:
        selected_paths_positive: Path to file with selected positive structures
        selected_paths_negative: Path to file with selected negative structures  
        structures_dir: Directory containing structure files
        output_tsv: Output TSV file
        analysis: Analysis identifier
        paramset: Parameter set identifier
    
    Returns:
        dict: Processing results
    """
    
    results = {
        'processed': 0,
        'failed': 0,
        'total_structures': 0,
        'skipped': 0
    }
    
    # Collect all selected FASTA paths
    selected_paths = []
    
    for paths_file in [selected_paths_positive, selected_paths_negative]:
        if Path(paths_file).exists():
            with open(paths_file, 'r') as f:
                paths = [line.strip() for line in f if line.strip()]
                selected_paths.extend(paths)
                print(f"Loaded {len(paths)} paths from {paths_file}")
    
    if not selected_paths:
        print("No selected paths found!")
        return results
    
    print(f"Processing {len(selected_paths)} selected 3D structures")
    print("=" * 50)
    
    tsv_data = []
    processed_structures = set()  # Track to avoid duplicates
    
    for fasta_path in selected_paths:
        if not fasta_path.strip():
            continue
            
        fasta_path = Path(fasta_path.strip())
        
        # Extract structure information with improved chain detection
        gene_name, structure_id, target_chain_id = extract_structure_info_from_fasta_path(fasta_path)
        
        # Check if we already processed this structure-chain combination
        structure_key = f"{structure_id}_{target_chain_id}"
        if structure_key in processed_structures:
            results['skipped'] += 1
            continue
        processed_structures.add(structure_key)
        
        # Find the actual structure file
        structure_file = find_structure_file(structure_id, structures_dir)
        if not structure_file:
            print(f"Warning: Structure file not found for {structure_id}")
            results['failed'] += 1
            continue
        
        # Check if FASTA file exists
        if not fasta_path.exists():
            print(f"Warning: FASTA file not found: {fasta_path}")
            results['failed'] += 1
            continue
        
        results['total_structures'] += 1
        
        print(f"Processing {gene_name}/{structure_id} (chain {target_chain_id})")
        
        # Extract chain information from structure file
        chain_info = extract_structure_chain_info(str(structure_file))
        
        if chain_info:
            # Get all chains from the FASTA header for this structure
            _, _, all_header_chains = extract_chain_info_from_fasta_header(fasta_path)
            
            # Find the first matching chain from header in the structure for residue numbering
            author_chain_id = None
            available_chains = list(chain_info.keys())
            
            # First round: Try to find any header chain in the structure
            for header_chain in all_header_chains:
                if header_chain in available_chains:
                    author_chain_id = header_chain
                    break
            
            # If no header chains found in structure, use first available
            if not author_chain_id and available_chains:
                author_chain_id = available_chains[0]
                print(f"  Warning: No header chains {all_header_chains} found in structure {structure_id}")
                # Limit the chain list display to avoid overwhelming output
                if len(available_chains) > 10:
                    chain_display = f"{', '.join(available_chains[:5])} ... {', '.join(available_chains[-3:])} ({len(available_chains)} total)"
                else:
                    chain_display = ', '.join(available_chains)
                print(f"  Available chains: {chain_display}")
                print(f"  Using fallback chain: {author_chain_id}")
            
            if author_chain_id:
                # Check if we got minimal residue range (indicating parsing issues)
                residue_range = chain_info[author_chain_id]['max_res'] - chain_info[author_chain_id]['min_res']
                
                # Second round: Try alternative chain search for problematic structures
                if residue_range <= 1:  # Range of 0 or 1 indicates parsing problem
                    print(f"  Minimal residue range detected ({residue_range}), trying alternative chain search...")
                    
                    # Try all available chains to find one with better residue range
                    best_chain = author_chain_id
                    best_range = residue_range
                    
                    for test_chain in available_chains:
                        test_range = chain_info[test_chain]['max_res'] - chain_info[test_chain]['min_res']
                        if test_range > best_range:
                            best_chain = test_chain
                            best_range = test_range
                            print(f"    Found better chain {test_chain} with range {test_range}")
                    
                    # If we found a better chain, use it
                    if best_range > residue_range:
                        author_chain_id = best_chain
                        print(f"  Using alternative chain {author_chain_id} with range {best_range}")
                    else:
                        # Try re-parsing the structure file with different approach
                        print(f"  Re-parsing structure file for better chain information...")
                        enhanced_chain_info = extract_enhanced_chain_info(str(structure_file))
                        
                        if enhanced_chain_info:
                            # Update chain_info with enhanced results
                            for enhanced_chain, enhanced_data in enhanced_chain_info.items():
                                enhanced_range = enhanced_data['max_res'] - enhanced_data['min_res']
                                if enhanced_range > best_range:
                                    chain_info[enhanced_chain] = enhanced_data
                                    author_chain_id = enhanced_chain
                                    best_range = enhanced_range
                                    print(f"    Enhanced parsing found chain {enhanced_chain} with range {enhanced_range}")
                
                results['processed'] += 1
                
                # Determine structure type
                structure_type = "computed" if structure_id.startswith("AF-") else "experimental"
                
                # Display all chains from header, use first found chain for numbering
                display_chains = ','.join(all_header_chains) if all_header_chains else author_chain_id
                
                tsv_data.append({
                    'Gene': gene_name,
                    'Structure_ID': structure_id,
                    'Chain': display_chains,
                    'Start': chain_info[author_chain_id]['min_res'],
                    'End': chain_info[author_chain_id]['max_res'],
                    'Range': f"{chain_info[author_chain_id]['min_res']}-{chain_info[author_chain_id]['max_res']}",
                    'FASTA_Path': str(fasta_path),
                    'Structure_Path': str(structure_file),
                    'Structure_Type': structure_type
                })
            else:
                print(f"  No chains found in structure!")
                results['failed'] += 1
                continue
        else:
            results['failed'] += 1
            print(f"  Failed to extract chain info for {structure_file}")
    
    # Save results
    if tsv_data:
        df = pd.DataFrame(tsv_data)
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"\nCreated selected PDB numbering TSV: {output_tsv}")
        print(f"Total entries: {len(tsv_data)}")
    else:
        print("No data to write to TSV file")
        # Create empty TSV with headers
        df = pd.DataFrame(columns=['Gene', 'Structure_ID', 'Chain', 'Start', 'End', 'Range', 'FASTA_Path', 'Structure_Path', 'Structure_Type'])
        df.to_csv(output_tsv, sep='\t', index=False)
    
    print("=" * 50)
    print(f"Selected structure processing complete!")
    print(f"Successfully processed: {results['processed']}")
    print(f"Failed: {results['failed']}")
    print(f"Skipped duplicates: {results['skipped']}")
    print(f"Total structures: {results['total_structures']}")
    
    return results

def create_enhanced_pdb_numbering_mapping(mapping_file, selected_paths_positive, selected_paths_negative, structures_dir, output_tsv, analysis, paramset):
    """
    Create PDB numbering mapping using the FASTA structure mapping file combined with selected paths
    
    Args:
        mapping_file: Path to the FASTA structure mapping TSV file
        selected_paths_positive: Path to file with selected positive structures
        selected_paths_negative: Path to file with selected negative structures  
        structures_dir: Directory containing structure files
        output_tsv: Output TSV file
        analysis: Analysis identifier
        paramset: Parameter set identifier
    
    Returns:
        dict: Processing results
    """
    
    results = {
        'processed': 0,
        'failed': 0,
        'total_structures': 0,
        'skipped': 0
    }
    
    # Load the FASTA structure mapping
    print(f"Loading FASTA structure mapping from: {mapping_file}")
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        print(f"Loaded {len(mapping_df)} entries from mapping file")
    except Exception as e:
        print(f"Error reading mapping file {mapping_file}: {e}")
        return results
    
    # Collect all selected FASTA paths
    selected_paths = []
    
    for paths_file in [selected_paths_positive, selected_paths_negative]:
        if Path(paths_file).exists():
            with open(paths_file, 'r') as f:
                paths = [line.strip() for line in f if line.strip()]
                selected_paths.extend(paths)
                print(f"Loaded {len(paths)} paths from {paths_file}")
    
    if not selected_paths:
        print("No selected paths found!")
        return results
    
    print(f"Processing {len(selected_paths)} selected 3D structures")
    print("=" * 50)
    
    # Filter mapping to only include selected FASTA files
    selected_mapping_entries = []
    
    for fasta_path in selected_paths:
        fasta_path = Path(fasta_path.strip())
        
        # Find matching entry in mapping file based on FASTA path
        matching_entries = mapping_df[mapping_df['fasta_path'] == str(fasta_path)]
        
        if len(matching_entries) > 0:
            selected_mapping_entries.extend(matching_entries.to_dict('records'))
        else:
            print(f"Warning: No mapping found for FASTA path: {fasta_path}")
    
    if not selected_mapping_entries:
        print("No mapping entries found for selected FASTA files!")
        return results
    
    print(f"Found {len(selected_mapping_entries)} matching mapping entries")
    
    tsv_data = []
    processed_structures = set()  # Track to avoid duplicates
    
    for entry in selected_mapping_entries:
        structure_id = entry['structure_id']
        structure_type = entry['structure_type']
        gene_name = entry['gene_name']
        protein_name = entry['protein_name']
        chain_info = entry['chain']
        species = entry['species']
        structure_path = entry['structure_path']
        fasta_path = entry['fasta_path']
        
        # Check if we already processed this structure-chain combination
        structure_key = f"{structure_id}_{chain_info}"
        if structure_key in processed_structures:
            results['skipped'] += 1
            continue
        processed_structures.add(structure_key)
        
        # Verify structure file exists
        structure_file = Path(structure_path)
        if not structure_file.exists():
            print(f"Warning: Structure file not found: {structure_file}")
            results['failed'] += 1
            continue
        
        # Verify FASTA file exists
        fasta_file = Path(fasta_path)
        if not fasta_file.exists():
            print(f"Warning: FASTA file not found: {fasta_file}")
            results['failed'] += 1
            continue
        
        results['total_structures'] += 1
        
        print(f"Processing {gene_name}/{structure_id} (chain {chain_info})")
        
        # Extract chain information from structure file
        chain_numbering = extract_structure_chain_info(str(structure_file))
        
        if chain_numbering:
            # Parse the chain info from mapping (could be multiple chains)
            chain_ids = [c.strip() for c in chain_info.split(',')]
            
            # Find the first available chain for numbering
            author_chain_id = None
            for chain_id in chain_ids:
                if chain_id in chain_numbering:
                    author_chain_id = chain_id
                    break
            
            # If no specific chain found, use first available
            if not author_chain_id and chain_numbering:
                available_chains = list(chain_numbering.keys())
                author_chain_id = available_chains[0]
                print(f"  Warning: Chains {chain_ids} not found in structure {structure_id}")
                print(f"  Available chains: {', '.join(available_chains)}")
                print(f"  Using fallback chain: {author_chain_id}")
            
            if author_chain_id:
                # Check if we got minimal residue range (indicating parsing issues)
                residue_range = chain_numbering[author_chain_id]['max_res'] - chain_numbering[author_chain_id]['min_res']
                
                # Try enhanced parsing for problematic structures
                if residue_range <= 1:
                    print(f"  Minimal residue range detected ({residue_range}), trying enhanced parsing...")
                    enhanced_chain_info = extract_enhanced_chain_info(str(structure_file))
                    
                    if enhanced_chain_info:
                        # Update with enhanced results if better
                        for enhanced_chain, enhanced_data in enhanced_chain_info.items():
                            enhanced_range = enhanced_data['max_res'] - enhanced_data['min_res']
                            if enhanced_range > residue_range:
                                chain_numbering[enhanced_chain] = enhanced_data
                                author_chain_id = enhanced_chain
                                residue_range = enhanced_range
                                print(f"    Enhanced parsing found chain {enhanced_chain} with range {enhanced_range}")
                
                results['processed'] += 1
                
                tsv_data.append({
                    'Gene': gene_name,
                    'Structure_ID': structure_id,
                    'Chain': chain_info,  # Use original chain info from mapping
                    'Start': chain_numbering[author_chain_id]['min_res'],
                    'End': chain_numbering[author_chain_id]['max_res'],
                    'Range': f"{chain_numbering[author_chain_id]['min_res']}-{chain_numbering[author_chain_id]['max_res']}",
                    'FASTA_Path': str(fasta_file),
                    'Structure_Path': str(structure_file),
                    'Structure_Type': structure_type
                })
            else:
                print(f"  No chains found in structure!")
                results['failed'] += 1
                continue
        else:
            results['failed'] += 1
            print(f"  Failed to extract chain info for {structure_file}")
    
    # Save results
    if tsv_data:
        df = pd.DataFrame(tsv_data)
        df.to_csv(output_tsv, sep='\t', index=False)
        print(f"\nCreated enhanced PDB numbering TSV: {output_tsv}")
        print(f"Total entries: {len(tsv_data)}")
    else:
        print("No data to write to TSV file")
        # Create empty TSV with headers
        df = pd.DataFrame(columns=['Gene', 'Structure_ID', 'Chain', 'Start', 'End', 'Range', 'FASTA_Path', 'Structure_Path', 'Structure_Type'])
        df.to_csv(output_tsv, sep='\t', index=False)
    
    print("=" * 50)
    print(f"Enhanced structure processing complete!")
    print(f"Successfully processed: {results['processed']}")
    print(f"Failed: {results['failed']}")
    print(f"Skipped duplicates: {results['skipped']}")
    print(f"Total structures: {results['total_structures']}")
    
    return results

def main():
    """Main function for Snakemake integration"""
    
    try:
        # Get inputs from Snakemake
        fasta_structure_mapping = snakemake.input.fasta_structure_mapping
        selected_paths_positive = snakemake.input.selected_3d_paths_positive
        selected_paths_negative = snakemake.input.selected_3d_paths_negative
        numbering_tsv = snakemake.output.numbering_tsv
        sentinel_file = snakemake.output.numbering_sentinel
        
        # Get parameters
        analysis = snakemake.params.analysis
        paramset = snakemake.params.paramset
        structures_dir = snakemake.params.structures_dir
        
        print(f"Selected PDB numbering extraction for {analysis}_{paramset}")
        print(f"FASTA structure mapping: {fasta_structure_mapping}")
        print(f"Positive paths: {selected_paths_positive}")
        print(f"Negative paths: {selected_paths_negative}")
        print(f"Structures directory: {structures_dir}")
        print(f"Output TSV: {numbering_tsv}")
        
    except NameError:
        # Test mode - command line arguments
        if len(sys.argv) != 7:
            print("Usage: python extract_pdb_numbering_selected.py <analysis> <paramset> <mapping_file> <positive_paths> <negative_paths> <structures_dir>")
            sys.exit(1)
        
        analysis = sys.argv[1]
        paramset = sys.argv[2]
        fasta_structure_mapping = sys.argv[3]
        selected_paths_positive = sys.argv[4]
        selected_paths_negative = sys.argv[5]
        structures_dir = sys.argv[6]
        
        numbering_tsv = f"results/{analysis}_{paramset}/protein_analysis/pdb_numbering_mapping.tsv"
        sentinel_file = f"results/{analysis}_{paramset}/protein_analysis/pdb_numbering_complete.sentinel"
    
    # Create output directory
    Path(numbering_tsv).parent.mkdir(parents=True, exist_ok=True)
    
    # Create the selected mapping using the FASTA structure mapping file
    results = create_enhanced_pdb_numbering_mapping(
        fasta_structure_mapping, selected_paths_positive, selected_paths_negative, 
        structures_dir, numbering_tsv, analysis, paramset
    )
    
    # Create sentinel file
    with open(sentinel_file, 'w') as f:
        f.write(f"Selected PDB numbering extraction completed for {analysis}_{paramset}\n")
        f.write(f"Successfully processed: {results.get('processed', 0)}\n")
        f.write(f"Failed: {results.get('failed', 0)}\n")
        f.write(f"Total structures: {results.get('total_structures', 0)}\n")
        f.write(f"Output TSV: {numbering_tsv}\n")
    
    print(f"Selected PDB numbering extraction completed. Sentinel: {sentinel_file}")

if __name__ == "__main__":
    main()