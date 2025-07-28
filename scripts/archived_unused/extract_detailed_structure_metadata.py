#!/usr/bin/env python3
"""
Detailed Structure Metadata Extractor
=====================================

Extracts comprehensive metadata for all protein structures associated with genes,
including experimental methods, resolution, chain information, and entity details.

This script provides detailed JSON information about structures per gene including:
- Method (X-RAY DIFFRACTION, NMR, CRYO-EM, etc.)
- Resolution (for applicable methods)
- Entity information (ID, molecule name, chains, sequence length)
- Organism details
- Additional structural and experimental details
"""

import requests
import pandas as pd
from pathlib import Path
import logging
import json
import time
from datetime import datetime
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def get_comprehensive_pdb_metadata(pdb_id):
    """
    Get comprehensive metadata for a PDB structure including entity details
    
    Args:
        pdb_id (str): PDB identifier
    
    Returns:
        dict: Comprehensive metadata dictionary
    """
    try:
        # Get entry-level metadata
        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        entry_response = requests.get(entry_url, timeout=30)
        entry_response.raise_for_status()
        entry_data = entry_response.json()
        
        # Get polymer entity information using different API endpoints
        polymer_entities = []
        
        # Try multiple approaches to get entity information
        # Approach 1: Direct polymer entity endpoint
        try:
            polymer_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}"
            polymer_response = requests.get(polymer_url, timeout=30)
            if polymer_response.status_code == 200:
                polymer_data = polymer_response.json()
                if isinstance(polymer_data, list):
                    polymer_entities.extend(polymer_data)
                elif polymer_data:
                    polymer_entities.append(polymer_data)
        except:
            pass
        
        # Approach 2: Get entity information from polymer entity instances
        if not polymer_entities:
            try:
                instances_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}"
                instances_response = requests.get(instances_url, timeout=30)
                if instances_response.status_code == 200:
                    instances_data = instances_response.json()
                    if isinstance(instances_data, list):
                        # Group instances by entity
                        entity_groups = {}
                        for instance in instances_data:
                            entity_id = instance.get("rcsb_polymer_entity_instance_container_identifiers", {}).get("entity_id")
                            if entity_id:
                                if entity_id not in entity_groups:
                                    entity_groups[entity_id] = instance
                        polymer_entities = list(entity_groups.values())
            except:
                pass
        
        # Get assembly information
        assembly_url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}"
        assembly_response = requests.get(assembly_url, timeout=30)
        assembly_data = []
        if assembly_response.status_code == 200:
            assembly_info = assembly_response.json()
            if isinstance(assembly_info, list):
                assembly_data = assembly_info
            elif assembly_info:
                assembly_data = [assembly_info]
        
        # Extract basic structure information
        metadata = {
            "pdb_id": pdb_id,
            "title": entry_data.get("struct", {}).get("title", "Unknown"),
            "classification": entry_data.get("struct_keywords", {}).get("pdbx_keywords", "Unknown"),
            "experimental_method": [],
            "resolution": None,
            "resolution_unit": "Å",
            "deposit_date": entry_data.get("rcsb_accession_info", {}).get("deposit_date", "Unknown"),
            "release_date": entry_data.get("rcsb_accession_info", {}).get("initial_release_date", "Unknown"),
            "revision_date": entry_data.get("rcsb_accession_info", {}).get("revision_date", "Unknown"),
            "authors": [],
            "journal": "Unknown",
            "doi": None,
            "space_group": entry_data.get("cell", {}).get("pdbx_space_group_name_H-M", "Unknown"),
            "crystal_system": None,
            "entities": [],
            "assembly_info": [],
            "total_chains": 0,
            "unique_chains": 0
        }
        
        # Get experimental methods
        exptl_methods = entry_data.get("exptl", [])
        for method in exptl_methods:
            method_name = method.get("method", "Unknown")
            metadata["experimental_method"].append(method_name)
        
        # Get resolution for applicable methods
        refine_info = entry_data.get("refine", [])
        if refine_info and any("DIFFRACTION" in method or "SCATTERING" in method for method in metadata["experimental_method"]):
            resolution = refine_info[0].get("ls_d_res_high")
            if resolution:
                metadata["resolution"] = float(resolution)
        
        # Also try getting resolution from rcsb_entry_info
        rcsb_entry_info = entry_data.get("rcsb_entry_info", {})
        if not metadata["resolution"] and rcsb_entry_info:
            resolution_combined = rcsb_entry_info.get("resolution_combined", [])
            if resolution_combined:
                metadata["resolution"] = float(resolution_combined[0])
        
        # Try diffraction resolution
        if not metadata["resolution"]:
            diffrn_detector = entry_data.get("diffrn_detector", [])
            for detector in diffrn_detector:
                resolution = detector.get("diffrn_radiation_wavelength")
                if resolution:
                    try:
                        metadata["resolution"] = float(resolution)
                        break
                    except:
                        continue
        
        # For electron microscopy, check different resolution fields
        if any("ELECTRON" in method for method in metadata["experimental_method"]):
            em_3d_reconstruction = entry_data.get("em_3d_reconstruction", [])
            if em_3d_reconstruction:
                em_resolution = em_3d_reconstruction[0].get("resolution")
                if em_resolution:
                    metadata["resolution"] = float(em_resolution)
        
        # Get crystal system
        symmetry = entry_data.get("symmetry", {})
        if symmetry:
            metadata["crystal_system"] = symmetry.get("space_group_name_H-M", "Unknown")
        
        # Get authors
        audit_author = entry_data.get("audit_author", [])
        metadata["authors"] = [author.get("name", "") for author in audit_author[:10]]  # Top 10 authors
        
        # Get journal and DOI information
        citation_info = entry_data.get("citation", [])
        if citation_info:
            citation = citation_info[0]
            metadata["journal"] = citation.get("journal_abbrev", "Unknown")
            metadata["doi"] = citation.get("pdbx_database_id_DOI")
        
        # Process polymer entities (protein chains)
        chain_count = 0
        unique_sequences = set()
        
        for entity in polymer_entities:
            entity_info = {
                "entity_id": entity.get("rcsb_id", "Unknown"),
                "entity_type": entity.get("entity_poly", {}).get("type", "Unknown"),
                "molecule_name": None,
                "chains": [],
                "sequence_length": 0,
                "molecular_weight": None,
                "organism_scientific": [],
                "organism_common": [],
                "gene_names": [],
                "uniprot_accessions": [],
                "sequence": None,
                "mutation_details": []
            }
            
            # Get molecule name from different possible sources
            polymer_comp = entity.get("entity_poly", {})
            if polymer_comp:
                entity_info["molecule_name"] = polymer_comp.get("pdbx_strand_id", "Unknown")
            
            # Try to get better molecule name from rcsb_polymer_entity
            rcsb_polymer = entity.get("rcsb_polymer_entity", {})
            if rcsb_polymer:
                formula_weight = rcsb_polymer.get("formula_weight")
                if formula_weight:
                    entity_info["molecular_weight"] = float(formula_weight)
            
            # Get chain information
            polymer_entity_instances = entity.get("rcsb_polymer_entity_container_identifiers", {})
            if polymer_entity_instances:
                auth_asym_ids = polymer_entity_instances.get("auth_asym_ids", [])
                entity_info["chains"] = auth_asym_ids
                chain_count += len(auth_asym_ids)
            
            # Get sequence information
            entity_poly = entity.get("entity_poly", {})
            if entity_poly:
                sequence = entity_poly.get("pdbx_seq_one_letter_code_can")
                if sequence:
                    entity_info["sequence"] = sequence.replace('\n', '').replace(' ', '')
                    entity_info["sequence_length"] = len(entity_info["sequence"])
                    unique_sequences.add(entity_info["sequence"])
            
            # Get organism information
            rcsb_entity_source_organism = entity.get("rcsb_entity_source_organism", [])
            for organism in rcsb_entity_source_organism:
                sci_name = organism.get("ncbi_scientific_name")
                common_name = organism.get("ncbi_common_name")
                if sci_name:
                    entity_info["organism_scientific"].append(sci_name)
                if common_name:
                    entity_info["organism_common"].append(common_name)
            
            # Get gene names and UniProt accessions
            rcsb_polymer_entity_annotation = entity.get("rcsb_polymer_entity_annotation", [])
            for annotation in rcsb_polymer_entity_annotation:
                if annotation.get("type") == "GO":
                    continue  # Skip GO annotations for now
                    
                annotation_id = annotation.get("annotation_id")
                if annotation_id and annotation_id.startswith("UniProt"):
                    entity_info["uniprot_accessions"].append(annotation_id.replace("UniProt:", ""))
            
            # Try to get gene names from entity_src_gen
            entity_src_gen = entity.get("entity_src_gen", [])
            for src in entity_src_gen:
                gene_name = src.get("gene_src_common_name")
                if gene_name:
                    entity_info["gene_names"].append(gene_name)
            
            metadata["entities"].append(entity_info)
        
        # Update chain counts
        metadata["total_chains"] = chain_count
        metadata["unique_chains"] = len(unique_sequences)
        
        # Process assembly information
        for assembly in assembly_data:
            assembly_info = {
                "assembly_id": assembly.get("rcsb_id", "Unknown"),
                "method_details": assembly.get("rcsb_assembly_info", {}).get("assembly_method_details", "Unknown"),
                "oligomeric_details": assembly.get("rcsb_assembly_info", {}).get("oligomeric_details", []),
                "polymer_entity_count": assembly.get("rcsb_assembly_info", {}).get("polymer_entity_count", 0),
                "polymer_entity_instance_count": assembly.get("rcsb_assembly_info", {}).get("polymer_entity_instance_count", 0)
            }
            metadata["assembly_info"].append(assembly_info)
        
        return metadata
        
    except Exception as e:
        logging.warning(f"Failed to get comprehensive metadata for PDB {pdb_id}: {e}")
        return {
            "pdb_id": pdb_id,
            "error": str(e),
            "title": "Unknown",
            "experimental_method": ["Unknown"],
            "resolution": None,
            "entities": [],
            "assembly_info": []
        }

def get_alphafold_comprehensive_metadata(alphafold_id):
    """
    Get comprehensive metadata for an AlphaFold model
    
    Args:
        alphafold_id (str): AlphaFold identifier (e.g., AF-P12345)
    
    Returns:
        dict: Comprehensive metadata dictionary
    """
    uniprot_accession = alphafold_id.replace('AF-', '')
    
    try:
        # Get AlphaFold metadata
        af_api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
        response = requests.get(af_api_url, timeout=30)
        response.raise_for_status()
        
        af_data = response.json()
        if not af_data:
            raise ValueError("No AlphaFold data found")
        
        af_entry = af_data[0] if isinstance(af_data, list) else af_data
        
        # Get additional UniProt metadata
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.json"
        uniprot_response = requests.get(uniprot_url, timeout=30)
        uniprot_data = {}
        if uniprot_response.status_code == 200:
            uniprot_data = uniprot_response.json()
        
        metadata = {
            "alphafold_id": alphafold_id,
            "uniprot_accession": uniprot_accession,
            "model_type": "AlphaFold Computed Model",
            "experimental_method": ["Computational Prediction"],
            "resolution": "N/A (Computed Model)",
            "model_confidence": af_entry.get('globalMetricValue', 'Unknown'),
            "model_version": af_entry.get('modelVersion', 'Unknown'),
            "sequence_length": af_entry.get('sequenceLength', 'Unknown'),
            "created_date": af_entry.get('modelCreatedDate', 'Unknown'),
            "gene_name": af_entry.get('gene', 'Unknown'),
            "organism_scientific": af_entry.get('organismScientificName', 'Unknown'),
            "protein_name": "Unknown",
            "entities": [{
                "entity_id": "1",
                "entity_type": "protein",
                "molecule_name": "Unknown",
                "chains": ["A"],
                "sequence_length": af_entry.get('sequenceLength', 0),
                "organism_scientific": [af_entry.get('organismScientificName', 'Unknown')],
                "gene_names": [af_entry.get('gene', 'Unknown')] if af_entry.get('gene') else [],
                "uniprot_accessions": [uniprot_accession],
                "confidence_score": af_entry.get('globalMetricValue', 'Unknown'),
                "model_version": af_entry.get('modelVersion', 'Unknown')
            }],
            "total_chains": 1,
            "unique_chains": 1
        }
        
        # Add protein name from UniProt if available
        if uniprot_data:
            protein_desc = uniprot_data.get('proteinDescription', {})
            recommended_name = protein_desc.get('recommendedName', {})
            if recommended_name:
                full_name = recommended_name.get('fullName', {})
                if full_name:
                    protein_name = full_name.get('value', 'Unknown')
                    metadata["protein_name"] = protein_name
                    metadata["entities"][0]["molecule_name"] = protein_name
        
        return metadata
        
    except Exception as e:
        logging.warning(f"Failed to get comprehensive AlphaFold metadata for {alphafold_id}: {e}")
        return {
            "alphafold_id": alphafold_id,
            "uniprot_accession": uniprot_accession,
            "model_type": "AlphaFold Computed Model",
            "experimental_method": ["Computational Prediction"],
            "resolution": "N/A (Computed Model)",
            "error": str(e),
            "entities": []
        }

def extract_structures_from_gene_directory(gene_dir):
    """
    Extract structure IDs from a gene directory
    
    Args:
        gene_dir (Path): Path to gene directory
    
    Returns:
        tuple: (pdb_ids, alphafold_ids)
    """
    pdb_ids = []
    alphafold_ids = []
    
    # Look for PDB structure files
    pdb_files = list(gene_dir.glob("*.pdb.gz"))
    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem.replace('.pdb', '')
        if not pdb_id.startswith('AF-'):
            pdb_ids.append(pdb_id)
    
    # Look for AlphaFold model files
    af_files = list(gene_dir.glob("AF-*.pdb"))
    for af_file in af_files:
        alphafold_id = af_file.stem
        alphafold_ids.append(alphafold_id)
    
    return pdb_ids, alphafold_ids

def process_gene_structures(gene_name, structures_base_dir):
    """
    Process all structures for a specific gene
    
    Args:
        gene_name (str): Gene name
        structures_base_dir (Path): Base directory containing gene structure folders
    
    Returns:
        dict: Gene structure metadata
    """
    gene_dir = structures_base_dir / gene_name
    
    if not gene_dir.exists():
        logging.warning(f"Gene directory not found: {gene_dir}")
        return {
            "gene_name": gene_name,
            "status": "directory_not_found",
            "structures": []
        }
    
    # Extract structure IDs
    pdb_ids, alphafold_ids = extract_structures_from_gene_directory(gene_dir)
    
    if not pdb_ids and not alphafold_ids:
        logging.info(f"No structures found for gene {gene_name}")
        return {
            "gene_name": gene_name,
            "status": "no_structures_found",
            "structures": []
        }
    
    logging.info(f"Processing {len(pdb_ids)} experimental + {len(alphafold_ids)} computed structures for gene {gene_name}")
    
    structures = []
    
    # Process experimental structures
    for pdb_id in pdb_ids:
        logging.info(f"Getting comprehensive metadata for PDB {pdb_id}")
        metadata = get_comprehensive_pdb_metadata(pdb_id)
        metadata["structure_type"] = "experimental"
        structures.append(metadata)
        time.sleep(0.5)  # Rate limiting
    
    # Process computed models
    for alphafold_id in alphafold_ids:
        logging.info(f"Getting comprehensive metadata for AlphaFold model {alphafold_id}")
        metadata = get_alphafold_comprehensive_metadata(alphafold_id)
        metadata["structure_type"] = "computed_model"
        structures.append(metadata)
        time.sleep(0.5)  # Rate limiting
    
    return {
        "gene_name": gene_name,
        "status": "completed",
        "experimental_structures_count": len(pdb_ids),
        "computed_models_count": len(alphafold_ids),
        "total_structures_count": len(structures),
        "structures": structures,
        "processing_date": datetime.now().isoformat()
    }

def main():
    """Main function for comprehensive structure metadata extraction"""
    
    parser = argparse.ArgumentParser(description='Extract comprehensive structure metadata for genes')
    parser.add_argument('--structures_dir', type=str, default='data/protein_structures',
                        help='Directory containing gene structure folders')
    parser.add_argument('--output_file', type=str, default='comprehensive_structure_metadata.json',
                        help='Output JSON file for comprehensive metadata')
    parser.add_argument('--gene', type=str, help='Process specific gene only')
    parser.add_argument('--summary_file', type=str, help='Summary TSV file with gene list')
    
    args = parser.parse_args()
    
    structures_dir = Path(args.structures_dir)
    output_file = Path(args.output_file)
    
    if not structures_dir.exists():
        logging.error(f"Structures directory not found: {structures_dir}")
        return
    
    # Determine which genes to process
    genes_to_process = []
    
    if args.gene:
        # Process single gene
        genes_to_process = [args.gene]
    elif args.summary_file:
        # Process genes from summary file
        summary_file = Path(args.summary_file)
        if summary_file.exists():
            df = pd.read_csv(summary_file, sep='\t')
            genes_to_process = df['gene_name'].unique().tolist()
            logging.info(f"Processing {len(genes_to_process)} genes from {summary_file}")
        else:
            logging.error(f"Summary file not found: {summary_file}")
            return
    else:
        # Process all gene directories
        gene_dirs = [d for d in structures_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
        genes_to_process = [d.name for d in gene_dirs]
        logging.info(f"Processing {len(genes_to_process)} genes from structure directories")
    
    if not genes_to_process:
        logging.error("No genes to process")
        return
    
    logging.info(f"Extracting comprehensive metadata for {len(genes_to_process)} genes")
    
    # Process all genes
    all_gene_data = {}
    
    for i, gene_name in enumerate(genes_to_process, 1):
        logging.info(f"Progress: {i}/{len(genes_to_process)} - Processing gene: {gene_name}")
        
        gene_data = process_gene_structures(gene_name, structures_dir)
        all_gene_data[gene_name] = gene_data
        
        # Save individual gene metadata
        gene_output_file = structures_dir / gene_name / "comprehensive_metadata.json"
        if gene_output_file.parent.exists():
            with open(gene_output_file, 'w') as f:
                json.dump(gene_data, f, indent=2)
            logging.info(f"Saved individual metadata: {gene_output_file}")
    
    # Create summary statistics
    total_genes = len(all_gene_data)
    genes_with_structures = len([g for g, data in all_gene_data.items() if data['status'] == 'completed'])
    total_experimental = sum(data.get('experimental_structures_count', 0) for data in all_gene_data.values())
    total_computed = sum(data.get('computed_models_count', 0) for data in all_gene_data.values())
    total_structures = sum(data.get('total_structures_count', 0) for data in all_gene_data.values())
    
    # Create final output
    final_output = {
        "extraction_info": {
            "extraction_date": datetime.now().isoformat(),
            "structures_directory": str(structures_dir),
            "total_genes_processed": total_genes,
            "genes_with_structures": genes_with_structures,
            "total_experimental_structures": total_experimental,
            "total_computed_models": total_computed,
            "total_structures": total_structures
        },
        "genes": all_gene_data
    }
    
    # Save comprehensive metadata
    with open(output_file, 'w') as f:
        json.dump(final_output, f, indent=2)
    
    logging.info(f"=== Comprehensive Metadata Extraction Complete ===")
    logging.info(f"Total genes processed: {total_genes}")
    logging.info(f"Genes with structures: {genes_with_structures}")
    logging.info(f"Total experimental structures: {total_experimental}")
    logging.info(f"Total computed models: {total_computed}")
    logging.info(f"Total structures: {total_structures}")
    logging.info(f"Comprehensive metadata saved to: {output_file}")

if __name__ == "__main__":
    main()