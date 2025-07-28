# Archived Unused Scripts

This directory contains scripts that are no longer used in the current TargSeek pipeline but are preserved for reference.

## Gene Selection Scripts (Download Pipeline)
- `extract_detailed_structure_metadata.py` - Extracted detailed metadata from protein structures (replaced by simplified version)
- `extract_simple_structure_metadata.py` - Simplified metadata extraction (functionality integrated into main pipeline)

## Protein Analysis Scripts (Analysis Pipeline)  
- `create_msa_fasta.py` - Created FASTA files for MSA (functionality integrated into other scripts)
- `map_fasta_to_pdb_residues_enhanced.py` - Enhanced PDB residue mapping (replaced by used_structures_mapping approach)
- `select_msa_proteins.py` - MSA protein selection (functionality integrated into sequence preparation)

## Epitope Analysis Utilities (Root Directory)
- `create_structure_filtered_epitopes.py` - Standalone utility to filter BepiPred epitope predictions based on 3D structure ranges
- `create_targeted_epitope_reports.py` - Standalone utility to create epitope reports for only the used structures  
- `debug_linear_epitopes.py` - Debugging tool for linear epitope detection analysis
- `generate_analysis_summary_report.py` - Standalone utility to generate summary reports of analyzed vs non-analyzed genes
- `test_bepipred_direct.py` - Testing tool for BepiPred direct execution

## Development and Testing Scripts
- `debug_uniprot_structure.py` - Debugging tool for UniProt structure API
- `gene_taxa_coverage_fast_test.py` - Fast testing version of coverage analysis  
- `select_proteins_to_study_original.py` - Original version of protein selection script
- `test_uniprot_api.py` - Testing tool for UniProt API functionality

## Archive Date
Scripts archived on: July 28, 2025

## Notes
These scripts may contain useful code patterns or approaches that could be referenced for future development, but they are not part of the active pipeline workflow.