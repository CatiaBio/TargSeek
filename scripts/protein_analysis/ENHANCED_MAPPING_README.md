# Enhanced PDB Mapping Tools

This directory contains enhanced tools for accurate PDB residue mapping using PDBML/XML format.

## New Enhanced Tools

### 1. `parse_pdbml_residue_numbering.py`
- **Purpose**: Direct parsing of PDBML/XML files for accurate residue numbering
- **Key Feature**: Extracts `auth_seq_id` values which represent original author numbering
- **Usage**: `python parse_pdbml_residue_numbering.py <pdbml_file.xml.gz> [--chain A]`

### 2. `map_fasta_to_pdb_residues_enhanced.py`
- **Purpose**: Enhanced epitope mapping using PDBML when available
- **Key Feature**: Prefers PDBML format over legacy PDB for accuracy
- **Usage**: `python map_fasta_to_pdb_residues_enhanced.py <bepipred_dir> [--prefer-pdbml]`
- **Output**: Files with `_enhanced_pdbml_mapping.tsv` suffix

## Why PDBML is Better

1. **Accurate Numbering**: Contains `auth_seq_id` which preserves original sequence numbering
2. **No Signal Peptide Issues**: Handles signal peptide removal correctly
3. **Insertion Codes**: Properly manages complex numbering schemes
4. **Example**: Position 423 in PDBML corresponds directly to author's intended numbering

## Migration from Old Tools

### Replaced Scripts:
- ❌ `map_fasta_to_pdb_residues.py` (old version)
- ❌ `analyze_bepipred_input_sequences.py` (diagnostic tool)
- ❌ `extract_pdb_residue_numbering.py` (basic version)

### Enhanced Scripts:
- ✅ `parse_pdbml_residue_numbering.py` (new)
- ✅ `map_fasta_to_pdb_residues_enhanced.py` (new)

## File Format Changes

### 3D Structure Downloads
The download script now automatically downloads both formats:
- `{pdb_id}.pdb.gz` (legacy format)
- `{pdb_id}.xml.gz` (PDBML format with accurate numbering)

### Enhanced Output Format
New mapping files include:
- `PDB_Start`, `PDB_End`: Accurate PDB residue numbers from PDBML
- `PDB_Residues`: Comma-separated list of all PDB residue numbers
- `Mapping_Method`: Indicates whether PDBML or legacy PDB was used

## Example PDBML Structure
```xml
<PDBx:atom_site id="3">
    <PDBx:auth_asym_id>A</PDBx:auth_asym_id>
    <PDBx:auth_atom_id>CA</PDBx:auth_atom_id>
    <PDBx:auth_comp_id>THR</PDBx:auth_comp_id>
    <PDBx:auth_seq_id>423</PDBx:auth_seq_id>
</PDBx:atom_site>
```

In this example, the `auth_seq_id` value of 423 represents the exact position numbering intended by the structure authors.