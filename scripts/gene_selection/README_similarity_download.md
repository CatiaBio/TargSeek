# Sequence Similarity-Based Structure Downloader

A simplified script that downloads protein structures based on sequence similarity using the PDB REST API.

## Features

- **Smart Reference Selection**: Automatically selects reference sequences from preferred organisms (E. coli, Bacillus, Salmonella, Pseudomonas)
- **PDB API Integration**: Uses the official PDB REST API for fast and accurate sequence similarity searches  
- **Experimental Structure Priority**: Prioritizes experimental structures over computed models
- **AlphaFold Fallback**: Includes AlphaFold structures if insufficient experimental structures are found
- **Concurrent Downloads**: Fast parallel downloading of structure files
- **Multiple Formats**: Downloads structures in PDB.gz, CIF.gz, or PDB formats

## Usage

### Basic Usage
```bash
python scripts/gene_selection/download_structures_by_similarity.py bamA
```

### Advanced Options
```bash
python scripts/gene_selection/download_structures_by_similarity.py bamA \
    --max-structures 10 \
    --identity-cutoff 0.4 \
    --evalue-cutoff 0.01 \
    --min-experimental 5
```

### Parameters

- `gene_name`: Gene name to search for (e.g., bamA, cls, lptE)
- `--max-structures`: Maximum number of structures to download (default: 20)
- `--identity-cutoff`: Minimum sequence identity (0.0-1.0, default: 0.3 = 30%)
- `--evalue-cutoff`: Maximum E-value cutoff (default: 0.1)
- `--min-experimental`: Minimum experimental structures before including AlphaFold (default: 5)
- `--no-alphafold`: Skip AlphaFold structures entirely
- `--output-dir`: Custom output directory

## How It Works

1. **Reference Selection**: Finds a reference sequence from your protein sequence data:
   - Prioritizes E. coli, Bacillus, Salmonella, Pseudomonas
   - Falls back to any available species
   - Uses first sequence from multi-sequence FASTA files

2. **PDB Search**: Searches PDB using sequence similarity:
   - Combines gene name search + sequence similarity
   - Filters by identity and E-value thresholds
   - Prioritizes experimental structures

3. **Structure Download**: Downloads matching structures:
   - Experimental structures first (PDB.gz, CIF.gz formats)
   - AlphaFold models if needed (PDB format)
   - Concurrent downloads for speed

4. **Results**: Creates organized output:
   - `{gene_name}/` directory with all structures
   - `download_metadata.json` with search details
   - Proper structure file naming

## Example Output

For `bamA` with identity cutoff 0.3:

```
data/protein_structures_similarity/bamA/
├── 6V05.pdb.gz          # Experimental structure
├── 7YE4.pdb.gz          # Experimental structure  
├── 7YE6.pdb.gz          # Experimental structure
├── 8BNZ.pdb.gz          # Experimental structure
├── 9A2O.cif.gz          # Experimental structure
└── download_metadata.json
```

## Metadata Format

The `download_metadata.json` contains:
- Gene name and reference sequence info
- Search parameters used
- Structure details with scores
- Download statistics

```json
{
  "gene_name": "bamA",
  "reference_organism": "Escherichia coli",
  "reference_sequence_length": 810,
  "identity_cutoff": 0.3,
  "experimental_structures_downloaded": 5,
  "computed_structures_downloaded": 0,
  "similar_structures": [...]
}
```

## Integration with Existing Pipeline

This script can be used alongside or instead of the main structure download pipeline:

- **Input**: Uses existing `data/protein_sequences/{gene}/` FASTA files
- **Output**: Creates `data/protein_structures_similarity/{gene}/` directories  
- **Independent**: Doesn't modify existing pipeline data or cache

## Benefits Over Complex Pipeline

- **Simpler**: Single script, no complex configuration
- **Faster**: Direct PDB API calls, no intermediate processing
- **More Accurate**: Sequence similarity ensures functional relevance
- **Flexible**: Easy to adjust parameters for different genes
- **Reliable**: Uses official PDB infrastructure