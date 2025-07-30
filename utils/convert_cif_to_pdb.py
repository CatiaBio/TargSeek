#!/usr/bin/env python3
"""
Convert CIF.gz files to PDB format for better PyMOL compatibility

This script converts compressed CIF (mmCIF) files to PDB format using PyMOL,
which provides more robust conversion than BioPython for problematic files.

Usage:
    python utils/convert_cif_to_pdb.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR] [--specific-file FILE]

Examples:
    # Convert all CIF.gz files in protein structures directory
    python utils/convert_cif_to_pdb.py
    
    # Convert specific file
    python utils/convert_cif_to_pdb.py --specific-file data/protein_structures/pdb_files/7BGL.cif.gz
    
    # Convert with custom directories
    python utils/convert_cif_to_pdb.py --input-dir data/protein_structures --output-dir data/protein_structures_pdb
"""

import gzip
import argparse
import logging
from pathlib import Path
from typing import Optional
import sys
import tempfile
import os

# Try PyMOL first (more robust for CIF conversion)
try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
    CONVERSION_METHOD = "pymol"
except ImportError:
    PYMOL_AVAILABLE = False
    CONVERSION_METHOD = "biopython"

# BioPython fallback
if not PYMOL_AVAILABLE:
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        BIOPYTHON_AVAILABLE = True
    except ImportError:
        print("Error: Neither PyMOL nor BioPython available.")
        print("Install PyMOL: conda install -c conda-forge pymol-open-source")
        print("Or BioPython: conda install -c conda-forge biopython")
        sys.exit(1)
else:
    BIOPYTHON_AVAILABLE = False

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class CIFtoPDBConverter:
    """Convert CIF files to PDB format"""
    
    def __init__(self):
        """Initialize the converter"""
        if PYMOL_AVAILABLE:
            # Initialize PyMOL in command line mode
            pymol.finish_launching(['pymol', '-cq'])
            logging.info(f"Using PyMOL for conversion")
        elif BIOPYTHON_AVAILABLE:
            self.mmcif_parser = MMCIFParser(QUIET=True)
            self.pdb_io = PDBIO()
            logging.info(f"Using BioPython for conversion")
        
    def convert_single_file(self, cif_file: Path, output_file: Path) -> bool:
        """
        Convert a single CIF file to PDB format
        
        Args:
            cif_file: Path to input CIF file (.cif or .cif.gz)
            output_file: Path to output PDB file
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Create output directory if it doesn't exist
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            if PYMOL_AVAILABLE:
                return self._convert_with_pymol(cif_file, output_file)
            elif BIOPYTHON_AVAILABLE:
                return self._convert_with_biopython(cif_file, output_file)
            else:
                logging.error("No conversion method available")
                return False
                
        except Exception as e:
            logging.error(f"✗ Failed to convert {cif_file}: {e}")
            return False
    
    def _convert_with_pymol(self, cif_file: Path, output_file: Path) -> bool:
        """Convert using PyMOL (more robust)"""
        try:
            # Clear PyMOL session
            cmd.delete("all")
            cmd.reinitialize()
            
            # Create temporary uncompressed file if needed
            temp_file = None
            if cif_file.suffix == '.gz':
                temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False)
                with gzip.open(cif_file, 'rt') as gz_f:
                    temp_file.write(gz_f.read())
                temp_file.close()
                load_file = temp_file.name
            else:
                load_file = str(cif_file)
            
            # Load structure
            cmd.load(load_file, "structure")
            
            # Save as PDB
            if output_file.suffix == '.gz':
                # Save to temporary file first, then compress
                temp_pdb = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
                temp_pdb.close()
                cmd.save(temp_pdb.name, "structure")
                
                # Compress the output
                with open(temp_pdb.name, 'r') as f_in:
                    with gzip.open(output_file, 'wt') as f_out:
                        f_out.write(f_in.read())
                
                # Clean up temp file
                os.unlink(temp_pdb.name)
            else:
                cmd.save(str(output_file), "structure")
            
            # Clean up temporary CIF file
            if temp_file:
                os.unlink(temp_file.name)
            
            logging.info(f"✓ Converted {cif_file.name} -> {output_file.name} (PyMOL)")
            return True
            
        except Exception as e:
            logging.error(f"✗ PyMOL conversion failed for {cif_file}: {e}")
            return False
    
    def _convert_with_biopython(self, cif_file: Path, output_file: Path) -> bool:
        """Convert using BioPython (fallback)"""
        try:
            # Handle compressed files
            if cif_file.suffix == '.gz':
                with gzip.open(cif_file, 'rt') as f:
                    structure = self.mmcif_parser.get_structure('structure', f)
            else:
                structure = self.mmcif_parser.get_structure('structure', str(cif_file))
            
            # Set structure and save as PDB
            self.pdb_io.set_structure(structure)
            
            # Save as compressed PDB if output filename ends with .gz
            if output_file.suffix == '.gz':
                with gzip.open(output_file, 'wt') as f:
                    self.pdb_io.save(f)
            else:
                self.pdb_io.save(str(output_file))
            
            logging.info(f"✓ Converted {cif_file.name} -> {output_file.name} (BioPython)")
            return True
            
        except Exception as e:
            logging.error(f"✗ BioPython conversion failed for {cif_file}: {e}")
            return False
    
    def convert_directory(self, input_dir: Path, output_dir: Path, pattern: str = "*.cif.gz") -> dict:
        """
        Convert all CIF files in a directory
        
        Args:
            input_dir: Directory containing CIF files
            output_dir: Directory to save converted PDB files
            pattern: File pattern to match (default: "*.cif.gz")
            
        Returns:
            Dictionary with conversion results
        """
        if not input_dir.exists():
            logging.error(f"Input directory does not exist: {input_dir}")
            return {}
        
        # Find all CIF files
        cif_files = list(input_dir.rglob(pattern))
        
        if not cif_files:
            logging.warning(f"No files matching pattern '{pattern}' found in {input_dir}")
            return {}
        
        logging.info(f"Found {len(cif_files)} CIF files to convert")
        
        results = {}
        successful = 0
        failed = 0
        
        for cif_file in cif_files:
            # Generate output filename
            # Convert: structure_id.cif.gz -> structure_id.pdb.gz
            if cif_file.name.endswith('.cif.gz'):
                output_name = cif_file.name.replace('.cif.gz', '.pdb.gz')
            elif cif_file.name.endswith('.cif'):
                output_name = cif_file.name.replace('.cif', '.pdb')
            else:
                output_name = f"{cif_file.stem}.pdb"
            
            # Maintain directory structure relative to input_dir
            relative_path = cif_file.relative_to(input_dir)
            output_file = output_dir / relative_path.parent / output_name
            
            # Skip if output already exists (unless we want to overwrite)
            if output_file.exists():
                logging.info(f"Skipping {cif_file.name} - output already exists")
                results[str(cif_file)] = "skipped"
                continue
            
            # Convert the file
            success = self.convert_single_file(cif_file, output_file)
            results[str(cif_file)] = "success" if success else "failed"
            
            if success:
                successful += 1
            else:
                failed += 1
        
        logging.info(f"\n=== Conversion Summary ===")
        logging.info(f"Successful: {successful}")
        logging.info(f"Failed: {failed}")
        logging.info(f"Skipped: {len([r for r in results.values() if r == 'skipped'])}")
        logging.info(f"Total processed: {len(results)}")
        
        return results
    
    def convert_in_place(self, input_dir: Path, pattern: str = "*.cif.gz") -> dict:
        """
        Convert CIF files to PDB in the same directory
        
        Args:
            input_dir: Directory containing CIF files
            pattern: File pattern to match
            
        Returns:
            Dictionary with conversion results
        """
        return self.convert_directory(input_dir, input_dir, pattern)

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Convert CIF files to PDB format')
    parser.add_argument('--input-dir', type=str, default='data/protein_structures',
                       help='Input directory containing CIF files (default: data/protein_structures)')
    parser.add_argument('--output-dir', type=str, 
                       help='Output directory for PDB files (default: same as input-dir)')
    parser.add_argument('--specific-file', type=str,
                       help='Convert specific file instead of directory')
    parser.add_argument('--pattern', type=str, default='*.cif.gz',
                       help='File pattern to match (default: *.cif.gz)')
    parser.add_argument('--overwrite', action='store_true',
                       help='Overwrite existing PDB files')
    
    args = parser.parse_args()
    
    if not PYMOL_AVAILABLE and not BIOPYTHON_AVAILABLE:
        logging.error("Neither PyMOL nor BioPython is available")
        sys.exit(1)
    
    # Initialize converter
    converter = CIFtoPDBConverter()
    
    # Convert specific file
    if args.specific_file:
        cif_file = Path(args.specific_file)
        if not cif_file.exists():
            logging.error(f"File does not exist: {cif_file}")
            sys.exit(1)
        
        # Generate output filename
        if args.output_dir:
            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            if cif_file.name.endswith('.cif.gz'):
                output_name = cif_file.name.replace('.cif.gz', '.pdb.gz')
            else:
                output_name = cif_file.name.replace('.cif', '.pdb')
            output_file = output_dir / output_name
        else:
            # Same directory as input
            if cif_file.name.endswith('.cif.gz'):
                output_file = cif_file.parent / cif_file.name.replace('.cif.gz', '.pdb.gz')
            else:
                output_file = cif_file.parent / cif_file.name.replace('.cif', '.pdb')
        
        success = converter.convert_single_file(cif_file, output_file)
        sys.exit(0 if success else 1)
    
    # Convert directory
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    
    results = converter.convert_directory(input_dir, output_dir, args.pattern)
    
    # Exit with appropriate code
    failed_count = len([r for r in results.values() if r == 'failed'])
    sys.exit(0 if failed_count == 0 else 1)

if __name__ == "__main__":
    main()