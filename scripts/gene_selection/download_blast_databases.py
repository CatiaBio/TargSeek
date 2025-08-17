#!/usr/bin/env python3
"""
Download and format BLAST databases for epitope cross-reactivity analysis.

This script selectively downloads protein databases based on user configuration.
It downloads proteomes from various organisms and formats them for BLAST searches.
Includes timestamp-based refresh mechanism to avoid unnecessary re-downloads.

Features:
- Selective database downloading based on configuration
- Timestamp tracking with automatic refresh after 6 months
- Complete database integrity checking
- Progress tracking and detailed logging
- Automatic cleanup of temporary files

Input:
- Configuration specifying which databases to download
- BLAST database URLs and metadata

Output:
- Downloaded and formatted BLAST databases
- Download log with status and statistics
- Timestamp files for refresh tracking
- Completion sentinel file
"""

import os
import sys
import json
import urllib.request
import gzip
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Any
import time
import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuration
DEFAULT_REFRESH_INTERVAL_MONTHS = 6  # Default refresh interval
TIMESTAMP_FORMAT = "%Y-%m-%d %H:%M:%S"

# Database definitions with URLs and metadata
DATABASE_CATALOG = {
    # Standard databases
    "human": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz",
        "description": "Human Proteome",
        "category": "standard",
        "size_mb": 25
    },
    "viral_proteins": {
        "url": "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz",
        "description": "Viral Proteins",
        "category": "standard",
        "size_mb": 15
    },
    "fungal_proteins": {
        "url": "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/fungi.1.protein.faa.gz",
        "description": "Fungal Proteins",
        "category": "standard",
        "size_mb": 80
    },
    
    # Farm animals & livestock
    "cow": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000009136/UP000009136_9913.fasta.gz",
        "description": "Cow Proteome",
        "category": "livestock",
        "size_mb": 22
    },
    "goat": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000291000/UP000291000_9925.fasta.gz",
        "description": "Goat Proteome",
        "category": "livestock",
        "size_mb": 20
    },
    "horse": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002281/UP000002281_9796.fasta.gz",
        "description": "Horse Proteome",
        "category": "livestock",
        "size_mb": 20
    },
    "pig": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000008227/UP000008227_9823.fasta.gz",
        "description": "Pig Proteome",
        "category": "livestock",
        "size_mb": 25
    },
    "chicken": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000539/UP000000539_9031.fasta.gz",
        "description": "Chicken Proteome",
        "category": "livestock",
        "size_mb": 17
    },
    "buffalo": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000261620/UP000261620_94237.fasta.gz",
        "description": "Water Buffalo Proteome",
        "category": "livestock",
        "size_mb": 20
    },
    "sheep": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002356/UP000002356_9940.fasta.gz",
        "description": "Sheep Proteome",
        "category": "livestock",
        "size_mb": 20
    },
    
    # Companion/domestic animals
    "cat": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000011712/UP000011712_9685.fasta.gz",
        "description": "Cat Proteome",
        "category": "companion",
        "size_mb": 20
    },
    # Note: Dog proteome (Canis lupus familiaris) UP000002254 is not available on UniProt FTP
    # Alternative: NCBI RefSeq or other sources can be added if needed
    
    # Research model organisms
    "mouse": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.fasta.gz",
        "description": "Mouse Proteome",
        "category": "model",
        "size_mb": 55
    },
    "rat": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002494/UP000002494_10116.fasta.gz",
        "description": "Rat Proteome",
        "category": "model",
        "size_mb": 30
    },
    "rabbit": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000001811/UP000001811_9986.fasta.gz",
        "description": "Rabbit Proteome",
        "category": "model",
        "size_mb": 22
    },
    
    # Aquatic food organisms
    "salmon": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000087266/UP000087266_8030.fasta.gz",
        "description": "Atlantic Salmon Proteome",
        "category": "aquatic",
        "size_mb": 35
    },
    "zebrafish": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000437/UP000000437_7955.fasta.gz",
        "description": "Zebrafish Proteome",
        "category": "aquatic",
        "size_mb": 45
    },
    
    # Plant crops
    "arabidopsis": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000006548/UP000006548_3702.fasta.gz",
        "description": "Arabidopsis Proteome",
        "category": "plant",
        "size_mb": 35
    },
    "rice": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000059680/UP000059680_39947.fasta.gz",
        "description": "Rice Proteome",
        "category": "plant",
        "size_mb": 40
    },
    "soybean": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000008827/UP000008827_3847.fasta.gz",
        "description": "Soybean Proteome",
        "category": "plant",
        "size_mb": 56
    },
    
    # Bacterial pathogens
    "ecoli": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000625/UP000000625_83333.fasta.gz",
        "description": "E. coli K-12 Proteome",
        "category": "pathogen",
        "size_mb": 1.3
    },
    "salmonella": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000000541/UP000000541_28901.fasta.gz",
        "description": "Salmonella enterica Proteome",
        "category": "pathogen",
        "size_mb": 1.5
    },
    "staph": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/UP000008816/UP000008816_93061.fasta.gz",
        "description": "Staphylococcus aureus Proteome",
        "category": "pathogen",
        "size_mb": 0.8
    }
}

def check_dependencies():
    """Check if required tools are available."""
    try:
        result = subprocess.run(['makeblastdb', '-help'], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            logger.error("makeblastdb not found. Please install BLAST+ tools.")
            return False
    except FileNotFoundError:
        logger.error("makeblastdb not found. Please install BLAST+ tools.")
        return False
    
    return True

def get_database_timestamp(base_dir: str, db_name: str) -> str:
    """Get the timestamp file path for a database."""
    return os.path.join(base_dir, f".{db_name}_last_updated.timestamp")

def read_timestamp(timestamp_file: str) -> datetime.datetime:
    """Read timestamp from file."""
    try:
        if os.path.exists(timestamp_file):
            with open(timestamp_file, 'r') as f:
                timestamp_str = f.read().strip()
                return datetime.datetime.strptime(timestamp_str, TIMESTAMP_FORMAT)
    except Exception as e:
        logger.debug(f"Could not read timestamp from {timestamp_file}: {e}")
    
    # Return epoch if no valid timestamp found
    return datetime.datetime(1970, 1, 1)

def write_timestamp(timestamp_file: str, timestamp: datetime.datetime = None):
    """Write current timestamp to file."""
    if timestamp is None:
        timestamp = datetime.datetime.now()
    
    try:
        os.makedirs(os.path.dirname(timestamp_file), exist_ok=True)
        with open(timestamp_file, 'w') as f:
            f.write(timestamp.strftime(TIMESTAMP_FORMAT))
        logger.debug(f"Updated timestamp: {timestamp_file}")
    except Exception as e:
        logger.warning(f"Could not write timestamp to {timestamp_file}: {e}")

def should_refresh_database(base_dir: str, db_name: str, refresh_months: int = DEFAULT_REFRESH_INTERVAL_MONTHS, force_refresh: bool = False) -> bool:
    """Check if database should be refreshed based on age."""
    if force_refresh:
        logger.info(f"Database {db_name} refresh forced by configuration")
        return True
    
    timestamp_file = get_database_timestamp(base_dir, db_name)
    last_updated = read_timestamp(timestamp_file)
    
    # Calculate refresh threshold
    refresh_threshold = datetime.datetime.now() - datetime.timedelta(days=30 * refresh_months)
    
    should_refresh = last_updated < refresh_threshold
    
    if should_refresh:
        age_days = (datetime.datetime.now() - last_updated).days
        if age_days > 0:
            logger.info(f"Database {db_name} is {age_days} days old, needs refresh (threshold: {refresh_months} months)")
        else:
            logger.info(f"Database {db_name} has no valid timestamp, needs download")
    else:
        age_days = (datetime.datetime.now() - last_updated).days
        logger.info(f"Database {db_name} is {age_days} days old, still fresh (threshold: {refresh_months} months)")
    
    return should_refresh

def database_exists_and_complete(fasta_path: str) -> bool:
    """Check if database exists and is properly formatted."""
    # Check if the main FASTA file exists
    if not os.path.exists(fasta_path):
        return False
    
    # Check if BLAST database files exist (.phr, .pin, .psq are essential)
    required_extensions = ['.phr', '.pin', '.psq']
    for ext in required_extensions:
        if not os.path.exists(fasta_path + ext):
            return False
    
    # Check if files are not empty
    try:
        if os.path.getsize(fasta_path) == 0:
            return False
        for ext in required_extensions:
            if os.path.getsize(fasta_path + ext) == 0:
                return False
    except OSError:
        return False
    
    return True

def download_file(url: str, output_path: str, description: str) -> bool:
    """Download a file from URL with progress tracking."""
    try:
        logger.info(f"🔽 Downloading {description}...")
        logger.info(f"  URL: {url}")
        logger.info(f"  Output: {output_path}")
        
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Download with progress
        with urllib.request.urlopen(url) as response:
            total_size = int(response.headers.get('Content-Length', 0))
            downloaded = 0
            
            with open(output_path, 'wb') as f:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        logger.debug(f"  Progress: {percent:.1f}% ({downloaded}/{total_size} bytes)")
        
        logger.info(f"  ✅ Downloaded {description}")
        return True
        
    except Exception as e:
        logger.error(f"  ❌ Failed to download {description}: {e}")
        return False

def decompress_file(input_path: str, output_path: str) -> bool:
    """Decompress a gzipped file."""
    try:
        logger.info(f"📦 Decompressing {os.path.basename(input_path)}...")
        
        with gzip.open(input_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                f_out.write(f_in.read())
        
        logger.info(f"  ✅ Decompressed to {os.path.basename(output_path)}")
        return True
        
    except Exception as e:
        logger.error(f"  ❌ Failed to decompress {input_path}: {e}")
        return False

def format_blast_database(fasta_path: str, description: str) -> bool:
    """Format a FASTA file as a BLAST database."""
    try:
        logger.info(f"📚 Formatting BLAST database for {description}...")
        
        cmd = [
            'makeblastdb',
            '-in', fasta_path,
            '-dbtype', 'prot',
            '-title', description,
            '-logfile', fasta_path + '.makeblastdb.log'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info(f"  ✅ BLAST database formatted for {description}")
            return True
        else:
            logger.error(f"  ❌ Failed to format BLAST database: {result.stderr}")
            return False
            
    except Exception as e:
        logger.error(f"  ❌ Failed to format BLAST database for {description}: {e}")
        return False

def process_database(db_name: str, db_info: Dict, base_dir: str, refresh_months: int = DEFAULT_REFRESH_INTERVAL_MONTHS, force_refresh: bool = False) -> Dict[str, Any]:
    """Download and process a single database with timestamp-based refresh."""
    start_time = time.time()
    
    # File paths
    gz_path = os.path.join(base_dir, f"{db_name}_proteome.fasta.gz")
    fasta_path = os.path.join(base_dir, f"{db_name}_proteome.fasta")
    timestamp_file = get_database_timestamp(base_dir, db_name)
    
    result = {
        "database": db_name,
        "description": db_info["description"],
        "category": db_info["category"],
        "expected_size_mb": db_info["size_mb"],
        "status": "failed",
        "download_time": 0,
        "file_size_mb": 0,
        "error": None,
        "last_updated": None,
        "age_days": None
    }
    
    try:
        # Check if database exists and is complete
        db_exists = database_exists_and_complete(fasta_path)
        
        # Get age information
        last_updated = read_timestamp(timestamp_file)
        age_days = None
        if last_updated.year > 1970:  # Valid timestamp
            age_days = (datetime.datetime.now() - last_updated).days
            result["last_updated"] = last_updated.strftime(TIMESTAMP_FORMAT)
            result["age_days"] = age_days
        
        # Check if we should refresh
        needs_refresh = should_refresh_database(base_dir, db_name, refresh_months, force_refresh)
        
        # Skip if already processed and fresh
        if db_exists and not needs_refresh:
            age_str = f"age: {age_days} days" if age_days is not None else "age: unknown"
            logger.info(f"✅ {db_info['description']} is up-to-date ({age_str})")
            result["status"] = "already_exists"
            result["file_size_mb"] = os.path.getsize(fasta_path) / (1024 * 1024)
            return result
        
        # If database exists but needs refresh, log the reason
        if db_exists and needs_refresh:
            age_str = f"age: {age_days} days" if age_days is not None else "age: unknown"
            logger.info(f"🔄 {db_info['description']} needs refresh ({age_str})")
        elif not db_exists:
            logger.info(f"📥 {db_info['description']} not found, downloading...")
        
        # Clean up old files if they exist
        for old_file in [fasta_path, gz_path]:
            if os.path.exists(old_file):
                os.remove(old_file)
                logger.debug(f"Removed old file: {old_file}")
        
        # Clean up old BLAST database files
        blast_extensions = ['.phr', '.pin', '.psq', '.pdb', '.pot', '.ptf', '.pto', '.makeblastdb.log']
        for ext in blast_extensions:
            old_blast_file = fasta_path + ext
            if os.path.exists(old_blast_file):
                os.remove(old_blast_file)
                logger.debug(f"Removed old BLAST file: {old_blast_file}")
        
        # Download compressed file
        if not download_file(db_info["url"], gz_path, db_info["description"]):
            result["error"] = "Download failed"
            return result
        
        # Decompress
        if not decompress_file(gz_path, fasta_path):
            result["error"] = "Decompression failed"
            return result
        
        # Format for BLAST
        if not format_blast_database(fasta_path, db_info["description"]):
            result["error"] = "BLAST formatting failed"
            return result
        
        # Update timestamp after successful processing
        current_time = datetime.datetime.now()
        write_timestamp(timestamp_file, current_time)
        
        # Success
        result["status"] = "success"
        result["download_time"] = time.time() - start_time
        result["file_size_mb"] = os.path.getsize(fasta_path) / (1024 * 1024)
        result["last_updated"] = current_time.strftime(TIMESTAMP_FORMAT)
        result["age_days"] = 0
        
        # Clean up compressed file to save space
        if os.path.exists(gz_path):
            os.remove(gz_path)
            logger.info(f"  🗑️ Removed compressed file {os.path.basename(gz_path)}")
        
        logger.info(f"  ✅ {db_info['description']} successfully processed and timestamped")
        
    except Exception as e:
        result["error"] = str(e)
        logger.error(f"❌ Failed to process {db_name}: {e}")
    
    return result

def main():
    """Main function to download selected BLAST databases."""
    
    # Get Snakemake parameters
    blast_databases_config = snakemake.config["blast_databases"]
    base_dir = snakemake.params.base_dir
    
    # Get refresh configuration
    refresh_config = snakemake.config.get("blast_database_refresh", {})
    refresh_months = refresh_config.get("refresh_interval_months", DEFAULT_REFRESH_INTERVAL_MONTHS)
    force_refresh = refresh_config.get("force_refresh_all", False)
    
    # Output files
    download_log = snakemake.output.download_log
    sentinel_file = snakemake.output.download_sentinel
    
    logger.info("Starting BLAST database downloads...")
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Create base directory
    os.makedirs(base_dir, exist_ok=True)
    
    # Filter enabled databases
    enabled_databases = {
        name: DATABASE_CATALOG[name] 
        for name, enabled in blast_databases_config.items() 
        if enabled and name in DATABASE_CATALOG
    }
    
    if not enabled_databases:
        logger.warning("No databases enabled in configuration")
        # Create empty log and sentinel
        with open(download_log, 'w') as f:
            json.dump({"databases": [], "summary": {"total": 0, "success": 0, "failed": 0}}, f, indent=2)
        Path(sentinel_file).touch()
        return
    
    logger.info(f"Processing {len(enabled_databases)} enabled databases:")
    if force_refresh:
        logger.info(f"Refresh policy: FORCE REFRESH ALL enabled (ignoring timestamps)")
    else:
        logger.info(f"Refresh policy: Re-download databases older than {refresh_months} months")
    
    # Calculate total expected download size
    total_size_mb = sum(db["size_mb"] for db in enabled_databases.values())
    logger.info(f"Expected total download size: ~{total_size_mb:.1f} MB (if all need downloading)")
    
    # Group by category for organized display
    by_category = {}
    for name, db_info in enabled_databases.items():
        category = db_info["category"]
        if category not in by_category:
            by_category[category] = []
        by_category[category].append(name)
    
    for category, db_names in by_category.items():
        logger.info(f"  📁 {category.title()}: {', '.join(db_names)}")
    
    logger.info("")
    
    # Process each database
    results = []
    for db_name, db_info in enabled_databases.items():
        logger.info(f"=== Processing {db_info['description']} ===")
        result = process_database(db_name, db_info, base_dir, refresh_months, force_refresh)
        results.append(result)
        logger.info("")
    
    # Generate summary
    successful = [r for r in results if r["status"] == "success"]
    already_exists = [r for r in results if r["status"] == "already_exists"]
    failed = [r for r in results if r["status"] == "failed"]
    
    # Calculate refresh statistics
    refresh_threshold_days = 30 * refresh_months
    refreshed_count = len([r for r in successful if r["age_days"] is not None])
    up_to_date_count = len([r for r in already_exists if r["age_days"] is not None and r["age_days"] < refresh_threshold_days])
    
    summary = {
        "total": len(results),
        "success": len(successful),
        "already_exists": len(already_exists),
        "failed": len(failed),
        "refreshed": refreshed_count,
        "up_to_date": up_to_date_count,
        "refresh_threshold_days": refresh_threshold_days,
        "total_download_time": sum(r["download_time"] for r in results),
        "total_size_mb": sum(r["file_size_mb"] for r in results if r["file_size_mb"] > 0)
    }
    
    # Save detailed log
    log_data = {
        "databases": results,
        "summary": summary,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "base_directory": base_dir
    }
    
    with open(download_log, 'w') as f:
        json.dump(log_data, f, indent=2)
    
    # Final status
    logger.info("="*60)
    logger.info("📊 BLAST DATABASE DOWNLOAD SUMMARY")
    logger.info("="*60)
    logger.info(f"Total databases processed: {summary['total']}")
    logger.info(f"Successfully downloaded/refreshed: {summary['success']}")
    logger.info(f"Already up-to-date: {summary['already_exists']}")
    logger.info(f"Failed downloads: {summary['failed']}")
    if summary['refreshed'] > 0:
        logger.info(f"Databases refreshed (>{refresh_months} months old): {summary['refreshed']}")
    if summary['up_to_date'] > 0:
        logger.info(f"Databases still fresh (<{refresh_months} months): {summary['up_to_date']}")
    logger.info(f"Refresh threshold: {summary['refresh_threshold_days']} days ({refresh_months} months)")
    logger.info(f"Total download time: {summary['total_download_time']:.1f} seconds")
    logger.info(f"Total data size: {summary['total_size_mb']:.1f} MB")
    logger.info(f"Database location: {base_dir}")
    
    if failed:
        logger.warning("Failed databases:")
        for result in failed:
            logger.warning(f"  - {result['database']}: {result['error']}")
    
    if summary['failed'] > 0:
        logger.error("Some database downloads failed!")
        sys.exit(1)
    
    # Create sentinel file
    Path(sentinel_file).touch()
    logger.info(f"✅ BLAST database setup complete! Sentinel: {sentinel_file}")

if __name__ == "__main__":
    main()