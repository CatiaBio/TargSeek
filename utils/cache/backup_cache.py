#!/usr/bin/env python3
"""
Cache Backup and Restore Utility
=================================

This script helps backup and restore cache files to prevent data loss.
"""

import json
import shutil
import logging
from pathlib import Path
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def backup_cache():
    """Backup all cache files with timestamp"""
    cache_dir = Path("cache")
    backup_dir = Path(f"cache_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    
    if not cache_dir.exists():
        logging.info("No cache directory found")
        return
    
    backup_dir.mkdir(exist_ok=True)
    
    # Copy all cache files
    for cache_file in cache_dir.rglob("*.json"):
        relative_path = cache_file.relative_to(cache_dir)
        backup_file = backup_dir / relative_path
        backup_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cache_file, backup_file)
        logging.info(f"Backed up: {cache_file} -> {backup_file}")
    
    logging.info(f"✅ Cache backup completed: {backup_dir}")

def restore_cache(backup_dir):
    """Restore cache from backup directory"""
    backup_path = Path(backup_dir)
    cache_dir = Path("cache")
    
    if not backup_path.exists():
        logging.error(f"Backup directory not found: {backup_dir}")
        return
    
    cache_dir.mkdir(exist_ok=True)
    
    # Copy files back
    for backup_file in backup_path.rglob("*.json"):
        relative_path = backup_file.relative_to(backup_path)
        cache_file = cache_dir / relative_path
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(backup_file, cache_file)
        logging.info(f"Restored: {backup_file} -> {cache_file}")
    
    logging.info(f"✅ Cache restored from: {backup_dir}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        if sys.argv[1] == "backup":
            backup_cache()
        elif sys.argv[1] == "restore" and len(sys.argv) > 2:
            restore_cache(sys.argv[2])
        else:
            print("Usage: python backup_cache.py [backup|restore <backup_dir>]")
    else:
        backup_cache()