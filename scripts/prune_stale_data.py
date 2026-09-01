#!/usr/bin/env python3
"""
scripts/prune_stale_data.py

Automated maintenance utility to delete old data files and conserve disk space.
Cleans up:
- User uploads in `uploads/` older than 7 days.
- Job status metadata files in `jobs/` older than 30 days.
- Job results directories (`Results_job_*` or `Results_*`) older than 30 days.
"""
import os
import shutil
import time
import logging
from pathlib import Path
from datetime import datetime, timedelta

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
    ]
)
logger = logging.getLogger("prune_stale_data")

# Determine project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent
UPLOADS_DIR = PROJECT_ROOT / "uploads"
JOBS_DIR = PROJECT_ROOT / "jobs"

# Retention limits
UPLOAD_RETENTION_DAYS = 7
RESULTS_RETENTION_DAYS = 30

def prune_directory_files(directory: Path, retention_days: int, dry_run: bool = False):
    """Prune all files in a directory that exceed the retention window."""
    if not directory.exists() or not directory.is_dir():
        logger.warning(f"Directory {directory} does not exist. Skipping.")
        return

    now = datetime.now()
    retention_cutoff = now - timedelta(days=retention_days)

    logger.info(f"Scanning {directory} for files older than {retention_days} days (cutoff: {retention_cutoff.strftime('%Y-%m-%d %H:%M:%S')})...")

    files_removed = 0
    bytes_freed = 0

    for file_path in directory.iterdir():
        if file_path.is_file():
            # Get file modification time
            mtime = datetime.fromtimestamp(file_path.stat().st_mtime)
            if mtime < retention_cutoff:
                file_size = file_path.stat().st_size
                logger.info(f"{'[DRY RUN] Would delete' if dry_run else 'Deleting'} file: {file_path.name} (Modified: {mtime.strftime('%Y-%m-%d %H:%M')}, Size: {file_size} bytes)")
                if not dry_run:
                    try:
                        file_path.unlink()
                    except Exception as e:
                        logger.error(f"Failed to delete {file_path.name}: {e}")
                        continue
                files_removed += 1
                bytes_freed += file_size

    logger.info(f"Pruned {files_removed} files from {directory.name}. Freed {bytes_freed} bytes.")

def prune_results_directories(project_root: Path, retention_days: int, dry_run: bool = False):
    """Prune results folders starting with 'Results_' older than retention window."""
    now = datetime.now()
    retention_cutoff = now - timedelta(days=retention_days)

    logger.info(f"Scanning {project_root} for results directories older than {retention_days} days...")

    dirs_removed = 0
    
    for item in project_root.iterdir():
        if item.is_dir() and item.name.startswith("Results_"):
            # Exclude base results templates or fixed tools if any
            mtime = datetime.fromtimestamp(item.stat().st_mtime)
            if mtime < retention_cutoff:
                logger.info(f"{'[DRY RUN] Would delete' if dry_run else 'Deleting'} directory: {item.name} (Modified: {mtime.strftime('%Y-%m-%d %H:%M')})")
                if not dry_run:
                    try:
                        shutil.rmtree(item)
                    except Exception as e:
                        logger.error(f"Failed to delete results directory {item.name}: {e}")
                        continue
                dirs_removed += 1

    logger.info(f"Pruned {dirs_removed} results directories.")

def run_prune(dry_run: bool = False):
    logger.info("Starting stale data pruning routine.")
    # Prune user uploads (7 days)
    prune_directory_files(UPLOADS_DIR, UPLOAD_RETENTION_DAYS, dry_run)
    # Prune job status JSONs (30 days)
    prune_directory_files(JOBS_DIR, RESULTS_RETENTION_DAYS, dry_run)
    # Prune Results directories (30 days)
    prune_results_directories(PROJECT_ROOT, RESULTS_RETENTION_DAYS, dry_run)
    logger.info("Pruning routine complete.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Prune stale EpitopePred files and directories.")
    parser.argument_default = False
    parser.add_argument("--dry-run", action="store_true", help="Perform a dry run without deleting files.")
    args = parser.parse_args()
    
    run_prune(dry_run=args.dry_run)
