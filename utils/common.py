import os
import re
from datetime import datetime
import logging
from Bio import SeqIO

# Configure logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def create_timestamped_dir(base_dir):
    """
    Creates a directory under base_dir with a timestamp in the format YYYY-MM-DD_HH-MM-SS.
    
    Args:
        base_dir (str): Base directory path (e.g., "Results").
    
    Returns:
        str: Path to the created timestamped directory.
    
    Raises:
        OSError: If directory creation fails.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    timestamp_dir = os.path.join(base_dir, timestamp)
    try:
        os.makedirs(timestamp_dir, exist_ok=True)
        logger.info(f"Created timestamped directory: {timestamp_dir}")
        return timestamp_dir
    except OSError as e:
        logger.error(f"Failed to create directory {timestamp_dir}: {e}")
        raise

def get_batch_number_from_filename(filename):
    """
    Extracts a batch number from the filename using a regex pattern.
    
    Args:
        filename (str): Path or name of the FASTA file.
    
    Returns:
        str: Extracted batch number, or "unknown" if no number is found.
    
    Raises:
        ValueError: If the filename is empty or invalid.
    """
    if not filename:
        logger.error("Filename is empty or None")
        raise ValueError("Filename cannot be empty")
    
    # Extract filename from path
    base_name = os.path.basename(filename)
    # Match digits in patterns like batch1, proteins_2, etc.
    match = re.search(r'\d+', base_name)
    if match:
        batch_no = match.group()
        logger.debug(f"Extracted batch number {batch_no} from {filename}")
        return batch_no
    else:
        logger.warning(f"No batch number found in {filename}, using 'unknown'")
        return "unknown"

def read_fasta(fasta_file):
    """
    Reads a FASTA file and returns a list of sequence records.
    
    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        list: List of Bio.SeqRecord.SeqRecord objects.
    
    Raises:
        FileNotFoundError: If the FASTA file does not exist.
        ValueError: If the FASTA file is invalid.
    """
    if not os.path.exists(fasta_file):
        logger.error(f"FASTA file not found: {fasta_file}")
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            logger.warning(f"FASTA file {fasta_file} is empty")
        else:
            logger.debug(f"Read {len(records)} sequences from {fasta_file}")
        return records
    except Exception as e:
        logger.error(f"Failed to parse FASTA file {fasta_file}: {e}")
        raise ValueError(f"Invalid FASTA file: {fasta_file}")
