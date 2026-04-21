
#!/usr/bin/env python3
## Importing the Libraries
import sys
import os
import requests
import logging
from pathlib import Path
from Bio import SeqIO
from io import StringIO
from datetime import datetime
import tempfile
import pandas as pd
from contextlib import contextmanager
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import shutil
import multiprocessing as mp
from functools import partial
import itertools
import atexit


## Importing tools
from modules.mhc_i import run_mhc1
from modules.mhc_ii import run_mhc2
from modules.netctl import run_netctl
from modules.netchop import run_netchop
from modules.bcell import run_bcell
from modules.pepmatch import run_pepmatch
from modules.psortb import run_psortb
from modules.iapred import run_iapred
from modules.algpred import run_algpred
from modules.instability import run_instability
from modules.molweight import run_molwt
from modules.toxicity import run_toxinpred
from modules.gutflora import run_gutflora
from modules.virulence import run_virulence
from modules.deeptmhmm import run_deeptmhmm
from modules.netsol import run_netsol
from modules.signalp import run_signalp
from modules.human_1 import run_human
from modules.clbtope import run_clbtope
from modules.toxin_epitope import run_toxinpred3
from modules.wolfpsort import run_wolf_psort
from modules.vaccine_assembly import run_assembly
from modules.vaccine_assembly import run_assembly, order_dict as assembly_orders
from modules.assembly_graph import plot_vaccine_architecture
from modules.esmfold import run_esmfold  # <--- ADD THIS
from modules.sasa_filter import run_sasa_analysis

# Define root directory
ROOT_DIR = os.environ.get("VAXELAN_ROOT", Path(__file__).resolve().parent)
ROOT_DIR = Path(ROOT_DIR)

# Define output directory with timestamp
TIMESTAMP = datetime.now().strftime('%Y%m%d_%H%M%S')
RESULTS_DIR = ROOT_DIR / f"Results_{TIMESTAMP}"

# Batch size for processing sequences
BATCH_SIZE = 100

def extract_uniprot_id(header):
    """Extract UniProt ID from FASTA header."""
    logger = logging.getLogger('VaxElan')
    try:
        if '|' in header:
            parts = header.split('|')
            for part in parts:
                part = part.strip()
                if len(part) >= 5 and all(c.isalnum() or c in {'-', '_'} for c in part):
                    return part
        # Fallback: Take the first word after '>'
        id_part = header.split()[0].replace(">", "").strip()
        if id_part:
            return id_part
        logger.warning(f"Could not extract UniProt ID from header: {header}")
        return "unknown_id"
    except Exception as e:
        logger.error(f"Error extracting UniProt ID from header '{header}': {e}")
        return "unknown_id"

def generate_tool_output_path(timestamp, strategy_no, batch_no, uniprot_id, tool_name, results_dir_base):
    """Generate the output file path for a given tool."""
    tool_extensions = {
        'netchop': 'txt',
        'instability': 'txt',
        'molwt': 'txt',
        'gutflora': 'txt',
        'virulence': 'txt',
        'deeptmhmm': 'csv',
        'signalp': 'csv',
    }
    extension = tool_extensions.get(tool_name, 'csv')
    tool_output_names = {
        'mhc1': 'mhci_out',
        'mhc2': 'mhcii_out',
        'wolfpsort': 'wolfpsort_out',
        'toxinpred3': 'toxinpred3_out',
    }
    output_filename = tool_output_names.get(tool_name, f"{tool_name}_out")
    
    if timestamp:
        # Legacy behavior (CLI)
        output_path = (
            results_dir_base /
            f"Results_{timestamp}" /
            f"strategy_{strategy_no}" /
            f"batch{batch_no}" /
            f"{uniprot_id}_output" /
            f"{output_filename}.{extension}"
        )
    else:
        # Direct behavior (Web or refactored CLI)
        output_path = (
            results_dir_base /
            f"batch{batch_no}" /
            f"{uniprot_id}_output" /
            f"{output_filename}.{extension}"
        )
    return output_path

def print_status(msg, status="info"):
    colors = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")
    logger = logging.getLogger('VaxElan')
    if status == "info":
        logger.info(msg)
    elif status == "warning":
        logger.warning(msg)
    elif status == "error":
        logger.error(msg)
    else:
        logger.debug(msg)

def setup_logging(output_dir):
    """Sets up logging to both file and console with proper flushing and isolation."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / f"vax_elan_{TIMESTAMP}.log"

    # Create logger
    logger = logging.getLogger('VaxElan')
    logger.setLevel(logging.DEBUG)
    logger.propagate = False  # Prevent double logging if root logger is used elsewhere

    # Clear old handlers (important in case of re-runs in Jupyter or similar)
    if logger.hasHandlers():
        logger.handlers.clear()

    # Log format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Flush logs on exit
    atexit.register(lambda: [h.flush() for h in logger.handlers if hasattr(h, 'flush')])

    # Confirmation message
    logger.info("Logging initialized for VaxElan")
    print("\033[94mLogging initialized for VaxElan\033[0m")

    return logger

def fetch_uniprot_sequence(uniprot_id, logger):
    logger.info(f"Fetching UniProt sequence for ID: {uniprot_id}")
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = session.get(url, timeout=10)
        response.raise_for_status()
        logger.debug(f"Successfully fetched UniProt sequence for {uniprot_id}")
        return response.text
    except requests.RequestException as e:
        logger.error(f"Failed to fetch UniProt ID: {uniprot_id}: {e}")
        print_status(f"Error: Failed to fetch UniProt ID: {uniprot_id}", "error")
        raise RuntimeError(f"Failed to fetch UniProt ID: {uniprot_id}: {e}")

def detect_pathogen_type_from_uniprot(uniprot_id):
    logger = logging.getLogger('VaxElan')
    logger.info(f"Detecting pathogen type from UniProt ID: {uniprot_id}")
    try:
        # Query UniProt for taxonomy information
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        session = requests.Session()
        retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
        session.mount("https://", HTTPAdapter(max_retries=retries))
        response = session.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()

        # Extract taxonomic lineage
        if 'organism' in data and 'lineage' in data['organism']:
            lineage = [taxon.lower() for taxon in data['organism']['lineage']]
            logger.debug(f"Taxonomic lineage for {uniprot_id}: {lineage}")

            # Classify based on lineage
            if "bacteria" in lineage:
                logger.info(f"Pathogen type detected: bacteria")
                return "bacteria"
            elif "viruses" in lineage:
                logger.info(f"Pathogen type detected: virus")
                return "virus"
            elif any(taxon in lineage for taxon in ["fungi", "ascomycota", "basidiomycota", "mucoromycota"]):
                logger.info(f"Pathogen type detected: fungi")
                return "fungi"
            elif any(taxon in lineage for taxon in ["apicomplexa", "kinetoplastida", "amoebozoa", "trichomonas"]):
                logger.info(f"Pathogen type detected: protozoa")
                return "protozoa"
            else:
                logger.warning(f"No recognizable pathogen type in lineage for UniProt ID: {uniprot_id}")
                return "unknown"
        else:
            logger.warning(f"No organism or lineage information found for UniProt ID: {uniprot_id}")
            return "unknown"
    except requests.RequestException as e:
        logger.error(f"Error fetching taxonomy data for {uniprot_id}: {e}")
        print_status(f"Error: Failed to fetch taxonomy data for {uniprot_id}: {e}", "error")
        return "unknown"

def detect_pathogen_type_from_fasta(input_file):
    logger = logging.getLogger('VaxElan')
    logger.info(f"Detecting pathogen type from FASTA file: {input_file}")
    try:
        input_file = Path(input_file)
        if not input_file.is_file():
            logger.error(f"FASTA file does not exist: {input_file}")
            print_status(f"Error: FASTA file does not exist: {input_file}", "error")
            return "unknown"

        with open(input_file, 'r') as f:
            raw_content = f.read().strip()
            if not raw_content:
                logger.error(f"FASTA file is empty: {input_file}")
                print_status(f"Error: FASTA file is empty: {input_file}", "error")
                return "unknown"
            if '>' not in raw_content:
                logger.error(f"No FASTA headers found in file: {input_file}")
                print_status(f"Error: No FASTA headers found in file: {input_file}", "error")
                return "unknown"
            cleaned_content = raw_content[raw_content.index('>'):]

        # Parse FASTA headers
        headers = []
        for record in SeqIO.parse(StringIO(cleaned_content), "fasta"):
            headers.append(record.description)
            logger.debug(f"Processing FASTA header: {record.description}")

        # Try to extract UniProt ID or organism name from the first header
        first_header = headers[0]
        uniprot_id = extract_uniprot_id(first_header)
        organism_name = None

        # If a UniProt ID is found, use it to fetch taxonomy
        if uniprot_id != "unknown_id":
            return detect_pathogen_type_from_uniprot(uniprot_id)
        else:
            # Try to extract organism name from header (common in FASTA format, e.g., after 'OS=')
            header_parts = first_header.split()
            for part in header_parts:
                if part.startswith("OS="):
                    organism_name = part[3:].strip()
                    break
            if not organism_name:
                # Fallback: take the second field as organism name if it looks reasonable
                if len(header_parts) > 1:
                    potential_name = header_parts[1]
                    if len(potential_name) > 3 and not potential_name.startswith(("sp", "tr", "|")):
                        organism_name = potential_name
                if not organism_name:
                    logger.warning(f"Could not extract organism name from header: {first_header}")
                    print_status(f"Warning: Could not extract organism name from header: {first_header}", "warning")
                    return "unknown"

        # Query UniProt Taxonomy API with organism name
        logger.info(f"Searching taxonomy for organism: {organism_name}")
        url = f"https://rest.uniprot.org/taxonomy/search?query={organism_name}&fields=scientific_name,lineage"
        try:
            session = requests.Session()
            retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
            session.mount("https://", HTTPAdapter(max_retries=retries))
            response = session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()

            if data.get("results"):
                lineage = [taxon.lower() for taxon in data["results"][0].get("lineage", [])]
                logger.debug(f"Taxonomic lineage for {organism_name}: {lineage}")

                # Classify based on lineage
                if "bacteria" in lineage:
                    logger.info(f"Pathogen type detected: bacteria")
                    return "bacteria"
                elif "viruses" in lineage:
                    logger.info(f"Pathogen type detected: virus")
                    return "virus"
                elif any(taxon in lineage for taxon in ["fungi", "ascomycota", "basidiomycota", "mucoromycota"]):
                    logger.info(f"Pathogen type detected: fungi")
                    return "fungi"
                elif any(taxon in lineage for taxon in ["apicomplexa", "kinetoplastida", "amoebozoa", "trichomonas"]):
                    logger.info(f"Pathogen type detected: protozoa")
                    return "protozoa"
                else:
                    logger.warning(f"No recognizable pathogen type in lineage for organism: {organism_name}")
                    return "unknown"
            else:
                logger.warning(f"No taxonomy results found for organism: {organism_name}")
                return "unknown"
        except requests.RequestException as e:
            logger.error(f"Error fetching taxonomy data for {organism_name}: {e}")
            print_status(f"Error: Failed to fetch taxonomy data for {organism_name}: {e}", "error")
            return "unknown"

    except Exception as e:
        logger.error(f"Error processing FASTA file {input_file}: {e}")
        print_status(f"Error: Error processing FASTA file: {e}", "error")
        return "unknown"

def validate_fasta_file(input_file):
    logger = logging.getLogger('VaxElan')
    logger.info(f"Validating FASTA file: {input_file}")
    try:
        input_file = Path(input_file)
        with open(input_file, 'r') as f:
            raw_content = f.read().strip()
            if not raw_content:
                logger.error(f"FASTA file is empty: {input_file}")
                print_status(f"Error: FASTA file is empty: {input_file}", "error")
                return False
            if '>' not in raw_content:
                logger.error(f"FASTA file lacks a header (no '>' found): {input_file}")
                print_status(f"Error: FASTA file lacks a header: {input_file}", "error")
                return False
            cleaned_content = raw_content[raw_content.index('>'):]
        records = list(SeqIO.parse(StringIO(cleaned_content), "fasta"))
        if not records:
            logger.error(f"No valid sequences found in FASTA file: {input_file}")
            print_status(f"Error: No valid sequences found in FASTA file: {input_file}", "error")
            return False
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        for record in records:
            seq = str(record.seq).upper()
            if not seq:
                logger.error(f"Empty sequence in FASTA file: {input_file}")
                print_status(f"Error: Empty sequence in FASTA file: {input_file}", "error")
                return False
            if not all(c in valid_aa for c in seq):
                logger.error(f"Invalid amino acids in sequence {record.id}: {seq}")
                print_status(f"Error: Invalid amino acids in sequence {record.id}", "error")
                return False
        logger.info(f"FASTA file validated successfully: {len(records)} sequence(s) found")
        print_status(f"FASTA file validated successfully: {len(records)} sequence(s)")
        return True
    except Exception as e:
        logger.error(f"Error validating FASTA file {input_file}: {e}")
        print_status(f"Error: Cannot validate FASTA file: {e}", "error")
        return False

def split_fasta_into_batches(input_file, batch_size, temp_dir):
    """Split a FASTA file into batches, generating batch FASTAs and per-sequence FASTAs."""
    logger = logging.getLogger('VaxElan')
    logger.info(f"Splitting FASTA file {input_file} into batches of {batch_size} sequences")
    input_file = Path(input_file)
    temp_dir = Path(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(input_file, "fasta"))
    total_sequences = len(records)
    if total_sequences == 0:
        logger.error(f"No sequences found in {input_file}")
        print_status(f"Error: No sequences found in {input_file}", "error")
        return [], [], {}
    
    batch_files = []
    batch_dirs = []
    uniprot_mapping = {}
    uniprot_id_counts = {}  # Track duplicates
    
    for batch_idx in range(0, total_sequences, batch_size):
        batch_records = records[batch_idx:batch_idx + batch_size]
        batch_no = (batch_idx // batch_size) + 1
        batch_dir = temp_dir / f"batch_{batch_no}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        batch_file = batch_dir / f"batch_{batch_no}.fasta"
        
        # Write batch FASTA
        try:
            SeqIO.write(batch_records, batch_file, "fasta")
            batch_files.append(batch_file)
            batch_dirs.append(batch_dir)
            logger.info(f"Created batch {batch_no} with {len(batch_records)} sequences: {batch_file}")
            print_status(f"Created batch {batch_no} with {len(batch_records)} sequences")
        except Exception as e:
            logger.error(f"Error writing batch FASTA {batch_file}: {e}")
            continue
        
        # Write per-sequence FASTAs
        for record in batch_records:
            uniprot_id = extract_uniprot_id(record.description)
            uniprot_id_counts[uniprot_id] = uniprot_id_counts.get(uniprot_id, 0) + 1
            if uniprot_id_counts[uniprot_id] > 1:
                unique_id = f"{uniprot_id}_{uniprot_id_counts[uniprot_id]}"
            else:
                unique_id = uniprot_id
            per_seq_path = batch_dir / f"{unique_id}.fasta"
            try:
                SeqIO.write([record], per_seq_path, "fasta")
                uniprot_mapping[unique_id] = {
                    'batch_file': batch_file,
                    'per_seq_file': per_seq_path,
                    'sequence': str(record.seq),
                    'header': record.description,
                    'batch_no': batch_no
                }
                logger.debug(f"Saved {unique_id} to {per_seq_path}")
            except Exception as e:
                logger.error(f"Error writing per-sequence FASTA {per_seq_path}: {e}")
                continue
    
    return batch_files, batch_dirs, uniprot_mapping

def check_dependencies():
    logger = logging.getLogger('VaxElan')
    logger.info("Checking dependencies...")
    try:
        import pandas
        import Bio.SeqIO
        import requests
        logger.info("Python dependencies (pandas, biopython, requests) are installed.")
    except ImportError as e:
        logger.error(f"Missing Python dependency: {e}")
        print_status(f"Error: Missing Python dependency: {e}. Install with 'conda install pandas biopython requests'", "error")
        raise RuntimeError(f"Failed to fetch UniProt ID: {uniprot_id}: {e}")
    if not shutil.which("python3"):
        logger.error("Python3 not found in system PATH")
        print_status("Error: Python3 not found in system PATH", "error")
        raise RuntimeError("Python3 not found in system PATH")
    blastp_path = shutil.which("blastp")
    if not blastp_path:
        blastp_path = os.environ.get("BLAST_PATH", ROOT_DIR / "tools" / "ncbi-blast" / "bin" / "blastp")
        if not Path(blastp_path).is_file():
            error_msg = (
                f"Error: BLASTP not found in system PATH or at {blastp_path}. "
                "Install with 'conda install -c bioconda blast' or download from "
                "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ and set BLAST_PATH environment variable."
            )
            logger.error(error_msg)
            print_status(error_msg, "error")
            raise RuntimeError(error_msg)
    logger.info(f"BLASTP found at {blastp_path}")
    print_status(f"BLASTP found at {blastp_path}", "success")
    clbtope_db = os.environ.get("CLBTOPE_DB", ROOT_DIR / "tools" / "clbtope"/ "clbtope" / "Database")
    if not Path(clbtope_db).is_dir():
        error_msg = (
            f"Error: ClbTope database not found at {clbtope_db}. "
            "Ensure the database is installed and set CLBTOPE_DB environment variable if located elsewhere."
        )
        logger.error(error_msg)
        print_status(error_msg, "error")
        raise RuntimeError(error_msg)
    logger.info(f"ClbTope database found at {clbtope_db}")
    print_status(f"ClbTope database found at {clbtope_db}", "success")
    logger.info("All dependencies verified.")
    print_status("All dependencies verified.", "success")

def select_pathogen_type(_choice=None):
    logger = logging.getLogger('VaxElan')
    types = {"1": "bacteria", "2": "virus", "3": "protozoa", "4": "fungi"}
    # Programmatic bypass
    if _choice is not None:
        # Accept either the key ("1") or the value ("bacteria")
        if _choice in types.values():
            logger.info(f"Programmatic pathogen type: {_choice}")
            return _choice
        selected_type = types.get(str(_choice))
        if selected_type:
            logger.info(f"Programmatic pathogen type: {selected_type}")
            return selected_type
        raise ValueError(f"Invalid pathogen type choice: {_choice}")
    # Interactive CLI
    print("\nSelect pathogen type:")
    print("1. Bacteria")
    print("2. Virus")
    print("3. Protozoa")
    print("4. Fungi")
    while True:
        choice = input("Enter number (1-4): ").strip()
        selected_type = types.get(choice)
        if selected_type:
            logger.info(f"User selected pathogen type: {selected_type}")
            return selected_type
        print_status("Invalid choice. Please enter 1, 2, 3, or 4.", "error")

def validate_pathogen_type(input_value, selected_type, is_uniprot=False):
    logger = logging.getLogger('VaxElan')
    logger.info(f"Validating pathogen type for input: {input_value}, selected: {selected_type}, is_uniprot: {is_uniprot}")
    detected_type = detect_pathogen_type_from_uniprot(input_value) if is_uniprot else detect_pathogen_type_from_fasta(input_value)
    logger.info(f"Detected pathogen type: {detected_type}")
    print_status(f"Detected pathogen type: {detected_type}")
    if detected_type == "unknown":
        logger.warning(f"Could not detect pathogen type from input. Using user-selected type: {selected_type}")
        print_status(f"Could not detect pathogen type from input. Using user-selected type: {selected_type}", "warning")
        return selected_type
    if detected_type != selected_type:
        types_to_choice = {"bacteria": "1", "virus": "2", "protozoa": "3", "fungi": "4"}
        correct_choice = types_to_choice.get(detected_type, "unknown")
        logger.error(f"Selected pathogen type ({selected_type}) does not match detected type ({detected_type})")
        raise ValueError(f"Selected pathogen type ({selected_type}) does not match detected type ({detected_type}).")
    logger.info(f"Pathogen type validated successfully: {selected_type}")
    print_status(f"Pathogen type validated: {selected_type}", "success")
    return selected_type

def select_mhci_method(_choice=None):
    MHCI_METHODS = {
        "a": ("ann", "Artificial Neural Network (ANN)"),
        "b": ("comblib_sidney2008", "Comblib Sidney 2008"),
        "c": ("consensus", "Consensus"),
        "d": ("netmhccons", "NetMHCcons"),
        "e": ("netmhcpan_ba", "NetMHCpan BA"),
        "f": ("netmhcpan_el", "NetMHCpan EL"),
        "g": ("netmhcstabpan", "NetMHCstabpan"),
        "h": ("pickpocket", "PickPocket"),
        "i": ("smm", "Stabilized Matrix Method (SMM)"),
        "j": ("smmpmbec", "SMM with Pseudo Matrix (SMMPMBEC)")
    }
    # Programmatic bypass
    if _choice is not None:
        choice = str(_choice).strip().lower()
        if choice in MHCI_METHODS:
            return choice, MHCI_METHODS[choice][1]
        raise ValueError(f"Invalid MHC-I method choice: {_choice}. Valid: {list(MHCI_METHODS.keys())}")
    # Interactive CLI
    print("\nSelect MHC-I prediction method:")
    for key, (_, name) in MHCI_METHODS.items():
        print(f"{key}. {name}")
    while True:
        choice = input("Enter your choice (a-j): ").strip().lower()
        if choice in MHCI_METHODS:
            print_status(f"{MHCI_METHODS[choice][1]} selected (mapped to \"{choice}\")")
            return choice, MHCI_METHODS[choice][1]
        print_status("Invalid choice. Please enter a-j.", "error")

def select_mhcii_method(_choice=None):
    MHCII_METHODS = {
         "1": ("nmel", "NetMHCIIpan_EL"),
         "2": ("nmba", "NetMHCIIpan_BA"),
         "3": ("cons", "Consensus3"),
         "4": ("nn", "NN_align"),
         "5": ("smm", "SMM_align"),
         "6": ("comblib", "Combinatorial library"),
         "7": ("sturniolo", "Sturniolo"),
         "8": ("nmel42", "NetMHCIIpan_EL-4.2"),
         "9": ("nmba42", "NetMHCIIpan_BA-4.2"),
        "10": ("nmel43", "NetMHCIIpan_EL-4.3"),
        "11": ("nmba43", "NetMHCIIpan_BA-4.3")
    }
    # Programmatic bypass — accept key ("1") or method code ("nmel")
    if _choice is not None:
        choice = str(_choice).strip()
        if choice in MHCII_METHODS:
            return MHCII_METHODS[choice][0], MHCII_METHODS[choice][1]
        # Try matching by method code
        for k, (code, name) in MHCII_METHODS.items():
            if code == choice:
                return code, name
        raise ValueError(f"Invalid MHC-II method choice: {_choice}. Valid keys: {list(MHCII_METHODS.keys())}")
    # Interactive CLI
    print("\nSelect MHC-II prediction method:")
    for key, (_, name) in MHCII_METHODS.items():
        print(f"{key}. {name}")
    while True:
        choice = input("Enter your choice (1-11): ").strip()
        if choice in MHCII_METHODS:
            print_status(f"{MHCII_METHODS[choice][1]} selected (mapped to \"{choice}\")")
            return MHCII_METHODS[choice][0], MHCII_METHODS[choice][1]
        print_status("Invalid selection. Please enter 1-11.", "error")

def run_tool_wrapper(tool_func, tool_name, *args, **kwargs):
    """Wrapper for running tools in multiprocessing, returning (tool_name, status, error)."""
    logger = logging.getLogger('VaxElan')
    try:
        status = tool_func(*args, **kwargs)
        return tool_name, status, None
    except Exception as e:
        logger.error(f"Error running {tool_name}: {e}")
        return tool_name, -1, str(e)

def execute_task(task):
    return task()

def move_misplaced_outputs(uniprot_id, batch_no, results_dir_base, timestamp, logger):
    """Move misplaced mhci_out and mhcii_out files to the correct output directory."""
    misplaced_files = {
        'mhc1': ROOT_DIR / 'mhci_out.csv',
        'mhc2': ROOT_DIR / 'mhcii_out.csv'
    }
    for tool_name, src_path in misplaced_files.items():
        if src_path.exists():
            dest_path = generate_tool_output_path(
                timestamp=timestamp,
                strategy_no=1,
                batch_no=batch_no,
                uniprot_id=uniprot_id,
                tool_name=tool_name,
                results_dir_base=results_dir_base
            )
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                shutil.move(src_path, dest_path)
                logger.info(f"Moved misplaced {tool_name} output from {src_path} to {dest_path}")
                print_status(f"Moved {tool_name} output to {dest_path}", "success")
            except Exception as e:
                logger.error(f"Failed to move {tool_name} output from {src_path} to {dest_path}: {e}")
                print_status(f"Error moving {tool_name} output: {e}", "error")

def strategy_1(input_file, results_dir, pathogen_type, logger, batch_files, batch_dirs, uniprot_mapping, failed_tools=None, selected_tools=None, _mhci_method=None, _mhcii_method=None):
    if failed_tools is None:
        failed_tools = []
    
    # Prompt for tool selection if none provided
    if selected_tools is None:
        print("\nChoose execution mode:")
        print("1. Run all tools (full strategy)")
        print("2. Select individual tools")
        choice = input("Enter your choice (1 or 2): ").strip()
        
        if choice == "1":
            selected_tools = ['mhc1', 'mhc2', 'netctl', 'netchop', 'bcell', 'psortb']
        elif choice == "2":
            print("\nSelect tools to run (enter numbers, separated by commas, e.g., 1,3,5):")
            tool_options = {
                '1': 'mhc1',
                '2': 'mhc2',
                '3': 'netctl',
                '4': 'netchop',
                '5': 'bcell',
                '6': 'psortb'
            }
            for num, tool in tool_options.items():
                tool_names = {
                    'mhc1': 'MHC-I (mhc_i.py)',
                    'mhc2': 'MHC-II (mhc_ii.py)',
                    'netctl': 'NetCTL (netctl.py)',
                    'netchop': 'Proteasomal Cleavage (netchop.py)',
                    'bcell': 'B-cell (IEDB) (bcell.py)',
                    'psortb': 'PSORTb (psortb.py)'
                }
                print(f"{num}. {tool_names[tool]}")
            selected = input("Enter tool numbers: ").strip().split(',')
            selected_tools = [tool_options[num.strip()] for num in selected if num.strip() in tool_options]
            if not selected_tools:
                logger.warning("No valid tools selected. Defaulting to all tools.")
                print_status("No valid tools selected. Defaulting to all tools.", "warning")
                selected_tools = ['mhc1', 'mhc2', 'netctl', 'netchop', 'bcell', 'psortb']
        else:
            logger.warning("Invalid choice. Defaulting to all tools.")
            print_status("Invalid choice. Defaulting to all tools.", "warning")
            selected_tools = ['mhc1', 'mhc2', 'netctl', 'netchop', 'bcell', 'psortb']

    logger.info("Executing Strategy 1: Epitope Prediction")
    print_status("\n## Executing Strategy 1: Epitope Prediction")
    print("Tools to be run per sequence:")
    for tool in selected_tools:
        tool_names = {
            'mhc1': 'MHC-I (mhc_i.py)',
            'mhc2': 'MHC-II (mhc_ii.py)',
            'netctl': 'NetCTL (netctl.py)',
            'netchop': 'Proteasomal Cleavage (netchop.py)',
            'bcell': 'B-cell (IEDB) (bcell.py)',
            'psortb': 'PSORTb (psortb.py)'
        }
        if tool in tool_names:
            print(f"- {tool_names[tool]}")
        else:
            logger.warning(f"Unknown tool selected: {tool}")
            print_status(f"Unknown tool selected: {tool}", "warning")

    mhci_key, mhcii_name = select_mhci_method(_choice=_mhci_method)
    mhcii_key, mhcii_name = select_mhcii_method(_choice=_mhcii_method)

    for batch_idx, batch_dir in enumerate(batch_dirs, 1):
        batch_no = batch_idx
        logger.info(f"Processing batch {batch_no}: {batch_dir}")
        print_status(f"\nProcessing batch {batch_no} in {batch_dir}")

        tasks = []
        for uniprot_id, info in uniprot_mapping.items():
            if info['batch_no'] != batch_no:
                continue
            tools = [
                ('mhc1', partial(run_mhc1, mhci_key, info['per_seq_file'])),
                ('mhc2', partial(run_mhc2, mhcii_key, info['per_seq_file'])),
                ('netctl', partial(run_netctl, info['per_seq_file'])),
                ('netchop', partial(run_netchop, info['per_seq_file'])),
                ('bcell', partial(run_bcell, info['per_seq_file'])),
                ('psortb', partial(run_psortb, info['per_seq_file'])),
            ]
            for tool_name, tool_func in tools:
                if tool_name not in selected_tools:
                    continue
                output_file = generate_tool_output_path(
                    timestamp=None,
                    strategy_no=1,
                    batch_no=batch_no,
                    uniprot_id=uniprot_id,
                    tool_name=tool_name,
                    results_dir_base=results_dir
                )
                output_file.parent.mkdir(parents=True, exist_ok=True)
                logger.debug(f"Task for {tool_name} (UniProt: {uniprot_id}): output_file={output_file}")
                tasks.append(
                    partial(run_tool_wrapper, tool_func, tool_name, output_file=output_file)
                )

        if not tasks:
            logger.warning(f"No tasks generated for batch {batch_no}")
            continue

        num_processes = min(mp.cpu_count(), len(tasks))
        logger.info(f"Running {len(tasks)} tasks in parallel with {num_processes} processes for batch {batch_no}")
        with mp.Pool(processes=num_processes) as pool:
            results = pool.map(execute_task, tasks)

        for tool_name, status, error in results:
            if status == 0:
                logger.info(f"{tool_name} completed successfully for batch {batch_no}")
                print_status(f"{tool_name} completed successfully for batch {batch_no}", "success")
                # Check for misplaced MHC-I and MHC-II outputs
                if tool_name in ['mhc1', 'mhc2']:
                    for uniprot_id, info in uniprot_mapping.items():
                        if info['batch_no'] == batch_no:
                            move_misplaced_outputs(uniprot_id, batch_no, results_dir, None, logger)
            else:
                logger.error(f"{tool_name} failed for batch {batch_no}: {error}")
                print_status(f"{tool_name} failed for batch {batch_no}: {error or 'Unknown error'}", "error")
                failed_tools.append(f"{tool_name}_batch_{batch_no}")

    logger.info("Strategy 1 completed.")
    print_status("Strategy 1 completed.", "success")
    return failed_tools

from functools import partial
import multiprocessing as mp
from pathlib import Path
from typing import List, Optional, Dict, Any

def strategy_2(
    input_file: Path,
    results_dir: Path,
    input_value: str,
    is_uniprot: bool,
    logger: Any,
    batch_files: List[Path],
    batch_dirs: List[Path],
    uniprot_mapping: Dict[str, Dict[str, Any]],
    failed_tools: Optional[List[str]] = None,
    selected_tools: Optional[List[str]] = None
) -> List[str]:
    if failed_tools is None:
        failed_tools = []

    available_tools = ['iapred', 'algpred', 'instability', 'molwt', 'wolfpsort']
    tool_func_map = {
        'iapred': run_iapred,
        'algpred': run_algpred,
        'instability': run_instability,
        'molwt': run_molwt,
        'wolfpsort': run_wolf_psort
    }

    print("\nDo you want to run:")
    print("1. Full strategy (multiple tools on all batches)")
    print("2. A single tool individually")
    choice = input("Enter 1 or 2: ").strip()

    if choice not in ['1', '2']:
        logger.error("Invalid choice. Please enter 1 or 2.")
        return failed_tools

    tools_to_run = []

    if choice == "2":
        print("Available tools: " + ", ".join(available_tools))
        tool_choice = input("Enter the tool name you want to run: ").strip().lower()
        if tool_choice not in available_tools:
            logger.error(f"Invalid tool: {tool_choice}")
            print_status(f"Invalid tool name. Available options: {', '.join(available_tools)}", "error")
            return failed_tools
        tools_to_run = [tool_choice]
    else:
        tools_to_run = available_tools
        print_status("\n## Running full strategy: all tools on all batches", "info")

    for tool_choice in tools_to_run:
        print_status(f"\n## Running tool: {tool_choice}")
        tool_output_dir = results_dir / tool_choice
        tool_output_dir.mkdir(parents=True, exist_ok=True)

        for batch_idx, batch_file in enumerate(batch_files, 1):
            print_status(f"\nProcessing batch {batch_idx} for {tool_choice}")
            batch_uniprot_ids = [uid for uid, info in uniprot_mapping.items() if info['batch_file'] == batch_file]

            tasks = []
            if tool_choice == 'wolfpsort':
                output_file = tool_output_dir / 'combined_results.txt'
                if output_file.exists():
                    logger.info(f"{tool_choice} already run for batch {batch_idx}. Skipping.")
                    continue
                func = partial(tool_func_map[tool_choice], input_fasta=batch_file, output_file=output_file)
                tasks.append(partial(run_tool_wrapper, func, tool_choice, output_file=output_file))
            else:
                for uniprot_id in batch_uniprot_ids:
                    output_file = tool_output_dir / f"{uniprot_id}_{tool_choice}.csv"
                    if output_file.exists():
                        logger.info(f"{tool_choice} already run for {uniprot_id}. Skipping.")
                        continue
                    if tool_choice == 'iapred':
                        func = partial(tool_func_map[tool_choice], batch_file, uniprot_id=input_value if is_uniprot else None, output_dir=tool_output_dir)
                    else:
                        func = partial(tool_func_map[tool_choice], input_fasta=batch_file, output_dir=tool_output_dir)
                    tasks.append(partial(run_tool_wrapper, func, tool_choice, output_file=output_file))

            if tasks:
                num_processes = min(mp.cpu_count(), len(tasks))
                with mp.Pool(processes=num_processes) as pool:
                    results = pool.map(execute_task, tasks)

                for _, status, error in results:
                    if status == 0:
                        print_status(f"{tool_choice} completed successfully for batch {batch_idx}", "success")
                    else:
                        print_status(f"{tool_choice} failed for batch {batch_idx}: {error}", "error")
                        failed_tools.append(f"{tool_choice}_batch_{batch_idx}")

    print_status(f"Tool execution completed: {', '.join(tools_to_run)}", "success")
    return failed_tools


def strategy_3(input_file, pathogen_type, results_dir, logger, batch_files, batch_dirs, uniprot_mapping, selected_tool=None, failed_tools=None):
    if failed_tools is None:
        failed_tools = []

    logger.info("Executing Strategy 3: Virulence, Toxicity, and Host Compatibility Assessment")
    print_status("\n## Executing Strategy 3: Virulence, Toxicity, and Host Compatibility Assessment")
    print(f"Pathogen type: {pathogen_type}")

    # Define available tools for Strategy 3
    all_tools = {
        'toxinpred': "ToxinPred (toxicity.py)",
        'deeptmhmm': "DeepTMHMM (deeptmhmm.py)",
        'netsol': "NetSol (netsol.py)",
        'signalp': "SignalP (signalp.py)",
        'human': "Human (human.py)"
    }
    if pathogen_type in ["virus", "bacteria"]:
        all_tools.update({
            'gutflora': "GutFlora (gutflora.py)",
            'virulence': "Virulence (virulence.py)"
        })

    # Interactive prompt for strategy or tool selection
    if selected_tool is None:
        print("\nDo you want to run:")
        print("1. Full strategy (multiple tools on all batches)")
        print("2. A single tool individually")
        user_choice = input("Enter 1 or 2: ").strip()
        
        if user_choice == "2":
            print("\nAvailable tools:")
            for tool_name, tool_desc in all_tools.items():
                print(f"- {tool_name}: {tool_desc}")
            print("\nEnter the name of the tool to run:")
            user_input = input().strip().lower()
            selected_tool = user_input if user_input else None
        elif user_choice != "1":
            logger.error(f"Invalid choice: {user_choice}. Please select 1 or 2.")
            print_status(f"Error: Invalid choice '{user_choice}'. Please select 1 or 2.", "error")
            return failed_tools

    if selected_tool and selected_tool not in all_tools:
        logger.error(f"Invalid tool selected: {selected_tool}. Available tools: {list(all_tools.keys())}")
        print_status(f"Error: Invalid tool '{selected_tool}'. Available tools: {list(all_tools.keys())}", "error")
        return failed_tools

    print(f"Selected tool: {all_tools[selected_tool]}" if selected_tool else "Running all tools per batch:")
    if not selected_tool:
        for tool_desc in all_tools.values():
            print(f"- {tool_desc}")

    # Use provided results_dir directly (do not nest strategy3 folder again)
    strategy_dir = Path(results_dir)
    strategy_dir.mkdir(parents=True, exist_ok=True)

    for batch_idx, batch_file in enumerate(batch_files, 1):
        logger.info(f"Processing batch {batch_idx}: {batch_file}")
        print_status(f"\nProcessing batch {batch_idx} with {batch_file}")

        batch_uniprot_ids = [uid for uid, info in uniprot_mapping.items() if info['batch_file'] == batch_file]

        tools = [
            ('toxinpred', partial(run_toxinpred, fasta_path=batch_file)),
            ('deeptmhmm', partial(run_deeptmhmm, input_fasta=batch_file)),
            ('netsol', partial(run_netsol, input_fasta=batch_file)),
            ('signalp', partial(run_signalp, input_fasta=batch_file)),
            ('human', partial(run_human, input_fasta=batch_file)),
        ]
        if pathogen_type in ["virus", "bacteria"]:
            tools.extend([
                ('gutflora', partial(run_gutflora, input_fasta=batch_file)),
                ('virulence', partial(run_virulence, input_fasta=batch_file)),
            ])

        if selected_tool:
            tools = [(name, func) for name, func in tools if name == selected_tool]

        tasks = []
        for tool_name, tool_func in tools:
            tool_dir = strategy_dir / tool_name
            tool_dir.mkdir(parents=True, exist_ok=True)

            if tool_name in ['deeptmhmm', 'signalp', 'netsol', 'virulence']:
                tasks.append(
                    partial(run_tool_wrapper, tool_func, tool_name, output_dir=tool_dir)
                )
            else:
                for uniprot_id in batch_uniprot_ids:
                    output_file = tool_dir / f"{uniprot_id}_{tool_name}.csv"
                    if tool_name == 'toxinpred':
                        tasks.append(
                            partial(run_tool_wrapper, partial(tool_func, output_dir=tool_dir, output_file=output_file), tool_name)
                        )
                    else:
                        tasks.append(
                            partial(run_tool_wrapper, tool_func, tool_name, output_file=output_file)
                        )

        if tasks:
            num_processes = min(mp.cpu_count(), len(tasks))
            logger.info(f"Running {len(tasks)} tasks in parallel with {num_processes} processes for batch {batch_idx}")
            with mp.Pool(processes=num_processes) as pool:
                results = pool.map(execute_task, tasks)

            for tool_name, status, error in results:
                if status == 0:
                    logger.info(f"{tool_name} completed successfully for batch {batch_idx}")
                    print_status(f"{tool_name} completed successfully for batch {batch_idx}", "success")
                else:
                    logger.error(f"{tool_name} failed for batch {batch_idx} with exit code {status}: {error}")
                    print_status(f"{tool_name} failed for batch {batch_idx}: {error or 'Unknown error'}", "error")
                    failed_tools.append(f"{tool_name}_batch_{batch_idx}")

    logger.info("Strategy 3 completed.")
    print_status("Strategy 3 completed.", "success")
    return failed_tools
import shutil
from pathlib import Path

def strategy_4(input_file, pathogen_type, results_dir, logger, batch_files, batch_dirs, uniprot_mapping, failed_tools=None):
    if failed_tools is None:
        failed_tools = []

    logger.info("Executing Strategy 4: Comprehensive Multi-Tool Analysis")
    print_status("\n## Executing Strategy 4: Comprehensive Multi-Tool Analysis")
    print_status(f"Pathogen type: {pathogen_type}")
    print_status("Running all tools from Strategies 1, 2, and 3:")

    parent_dir = results_dir.parent
    results_dir.mkdir(parents=True, exist_ok=True)

    def run_and_move(strategy_func, strategy_no, *args):
        # Run in default top-level dir
        src_dir = parent_dir / f"strategy_{strategy_no}"
        src_dir.mkdir(parents=True, exist_ok=True)
        dest_dir = results_dir / f"strategy_{strategy_no}"

        try:
            # Run strategy
            if strategy_no == 1:
                strategy_func(input_file, src_dir, pathogen_type, logger, batch_files, batch_dirs, uniprot_mapping, failed_tools)
            elif strategy_no == 2:
                strategy_func(input_file, src_dir, input_file.name, not input_file.is_file(), logger, batch_files, batch_dirs, uniprot_mapping, failed_tools)
            elif strategy_no == 3:
                strategy_func(input_file, pathogen_type, src_dir, logger, batch_files, batch_dirs, uniprot_mapping, failed_tools)

            # Move it
            if src_dir.exists():
                shutil.move(str(src_dir), str(dest_dir))
                logger.info(f"Moved {src_dir} → {dest_dir}")
                print_status(f"Moved strategy_{strategy_no} results inside strategy_4", "success")
            else:
                logger.warning(f"{src_dir} does not exist, skipping move")
                print_status(f"Warning: strategy_{strategy_no} output not found for move", "warning")

        except Exception as e:
            logger.error(f"Error in strategy_{strategy_no}: {e}")
            print_status(f"Error running or moving strategy_{strategy_no}: {e}", "error")
            failed_tools.append(f"strategy_{strategy_no}")

    # Run and move all 3 strategies unconditionally
    run_and_move(strategy_1, 1)
    run_and_move(strategy_2, 2)
    run_and_move(strategy_3, 3)

    logger.info("Strategy 4 completed.")
    print_status("Strategy 4 completed.", "success")
    return failed_tools



import pandas as pd
import shutil
import os

try:
    from modules.algpred import run_algpred
    from modules.iapred import run_iapred
    from modules.toxicity import run_toxinpred
    from modules.toxin_epitope import run_toxinpred3
    from modules.algpred_down import run_algpred_down
    from modules.iapred_down import run_iapred_down
    from modules.pepmatch import run_pepmatch
except ImportError:
    from algpred import run_algpred
    from iapred import run_iapred
    from toxicity import run_toxinpred
    from toxin_epitope import run_toxinpred3
    from algpred_down import run_algpred_down
    from iapred_down import run_iapred_down
    from pepmatch import run_pepmatch
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from modules.graph import plot_epitope_analysis
import os
import shutil
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from datetime import datetime

def create_timestamped_dir(base_dir):
    """Create a directory with a timestamp in its name."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    timestamp_dir = Path(base_dir) / f"results_{timestamp}"
    timestamp_dir.mkdir(parents=True, exist_ok=True)
    return timestamp_dir
def csv_to_fasta(csv_file, fasta_file, sequence_column="Peptide", uniprot_id=None):
    """Convert CSV to FASTA, writing unique sequences with UniProt ID in header."""
    try:
        df = pd.read_csv(csv_file)
        if sequence_column not in df.columns:
            print_status(f"Warning: Column '{sequence_column}' not found in {csv_file}.", "warning")
            return False
        sequences = df[sequence_column][df.get('Score', 0) > 0.5].dropna().unique()[:100]
        if not sequences.size:
            print_status(f"Warning: No valid sequences after filtering in {csv_file}.", "warning")
            return False
        with open(fasta_file, "w") as f:
            for i, seq in enumerate(sequences):
                if seq.strip():
                    header = f">{uniprot_id}_Epitope_{i+1}" if uniprot_id else f">Epitope_{i+1}"
                    f.write(f"{header}\n{seq}\n")
        print_status(f"Wrote {len(sequences)} sequences to {fasta_file}")
        return True
    except Exception as e:
        print_status(f"Error converting {csv_file} to FASTA: {e}", "error")
        return False

def is_csv_valid(csv_path):
    """Check if a CSV file exists, is non-empty, and contains valid data."""
    if not csv_path.exists():
        return False
    try:
        df = pd.read_csv(csv_path)
        return not df.empty and df.shape[0] > 0
    except Exception:
        return False

import pandas as pd
import shutil
try:
    from modules.algpred import run_algpred
    from modules.iapred import run_iapred
    from modules.toxicity import run_toxinpred
    from modules.toxin_epitope import run_toxinpred3
    from modules.algpred_down import run_algpred_down
    from modules.iapred_down import run_iapred_down
    from modules.pepmatch import run_pepmatch
except ImportError:
    from algpred import run_algpred
    from iapred import run_iapred
    from toxicity import run_toxinpred
    from toxin_epitope import run_toxinpred3
    from algpred_down import run_algpred_down
    from iapred_down import run_iapred_down
    from pepmatch import run_pepmatch
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from modules.graph import plot_epitope_analysis

def strategy_5(input_file, pathogen_type, results_dir, logger, batch_files, batch_dirs, uniprot_mapping, uniprot_id_to_extract=None, failed_tools=None, _pre_predicted_fastas=None, _mhci_method=None, _mhcii_method=None):
    if failed_tools is None:
        failed_tools = []
    logger.info("Executing Strategy 5: Sequential Epitope Analysis with Pre-Filtering")
    print_status("\n## Executing Strategy 5: Sequential Epitope Analysis with Pre-Filtering")
    print_status(f"Pathogen type: {pathogen_type}")
    print_status("Workflow:")
    print_status("- Check for pre-predicted epitope FASTA files")
    print_status("- Pre-analysis: AlgPred, IAPred, ToxinPred")
    print_status("- Filter sequences (non-allergen, antigenic, non-toxic)")
    print_status("- B-cell epitope prediction (bcell.py) → ClbTope (clbtope.py) → Downstream analysis")
    print_status("- MHC-I prediction (mhc_i.py) → Downstream analysis")
    print_status("- MHC-II prediction (mhc_ii.py) → Downstream analysis")
    print_status("Downstream tools: Allergenicity (algpred.py), Antigenicity (iapred.py), Toxicity (toxinpred3.py), Peptide Matching (pepmatch.py)")

    check_dependencies()

    timestamp_dir = create_timestamped_dir(str(results_dir))
    logger.debug(f"Created directory: {timestamp_dir}")

    # Programmatic bypass for pre-predicted FASTA paths
    if _pre_predicted_fastas is not None:
        have_predicted = "yes" if _pre_predicted_fastas else "no"
        bcell_fasta = _pre_predicted_fastas.get("bcell", "")
        mhci_fasta = _pre_predicted_fastas.get("mhci", "")
        mhcii_fasta = _pre_predicted_fastas.get("mhcii", "")
    else:
        have_predicted = input("Do you have pre-predicted epitope FASTA files for B-cell, MHC-I, and MHC-II? (yes/no): ").strip().lower()
        while have_predicted not in ["yes", "no"]:
            print_status("Please enter 'yes' or 'no'.", "error")
            have_predicted = input("Do you have pre-predicted FASTA files for B-cell, MHC-I, and MHC-II? ").strip().lower()

    epitope_types = ["bcell", "mhci", "mhcii"]

    if have_predicted == "yes":
        if _pre_predicted_fastas is None:
            bcell_fasta = input("Enter the path to B-cell epitope FASTA file: ").strip()
            mhci_fasta = input("Enter the path to MHC-I epitope FASTA file: ").strip()
            mhcii_fasta = input("Enter the path to MHC-II epitope FASTA file: ").strip()

        epitope_fastas = []
        for fasta, name in [(bcell_fasta, "bcell"), (mhci_fasta, "mhci"), (mhcii_fasta, "mhcii")]:
            try:
                validate_fasta_file(fasta)
                logger.debug(f"Validated FASTA file: {fasta}")
                epitope_fastas.append((fasta, name))
            except Exception as e:
                logger.error(f"FASTA file validation failed: {fasta}: {e}")
                print_status(f"Error: Invalid FASTA file {fasta}: {e}", "error")
                failed_tools.append(f"Validation_{name}")
                continue

        print_status(f"\n[Downstream Analysis] Processing pre-predicted epitope FASTA files...")
        for epitope_fasta, epitope_type in epitope_fastas:
            uniprot_id = uniprot_id_to_extract or "unknown"
            try:
                output_base = timestamp_dir
                output_base.mkdir(parents=True, exist_ok=True)

                tool_outputs = {
                    "ToxinPred3": output_base / f"{epitope_type}_Combined_ToxinPred3_1.csv",
                    "AlgPred": output_base / f"{epitope_type}_Combined_AlgPred_1.csv",
                    "IAPred": output_base / f"{epitope_type}_combined_iapred.csv",
                    "PepMatch": output_base / f"{epitope_type}_Combined_PepMatch_1.csv"
                }

                run_toxinpred3(epitope_fasta, tool_outputs["ToxinPred3"])
                run_algpred_down(epitope_fasta, output_base, output_file=tool_outputs["AlgPred"].name)
                run_iapred_down(epitope_fasta, output_base, output_file=tool_outputs["IAPred"].name)
                run_pepmatch(epitope_fasta, output_base, output_file=tool_outputs["PepMatch"])

                for tool, out_file in tool_outputs.items():
                    if out_file.exists():
                        logger.debug(f"{tool} output saved to {out_file}")
                        print_status(f"[✓] {tool} analysis for {epitope_type} completed successfully", "success")
                    else:
                        logger.warning(f"{tool} output file {out_file} not found")
                        print_status(f"Error: {tool} failed for {epitope_type}: Output file not found", "error")
                        failed_tools.append(f"{tool}_{epitope_type}")
            except Exception as e:
                logger.error(f"{uniprot_id} - Downstream failed for {epitope_type}: {e}", exc_info=True)
                print_status(f"Error: Downstream failed for {epitope_type}: {e}", "error")
                failed_tools.append(f"{uniprot_id}_{epitope_type}_downstream")

        plot_epitope_analysis(
            results_dir=timestamp_dir,
            uniprot_ids=[uniprot_id_to_extract] if uniprot_id_to_extract else ["unknown"],
            epitope_types=epitope_types,
            timestamp_dir=timestamp_dir,
            batch_no=1
        )
        return failed_tools

    for batch_no, (batch_file, pre_dir) in enumerate(zip(batch_files, batch_dirs), 1):
        logger.debug(f"Processing batch {batch_no}: file={batch_file}, pre_dir={pre_dir}")
        try:
            validate_fasta_file(batch_file)
            logger.debug(f"Validated batch file: {batch_file}")
        except Exception as e:
            logger.error(f"Batch file validation failed: {batch_file}: {e}")
            print_status(f"Error: Invalid batch file {batch_file}: {e}", "error")
            failed_tools.append(f"Batch_validation_batch{batch_no}")
            continue

        print_status(f"\n[Pre-Analysis] Running AlgPred, IAPred, and ToxinPred for batch {batch_no}...")
        pre_analysis_success = True
        for tool, func in [
            ("AlgPred", run_algpred),
            ("IAPred", run_iapred),
            ("ToxinPred", run_toxinpred)
        ]:
            try:
                output_file = func(str(batch_file), timestamp_dir)
                logger.debug(f"{tool} output written to: {output_file}")
            except Exception as e:
                logger.error(f"Pre-analysis tool {tool} failed for batch {batch_no}: {e}")
                print_status(f"Error: {tool} failed for batch {batch_no}: {e}", "error")
                failed_tools.append(f"{tool}_batch{batch_no}")
                pre_analysis_success = False

        if not pre_analysis_success:
            logger.warning(f"Skipping batch {batch_no} due to pre-analysis failures")
            print_status(f"Warning: Skipping batch {batch_no} due to pre-analysis failures", "warning")
            continue

        tox_subdir = timestamp_dir / "toxinpred"
        tox_csv_expected = timestamp_dir / f"combined_toxinpred_batch_{batch_no}.csv"
        tox_csv_source = tox_subdir / "combined.csv"
        if tox_subdir.exists() and tox_csv_source.exists():
            try:
                shutil.copy(tox_csv_source, tox_csv_expected)
                logger.debug(f"Copied ToxinPred output from {tox_csv_source} to {tox_csv_expected}")
            except Exception as e:
                logger.error(f"Failed to copy ToxinPred CSV from {tox_csv_source} to {tox_csv_expected}: {e}")
                print_status(f"Error: Failed to copy ToxinPred CSV: {e}", "error")
                failed_tools.append(f"CSV_copy_ToxinPred_batch{batch_no}")
        else:
            logger.error(f"ToxinPred output missing: {tox_csv_source}")
            print_status(f"Error: ToxinPred output missing: {tox_csv_source}", "error")
            failed_tools.append(f"ToxinPred_output_missing_batch{batch_no}")

        alg_csv = timestamp_dir / f"combined_algpred_batch_{batch_no}.csv"
        iap_csv = timestamp_dir / f"combined_iapred_batch_{batch_no}.csv"
        tox_csv = tox_csv_expected

        all_csvs_present = True
        for csv_file, tool_name in [(alg_csv, "AlgPred"), (iap_csv, "IAPred"), (tox_csv, "ToxinPred")]:
            if not csv_file.exists():
                logger.error(f"Pre-analysis CSV not found: {csv_file}")
                print_status(f"Error: Pre-analysis CSV not found: {csv_file}", "error")
                failed_tools.append(f"CSV_missing_{tool_name}_batch{batch_no}")
                all_csvs_present = False

        if not all_csvs_present:
            logger.warning(f"Skipping batch {batch_no} due to missing CSV files")
            print_status(f"Warning: Skipping batch {batch_no} due to missing CSV files", "warning")
            continue

        input_records = list(SeqIO.parse(batch_file, "fasta"))
        seq_dict = {}
        for rec in input_records:
            seq_dict[rec.id] = str(rec.seq)
            description_parts = rec.description.split("|")
            for part in description_parts:
                if part and part not in seq_dict:
                    seq_dict[part] = str(rec.seq)
            if len(description_parts) > 1 and description_parts[0].startswith("sp"):
                uniprot_id = description_parts[1]
                seq_dict[uniprot_id] = str(rec.seq)
        logger.debug(f"Sequence dictionary keys: {list(seq_dict.keys())}")

        try:
            dfs = []
            for csv_file, tool_name, out_col in [
                (alg_csv, "AlgPred", "AlgPred_Prediction"),
                (iap_csv, "IAPred", "Antigenicity_Category"),
                (tox_csv, "ToxinPred", "ToxinPred_Prediction")
            ]:
                try:
                    df = pd.read_csv(csv_file)
                    if df.empty:
                        logger.error(f"CSV file is empty: {csv_file}")
                        print_status(f"Error: CSV file is empty: {csv_file}", "error")
                        failed_tools.append(f"CSV_empty_{tool_name}_batch{batch_no}")
                        all_csvs_present = False
                        continue
                    id_col = next((col for col in df.columns if col.lower() in ["seq_id", "id", "sequence_id", "sequence", "identifier", "header", "subject"]), None)
                    if not id_col:
                        for col in df.columns:
                            if df[col].dtype == "object" and df[col].str.match(r"^[A-Za-z0-9_|\-]+$").all():
                                id_col = col
                                logger.debug(f"Inferred ID column for {tool_name}: {id_col}")
                                break
                        if not id_col:
                            logger.error(f"No valid ID column found in {tool_name} CSV: {list(df.columns)}")
                            print_status(f"Error: No valid ID column found in {tool_name} CSV: {list(df.columns)}", "error")
                            failed_tools.append(f"{tool_name}_id_column_batch{batch_no}")
                            all_csvs_present = False
                            continue
                    pred_col = "Prediction" if tool_name == "AlgPred" else next((col for col in df.columns if col.lower() in ["prediction", "result", "output", "label", "class", "antigenicity_category"]), None)
                    if not pred_col:
                        logger.error(f"Missing prediction column in {tool_name} CSV: {list(df.columns)}")
                        print_status(f"Error: Missing prediction column in {tool_name} CSV: {list(df.columns)}", "error")
                        failed_tools.append(f"{tool_name}_columns_batch{batch_no}")
                        all_csvs_present = False
                        continue
                    df = df[[id_col, pred_col]].rename(columns={id_col: "seq_id", pred_col: out_col})
                    if tool_name in ["IAPred", "ToxinPred"]:
                        df["seq_id"] = df["seq_id"].apply(lambda x: x.split("|")[1] if isinstance(x, str) and "|" in x else x)
                    if tool_name == "AlgPred":
                        df["sequence"] = df["seq_id"].map(seq_dict)
                    dfs.append((df, tool_name))
                except Exception as e:
                    logger.error(f"Failed to read {tool_name} CSV {csv_file}: {e}")
                    print_status(f"Error: Failed to read {tool_name} CSV {csv_file}: {e}", "error")
                    failed_tools.append(f"{tool_name}_read_batch{batch_no}")
                    all_csvs_present = False

            if not all_csvs_present or not dfs:
                logger.warning(f"Skipping batch {batch_no} due to CSV reading or column issues")
                print_status(f"Warning: Skipping batch {batch_no} due to CSV reading or column issues", "warning")
                continue

            alg_df, iap_df, tox_df = [df for df, _ in dfs]
            merged_df = alg_df.merge(iap_df, on="seq_id", how="inner").merge(tox_df, on="seq_id", how="inner")
            if merged_df.empty:
                logger.warning(f"Merged DataFrame is empty for batch {batch_no}")
                print_status(f"Warning: Merged DataFrame is empty for batch {batch_no}", "warning")
                continue
            check_csv_path = timestamp_dir / f"pre_analysis_result_batch_{batch_no}.csv"
            merged_df.to_csv(check_csv_path, index=False)
            print_status(f"[✓] Pre-analysis results saved to {check_csv_path}", "success")

            print_status(f"\n[Filtering] Selecting sequences that are non-allergen, antigenic, and non-toxic for batch {batch_no}...")
            filtered_ids = merged_df[
                (merged_df["AlgPred_Prediction"] == "Non-Allergen") &
                (merged_df["Antigenicity_Category"] == "High") &
                (merged_df["ToxinPred_Prediction"] == "Non-Toxin")
            ]["seq_id"].tolist()
            if not filtered_ids:
                logger.warning(f"No sequences passed strict filters for batch {batch_no}, using all sequences from merged_df")
                print_status(f"Warning: No sequences passed strict filters, using all sequences", "warning")
                filtered_ids = merged_df["seq_id"].tolist()

            filtered_fasta = timestamp_dir / f"filtered_sequences_batch_{batch_no}.fasta"
            written_sequences = 0
            with open(filtered_fasta, "w") as f:
                for _, row in merged_df[merged_df["seq_id"].isin(filtered_ids)].iterrows():
                    seq_id = row["seq_id"]
                    sequence = row.get("sequence", None)
                    if pd.notna(sequence):
                        f.write(f">{seq_id}\n{sequence}\n")
                        written_sequences += 1
            logger.debug(f"Filtered FASTA file created: {filtered_fasta} with {written_sequences} sequences")
            print_status(f"[✓] Successfully filtered {written_sequences} sequences to {filtered_fasta}", "success")
            if not written_sequences:
                logger.warning(f"No sequences written to filtered FASTA file: {filtered_fasta}")
                print_status(f"Warning: No sequences written to filtered FASTA file: {filtered_fasta}", "warning")
                continue
        except Exception as e:
            logger.error(f"Batch {batch_no} - Failed to process CSVs or filter sequences: {e}")
            print_status(f"Error: Batch processing failed for batch {batch_no}: {e}", "error")
            failed_tools.append(f"Batch_processing_batch{batch_no}")
            continue

        mhci_key, mhci_name = select_mhci_method(_choice=_mhci_method)
        mhcii_key, mhcii_name = select_mhcii_method(_choice=_mhcii_method)

        batch_dir = timestamp_dir / f"batch{batch_no}"
        batch_dir.mkdir(parents=True, exist_ok=True)

        batch_uniprot_ids = [uid for uid, info in uniprot_mapping.items() if info['batch_file'] == batch_file]
        strategy5_fasta = filtered_fasta

        bcell_outputs = []
        for uniprot_id in batch_uniprot_ids:
            bcell_csv = timestamp_dir / f"{uniprot_id}_bcell_out.csv"
            try:
                if not strategy5_fasta.exists():
                    raise FileNotFoundError(f"Input FASTA not found: {strategy5_fasta}")
                logger.debug(f"Running B-cell prediction for {strategy5_fasta} -> {bcell_csv}")
                run_bcell(strategy5_fasta, output_file=bcell_csv)
                if not bcell_csv.exists():
                    raise FileNotFoundError(f"B-cell output file not created: {bcell_csv}")
                logger.debug(f"B-cell output generated: {bcell_csv}")
                bcell_outputs.append(bcell_csv)
            except Exception as e:
                logger.error(f"{uniprot_id} - B-cell prediction failed: {e}", exc_info=True)
                failed_tools.append(f"{uniprot_id}_bcell")
                continue

            bcell_input_fasta = timestamp_dir / f"{uniprot_id}_bcell_input.fasta"
            clbtope_csv = timestamp_dir / f"{uniprot_id}_clbtope_out.csv"
            bcell_epitope_fasta = timestamp_dir / f"{uniprot_id}_bcell_epitope.fasta"

            try:
                if not bcell_csv.exists():
                    raise FileNotFoundError(f"B-cell CSV not found: {bcell_csv}")
                df_b = pd.read_csv(bcell_csv)
                logger.debug(f"B-cell CSV columns: {df_b.columns}")
                if 'Peptipe' in df_b.columns:
                    logger.debug(f"Renaming 'Peptipe' to 'Peptide' in {bcell_csv}")
                    df_b.rename(columns={'Peptipe': 'Peptide'}, inplace=True)
                if 'Peptide' not in df_b.columns:
                    raise ValueError(f"'Peptide' column not found in {bcell_csv}")
                peptides = df_b['Peptide'].dropna().unique()
                logger.debug(f"Found {len(peptides)} unique peptides: {peptides}")
                if not peptides.size:
                    raise ValueError(f"No valid peptides found in {bcell_csv}")
                logger.debug(f"Writing {len(peptides)} peptides to {bcell_input_fasta}")
                SeqIO.write(
                    [SeqRecord(Seq(p), id=f"{uniprot_id}_Bcell_{i+1}", description="") for i, p in enumerate(peptides)],
                    bcell_input_fasta, "fasta"
                )
                if not bcell_input_fasta.exists():
                    raise FileNotFoundError(f"B-cell input FASTA not created: {bcell_input_fasta}")
                logger.debug(f"Running ClbTope on {bcell_input_fasta}")
                run_clbtope(bcell_input_fasta, clbtope_csv)
                if not clbtope_csv.exists():
                    raise FileNotFoundError(f"ClbTope output not created: {clbtope_csv}")
                df_clb = pd.read_csv(clbtope_csv)
                logger.debug(f"ClbTope CSV columns: {df_clb.columns}")
                if 'Sequence' not in df_clb.columns:
                    raise ValueError(f"'Sequence' column not found in {clbtope_csv}")
                final_peptides = df_clb['Sequence'].dropna().unique()
                logger.debug(f"Found {len(final_peptides)} final peptides: {final_peptides}")
                logger.debug(f"Writing {len(final_peptides)} final peptides to {bcell_epitope_fasta}")
                SeqIO.write(
                    [SeqRecord(Seq(p), id=f"{uniprot_id}_BcellEpitope_{i+1}", description="") for i, p in enumerate(final_peptides)],
                    bcell_epitope_fasta, "fasta"
                )
                if not bcell_epitope_fasta.exists():
                    raise FileNotFoundError(f"B-cell epitope FASTA not created: {bcell_epitope_fasta}")
            except Exception as e:
                logger.error(f"{uniprot_id} - ClbTope or FASTA generation failed: {e}", exc_info=True)
                failed_tools.append(f"{uniprot_id}_clbtope")
                bcell_epitope_fasta = None

            mhci_csv = timestamp_dir / f"{uniprot_id}_mhci_out.csv"
            mhcii_csv = timestamp_dir / f"{uniprot_id}_mhcii_out.csv"

            try:
                run_mhc1(mhci_key, strategy5_fasta, output_file=mhci_csv)
                run_mhc2(mhcii_key, strategy5_fasta, output_file=mhcii_csv)
            except Exception as e:
                logger.error(f"{uniprot_id} - MHC prediction failed: {e}", exc_info=True)
                failed_tools.append(f"{uniprot_id}_mhc_prediction")

            def extract_fasta(csv_path, col_name, output_fasta, epitope_tag):
                try:
                    df = pd.read_csv(csv_path)
                    peptides = df[col_name].dropna().unique()
                    SeqIO.write(
                        [SeqRecord(Seq(p), id=f"{uniprot_id}_{epitope_tag}_{i+1}", description="") for i, p in enumerate(peptides)],
                        output_fasta, "fasta"
                    )
                    return output_fasta
                except Exception as e:
                    logger.error(f"{uniprot_id} - Failed to extract {epitope_tag} FASTA: {e}", exc_info=True)
                    failed_tools.append(f"{uniprot_id}_{epitope_tag}_fasta")
                    return None

            mhci_fasta = extract_fasta(mhci_csv, "peptide", timestamp_dir / f"{uniprot_id}_mhci_epitope.fasta", "MHCI")
            mhcii_fasta = extract_fasta(mhcii_csv, "peptide", timestamp_dir / f"{uniprot_id}_mhcii_epitope.fasta", "MHCII")

            for epitope_type, fasta_file in [
                ("bcell", bcell_epitope_fasta),
                ("mhci", mhci_fasta),
                ("mhcii", mhcii_fasta)
            ]:
                if not fasta_file or not Path(fasta_file).exists():
                    logger.warning(f"{uniprot_id} - Missing FASTA for {epitope_type}, skipping downstream.")
                    continue
                try:
                    output_base = timestamp_dir
                    output_base.mkdir(parents=True, exist_ok=True)

                    tool_outputs = {
                        "ToxinPred3": output_base / f"{epitope_type}_Combined_ToxinPred3_{batch_no}.csv",
                        "AlgPred": output_base / f"{epitope_type}_Combined_AlgPred_{batch_no}.csv",
                        "IAPred": output_base / f"{epitope_type}_combined_iapred.csv",
                        "PepMatch": output_base / f"{epitope_type}_Combined_PepMatch_{batch_no}.csv"
                    }

                    run_toxinpred3(fasta_file, tool_outputs["ToxinPred3"])
                    run_algpred_down(fasta_file, output_base, output_file=tool_outputs["AlgPred"].name)
                    run_iapred_down(fasta_file, output_base, output_file=tool_outputs["IAPred"].name)
                    run_pepmatch(fasta_file, output_base, output_file=tool_outputs["PepMatch"])

                    for tool, out_file in tool_outputs.items():
                        if out_file.exists():
                            logger.debug(f"{tool} output saved to {out_file}")
                        else:
                            logger.warning(f"{tool} output file {out_file} not found")
                            failed_tools.append(f"{uniprot_id}_{epitope_type}_{tool}")
                except Exception as e:
                    logger.error(f"{uniprot_id} - Downstream failed for {epitope_type}: {e}", exc_info=True)
                    failed_tools.append(f"{uniprot_id}_{epitope_type}_downstream")

        plot_epitope_analysis(
            results_dir=timestamp_dir,
            uniprot_ids=batch_uniprot_ids,
            epitope_types=epitope_types,
            timestamp_dir=timestamp_dir,
            batch_no=batch_no
        )

    logger.info("Strategy 5 completed.")
    print_status("Strategy 5 completed.", "success")

    if failed_tools:
        logger.warning(f"Failed tools: {', '.join(failed_tools)}")
        print_status(f"Warning: The following tools failed: {', '.join(failed_tools)}", "warning")
    return failed_tools

import os
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from io import StringIO
import shutil  # only used for moving mhcii_out.csv when found

# ───────────────────────────────────────────────────────────────
# Import your REAL prediction functions
# Adjust these import paths to match your project structure
# ───────────────────────────────────────────────────────────────
from modules.mhc_i import run_mhc1
from modules.mhc_ii import run_mhc2
from modules.bcell import run_bcell

def extract_positions_from_prediction_csv(csv_path):
    """
    Extract (Start, End) from prediction CSVs (MHC-I, MHC-II, B-cell).
    Tries multiple common column patterns from IEDB/NetMHC-style outputs.
    """
    csv_path = str(csv_path)  # ensure string
    if not os.path.exists(csv_path):
        print(f"CSV path does not exist: {csv_path}")
        return []
    if os.path.isdir(csv_path):
        print(f"ERROR: {csv_path} is a DIRECTORY, not a CSV file! Skipping (tool bug).")
        return []
    if not os.path.isfile(csv_path):
        print(f"ERROR: {csv_path} exists but is not a regular file! Skipping.")
        return []
    try:
        df = pd.read_csv(csv_path)
        df.columns = [col.strip().lower() for col in df.columns]  # normalize case and whitespace
        # Preferred: explicit Start/End
        if 'start' in df.columns and 'end' in df.columns:
            return list(zip(df['start'].astype(int), df['end'].astype(int)))
        # Many IEDB/NetMHC tools use 'Pos' or 'Start Pos' for starting position
        elif 'pos' in df.columns or 'start pos' in df.columns:
            pos_col = 'pos' if 'pos' in df.columns else 'start pos'
            return [(int(row[pos_col]), int(row[pos_col]) + len(row.get('peptide', '')) - 1)
                    for _, row in df.iterrows() if pd.notna(row.get('peptide'))]
        # Fallback: peptide column + infer length
        elif 'peptide' in df.columns:
            return [(i+1, i + len(p)) for i, p in enumerate(df['peptide'].dropna())]
        print(f"Warning: No usable position columns found in {csv_path}")
        print(f"Available columns: {list(df.columns)}")
        return []
    except Exception as e:
        print(f"Error parsing prediction CSV {csv_path}: {e}")
        return []

def strategy_6(results_dir, logger, _mode=None, _assembly_config=None):
    """Strategy 6: Multi-Epitope Vaccine Assembly & 3D Structural Validation"""
    print_status("\n## Strategy 6: Multi-Epitope Vaccine Assembly & 3D Modeling", "info")
 
    # 1. INPUT MODE SELECTION
    if _mode is not None:
        mode_choice = "2" if _mode == "custom" else "1"
    else:
        print("\nHow would you like to provide the vaccine sequence?")
        print("1. Assemble from epitope CSV files (Strategy 5 outputs)")
        print("2. Upload/Paste a custom pre-assembled vaccine sequence")
        mode_choice = input("Enter choice (1 or 2): ").strip()
 
    final_vaccine_sequence = ""
 
    if mode_choice == "2":
        print("\n[Custom Input Mode]")
        print("Please provide the full vaccine construct sequence for 3D modeling.")
        if _assembly_config and _assembly_config.custom_fasta_path:
            path_or_seq = _assembly_config.custom_fasta_path
        elif _assembly_config and _assembly_config.custom_sequence:
            path_or_seq = _assembly_config.custom_sequence
        else:
            path_or_seq = input("Enter path to FASTA/TXT file OR paste the raw sequence: ").strip()
     
        path_or_seq = path_or_seq.replace("'", "").replace('"', "")
     
        if os.path.isfile(path_or_seq):
            try:
                with open(path_or_seq, 'r', encoding='utf-8-sig') as f:
                    content = f.read().strip()
             
                if content.startswith(">"):
                    records = list(SeqIO.parse(StringIO(content), "fasta"))
                    if records:
                        final_vaccine_sequence = str(records[0].seq).upper()
                    else:
                        raise ValueError("FASTA file appears empty or malformed.")
                else:
                    final_vaccine_sequence = "".join(content.split()).upper()
            except Exception as e:
                print_status(f"Error reading file: {e}", "error")
                logger.error(f"Error reading custom sequence file: {e}")
                return ["Custom_Sequence_Read_Error"]
        else:
            final_vaccine_sequence = "".join(path_or_seq.split()).upper()
     
        if not final_vaccine_sequence:
            print_status("Error: No sequence detected!", "error")
            return ["Empty_Sequence_Error"]
        
        print_status(f"Loaded sequence of length {len(final_vaccine_sequence)} aa.", "success")
      
        # ───────────────────────────────────────────────────────────────
        # Run 3D modeling FIRST
        # ───────────────────────────────────────────────────────────────
        print("\n" + "="*55)
        print(" PHASE 2: Automatic 3D Structural Validation (ESMFold)")
        print("="*55)
      
        from modules.esmfold import run_esmfold
      
        res = run_esmfold(final_vaccine_sequence, results_dir, logger)
      
        if res['status'] == 0:
            print_status("\n[✓] 3D Structural Validation Complete!", "success")
         
            conf = res.get('avg_confidence', 0)
            if isinstance(conf, (int, float)):
                print(f"Overall Model Confidence (pLDDT): {conf:.2f}%")
            else:
                print(f"Overall Model Confidence: {conf}")
             
            print(f"PDB Structure: {res['pdb_file']}")
            if 'report_file' in res:
                print(f"Confidence Report: {res['report_file']}")
         
            print_status("\nTip: PyMOL has been launched to visualize the 3D structure.", "info")
        else:
            print_status(f"3D Modeling skipped or failed: {res.get('error')}", "error")
            return ["3D_Modeling_Failure"]
      
        # ───────────────────────────────────────────────────────────────
        # Ask about SASA AFTER showing 3D results
        # ───────────────────────────────────────────────────────────────
        if _assembly_config is not None:
            do_sasa = 'y' if _assembly_config.run_sasa else 'n'
        elif _mode is not None:
            do_sasa = 'n'
        else:
            do_sasa = input("\nDo you want to perform SASA exposure analysis on the predicted structure? (y/n): ").strip().lower()
      
        if do_sasa == 'y':
            if _assembly_config is not None and getattr(_assembly_config, "sasa_csv_path", None):
                has_csv = 'y'
            elif _mode is not None:
                has_csv = 'n'
            else:
                has_csv = input("Do you have a CSV file with epitope positions? (y/n): ").strip().lower()
            epitope_positions = {}
         
            if has_csv == 'y':
                if _assembly_config and _assembly_config.sasa_csv_path:
                    csv_path = _assembly_config.sasa_csv_path
                else:
                    csv_path = input("Enter path to epitope positions CSV (columns: Type,Start,End): ").strip()
                if os.path.exists(csv_path):
                    try:
                        df = pd.read_csv(csv_path)
                        for key in ['bcell', 'mhci', 'mhcii', 'Bcell', 'CTL', 'HTL']:
                            sub = df[df['Type'].str.contains(key, case=False, na=False)]
                            if not sub.empty:
                                norm_key = 'bcell' if 'b' in key.lower() else 'mhci' if 'c' in key.lower() or 'ctl' in key.lower() else 'mhcii'
                                epitope_positions[norm_key] = list(zip(sub['Start'].astype(int), sub['End'].astype(int)))
                    except Exception as e:
                        print_status(f"Error loading CSV: {e}", "warning")
                else:
                    print_status("CSV not found.", "warning")
            else:
                print("\nNo position CSV provided → predicting epitopes on the vaccine sequence.")
                if _assembly_config is None:
                    vaccine_fasta_input = input("Path to vaccine FASTA file (Enter to use loaded sequence): ").strip()
                else:
                    vaccine_fasta_input = _assembly_config.custom_fasta_path if _assembly_config.custom_fasta_path else ""
            
                if vaccine_fasta_input and os.path.exists(vaccine_fasta_input):
                    try:
                        rec = next(SeqIO.parse(vaccine_fasta_input, "fasta"))
                        final_vaccine_sequence = str(rec.seq).upper()
                    except:
                        print_status("Could not read FASTA → using loaded sequence.", "warning")
            
                temp_fasta = results_dir / "vaccine_for_prediction.fasta"
                with open(temp_fasta, "w") as f:
                    f.write(f">Vaccine_construct\n{final_vaccine_sequence}\n")
            
                # Output files
                bcell_out = results_dir / "predicted_bcell.csv"
                mhci_out = results_dir / "predicted_mhci.csv"
                mhcii_out = results_dir / "predicted_mhcii.csv"
            
                print_status("Predicting B-cell, MHC-I and MHC-II epitopes using your tools...", "info")
            
                try:
                    # B-cell prediction
                    print(f"DEBUG: Calling B-cell tool → output path: {bcell_out}")
                    run_bcell(str(temp_fasta), str(bcell_out))
                    print(f"DEBUG B-cell finished → exists: {os.path.exists(bcell_out)} | is_file: {os.path.isfile(bcell_out)} | is_dir: {os.path.isdir(bcell_out)}")
                 
                    # MHC-I prediction
                    print(f"DEBUG: Calling MHC-I tool → output path: {mhci_out}")
                    run_mhc1("f", str(temp_fasta), str(mhci_out))
                    print(f"DEBUG MHC-I finished → exists: {os.path.exists(mhci_out)} | is_file: {os.path.isfile(mhci_out)} | is_dir: {os.path.isdir(mhci_out)}")
                 
                    # MHC-II prediction — FIXED CALLING CONVENTION
                    print(f"DEBUG: Calling MHC-II tool with 'nmel'...")
                    # Use keyword arguments to match your function signature
                    run_mhc2(
                        method_code="nmel",
                        fasta_file=str(temp_fasta),
                        output_file=str(mhcii_out)
                    )
                    print(f"DEBUG: MHC-II call completed — checking output...")
                   
                    # ROBUST FILE DETECTION for MHC-II
                    mhcii_effective_path = None
                   
                    # Path A: The one we requested
                    if mhcii_out.exists() and mhcii_out.is_file():
                        mhcii_effective_path = str(mhcii_out)
                        print_status(f"MHC-II output found at requested path: {mhcii_effective_path}", "success")
                   
                    # Path B: Check current working directory for default name
                    elif Path("mhcii_out.csv").exists() and Path("mhcii_out.csv").is_file():
                        shutil.move("mhcii_out.csv", str(mhcii_out))
                        mhcii_effective_path = str(mhcii_out)
                        print_status("Recovered MHC-II output from current directory → moved to expected path", "success")
                   
                    # Path C: Check results_dir for default name
                    elif (results_dir / "mhcii_out.csv").exists() and (results_dir / "mhcii_out.csv").is_file():
                        shutil.move(str(results_dir / "mhcii_out.csv"), str(mhcii_out))
                        mhcii_effective_path = str(mhcii_out)
                        print_status("Recovered MHC-II output from results_dir → moved to expected path", "success")
                   
                    else:
                        print_status("MHC-II output could not be located — skipping MHC-II", "warning")
                        print("→ Tool probably saved 'mhcii_out.csv' somewhere else.")
                        print("→ Run: find . -name 'mhcii_out.csv' 2>/dev/null to locate it")
                except Exception as e:
                    print_status(f"Error during epitope prediction: {e}", "error")
                    logger.error(f"Prediction calls failed: {e}")
                    mhcii_effective_path = None
         
            # Handle B-cell directory bug (already working in your previous run)
            bcell_effective_path = str(bcell_out)
            if os.path.isdir(bcell_out):
                csv_files = list(Path(bcell_out).glob("*.csv"))
                if csv_files:
                    bcell_effective_path = str(csv_files[0])
                    print_status(f"Found B-cell CSV inside directory: {bcell_effective_path}", "success")
                else:
                    bcell_effective_path = None
                    print_status("No .csv found inside B-cell directory — skipping B-cell", "warning")
            elif not os.path.isfile(bcell_out):
                bcell_effective_path = None
                print_status("B-cell output missing or not a file — skipping B-cell", "warning")
             
            # Extract positions
            epitope_positions = {}
            check_list = [
                ('bcell', bcell_effective_path),
                ('mhci', str(mhci_out)),
                ('mhcii', mhcii_effective_path)
            ]
           
            for typ, out_path_str in check_list:
                if out_path_str and os.path.exists(out_path_str) and os.path.isfile(out_path_str):
                    positions = extract_positions_from_prediction_csv(out_path_str)
                    if positions:
                        epitope_positions[typ] = positions
                        print_status(f"Success: {len(positions)} {typ.upper()} positions extracted", "success")
                    else:
                        print_status(f"No valid positions in {typ.upper()} CSV", "warning")
                else:
                    print_status(f"{typ.upper()} results not found, skipping", "warning")
        
            # Run SASA
            if epitope_positions and 'pdb_file' in res and Path(res['pdb_file']).exists():
                print("\n→ Running SASA exposure analysis on the 3D model...")
                pdb_path = res['pdb_file']
             
                for etype, pos_list in epitope_positions.items():
                    if not pos_list:
                        print_status(f"Skipping SASA for {etype.upper()}: no positions", "warning")
                        continue
                    temp_csv = results_dir / f"sasa_{etype}_input.csv"
                    pd.DataFrame(pos_list, columns=['Start', 'End']).to_csv(temp_csv, index=False)
                 
                    sasa_out = results_dir / f"{etype}_sasa_results.csv"
                 
                    result = run_sasa_analysis(
                        pdb_path=str(pdb_path),
                        epitope_csv=str(temp_csv),
                        output_csv=str(sasa_out),
                        sasa_threshold=30.0,
                        fraction_threshold=0.3
                    )
                 
                    if result["status"] == "success":
                        print_status(f" ✓ {etype.upper()}: {result['exposed_count']}/{result['total_count']} epitopes exposed", "success")
                    else:
                        print_status(f" ✗ SASA failed for {etype.upper()}: {result.get('message', 'unknown error')}", "error")
            elif do_sasa == 'y':
                print_status("Cannot run SASA: missing positions or PDB file", "warning")
 
    else:
        # Mode 1: Assembly logic — NOW WITH FULL ORDER DESCRIPTIONS & ALL 6 ORDERS
        print("\nPlease provide paths to the selected epitope CSV files.")
        if _assembly_config is not None:
            b_path = _assembly_config.bcell_csv_path
            c_path = _assembly_config.ctl_csv_path
            h_path = _assembly_config.htl_csv_path
        else:
            b_path = input("B-cell Epitopes CSV Path: ").strip()
            c_path = input("CTL (MHC-I) Epitopes CSV Path: ").strip()
            h_path = input("HTL (MHC-II) Epitopes CSV Path: ").strip()
   
        def load_csv(p):
            if p and os.path.exists(p): return pd.read_csv(p)
            return pd.DataFrame()
   
        dfs = [load_csv(b_path), load_csv(c_path), load_csv(h_path)]
   
        print("\n[Configuration]")
        n_term = ""
        if input("Add L7/L12 Adjuvant? (y/n): ").lower() == 'y':
            n_term = "MAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGAAVEAAEEQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAKAKLEAAGATVTVKEAAAK"
   
        c_term = "HHHHHH" if input("Add 6x-Histidine purification tag? (y/n): ").lower() == 'y' else ""
   
        try:
            if _assembly_config is not None:
                b_num = _assembly_config.bcell_count
                c_num = _assembly_config.ctl_count
                h_num = _assembly_config.htl_count
                order_opt = _assembly_config.assembly_order or "1"
            else:
                b_num = int(input("No. of B-cell epitopes: ") or 0)
                c_num = int(input("No. of CTL epitopes: ") or 0)
                h_num = int(input("No. of HTL epitopes: ") or 0)
              
                # ───────────────────────────────────────────────────────────────
                # Show user all order options BEFORE they choose
                # ───────────────────────────────────────────────────────────────
                print("\nAvailable Assembly Orders (choose 1-6):")
                print("1: B-cell → MHC-I → MHC-II")
                print("2: MHC-I → MHC-II → B-cell")
                print("3: MHC-II → B-cell → MHC-I")
                print("4: B-cell → MHC-II → MHC-I")
                print("5: MHC-I → B-cell → MHC-II")
                print("6: MHC-II → MHC-I → B-cell")
                print("Default: 1")
                order_opt = input("Assembly Order (1-6): ").strip() or "1"
     
            print_status("Assembling vaccine candidates...", "info")
            final_df = run_assembly(dfs, [b_num, c_num, h_num], ["KK", "AAY", "GPGPG"], order_opt, n_term, c_term)
    
            # ───────────────────────────────────────────────────────────────
            # Debug output
            # ───────────────────────────────────────────────────────────────
            print_status("Debug: Columns in final_df after assembly:", "info")
            print(final_df.columns.tolist())
         
            if final_df.empty:
                print_status("Error: No vaccine candidates assembled (empty DataFrame)", "error")
                return ["Vaccine_Assembly_Empty"]
        
            # ───────────────────────────────────────────────────────────────
            # MANUAL SEQUENCE BUILDING FOR ALL 6 ORDERS
            # ───────────────────────────────────────────────────────────────
            candidate = final_df.iloc[0]
        
            # Extract epitope lists
            b_epitopes = candidate.get('B_cell_Epitopes', '')
            ctl_epitopes = candidate.get('CTL_Epitopes', '')
            htl_epitopes = candidate.get('HTL_Epitopes', '')
        
            # Convert to lists if comma-separated strings
            if isinstance(b_epitopes, str):
                b_epitopes = [e.strip() for e in b_epitopes.split(',') if e.strip()]
            if isinstance(ctl_epitopes, str):
                ctl_epitopes = [e.strip() for e in ctl_epitopes.split(',') if e.strip()]
            if isinstance(htl_epitopes, str):
                htl_epitopes = [e.strip() for e in htl_epitopes.split(',') if e.strip()]
        
            # Linkers per type
            linkers = {"bcell": "KK", "mhci": "AAY", "mhcii": "GPGPG"}
        
            # Define all 6 orders with descriptions
            order_patterns = {
                1: (["bcell", "mhci", "mhcii"], "B-cell → MHC-I → MHC-II"),
                2: (["mhci", "mhcii", "bcell"], "MHC-I → MHC-II → B-cell"),
                3: (["mhcii", "bcell", "mhci"], "MHC-II → B-cell → MHC-I"),
                4: (["bcell", "mhcii", "mhci"], "B-cell → MHC-II → MHC-I"),
                5: (["mhci", "bcell", "mhcii"], "MHC-I → B-cell → MHC-II"),
                6: (["mhcii", "mhci", "bcell"], "MHC-II → MHC-I → B-cell"),
            }
        
            order_num = int(order_opt)
            pattern, order_desc = order_patterns.get(order_num, (["bcell", "mhci", "mhcii"], "Default: B-cell → MHC-I → MHC-II"))
        
            print_status(f"Building vaccine sequence using order {order_num}: {order_desc}", "info")
        
            sequence_parts = [n_term] if n_term else []
        
            for typ in pattern:
                epitopes = b_epitopes if typ == "bcell" else ctl_epitopes if typ == "mhci" else htl_epitopes
                linker = linkers[typ]
                sequence_parts += [pep + linker for pep in epitopes]
        
            if c_term:
                sequence_parts.append(c_term)
        
            final_vaccine_sequence = "".join(sequence_parts)
        
            print_status(f"Assembled vaccine sequence (order {order_num}): length {len(final_vaccine_sequence)} aa", "success")
            print(f"Sequence preview: {final_vaccine_sequence[:100]}...")
            if len(final_vaccine_sequence) > 100:
                print(f"... (total {len(final_vaccine_sequence)} aa)")
        
            output_file = results_dir / f"Assembled_Vaccines_{TIMESTAMP}.csv"
            final_df.to_csv(output_file, index=False)
            print_status(f"[✓] Vaccine candidates generated: {output_file}", "success")
      
            graph_file = results_dir / f"Vaccine_Architecture_{TIMESTAMP}.png"
            order_list = assembly_orders.get(order_num, assembly_orders[1])
            plot_vaccine_architecture(n_term, [], order_list, ["KK", "AAY", "GPGPG"], c_term, graph_file)
    
        except Exception as e:
            logger.error(f"Assembly failed: {e}", exc_info=True)
            print_status(f"Error during assembly: {e}", "error")
            return ["Vaccine_Assembly"]
    return []
# 3. UPDATED MAIN FUNCTION
def main():
    logger = setup_logging(RESULTS_DIR)
    check_dependencies()
    
    while True:
        # 1. ASK STRATEGY FIRST
        print("\n" + "="*30)
        print(" VAXELAN MAIN MENU")
        print("="*30)
        print("1 - Epitope Prediction")
        print("2 - Protein Prioritization")
        print("3 - Virulence, Toxicity, and Host Compatibility Assessment")
        print("4 - Comprehensive Multi-Tool Analysis")
        print("5 - Predicted Epitope Analysis")
        print("6 - Multi-Epitope Vaccine Assembly & 3D Modeling")
        
        strategy = input("\nSelect strategy (1-6) or 'q' to quit: ").strip().lower()
        
        if strategy == 'q':
            print_status("Exiting program. Have a good day!", "success")
            break
            
        if strategy not in ["1", "2", "3", "4", "5", "6"]:
            print_status("Invalid strategy. Please choose 1, 2, 3, 4, 5 or 6.", "error")
            continue

        results_dir = RESULTS_DIR / f"strategy_{strategy}"
        results_dir.mkdir(parents=True, exist_ok=True)
        failed_tools = []

        # 2. CONDITIONAL INPUT FLOW
        if strategy == "6":
            # Strategy 6 handles its own sequence input (custom or CSV)
            failed_tools = strategy_6(results_dir, logger)
        
        else:
            # STRATEGIES 1-5: Need Sequence and Pathogen Type
            input_value = input("\nEnter the input (FASTA file path or UniProt ID): ").strip()
            is_uniprot = not os.path.isfile(input_value)
            temp_dir = None
            
            try:
                # Handle Sequence Input
                if is_uniprot:
                    fasta_content = fetch_uniprot_sequence(input_value, logger)
                    temp_dir = Path(tempfile.mkdtemp(prefix="vaxelan_"))
                    input_file = temp_dir / f"{input_value}.fasta"
                    with open(input_file, "w") as f:
                        f.write(fasta_content)
                else:
                    input_file = Path(input_value)
                    if not input_file.is_file() or not validate_fasta_file(input_file):
                        print_status("Error: Invalid FASTA file or Path. Please try again.", "error")
                        continue
                    temp_dir = Path(tempfile.mkdtemp(prefix="vaxelan_"))

                # Split Batches
                batch_files, batch_dirs, uniprot_mapping = split_fasta_into_batches(input_file, BATCH_SIZE, temp_dir)
                if not batch_files:
                    print_status("Error: No valid sequences found.", "error")
                    continue

                # Ask for Pathogen Type
                selected_pathogen_type = select_pathogen_type()
                confirmed_pathogen_type = validate_pathogen_type(input_value, selected_pathogen_type, is_uniprot)

                # Execute selected Strategy 1-5
                if strategy == "1":
                    print_status("Executing Strategy 1: Epitope Prediction", "info")
                    failed_tools = strategy_1(input_file, results_dir, confirmed_pathogen_type, logger, batch_files, batch_dirs, uniprot_mapping)
                elif strategy == "2":
                    print_status("Executing Strategy 2: Protein Prioritization", "info")
                    failed_tools = strategy_2(input_file, results_dir, input_value, is_uniprot, logger, batch_files, batch_dirs, uniprot_mapping)
                elif strategy == "3":
                    print_status("Executing Strategy 3: Virulence Assessment", "info")
                    failed_tools = strategy_3(input_file, confirmed_pathogen_type, results_dir, logger, batch_files, batch_dirs, uniprot_mapping)
                elif strategy == "4":
                    print_status("Executing Strategy 4: Comprehensive Analysis", "info")
                    failed_tools = strategy_4(input_file, confirmed_pathogen_type, results_dir, logger, batch_files, batch_dirs, uniprot_mapping)
                elif strategy == "5":
                    print_status("Executing Strategy 5: Predicted Epitope Analysis", "info")
                    failed_tools = strategy_5(input_file, confirmed_pathogen_type, results_dir, logger, batch_files, batch_dirs, uniprot_mapping)

            except ValueError as e:
                print_status(f"Error: {e}", "error")
                continue
            except Exception as e:
                print_status(f"Unexpected error: {e}", "error")
                logger.error(f"Unexpected error in main loop: {e}", exc_info=True)
                continue
            finally:
                if temp_dir and temp_dir.exists():
                    shutil.rmtree(temp_dir, ignore_errors=True)
                    logger.info(f"Cleaned up temporary directory: {temp_dir}")

        # 3. FINAL STATUS & LOOP BACK
        if failed_tools:
            logger.warning(f"Failed tools: {', '.join(failed_tools)}")
            print_status(f"Warning: The following tools failed: {', '.join(failed_tools)}", "warning")

        run_again = input("\nAction complete. Back to Main Menu? (y/n): ").strip().lower()
        if run_again != 'y':
            print_status("Have a good day!", "success")
            break

if __name__ == "__main__":
    main()
