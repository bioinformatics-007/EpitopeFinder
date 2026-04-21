#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
import argparse
from io import StringIO
from pathlib import Path
import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import psutil
import shutil

# Set up logging
logger = logging.getLogger('VaxElan')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Dynamically resolve root directory up to Vaxelan_2_0
current_path = Path(__file__).resolve()
while current_path.name != "Vaxelan_2_0":
    current_path = current_path.parent
    if current_path == current_path.parent:
        logger.error("Could not locate 'Vaxelan_2_0' in path hierarchy.")
        raise RuntimeError("Could not locate 'Vaxelan_2_0' in path hierarchy.")
root_dir = current_path

# MHC-II tool path
mhcii_tool = os.path.join(root_dir, "tools", "IEDB_MHC_II-3.1.12", "mhc_ii", "mhc_II_binding.py")

# Complete allele list (optimized for consensus3 and netmhciipan methods)
ALLELES_LIST = [
    # HLA-DRB alleles (commonly supported by consensus3)
    'HLA-DRB1*01:01', 'HLA-DRB1*01:02', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01',
    'HLA-DRB1*04:04', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:01',
    'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01', 'HLA-DRB1*11:04',
    'HLA-DRB1*12:01', 'HLA-DRB1*13:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01',
    'HLA-DRB1*15:02', 'HLA-DRB3*01:01', 'HLA-DRB3*02:01', 'HLA-DRB4*01:01',
    'HLA-DRB5*01:01',
    # DQA1/DQB1 alleles (included for netmhciipan methods)
    'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01',
    'HLA-DQA1*03:01/DQB1*03:02', 'HLA-DQA1*04:01/DQB1*04:02',
    'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02',
    # DPA1/DPB1 alleles (included for netmhciipan methods)
    'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DPA1*01:03/DPB1*02:01',
    'HLA-DPA1*03:01/DPB1*04:02', 'HLA-DPA1*02:01/DPB1*05:01',
]
UNSUPPORTED_ALLELES = ['HLA-DRB3*02:02', 'HLA-DPA1*01:03/DPB1*04:01', 'HLA-DPA1*02:01/DPB1*14:01']

def find_python_interpreter():
    """Find the current Python 3 interpreter."""
    python_cmd = shutil.which("python3")
    if not python_cmd:
        logger.error("No Python 3 interpreter found. Ensure Python 3 is installed.")
        return None
    try:
        result = subprocess.run(
            [python_cmd, "--version"],
            capture_output=True, text=True, check=True
        )
        version_output = result.stderr or result.stdout
        if "Python 3" in version_output:
            logger.info(f"Using Python 3 interpreter: {python_cmd}")
            return python_cmd
        else:
            logger.error(f"Interpreter {python_cmd} is not Python 3: {version_output.strip()}")
            return None
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Error checking interpreter {python_cmd}: {e}")
        return None

def print_status(msg, status="info"):
    """Print colored status messages and log them."""
    colors = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")
    if status == "info":
        logger.info(msg)
    elif status == "warning":
        logger.warning(msg)
    elif status == "error":
        logger.error(msg)
    else:
        logger.debug(msg)

def log_memory_usage():
    """Log current memory usage in MB."""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / 1024 / 1024

def validate_fasta(fasta_file):
    """Validate FASTA file contains amino acid sequences."""
    try:
        with open(fasta_file, 'r') as f:
            sequence = ""
            header_count = 0
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header_count += 1
                elif line:
                    sequence += line
            if not sequence:
                print_status(f"Input FASTA {fasta_file} is empty.", "error")
                logger.error(f"Input FASTA {fasta_file} is empty.")
                return False
            if header_count == 0:
                print_status(f"Input FASTA {fasta_file} has no headers.", "error")
                logger.error(f"Input FASTA {fasta_file} has no headers.")
                return False
        # Check for DNA-like sequences
        if re.match(r'^[ACGT]+$', sequence.upper()):
            print_status(f"Input FASTA {fasta_file} contains DNA-like sequences. MHC-II predictions require amino acid sequences.", "error")
            logger.error(f"Input FASTA {fasta_file} contains DNA-like sequences.")
            return False
        # Basic amino acid check
        if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence.upper()):
            print_status(f"Input FASTA {fasta_file} contains invalid amino acid characters.", "error")
            logger.error(f"Input FASTA {fasta_file} contains invalid amino acid characters.")
            return False
        return True
    except Exception as e:
        print_status(f"Error validating FASTA file {fasta_file}: {e}", "error")
        logger.error(f"Error validating FASTA file {fasta_file}: {e}")
        return False

def get_valid_alleles(method):
    """Get list of valid alleles for a given method."""
    python_cmd = find_python_interpreter()
    if not python_cmd:
        return []
    try:
        result = subprocess.run(
            [python_cmd, mhcii_tool, "allele"],
            capture_output=True,
            text=True,
            check=True,
        )
        logger.debug(f"Raw allele output for method {method}:\n{result.stdout}")
        alleles = []
        lines = result.stdout.splitlines()
        for line in lines:
            line = line.strip()
            if not line or line.startswith(("Available", "----", "Method")):
                continue
            # Extract allele names, handling various formats
            allele = line.split('\t')[0].strip()
            if allele and re.match(r'^HLA-[A-Z0-9*:/]+$', allele):
                alleles.append(allele)
        if not alleles:
            logger.warning(f"No alleles parsed for method {method}. Falling back to ALLELES_LIST.")
            print_status(f"Warning: No alleles parsed for {method}. Using predefined ALLELES_LIST.", "warning")
            return ALLELES_LIST  # Fallback to predefined list
        return alleles
    except subprocess.CalledProcessError as e:
        print_status(f"Error fetching valid alleles for {method}: {e.stderr}", "error")
        logger.error(f"Error fetching valid alleles for {method}: {e.stderr}")
        return []

def suggest_alternative_methods():
    """Suggest alternative methods with supported alleles."""
    python_cmd = find_python_interpreter()
    if not python_cmd:
        return "No Python interpreter found."
    try:
        result = subprocess.run(
            [python_cmd, mhcii_tool, "method"],
            capture_output=True,
            text=True,
            check=True,
        )
        methods = [line.strip() for line in result.stdout.splitlines() if line.strip() and not line.startswith(("Available", "----"))]
        suggestion = "Try alternative methods:\n"
        for method in methods:
            alleles = get_valid_alleles(method)
            if alleles:
                suggestion += f"- {method}: Supports alleles like {', '.join(alleles[:3])}{', ...' if len(alleles) > 3 else ''}\n"
            else:
                suggestion += f"- {method}: No supported alleles found\n"
        return suggestion
    except subprocess.CalledProcessError as e:
        print_status(f"Error fetching alternative methods: {e.stderr}", "error")
        logger.error(f"Error fetching alternative methods: {e.stderr}")
        return "Unable to fetch alternative methods."

def run_mhc2_prediction(method, fasta_file, allele, output_dir, timeout=3600):
    """Run MHC-II prediction for a single allele."""
    logger.info(f"Running prediction for allele: {allele}, method: {method}")
    
    # Validate tool and input
    if not os.path.isfile(mhcii_tool):
        print_status(f"MHC-II tool '{mhcii_tool}' not found.", "error")
        logger.error(f"MHC-II tool '{mhcii_tool}' not found.")
        return None
    if not os.access(mhcii_tool, os.X_OK):
        print_status(f"MHC-II tool '{mhcii_tool}' is not executable. Run 'chmod +x {mhcii_tool}'.", "error")
        logger.error(f"MHC-II tool '{mhcii_tool}' is not executable.")
        return None
    if not os.path.isfile(fasta_file):
        print_status(f"Input file '{fasta_file}' not found.", "error")
        logger.error(f"Input file '{fasta_file}' not found.")
        return None
    if not os.access(fasta_file, os.R_OK):
        print_status(f"Input file '{fasta_file}' is not readable.", "error")
        logger.error(f"Input file '{fasta_file}' is not readable.")
        return None

    # Find Python interpreter
    python_cmd = find_python_interpreter()
    if not python_cmd:
        print_status("No Python 3 interpreter found.", "error")
        logger.error("No Python 3 interpreter found.")
        return None

    # Ensure debug directory exists
    debug_dir = os.path.join(output_dir, "mhcii_debug")
    os.makedirs(debug_dir, exist_ok=True)
    logger.info(f"Created MHC-II debug directory: {debug_dir}")

    # Build command
    command = [python_cmd, str(mhcii_tool), method, allele, str(fasta_file)]
    print_status(f"→ Running prediction: {method} with allele: {allele}", "info")
    logger.info(f"Executing command: {' '.join(command)}")

    # Save debug output
    safe_allele = allele.replace('*', '_').replace('/', '_').replace(':', '_')
    debug_output = os.path.join(debug_dir, f"mhc2_{safe_allele}.txt")

    try:
        mem_usage = log_memory_usage()
        total_mem = psutil.virtual_memory().total / 1024 / 1024
        if mem_usage > 0.8 * total_mem:
            print_status(f"Warning: High memory usage ({mem_usage:.2f}/{total_mem:.2f} MB).", "warning")
            logger.warning(f"High memory usage: {mem_usage:.2f}/{total_mem:.2f} MB")

        start = time.time()
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout,
            stdin=subprocess.DEVNULL,
        )
        elapsed = time.time() - start
        logger.info(f"Prediction for {allele} took {elapsed:.2f} seconds")
        
        if "Error" in result.stdout or "Exception" in result.stdout:
            print_status(f"Error in prediction output for {allele}:\n{result.stdout}", "error")
            logger.error(f"Error in prediction output for {allele}: {result.stdout}")
            with open(debug_output + ".error", "w") as f:
                f.write(result.stdout)
            return None
        if result.stderr:
            if "Selected prediction method" in result.stderr and "does not exist" in result.stderr:
                print_status(f"Invalid method '{method}' for {allele}. Run 'python3 {mhcii_tool} method' to see valid methods.", "error")
                logger.error(f"Invalid method '{method}' for {allele}: {result.stderr}")
                try:
                    valid_methods = subprocess.run(
                        [python_cmd, mhcii_tool, "method"],
                        capture_output=True,
                        text=True,
                        check=True,
                    ).stdout
                    print_status(f"Available methods:\n{valid_methods}", "info")
                except subprocess.CalledProcessError as e:
                    print_status(f"Could not fetch valid methods: {e.stderr}", "error")
            else:
                print_status(f"Prediction stderr for {allele}: {result.stderr}", "warning")
                logger.warning(f"Prediction stderr for {allele}: {result.stderr}")
            with open(debug_output + ".error", "w") as f:
                f.write(result.stderr)
            return None
        with open(debug_output, "w") as f:
            f.write(result.stdout)
        return result.stdout
    except subprocess.TimeoutExpired as e:
        print_status(f"Timeout after {timeout}s for {allele}: {e}", "error")
        logger.error(f"Timeout after {timeout}s for {allele}: {e}")
        with open(debug_output + ".error", "w") as f:
            f.write(str(e))
        return None
    except subprocess.CalledProcessError as e:
        if "Selected prediction method" in e.stderr and "does not exist" in e.stderr:
            print_status(f"Invalid method '{method}' for {allele}. Run 'python3 {mhcii_tool} method' to see valid methods.", "error")
            logger.error(f"Invalid method '{method}' for {allele}: {e.stderr}")
            try:
                valid_methods = subprocess.run(
                    [python_cmd, mhcii_tool, "method"],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout
                print_status(f"Available methods:\n{valid_methods}", "info")
            except subprocess.CalledProcessError as e:
                print_status(f"Could not fetch valid methods: {e.stderr}", "error")
        else:
            print_status(f"Failed to run prediction for {allele}: {e.stderr}", "error")
            logger.error(f"Failed to run prediction for {allele}: {e.stderr}")
        with open(debug_output + ".error", "w") as f:
            f.write(e.stderr)
        return None
    except Exception as e:
        print_status(f"Unexpected error for {allele}: {e}", "error")
        logger.error(f"Unexpected error for {allele}: {e}")
        with open(debug_output + ".error", "w") as f:
            f.write(str(e))
        return None

def parse_results(response_text, method, allele, score_threshold):
    """Parse prediction results into a DataFrame."""
    if not response_text or not response_text.strip():
        print_status(f"No results for {allele}", "warning")
        logger.warning(f"No results for {allele}")
        return pd.DataFrame()
    try:
        logger.debug(f"Raw response for {allele}:\n{response_text[:500]}...")  # Log first 500 chars
        print_status(f"Parsing results for {allele}", "info")
        data_io = StringIO(response_text)
        df = pd.read_csv(data_io, sep=r'\s*,\s*|\t|\s+', engine="python", skipinitialspace=True)
        # Normalize column names: strip whitespace, convert to lowercase for matching
        df.columns = [col.strip().lower() for col in df.columns]
        print_status(f"Columns found: {df.columns.tolist()}", "info")
        logger.debug(f"Columns found for {allele}: {df.columns.tolist()}")

        score_col = None
        # Prioritize consensus_percentile_rank for consensus3 method
        if method == "consensus3":
            for col in ['consensus_percentile_rank', 'smm_align_ic50', 'nn_align_ic50']:
                if col in df.columns:
                    score_col = col
                    break
        # Fallback for other methods
        if not score_col:
            for col in ['score', 'ic50', 'affinity', 'smm_align_ic50', 'nn_align_ic50']:
                if col in df.columns:
                    score_col = col
                    break

        if score_col:
            logger.info(f"Selected score column for {allele}: {score_col}")
            print_status(f"Using score column: {score_col}", "info")
            df['score'] = pd.to_numeric(df[score_col], errors='coerce')
            # For percentile ranks (lower is better), use <= threshold
            if 'percentile_rank' in score_col:
                df = df[df['score'] <= score_threshold]
            # For IC50 (lower is better), use <= threshold
            elif 'ic50' in score_col:
                df = df[df['score'] <= score_threshold]
            # For other scores (e.g., Affinity), assume lower is better
            else:
                df = df[df['score'] <= score_threshold]
            if df.empty:
                print_status(f"No peptides with score <= {score_threshold} for {allele} in column {score_col}", "warning")
                logger.warning(f"No peptides with score <= {score_threshold} for {allele} in column {score_col}")
        else:
            print_status(f"No score column found for {allele}. Keeping all results.", "warning")
            logger.warning(f"No score column found for {allele}. Keeping all results.")
            # Keep all results if no score column is found
            df['score'] = pd.NA

        df.insert(0, "Allele", allele)
        df.insert(0, "Method", method)
        return df
    except Exception as e:
        print_status(f"Error parsing results for {allele}: {e}", "error")
        logger.error(f"Error parsing results for {allele}: {e}")
        return pd.DataFrame()

def run_mhc2(method_code, fasta_file, output_dir=".", output_file="mhcii_out.csv", score_threshold=50):
    """Run MHC-II predictions for all alleles and save to a single CSV."""
    logger.info(f"Starting run_mhc2 with method_code: {method_code}, input_file: {fasta_file}, output_file: {output_file}")
    start_time = time.time()
    methods_ls = {
        "nmel": "netmhciipan_el",
        "nmba": "netmhciipan_ba",
        "cons": "consensus3",
        "nn": "nn_align",
        "smm": "smm_align",
        "comblib": "comblib",
        "sturniolo": "sturniolo",
        "nmel42": "netmhciipan_el-4.2",
        "nmba42": "netmhciipan_ba-4.2",
        "nmel43": "netmhciipan_el-4.3",
        "nmba43": "netmhciipan_ba-4.3",
    }
    if method_code not in methods_ls:
        valid_methods = ', '.join(methods_ls.keys())
        print_status(f"Invalid method code: {method_code}. Choose from: {valid_methods}", "error")
        logger.error(f"Invalid method code: {method_code}. Valid options: {valid_methods}")
        return 1
    method = methods_ls[method_code]
    logger.info(f"Resolved method: {method}")

    # Validate FASTA file
    if not validate_fasta(fasta_file):
        return 1

    # Validate output path
    output_path = Path(output_file)
    output_file_str = str(output_file)  # Convert to string for string operations
    if output_path.is_dir():
        print_status(f"Output path {output_file} is a directory, not a file.", "error")
        logger.error(f"Output path {output_file} is a directory, not a file.")
        return 1
    if not output_file_str.lower().endswith('.csv'):
        print_status(f"Output file {output_file} should have a .csv extension.", "error")
        logger.error(f"Output file {output_file} should have a .csv extension.")
        return 1
    output_dir = str(output_path.parent)
    if not os.access(output_dir, os.W_OK):
        print_status(f"Output directory {output_dir} is not writable.", "error")
        logger.error(f"Output directory {output_dir} is not writable.")
        return 1
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Ensuring output directory exists: {output_dir}")

    # Validate alleles for the method
    valid_alleles = get_valid_alleles(method)
    filtered_alleles = [allele for allele in ALLELES_LIST if allele in valid_alleles]
    if not filtered_alleles:
        print_status(f"No valid alleles for method {method}. Check supported alleles with: python3 {mhcii_tool} allele", "warning")
        logger.warning(f"No valid alleles for method {method}. Using predefined ALLELES_LIST.")
        filtered_alleles = [allele for allele in ALLELES_LIST if allele not in UNSUPPORTED_ALLELES]
        if not filtered_alleles:
            print_status("No alleles available after excluding unsupported alleles.", "error")
            logger.error("No alleles available after excluding unsupported alleles.")
            suggestion = suggest_alternative_methods()
            print_status(suggestion, "info")
            return 1

    if len(filtered_alleles) < len(ALLELES_LIST):
        unsupported = set(ALLELES_LIST) - set(filtered_alleles)
        print_status(f"Unsupported alleles for {method}: {', '.join(unsupported)}", "warning")
        logger.warning(f"Unsupported alleles for {method}: {', '.join(unsupported)}")

    # Log memory usage at start
    mem_usage = log_memory_usage()
    logger.info(f"Initial memory usage: {mem_usage:.2f} MB")

    all_results = []
    failed_alleles = []
    print_status("\n**Starting MHC Class II Prediction for ALL Alleles...**", "info")
    logger.info(f"Skipping unsupported alleles: {', '.join(UNSUPPORTED_ALLELES)}")

    max_workers = min(len(filtered_alleles), max(1, os.cpu_count() // 2))
    logger.info(f"Using {max_workers} workers for parallel predictions")

    def process_allele(allele):
        """Process a single allele's prediction."""
        try:
            logger.info(f"Processing allele: {allele}")
            response_text = run_mhc2_prediction(method, fasta_file, allele, output_dir)
            if response_text:
                df = parse_results(response_text, method, allele, score_threshold)
                return df, None
            return pd.DataFrame(), allele
        except Exception as e:
            logger.error(f"Error in process_allele for {allele}: {e}")
            return pd.DataFrame(), allele

    try:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_allele = {executor.submit(process_allele, allele): allele for allele in filtered_alleles}
            for future in as_completed(future_to_allele):
                try:
                    df, failed = future.result()
                    if not df.empty:
                        all_results.append(df)
                    if failed:
                        failed_alleles.append(failed)
                except Exception as e:
                    allele = future_to_allele[future]
                    print_status(f"Error processing {allele}: {e}", "error")
                    logger.error(f"Error processing {allele}: {e}")
                    failed_alleles.append(allele)
    except Exception as e:
        print_status(f"Error in parallel execution: {e}", "error")
        logger.error(f"Error in parallel execution: {e}")
        return 1

    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        try:
            final_df.to_csv(output_path, index=False, mode='w')
            elapsed_time = time.time() - start_time
            print_status(f"\n✔ Predictions complete in {elapsed_time:.2f} seconds. Results saved to: {output_path}", "success")
            logger.info(f"Predictions complete in {elapsed_time:.2f} seconds. Results saved to: {output_path}")
        except Exception as e:
            print_status(f"Error saving results to {output_path}: {e}", "error")
            logger.error(f"Error saving results to {output_path}: {e}")
            return 1
        if failed_alleles:
            print_status(f"Warning: Predictions failed for: {', '.join(failed_alleles)}", "warning")
            logger.warning(f"Predictions failed for: {', '.join(failed_alleles)}")
        return 0
    else:
        elapsed_time = time.time() - start_time
        print_status(f"No valid predictions found in {elapsed_time:.2f} seconds. No output saved.", "warning")
        logger.warning(f"No valid predictions found in {elapsed_time:.2f} seconds. No output saved.")
        return 2

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MHC-II prediction for all alleles.")
    parser.add_argument("-m", "--method", required=True, help="Method code (e.g., cons, nmel, nmba, nmel42, nmba42, nmel43, nmba43, nn, smm, comblib, sturniolo).")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file with amino acid sequences.")
    parser.add_argument("-o", "--output", default="mhcii_out.csv", help="Output CSV file for results.")
    parser.add_argument("-d", "--output_dir", default=".", help="Output directory for debug files.")
    parser.add_argument("-s", "--score_threshold", type=float, default=50, help="Score threshold for filtering results (default: 50).")
    args = parser.parse_args()
    logger.info(f"Running mhc_ii.py with args: {vars(args)}")
    status = run_mhc2(
        method_code=args.method,
        fasta_file=args.input,
        output_dir=args.output_dir,
        output_file=args.output,
        score_threshold=args.score_threshold
    )
    exit(status)
