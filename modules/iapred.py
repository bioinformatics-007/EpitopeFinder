
import os
import subprocess
import shutil
from pathlib import Path
import tempfile
import pandas as pd
import argparse
import logging
from itertools import cycle
import uuid
import warnings

# Suppress scikit-learn version mismatch warnings
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn.base")

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Dynamically set IAPRED_TOOL_PATH relative to this script's parent directory
BASE_DIR = Path(__file__).resolve().parent.parent
IAPRED_TOOL_PATH = BASE_DIR / "tools" / "IApred"

def print_status(msg, status="info"):
    colors = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")
    logging.log(logging.INFO if status != "error" else logging.ERROR, msg)

def parse_fasta(input_fasta):
    header = None
    sequence = []
    sequences_found = False
    sequence_count = 0
    try:
        with open(input_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if header and sequence:
                        sequences_found = True
                        full_header = header
                        base_header = header.split('|')[1] if '|' in header else header
                        sequence_count += 1
                        yield full_header, base_header, ''.join(sequence)
                    header = line[1:].split()[0]
                    sequence = []
                else:
                    sequence.append(line)
            if header and sequence:
                sequences_found = True
                full_header = header
                base_header = header.split('|')[1] if '|' in header else header
                sequence_count += 1
                yield full_header, base_header, ''.join(sequence)
        if not sequences_found:
            raise ValueError(f"No valid sequences found in FASTA file: {input_fasta}")
        logging.debug(f"Parsed {sequence_count} sequences from FASTA")
    except Exception as e:
        logging.error(f"Error parsing FASTA file {input_fasta}: {str(e)}")
        raise

def run_iapred(input_fasta, output_dir, batch_idx=1, output_file=None, uniprot_id=None):
    """
    Runs IAPred to predict immunogenic antigens, producing individual and combined CSV outputs.

    Args:
        input_fasta (str or Path): Path to the input FASTA file.
        output_dir (str or Path): Directory where output should be saved.
        batch_idx (int): Batch index for naming combined output (default: 1).
        output_file (str, optional): Path to save output for a single sequence (default: None).
        uniprot_id (str, optional): UniProt ID for specific sequence output (default: None).

    Returns:
        int: Exit status code (0 for success, 1 for error).
    """
    input_fasta = str(input_fasta)
    # Use the provided output_dir directly without appending 'iapred'
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not os.access(output_dir, os.W_OK):
        logging.error(f"Output directory is not writable: {output_dir}")
        print_status(f"Output directory is not writable: {output_dir}", "error")
        return 1

    combined_output_path = output_dir / f"combined_iapred_batch_{batch_idx}.csv"
    iapred_script = IAPRED_TOOL_PATH / "IApred.py"

    if not iapred_script.is_file():
        logging.error(f"IAPred script not found at {iapred_script}")
        print_status(f"IAPred script not found at {iapred_script}", "error")
        return 1

    models_path = IAPRED_TOOL_PATH / "models"
    if not models_path.is_dir():
        logging.error(f"'models' directory not found in {models_path}")
        print_status(f"'models' directory not found in {models_path}", "error")
        return 1

    # Copy models directory to current working directory
    try:
        shutil.copytree(models_path, "models", dirs_exist_ok=True)
        logging.debug(f"Copied models directory to {os.getcwd()}/models")
    except Exception as e:
        logging.error(f"Failed to copy models directory: {str(e)}")
        print_status(f"Failed to copy models directory: {str(e)}", "error")
        return 1

    # Use unique temporary file names to avoid conflicts
    batch_fasta = os.path.join(tempfile.gettempdir(), f"batch_input_{uuid.uuid4().hex}.fasta")
    batch_output = os.path.join(tempfile.gettempdir(), f"batch_iapred_out_{uuid.uuid4().hex}.csv")

    try:
        sequence_dict = {}
        base_to_full = {}
        sequences = []
        for full_header, base_header, seq in parse_fasta(input_fasta):
            sequence_dict[full_header] = seq
            base_to_full.setdefault(base_header, []).append(full_header)
            sequences.append(seq)
        with open(batch_fasta, 'w') as f:
            for full_header, seq in sequence_dict.items():
                f.write(f">{full_header}\n{seq}\n")
        logging.debug(f"Created batch FASTA file: {batch_fasta}")
    except Exception as e:
        logging.error(f"Failed to process FASTA file: {str(e)}")
        print_status(f"Failed to process FASTA file: {str(e)}", "error")
        return 1

    command = ["python3", str(iapred_script), batch_fasta, batch_output, "-v"]

    print_status("Running IAPred batch analysis...", "info")
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        logging.debug(f"IAPred stdout: {result.stdout}")
        if result.stderr:
            logging.warning(f"IAPred stderr: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logging.error(f"IAPred failed: {e.stderr}")
        print_status(f"IAPred failed: {e.stderr}", "error")
        return 1
    except Exception as e:
        logging.error(f"Unexpected error running IAPred: {str(e)}")
        print_status(f"Unexpected error running IAPred: {str(e)}", "error")
        return 1

    if not os.path.isfile(batch_output):
        logging.error(f"IAPred did not produce output file: {batch_output}")
        print_status(f"IAPred did not produce output file: {batch_output}", "error")
        return 1

    try:
        df = pd.read_csv(batch_output)
        logging.debug(f"IAPred output CSV columns: {df.columns.tolist()}")
        if df.empty:
            logging.warning("IAPred output CSV is empty")
            print_status("IAPred output CSV is empty", "warning")
    except Exception as e:
        logging.error(f"Failed to read IAPred output CSV: {str(e)}")
        print_status(f"Failed to read IAPred output CSV: {str(e)}", "error")
        return 1

    if 'Header' not in df.columns:
        logging.error("'Header' column missing in IAPred output.")
        print_status("'Header' column missing in IAPred output.", "error")
        return 1

    df['UniProt_ID'] = df['Header'].apply(lambda x: x.split('|')[1] if '|' in x else x)
    df['Header'] = df['Header'].str.rstrip('_')

    def map_sequence(row, sequence_iterator):
        header = row['Header']
        if header in sequence_dict:
            return sequence_dict[header]
        for base_header, full_headers in base_to_full.items():
            if header.lower() == base_header.lower() or header in base_header or base_header in header:
                return next(sequence_iterator)
        logging.warning(f"No sequence found for header: {header}")
        return "UNMAPPED"

    sequence_iterator = cycle(sequences)
    df["Epitope_Sequence"] = df.apply(lambda row: map_sequence(row, sequence_iterator), axis=1)

    if 'Start' in df.columns and 'End' in df.columns:
        def extract_epitope(row):
            if row['Epitope_Sequence'] == "UNMAPPED" or pd.isna(row['Start']) or pd.isna(row['End']):
                return row['Epitope_Sequence']
            start = int(row['Start']) - 1
            end = int(row['End'])
            try:
                return row['Epitope_Sequence'][start:end]
            except Exception as e:
                logging.warning(f"Failed to extract epitope for {row['Header']}: {str(e)}")
                return "ERROR"
        df['Epitope_Sequence'] = df.apply(extract_epitope, axis=1)
    else:
        logging.warning("Start and End columns missing; using full sequence as epitope.")

    try:
        df.to_csv(combined_output_path, index=False)
        print_status(f"Combined IAPred results saved to: {combined_output_path} with {len(df)} records", "success")
        logging.info(f"Saved combined output CSV: {combined_output_path}")
    except Exception as e:
        logging.error(f"Failed to save combined output CSV: {str(e)}")
        print_status(f"Failed to save combined output CSV: {str(e)}", "error")
        return 1

    try:
        unique_uniprot_ids = df['UniProt_ID'].unique()
        logging.debug(f"Unique UniProt IDs: {unique_uniprot_ids}")
        for uniprot_id in unique_uniprot_ids:
            indiv_df = df[df['UniProt_ID'] == uniprot_id]
            indiv_csv_path = Path(output_file) if output_file and uniprot_id in str(output_file) else output_dir / f"{uniprot_id}_iapred.csv"
            indiv_csv_path.parent.mkdir(parents=True, exist_ok=True)
            indiv_df.to_csv(indiv_csv_path, index=False)
            print_status(f"Individual CSV saved to: {indiv_csv_path} with {len(indiv_df)} records", "success")
            logging.info(f"Saved individual CSV: {indiv_csv_path}")
    except Exception as e:
        logging.error(f"Failed to save individual CSV files: {str(e)}")
        print_status(f"Failed to save individual CSV files: {str(e)}", "error")
        return 1

    try:
        os.remove(batch_fasta)
        logging.debug(f"Removed temporary FASTA file: {batch_fasta}")
    except FileNotFoundError:
        logging.debug(f"Temporary FASTA file not found for removal: {batch_fasta}")
    shutil.rmtree("models", ignore_errors=True)
    logging.debug("Removed temporary models directory")

    print_status("IAPred processing completed successfully.", "success")
    return 0

def main():
    parser = argparse.ArgumentParser(description="Run IAPred for immunogenic antigen prediction.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-d", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--batch_idx", type=int, default=1, help="Batch index")
    parser.add_argument("--output_file", help="Optional path for output file")
    args = parser.parse_args()

    logging.info(f"Input arguments: input={args.input}, output_dir={args.output_dir}, batch_idx={args.batch_idx}, output_file={args.output_file}")
    result = run_iapred(args.input, args.output_dir, args.batch_idx, args.output_file)
    exit(result)

if __name__ == "__main__":
    main()

