import os
import subprocess
import shutil
from pathlib import Path
import tempfile
import pandas as pd
import argparse
import logging
import uuid
import warnings
import sys
import sklearn

# Suppress scikit-learn version mismatch warnings
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn.base")

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Define tool path
try:
    BASE_DIR = Path(__file__).resolve().parent.parent
except NameError:
    BASE_DIR = Path.cwd().parent
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
        if not os.path.getsize(input_fasta):
            raise ValueError(f"Input FASTA file is empty: {input_fasta}")
        with open(input_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if header and sequence:
                        sequences_found = True
                        full_header = header
                        base_header = header.split('|')[1] if '|' in header else header
                        sequence_count += 1
                        logging.debug(f"Parsed sequence with full header: {full_header}, base header: {base_header}")
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
                logging.debug(f"Parsed sequence with full header: {full_header}, base header: {base_header}")
                yield full_header, base_header, ''.join(sequence)
        if not sequences_found:
            raise ValueError(f"No valid sequences found in FASTA file: {input_fasta}")
        logging.debug(f"Parsed {sequence_count} sequences from FASTA")
    except Exception as e:
        logging.error(f"Error parsing FASTA file {input_fasta}: {str(e)}")
        raise

def run_iapred_down(input_fasta, output_dir, output_file=None, batch_idx=1, uniprot_id=None, epitope_type=None):
    input_fasta = str(input_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not os.access(output_dir, os.W_OK):
        print_status(f"Output directory is not writable: {output_dir}", "error")
        return 1

    # Use output_file if provided, otherwise default to epitope_type-based name
    if output_file:
        combined_output_path = output_dir / output_file
    else:
        combined_output_path = output_dir / f"{epitope_type}_combined_iapred.csv"

    IApred_script = IAPRED_TOOL_PATH / "IApred.py"
    models_path = IAPRED_TOOL_PATH / "models"

    # Check scikit-learn version
    sklearn_version = sklearn.__version__
    if not sklearn_version.startswith("1.5"):
        print_status(f"Warning: scikit-learn version {sklearn_version} detected. Models may require version 1.5.x.", "warning")
        logging.warning(f"scikit-learn version {sklearn_version} may be incompatible with models expecting 1.5.x")

    if not IApred_script.is_file():
        print_status(f"IAPred script not found at {IApred_script}", "error")
        return 1
    if not models_path.is_dir():
        print_status(f"'models' directory not found at {models_path}", "error")
        return 1

    try:
        shutil.copytree(models_path, "models", dirs_exist_ok=True)
    except Exception as e:
        print_status(f"Failed to copy models directory: {str(e)}", "error")
        return 1

    # Use a temporary subdirectory in output_dir
    temp_subdir = output_dir / f"tmp_iapred_{uuid.uuid4().hex}"
    temp_subdir.mkdir(parents=True, exist_ok=True)
    batch_fasta = temp_subdir / f"batch_input_{uuid.uuid4().hex}.fasta"
    batch_output = temp_subdir / f"batch_iapred_out_{uuid.uuid4().hex}.csv"

    try:
        sequence_dict = {}
        for full_header, base_header, seq in parse_fasta(input_fasta):
            sequence_dict[full_header] = seq
        with open(batch_fasta, 'w') as f:
            for full_header, seq in sequence_dict.items():
                f.write(f">{full_header}\n{seq}\n")

        print_status(f"Running IAPred for {epitope_type} analysis...", "info")
        try:
            result = subprocess.run(
                ["python3", str(IApred_script), str(batch_fasta), str(batch_output), "-v"],
                check=True,
                capture_output=True,
                text=True
            )
            logging.debug(f"IAPred stdout: {result.stdout}")
            if result.stderr:
                logging.warning(f"IAPred stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            print_status(f"IAPred failed: {e.stderr}", "error")
            return 1
        except Exception as e:
            print_status(f"Unexpected error running IAPred: {str(e)}", "error")
            return 1

        if not batch_output.is_file():
            print_status(f"IAPred did not produce output file: {batch_output}", "error")
            return 1

        try:
            df = pd.read_csv(batch_output)
            if df.empty:
                print_status("IAPred output CSV is empty", "warning")
                return 1
        except Exception as e:
            print_status(f"Failed to read IAPred output CSV: {str(e)}", "error")
            return 1

        # Log CSV content for debugging
        logging.debug(f"IAPred output CSV columns: {df.columns.tolist()}")
        logging.debug(f"IAPred output CSV head:\n{df.head().to_string()}")

        # Flexible column handling
        score_col = next((col for col in df.columns if col.lower() in ['intrinsic_antigenicity_score', 'iascore', 'score']), None)
        pred_col = next((col for col in df.columns if col.lower() in ['antigenicity_category', 'prediction', 'category']), None)
        id_col = next((col for col in df.columns if col.lower() in ['header', 'id', 'sequence_id']), None)

        if not all([id_col, score_col, pred_col]):
            missing_cols = [col for col, found in [('Header', id_col), ('IAscore', score_col), ('Antigenicity_Category', pred_col)] if not found]
            print_status(f"IAPred output is missing required columns: {missing_cols}", "error")
            return 1

        df = df.rename(columns={id_col: 'Header', score_col: 'IAscore', pred_col: 'Antigenicity_Category'})

        df['Epitope_Type'] = epitope_type
        df['Header'] = df['Header'].str.rstrip('_')
        df['Epitope_Sequence'] = df['Header'].map(sequence_dict).fillna('UNMAPPED')
        df['Sequence_Length'] = df['Epitope_Sequence'].apply(lambda x: len(x) if x != 'UNMAPPED' else 0)

        df = df[['Header', 'Epitope_Type', 'Sequence_Length', 'IAscore', 'Antigenicity_Category']]

        try:
            if combined_output_path.is_file():
                df.to_csv(combined_output_path, mode='a', header=False, index=False)
            else:
                df.to_csv(combined_output_path, mode='w', header=True, index=False)
            print_status(f"Combined IAPred results saved to: {combined_output_path} with {len(df)} records", "success")
        except Exception as e:
            print_status(f"Failed to save combined output CSV: {str(e)}", "error")
            return 1

    finally:
        shutil.rmtree(temp_subdir, ignore_errors=True)
        shutil.rmtree("models", ignore_errors=True)

    print_status(f"IAPred processing completed successfully for {epitope_type}.", "success")
    return 0

def main():
    parser = argparse.ArgumentParser(description="Run IAPred for immunogenic antigen prediction.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-d", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--batch", type=int, default=1, help="Batch index")
    parser.add_argument("--output-file", help="Optional output file name")
    parser.add_argument("--eid", help="UniProt ID")
    parser.add_argument("--epitope-type", required=True, choices=["bcell", "mhci", "mhcii"], help="Epitope type")
    args = parser.parse_args()

    logging.info(f"Input arguments: input={args.input}, output_dir={args.output_dir}, batch={args.batch}, "
                 f"output_file={args.output_file}, uniprot_id={args.eid}, epitope_type={args.epitope_type}")
    result = run_iapred_down(args.input, args.output_dir, args.output_file, args.batch, args.eid, args.epitope_type)
    sys.exit(result)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("Error: No arguments provided. Use -h or --help for usage information.")
        sys.exit(1)
    main()
