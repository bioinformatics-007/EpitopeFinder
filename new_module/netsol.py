import os
import subprocess
import argparse
import traceback
from pathlib import Path
import pandas as pd
from Bio import SeqIO

# Define root directory
ROOT_DIR = Path("/home/yuktika/Downloads/Vaxelan_2_0")  # Hardcoded to ensure correct base directory

# Define default paths
NETSOLP_SCRIPT = ROOT_DIR / "tools/netsolp-1.0.ALL/predict.py"  # Points to predict.py
NETSOLP_MODELS = ROOT_DIR / "tools/netsolp-1.0.ALL/models"  # Points to models directory
NETSOLP_OUTPUT_DIR = ROOT_DIR / "results" / "NetSolP"  # Dedicated NetSolP output directory

def print_status(msg, status="info"):
    """Print colored status messages."""
    colors = {
        "info": "\033[94m",     # Blue
        "success": "\033[92m",  # Green
        "warning": "\033[93m",  # Yellow
        "error": "\033[91m"     # Red
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")

def validate_netsolp_arguments(models_path, model_type, prediction_type, netsolp_script):
    """
    Validate NetSolP arguments and script compatibility.

    Args:
        models_path (str or Path): Path to the NetSolP models directory.
        model_type (str): Type of model to use (e.g., ESM1b, ESM2).
        prediction_type (str): Type of prediction (S, U, SU).
        netsolp_script (str or Path): Path to the NetSolP predict.py script.

    Returns:
        bool: True if valid, False otherwise.
    """
    # Convert models_path and netsolp_script to Path objects
    models_path = Path(models_path)
    netsolp_script = Path(netsolp_script)

    # Check if predict.py exists
    if not netsolp_script.is_file():
        print_status(f"NetSolP script not found: {netsolp_script}", "error")
        return False

    # Check if models directory exists
    if not models_path.is_dir():
        print_status(f"Models directory not found: {models_path}", "error")
        return False

    # Check for ESM12_alphabet.pkl (commonly required by NetSolP)
    alphabet_file = models_path / "ESM12_alphabet.pkl"
    if not alphabet_file.is_file():
        print_status(f"Required model file missing: {alphabet_file}", "error")
        return False

    # Validate model_type
    valid_model_types = ["ESM1b", "ESM2"]
    if model_type not in valid_model_types:
        print_status(f"Invalid model type: {model_type}. Supported: {valid_model_types}", "error")
        return False

    # Validate prediction_type
    valid_prediction_types = ["S", "U", "SU"]
    if prediction_type not in valid_prediction_types:
        print_status(f"Invalid prediction type: {prediction_type}. Supported: {valid_prediction_types}", "error")
        return False

    return True

def reformat_output(output_file, input_fasta, output_dir):
    """
    Reformat NetSolP output to create individual CSV files for each protein sequence.

    Args:
        output_file (Path): Path to the NetSolP output CSV file.
        input_fasta (Path): Path to the input FASTA file.
        output_dir (Path): Directory to save the individual CSV files.

    Returns:
        pd.DataFrame: DataFrame with all results (for potential further use), or empty DataFrame on failure.
    """
    try:
        df = pd.read_csv(output_file, sep=',')
        print_status(f"Loaded output CSV: {output_file}", "info")
        print_status(f"Columns: {df.columns}", "info")

        # Map sequences from FASTA file
        sequences = {record.id.split("|")[1] if "|" in record.id else record.id: str(record.seq) 
                     for record in SeqIO.parse(input_fasta, "fasta")}
        
        # Ensure 'sid' exists; if not, use the first column as identifier
        id_column = 'sid' if 'sid' in df.columns else df.columns[0]
        df['Sequence'] = df[id_column].map(sequences)
        if df['Sequence'].isna().any():
            print_status(f"Some sequence IDs not found in FASTA file: {input_fasta}", "warning")

        # Select relevant columns for output
        desired_columns = [
            'sid', 'fasta', 'Sequence', 'predicted_solubility', 'predicted_usability',
            'predicted_solubility_model_0', 'predicted_solubility_model_1',
            'predicted_solubility_model_2', 'predicted_solubility_model_3',
            'predicted_solubility_model_4', 'predicted_usability_model_0',
            'predicted_usability_model_1', 'predicted_usability_model_2',
            'predicted_usability_model_3', 'predicted_usability_model_4'
        ]
        available_columns = [col for col in desired_columns if col in df.columns]
        if not available_columns:
            raise ValueError(f"No relevant columns found in {df.columns}")
        
        results_df = df[available_columns].copy()

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Generate individual CSV files for each sequence
        for _, row in results_df.iterrows():
            sequence_id = row[id_column]
            # Sanitize sequence_id to create a valid filename
            safe_sequence_id = "".join(c for c in str(sequence_id) if c.isalnum() or c in ('.', '_', '-'))
            individual_output = output_dir / f"{safe_sequence_id}.csv"
            # Create a single-row DataFrame
            single_row_df = pd.DataFrame([row], columns=available_columns)
            single_row_df.to_csv(individual_output, index=False, sep=',')
            print_status(f"Saved individual results to: {individual_output}", "success")

        return results_df

    except Exception as e:
        print_status(f"Error reformatting output: {e}", "error")
        return pd.DataFrame()

def run_netsol(input_fasta, output_dir, models_path=NETSOLP_MODELS, model_type="ESM1b", prediction_type="SU", netsolp_script=NETSOLP_SCRIPT):
    """
    Run NetSolP-1.0 to predict protein solubility and usability.

    Args:
        input_fasta (str or Path): Path to the input FASTA file.
        output_dir (str or Path): Directory to save output CSV files (individual files per sequence and combined_result.csv).
        models_path (str or Path): Path to the NetSolP models directory.
        model_type (str): Type of model to use (e.g., ESM1b, ESM2).
        prediction_type (str): Type of prediction (S for solubility, U for usability, SU for both).
        netsolp_script (str or Path): Path to the NetSolP predict.py script.

    Returns:
        int: Exit status code (0 if successful, 1 if failed).
    """
    try:
        # Convert paths to Path objects
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir) / "NetSolP"  # Append NetSolP subdirectory
        models_path = Path(models_path)
        netsolp_script = Path(netsolp_script)

        # Define output file path for NetSolP temporary output
        temp_output = output_dir / "combined_result.csv"

        # Validate inputs
        if not input_fasta.is_file():
            print_status(f"Input FASTA file not found: {input_fasta}", "error")
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

        # Validate FASTA file content
        records = list(SeqIO.parse(input_fasta, "fasta"))
        if not records:
            print_status(f"FASTA file is empty or invalid: {input_fasta}", "error")
            raise ValueError(f"FASTA file is empty or invalid: {input_fasta}")

        # Validate NetSolP arguments
        if not validate_netsolp_arguments(models_path, model_type, prediction_type, netsolp_script):
            raise ValueError("NetSolP argument validation failed")

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Construct the command
        command = [
            "python3",
            str(netsolp_script),
            "--FASTA_PATH", str(input_fasta),
            "--OUTPUT_PATH", str(temp_output),
            "--MODELS_PATH", str(models_path),
            "--MODEL_TYPE", model_type,
            "--PREDICTION_TYPE", prediction_type
        ]
        print_status(f"Running NetSolP: {' '.join(command)}", "info")

        # Run NetSolP
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False
        )

        # Log output
        print_status(f"STDOUT: {result.stdout}", "info")
        if result.stderr:
            print_status(f"STDERR: {result.stderr}", "warning")

        # Check exit code
        if result.returncode != 0:
            print_status(f"NetSolP failed with exit code {result.returncode}", "error")
            print_status(f"Command: {' '.join(command)}", "error")
            return 1

        # Verify output file
        if not temp_output.is_file() or temp_output.stat().st_size == 0:
            print_status(f"Output file not created or empty: {temp_output}", "error")
            return 1

        # Reformat output to generate individual CSV files
        results_df = reformat_output(temp_output, input_fasta, output_dir)
        if results_df.empty:
            print_status(f"No valid results generated from: {temp_output}", "error")
            return 1

        return 0

    except subprocess.CalledProcessError as e:
        print_status(f"NetSolP subprocess failed: {e}", "error")
        return 1
    except Exception as e:
        print_status(f"Unexpected error running NetSolP: {e}", "error")
        return 1

def main():
    """Command-line interface for NetSolP-1.0 prediction."""
    parser = argparse.ArgumentParser(description="Run NetSolP-1.0 on a FASTA file.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to the input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory")
    parser.add_argument("--models", default=str(NETSOLP_MODELS), help="Path to NetSolP models directory")
    parser.add_argument("--model-type", default="ESM1b", help="Model type (e.g., ESM1b, ESM2)")
    parser.add_argument("--prediction-type", default="SU", choices=["S", "U", "SU"],
                        help="Prediction type (S for solubility, U for usability, SU for both)")
    parser.add_argument("--script", default=str(NETSOLP_SCRIPT), help="Path to NetSolP predict.py script")

    args = parser.parse_args()
    return run_netsol(args.fasta, args.output, args.models, args.model_type, args.prediction_type, args.script)

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)

