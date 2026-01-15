import os
import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import sys
import argparse
import tempfile
import shutil
import uuid

# Define root directory
ROOT_DIR = Path.cwd().parent

def print_status(msg, status="info"):
    """Print colored status messages."""
    colors = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")

def is_csv_empty(csv_path):
    """Check if a CSV file is empty or contains only headers."""
    if not os.path.isfile(csv_path):
        return True
    try:
        df = pd.read_csv(csv_path)
        return df.empty or len(df) == 0
    except Exception:
        return True

def get_fasta_sequences(fasta_file):
    """Read all sequences from a FASTA file with UniProt IDs."""
    try:
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            seq_id = record.id.split("|")[1] if "|" in record.id else record.id  # Extract UniProt ID
            sequences.append((seq_id, seq))
        print_status(f"Extracted {len(sequences)} sequences from {fasta_file}", "info")
        return sequences
    except Exception as e:
        print_status(f"Error reading FASTA file {fasta_file}: {e}", "error")
        return []

def validate_fasta_file(fasta_file):
    """Validate and clean the FASTA file for AlgPred2 compatibility."""
    try:
        # Check if input file exists and is readable
        if not os.path.isfile(fasta_file):
            return False, None, None, f"Input FASTA file not found: {fasta_file}"
        
        with open(fasta_file, 'r') as f:
            content = f.read().strip()
            if not content:
                return False, None, None, "FASTA file is empty"
            if not content.startswith('>'):
                return False, None, None, "FASTA file does not start with '>'"

        valid_sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq).upper()
            seq_id = record.id.split("|")[1] if "|" in record.id else record.id  # Extract UniProt ID
            # Enhanced validation: check sequence length and characters
            if not sequence:
                print_status(f"Empty sequence for ID {seq_id}. Skipping.", "warning")
                continue
            if len(sequence) < 5:
                print_status(f"Sequence too short for ID {seq_id} (length {len(sequence)}). Skipping.", "warning")
                continue
            if not all(c in "ACDEFGHIKLMNPQRSTVWY" for c in sequence):
                print_status(f"Invalid characters in sequence {seq_id}: {sequence}. Skipping.", "warning")
                continue
            valid_sequences.append((seq_id, sequence))

        if not valid_sequences:
            return False, None, None, "No valid sequences found in FASTA file"

        # Create temporary directory
        temp_dir = tempfile.mkdtemp(prefix="vaxelan_")
        if not os.path.exists(temp_dir):
            return False, None, None, f"Failed to create temporary directory: {temp_dir}"

        # Write cleaned FASTA file to temporary directory
        clean_fasta = os.path.join(temp_dir, f"input_{uuid.uuid4().hex}.fasta")
        with open(clean_fasta, 'w') as f:
            for seq_id, sequence in valid_sequences:
                f.write(f">{seq_id}\n{sequence}\n")

        # Verify cleaned FASTA file was created
        if not os.path.isfile(clean_fasta):
            return False, None, temp_dir, f"Failed to create cleaned FASTA file: {clean_fasta}"

        print_status(f"FASTA file validated and cleaned: {len(valid_sequences)} sequences saved to {clean_fasta}", "success")
        return True, clean_fasta, temp_dir, None
    except Exception as e:
        return False, None, None, f"Error validating FASTA file: {e}"

def verify_output_files(output_files):
    """Verify that output files exist and are non-empty."""
    for file_path in output_files:
        if file_path.exists():
            file_size = file_path.stat().st_size
            if file_size > 0:
                print_status(f"Output file found and non-empty: {file_path} ({file_size} bytes)", "success")
            else:
                print_status(f"Output file is empty: {file_path}", "warning")
        else:
            print_status(f"Output file not found: {file_path}", "error")

def run_algpred(input_fasta, output_dir, batch_idx=1, output_file=None, model="1"):
    """
    Runs AlgPred2 via Conda environment to predict allergenicity, ensuring both allergens and non-allergens are output.

    Args:
        input_fasta (str or Path): Path to the input FASTA file.
        output_dir (str or Path): Directory where output should be saved.
        batch_idx (int): Batch index for naming combined output (default: 1).
        output_file (str, optional): Path to save output for a single sequence (default: None).
        model (str): AlgPred2 model (1 or 2, default 1).

    Returns:
        int: Exit status code (0 for success, 1 for error, 2 for no valid output).
    """
    temp_dir = None

    try:
        if not output_dir:
            raise ValueError("output_dir must be provided and cannot be None or empty")

        # Use the provided output_dir directly without appending 'algpred'
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        if not os.access(output_dir, os.W_OK):
            raise PermissionError(f"Output directory is not writable: {output_dir}")

        input_fasta = str(input_fasta)

        is_valid, clean_fasta, temp_dir, error_msg = validate_fasta_file(input_fasta)
        if not is_valid:
            print_status(f"Invalid FASTA file {input_fasta}: {error_msg}", "error")
            return 1

        if not os.path.isfile(clean_fasta):
            print_status(f"Cleaned FASTA file not found: {clean_fasta}", "error")
            return 1

        input_sequences = get_fasta_sequences(clean_fasta)
        if not input_sequences:
            print_status(f"No valid sequences found in cleaned FASTA file: {clean_fasta}", "error")
            return 1

        combined_results = []
        expected_columns = ["Sequence_ID", "Sequence", "Score", "Prediction"]

        thresholds = ["0.3", "0.5", "0.7"]
        success = False

        for thresh in thresholds:
            print_status(f"Running AlgPred2 on: {clean_fasta} (Model {model}, Threshold {thresh})", "info")
            temp_output = os.path.join(temp_dir, f"algpred_temp_{uuid.uuid4().hex}.csv")

            if not os.path.exists(temp_dir):
                print_status(f"Temporary directory not found: {temp_dir}", "error")
                return 1

            command = [
                "conda", "run", "-n", "vaxelan_new",
                "algpred2",
                "-i", clean_fasta,
                "-o", temp_output,
                "-m", model,
                "-t", thresh
            ]
            command_str = " ".join(command)
            print_status(f"Running command: {command_str}", "info")

            try:
                result = subprocess.run(
                    command_str,
                    shell=True,
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=600,
                    cwd=temp_dir
                )
                print_status(f"AlgPred2 stdout: {result.stdout}", "info")
                if result.stderr:
                    print_status(f"AlgPred2 stderr: {result.stderr}", "warning")

                if os.path.isfile(temp_output) and not is_csv_empty(temp_output):
                    try:
                        df = pd.read_csv(temp_output)
                        sequence_col = next((col for col in df.columns if col.lower() in ['sequence', 'peptide']), None)
                        score_col = next((col for col in df.columns if col.lower() in ['score', 'ml_score', 'prediction_score']), None)
                        pred_col = next((col for col in df.columns if col.lower() in ['prediction', 'class']), None)

                        if not sequence_col or not score_col or not pred_col:
                            print_status(f"Expected columns not found in {temp_output}. Columns: {df.columns.tolist()}", "warning")
                            continue  # Try next threshold

                        df = df.rename(columns={sequence_col: "Sequence", score_col: "Score", pred_col: "Prediction"})
                        df["Sequence_ID"] = df["Sequence"].map(
                            {seq: seq_id for seq_id, seq in input_sequences}
                        )
                        cols = ["Sequence_ID"] + [col for col in df.columns if col != "Sequence_ID"]
                        df = df[cols]

                        output_sequences = set(df["Sequence"].str.strip())
                        missing_sequences = [(seq_id, seq) for seq_id, seq in input_sequences if seq not in output_sequences]

                        if missing_sequences:
                            non_allergen_df = pd.DataFrame({
                                "Sequence_ID": [seq_id for seq_id, _ in missing_sequences],
                                "Sequence": [seq for _, seq in missing_sequences],
                                "Score": [0.0] * len(missing_sequences),
                                "Prediction": ["Non-Allergen"] * len(missing_sequences)
                            })
                            df = pd.concat([df, non_allergen_df], ignore_index=True)
                            print_status(f"Added {len(missing_sequences)} non-allergen sequences to output", "info")

                        for seq_id, seq in input_sequences:
                            seq_df = df[df["Sequence_ID"] == seq_id]
                            if not seq_df.empty:
                                output_csv_file = Path(output_file) if output_file and seq_id in str(output_file) else output_dir / f"{seq_id}_algpred.csv"
                                output_csv_file.parent.mkdir(parents=True, exist_ok=True)
                                seq_df.to_csv(output_csv_file, index=False)
                                combined_results.append(seq_df)
                                print_status(f"Individual CSV saved to: {output_csv_file} with {len(seq_df)} records", "success")
                            else:
                                print_status(f"No results for {seq_id}; creating non-allergen CSV", "warning")
                                non_allergen_df = pd.DataFrame({
                                    "Sequence_ID": [seq_id],
                                    "Sequence": [seq],
                                    "Score": [0.0],
                                    "Prediction": ["Non-Allergen"]
                                })
                                output_csv_file = Path(output_file) if output_file and seq_id in str(output_file) else output_dir / f"{seq_id}_algpred.csv"
                                output_csv_file.parent.mkdir(parents=True, exist_ok=True)
                                non_allergen_df.to_csv(output_csv_file, index=False)
                                combined_results.append(non_allergen_df)
                                print_status(f"Non-allergen CSV saved to: {output_csv_file}", "success")

                        print_status(f"AlgPred2 completed (Model {model}, Threshold {thresh}). Processed {len(df)} sequences", "success")
                        success = True
                        break  # Exit threshold loop on success
                    except Exception as e:
                        print_status(f"Error processing output CSV {temp_output}: {type(e).__name__}: {e}", "error")
                        continue  # Try next threshold
                else:
                    print_status(f"Output file empty or not created: {temp_output}", "warning")

                    for seq_id, seq in input_sequences:
                        non_allergen_df = pd.DataFrame({
                            "Sequence_ID": [seq_id],
                            "Sequence": [seq],
                            "Score": [0.0],
                            "Prediction": ["Non-Allergen"]
                        })
                        output_csv_file = Path(output_file) if output_file and seq_id in str(output_file) else output_dir / f"{seq_id}_algpred.csv"
                        output_csv_file.parent.mkdir(parents=True, exist_ok=True)
                        non_allergen_df.to_csv(output_csv_file, index=False)
                        combined_results.append(non_allergen_df)
                        print_status(f"No allergens detected for {seq_id}. Saved non-allergen CSV to: {output_csv_file}", "success")

                    print_status(f"No allergens detected. Saved {len(input_sequences)} non-allergen sequences", "success")
                    success = True
                    break  # Exit threshold loop

            except subprocess.CalledProcessError as e:
                error_msg = e.stderr or "No stderr output"
                print_status(f"AlgPred2 execution failed (Model {model}, Threshold {thresh}): {error_msg} (Return code: {e.returncode})", "error")
                continue  # Try next threshold
            except subprocess.TimeoutExpired as e:
                print_status(f"AlgPred2 timed out after 10 minutes (Model {model}, Threshold {thresh}): {e}", "error")
                continue  # Try next threshold
            except Exception as e:
                print_status(f"Unexpected error running AlgPred2 (Model {model}, Threshold {thresh}): {type(e).__name__}: {e}", "error")
                continue  # Try next threshold

        combined_file = output_dir / f"combined_algpred_batch_{batch_idx}.csv"
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            if combined_results:
                combined_df = pd.concat(combined_results, ignore_index=True)
                for col in expected_columns:
                    if col not in combined_df.columns:
                        combined_df[col] = pd.NA
                combined_df = combined_df[expected_columns]
                combined_df.to_csv(combined_file, index=False)
                print_status(f"Combined output saved to: {combined_file} with {len(combined_df)} records", "success")
            else:
                print_status(f"No valid sequences processed; creating empty combined output: {combined_file}", "warning")
                pd.DataFrame(columns=expected_columns).to_csv(combined_file, index=False)
        except Exception as e:
            print_status(f"Failed to save combined output {combined_file}: {e}", "error")
            return 1

        output_files = [output_dir / f"{seq_id}_algpred.csv" for seq_id, _ in input_sequences] + [combined_file]
        verify_output_files(output_files)

        if success:
            print_status("AlgPred2 processing completed.", "success")
            return 0
        else:
            if model == "1":
                print_status("All thresholds failed with Model 1. Retrying with Model 2...", "info")
                return run_algpred(input_fasta, output_dir, batch_idx, output_file, model="2")
            print_status("AlgPred2 did not produce usable output after all attempts.", "error")
            return 2

    except Exception as e:
        print_status(f"Unexpected error in run_algpred: {type(e).__name__}: {e}", "error")
        return 1

    finally:
        if temp_dir and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)
            print_status(f"Cleaned up temporary directory: {temp_dir}", "info")

def main():
    """Command-line interface for running AlgPred2."""
    parser = argparse.ArgumentParser(description="Run AlgPred2 for allergenicity prediction.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-d", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--batch_idx", type=int, default=1, help="Batch index for naming combined output (default: 1)")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print_status(f"Input FASTA file does not exist: {args.input}", "error")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            print_status(f"Failed to create output directory {args.output_dir}: {e}", "error")
            sys.exit(1)

    if not os.access(args.output_dir, os.W_OK):
        print_status(f"Output directory is not writable: {args.output_dir}", "error")
        sys.exit(1)

    result = run_algpred(args.input, args.output_dir, args.batch_idx)
    sys.exit(result)

if __name__ == "__main__":
    main()
