import os
import subprocess
import pandas as pd
from pathlib import Path
import sys
import argparse
import shutil
import uuid
from Bio import SeqIO

def print_status(msg, status="info"):
    colors = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")

def is_csv_empty(csv_path):
    if not csv_path.exists():
        return True
    try:
        df = pd.read_csv(csv_path)
        return df.empty or len(df) == 0
    except Exception:
        return True

def get_fasta_sequences(fasta_file):
    try:
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            seq_id = record.id.split("|")[0] if "|" in record.id else record.id
            sequences.append((seq_id, seq))
        return sequences
    except Exception as e:
        print_status(f"Error reading FASTA file {fasta_file}: {e}", "error")
        return []

def validate_fasta_file(fasta_file):
    try:
        if not fasta_file.exists():
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
            seq_id = record.id.split("|")[0] if "|" in record.id else record.id
            if not sequence:
                print_status(f"Empty sequence for ID {seq_id}. Skipping.", "warning")
                continue
            if len(sequence) < 5:
                print_status(f"Sequence too short for ID {seq_id} (length {len(sequence)}). Skipping.", "warning")
                continue
            if not all(c in "ACDEFGHIKLMNPQRSTVWY" for c in sequence):
                print_status(f"Invalid characters in sequence {seq_id}. Skipping.", "warning")
                continue
            valid_sequences.append((seq_id, sequence))

        if not valid_sequences:
            return False, None, None, "No valid sequences found in FASTA file"

        temp_dir = fasta_file.parent / f"tmp_dir_{uuid.uuid4().hex}"
        temp_dir.mkdir(parents=True, exist_ok=True)

        clean_fasta = temp_dir / f"input_{uuid.uuid4().hex}.fasta"
        with open(clean_fasta, 'w') as f:
            for seq_id, sequence in valid_sequences:
                f.write(f">{seq_id}\n{sequence}\n")

        if not clean_fasta.exists():
            return False, None, None, f"Failed to create cleaned FASTA file: {clean_fasta}"

        print_status(f"FASTA file validated and cleaned: {len(valid_sequences)} sequences saved to {clean_fasta}", "success")
        return True, clean_fasta, temp_dir, None
    except Exception as e:
        return False, None, None, f"Error validating FASTA file: {e}"

def verify_output_files(output_files):
    for file_path in output_files:
        if file_path.exists():
            file_size = file_path.stat().st_size
            if file_size > 0:
                print_status(f"Output file found and non-empty: {file_path} ({file_size} bytes)", "success")
            else:
                print_status(f"Output file is empty: {file_path}", "warning")
        else:
            print_status(f"Output file not found: {file_path}", "error")

def run_algpred_down(input_fasta, output_dir, batch_idx=1, output_file=None, model="1", epitope_type=None):
    input_fasta = Path(input_fasta)
    output_dir = Path(output_dir)
    temp_dir = None

    try:
        if not output_dir:
            raise ValueError("output_dir must be provided and cannot be None or empty")

        output_dir.mkdir(parents=True, exist_ok=True)
        if not os.access(output_dir, os.W_OK):
            raise PermissionError(f"Output directory is not writable: {output_dir}")

        is_valid, clean_fasta, temp_dir, error_msg = validate_fasta_file(input_fasta)
        if not is_valid:
            print_status(f"{error_msg}", "error")
            return 1

        if not clean_fasta.exists():
            print_status(f"No cleaned FASTA file: {clean_fasta}", "error")
            return 1

        input_sequences = get_fasta_sequences(clean_fasta)
        if not input_sequences:
            print_status(f"No valid sequences found in {clean_fasta}", "error")
            return 1

        combined_results = []
        expected_columns = ["Sequence_ID", "Epitope_Type", "Sequence", "Sequence_Length", "Score", "Prediction"]

        thresholds = ["0.3", "0.5", "0.7"]
        success = False

        for thresh in thresholds:
            print_status(f"Running AlgPred2 on batch {batch_idx}, threshold {thresh} (Model {model})", "info")
            temp_output = temp_dir / f"algpred_temp_{batch_idx}_{uuid.uuid4().hex}.csv"
            debug_output = output_dir / f"algpred_model_{model}_thresh_{thresh}.txt"
            debug_error = output_dir / f"algpred_model_{model}_thresh_{thresh}.error"

            if not temp_dir.exists():
                print_status(f"Temp directory not found: {temp_dir}", "error")
                return 1

            command = [
                "conda", "run", "-n", "vaxelan_new", "algpred2",
                "-i", str(clean_fasta),
                "-o", str(temp_output),
                "-m", model,
                "-t", thresh
            ]
            command_str = " ".join(command)
            print_status(f"Executing: {command_str}", "info")
            try:
                result = subprocess.run(
                    command_str,
                    shell=True,
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=600,
                    cwd=str(temp_dir)
                )

                with open(debug_output, "w") as f:
                    f.write(result.stdout)
                with open(debug_error, "w") as f:
                    f.write(result.stderr or "No stderr output")

                if result.stderr:
                    print_status(f"AlgPred2 stderr: {result.stderr}", "warning")

                if temp_output.exists() and not is_csv_empty(temp_output):
                    try:
                        df = pd.read_csv(temp_output)
                        sequence_col = next((col for col in df.columns if col.lower() in ['sequence', 'peptide']), None)
                        score_col = next((col for col in df.columns if col.lower() in ['score', 'ml_score', 'prediction_score']), None)
                        pred_col = next((col for col in df.columns if col.lower() in ['prediction', 'class']), None)

                        if not sequence_col or not score_col or not pred_col:
                            print_status(f"Warning: Expected columns not found in {temp_output}. Columns: {df.columns.tolist()}", "warning")
                            continue

                        df = df.rename(columns={sequence_col: "Sequence", score_col: "Score", pred_col: "Prediction"})
                        df["Sequence_ID"] = df["Sequence"].map(
                            {seq: seq_id for seq_id, seq in input_sequences}
                        )
                        df["Epitope_Type"] = str(epitope_type) if epitope_type else "NA"
                        df["Sequence_Length"] = df["Sequence"].apply(len)
                        cols = ["Sequence_ID", "Epitope_Type", "Sequence", "Sequence_Length", "Score", "Prediction"]
                        df = df[cols]

                        output_sequences = set(df["Sequence"].str.strip())
                        missing_sequences = [(seq_id, seq) for seq_id, seq in input_sequences if seq not in output_sequences]

                        if missing_sequences:
                            non_allergen_df = pd.DataFrame({
                                "Sequence_ID": [seq_id for seq_id, _ in missing_sequences],
                                "Epitope_Type": [epitope_type if epitope_type else "NA" for _ in missing_sequences],
                                "Sequence": [seq for _, seq in missing_sequences],
                                "Sequence_Length": [len(seq) for _, seq in missing_sequences],
                                "Score": [0.0] * len(missing_sequences),
                                "Prediction": ["Non-Allergen"] * len(missing_sequences)
                            })
                            df = pd.concat([df, non_allergen_df], ignore_index=True)

                        combined_results.append(df)
                        print_status(f"Success (Model {model}, Threshold {thresh}): {len(df)} sequences", "success")
                        success = True
                        break
                    except Exception as e:
                        print_status(f"Error processing CSV {temp_output}: {str(e)}", "error")
                        continue
                else:
                    print_status(f"Output empty or missing: {temp_output}", "warning")
                    non_allergen_df = pd.DataFrame({
                        "Sequence_ID": [seq_id for seq_id, _ in input_sequences],
                        "Epitope_Type": [epitope_type if epitope_type else "NA" for _ in input_sequences],
                        "Sequence": [seq for _, seq in input_sequences],
                        "Sequence_Length": [len(seq) for _, seq in input_sequences],
                        "Score": [0.0] * len(input_sequences),
                        "Prediction": ["Non-Allergen"] * len(input_sequences)
                    })
                    combined_results.append(non_allergen_df)
                    print_status(f"No allergens: {len(input_sequences)} non-allergen sequences", "success")
                    success = True
                    break

            except subprocess.CalledProcessError as e:
                print_status(f"AlgPred2 failed (Threshold {thresh}): {e.stderr}", "error")
                continue
            except Exception as e:
                print_status(f"Unexpected error (Threshold {thresh}): {str(e)}", "error")
                continue

        combined_file = output_dir / output_file if output_file else output_dir / f"{epitope_type}_Combined_AlgPred_{batch_idx}.csv" if epitope_type else output_dir / f"combined_algpred_batch_{batch_idx}.csv"
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            if combined_results:
                combined_df = pd.concat(combined_results, ignore_index=True)
                for col in expected_columns:
                    if col not in combined_df.columns:
                        combined_df[col] = pd.NA
                combined_df = combined_df[expected_columns]
                combined_df.to_csv(combined_file, index=False)
                print_status(f"Saved {len(combined_df)} to {combined_file}", "success")
            else:
                print_status(f"No valid sequences processed; empty output: {combined_file}", "warning")
                pd.DataFrame(columns=expected_columns).to_csv(combined_file, index=False)
        except Exception as e:
            print_status(f"Failed to save {combined_file}: {str(e)}", "error")
            return 1

        output_files = [combined_file]
        verify_output_files(output_files)

        if success:
            print_status(f"AlgPred2 completed for {epitope_type if epitope_type else 'batch'}.", "success")
            return 0
        else:
            if model == "1":
                print_status("Retrying with Model 2...", "info")
                return run_algpred_down(input_fasta, output_dir, batch_idx, output_file=output_file, model="2", epitope_type=epitope_type)
            print_status("AlgPred2 did not produce output.", "error")
            return 2

    except Exception as e:
        print_status(f"Error in run_algpred_down: {str(e)}", "error")
        return 1

    finally:
        if temp_dir and temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)

def main():
    parser = argparse.ArgumentParser(description="Run AlgPred2 for allergenicity.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA")
    parser.add_argument("-d", "--output-dir", required=True, help="Output dir")
    parser.add_argument("--batch_idx", type=int, default=1, help="Batch index")
    parser.add_argument("--output-file", help="Output file name")
    parser.add_argument("--epitope-type", help="Epitope type")

    args = parser.parse_args()

    if not Path(args.input).exists():
        print_status(f"Error: Input FASTA not found: {args.input}", "error")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            print_status(f"Error: Failed to create dir {args.output_dir}: {str(e)}", "error")
            sys.exit(1)

    result = run_algpred_down(args.input, args.output_dir, args.batch_idx, args.output_file, epitope_type=args.epitope_type)
    sys.exit(result)

if __name__ == "__main__":
    main()


