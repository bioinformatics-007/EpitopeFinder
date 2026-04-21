import os
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

# Define root directory
ROOT_DIR = Path(__file__).resolve().parent.parent

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

def is_valid_sequence(seq: str) -> bool:
    """Validate that a sequence contains only standard amino acids."""
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.upper()
    return all(c in valid_amino_acids for c in seq)

def protparm(seq: str) -> tuple:
    """Calculate the instability index and determine stability."""
    sequence = str(seq).upper()
    try:
        if not sequence:
            raise ValueError("Empty sequence provided")
        if not is_valid_sequence(sequence):
            raise ValueError("Sequence contains invalid amino acid characters")
        analysis = ProteinAnalysis(sequence)
        ii = round(analysis.instability_index(), 2)
        stab = "stable" if ii < 40 else "unstable"
        stab_coff = 1 if ii < 40 else 0
        return ii, stab, stab_coff
    except Exception as e:
        print_status(f"Error computing instability index: {e}", "error")
        raise

def validate_paths(input_fasta, output_dir):
    """Validate input and output directory paths."""
    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()

    # Check input FASTA
    if not input_fasta.is_file():
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    if not os.access(input_fasta, os.R_OK):
        raise PermissionError(f"No read permission for FASTA file: {input_fasta}")

    # Check output directory
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise PermissionError(f"Cannot create output directory {output_dir}: {e}")
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"No write permission for output directory: {output_dir}")

    return input_fasta, output_dir

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

def run_instability(input_fasta, output_dir, batch_idx=1, **kwargs):
    """
    Compute instability indices for each sequence in the input FASTA file.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory to save output CSV files.
        batch_idx (int): Batch index for naming combined output (default: 1).
        **kwargs: Handle deprecated arguments gracefully.

    Returns:
        int: Exit status code (0 for success, 1 for error).
    """
    # Handle deprecated arguments
    if kwargs.get('binary_output_file') or kwargs.get('value_output_file'):
        print_status("Warning: 'binary_output_file' and 'value_output_file' are deprecated and ignored. Use CSV outputs instead.", "warning")

    try:
        # Validate paths
        input_fasta, output_dir = validate_paths(input_fasta, output_dir)

        # Define combined CSV path
        combined_csv_path = output_dir / f"combined_instability_batch_{batch_idx}.csv"

        # Read and validate FASTA file
        sequences = []
        try:
            with open(input_fasta, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    seq_id = record.id.split("|")[1] if "|" in record.id else record.id  # Extract UniProt ID
                    sequence = str(record.seq).upper()
                    if not sequence:
                        print_status(f"Empty sequence for {seq_id}; skipping", "warning")
                        continue
                    if not is_valid_sequence(sequence):
                        print_status(f"Invalid amino acid characters in {seq_id}; skipping", "warning")
                        continue
                    sequences.append((seq_id, sequence, record.description))
        except Exception as e:
            print_status(f"Error parsing FASTA file {input_fasta}: {e}", "error")
            return 1

        if not sequences:
            print_status(f"No valid sequences found in FASTA file: {input_fasta}", "error")
            return 1

        # Initialize list for combined results
        combined_results = []
        sequence_count = 0
        overall_success = False

        # Process sequences
        for seq_id, sequence, description in sequences:
            print_status(f"Processing sequence: {seq_id}", "info")
            individual_csv_path = output_dir / f"{seq_id}_instability.csv"

            try:
                # Compute instability index
                ii, stability, stability_score = protparm(sequence)
                sequence_count += 1

                # Save individual CSV
                individual_df = pd.DataFrame({
                    "Sequence_ID": [seq_id],
                    "Description": [description],
                    "Sequence": [sequence],
                    "Instability_Index": [ii],
                    "Stability": [stability],
                    "Stability_Score": [stability_score]
                })
                individual_df.to_csv(individual_csv_path, index=False)

                # Add to combined results
                combined_results.append(individual_df)

                # Verify individual output
                verify_output_files([individual_csv_path])

                print_status(f"Results for {seq_id} saved to {individual_csv_path}", "success")
                overall_success = True

            except Exception as e:
                print_status(f"Error processing sequence {seq_id}: {e}", "error")
                continue

        # Save combined CSV
        try:
            if combined_results:
                combined_df = pd.concat(combined_results, ignore_index=True)
                combined_df.to_csv(combined_csv_path, index=False)
                print_status(f"Combined output saved to: {combined_csv_path} with {len(combined_df)} records", "success")
            else:
                print_status(f"No valid sequences processed; creating empty combined output: {combined_csv_path}", "warning")
                pd.DataFrame(columns=["Sequence_ID", "Description", "Sequence", "Instability_Index", "Stability", "Stability_Score"]).to_csv(combined_csv_path, index=False)
        except Exception as e:
            print_status(f"Failed to save combined output {combined_csv_path}: {e}", "error")
            return 1

        # Verify combined output
        verify_output_files([combined_csv_path])

        if sequence_count == 0:
            print_status("No valid sequences processed", "warning")
            return 1

        print_status(f"Instability analysis completed. Processed {sequence_count} sequences.", "success")
        return 0 if overall_success else 1

    except Exception as e:
        print_status(f"General error during instability analysis: {e}", "error")
        return 1

def main():
    """Command-line interface for instability analysis."""
    parser = argparse.ArgumentParser(description="Compute protein stability from a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-d", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--batch_idx", type=int, default=1, help="Batch index for naming combined output (default: 1)")

    args = parser.parse_args()
    return run_instability(args.input, args.output_dir, args.batch_idx)

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)

