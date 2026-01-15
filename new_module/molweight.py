import os
import argparse
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd



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

def is_valid_sequence(seq: str) -> bool:
    """Validate that a sequence contains only standard amino acids."""
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.upper()
    return all(c in valid_amino_acids for c in seq)

def protparm(seq: str) -> float:
    """Compute molecular weight of a protein sequence."""
    try:
        if not is_valid_sequence(seq):
            raise ValueError("Sequence contains invalid amino acid characters")
        analysis = ProteinAnalysis(seq)
        mol_weight = round(analysis.molecular_weight(), 2)

        return mol_weight
    except Exception as e:
        print_status(f"Error computing molecular weight: {e}", "error")
        
        raise

def validate_paths(input_fasta, output_dir, output_file=None):
    """Validate input, output directory, and optional output file paths."""
    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()

    # Handle output_dir: ensure it's a directory, not a file
    if output_dir.suffix and output_dir.suffix != ".":
        output_dir = output_dir.parent
        
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
    

    # Validate output_file if provided
    if output_file:
        output_file = Path(output_file).resolve()
        if not output_file.parent.exists():
            try:
                output_file.parent.mkdir(parents=True, exist_ok=True)
                
            except Exception as e:
                raise PermissionError(f"Cannot create output file directory {output_file.parent}: {e}")
        if not os.access(output_file.parent, os.W_OK):
            raise PermissionError(f"No write permission for output file directory: {output_file.parent}")

    return input_fasta, output_dir, output_file

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

def run_molwt(input_fasta, output_dir, batch_idx=1, output_file=None):
    """
    Compute molecular weights for each sequence in the input FASTA file.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory to save output files (CSV) for each sequence.
        batch_idx (int): Batch index for naming combined output (default: 1).
        output_file (str, optional): Path to save output for a single sequence (default: None).

    Returns:
        int: Exit status code (0 if at least one sequence processed successfully, 1 otherwise).
    """
    try:
        # Validate paths
        input_fasta, output_dir, output_file = validate_paths(input_fasta, output_dir, output_file)

        # Read FASTA file
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

        # Initialize list to collect results for combined output
        combined_results = []
        overall_success = False

        # Process sequences
        for seq_id, sequence, description in sequences:
            print_status(f"Processing sequence: {seq_id}", "info")
            # Use output_file if provided, else default to standard naming
            output_csv_file = Path(output_file) if output_file and seq_id in str(output_file) else output_dir / f"{seq_id}_molwt.csv"

            try:
                # Compute molecular weight
                molweight = protparm(sequence)

                # Save individual CSV file
                output_csv_file.parent.mkdir(parents=True, exist_ok=True)
                df = pd.DataFrame({
                    "Sequence_ID": [seq_id],
                    "Description": [description],
                    "Sequence": [sequence],
                    "Molecular_Weight": [molweight]
                })
                df.to_csv(output_csv_file, index=False)
                

                # Add to combined results
                combined_results.append(df)

                # Verify individual output
                verify_output_files([output_csv_file])

                print_status(f"Results for {seq_id} saved to {output_csv_file}", "success")
                overall_success = True

            except Exception as e:
                print_status(f"Error processing sequence {seq_id}: {e}", "error")
                
                continue

        # Save combined results
        combined_file = output_dir / f"combined_molwt_batch_{batch_idx}.csv"
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            if combined_results:
                combined_df = pd.concat(combined_results, ignore_index=True)
                combined_df.to_csv(combined_file, index=False)
                print_status(f"Combined output saved to: {combined_file} with {len(combined_df)} records", "success")
                
            else:
                print_status(f"No valid sequences processed; creating empty combined output: {combined_file}", "warning")
                
                pd.DataFrame(columns=["Sequence_ID", "Description", "Sequence", "Molecular_Weight"]).to_csv(combined_file, index=False)
        except Exception as e:
            print_status(f"Failed to save combined output {combined_file}: {e}", "error")
            
            return 1

        # Verify combined file
        verify_output_files([combined_file])

        print_status("Molecular weight processing completed.", "success")
        
        return 0 if overall_success else 1

    except Exception as e:
        print_status(f"General error: {e}", "error")
        
        return 1

def main():
    """Command-line interface for computing molecular weights."""
    parser = argparse.ArgumentParser(description="Compute molecular weights from a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-d", "--output-dir", required=True, help="Directory to store results")
    parser.add_argument("--batch_idx", type=int, default=1, help="Batch index for naming combined output (default: 1)")

    args = parser.parse_args()
    return run_molwt(args.input, args.output_dir, args.batch_idx)

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)

