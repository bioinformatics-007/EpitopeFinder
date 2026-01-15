import subprocess
import os
import sys
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Define root directory
ROOT_DIR = Path(__file__).resolve().parent.parent

# Path to NetChop executable (relative to project root)
NETCHOP_SCRIPT_PATH = os.path.join(
    ROOT_DIR,
    "tools/IEDB_NetChop-3.0/netchop/predict.py"
)

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

def get_sequence_id(input_fasta):
    """
    Extract the first sequence ID from the FASTA file for plot naming.

    Args:
        input_fasta (str): Path to the input FASTA file.

    Returns:
        str: Formatted sequence ID or 'unknown' if not found.
    """
    try:
        with open(input_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    seq_id = line.strip()[1:].replace(' ', '_').replace('=', '_').replace(':', '_')[:100]  # Limit length, minimal replacements
                    return seq_id
        return 'unknown'
    except Exception as e:
        print_status(f"Error reading FASTA for sequence ID: {e}", "warning")
        return 'unknown'

def plot_results(csv_file, output_dir, method, seq_id):
    """
    Generate a bar plot for the top 10 peptides/amino acids by score from a NetChop CSV file.

    Args:
        csv_file (str): Path to the CSV file.
        output_dir (str): Directory to save the plot.
        method (str): Method name (e.g., 'netchop', 'netctlpan') for naming the plot.
        seq_id (str): Sequence ID for plot naming.
    """
    try:
        print_status(f"Using output_dir for plot: {output_dir}", "info")
        # Read CSV file
        df = pd.read_csv(csv_file, header=None)
        
        # Find the header row (starts with '#')
        header_row = df[df[0].str.startswith('#')].iloc[0][0]
        headers = header_row.lstrip('#').split('\t')
        
        # Filter data rows (exclude header and metadata)
        data_rows = df[~df[0].str.startswith(('#', 'tr|'))].copy()
        
        # Parse tab-separated data
        data_rows['split'] = data_rows[0].str.split('\t')
        data_rows = pd.DataFrame(data_rows['split'].tolist(), columns=headers, index=data_rows.index)
        
        print_status(f"Columns in {csv_file}: {headers}", "info")
        print_status(f"First few rows:\n{data_rows.head().to_string()}", "info")
        
        if method == 'netchop':
            peptide_col = 'amino_acid'
            score_col = 'prediction_score'
        else:  # netctlpan
            peptide_col = 'peptide'
            score_col = 'combined_prediction_score'
        
        if peptide_col not in headers or score_col not in headers:
            print_status(f"CSV file {csv_file} missing {peptide_col} or {score_col}. Found columns: {headers}", "error")
            return

        # Convert score column to numeric
        data_rows[score_col] = pd.to_numeric(data_rows[score_col], errors='coerce')
        
        # Sort by score and take top 10
        df = data_rows.sort_values(by=score_col, ascending=False).head(10)
        peptides = df[peptide_col]
        scores = df[score_col]

        # Create bar plot
        plt.figure(figsize=(10, 6))
        plt.bar(peptides, scores, color='skyblue')
        plt.xlabel('Peptide' if method == 'netctlpan' else 'Amino Acid')
        plt.ylabel('Prediction Score')
        plt.title(f'Top 10 {"Peptides" if method == "netctlpan" else "Amino Acids"} - {method.capitalize()}')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        # Save plot
        plot_path = os.path.join(output_dir, f"plot_sequence_1_{seq_id}_{method}.png")
        plt.savefig(plot_path)
        plt.close()
        print_status(f"Plot saved to: {plot_path}", "success")
    except Exception as e:
        print_status(f"Error generating plot for {method}: {e}", "error")

def run_netchop(input_fasta, output_dir=".", method="netchop", threshold=None, output_file=None):
    """
    Runs the NetChop prediction tool for a single method.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory where output should be saved (default: current directory).
        method (str): Prediction method (e.g., 'netchop', 'netctlpan').
        threshold (float): Optional threshold value for predictions (default: None).
        output_file (str, optional): Custom output file name (default: None, uses method_out.csv).

    Returns:
        int: Exit status code (0 if successful).
    """
    if not os.path.isfile(input_fasta):
        print_status(f"Input FASTA file not found: {input_fasta}", "error")
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    if not os.path.isfile(NETCHOP_SCRIPT_PATH):
        print_status(f"NetChop executable not found: {NETCHOP_SCRIPT_PATH}", "error")
        raise FileNotFoundError(f"NetChop executable not found: {NETCHOP_SCRIPT_PATH}")

    os.makedirs(output_dir, exist_ok=True)
    # Use provided output_file with method suffix or default to method_out.csv
    if output_file:
        output_file = f"{output_file}_{method}.csv"
    else:
        output_file = f"{method}_out.csv"
    output_path = os.path.join(output_dir, output_file)

    # Get sequence ID for plot naming
    seq_id = get_sequence_id(input_fasta)

    print_status(f"Running {method} on: {input_fasta}", "info")
    print_status(f"Output will be saved to: {output_path}", "info")
    print_status(f"Executing command: {sys.executable} {NETCHOP_SCRIPT_PATH} -m {method} {input_fasta}", "info")

    # Build the command with input FASTA directly
    command = [sys.executable, NETCHOP_SCRIPT_PATH, "-m", method, input_fasta]
    if threshold is not None:
        command.extend(["--threshold", str(threshold)])

    try:
        # Run the command and capture stdout
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # Write the captured stdout to the output file
        with open(output_path, "w") as out:
            out.write(result.stdout)
        print_status(f"{method} ran successfully. Output saved to: {output_path}", "success")
        # Generate plot for the results
        plot_results(output_path, output_dir, method, seq_id)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print_status(f"Error occurred for {method}: {e.stderr}", "error")
        raise RuntimeError(f"{method} execution failed: {e.stderr}")
    except Exception as e:
        print_status(f"Unexpected error for {method}: {e}", "error")
        raise

def run_all_methods(input_fasta, output_dir=".", threshold=None, output_file=None):
    """
    Runs the NetChop prediction tool for all specified methods.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory where output should be saved (default: current directory).
        threshold (float): Optional threshold value for predictions (default: None).
        output_file (str, optional): Custom output file name prefix (default: None).

    Returns:
        dict: Mapping of method names to their exit status codes.
    """
    methods = ["netchop", "netctlpan"]
    results = {}

    for method in methods:
        try:
            returncode = run_netchop(input_fasta, output_dir, method, threshold, output_file)
            results[method] = returncode
        except Exception as e:
            results[method] = -1  # Indicate failure
            print_status(f"Skipping {method} due to error.", "warning")

    return results

def main():
    """Command-line interface for NetChop predictions."""
    parser = argparse.ArgumentParser(description="Run NetChop on a FASTA file for multiple methods.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-d", "--output_dir", default=".", help="Output directory")
    parser.add_argument("-t", "--threshold", type=float, help="Prediction threshold (optional)")
    parser.add_argument("-o", "--output_file", help="Custom output file name prefix (optional)")

    args = parser.parse_args()
    print_status(f"Output directory set to: {args.output_dir}", "info")
    results = run_all_methods(args.input, args.output_dir, args.threshold, args.output_file)

    # Summarize results
    seq_id = get_sequence_id(args.input)
    for method, returncode in results.items():
        if returncode == 0:
            print_status(f"Completed {method} successfully. Output in {os.path.join(args.output_dir, f'{method}_out.csv' if not args.output_file else f'{args.output_file}_{method}.csv')}", "success")
            print_status(f"Plot in {os.path.join(args.output_dir, f'plot_sequence_1_{seq_id}_{method}.png')}", "success")
        else:
            print_status(f"Failed to complete {method}.", "error")

if __name__ == "__main__":
    main()
