import os
import subprocess
import argparse
from pathlib import Path
import shutil
import sys
import psutil
import csv

# Define root directory dynamically (Vaxelan_2_0 directory)
SCRIPT_DIR = Path(__file__).resolve().parent  # modules/
ROOT_DIR = SCRIPT_DIR.parent  # Vaxelan_2_0/

# Define path to SignalP model directory
MODEL_DIR = ROOT_DIR / "tools/signalp-6.0h.fast/signalp6_fast/signalp-6-package/signalp/model_weights"

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

def find_model_dir(model_path):
    """
    Search for the model directory containing distilled_model_signalp6.pt.

    Args:
        model_path (Path): Primary path to the model directory.

    Returns:
        Path or None: Path to the model directory if found, else None.
    """
    print_status(f"Checking primary model path: {model_path}", "info")
    if model_path.is_dir() and (model_path / "distilled_model_signalp6.pt").is_file():
        print_status(f"Found model file in primary path: {model_path}", "success")
        return model_path

    # Alternative predefined paths
    alternative_paths = [
        ROOT_DIR / "tools/signalp-6.0h.fast/signalp6_fast/signalp-6-package/models",
        ROOT_DIR / "tools/signalp-6.0h.fast/models",
        Path("/home/hp-lapi/Downloads/tools/signalp-6.0h.fast/signalp6_fast/signalp-6-package/signalp/model_weights")
    ]

    for alt_path in alternative_paths:
        print_status(f"Checking alternative path: {alt_path}", "info")
        if alt_path.is_dir() and (alt_path / "distilled_model_signalp6.pt").is_file():
            print_status(f"Found model file at alternative location: {alt_path}", "success")
            return alt_path

    # Recursive search in tools directory
    print_status("Performing recursive search for distilled_model_signalp6.pt in tools directory...", "info")
    tools_dir = ROOT_DIR / "tools"
    for root, _, files in os.walk(tools_dir):
        if "distilled_model_signalp6.pt" in files:
            model_dir = Path(root)
            print_status(f"Found model file in recursive search: {model_dir}", "success")
            return model_dir

    print_status("Model file distilled_model_signalp6.pt not found in any searched paths.", "error")
    return None

def validate_fasta(fasta_file):
    """
    Validate FASTA file format and content.

    Args:
        fasta_file (Path): Path to the FASTA file.

    Returns:
        bool: True if valid, False otherwise.
    """
    try:
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
            if not lines:
                print_status("FASTA file is empty.", "error")
                return False
            if not lines[0].startswith('>'):
                print_status("FASTA file does not start with a header.", "error")
                return False
            sequence_found = False
            for line in lines[1:]:
                if line.strip() and not line.startswith('>'):
                    sequence_found = True
                    break
            if not sequence_found:
                print_status("FASTA file contains no valid sequences.", "error")
                return False
        return True
    except Exception as e:
        print_status(f"Error validating FASTA file: {e}", "error")
        return False

def convert_to_csv(input_file, output_csv):
    """
    Convert SignalP prediction_results.txt to CSV format.

    Args:
        input_file (Path): Path to the SignalP prediction_results.txt file.
        output_csv (Path): Path to the output CSV file.

    Returns:
        bool: True if successful, False otherwise.
    """
    try:
        if not input_file.is_file() or input_file.stat().st_size == 0:
            print_status(f"Input file {input_file} is missing or empty.", "error")
            return False

        # Read the input file
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Skip header lines (lines starting with '#' or empty)
        data_lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
        if not data_lines:
            print_status("No data found in prediction_results.txt.", "error")
            return False

        # Assume the first non-comment line is the header
        header = data_lines[0].split('\t')
        data = [line.split('\t') for line in data_lines[1:]]

        # Write to CSV
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)  # Write header
            writer.writerows(data)  # Write data rows

        print_status(f"Successfully converted {input_file} to {output_csv}", "success")
        return True

    except Exception as e:
        print_status(f"Error converting to CSV: {e}", "error")
        return False

def run_signalp(input_fasta, output_dir, model_dir=MODEL_DIR):
    """
    Runs SignalP-6.0 to predict signal peptides using the signalp6 executable.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Path to the output directory.
        model_dir (str): Path to the SignalP model directory.

    Returns:
        int: Exit status code (0 if successful).
    """
    try:
        # Convert paths to Path objects
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)
        model_dir = Path(model_dir)

        # Log system resources
        print_status(f"Available memory: {psutil.virtual_memory().available / (1024**3):.2f} GB", "info")
        print_status(f"CPU count: {psutil.cpu_count()}", "info")

        # Locate the signalp6 executable
        signalp_executable = shutil.which("signalp6")
        if not signalp_executable:
            print_status("Error: 'signalp6' executable not found in PATH. Ensure it is installed in your Conda environment (vaxelan_new).", "error")
            raise FileNotFoundError("signalp6 executable not found in PATH")

        # Validate inputs
        if not input_fasta.is_file():
            print_status(f"Input FASTA file not found: {input_fasta}", "error")
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

        if not validate_fasta(input_fasta):
            print_status("Invalid FASTA file format.", "error")
            return 1

        # Check for model directory with fallback search
        model_dir = find_model_dir(model_dir)
        if not model_dir:
            print_status(
                f"Model directory not found or missing distilled_model_signalp6.pt. "
                f"Ensure the model weights are in {MODEL_DIR} or specify the correct path with --model-dir. "
                f"Check Vaxelan_2_0 documentation or SignalP provider for model weights.",
                "error"
            )
            raise FileNotFoundError("Model directory not found or missing model file")

        # Check for the specific model file
        model_file = model_dir / "distilled_model_signalp6.pt"
        if not model_file.is_file():
            print_status(
                f"Required model file not found: {model_file}. "
                f"Please ensure the SignalP model weights are installed correctly. "
                f"Refer to Vaxelan_2_0 or SignalP documentation for obtaining distilled_model_signalp6.pt.",
                "error"
            )
            raise FileNotFoundError(f"Model file {model_file} not found")

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Set environment variable for model weights
        os.environ["SIGNARP_MODEL_DIR"] = str(model_dir)

        # Build the SignalP command
        command = [
            signalp_executable,
            "--fastafile", str(input_fasta),
            "--organism", "other",
            "--output_dir", str(output_dir),
            "--format", "txt",
            "--mode", "fast",
            "--bsize", "10",  # Added to reduce memory usage
            "--torch_num_threads", str(int(psutil.cpu_count()))  # Fixed: Explicitly cast to int, then to str
        ]

        # Log file size and command
        print_status(f"FASTA file size: {input_fasta.stat().st_size / (1024**2):.2f} MB", "info")
        command_str = " ".join(str(x) for x in command)  # Ensure all arguments are strings for logging
        print_status(f"Debug command: {command_str}", "info")  # Added debug logging

        # Run SignalP
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False
        )

        print_status(f"STDOUT: {result.stdout}", "info")
        if result.stderr:
            print_status(f"STDERR: {result.stderr}", "warning")
        print_status(f"SignalP returncode: {result.returncode}", "info")

        # Check for output file
        output_summary = output_dir / "prediction_results.txt"  # SignalP 6.0 default output file
        output_csv = output_dir / "prediction_results.csv"  # Output CSV file
        if result.returncode == 0:
            print_status("SignalP ran successfully.", "success")
            if output_summary.is_file() and output_summary.stat().st_size > 0:
                print_status(f"Output file created: {output_summary}", "success")
                # Convert to CSV
                if convert_to_csv(output_summary, output_csv):
                    print_status(f"CSV output created: {output_csv}", "success")
                else:
                    print_status("Failed to create CSV output.", "error")
                    return 1
            else:
                print_status(f"Output file '{output_summary}' is empty or not created.", "error")
                return 1
            return 0
        else:
            print_status(f"SignalP failed with return code {result.returncode}", "error")
            print_status(f"STDERR: {result.stderr}", "error")
            return 1

    except FileNotFoundError as e:
        print_status(f"File not found error: {e}", "error")
        return 1
    except subprocess.SubprocessError as e:
        print_status(f"Subprocess error running SignalP: {e}", "error")
        return 1
    except Exception as e:
        print_status(f"Unexpected error running SignalP: {e}", "error")
        return 1

def main():
    """Command-line interface for SignalP-6.0 prediction."""
    parser = argparse.ArgumentParser(description="Run SignalP-6.0 for signal peptide prediction.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--model-dir", default=str(MODEL_DIR), help="Path to SignalP model directory")

    args = parser.parse_args()
    return run_signalp(args.fasta, args.output, args.model_dir)

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)

