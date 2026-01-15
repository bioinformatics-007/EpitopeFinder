import subprocess
import os
import tempfile
import sys
import argparse
import csv
from pathlib import Path
import psutil  # For memory monitoring

# Define root directory
ROOT_DIR = Path(__file__).resolve().parent.parent

# Path to IEDB B-cell prediction script (relative to project root)
IEDB_TOOL_PATH = os.path.join(
    ROOT_DIR,
    "tools/IEDB_BCell-3.1/bcell_standalone/predict_antibody_epitope.py"
)

# Maximum sequence length warning threshold (adjust as needed)
MAX_SEQUENCE_LENGTH = 100000  # Warn for sequences >100,000 amino acids

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

def check_system_resources():
    """Check available memory to warn about potential issues."""
    available_memory = psutil.virtual_memory().available / (1024 ** 2)  # MB
    if available_memory < 1024:  # Less than 1GB
        print_status("Low system memory detected (<1GB). Long sequences may cause issues.", "warning")
    return available_memory

def read_fasta(fasta_path):
    """Reads the first sequence from a FASTA file, streaming to avoid memory overload."""
    try:
        sequence = []
        with open(fasta_path, "r") as file:
            for line in file:
                line = line.strip()
                if not line.startswith(">"):
                    sequence.append(line)
        sequence = "".join(sequence)
        if not sequence:
            raise ValueError("FASTA file is empty or contains no valid sequence")
        seq_length = len(sequence)
        print_status(f"Read sequence of length {seq_length}", "info")
        if seq_length > MAX_SEQUENCE_LENGTH:
            print_status(
                f"Sequence length ({seq_length} amino acids) exceeds recommended limit ({MAX_SEQUENCE_LENGTH}). "
                "Processing may be slow or memory-intensive.",
                "warning"
            )
        return sequence
    except FileNotFoundError:
        print_status(f"Input FASTA file not found: {fasta_path}", "error")
        raise
    except Exception as e:
        print_status(f"Error reading FASTA file: {e}", "error")
        raise

def query_iedb(sequence, method="Emini", chunk_size=50000, window_size=10):
    """Run IEDB prediction tool with the given sequence, method, and window size, processing in chunks if needed."""
    if not os.path.isfile(IEDB_TOOL_PATH):
        print_status(f"IEDB tool not found: {IEDB_TOOL_PATH}", "error")
        raise FileNotFoundError(f"IEDB tool not found: {IEDB_TOOL_PATH}")

    seq_length = len(sequence)
    if seq_length > chunk_size and method in ["Bepipred", "Bepipred2"]:
        print_status(
            f"Sequence length ({seq_length}) is large for {method}. Processing in chunks of {chunk_size} residues.",
            "info"
        )
        return process_in_chunks(sequence, method, chunk_size, window_size)

    # Single chunk processing for smaller sequences or non-intensive methods
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as temp_file:
        temp_file.write(f">sequence\n{sequence}\n")
        temp_file_path = temp_file.name

    command = ['python3', IEDB_TOOL_PATH, '-f', temp_file_path, '-m', method, '-w', str(window_size)]
    print_status(f"Running command: {' '.join(command)}", "info")

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        os.unlink(temp_file_path)  # Clean up temporary file
        if not result.stdout.strip():
            print_status("IEDB tool returned empty output", "warning")
        else:
            print_status("IEDB tool output received. Writing raw output to debug.log", "info")
            with open("debug.log", "w") as f:
                f.write(result.stdout)
        return result.stdout
    except subprocess.CalledProcessError as e:
        os.unlink(temp_file_path)
        print_status(f"IEDB command failed: {e.stderr}", "error")
        raise RuntimeError(f"Error querying IEDB: {e.stderr}")
    except Exception as e:
        os.unlink(temp_file_path)
        print_status(f"Unexpected error: {e}", "error")
        raise

def process_in_chunks(sequence, method, chunk_size, window_size):
    """Process long sequences in chunks to reduce memory usage."""
    results = []
    seq_length = len(sequence)
    for start in range(0, seq_length, chunk_size):
        end = min(start + chunk_size, seq_length)
        chunk = sequence[start:end]
        print_status(f"Processing chunk {start+1}-{end} of {seq_length} residues", "info")
        chunk_output = query_iedb(chunk, method, chunk_size, window_size)
        results.append(chunk_output)
    # Combine chunk outputs (simplified; assumes concatenation of text output)
    return "\n".join(results)

def parse_output(output, min_peptide_length=10):
    """Parse IEDB output to extract predicted peptides and epitopes, filtering peptides by minimum length."""
    peptides = []
    epitopes = []
    current_section = None

    print_status("Parsing IEDB output", "info")
    # Process output line by line to reduce memory usage
    for line in output.strip().split('\n'):
        line = line.strip()
        print_status(f"Parsing line: {line}", "info")
        if line.startswith("Predicted peptides"):
            current_section = "peptides"
            print_status("Entered peptides section", "info")
            continue
        elif line.startswith("Position"):
            current_section = "epitopes"
            print_status("Entered epitopes section", "info")
            continue
        elif not line or line.startswith("No"):
            print_status(f"Skipping line: {line}", "info")
            continue

        if current_section == "peptides":
            parts = line.split()
            if len(parts) == 5:
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    peptide = parts[2]
                    length = float(parts[3])
                    score = float(parts[4])
                    if length >= min_peptide_length:
                        peptides.append((start, end, peptide, length, score))
                        print_status(f"Added peptide: {peptide} (length: {length})", "info")
                    else:
                        print_status(f"Skipping peptide {peptide} (length {length} < {min_peptide_length})", "info")
                except ValueError:
                    print_status(f"Skipping peptide line due to parsing error: {line}", "warning")
                    continue

        elif current_section == "epitopes":
            parts = line.split()
            if len(parts) == 6:
                try:
                    position = int(parts[0])
                    residue = parts[1]
                    start = int(parts[2])
                    end = int(parts[3])
                    peptide = parts[4]
                    score = float(parts[5])
                    if (end - start + 1) >= min_peptide_length:
                        epitopes.append((position, residue, start, end, peptide, score))
                        print_status(f"Added epitope: {peptide} (length: {end - start + 1})", "info")
                    else:
                        print_status(f"Skipping epitope {peptide} (length {end - start + 1} < {min_peptide_length})", "info")
                except ValueError:
                    print_status(f"Skipping epitope line due to parsing error: {line}", "warning")
                    continue

    print_status(f"Parsed {len(peptides)} peptides and {len(epitopes)} epitopes", "info")
    return peptides, epitopes

def run_bcell(input_fasta, output_dir=".", output_file="bcell_out.csv", min_peptide_length=10):
    """
    Runs the IEDB B-cell epitope prediction tool, optimized for long sequences.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory where output should be saved (default: current directory).
        output_file (str): Name of the output CSV file (default: 'bcell_out.csv').
        min_peptide_length (int): Minimum peptide length to include in output (default: 10).

    Returns:
        int: Exit status code (0 if successful).
    """
    # Hardcode method and window size
    method = "Emini"
    window_size = 10

    # Use output_file path directly if it includes a directory, otherwise combine with output_dir
    output_path = output_file if os.path.dirname(output_file) else os.path.join(output_dir, output_file)
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    print_status(f"Running B-cell prediction on: {input_fasta}", "info")
    print_status(f"Output will be saved to: {output_path}", "info")
    print_status(f"Using method: {method}", "info")
    print_status(f"Minimum peptide length: {min_peptide_length}", "info")
    print_status(f"Using window size: {window_size}", "info")

    # Check system resources
    check_system_resources()

    try:
        # Read sequence from FASTA
        sequence = read_fasta(input_fasta)
        # Run IEDB prediction
        epitopes_output = query_iedb(sequence, method, window_size=window_size)
        if not epitopes_output.strip():
            print_status("No output from IEDB tool. CSV will be empty.", "warning")
        # Parse output with minimum peptide length filter
        peptides, epitopes = parse_output(epitopes_output, min_peptide_length)
        # Save to CSV incrementally
        with open(output_path, "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Type", "Position", "Residue", "Start", "End", "Peptide", "Length", "Score"])
            for start, end, peptide, length, score in peptides:
                writer.writerow(["Peptide", "", "", start, end, peptide, length, score])
                csvfile.flush()  # Write immediately to reduce memory usage
                print_status(f"Wrote peptide to CSV: {peptide}", "info")
            for position, residue, start, end, peptide, score in epitopes:
                writer.writerow(["Epitope", position, residue, start, end, peptide, "", score])
                csvfile.flush()
                print_status(f"Wrote epitope to CSV: {peptide}", "info")
        if not peptides and not epitopes:
            print_status("No peptides or epitopes meet the criteria. CSV may be empty except for header.", "warning")
        else:
            print_status(f"B-cell prediction completed. Results saved to: {output_path}", "success")
        return 0
    except Exception as e:
        print_status(f"Error running B-cell prediction: {e}", "error")
        raise

def main():
    """Command-line interface for B-cell epitope predictions."""
    parser = argparse.ArgumentParser(description="Run IEDB B-cell epitope prediction on a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", default="bcell_out.csv", help="Output CSV file (full path or filename)")
    parser.add_argument("-d", "--output_dir", default=".", help="Output directory (used if output is filename only)")
    parser.add_argument("-l", "--min_peptide_length", type=int, default=10, help="Minimum peptide length (default: 10)")

    args = parser.parse_args()
    run_bcell(args.input, args.output_dir, args.output, args.min_peptide_length)

if __name__ == "__main__":
    main()
