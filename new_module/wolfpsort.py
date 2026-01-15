import os
import argparse
import subprocess
from pathlib import Path

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

def validate_file(file_path, description):
    """Validate that a file exists and is readable."""
    file_path = Path(file_path)
    if not file_path.is_file() or not os.access(file_path, os.R_OK):
        raise FileNotFoundError(f"{description} not found or not readable: {file_path}")
    return file_path

def validate_directory(dir_path, description):
    """Validate that a directory exists or can be created."""
    dir_path = Path(dir_path).resolve()
    try:
        dir_path.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise OSError(f"Cannot create {description}: {dir_path} ({e})")
    return dir_path

def run_wolf_psort(input_fasta, output_file, wolf_psort_path=ROOT_DIR / "tools" / "WoLFPSort-master",
                   organism_type="fungi", html_output=False):
    """
    Run WoLF PSORT command and save results to a file and individual files by UniProt IDs.

    Args:
        input_fasta (str or Path): Path to input FASTA file.
        output_file (str or Path): Path to output text file.
        wolf_psort_path (str or Path): Path to WoLF PSORT installation directory.
        organism_type (str): Organism type ("fungi", "plant", "animal").
        html_output (bool): If True, generate HTML output instead of summary.

    Returns:
        int: Exit status code.
    """
    # Validate inputs
    wolf_psort_path = Path(wolf_psort_path).resolve()
    bin_dir = wolf_psort_path / "bin"
    if not bin_dir.is_dir():
        raise NotADirectoryError(f"WoLF PSORT bin directory not found: {bin_dir}")

    script_name = "runWolfPsortHtmlTables" if html_output else "runWolfPsortSummary"
    script_path = bin_dir / script_name
    validate_file(script_path, f"WoLF PSORT script ({script_name})")

    valid_organisms = {"fungi", "plant", "animal"}
    if organism_type not in valid_organisms:
        raise ValueError(f"Organism type must be one of {valid_organisms}, got: {organism_type}")

    input_fasta = validate_file(input_fasta, "Input FASTA file")
    output_path = Path(output_file)
    output_dir = validate_directory(output_path.parent, "Output directory")

    print_status(f"Running WoLF PSORT analysis on: {input_fasta}", "info")
    print_status(f"Output will be saved to: {output_path}", "info")

    original_dir = Path.cwd()
    try:
        os.chdir(bin_dir)

        if html_output:
            base_name = input_fasta.stem
            command = [f"./{script_name}", organism_type, str(output_dir), base_name]
            with open(input_fasta, "r") as fasta_input:
                process = subprocess.run(
                    command,
                    stdin=fasta_input,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True)
        else:
            command = [f"./{script_name}", organism_type]
            with open(input_fasta, "r") as fasta_input, open(output_path, "w") as out_file:
                process = subprocess.run(
                    command,
                    stdin=fasta_input,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True)
                out_file.write(process.stdout)
                print(process.stdout)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(
                process.returncode,
                command,
                output=process.stdout,
                stderr=process.stderr)

        if process.stderr:
            print_status(f"Warnings/Errors: {process.stderr}", "warning")

        print_status("WoLF PSORT analysis completed.", "success")

        # Parse output and create individual UniProt ID files + combined file
        output_lines = process.stdout.strip().split("\n")

        # Try to detect header line (contains "UniProt" or "ID")
        header = None
        data_lines = []
        for line in output_lines:
            if not line.strip():
                continue
            if header is None and ("UniProt" in line or "ID" in line):
                header = line.strip()
            else:
                data_lines.append(line.strip())

        combined_file_path = output_dir / "combined_results.txt"
        with open(combined_file_path, "w") as combined_file:
            if header:
                combined_file.write(header + "\n")
            combined_file.write("\n".join(data_lines))

        for line in data_lines:
            columns = line.split()
            if not columns:
                continue
            uniprot_id = columns[0]
            individual_file = output_dir / f"{uniprot_id}.txt"
            with open(individual_file, "w") as indfile:
                if header:
                    indfile.write(header + "\n")
                indfile.write(line + "\n")

        print_status(f"Created individual UniProt ID files and combined_results.txt in {output_dir}", "success")

        return 0

    except Exception as e:
        print_status(f"Error during WoLF PSORT analysis: {e}", "error")
        return 1

    finally:
        os.chdir(original_dir)

def main():
    parser = argparse.ArgumentParser(
        description="Wrapper script for running WoLF PSORT protein localization prediction."
    )
    parser.add_argument(
        "--wolf-path",
        default=ROOT_DIR / "tools" / "WoLFPSort-master",
        help="Path to WoLF PSORT installation directory (default: %(default)s)"
    )
    parser.add_argument(
        "--organism",
        default="fungi",
        help="Organism type for prediction (fungi, plant, animal; default: %(default)s)"
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to input FASTA file"
    )
    parser.add_argument(
        "--output-dir",
        default="wolf_psort_output",
        help="Directory for output files (default: %(default)s)"
    )
    parser.add_argument(
        "--output-file",
        default="wolf_psort_out.txt",
        help="Output text file (default: %(default)s)"
    )
    parser.add_argument(
        "--html",
        action="store_true",
        help="Generate HTML output instead of summary"
    )

    args = parser.parse_args()

    # Construct output file path from directory and filename
    output_file_path = Path(args.output_dir) / args.output_file

    exit_code = run_wolf_psort(
        input_fasta=args.fasta,
        output_file=output_file_path,
        wolf_psort_path=args.wolf_path,
        organism_type=args.organism,
        html_output=args.html
    )
    return exit_code

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
