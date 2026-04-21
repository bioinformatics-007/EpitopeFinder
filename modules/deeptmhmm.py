import os
import subprocess
import argparse
import shutil
from pathlib import Path
import traceback
import csv
from Bio import SeqIO



# Define root directory
ROOT_DIR = Path(__file__).resolve().parent.parent  # Points to Vaxelan_2_0/vax_elan/

# Define default biolib path
BIOLIB_PATH = os.environ.get("BIOLIB_PATH", "/home/yuktika/anaconda3/envs/vaxelan_new/bin/biolib")

def print_status(msg, status="info"):
    """Print colored status messages."""
    colors = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m",
        "debug": "\033[95m"
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")
    
def validate_paths(input_fasta, output_dir, biolib_path=None):
    """Validate input paths, FASTA content, and biolib executable."""
    input_fasta = Path(input_fasta).resolve()
    output_dir = Path(output_dir).resolve()
    if biolib_path:
        biolib_path = Path(biolib_path).resolve()

    # Check input FASTA
    if not input_fasta.is_file():
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    if not os.access(input_fasta, os.R_OK):
        raise PermissionError(f"No read permission for FASTA file: {input_fasta}")

    # Check biolib executable if provided
    if biolib_path and not biolib_path.is_file():
        raise FileNotFoundError(f"Biolib executable not found: {biolib_path}")
    if biolib_path and not os.access(biolib_path, os.X_OK):
        raise PermissionError(f"No execute permission for biolib: {biolib_path}")

    # Check output directory
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise PermissionError(f"Cannot create output directory {output_dir}: {e}")
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"No write permission for output directory: {output_dir}")

    # Validate FASTA content
    try:
        sequences = list(SeqIO.parse(input_fasta, "fasta"))
        if not sequences:
            raise ValueError(f"Input FASTA file is empty or invalid: {input_fasta}")
        print_status(f"Found {len(sequences)} sequences in {input_fasta}", "info")
        seq_ids = [seq.id for seq in sequences]
        if len(seq_ids) != len(set(seq_ids)):
            duplicates = [sid for sid in set(seq_ids) if seq_ids.count(sid) > 1]
            raise ValueError(f"Duplicate sequence IDs found in {input_fasta}: {duplicates}")
        for seq in sequences:
            if not seq.seq:
                raise ValueError(f"Sequence {seq.id} is empty in {input_fasta}")
            if not all(c in "ACDEFGHIKLMNPQRSTVWY" for c in seq.seq.upper()):
                print_status(f"Warning: Sequence {seq.id} contains non-standard amino acids", "warning")
        return input_fasta, output_dir, sequences, biolib_path
    except Exception as e:
        raise ValueError(f"Invalid FASTA format in {input_fasta}: {e}")

def create_gff_file(gff_file):
    """Create TMRs.gff3 with placeholder content if it doesn't exist."""
    gff_file = Path(gff_file)
    gff_content = """# tr|P88147|P88147_9HIV2 Number of predicted TMRs: 1
tr|P88147|P88147_9HIV2	signal	1	22
tr|P88147|P88147_9HIV2	outside	23	686
tr|P88147|P88147_9HIV2	TMhelix	687	701
tr|P88147|P88147_9HIV2	inside	702	865
# sp|Q05322|VP24_EBOZM Length: 251
# sp|Q05322|VP24_EBOZM Number of predicted TMRs: 0
sp|Q05322|VP24_EBOZM	inside	1	251
"""
    try:
        gff_file.parent.mkdir(parents=True, exist_ok=True)
        with open(gff_file, "w", encoding="utf-8") as f:
            f.write(gff_content)
        print_status(f"Created GFF file: {gff_file}", "success")
        return True
    except Exception as e:
        print_status(f"Failed to create GFF file {gff_file}: {e}", "error")
        return False

def parse_output_file(file_path, fasta_sequences):
    """
    Parse DeepTMHMM GFF-like output file (TMRs.gff3).

    Args:
        file_path: Path to TMRs.gff3.
        fasta_sequences: List of SeqRecord objects from input FASTA.

    Returns:
        tuple: (gff_data, header_data)
        - gff_data: list of dicts with 'Sequence_ID', 'Type', 'Start', 'End', 'TMR_Count'
        - header_data: list of dicts with 'Sequence_ID', 'Length', 'TMR_Count'
    """
    gff_data = []
    header_data = []
    file_path = Path(file_path)
    print_status(f"DEBUG: Parsing file {file_path}", "debug")

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            if not lines:
                print_status(f"DEBUG: File {file_path} is empty", "warning")
                return [], []

            print_status(f"DEBUG: File {file_path} contents:\n{' '.join(line.strip() for line in lines)}", "debug")

            seq_id_to_tmr_count = {seq.id: 0 for seq in fasta_sequences}
            seq_id_to_length = {seq.id: str(len(seq.seq)) for seq in fasta_sequences}
            current_seq_id = None
            length = None
            tmr_count = 0

            for line in lines:
                line = line.strip()
                print_status(f"DEBUG: Processing line: {line}", "debug")
                if not line or line.startswith("##") or line == "//":
                    continue
                if line.startswith("#"):
                    if "Length:" in line:
                        length = line.split("Length:")[1].strip()
                    elif "Number of predicted TMRs:" in line:
                        tmr_count = int(line.split("Number of predicted TMRs:")[1].strip())
                        if current_seq_id:
                            seq_id_to_tmr_count[current_seq_id] = tmr_count
                            header_data.append({
                                "Sequence_ID": current_seq_id,
                                "Length": length or seq_id_to_length.get(current_seq_id, "Unknown"),
                                "TMR_Count": str(tmr_count)
                            })
                            if tmr_count == 0:
                                gff_data.append({
                                    "Sequence_ID": current_seq_id,
                                    "Type": "None",
                                    "Start": "0",
                                    "End": "0",
                                    "TMR_Count": "0"
                                })
                    elif not line.startswith("##"):
                        current_seq_id = line[1:].split()[0]
                        if current_seq_id not in seq_id_to_length:
                            print_status(f"Warning: Sequence ID {current_seq_id} not found in input FASTA", "warning")
                            current_seq_id = None
                else:
                    fields = line.split("\t")
                    if len(fields) >= 4 and current_seq_id:
                        seq_id = fields[0]
                        feature_type = fields[1]
                        start = fields[2]
                        end = fields[3]
                        if seq_id == current_seq_id and feature_type == "TMhelix":
                            gff_data.append({
                                "Sequence_ID": seq_id,
                                "Type": feature_type,
                                "Start": start,
                                "End": end,
                                "TMR_Count": str(seq_id_to_tmr_count.get(seq_id, 0))
                            })

            # Ensure all FASTA sequences are represented
            for seq in fasta_sequences:
                if seq.id not in seq_id_to_tmr_count or seq.id not in {entry["Sequence_ID"] for entry in header_data}:
                    seq_id_to_tmr_count[seq.id] = 0
                    header_data.append({
                        "Sequence_ID": seq.id,
                        "Length": str(len(seq.seq)),
                        "TMR_Count": "0"
                    })
                    gff_data.append({
                        "Sequence_ID": seq.id,
                        "Type": "None",
                        "Start": "0",
                        "End": "0",
                        "TMR_Count": "0"
                    })

        if not gff_data and not header_data:
            print_status(f"No valid data parsed from {file_path}", "warning")
        return gff_data, header_data

    except Exception as e:
        print_status(f"Failed to parse output file {file_path}: {e}", "warning")
        print_status(traceback.format_exc(), "debug")
        return [], []

def run_deeptmhmm(input_fasta, output_dir, biolib_path=BIOLIB_PATH):
    """
    Run DeepTMHMM via biolib and preserve only the biolib_results folder.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory to save biolib results.
        biolib_path (str): Path to biolib executable.

    Returns:
        int: Exit code 0 success, 1 failure.
    """
    try:
        # Validate input and environment
        input_fasta, output_dir, _, biolib_path = validate_paths(input_fasta, output_dir, biolib_path)

        # Run biolib command
        command = [
            str(biolib_path),
            "run",
            "DTU/DeepTMHMM",
            "--fasta", str(input_fasta)
        ]
        print_status(f"Running DeepTMHMM: {' '.join(command)}", "info")
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        print_status(f"STDOUT: {result.stdout}", "info")
        if result.stderr:
            print_status(f"STDERR: {result.stderr}", "warning")

        # Move biolib_results to output directory
        biolib_results_dir = Path.cwd() / "biolib_results"
        if biolib_results_dir.exists():
            target_dir = output_dir / "biolib_results"
            shutil.move(str(biolib_results_dir), str(target_dir))
            print_status(f"Results moved to: {target_dir}", "success")
        else:
            print_status(f"'biolib_results' directory not found in {Path.cwd()}", "error")
            return 1

        return 0

    except subprocess.CalledProcessError as e:
        print_status(f"Error running DeepTMHMM: {e}", "error")
        print_status(f"Command: {' '.join(e.cmd)}", "error")
        print_status(f"STDOUT: {e.stdout}", "error")
        print_status(f"STDERR: {e.stderr}", "error")
        return 1
    except Exception as e:
        print_status(f"Unexpected error: {e}", "error")
        print_status(f"Traceback: {traceback.format_exc()}", "error")
        return 1


def main():
    """Command-line interface for DeepTMHMM prediction and CSV generation."""
    parser = argparse.ArgumentParser(description="Run DeepTMHMM using biolib and generate combined CSV output.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to the input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Directory to save results")
    parser.add_argument("--gff_file", help="Path to TMRs.gff3 file (default: output_dir/deeptmhmm_results/TMRs.gff3)")
    parser.add_argument("--output_file", help="Filename for combined CSV result")
    parser.add_argument("--biolib", default=BIOLIB_PATH, help="Path to biolib executable")

    args = parser.parse_args()
    return run_deeptmhmm(args.fasta, args.output, args.gff_file, args.output_file, args.biolib)

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
