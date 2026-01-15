import os
import subprocess
import argparse
import csv
import glob
from pathlib import Path
import shutil
import re

# ... [Keep Search Paths and setup functions exactly the same] ...
ROOT_DIR = Path(__file__).resolve().parent.parent
PSORTB_SEARCH_PATHS = [
    os.path.join(ROOT_DIR, "tools/psortb/psortb"),
    os.path.join(ROOT_DIR, "vax_elan/tools/psortb/psortb"),
    shutil.which("psortb")
]
BIOPERL_DIRS = [
    os.path.join(ROOT_DIR, "tools/psortb/bioperl-live"),
    os.path.join(ROOT_DIR, "tools/psortb/bioperl-run")
]

def print_status(msg, status="info"):
    colors = {"info": "\033[94m", "success": "\033[92m", "warning": "\033[93m", "error": "\033[91m"}
    print(f"{colors.get(status, '')}{msg}\033[0m")

def find_psortb_executable():
    for path in PSORTB_SEARCH_PATHS:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None

PSORTB_EXECUTABLE = find_psortb_executable()

def setup_bioperl():
    perl5lib = os.environ.get("PERL5LIB", "")
    for bioperl_dir in BIOPERL_DIRS:
        if os.path.isdir(bioperl_dir):
            perl5lib = f"{bioperl_dir}:{perl5lib}"
    os.environ["PERL5LIB"] = perl5lib

def check_blast():
    return shutil.which("blastp") is not None

def read_fasta(fasta_path):
    sequences = {}
    current_id = None
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            elif current_id:
                sequences[current_id].append(line)
    return {seq_id: ''.join(seq) for seq_id, seq in sequences.items()}

# =========================================================================
# UPDATED: Changed "-o long" to "-o terse" for compatibility with CSV parsing
# =========================================================================
def run_psortb_tool(fasta_file, output_dir, organism_type="n"):
    if not os.path.isfile(fasta_file) or not PSORTB_EXECUTABLE:
        raise FileNotFoundError("Input or Executable missing")
 
    if not check_blast():
        raise RuntimeError("BLAST+ is required for PSORTb")
    
    setup_bioperl()
    os.makedirs(output_dir, exist_ok=True)
   
    organism_flag = {"n": "-n", "p": "-p", "a": "-a"}.get(organism_type, "-n")
    
    command = [
        PSORTB_EXECUTABLE,
        "-i", str(fasta_file),
        "-r", str(output_dir),
        organism_flag,
        "-o", "terse"  # CHANGED FROM "long" TO "terse"
    ]
   
    try:
        result = subprocess.run(
            command,
            check=True,
            cwd=output_dir, # Run within the target folder
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print_status(f"PSORTb failed: {e.stderr.strip()}", "error")
        raise

def find_latest_output(output_dir):
    output_dir = Path(output_dir)
    candidates = list(output_dir.glob("*psortb*.txt"))
    if not candidates:
        candidates = [p for p in output_dir.glob("*.txt") if "psortb" in p.name.lower()]
    if not candidates:
        raise FileNotFoundError(f"No PSORTb output in {output_dir}")
    return str(max(candidates, key=os.path.getctime))

# =========================================================================
# UPDATED: Robust parsing for the tab-separated "terse" output
# =========================================================================
def parse_psortb_output(output_file):
    try:
        with open(output_file, 'r') as f:
            # Skip empty lines or comment lines
            lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
            if not lines:
                return [], []
            
            header = lines[0].split("\t")
            predictions = []
            for line in lines[1:]:
                data = line.split("\t")
                data += ["N/A"] * (len(header) - len(data))
                predictions.append(dict(zip(header, data)))
            return header, predictions
    except Exception as e:
        print_status(f"Error parsing PSORTb output: {e}", "error")
        raise

# =========================================================================
# UPDATED: Logic to handle path splitting correctly
# =========================================================================
def run_psortb(input_fasta, output_dir=".", output_file="psortb_out.csv", organism_type="n"):
    print_status("Starting PSORTb analysis", "info")
    
    # If the output_file is passed as a full path, split it
    if os.path.isabs(output_file) or "/" in output_file or "\\" in output_file:
        target_path = Path(output_file)
        output_dir = str(target_path.parent)
        output_file = target_path.name
    
    os.makedirs(output_dir, exist_ok=True)
    output_csv = os.path.join(output_dir, output_file)
   
    try:
        read_fasta(input_fasta) 
        run_psortb_tool(input_fasta, output_dir, organism_type)
        
        raw_file = find_latest_output(output_dir)
        header, predictions = parse_psortb_output(raw_file)
       
        if header:
            with open(output_csv, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=header)
                writer.writeheader()
                writer.writerows(predictions)
            print_status(f"PSORTb completed → {output_csv}", "success")
        else:
            print_status("PSORTb finished but no data was parsed.", "warning")
            
        return 0
    except Exception as e:
        print_status(f"PSORTb failed: {e}", "error")
        return 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", default="psortb_out.csv")
    parser.add_argument("-d", "--output_dir", default=".")
    parser.add_argument("-t", "--type", default="n", choices=["n", "p", "a"])
    args = parser.parse_args()
    run_psortb(args.input, args.output_dir, args.output, args.type)

if __name__ == "__main__":
    main()
