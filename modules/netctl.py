import subprocess
import os
import sys
import argparse
from pathlib import Path

# Dynamically Resolve root directory
ROOT_DIR = Path(__file__).resolve().parent.parent

# Path to NetCTL executable (relative to project root)
NETCTL_SCRIPT_PATH = os.path.join(ROOT_DIR, "tools/netCTL-1.2b/netCTL")

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

def run_netctl(input_fasta, output_dir=".", output_file="netctl_out.txt", netctl_script_path=NETCTL_SCRIPT_PATH, supertype="A2"):
    """
    Runs the NetCTL prediction tool.

    Args:
        input_fasta (str): Path to the input FASTA file or 'stdin' for standard input.
        output_dir (str): Directory where output should be saved (default: current directory).
        output_file (str): Name of the output file or 'stdout' for console output (default: 'netctl_out.txt').
        netctl_script_path (str): Path to the NetCTL executable script (default: NETCTL_SCRIPT_PATH).
        supertype (str): MHC supertype for prediction (default: 'A2').

    Returns:
        int: Exit status code (0 if successful).
    """
    # Convert paths to strings to avoid PosixPath issues
    input_fasta = str(input_fasta)
    output_dir = str(output_dir)
    netctl_script_path = str(netctl_script_path)

    # Validate inputs
    if input_fasta != 'stdin' and not os.path.isfile(input_fasta):
        print_status(f"Input FASTA file not found: {input_fasta}", "error")
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    if not os.path.isfile(netctl_script_path):
        print_status(f"NetCTL script not found: {netctl_script_path}", "error")
        raise FileNotFoundError(f"NetCTL script not found: {netctl_script_path}")

    # Debug: Print NETCTLPATH and binary directory
    netctl_base_dir = os.path.dirname(netctl_script_path)
    binary_dir = os.path.join(netctl_base_dir, "Linux_x86_64")
    print_status(f"NETCTLPATH: {os.environ.get('NETCTLPATH', 'Not set')}", "info")
    print_status(f"Checking binaries in: {binary_dir}", "info")
    if os.path.isdir(binary_dir):
        print_status(f"Binary directory contents: {os.listdir(binary_dir)}", "info")
    else:
        print_status(f"Binary directory not found: {binary_dir}", "error")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)

    # Sanitize FASTA input to temporary file to avoid header parsing issues in NetCTL binary
    import tempfile, csv
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    with open(input_fasta, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                header_id = line[1:].strip().split()[0].replace('|', '_').replace('=', '_')[:30]
                temp_fasta.write(f">{header_id}\n")
            else:
                temp_fasta.write(line)
    temp_fasta.close()
    clean_input_fasta = temp_fasta.name

    # Set CTLHOME env variable so tcsh script can find its binaries relative to itself
    env = os.environ.copy()
    env["CTLHOME"] = netctl_base_dir

    try:
        result = subprocess.run(
            ["tcsh", netctl_script_path, "-f", clean_input_fasta, "-s", supertype],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            stdin=subprocess.DEVNULL,
            env=env
        )

        stdout_text = result.stdout
        # Parse NetCTL output lines into structured CSV format
        rows = []
        for line in stdout_text.splitlines():
            line_str = line.strip()
            if "ID " in line_str and " pep " in line_str and " COMB " in line_str:
                is_epitope = "<-E" in line_str
                parts = line_str.replace("<-E", "").split()
                try:
                    seq_id = parts[2]
                    peptide = parts[4]
                    aff = parts[6]
                    aff_rescale = parts[8]
                    cle = parts[10]
                    tap = parts[12]
                    comb = parts[14]
                    rows.append({
                        "SeqID": seq_id,
                        "Peptide": peptide,
                        "MHC_Affinity": aff,
                        "Rescaled_Affinity": aff_rescale,
                        "Cleavage_Score": cle,
                        "TAP_Score": tap,
                        "Combined_Score": comb,
                        "Is_Epitope": "YES" if is_epitope else "NO"
                    })
                except Exception:
                    continue

        if output_file == 'stdout':
            if rows:
                writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)
            else:
                sys.stdout.write(stdout_text)
        else:
            with open(output_path, "w", newline="") as out:
                if rows:
                    writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()))
                    writer.writeheader()
                    writer.writerows(rows)
                else:
                    out.write(stdout_text)

        print_status("NetCTL ran successfully.", "success")
        return result.returncode
    except subprocess.CalledProcessError as e:
        print_status(f"NetCTL error: {e.stderr}", "error")
        raise RuntimeError(f"NetCTL execution failed: {e.stderr}")
    except Exception as e:
        print_status(f"Unexpected error: {e}", "error")
        raise
    finally:
        if os.path.exists(clean_input_fasta):
            os.remove(clean_input_fasta)

def main():
    """Command-line interface for NetCTL predictions."""
    parser = argparse.ArgumentParser(description="Run NetCTL on a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file or 'stdin' for standard input")
    parser.add_argument("-o", "--output", default="netctl_out.txt", help="Output result file or 'stdout' for console output")
    parser.add_argument("-d", "--output_dir", default=".", help="Output directory")
    parser.add_argument("-s", "--script", default=NETCTL_SCRIPT_PATH, help="Path to NetCTL script (default: derived from ROOT_DIR)")
    parser.add_argument("--supertype", default="A2", help="MHC supertype (e.g., A2, A1, B7)")

    args = parser.parse_args()
    print_status(f"NETCTL_SCRIPT_PATH: {NETCTL_SCRIPT_PATH}", "info")
    run_netctl(args.input, args.output_dir, args.output, args.script, args.supertype)

if __name__ == "__main__":
    main()
