import subprocess
import os
import sys
import argparse
from pathlib import Path

# Dynamically Resolve root directory
current_path = Path(__file__).resolve()
while current_path.name != "Vaxelan_2_0" and current_path != current_path.parent:
    current_path = current_path.parent
ROOT_DIR = current_path

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

    # Validate output path
    if output_file != 'stdout' and os.path.isdir(output_path):
        print_status(f"Output file path is a directory: {output_path}", "error")
        raise IsADirectoryError(f"Output file path is a directory: {output_path}")

    print_status(f"Running NetCTL on: {input_fasta}", "info")
    print_status(f"Using NetCTL script: {netctl_script_path}", "info")
    print_status(f"MHC supertype: {supertype}", "info")
    if output_file != 'stdout':
        print_status(f"Output will be saved to: {output_path}", "info")

    # Determine output destination
    out = sys.stdout if output_file == 'stdout' else open(output_path, "w")

    try:
        result = subprocess.run(
            ["tcsh", netctl_script_path, "-f", input_fasta, "-s", supertype],
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            stdin=subprocess.DEVNULL
        )
        print_status("NetCTL ran successfully.", "success")
        return result.returncode
    except subprocess.CalledProcessError as e:
        print_status(f"NetCTL error: {e.stderr}", "error")
        raise RuntimeError(f"NetCTL execution failed: {e.stderr}")
    except Exception as e:
        print_status(f"Unexpected error: {e}", "error")
        raise
    finally:
        if output_file != 'stdout':
            out.close()

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
