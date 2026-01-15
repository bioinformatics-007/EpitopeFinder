import os
import sys
import json
import subprocess
import argparse
import pandas as pd
from pathlib import Path
from io import StringIO
import logging
import time
import tempfile

# Set up logging
logger = logging.getLogger('VaxElan')
logging.basicConfig(level=logging.INFO)

# Dynamically resolve root directory
current_path = Path(__file__).resolve()
while current_path.name != "Vaxelan_2_0" and current_path != current_path.parent:
    current_path = current_path.parent
root_dir = current_path

MHC_I_TOOL = os.path.join(root_dir, "tools/IEDB_NG_TC1-0.1.2-beta/src/tcell_mhci.py")

METHOD_CODES = {
    "a": "ann",
    "b": "comblib_sidney2008",
    "c": "consensus",
    "d": "netmhcpan_ba",
    "e": "netmhcpan_ba",
    "f": "netmhcpan_el",
    "g": "netmhcpan_el",
    "h": "netmhcpan_ba",
    "i": "smm",
    "j": "smmpmbec"
}

SUPPORTED_ALLELES = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01",
    "HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01",
    "HLA-B*53:01", "HLA-B*57:01", "HLA-B*58:01"
]

def create_iedb_json(input_fasta, method, alleles):
    data = {
        "input_sequence_text_file_path": os.path.abspath(input_fasta),
        "peptide_length_range": [8, 11],
        "alleles": ",".join(alleles),
        "predictors": [{"type": "binding", "method": method}]
    }
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".json", mode='w')
    json.dump(data, tmp)
    tmp.close()
    return tmp.name

def run_mhc1(method_code, input_file, output_file="mhci_out.csv", score_threshold=5000):
    start_time = time.time()
    method = METHOD_CODES.get(method_code)
    
    if not method:
        logger.error(f"Invalid method code: {method_code}")
        return 1

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    json_path = create_iedb_json(input_file, method, SUPPORTED_ALLELES)
    
    env = os.environ.copy()
    env["PYTHONWARNINGS"] = "ignore"
    
    command = [sys.executable, MHC_I_TOOL, "-j", json_path, "-f", "tsv"]
    
    logger.info(f"Running {method} via IEDB Next-Gen Tools...")
    
    try:
        result = subprocess.run(command, capture_output=True, text=True, env=env)
        
        # FIX 1: If the tool fails, return exit code 1 immediately
        if result.returncode != 0:
            logger.error(f"IEDB Tool Error: {result.stderr}")
            return 1

        if not result.stdout.strip():
            # FIX 2: Ensure column names match 'peptide' (lowercase)
            pd.DataFrame(columns=['Method', 'allele', 'peptide', 'score']).to_csv(output_path, index=False)
            logger.warning("No data returned. Created empty output file.")
            return 0

        df = pd.read_csv(StringIO(result.stdout), sep='\t')
        
        if 'ic50' in df.columns:
            df['score'] = pd.to_numeric(df['ic50'], errors='coerce')
        
        df = df[df['score'] <= score_threshold].copy()
        
        # FIX 3: DO NOT rename 'peptide' to 'Peptide'. 
        # main.py expects lowercase 'peptide'.
        df.insert(0, "Method", method)

        df.to_csv(output_path, index=False)
        
        elapsed = time.time() - start_time
        logger.info(f"✔ Completed in {elapsed:.1f}s. Saved to {output_file}")
        return 0

    except Exception as e:
        logger.error(f"Unexpected Script Error: {str(e)}")
        return 1
    finally:
        if os.path.exists(json_path):
            os.remove(json_path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--method", choices=METHOD_CODES.keys(), required=True)
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", default="mhci_out.csv")
    parser.add_argument("-s", "--score_threshold", type=float, default=5000)
    args = parser.parse_args()
    
    exit_code = run_mhc1(args.method, args.input_file, args.output_file, args.score_threshold)
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
