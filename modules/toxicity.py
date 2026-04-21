import os
import subprocess
import argparse
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import traceback
import platform
import tempfile
import shutil
import numpy as np


# Paths
# Dynamically Resolve root directory
current_path = Path(__file__).resolve()
while current_path.name != "Vaxelan_2_0" and current_path != current_path.parent:
    current_path = current_path.parent
ROOT_DIR = current_path

TOXINPRED2_SCRIPT = ROOT_DIR / "tools/toxinpred2/toxinpred2.py"
MODEL_FILE = ROOT_DIR / "tools/toxinpred2/RF_model"
BLASTDB = ROOT_DIR / "tools/Database/data"
MERCI = ROOT_DIR / "tools/progs/MERCI_motif_locator.pl"
MOTIFS = ROOT_DIR / "tools/Database/pos_motif.txt"

# Define BLAST binary paths for different operating systems
BLAST_BINARIES = {
    "Linux": "/../tools/ncbi-blast-2.16.0+-aarch64-linux/ncbi-blast-2.16.0+/bin/blastp",
    "Darwin": "/../blast_binaries/mac/blastp",
    "Windows": "/../blast_binaries/windows/blastp.exe",
}

def print_status(message, status="info"):
    color_map = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    print(f"{color_map.get(status, '')}{message}\033[0m")


def find_model_file(path):
    path = Path(path) if path else None
    if path and path.is_file():
        print_status(f"Selected model file: {path}", "info")
        return path
    alternatives = [
        ROOT_DIR / "tools/toxinpred2/RF_model",
        ROOT_DIR / "tools/RF_model",
        Path("/home/yuktika/Downloads/Vaxelan_2_0/tools/toxinpred2/RF_model")
    ]
    for alt in alternatives:
        if alt.is_file():
            print_status(f"Using alternative model file: {alt}", "info")
            return alt
    print_status("No model file found in specified or alternative paths", "error")
    return None

def get_blastp_path():
    """Get the appropriate BLASTP binary path based on the operating system."""
    nf_path = TOXINPRED2_SCRIPT.parent
    operating_system = platform.system()
    blastp_path = str(nf_path) + BLAST_BINARIES.get(operating_system)
    if not blastp_path:
        raise ValueError(f"Unsupported operating system: {operating_system}")
    if not os.path.exists(blastp_path):
        raise FileNotFoundError(f"BLASTP binary not found: {blastp_path}")
    print_status(f"BLASTP binary path: {blastp_path}", "info")
    return blastp_path

def validate_fasta_file(fasta_path):
    """Validate FASTA file and return sequence IDs."""
    try:
        sequences = list(SeqIO.parse(fasta_path, "fasta"))
        if not sequences:
            print_status(f"No valid sequences found in {fasta_path}", "error")
            return None
        seq_ids = [record.id for record in sequences]
        print_status(f"Validated {len(seq_ids)} sequences in {fasta_path}", "info")
        return seq_ids
    except Exception as e:
        print_status(f"Failed to parse FASTA file {fasta_path}: {e}", "error")
        return None

def reformat_output(output_file, fasta_file, output_dir, display=2):
    try:
        df = pd.read_csv(output_file)
        required = ['Subject', 'ML Score', 'Prediction']
        if not all(col in df.columns for col in required):
            raise ValueError(f"Missing required columns: {required}")

        # Check for borderline scores
        for _, row in df.iterrows():
            score = row['ML Score']
            if 0.45 <= score <= 0.55:
                print_status(
                    f"Borderline ML Score for {row['Subject']}: {score}. Prediction ({row['Prediction']}) may be sensitive to threshold.",
                    "warning"
                )

        seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
        df['Sequence'] = df['Subject'].map(seq_dict)
        df = df[['Subject', 'ML Score', 'Sequence', 'Prediction']]

        if display == 1:
            df = df[df['Prediction'] == 'Toxin']

        output_dir = Path(output_dir) / "toxinpred"
        output_dir.mkdir(parents=True, exist_ok=True)

        combined_output_path = output_dir / "combined.csv"
        df.to_csv(combined_output_path, index=False)
        print_status(f"Saved combined output to: {combined_output_path}", "success")

        for _, row in df.iterrows():
            safe_subject = str(row['Subject']).replace("|", "_").replace("/", "_")
            out_path = output_dir / f"{safe_subject}_toxinpred.csv"
            pd.DataFrame([row]).to_csv(out_path, index=False)
            print_status(f"Saved individual output: {out_path}", "success")

        return 0
    except Exception as e:
        print_status(f"Failed to reformat output: {e}", "error")
        return 1

def MERCI_Processor(merci_file, merci_processed, seq_names):
    """Process MERCI output to classify sequences."""
    hh = []
    check = '>'
    with open(merci_file) as f:
        l = []
        for line in f:
            if line.strip():
                l.append(line)
            if 'COVERAGE' in line:
                for item in l:
                    if item.startswith(check):
                        hh.append(item)
                l = []
    
    jj, qq, kk = [], [], []
    if not hh:
        ff = [w.replace('>', '') for w in seq_names]
        for a in ff:
            jj.append(a)
            qq.append(np.array(['0']))
            kk.append('Non-Toxin')
    else:
        ff = [w.replace('\n', '').replace('>', '') for w in hh]
        rr = [w.replace('>', '') for w in seq_names]
        oo = np.unique(ff + rr)
        df1 = pd.DataFrame([x.strip() for x in l[1:]], columns=['Name'])
        df1['Name'] = df1['Name'].str.strip('(')
        df1[['Seq', 'Hits']] = df1['Name'].str.split("(", expand=True)
        df1[['Seq', 'Hits']] = df1[['Seq', 'Hits']].replace(to_replace=r"[)\s]|motifs match", value='', regex=True)
        for j in oo:
            if j in df1['Seq'].values:
                jj.append(j)
                qq.append(df1.loc[df1['Seq'] == j, 'Hits'].values)
                kk.append('Toxin')
            else:
                jj.append(j)
                qq.append(np.array(['0']))
                kk.append('Non-Toxin')
    
    df3 = pd.DataFrame({'Name': jj, 'Hits': qq, 'Prediction': kk})
    if df3.empty:
        raise ValueError(f"MERCI output is empty: {merci_processed}")
    df3.to_csv(merci_processed, index=None)
    if not os.path.exists(merci_processed):
        raise RuntimeError(f"Failed to write MERCI output: {merci_processed}")
    if os.stat(merci_processed).st_size == 0:
        raise ValueError(f"MERCI output is empty: {merci_processed}")

def Merci_after_processing(merci_processed, final_merci):
    """Convert MERCI output to scores."""
    df = pd.read_csv(merci_processed)[['Name', 'Hits']]
    df['MERCI Score'] = df['Hits'].apply(lambda x: 0.5 if x > 0 else 0)
    df[['Name', 'MERCI Score']].to_csv(final_merci, index=None, header=['Subject', 'MERCI Score'])
    if not os.path.exists(final_merci):
        raise RuntimeError(f"Failed to write final MERCI output: {final_merci}")
    if os.stat(final_merci).st_size == 0:
        raise ValueError(f"Final MERCI output is empty: {final_merci}")

def BLAST_processor(blast_result, blast_processed, seq_names):
    """Process BLAST output to assign scores."""
    if os.stat(blast_result).st_size != 0:
        df1 = pd.read_csv(blast_result, sep="\t", header=None)
        df2 = df1.iloc[:, :2]
        df2.columns = ['Subject', 'Query']
        df3 = pd.DataFrame()
        for i in df2['Subject'].unique():
            df3 = pd.concat([df3, df2[df2['Subject'] == i].iloc[:5]]).reset_index(drop=True)
        df3['label'] = df3['Query'].apply(lambda x: x.split("_")[0])
        df3['vote'] = df3['label'].apply(lambda x: 1 if x == 'P' else -1)
        df4 = df3.groupby('Subject')['vote'].sum().reset_index()
        df4['BLAST Score'] = df4['vote'].apply(lambda x: 0.5 if x > 0 else (-0.5 if x < 0 else 0))
        df4 = df4[['Subject', 'BLAST Score']]
    else:
        df4 = pd.DataFrame({'Subject': [w.replace('>', '') for w in seq_names], 'BLAST Score': [0] * len(seq_names)})
    df4.to_csv(blast_processed, index=None)
    if not os.path.exists(blast_processed):
        raise RuntimeError(f"Failed to write BLAST output: {blast_processed}")
    if os.stat(blast_processed).st_size == 0:
        raise ValueError(f"BLAST output is empty: {blast_processed}")

def hybrid(ml_output, seq_names, merci_output, blast_output, threshold, final_output):
    """Combine ML, MERCI, and BLAST scores for hybrid prediction."""
    df_ml = pd.read_csv(ml_output, header=None, names=['ML Score'])
    df_seq = pd.DataFrame({'Subject': [w.replace('>', '') for w in seq_names]})
    df_ml = pd.concat([df_seq, df_ml], axis=1)
    df_merci = pd.read_csv(merci_output)
    df_blast = pd.read_csv(blast_output)
    df = pd.merge(df_ml, df_merci, how='outer', on='Subject')
    df = pd.merge(df, df_blast, how='outer', on='Subject').fillna(0)
    df['Hybrid Score'] = df[['ML Score', 'MERCI Score', 'BLAST Score']].sum(axis=1)
    df['Prediction'] = df['Hybrid Score'].apply(lambda x: 'Toxin' if x > float(threshold) else 'Non-Toxin')
    df = df.round(3)
    df.to_csv(final_output, index=None)
    if not os.path.exists(final_output):
        raise RuntimeError(f"Failed to write hybrid output: {final_output}")
    if os.stat(final_output).st_size == 0:
        raise ValueError(f"Hybrid output is empty: {final_output}")

def run_toxinpred(fasta_path, output_dir, output_file=None, model_file=None, threshold=0.6, model=1, display=2, batch_idx=1):
    try:
        fasta_path = Path(fasta_path)
        output_dir = Path(output_dir)
        model_file = find_model_file(model_file)
        if not model_file:
            raise FileNotFoundError(f"Model file not found: {model_file}")

        if not TOXINPRED2_SCRIPT.exists():
            raise FileNotFoundError(f"ToxinPred2 script missing: {TOXINPRED2_SCRIPT}")

        output_dir.mkdir(parents=True, exist_ok=True)

        if fasta_path.is_file():
            fasta_files = [fasta_path]
        elif fasta_path.is_dir():
            fasta_files = list(fasta_path.glob("*.fa")) + list(fasta_path.glob("*.fasta"))
            if not fasta_files:
                raise FileNotFoundError("No FASTA files found in directory.")
        else:
            raise FileNotFoundError(f"Invalid input path: {fasta_path}")

        # Validate FASTA files
        valid_fasta_files = []
        for fasta in fasta_files:
            seq_ids = validate_fasta_file(fasta)
            if seq_ids:
                valid_fasta_files.append((fasta, seq_ids))
            else:
                print_status(f"Skipping invalid FASTA file: {fasta}", "warning")
                continue

        if not valid_fasta_files:
            raise ValueError("No valid FASTA files to process.")

        # Create temporary directory
        temp_dir = tempfile.mkdtemp(prefix=f"toxinpred2_batch_{batch_idx}_")
        print_status(f"Created temporary directory: {temp_dir}", "info")

        for fasta, seq_ids in valid_fasta_files:
            output_csv = output_dir / "toxinpred" / f"temp_output_batch_{batch_idx}.csv"
            sequence_file = Path(temp_dir) / f"Sequence_{batch_idx}"
            with open(sequence_file, 'w') as f:
                for record in SeqIO.parse(fasta, "fasta"):
                    f.write(f">{record.id}\n{record.seq}\n")

            # Ensure ToxinPred2 script is executable
            if not os.access(TOXINPRED2_SCRIPT, os.X_OK):
                os.chmod(TOXINPRED2_SCRIPT, 0o755)
                print_status(f"Set executable permissions for {TOXINPRED2_SCRIPT}", "info")

            command = [
                "python3", str(TOXINPRED2_SCRIPT),
                "-i", str(sequence_file),
                "-m", str(model),
                "-o", str(output_csv),
                "-t", str(threshold),
                "-d", str(display)
            ]

            print_status(f"Executing command for batch {batch_idx}: {' '.join(command)}", "info")
            os.chdir(TOXINPRED2_SCRIPT.parent)
            try:
                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    check=True,  # Raises CalledProcessError on non-zero exit code
                    stdin=subprocess.DEVNULL
                )
                print_status(f"ToxinPred2 stdout for batch {batch_idx}: {result.stdout}", "info")
                if result.stderr:
                    print_status(f"ToxinPred2 stderr for batch {batch_idx}: {result.stderr}", "warning")
            except subprocess.CalledProcessError as e:
                print_status(f"ToxinPred2 failed for batch {batch_idx} with exit code {e.returncode}: {e.stderr}", "error")
                return 1

            if output_csv.exists() and output_csv.stat().st_size > 0:
                with open(output_csv, 'r') as f:
                    raw_content = f.read()
                print_status(f"Raw ToxinPred2 output for batch {batch_idx}:\n{raw_content}", "info")
                
                if model == 2:
                    # Process hybrid model output
                    try:
                        df = pd.read_csv(output_csv)
                        if 'Hybrid Score' in df.columns:
                            df['ML Score'] = df['Hybrid Score']  # Use Hybrid Score as ML Score for consistency
                        status = reformat_output(output_csv, fasta, output_dir, display)
                        if status != 0:
                            print_status(f"Failed to reformat output for batch {batch_idx}", "error")
                            return 1
                    except Exception as e:
                        print_status(f"Error processing hybrid output for batch {batch_idx}: {e}", "error")
                        return 1
                else:
                    status = reformat_output(output_csv, fasta, output_dir, display)
                    if status != 0:
                        print_status(f"Failed to reformat output for batch {batch_idx}", "error")
                        return 1
            else:
                print_status(f"Empty or missing output for batch {batch_idx}: {output_csv}", "error")
                return 1

        # Clean up temporary directory
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                print_status(f"Removed temporary directory: {temp_dir}", "info")
            except Exception as e:
                print_status(f"Failed to remove temporary directory {temp_dir}: {e}", "warning")

        print_status(f"toxinpred completed successfully for batch {batch_idx}", "success")
        return 0
    except Exception as e:
        print_status(f"Execution failed for batch {batch_idx}: {e}", "error")
        return 1

def main():
    parser = argparse.ArgumentParser(description="Run ToxinPred2 for protein toxicity prediction.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file or directory")
    parser.add_argument("-o", "--output", default=str(ROOT_DIR / "Results"), help="Output directory (default: Results)")
    parser.add_argument("-f", "--output_file", help="Optional: output CSV filename (single FASTA only)")
    parser.add_argument("-m", "--model_file", default=str(MODEL_FILE), help="Path to model file (default: tools/toxinpred2/RF_model)")
    parser.add_argument("-a", "--approach", type=int, choices=[1, 2], default=1, help="Model: 1: AAC-based RF, 2: Hybrid (default: 1)")
    parser.add_argument("-t", "--threshold", type=float, default=0.6, help="Threshold: Value between 0 to 1 (default: 0.6)")
    parser.add_argument("-d", "--display", type=int, choices=[1, 2], default=2, help="Display: 1: Toxin only, 2: All peptides (default: 2)")
    parser.add_argument("--batch_idx", type=int, default=1, help="Batch index (for naming output)")

    args = parser.parse_args()

    if args.output_file and not args.output:
        args.output = str(Path(args.output_file).parent)

    status = run_toxinpred(args.input, args.output, args.output_file, args.model_file, args.threshold, args.approach, args.display, args.batch_idx)
    exit(status)

if __name__ == "__main__":
    main()
