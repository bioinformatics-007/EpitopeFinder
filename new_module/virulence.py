import os
import subprocess
import argparse
from pathlib import Path
import traceback
import shutil
import warnings
from Bio import SeqIO
import tempfile
import pandas as pd

# Suppress warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="scipy")
warnings.filterwarnings("ignore", category=UserWarning, module="autogluon")           # pkg_resources deprecation
warnings.filterwarnings("ignore", category=UserWarning, module="pkg_resources")      # extra safety

# Root directory
ROOT_DIR = Path(__file__).resolve().parent.parent

# Paths
VIRULENTPRED_SCRIPT = os.environ.get(
    "VIRULENTPRED_SCRIPT",
    ROOT_DIR / "tools" / "virulentpred_2_0" / "predict.pl"
)
WORKING_DIR = os.environ.get(
    "VIRULENTPRED_WORKING_DIR",
    ROOT_DIR / "tools" / "virulentpred_2_0"
)

# Conda environment where Perl + VirulentPred dependencies live
VIRULENTPRED_ENV = "autogluon_env"


def print_status(msg, status="info"):
    colors = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")


def validate_paths(input_fasta, output_dir, virulentpred_script, working_dir):
    input_fasta = Path(input_fasta).resolve()
    virulentpred_script = Path(virulentpred_script).resolve()
    working_dir = Path(working_dir).resolve()

    output_dir = Path(output_dir).resolve()
    if output_dir.suffix and output_dir.suffix != ".":
        output_dir = output_dir.parent
    output_dir = output_dir / "virulence"

    if not input_fasta.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")
    if not virulentpred_script.is_file():
        raise FileNotFoundError(f"VirulentPred script not found: {virulentpred_script}")
    if not working_dir.is_dir():
        raise FileNotFoundError(f"Working directory not found: {working_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    return input_fasta, output_dir, virulentpred_script, working_dir


def run_virulence(
    input_fasta,
    output_dir,
    virulentpred_script=str(VIRULENTPRED_SCRIPT),
    working_dir=str(WORKING_DIR),
    batch_idx=1
):
    try:
        input_fasta, output_dir, virulentpred_script, working_dir = validate_paths(
            input_fasta, output_dir, virulentpred_script, working_dir
        )

        sequences = []
        with open(input_fasta, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_id = record.id.split("|")[1] if "|" in record.id else record.id
                sequences.append((seq_id, str(record.seq)))

        if not sequences:
            print_status(f"No sequences found in FASTA: {input_fasta}", "error")
            return 1

        overall_success = 0
        failures = 0
        original_dir = os.getcwd()
        os.chdir(working_dir)

        combined_results = []
        expected_columns = None

        for seq_id, sequence in sequences:
            print_status(f"Processing sequence: {seq_id}", "info")

            temp_fasta_path = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode="w", suffix=".fasta", delete=False, dir=working_dir
                ) as temp_fasta:
                    temp_fasta.write(f">{seq_id}\n{sequence}\n")
                    temp_fasta_path = temp_fasta.name

                output_text_file = output_dir / f"{seq_id}.txt"

                perl_cmd = f"perl {virulentpred_script} {temp_fasta_path}"

                cmd = [
                    "conda", "run",
                    "--no-capture-output",
                    "-n", VIRULENTPRED_ENV,
                    "bash", "-c", perl_cmd
                ]

                print_status(f"Running in {VIRULENTPRED_ENV}: {perl_cmd}", "info")

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=300,  # 5 min timeout to prevent hangs
                    check=True
                )

                with open(output_text_file, "w") as f:
                    f.write(f"=== VirulentPred 2.0 Output for {seq_id} ===\n")
                    f.write(result.stdout)
                    if result.stderr:
                        f.write("\n=== Errors/Warnings ===\n")
                        f.write(result.stderr)

                output_files = move_output_files(working_dir, output_dir, seq_id)

                csv_file = next((f for f in output_files if f.suffix == ".csv"), None)
                if csv_file and csv_file.exists():
                    try:
                        df = pd.read_csv(csv_file)
                        df["Sequence_ID"] = seq_id
                        cols = ["Sequence_ID"] + [c for c in df.columns if c != "Sequence_ID"]
                        df = df[cols]
                        combined_results.append(df)
                        if expected_columns is None:
                            expected_columns = df.columns.tolist()
                    except Exception as e:
                        print_status(f"Failed to parse CSV for {seq_id}: {e}", "warning")

                verify_output_files(output_files, output_text_file)
                print_status(f"Results saved for {seq_id}", "success")
                overall_success += 1

            except subprocess.TimeoutExpired:
                print_status(f"Timeout processing {seq_id} (took >5 min)", "error")
                failures += 1
            except subprocess.CalledProcessError as e:
                print_status(f"VirulentPred failed for {seq_id} (code {e.returncode})", "error")
                if e.stderr:
                    print_status(e.stderr.strip(), "error")
                failures += 1
            except Exception as e:
                print_status(f"Unexpected error for {seq_id}: {e}", "error")
                print_status(traceback.format_exc(), "error")
                failures += 1
            finally:
                if temp_fasta_path and os.path.exists(temp_fasta_path):
                    try:
                        os.remove(temp_fasta_path)
                    except:
                        pass

        # Save combined results
        combined_file = output_dir / f"combined_virulence_batch_{batch_idx}.csv"
        try:
            if combined_results:
                combined_df = pd.concat(combined_results, ignore_index=True)
                if expected_columns:
                    for col in expected_columns:
                        if col not in combined_df.columns:
                            combined_df[col] = pd.NA
                    combined_df = combined_df[expected_columns]
                combined_df.to_csv(combined_file, index=False)
                print_status(f"Combined output saved: {combined_file} ({len(combined_df)} sequences)", "success")
            else:
                combined_file.write_text("No valid predictions generated\n")
                print_status(f"Combined output saved (empty): {combined_file}", "warning")
        except Exception as e:
            print_status(f"Failed to save combined output: {e}", "error")

        # Final summary
        print_status(
            f"VirulentPred completed: {overall_success}/{len(sequences)} sequences successful, {failures} failed",
            "success" if failures == 0 else "warning"
        )

        return 0 if failures == 0 else 1

    except Exception as e:
        print_status(f"Critical pipeline error: {e}", "error")
        print_status(traceback.format_exc(), "error")
        return 1

    finally:
        if 'original_dir' in locals():
            os.chdir(original_dir)


def move_output_files(working_dir, output_dir, sequence_id, prefix="VP_2.0_pred_res"):
    output_files = []
    for ext in [".xls", ".xlsx", ".csv"]:
        files = list(Path(working_dir).glob(f"{prefix}*{ext}"))
        for src in files:
            dest = output_dir / f"{sequence_id}{ext}"
            try:
                shutil.move(src, dest)
                output_files.append(dest)
                print_status(f"Moved: {dest}", "success")
            except Exception as e:
                print_status(f"Failed to move {src}: {e}", "warning")
    return output_files


def verify_output_files(output_files, text_file):
    for p in output_files + [text_file]:
        if p.exists():
            size = p.stat().st_size
            msg = f"Output found: {p} ({size} bytes)"
            print_status(msg, "success" if size > 100 else "warning")  # flag very small files
        else:
            print_status(f"Output missing: {p}", "error")


def main():
    parser = argparse.ArgumentParser(description="Run VirulentPred 2.0 via conda env")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--script", default=str(VIRULENTPRED_SCRIPT), help="Path to predict.pl")
    parser.add_argument("--working-dir", default=str(WORKING_DIR), help="Working dir")
    parser.add_argument("--batch_idx", type=int, default=1, help="Batch index")
    args = parser.parse_args()

    return run_virulence(
        args.input,
        args.output_dir,
        args.script,
        args.working_dir,
        args.batch_idx
    )


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
