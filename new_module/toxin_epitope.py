import subprocess
import os
import shutil
import zipfile
import argparse
from pathlib import Path

def run_toxinpred3(input_file, output_file):
    # Determine the root directory two levels above this script
    ROOT_DIR = Path(__file__).resolve().parents[1]

    # Define paths properly using os.path.join
    toxinpred_script = os.path.join(ROOT_DIR, "tools", "toxinpred3", "toxinpred3.py")
    model_base_path = os.path.join(ROOT_DIR, "tools", "toxinpred3")
    model_path = os.path.join(model_base_path, "toxinpred3.0_model.pkl")
    nested_model_path = os.path.join(model_base_path, "model", "toxinpred3.0_model.pkl")
    zip_model_path = os.path.join(model_base_path, "model", "toxinpred3.0_model.pkl.zip")
    working_directory = os.path.join(ROOT_DIR, "tools", "toxinpred3")

    # Convert output_file to Path object and ensure it has a file name
    output_file = Path(output_file)
    if output_file.is_dir() or not output_file.suffix:
        output_file = output_file / "toxinpred3_out.csv"
    
    # Ensure output directory exists
    output_dir = output_file.parent
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Created output directory: {output_dir}")

    def unzip_model(zip_path, extract_to):
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_to)
            print(f"Unzipped {zip_path} to {extract_to}")
            return True
        except zipfile.BadZipFile:
            print(f"Error: {zip_path} is not a valid zip file.")
            return False
        except Exception as e:
            print(f"Error unzipping {zip_path}: {e}")
            return False

    # Locate model file
    model_found = False
    if os.path.isfile(model_path):
        print(f"Model file found at {model_path}")
        model_found = True
    elif os.path.isfile(nested_model_path):
        print(f"Model file found at nested path {nested_model_path}. Copying to {model_path}")
        shutil.copy(nested_model_path, model_path)
        model_found = True
    elif os.path.exists(zip_model_path):
        print(f"Found zipped model file at {zip_model_path}. Attempting to unzip.")
        unzip_dir = os.path.join(model_base_path, "model")
        if unzip_model(zip_model_path, unzip_dir):
            unzipped_model = os.path.join(unzip_dir, "toxinpred3.0_model.pkl")
            if os.path.isfile(unzipped_model):
                print(f"Unzipped model file to {unzipped_model}. Copying to {model_path}")
                shutil.copy(unzipped_model, model_path)
                model_found = True
            elif os.path.isfile(nested_model_path):
                print(f"Unzipped model file found at {nested_model_path}. Copying to {model_path}")
                shutil.copy(nested_model_path, model_path)
                model_found = True

    # Alternative model file without .pkl extension
    if not model_found:
        alt_model_path = os.path.join(model_base_path, "toxinpred3.0_model")
        if os.path.isfile(alt_model_path):
            print(f"Found model file without .pkl extension at {alt_model_path}. Renaming to {model_path}")
            os.rename(alt_model_path, model_path)
            model_found = True

    if not model_found:
        print("Error: Model file not found at any expected location")
        print(f"Checked: {model_path}, {nested_model_path}, {zip_model_path}, {alt_model_path}")
        print("Please ensure the model file (toxinpred3.0_model.pkl or toxinpred3.0_model) is available.")
        exit(1)

    # Verify model path is a regular file
    if not os.path.isfile(model_path):
        print(f"Error: Model path {model_path} is not a regular file (might be a directory).")
        exit(1)

    # Check toxinpred script exists
    if not os.path.isfile(toxinpred_script):
        print(f"Error: ToxinPred3.0 script not found at {toxinpred_script}")
        exit(1)

    # Check input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file not found at {input_file}")
        exit(1)

    # Validate input FASTA file format
    try:
        with open(input_file, 'r') as f:
            content = f.read().strip()
            if not content:
                print("Error: Input file is empty")
                exit(1)
            if not content.startswith('>'):
                print("Error: Input file does not appear to be in FASTA format (missing '>' header)")
                exit(1)
            seq_count = content.count('>')
            print(f"Input file contains {seq_count} sequence(s)")
    except Exception as e:
        print(f"Error reading input file: {e}")
        exit(1)

    # Parameters for ToxinPred3.0
    threshold = 0.38
    model = 2
    display = 2

    # Ensure model file symlink or copy in expected location
    expected_model_dir = os.path.join(working_directory, "model")
    expected_model_path = os.path.join(expected_model_dir, "toxinpred3.0_model.pkl")

    if os.path.exists(expected_model_path):
        if os.path.isdir(expected_model_path):
            print(f"Removing invalid model path (directory): {expected_model_path}")
            shutil.rmtree(expected_model_path)
        elif not os.path.isfile(expected_model_path):
            print(f"Removing invalid model path: {expected_model_path}")
            os.remove(expected_model_path)
        else:
            # Check if symlink points correctly
            if os.path.islink(expected_model_path):
                target = os.readlink(expected_model_path)
                if os.path.abspath(target) == os.path.abspath(model_path):
                    print("Model symlink is correct.")
                else:
                    print("Replacing incorrect symlink at model path.")
                    os.remove(expected_model_path)
                    os.makedirs(expected_model_dir, exist_ok=True)
                    os.symlink(model_path, expected_model_path)
                    print(f"Created symlink for model file at {expected_model_path}")
            else:
                print("Replacing existing model file with symlink.")
                os.remove(expected_model_path)
                os.makedirs(expected_model_dir, exist_ok=True)
                os.symlink(model_path, expected_model_path)
                print(f"Created symlink for model file at {expected_model_path}")
    else:
        os.makedirs(expected_model_dir, exist_ok=True)
        try:
            os.symlink(model_path, expected_model_path)
            print(f"Created symlink for model file at {expected_model_path}")
        except OSError as e:
            print(f"Error creating symlink at {expected_model_path}: {e}")
            exit(1)

    # Execute ToxinPred3.0 script
    try:
        cmd = [
            "python3",
            toxinpred_script,
            "-i", str(input_file),
            "-o", str(output_file),
            "-t", str(threshold),
            "-m", str(model),
            "-d", str(display)
        ]
        print(f"Executing ToxinPred3.0 with command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            cwd=working_directory,
            capture_output=True,
            text=True,
            check=True
        )
        print(f"ToxinPred3.0 stdout: {result.stdout}")
        if result.stderr:
            print(f"ToxinPred3.0 stderr: {result.stderr}")

        # Verify output file exists and is non-empty
        if not os.path.isfile(output_file):
            raise RuntimeError(f"ToxinPred3.0 output file not found at {output_file}")
        if os.path.getsize(output_file) == 0:
            raise RuntimeError(f"ToxinPred3.0 output file is empty at {output_file}")

        print(f"[✓] ToxinPred3.0 completed successfully, output saved to {output_file}")
        return output_file

    except subprocess.CalledProcessError as e:
        print(f"Error executing ToxinPred3.0: {e}")
        print(f"Command: {' '.join(e.cmd)}")
        print(f"Return code: {e.returncode}")
        print(f"Stderr: {e.stderr}")
        raise RuntimeError(f"ToxinPred3.0 execution failed: {e}")
    except Exception as e:
        print(f"Unexpected error during ToxinPred3.0 execution: {e}")
        raise RuntimeError(f"ToxinPred3.0 execution failed: {e}")

def main():
    """Parse command-line arguments and run ToxinPred3.0."""
    parser = argparse.ArgumentParser(description="Run ToxinPred3.0 for toxin epitope prediction")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    run_toxinpred3(args.input, args.output)

if __name__ == "__main__":
    main()
