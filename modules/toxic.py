import subprocess
import os
import shutil
import zipfile
from pathlib import Path
import pandas as pd
import tempfile
import re

def run_toxinpred3(input_file, output_dir):
    # Configuration
    CONFIG = {
        "batch_size": 10000,
        "temp_base_path": tempfile.gettempdir()
    }
   
    # Determine the root directory two levels above this script
    ROOT_DIR = Path(__file__).resolve().parents[1]
    # Define paths
    toxinpred_script = os.path.join(ROOT_DIR, "tools", "toxinpred3", "toxinpred3.py")
    model_base_path = os.path.join(ROOT_DIR, "tools", "toxinpred3")
    model_path = os.path.join(model_base_path, "toxinpred3.0_model.pkl")
    nested_model_path = os.path.join(model_base_path, "model", "toxinpred3.0_model.pkl")
    zip_model_path = os.path.join(model_base_path, "model", "toxinpred3.0_model.pkl.zip")
    working_directory = os.path.join(ROOT_DIR, "tools", "toxinpred3")
    # Convert output_dir to Path object and define output file
    output_dir = Path(output_dir)
    output_file = output_dir / "toxinpred3_output.csv"
   
    # Ensure output directory exists
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Created output directory: {output_dir}")
   
    # Ensure temporary directory exists
    temp_base_dir = Path(CONFIG["temp_base_path"])
    if not temp_base_dir.exists():
        temp_base_dir.mkdir(parents=True, exist_ok=True)
        print(f"Created temporary directory: {temp_base_dir}")
   
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
   
    # Read FASTA file and convert to DataFrame
    try:
        sequences = []
        with open(input_file, 'r') as f:
            lines = f.readlines()
            current_id = None
            current_seq = []
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        # Split header on '|' to extract Master_Epitope_ID and Sequence_Clade_ID
                        parts = current_id.split('|')
                        if len(parts) == 2:
                            master_id, clade_id = parts
                        else:
                            master_id = current_id
                            clade_id = current_id
                        sequences.append({
                            'Master_Epitope_ID': master_id,
                            'Sequence_Clade_ID': clade_id,
                            'Mutated_Seq': ''.join(current_seq),
                            'Sequence_ID': current_id  # Store the full header as Sequence_ID
                        })
                    current_id = line[1:]  # Remove '>'
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id is not None:
                parts = current_id.split('|')
                if len(parts) == 2:
                    master_id, clade_id = parts
                else:
                    master_id = current_id
                    clade_id = current_id
                sequences.append({
                    'Master_Epitope_ID': master_id,
                    'Sequence_Clade_ID': clade_id,
                    'Mutated_Seq': ''.join(current_seq),
                    'Sequence_ID': current_id
                })
        df = pd.DataFrame(sequences)
        if df.empty:
            print("Error: No sequences found in FASTA file")
            exit(1)
        if df['Mutated_Seq'].isnull().any():
            print("Error: Input contains null values in sequences")
            exit(1)
        print(f"Input FASTA contains {len(df)} sequence(s)")
    except Exception as e:
        print(f"Error reading input FASTA file: {e}")
        exit(1)
   
    # Process sequences in chunks
    batch_size = CONFIG["batch_size"]
    result_dfs = []
   
    for start_idx in range(0, len(df), batch_size):
        end_idx = min(start_idx + batch_size, len(df))
        chunk_df = df[start_idx:end_idx]
        print(f"Processing chunk {start_idx // batch_size + 1} ({start_idx} to {end_idx})")
       
        # Convert chunk to FASTA format
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', dir=CONFIG["temp_base_path"], delete=False) as temp_fasta:
            for idx, row in chunk_df.iterrows():
                temp_fasta.write(f">{row['Sequence_ID']}\n")
                temp_fasta.write(f"{row['Mutated_Seq']}\n")
            temp_fasta_path = temp_fasta.name
       
        try:
            # Validate FASTA file format
            with open(temp_fasta_path, 'r') as f:
                content = f.read().strip()
                if not content:
                    print("Error: Generated FASTA file is empty")
                    exit(1)
                if not content.startswith('>'):
                    print("Error: Generated FASTA file does not start with '>' header")
                    exit(1)
                seq_count = content.count('>')
                print(f"Generated FASTA file contains {seq_count} sequence(s) for chunk")
           
            # Parameters for ToxinPred3.0
            threshold = 0.5
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
            temp_output = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', dir=CONFIG["temp_base_path"], delete=False)
            temp_output_path = temp_output.name
            temp_output.close()
           
            try:
                cmd = [
                    "python3",
                    toxinpred_script,
                    "-i", temp_fasta_path,
                    "-o", temp_output_path,
                    "-t", str(threshold),
                    "-m", str(model),
                    "-d", str(display)
                ]
                print(f"Executing ToxinPred3.0 with command: {' '.join(cmd)}")
                result = subprocess.run(
                    cmd,
                    cwd=working_directory,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    check=True,
                    stdin=subprocess.DEVNULL
                )
                print(f"ToxinPred3.0 stdout: {result.stdout}")
                if result.stderr:
                    print(f"ToxinPred3.0 stderr: {result.stderr}")
               
                # Verify output file exists and is non-empty
                if not os.path.isfile(temp_output_path):
                    raise RuntimeError(f"ToxinPred3.0 output file not found at {temp_output_path}")
                if os.path.getsize(temp_output_path) == 0:
                    raise RuntimeError(f"ToxinPred3.0 output file is empty at {temp_output_path}")
               
                # Read ToxinPred3.0 output and merge with chunk DataFrame
                toxinpred_df = pd.read_csv(temp_output_path)
                print("ToxinPred3.0 output columns:", toxinpred_df.columns.tolist())
                if 'Prediction' not in toxinpred_df.columns:
                    print("Error: ToxinPred3.0 output CSV missing 'Prediction' column")
                    exit(1)
               
                # Map Prediction to Is_Non_Toxin (0 for Toxin, 1 for Non-Toxin)
                toxinpred_df['Is_Non_Toxin'] = toxinpred_df['Prediction'].apply(lambda x: 1 if x == 'Non-Toxin' else 0)
               
                # Use 'Subject' as the sequence ID column
                id_column = 'Subject'
                if id_column not in toxinpred_df.columns:
                    print(f"Error: Sequence ID column '{id_column}' not found in ToxinPred3.0 output. Available columns: {toxinpred_df.columns.tolist()}")
                    exit(1)
               
                # Extract sequence ID
                toxinpred_df['Sequence_ID'] = toxinpred_df[id_column].str.lstrip('>')
               
                # Prepare input DataFrame
                input_df = chunk_df.copy()
                input_df['Epitope_ID'] = input_df['Master_Epitope_ID']
                input_df['Sequence_Clade_ID'] = input_df['Sequence_Clade_ID']
                input_df['Mutated_Seq'] = input_df['Mutated_Seq']
               
                # Merge results
                output_columns = ['Epitope_ID', 'Sequence_Clade_ID', 'Mutated_Seq', 'Prediction']
                score_column = None
                if 'Hybrid Score' in toxinpred_df.columns:
                    score_column = 'Hybrid Score'
                    output_columns.append('Hybrid Score')
                elif 'ML Score' in toxinpred_df.columns:
                    score_column = 'ML Score'
                    output_columns.append('ML Score')
                output_columns.append('Is_Non_Toxin')
               
                merge_columns = ['Sequence_ID', 'Prediction', 'Is_Non_Toxin']
                if score_column:
                    merge_columns.append(score_column)
                result_df = input_df.merge(toxinpred_df[merge_columns], on='Sequence_ID', how='left')
               
                # Select output columns
                result_df = result_df[output_columns]
               
                # Append to results
                result_dfs.append(result_df)
               
            except subprocess.CalledProcessError as e:
                print(f"Error executing ToxinPred3.0 for {input_file} chunk {start_idx // batch_size + 1}: {e}")
                print(f"Command: {' '.join(e.cmd)}")
                print(f"Return code: {e.returncode}")
                print(f"Stderr: {e.stderr}")
                raise RuntimeError(f"ToxinPred3.0 execution failed: {e}")
            except Exception as e:
                print(f"Unexpected error during ToxinPred3.0 execution for {input_file} chunk {start_idx // batch_size + 1}: {e}")
                raise RuntimeError(f"ToxinPred3.0 execution failed: {e}")
            finally:
                # Clean up temporary files
                if os.path.exists(temp_fasta_path):
                    os.remove(temp_fasta_path)
                if os.path.exists(temp_output_path):
                    os.remove(temp_output_path)
       
        except Exception as e:
            print(f"Error processing temporary FASTA file for {input_file} chunk {start_idx // batch_size + 1}: {e}")
            if os.path.exists(temp_fasta_path):
                os.remove(temp_fasta_path)
            raise RuntimeError(f"Processing failed: {e}")
   
    # Combine all chunk results
    final_df = pd.concat(result_dfs, ignore_index=True)
   
    # Save final output
    final_df.to_csv(output_file, index=False)
    print(f"[✓] ToxinPred3.0 completed successfully for {input_file}, output saved to {output_file}")
    return output_file

def main():
    # Define input files and corresponding output directories
    input_output_pairs = [
        (
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/Fasta_folder/H3_envelope_glycoprotein_sequences_Ia.fasta",
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/H3_envelope_glycoprotein_sequences_Ia"
        ),
        (
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/Fasta_folder/H3_envelope_glycoprotein_sequences_Ib.fasta",
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/H3_envelope_glycoprotein_sequences_Ib"
        ),
        (
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/Fasta_folder/H3_envelope_glycoprotein_sequences_IIa.fasta",
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/H3_envelope_glycoprotein_sequences_IIa"
        ),
        (
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/Fasta_folder/H3_envelope_glycoprotein_sequences_IIb.fasta",
            "/home/yuktika/Downloads/epitope_mapping/H3L/epitope_mutation_A0A7H0DN80 (H3L)/H3_envelope_glycoprotein_sequences_IIb"
        )
    ]
   
    # Process each input file
    for input_file, output_dir in input_output_pairs:
        print(f"Processing {input_file}...")
        run_toxinpred3(input_file, output_dir)

if __name__ == "__main__":
    main()
