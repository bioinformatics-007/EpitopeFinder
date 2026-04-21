import os
import subprocess
import pandas as pd
from pathlib import Path
from io import StringIO
import logging
import time
import psutil
import signal
import tempfile
import re
import shutil

# Set up logging
logger = logging.getLogger('VaxElan')
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Dynamically resolve root directory up to Vaxelan_2_0
current_path = Path(__file__).resolve()
while current_path.name != "Vaxelan_2_0":
    current_path = current_path.parent
    if current_path == current_path.parent:
        logger.error("Could not locate 'Vaxelan_2_0' in path hierarchy.")
        raise RuntimeError("Could not locate 'Vaxelan_2_0' in path hierarchy.")
root_dir = current_path

# MHC-I tool path
MHC_I_TOOL = os.environ.get('MHC_I_TOOL_PATH', os.path.join(root_dir, "tools/mhc_i/src/predict_binding.py"))

# Path to required pickle file
REQUIRED_PICKLE_FILE = os.environ.get('MHC_I_PICKLE_FILE', os.path.join(root_dir, "tools/mhc_i/method/allele-info/allele_info/pickles/mhci_info_dict.p"))

# NetMHCpan executable path
NETMHCPAN_PATH = os.environ.get('NETMHCPAN', os.path.join(root_dir, "tools/mhc_i/method/netmhcpan-4.1-executable/netmhcpan_4_1_executable/Linux_x86_64"))

# Hard-coded prediction method
PREDICTION_METHOD_CODE = "f"  # netmhcpan_el

# Method code mapping
METHOD_CODES = {
    "a": "ann",
    "b": "comblib_sidney2008",
    "c": "consensus",
    "d": "netmhccons",
    "e": "netmhcpan_ba",
    "f": "netmhcpan_el",
    "g": "netmhcstabpan",
    "h": "pickpocket",
    "i": "smm",
    "j": "smmpmbec"
}

# HLA alleles mapping based on Clade_Name
HLA_ALLELES = {
    "A1": "HLA-A*01:01",
    "A2.1": "HLA-A*02:01",
    "A2.2.1": "HLA-A*02:01",
    "A2.2.2": "HLA-A*02:01",
    "B1": "HLA-B*07:02",
    "B2": "HLA-B*08:01",
    "Ia": "HLA-A*01:01",
    "Ib": "HLA-B*07:02",
    "IIa": "HLA-A*02:01",
    "IIb": "HLA-B*08:01"
}

# Allowed peptide lengths for MHC-I
ALLOWED_LENGTHS = [8, 9, 10, 11]

def print_status(msg, status="info"):
    """Print colored status messages."""
    colors = {
        "info": "\033[94m",  # Blue
        "success": "\033[92m",  # Green
        "warning": "\033[93m",  # Yellow
        "error": "\033[91m"  # Red
    }
    endc = "\033[0m"
    print(f"{colors.get(status, '')}{msg}{endc}")
    if status == "info":
        logger.info(msg)
    elif status == "warning":
        logger.warning(msg)
    else:
        logger.error(msg)

def log_memory_usage():
    """Log current memory usage of the process."""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    mem_usage = mem_info.rss / 1024 / 1024
    logger.info(f"Memory usage: RSS={mem_usage:.2f} MB, VMS={mem_info.vms / 1024 / 1024:.2f} MB")
    return mem_usage

def is_valid_sequence(seq):
    """Check if a sequence is valid (contains only amino acids and is non-empty)."""
    if not isinstance(seq, str):
        logger.error(f"Invalid sequence type: {type(seq)}")
        return False
    seq = seq.strip()
    amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    is_valid = bool(seq) and all(c in amino_acids for c in seq.upper())
    if not is_valid:
        logger.error(f"Invalid sequence: {seq} (contains non-amino acid characters or empty)")
    return is_valid

def pad_sequence(seq, target_length=9):
    """Preserve original sequence without padding or trimming."""
    seq = seq.strip()
    if not is_valid_sequence(seq):
        logger.error(f"Cannot process invalid sequence: {seq}")
        return None
    if len(seq) not in ALLOWED_LENGTHS:
        logger.warning(f"Sequence '{seq}' (length={len(seq)}) not in allowed lengths {ALLOWED_LENGTHS}")
        return None
    logger.info(f"Using sequence '{seq}' (length={len(seq)}) as is")
    return seq

def generate_peptides(sequence, seq_id):
    """Generate peptides, keeping original sequence if within allowed lengths."""
    if not is_valid_sequence(sequence):
        logger.error(f"Cannot generate peptides for invalid sequence: {sequence}")
        return []
    seq_len = len(sequence)
    peptides = []
    if seq_len in ALLOWED_LENGTHS:
        peptides.append((sequence, f"{seq_id}_len{seq_len}"))
    else:
        logger.warning(f"Sequence '{sequence}' (length={seq_len}) not in allowed lengths {ALLOWED_LENGTHS}")
    logger.info(f"Generated {len(peptides)} peptides for sequence ID {seq_id}: {[p[0] for p in peptides]}")
    return peptides

def create_fasta_file(sequences, sequence_ids, output_file):
    """Create a FASTA file with original sequences, no length enforcement."""
    valid_sequences = []
    valid_ids = []
    for seq_id, seq in zip(sequence_ids, sequences):
        if not is_valid_sequence(seq):
            logger.error(f"Skipping invalid sequence: '{seq}'")
            continue
        if len(seq) not in ALLOWED_LENGTHS:
            logger.warning(f"Skipping sequence '{seq}' (length={len(seq)}) not in allowed lengths {ALLOWED_LENGTHS}")
            continue
        valid_sequences.append(seq)
        valid_ids.append(seq_id)
        logger.debug(f"Added sequence '{seq}' (length={len(seq)}) with ID {seq_id} to {output_file}")
    if not valid_sequences:
        logger.warning(f"No valid sequences for FASTA file {output_file}")
        return False
    with open(output_file, 'w') as f:
        for seq_id, seq in zip(valid_ids, valid_sequences):
            f.write(f">{seq_id}\n{seq}\n")
    # Validate FASTA file
    with open(output_file, 'r') as f:
        content = f.read()
    if not content.strip():
        logger.error(f"FASTA file {output_file} is empty.")
        return False
    lines = content.splitlines()
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        if line.startswith('>'):
            continue
        seq = line.strip()
        if not re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', seq):
            logger.error(f"FASTA file {output_file} contains invalid sequence at line {i+1}: {seq}")
            return False
        if len(seq) not in ALLOWED_LENGTHS:
            logger.error(f"FASTA file {output_file} contains sequence of incorrect length at line {i+1}: '{seq}' (length={len(seq)}, expected={ALLOWED_LENGTHS})")
            return False
    logger.info(f"Created and validated FASTA file: {output_file} with {len(valid_sequences)} sequences")
    with open(output_file, 'r') as f:
        content = f.read()
    logger.debug(f"FASTA file content:\n{content}")
    return True

def find_python_interpreter():
    """Find the current Python 3 interpreter."""
    python_cmd = shutil.which("python3")
    if not python_cmd:
        logger.error("No Python 3 interpreter found. Ensure Python 3 is installed.")
        return None
    try:
        result = subprocess.run(
            [python_cmd, "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        version_output = result.stderr or result.stdout
        if "Python 3" in version_output:
            logger.info(f"Using Python 3 interpreter: {python_cmd}")
            return python_cmd
        else:
            logger.error(f"Interpreter {python_cmd} is not Python 3: {version_output.strip()}")
            return None
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.error(f"Error checking interpreter {python_cmd}: {e}")
        return None

def extract_clade_identifier(clade_name):
    """Extract the clade identifier (e.g., 'IIb' from 'Envelop protein A28_IIb')."""
    if pd.isna(clade_name) or not isinstance(clade_name, str):
        logger.error(f"Invalid Clade_Name: {clade_name}")
        return None
    cleaned_clade = clade_name.replace('__', '_')
    clade_identifier = cleaned_clade.split('_')[-1]
    logger.debug(f"Extracted clade identifier '{clade_identifier}' from Clade_Name '{clade_name}'")
    if clade_identifier not in HLA_ALLELES:
        logger.error(f"Clade identifier '{clade_identifier}' not found in HLA_ALLELES: {list(HLA_ALLELES.keys())}")
        return None
    return clade_identifier

def validate_netmhcpan_path():
    """Validate NetMHCpan executable and required files."""
    if not os.path.exists(NETMHCPAN_PATH):
        print_status(f"NetMHCpan executable not found at: {NETMHCPAN_PATH}", "error")
        return False
    version_file = os.path.join(NETMHCPAN_PATH, "data/version")
    if not os.path.isfile(version_file):
        print_status(f"NetMHCpan version file not found at: {version_file}", "error")
        return False
    return True

def run_mhc_predictor(method, allele, input_file, output_dir, temp_dir, timeout=3600):
    """Run the MHC-I prediction tool with retry mechanism for sequence length issues."""
    if not os.path.isfile(MHC_I_TOOL):
        print_status(f"MHC-I tool '{MHC_I_TOOL}' not found.", "error")
        return None, "MHC-I tool not found"
    if not os.path.isfile(input_file):
        print_status(f"Input file '{input_file}' not found.", "error")
        return None, "Input file not found"
    if not os.path.isfile(REQUIRED_PICKLE_FILE):
        print_status(f"Required pickle file not found at: {REQUIRED_PICKLE_FILE}", "error")
        return None, "Required pickle file not found"
    if not validate_netmhcpan_path():
        return None, "NetMHCpan validation failed"
    # Read and log FASTA file content
    try:
        with open(input_file, 'r') as f:
            fasta_content = f.read()
        logger.debug(f"FASTA file {input_file} content:\n{fasta_content}")
        sequences = []
        sequence_ids = []
        peptide_lengths = set()
        current_id = None
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    current_id = line.strip()[1:]
                    continue
                seq = line.strip()
                if seq and is_valid_sequence(seq) and current_id:
                    if len(seq) in ALLOWED_LENGTHS:
                        peptide_lengths.add(len(seq))
                        sequences.append(seq)
                        sequence_ids.append(current_id)
                        logger.debug(f"Validated sequence '{seq}' (length={len(seq)}) with ID {current_id}")
                    else:
                        logger.warning(f"Skipping sequence '{seq}' (length={len(seq)}) not in allowed lengths {ALLOWED_LENGTHS}")
                elif not seq and current_id:
                    logger.warning(f"Empty sequence for ID {current_id} in {input_file}")
    except Exception as e:
        print_status(f"Error reading FASTA file {input_file}: {e}", "error")
        return None, f"Error reading FASTA file: {e}"
    if not peptide_lengths:
        print_status(f"No peptides with valid lengths {ALLOWED_LENGTHS} in {input_file}", "error")
        return None, f"No peptides with valid lengths {ALLOWED_LENGTHS}"
    python_cmd = find_python_interpreter()
    if not python_cmd:
        print_status("No Python 3 interpreter found.", "error")
        return None, "No Python 3 interpreter found"
    debug_dir = os.path.join(output_dir, "mhci_debug")
    os.makedirs(debug_dir, exist_ok=True)
    safe_allele = allele.replace('*', '_').replace(':', '_')
    combined_output = []
    for length in sorted(peptide_lengths):
        # Create length-specific FASTA file
        length_fasta = os.path.join(temp_dir, f"epitopes_len{length}.fasta")
        if not create_fasta_file(
            [seq for seq in sequences if len(seq) == length],
            [seq_id for seq, seq_id in zip(sequences, sequence_ids) if len(seq) == length],
            length_fasta
        ):
            print_status(f"No valid sequences for length {length}, skipping", "warning")
            continue
        # Copy length-specific FASTA to debug directory
        debug_fasta = os.path.join(debug_dir, f"epitopes_len{length}.fasta")
        shutil.copy(length_fasta, debug_fasta)
        # Log FASTA content for this length
        with open(length_fasta, 'r') as f:
            length_fasta_content = f.read()
        logger.debug(f"FASTA file {length_fasta} content for length {length}:\n{length_fasta_content}")
        # Try prediction with specified length
        command = [python_cmd, MHC_I_TOOL, method, allele, str(length), length_fasta]
        command_str = ' '.join(map(str, command))
        print_status(f"Running command for length {length}: {command_str}", "info")
        debug_output = os.path.join(debug_dir, f"mhc1_{safe_allele}_len{length}.txt")
        debug_error = debug_output + ".error"
        try:
            mem_usage = log_memory_usage()
            if mem_usage > 0.8 * psutil.virtual_memory().total / 1024 / 1024:
                print_status(f"Warning: High memory usage ({mem_usage:.2f} MB).", "warning")
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout,
                preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL)
            )
            with open(debug_output, "w") as f:
                f.write(result.stdout)
            if result.stderr:
                print_status(f"Error output for length {length}: {result.stderr}", "warning")
                with open(debug_error, "w") as f:
                    f.write(result.stderr)
            if result.stdout.strip():
                combined_output.append(result.stdout)
                logger.debug(f"Output for {allele} length {length}:\n{result.stdout[:1000]}...")
            else:
                print_status(f"No output from MHC-I tool for {allele} at length {length}", "warning")
                if "too short" in result.stderr.lower():
                    logger.warning(f"MHC-I tool rejected length {length} sequences as too short. Retrying with individual sequences.")
                    # Retry each sequence individually
                    for seq, seq_id in zip(sequences, sequence_ids):
                        if len(seq) != length:
                            continue
                        single_fasta = os.path.join(temp_dir, f"single_seq_{seq_id}.fasta")
                        with open(single_fasta, 'w') as f:
                            f.write(f">{seq_id}\n{seq}\n")
                        single_command = [python_cmd, MHC_I_TOOL, method, allele, str(length), single_fasta]
                        single_command_str = ' '.join(map(str, single_command))
                        logger.debug(f"Retrying single sequence '{seq}' with command: {single_command_str}")
                        try:
                            single_result = subprocess.run(
                                single_command,
                                capture_output=True,
                                text=True,
                                timeout=timeout,
                                preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL)
                            )
                            single_debug_output = os.path.join(debug_dir, f"mhc1_{safe_allele}_single_{seq_id}.txt")
                            single_debug_error = single_debug_output + ".error"
                            with open(single_debug_output, "w") as f:
                                f.write(single_result.stdout)
                            if single_result.stderr:
                                with open(single_debug_error, "w") as f:
                                    f.write(single_result.stderr)
                            if single_result.stdout.strip():
                                combined_output.append(single_result.stdout)
                                logger.debug(f"Successful retry for sequence '{seq}' with output:\n{single_result.stdout[:1000]}...")
                            else:
                                logger.warning(f"Retry failed for sequence '{seq}': stderr={single_result.stderr}")
                        except Exception as e:
                            logger.error(f"Retry failed for sequence '{seq}': {e}")
                else:
                    return None, f"No output from MHC-I tool for length {length}: stderr={result.stderr}"
        except subprocess.TimeoutExpired as e:
            print_status(f"Timeout after {timeout}s for {allele} at length {length}: {e}", "error")
            return None, f"Timeout after {timeout}s for length {length}: {e}"
        except subprocess.CalledProcessError as e:
            print_status(f"Command failed for {allele} at length {length}: {e.stderr}", "error")
            return None, f"Command failed for length {length}: {e.stderr}"
        except Exception as e:
            print_status(f"Unexpected error for {allele} at length {length}: {e}", "error")
            return None, f"Unexpected error for length {length}: {e}"
    if not combined_output:
        print_status(f"No output from MHC-I tool for {allele} across all lengths", "error")
        return None, f"No output from MHC-I tool for {allele} across all lengths"
    return "\n".join(combined_output), None

def parse_results(response_text, method, allele, score_threshold=100000):
    """Parse prediction results into a DataFrame with enhanced error handling."""
    if response_text is None or not isinstance(response_text, str) or not response_text.strip():
        print_status(f"No results for {allele}", "warning")
        logger.debug(f"No response text for {allele}")
        return pd.DataFrame(), "No response text"
    logger.debug(f"Raw MHC-I output for {allele}:\n{response_text[:1000]}...")
    try:
        data_io = StringIO(response_text)
        df = pd.read_csv(data_io, sep=r'\s*,\s*|\t|\s+', engine="python", skipinitialspace=True)
        if df.empty:
            print_status(f"Empty DataFrame after parsing for {allele}", "warning")
            logger.debug(f"Empty DataFrame after parsing for {allele}")
            return pd.DataFrame(), "Empty DataFrame after parsing"
        df.columns = [col.strip().lower().replace(' ', '_') for col in df.columns]
        logger.info(f"Parsed columns for {allele}: {list(df.columns)}")
        logger.debug(f"First few rows of parsed DataFrame:\n{df.head().to_string()}")
        peptide_col = next((col for col in df.columns if col.lower() in ['peptide', 'seq', 'sequence', 'allele_peptide', 'pep']), None)
        score_col = next((col for col in df.columns if 'affinity' in col.lower() or 'ic50' in col.lower() or 'score' in col.lower()), None)
        if not peptide_col or not score_col:
            print_status(f"No peptide or score column found in results for {allele}. Available columns: {list(df.columns)}", "error")
            logger.debug(f"No peptide or score column for {allele}. Columns: {list(df.columns)}")
            return pd.DataFrame(), f"No peptide or score column found. Available columns: {list(df.columns)}"
        df = df[[peptide_col, score_col]].rename(columns={peptide_col: 'Sequence', score_col: 'Binding_Affinity'})
        df['Binding_Affinity'] = pd.to_numeric(df['Binding_Affinity'], errors='coerce')
        invalid_rows = df[df['Binding_Affinity'].isna()]
        if not invalid_rows.empty:
            logger.warning(f"Rows with invalid binding affinities for {allele}:\n{invalid_rows.to_string()}")
        df = df[df['Binding_Affinity'].notna() & (df['Binding_Affinity'] <= score_threshold)]
        if df.empty:
            print_status(f"No peptides with valid scores <= {score_threshold} for {allele}", "warning")
            return pd.DataFrame(), f"No peptides with valid scores <= {score_threshold}"
        print_status(f"Parsed {len(df)} valid peptides for {allele}", "info")
        df['Allele'] = allele
        return df, None
    except Exception as e:
        print_status(f"Error parsing results for {allele}: {str(e)}", "error")
        logger.debug(f"Error parsing results for {allele}: {str(e)}")
        logger.debug(f"Raw response text: {response_text[:1000]}...")
        return pd.DataFrame(), f"Error parsing results: {str(e)}"

def process_input_file(input_file, method, output_dir, score_threshold=100000):
    """Process the input CSV file and save results for all HLA alleles."""
    output_csv = os.path.join(output_dir, "binding_affinity_results_mhci.csv")
    skipped_csv = os.path.join(output_dir, "skipped_rows_mhci.csv")
    try:
        df = pd.read_csv(input_file)
        print_status(f"Successfully read file: {input_file}", "info")
    except Exception as e:
        print_status(f"Error reading {input_file}: {e}", "error")
        return
    # Determine Clade_Name from input file name
    file_name = os.path.basename(input_file).lower()
    if 'ia' in file_name:
        clade_name = 'Envelop protein A28_Ia'
    elif 'ib' in file_name:
        clade_name = 'Envelop protein A28_Ib'
    elif 'iia' in file_name:
        clade_name = 'Envelop protein A28_IIa'
    elif 'iib' in file_name:
        clade_name = 'Envelop protein A28_IIb'
    else:
        clade_name = 'Envelop protein A28_Unknown'
    # Update Epitope_Type and Clade_Name
    df['Epitope_Type'] = 'Mhci'
    df['Clade_Name'] = clade_name
    logger.info(f"Assigned Clade_Name: {clade_name} for file: {input_file}")
    logger.info(f"Filtered {len(df)} MHC-I epitopes")
    with tempfile.TemporaryDirectory() as temp_dir:
        ref_sequences = df['Master_Sequence'].dropna().astype(str).str.strip().tolist()
        mut_sequences = df['Clade_Sequence'].dropna().astype(str).str.strip().tolist()
        valid_sequences = []
        valid_sequence_ids = []
        invalid_sequences = []
        sequence_to_id = {}
        subsequence_to_original = {}
        for i, seq in enumerate(set(ref_sequences + mut_sequences)):
            seq_id = f"seq_{i}"
            if not is_valid_sequence(seq):
                invalid_sequences.append((seq, "invalid amino acids"))
                continue
            # Generate peptides
            peptides = generate_peptides(seq, seq_id)
            if not peptides:
                invalid_sequences.append((seq, "failed peptide generation"))
                continue
            sequence_to_id[seq] = seq_id
            for sub_seq, sub_seq_id in peptides:
                valid_sequences.append(sub_seq)
                valid_sequence_ids.append(sub_seq_id)
                subsequence_to_original[sub_seq] = seq
        if invalid_sequences:
            print_status(f"Invalid sequences: {[(seq[:10] + '...' if isinstance(seq, str) and len(seq) > 10 else seq, reason) for seq, reason in invalid_sequences]}", "warning")
        if not valid_sequences:
            print_status("No valid sequences found", "error")
            return
        print_status(f"Valid peptides: {[seq for seq in valid_sequences]}", "info")
        fasta_file = os.path.join(temp_dir, "epitopes.fasta")
        if not create_fasta_file(valid_sequences, valid_sequence_ids, fasta_file):
            print_status("Failed to create valid FASTA file", "error")
            return
        debug_fasta = os.path.join(output_dir, "mhci_debug", "epitopes.fasta")
        os.makedirs(os.path.dirname(debug_fasta), exist_ok=True)
        create_fasta_file(valid_sequences, valid_sequence_ids, debug_fasta)
        all_results = []
        skipped_rows = []
        processed_alleles = set()
        for _, row in df.iterrows():
            epitope_id = row['Epitope_ID']
            seq_id = row['Sequence_ID_Clade']
            ref_seq = str(row['Master_Sequence']).strip()
            mut_seq = str(row['Clade_Sequence']).strip()
            current_clade_name = row['Clade_Name']
            epitope_type = row['Epitope_Type']
            clade_identifier = extract_clade_identifier(current_clade_name)
            if not clade_identifier or clade_identifier not in HLA_ALLELES:
                print_status(f"Invalid or missing Clade_Name: {current_clade_name} (extracted: {clade_identifier}) for {epitope_id}", "warning")
                skipped_rows.append({
                    'Epitope_ID': epitope_id,
                    'Sequence_ID_Clade': seq_id,
                    'Allele': None,
                    'Master_Sequence': ref_seq,
                    'Clade_Sequence': mut_seq,
                    'Clade_Name': current_clade_name,
                    'Epitope_Type': epitope_type,
                    'Reason': f"Invalid or missing Clade_Name: {current_clade_name} (extracted: {clade_identifier})"
                })
                continue
            hla_allele = HLA_ALLELES[clade_identifier]
            print_status(f"Processing {epitope_id} for {hla_allele} (Clade_Name: {current_clade_name}, extracted: {clade_identifier})", "info")
            if hla_allele not in processed_alleles:
                response_text, error_reason = run_mhc_predictor(method, hla_allele, fasta_file, output_dir, temp_dir)
                if response_text is None:
                    print_status(f"Failed to obtain binding scores for {hla_allele}: {error_reason}", "error")
                    skipped_rows.append({
                        'Epitope_ID': epitope_id,
                        'Sequence_ID_Clade': seq_id,
                        'Allele': hla_allele,
                        'Master_Sequence': ref_seq,
                        'Clade_Sequence': mut_seq,
                        'Clade_Name': current_clade_name,
                        'Epitope_Type': epitope_type,
                        'Reason': f"Failed to obtain binding scores: {error_reason}"
                    })
                    continue
                binding_df, parse_error = parse_results(response_text, method, hla_allele, score_threshold)
                if binding_df.empty:
                    print_status(f"No valid binding scores for {hla_allele}: {parse_error}", "warning")
                    skipped_rows.append({
                        'Epitope_ID': epitope_id,
                        'Sequence_ID_Clade': seq_id,
                        'Allele': hla_allele,
                        'Master_Sequence': ref_seq,
                        'Clade_Sequence': mut_seq,
                        'Clade_Name': current_clade_name,
                        'Epitope_Type': epitope_type,
                        'Reason': f"No valid binding scores: {parse_error}"
                    })
                    continue
                all_results.append(binding_df)
                processed_alleles.add(hla_allele)
        if not all_results:
            print_status("No valid results for any allele", "error")
            if skipped_rows:
                pd.DataFrame(skipped_rows).to_csv(skipped_csv, index=False)
                print_status(f"Skipped rows saved to {skipped_csv}", "info")
            return
        combined_binding_df = pd.concat(all_results, ignore_index=True)
        binding_scores = {}
        for _, row in combined_binding_df.iterrows():
            sub_seq = row['Sequence']
            score = row['Binding_Affinity']
            allele = row['Allele']
            original_seq = subsequence_to_original.get(sub_seq)
            if original_seq:
                key = (original_seq, allele)
                if key not in binding_scores or score < binding_scores[key]:
                    binding_scores[key] = score
        clade_results = []
        for _, row in df.iterrows():
            epitope_id = row['Epitope_ID']
            seq_id = row['Sequence_ID_Clade']
            ref_seq = str(row['Master_Sequence']).strip()
            mut_seq = str(row['Clade_Sequence']).strip()
            current_clade_name = row['Clade_Name']
            epitope_type = row['Epitope_Type']
            clade_identifier = extract_clade_identifier(current_clade_name)
            if not clade_identifier or clade_identifier not in HLA_ALLELES:
                continue
            hla_allele = HLA_ALLELES[clade_identifier]
            ref_score = binding_scores.get((ref_seq, hla_allele), None)
            mut_score = binding_scores.get((mut_seq, hla_allele), None)
            if ref_score is None or mut_score is None:
                reason = f"Missing binding score for {hla_allele}: ref_seq={ref_seq[:10]}... or mut_seq={mut_seq[:10]}..."
                skipped_rows.append({
                    'Epitope_ID': epitope_id,
                    'Sequence_ID_Clade': seq_id,
                    'Allele': hla_allele,
                    'Master_Sequence': ref_seq,
                    'Clade_Sequence': mut_seq,
                    'Clade_Name': current_clade_name,
                    'Epitope_Type': epitope_type,
                    'Reason': reason
                })
                continue
            affinity_diff = ref_score - mut_score
            clade_results.append({
                'Clade_Name': current_clade_name,
                'Epitope_Type': epitope_type,
                'Epitope_ID': epitope_id,
                'Sequence_ID_Clade': seq_id,
                'Master_Sequence': ref_seq,
                'Clade_Sequence': mut_seq,
                'Allele': hla_allele,
                'Master_Binding_Affinity_nM': ref_score,
                'Mutated_Binding_Affinity_nM': mut_score,
                'Affinity_Difference_nM': affinity_diff
            })
        if clade_results:
            results_df = pd.DataFrame(clade_results)
            results_df = results_df[[
                'Clade_Name', 'Epitope_Type', 'Epitope_ID', 'Sequence_ID_Clade',
                'Master_Sequence', 'Clade_Sequence', 'Allele', 'Master_Binding_Affinity_nM',
                'Mutated_Binding_Affinity_nM', 'Affinity_Difference_nM'
            ]]
            results_df.to_csv(output_csv, index=False)
            print_status(f"Results saved to {output_csv}", "success")
        else:
            print_status("No valid results to save", "warning")
        if skipped_rows:
            pd.DataFrame(skipped_rows).to_csv(skipped_csv, index=False)
            print_status(f"Skipped rows saved to {skipped_csv}", "info")

def run_mhc1():
    """Run MHC-I prediction for the specified input files."""
    input_files = [
        "/home/yuktika/Downloads/A_28/A28/A28_entry_fusion_protein_sequences_Ia/mhci_mutations_A28_entry_fusion_protein_sequences_Ia.csv",
        "/home/yuktika/Downloads/A_28/A28/A28_entry_fusion_protein_sequences_Ib/mhci_mutations_A28_entry_fusion_protein_sequences_Ib.csv",
        "/home/yuktika/Downloads/A_28/A28/A28_entry_fusion_protein_sequences_IIa/mhci_mutations_A28_entry_fusion_protein_sequences_IIa.csv",
        "/home/yuktika/Downloads/A_28/A28/A28_entry_fusion_protein_sequences_IIb/mhci_mutations_A28_entry_fusion_protein_sequences_IIb.csv"
    ]
    method = METHOD_CODES[PREDICTION_METHOD_CODE]
    score_threshold = 100000
    for input_file in input_files:
        if not os.path.isfile(input_file):
            print_status(f"Input file '{input_file}' does not exist.", "error")
            continue
        if not input_file.lower().endswith('.csv'):
            print_status(f"Input file '{input_file}' is not a CSV file.", "error")
            continue
        output_dir = os.path.dirname(input_file)
        os.makedirs(output_dir, exist_ok=True)
        start_time = time.time()
        print_status(f"\n**Starting MHC Class I Prediction with {method} for {input_file}**", "info")
        process_input_file(input_file, method, output_dir, score_threshold)
        output_csv = os.path.join(output_dir, "binding_affinity_results_mhci.csv")
        if not os.path.exists(output_csv):
            print_status(f"Processing failed for {input_file}: No output file generated", "error")
            continue
        elapsed_time = time.time() - start_time
        print_status(f"\n✔ Prediction complete for {input_file} in {elapsed_time:.2f} seconds.", "success")

if __name__ == "__main__":
    run_mhc1()
