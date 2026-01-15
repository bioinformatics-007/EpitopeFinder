import os
import subprocess
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
import re
import tempfile
import pandas as pd

# Define root directory
ROOT_DIR = Path("/home/yuktika/Downloads/Vaxelan_2_0")

# Define default database path
GUT_DB = ROOT_DIR / "tools/databases/gut_microbiome_db"

def print_status(message, status="info"):
    color_map = {
        "info": "\033[94m",
        "success": "\033[92m",
        "warning": "\033[93m",
        "error": "\033[91m"
    }
    print(f"{color_map.get(status, '')}{message}\033[0m")

def parse_fasta(file_path):
    """
    Parse all sequences from a FASTA file to extract protein ID, description, and sequence.
    """
    try:
        with open(file_path, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
            
            if not lines:
                print_status(f"FASTA file is empty: {file_path}", "error")
                raise ValueError("Empty FASTA file")

            sequences = []
            current_header = None
            current_seq_lines = []
            header_index = -1

            for i, line in enumerate(lines):
                if line.startswith(">"):
                    if current_header is not None:
                        seq = ''.join(current_seq_lines).replace(" ", "")
                        if not seq:
                            print_status(f"No sequence found for header: {current_header} in {file_path}", "error")
                            raise ValueError(f"No sequence for header: {current_header}")

                        valid_aa_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY*]+$', re.IGNORECASE)
                        if not valid_aa_pattern.match(seq):
                            print_status(f"Invalid sequence for {current_header} in {file_path}", "error")
                            raise ValueError(f"Invalid sequence: {seq[:50]}...")

                        header_parts = current_header[1:].split()
                        if not header_parts:
                            print_status(f"FASTA header lacks protein ID at line {header_index+1} in {file_path}", "error")
                            raise ValueError("FASTA header lacks protein ID")

                        prot_id = header_parts[0]
                        prot_desc = ' '.join(header_parts[1:]) if len(header_parts) > 1 else ""
                        sequences.append((prot_id, prot_desc, seq))
                       

                    current_header = line
                    current_seq_lines = []
                    header_index = i
                else:
                    if current_header is None:
                        print_status(f"Sequence data found before header in {file_path}", "error")
                        raise ValueError("Sequence data found before header")
                    current_seq_lines.append(line)

            if current_header is not None:
                seq = ''.join(current_seq_lines).replace(" ", "")
                if not seq:
                    print_status(f"No sequence found for header: {current_header} in {file_path}", "error")
                    raise ValueError(f"No sequence for header: {current_header}")

                valid_aa_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY*]+$', re.IGNORECASE)
                if not valid_aa_pattern.match(seq):
                    print_status(f"Invalid sequence for {current_header} in {file_path}", "error")
                    raise ValueError(f"Invalid sequence: {seq[:50]}...")

                header_parts = current_header[1:].split()
                if not header_parts:
                    print_status(f"FASTA header lacks protein ID at line {header_index+1} in {file_path}", "error")
                    raise ValueError("FASTA header lacks protein ID")

                prot_id = header_parts[0]
                prot_desc = ' '.join(header_parts[1:]) if len(header_parts) > 1 else ""
                sequences.append((prot_id, prot_desc, seq))
                

            if not sequences:
                print_status(f"No valid sequences found in {file_path}", "error")
                raise ValueError("No valid sequences found")

            return sequences
    except Exception as e:
        print_status(f"Error parsing FASTA file {file_path}: {e}", "error")
        raise

def write_temp_fasta(prot_id, prot_desc, prot_seq, temp_file):
    """Write a temporary FASTA file for BLASTP query."""
    try:
        with open(temp_file, "w") as f:
            f.write(f">{prot_id} {prot_desc}\n{prot_seq}\n")
        print_status(f"Created temporary FASTA file: {temp_file}", "info")
    except Exception as e:
        print_status(f"Error writing temporary FASTA file {temp_file}: {e}", "error")
        raise

def run_blastp(db_path, query_file, out_xml, temp_dir):
    """Run BLASTP against a database with relaxed parameters."""
    db_path_str = str(db_path)
    if not os.path.exists(db_path_str + ".phr"):
        print_status(f"BLAST database not found: {db_path}", "error")
        raise FileNotFoundError(f"BLAST database not found: {db_path}")

    query_file = str(Path(temp_dir) / query_file)
    out_xml = str(Path(temp_dir) / out_xml)

    cmd = [
        "blastp",
        "-query", query_file,
        "-db", db_path_str,
        "-out", out_xml,
        "-outfmt", "5",
        "-evalue", "1e-2",
        "-max_target_seqs", "50",
        "-num_threads", "4"
    ]
    print_status(f"Running BLASTP: {' '.join(cmd)}", "info")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=temp_dir)
        if result.stdout:
            print_status(f"BLASTP STDOUT: {result.stdout[:500]}...", "info")
        if result.stderr:
            print_status(f"BLASTP STDERR: {result.stderr[:500]}...", "warning")
        return 0
    except subprocess.CalledProcessError as e:
        print_status(f"BLASTP failed: {e}", "error")
        print_status(f"STDERR: {e.stderr}", "error")
        return 1
    except Exception as e:
        print_status(f"Unexpected error running BLASTP: {e}", "error")
        return 1

def parse_blast_output(xml_path):
    """Parse BLASTP XML output to extract matches."""
    try:
        if not os.path.exists(xml_path):
            print_status(f"BLAST output file not found: {xml_path}", "error")
            return []
        tree = ET.parse(xml_path)
        root = tree.getroot()
        matches = []

        for hit in root.findall(".//Hit"):
            hit_id = hit.findtext("Hit_id") or "Unknown"
            hit_def = hit.findtext("Hit_def") or "No description"
            hsp = hit.find(".//Hsp")
            if hsp is None:
                continue
            identity = int(hsp.findtext("Hsp_identity") or 0)
            align_len = int(hsp.findtext("Hsp_align-len") or 1)
            evalue = float(hsp.findtext("Hsp_evalue") or 1)
            bit_score = float(hsp.findtext("Hsp_bit-score") or 0)
            identity_percent = round((identity / align_len) * 100, 2)

            if evalue <= 1e-2 and identity_percent >= 30:
                matches.append({
                    "subject_id": hit_id,
                    "subject_def": hit_def,
                    "identity_percent": identity_percent,
                    "evalue": evalue,
                    "bit_score": bit_score
                })

        print_status(f"Parsed {len(matches)} matches from {xml_path}", "info")
        return matches
    except Exception as e:
        print_status(f"Error parsing BLAST output {xml_path}: {e}", "error")
        return []

def save_output(output_file, results):
    """Save cross-reactivity results to a CSV file."""
    try:
        os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
        data = []
        for protein_id, protein_name, gut_matches in results:
            matches_str = ""
            if gut_matches:
                matches_str = ";".join(
                    f"{m['subject_id']} ({m['subject_def']}) - Identity: {m['identity_percent']}% | E-value: {m['evalue']} | Bit-score: {m['bit_score']}"
                    for m in gut_matches
                )
            data.append({
                "Sequence_ID": protein_id,
                "Description": protein_name,
                "Cross_Reactivity_Detected": bool(gut_matches),
                "Matches": matches_str
            })
        df = pd.DataFrame(data)
        df.to_csv(output_file, index=False)
        print_status(f"Results saved to {output_file}", "success")
    except Exception as e:
        print_status(f"Error saving output to {output_file}: {e}", "error")
        raise

def cleanup_files(temp_files):
    """Clean up temporary files."""
    for file_path in temp_files:
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                print_status(f"Cleaned up {file_path}", "info")
            except Exception as e:
                print_status(f"Error cleaning up {file_path}: {e}", "warning")

def run_gutflora(input_fasta, output_file, gut_db=GUT_DB):
    """
    Check cross-reactivity with gut flora proteins using BLASTP.

    Args:
        input_fasta (str): Path to input FASTA file.
        output_file (str): Path to output results file or directory.
        gut_db (Path): Path to gut microbiome BLAST database.

    Returns:
        int: Exit status code (0 if successful).
    """
    with tempfile.TemporaryDirectory(dir=ROOT_DIR) as temp_dir:
        temp_dir = Path(temp_dir)
        temp_files = []
        results = []

        try:
            # Validate input
            if not os.path.isfile(input_fasta):
                print_status(f"Input FASTA file not found: {input_fasta}", "error")
                raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

            # Create gutflora subfolder
            output_path = Path(output_file)
            tool_dir = output_path.parent / "gutflora"
            tool_dir.mkdir(parents=True, exist_ok=True)

            # Set combined output file
            combined_output = tool_dir / "combined.csv"

            # Verify database
            if not os.path.exists(gut_db.with_suffix(".phr")):
                print_status(f"Gut microbiome database not found: {gut_db}", "error")
                return 1

            # Parse all sequences
            print_status(f"Parsing FASTA file: {input_fasta}", "info")
            sequences = parse_fasta(input_fasta)
            print_status(f"Found {len(sequences)} sequences in {input_fasta}", "info")

            # Process each sequence
            for idx, (prot_id, prot_desc, prot_seq) in enumerate(sequences, 1):
                print_status(f"\nProcessing sequence {idx}: {prot_id}", "info")

                # Create unique temporary files
                temp_query = f"temp_query_{idx}.fasta"
                blast_gut_output = f"blast_result_gut_{idx}.xml"
                temp_files.extend([str(temp_dir / temp_query), str(temp_dir / blast_gut_output)])

                # Write temporary FASTA
                write_temp_fasta(prot_id, prot_desc, prot_seq, temp_dir / temp_query)

                # Verify file exists
                if not os.path.exists(temp_dir / temp_query):
                    print_status(f"Temporary FASTA file not created: {temp_dir / temp_query}", "error")
                    return 1

                # Run BLASTP
                print_status("Running BLASTP against gut microbiome database...", "info")
                status = run_blastp(gut_db, temp_query, blast_gut_output, temp_dir)
                if status != 0:
                    print_status(f"Gut microbiome BLASTP failed for sequence {idx}: {prot_id}", "error")
                    return status

                # Parse BLAST output
                gut_matches = parse_blast_output(temp_dir / blast_gut_output)
                if gut_matches:
                    print_status(f"Found {len(gut_matches)} matches for {prot_id}: {gut_matches}", "info")
                else:
                    print_status(f"No matches found for {prot_id}", "warning")

                # Store results
                results.append((prot_id, prot_desc, gut_matches))

                # Save individual CSV
                sanitized_id = "".join(c for c in prot_id if c.isalnum() or c in ('.', '_', '-'))
                individual_output = tool_dir / f"{sanitized_id}_gutflora.csv"
                save_output(individual_output, [(prot_id, prot_desc, gut_matches)])

            # Save combined results
            save_output(combined_output, results)

            # Check for cross-reactivity
            if any(gut_matches for _, _, gut_matches in results):
                print_status(f"Cross-reactivity detected. Check: {combined_output}", "warning")
            else:
                print_status("No significant cross-reactivity found.", "info")

            return 0

        except Exception as e:
            print_status(f"Error in gutflora analysis: {e}", "error")
            return 1
        finally:
            cleanup_files(temp_files)

def main():
    """Command-line interface for gut flora cross-reactivity analysis."""
    parser = argparse.ArgumentParser(description="Check cross-reactivity with gut flora proteins.")
    parser.add_argument("-i", "--input", required=True, help="Path to input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Path to output results file or directory")
    parser.add_argument("--gut-db", default=str(GUT_DB), help="Path to gut microbiome BLAST database")

    args = parser.parse_args()
    status = run_gutflora(args.input, args.output, Path(args.gut_db))
    exit(status)

if __name__ == "__main__":
    main()
