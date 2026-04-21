import os
import subprocess
import argparse
import tempfile
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML


# Define default database path (absolute path to your DB folder)
SWISSPROT_DB = Path(__file__).resolve().parent / "human_db" / "human_protein_db"

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

def validate_fasta_file(fasta_path):
    """Validate FASTA file format and content."""
    try:
        with open(fasta_path, "r") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
        if not sequences:
            raise ValueError("FASTA file is empty or contains no valid sequences")
        return sequences
    except Exception as e:
        raise ValueError(f"Invalid FASTA file {fasta_path}: {e}")

def validate_paths(input_fasta, output_file, db_path, output_dir):
    """Validate input, output, and database paths."""
    input_fasta = Path(input_fasta).resolve()
    output_file = Path(output_file).resolve() if output_file else None
    db_path = Path(db_path).resolve()

    # Set default output_dir if not provided
    if output_dir is None:
        output_dir = output_file.parent if output_file else Path.cwd()
    else:
        output_dir = Path(output_dir).resolve()

    # Always use human_blast subdirectory
    output_dir = output_dir / "human_blast"

    # Set combined output file to combined.txt inside human_blast
    output_file = output_dir / "combined.txt"

    # Check input FASTA
    if not input_fasta.is_file():
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    if not os.access(input_fasta, os.R_OK):
        raise PermissionError(f"No read permission for FASTA file: {input_fasta}")
    validate_fasta_file(input_fasta)
    
    # Check if required BLAST DB files exist
    required_exts = [".phr", ".pin", ".psq"]
    missing_files = [str(db_path.with_suffix(ext)) for ext in required_exts if not db_path.with_suffix(ext).is_file()]
    if missing_files:
        raise FileNotFoundError(f"Missing BLAST DB files: {missing_files}")
        
    if not os.access(db_path.parent, os.R_OK):
        raise PermissionError(f"No read permission for database directory: {db_path.parent}")

    # Check output directory
    if not output_dir.exists():
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise PermissionError(f"Cannot create output directory {output_dir}: {e}")
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"No write permission for output directory: {output_dir}")
    return input_fasta, output_file, db_path, output_dir

def check_blastp():
    """Check if blastp is available and executable."""
    try:
        blastp_path = subprocess.run(
            ["which", "blastp"], capture_output=True, text=True, check=True
        ).stdout.strip()
        if not blastp_path:
            raise FileNotFoundError("blastp executable not found in PATH")
        if not os.access(blastp_path, os.X_OK):
            raise PermissionError(f"blastp executable is not executable: {blastp_path}")
        return blastp_path
    except subprocess.CalledProcessError as e:
        raise FileNotFoundError(f"blastp not found: {e}")

def parse_blast_xml(blast_out, identity_thresh):
    """Parse BLAST XML output."""
    try:
        with open(blast_out, "r") as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            results = []
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        identity = (hsp.identities / hsp.align_length) * 100
                        if identity >= identity_thresh:
                            results.append({
                                "qseqid": record.query,
                                "sseqid": alignment.hit_id,
                                "identity": identity,
                                "length": hsp.align_length,
                                "mismatch": hsp.align_length - hsp.identities,
                                "gapopen": hsp.gaps,
                                "qstart": hsp.query_start,
                                "qend": hsp.query_end,
                                "sstart": hsp.sbjct_start,
                                "send": hsp.sbjct_end,
                                "evalue": hsp.expect,
                                "bitscore": hsp.score
                            })
        df = pd.DataFrame(results)
        return df
    except Exception as e:
        print_status(f"Failed to parse BLAST XML output: {e}", "error")
        raise

def run_human(fasta_file=None, input_fasta=None, evalue_thresh=1e-5, identity_thresh=30.0, output_file=None, outfmt=6, output_dir=None, batch_idx=1, combined_output_name=None, individual_suffix="_human_blast"):
    """
    Run BLASTP against a database to identify protein homology.

    Args:
        fasta_file (str, optional): Path to the input FASTA file.
        input_fasta (str, optional): Deprecated alias for fasta_file.
        evalue_thresh (float): E-value threshold for BLASTP.
        identity_thresh (float): Minimum percent identity for filtering BLAST hits.
        output_file (str, optional): Path to the combined output TXT file (will be overridden to combined.txt in human_blast).
        outfmt (int): BLAST output format (6 for tabular, 5 for XML).
        output_dir (str, optional): Directory for output files; defaults to output_file's parent or cwd.
        batch_idx (int): Batch index (unused since combined output name is fixed).
        combined_output_name (str, optional): Ignored; combined output is always combined.txt.
        individual_suffix (str, optional): Suffix for individual output files (default: _human_blast).

    Returns:
        int: Exit status code (0 if successful, 1 if error).
    """
    try:
        # Handle deprecated input_fasta
        if input_fasta is not None:
            fasta_file = input_fasta
        if fasta_file is None:
            raise ValueError("No FASTA file provided (expected 'fasta_file' or 'input_fasta')")

        # Validate paths and FASTA content
        input_fasta, output_file, db_path, output_dir = validate_paths(fasta_file, output_file, SWISSPROT_DB, output_dir)

        # Check blastp availability
        blastp_path = check_blastp()

        # Parse and validate FASTA
        print_status("Parsing FASTA...", "info")
        sequences = validate_fasta_file(input_fasta)
        print_status(f"Total sequences found: {len(sequences)}", "info")

        if not sequences:
            print_status("No sequences found in the FASTA file", "warning")
            return 0

        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_input = Path(temp_dir) / "temp_query.fasta"
            blast_out = Path(temp_dir) / ("blast_results.xml" if outfmt == 5 else "blast_results.tsv")

            # Write temporary FASTA
            print_status(f"Writing temporary FASTA to {temp_input}...", "info")
            SeqIO.write(sequences, temp_input, "fasta")
           
            # Construct and run BLASTP command
            blast_cmd = [
                blastp_path,
                "-query", str(temp_input),
                "-db", str(db_path),
                "-outfmt", str(outfmt),
                "-evalue", str(evalue_thresh),
                "-out", str(blast_out)
            ]
            if outfmt == 5:
                blast_cmd.extend(["-max_target_seqs", "10"])
            print_status(f"Running BLASTP: {' '.join(blast_cmd)}", "info")

            try:
                result = subprocess.run(
                    blast_cmd,
                    check=True,
                    capture_output=True,
                    text=True
                )
                if result.stderr:
                    print_status(f"BLASTP stderr: {result.stderr}", "warning")
                    # Optionally consider this a warning, not an error
            except subprocess.CalledProcessError as e:
                print_status(f"BLASTP failed with exit code {e.returncode}: {e.stderr}", "error")
                return 1

            # Read and filter BLAST output
            print_status("Extracting BLAST results...", "info")
            try:
                if outfmt == 5:
                    df = parse_blast_xml(blast_out, identity_thresh)
                else:
                    col_names = [
                        "qseqid", "sseqid", "identity", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
                    ]
                    df = pd.read_csv(blast_out, sep="\t", names=col_names)
                    df = df[df["identity"] >= identity_thresh]
                
            except Exception as e:
                print_status(f"Failed to process BLAST output: {e}", "error")
                
                return 1

            print_status(f"Found {len(df)} hits meeting identity threshold (>= {identity_thresh}%)", "info")
            

            # Save results
            os.makedirs(output_dir, exist_ok=True)
            excel_file = output_dir / "combined.xlsx"

            try:
                # Save combined TXT (tab-separated)
                if not df.empty:
                    df.to_csv(output_file, sep="\t", index=False)
                    print_status(f"Combined output saved to: {output_file} with {len(df)} records", "success")
                    
                else:
                    with open(output_file, "w") as f:
                        f.write("No BLAST hits found for any sequence\n")
                    print_status(f"Combined output saved to: {output_file} with no hits", "warning")

                # Save combined Excel
                if not df.empty:
                    df.to_excel(excel_file, index=False)
                    print_status(f"Combined Excel output saved to: {excel_file}", "success")

                else:
                    pd.DataFrame(columns=["qseqid", "sseqid", "identity", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]).to_excel(excel_file, index=False)
                    print_status(f"Combined Excel output saved to: {excel_file} with no hits", "warning")

                # Save individual TXT files for each sequence
                sequence_ids = [seq.id for seq in sequences]
                for qseqid in sequence_ids:
                    safe_qseqid = "".join(c for c in qseqid if c.isalnum() or c in "_-")
                    individual_file = output_dir / f"{safe_qseqid}{individual_suffix}.txt"
                    individual_df = df[df["qseqid"] == qseqid] if not df.empty else pd.DataFrame()
                    try:
                        if not individual_df.empty:
                            individual_df.to_csv(individual_file, sep="\t", index=False)
                            print_status(f"Individual output saved to: {individual_file} with {len(individual_df)} records", "success")
                            
                        else:
                            with open(individual_file, "w") as f:
                                f.write(f"No BLAST hits found for sequence {qseqid}\n")
                            print_status(f"Individual output saved to: {individual_file} with no hits", "warning")
                            
                    except Exception as e:
                        print_status(f"Failed to save individual file {individual_file}: {e}", "error")
                        
                        continue

            except Exception as e:
                print_status(f"Failed to save results: {e}", "error")
            
                return 1

            # Verify output files
            all_files = [output_file, excel_file] + [output_dir / f"{''.join(c for c in seq.id if c.isalnum() or c in '_-')}{individual_suffix}.txt" for seq in sequences]
            for file_path in all_files:
                if file_path.exists() and file_path.stat().st_size > 0:
                    print_status(f"Output file verified: {file_path} ({file_path.stat().st_size} bytes)", "success")
                else:
                    print_status(f"Output file is empty or not found: {file_path}", "warning")
                    
            return 0

    except Exception as e:
        print_status(f"Unexpected error running BLASTP: {e}", "error")
        
        return 1


