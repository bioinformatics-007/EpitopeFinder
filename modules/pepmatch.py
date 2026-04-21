import subprocess
import os
import pandas as pd
import json
import sys
import argparse
from pathlib import Path
from Bio import SeqIO

# --- CONFIGURATION ---
current_path = Path(__file__).resolve()
while current_path.name != "Vaxelan_2_0" and current_path != current_path.parent:
    current_path = current_path.parent
ROOT_DIR = current_path

PEPMATCH_EXECUTABLE = ROOT_DIR / "tools" / "IEDB_NG_PEPMATCH-0.1.2-beta" / "ng_pepmatch-0.1.2-beta" / "src" / "match.py"
PEPMATCH_PROTEOMES_PATH = ROOT_DIR / "tools" / "iedb_proteomes" / "proteomes"
proteome_keywords = {'human': '9606', 'mouse': '10090'}

def print_status(msg, status="info"):
    colors = {"info": "\033[94m", "success": "\033[92m", "warning": "\033[93m", "error": "\033[91m"}
    print(f"{colors.get(status, '')}{msg}\033[0m")

def tile_sequence(input_fasta, tiled_fasta, window_size=9):
    tiled_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq).upper()
        for i in range(len(seq) - window_size + 1):
            sub_seq = seq[i:i + window_size]
            # Use the sequence as the ID so PepMatch reports it in 'query_sequence'
            new_id = f"{sub_seq}" 
            tiled_records.append(f">{new_id}\n{sub_seq}")
    with open(tiled_fasta, "w") as f:
        f.write("\n".join(tiled_records))
    return len(tiled_records)

def run_pepmatch_organism(fasta_file, output_dir, organism, taxid, mismatch):
    output_prefix = os.path.join(output_dir, f"pepmatch_{organism}")
    json_out = f"{output_prefix}.json"
   
    command = [
        'python', str(PEPMATCH_EXECUTABLE),
        '--fasta', str(fasta_file),
        '--proteome', str(taxid),
        '--output-prefix', output_prefix,
        '--mismatch', str(mismatch)
    ]
    env = os.environ.copy()
    env['PEPMATCH_PROTEOMES_PATH'] = str(PEPMATCH_PROTEOMES_PATH)
    
    try:
        subprocess.run(command, env=env, capture_output=True, text=True, timeout=600)
        if not os.path.exists(json_out): 
            return False, None
            
        with open(json_out, 'r') as f:
            data = json.load(f)
        
        if 'table_data' in data and data['table_data']:
            df = pd.DataFrame(data['table_data'], columns=data['table_columns'])
            
            # --- IMPROVED MAPPING: Matching the fields in your screenshot ---
            rename_map = {
                'query_id': 'input sequence',
                'peptide': 'matched sequence',
                'protein_id': 'protein id',
                'protein_name': 'protein name',
                'gene': 'gene',
                'mismatches': 'mismatches',
                'start_pos': 'start',
                'end_pos': 'end',
                'taxid': 'taxid'
            }
            df = df.rename(columns=rename_map)
            return True, df
        return True, None
    except Exception as e:
        print_status(f"Error: {e}", "error")
        return False, None

def run_pepmatch(input_fasta, output_dir=".", output_file=None, seq_id=None, mismatch=1, length=9):
    os.makedirs(output_dir, exist_ok=True)
   
    tiled_fasta = os.path.join(output_dir, "tiled_temp.fasta")
    tile_count = tile_sequence(input_fasta, tiled_fasta, window_size=length)
    print_status(f"Generated {tile_count} peptides for matching.")
    
    combined_dfs = []
    for org, taxid in proteome_keywords.items():
        db_path = PEPMATCH_PROTEOMES_PATH / f"{taxid}.db"
        if not db_path.exists(): 
            continue
       
        success, df = run_pepmatch_organism(tiled_fasta, output_dir, org, taxid, mismatch)
        if success and df is not None:
            df['organism'] = 'Homo sapiens' if taxid == '9606' else 'Mus musculus'
            combined_dfs.append(df)
    
    if combined_dfs:
        final_df = pd.concat(combined_dfs, ignore_index=True)
        
        # --- EXACT COLUMN ORDERING AS PER YOUR IMAGE ---
        desired_order = [
            'input sequence',
            'matched sequence',
            'protein id',
            'protein name',
            'organism',
            'taxid',
            'gene',
            'start',
            'end',
            'mismatches'
        ]
        
        # Filter only existing columns and maintain order
        existing_cols = [col for col in desired_order if col in final_df.columns]
        final_df = final_df[existing_cols]
        
        if output_file:
            final_df.to_csv(output_file, index=False)
            print_status(f"Success! Results saved to: {output_file}", "success")
        else:
            print(final_df.to_string(index=False))
    else:
        print_status("No matches found.", "warning")
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-d", "--output_dir", default=".")
    parser.add_argument("-o", "--output_file", default="pepmatch_results.csv")
    parser.add_argument("--mismatch", type=int, default=1)
    parser.add_argument("--length", type=int, default=9)
    args = parser.parse_args()
   
    sys.exit(run_pepmatch(args.input, args.output_dir, args.output_file, 
                          mismatch=args.mismatch, length=args.length))

if __name__ == "__main__":
    main()
