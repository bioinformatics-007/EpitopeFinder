def write_fasta(sequences, output_file):
    try:
        with open(output_file, 'w') as f:
            for seq_id, desc, seq in sequences:
                f.write(f">{seq_id} {desc}\n{seq}\n")
        return True
    except Exception as e:
        print(f"Error writing file: {e}")
        return False
