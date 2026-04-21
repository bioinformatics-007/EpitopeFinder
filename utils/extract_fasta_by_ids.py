def extract_sequences_by_ids(fasta_file, id_list):
    """
    Extract sequences from a FASTA file for given sequence IDs.
    Args:
        fasta_file (str): Path to the FASTA file.
        id_list (list): List of sequence IDs to find.
    Returns:
        list: List of (ID, description, sequence) for matching sequences.
    """
    sequences = []
    current_id = None
    current_desc = None
    current_seq = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_id in id_list:
                        sequences.append((current_id, current_desc, ''.join(current_seq)))
                    header = line[1:].split(maxsplit=1)
                    current_id = header[0]
                    current_desc = header[1] if len(header) > 1 else ''
                    current_seq = []
                else:
                    current_seq.append(line)
        
        if current_id and current_id in id_list:
            sequences.append((current_id, current_desc, ''.join(current_seq)))
        
        return sequences
    except FileNotFoundError:
        print(f"Error: FASTA file {fasta_file} not found")
        return []
    except Exception as e:
        print(f"Error parsing FASTA file: {e}")
        return []
