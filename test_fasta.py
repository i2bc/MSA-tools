def check_fasta_lengths(file_path):
    sequences = {}
    sequence_id = ""
    sequence = ""
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith(">"):
                if sequence_id and sequence:
                    sequences[sequence_id] = len(sequence)
                sequence_id = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        if sequence_id and sequence:
            sequences[sequence_id] = len(sequence)

    for seq_id, seq_len in sequences.items():
        print(f"{seq_id}: {seq_len} characters")

    lengths = list(sequences.values())
    if len(set(lengths)) == 1:
        print("All sequences have the same length.")
    else:
        print("Sequences have different lengths.")
        print(f"Lengths: {lengths}")


check_fasta_lengths("./split_1.fasta")
