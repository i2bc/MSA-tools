def parse_fasta(file):
    """Parse the FASTA file and return a list of dictionaries with headers and sequences."""
    sequences = []
    with open(file, "r") as f:
        header = ""
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and sequence:
                    sequences.append({"header": header, "sequence": sequence})
                header = line
                sequence = ""
            else:
                sequence += line
        if header and sequence:
            sequences.append({"header": header, "sequence": sequence})
    return sequences


def sort_by_best_match(sequences):
    """Get the best match compared to the first sequence in the list."""
    if not sequences:
        return None

    main_header = sequences[0]["header"]
    main_sequence = sequences[0]["sequence"]

    for i in range(1, len(sequences)):
        other_sequence = sequences[i]["sequence"]
        match_count = sum(1 for a, b in zip(main_sequence, other_sequence) if a == b)
        sequences[i]["matches"] = match_count
    sorted_by_match = sorted(sequences[1:], reverse=True, key=lambda x: x["matches"])

    sorted_by_match.insert(0, {"header": main_header, "sequence": main_sequence})

    return sorted_by_match


def write_best_match(sequences, output_file, number_of_sequences=10):
    with open(output_file, "w") as f:
        # Start at 0 to account for main sequence
        counter = 0
        for sequence in sequences:
            f.write(sequence["header"] + "\n")
            f.write(sequence["sequence"] + "\n")
            if counter == number_of_sequences:
                break
            counter += 1


# Example usage:
fasta_file = "./msa_hhfilter.fasta"
sequences = parse_fasta(fasta_file)
write_best_match(sort_by_best_match(sequences), "test.fasta")
