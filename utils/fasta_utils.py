def parse_fasta(fasta_str):
    sequences = []
    entries = fasta_str.split(">")
    for entry in entries:
        if not entry.strip():
            continue
        parts = entry.strip().split("\n", 1)
        header = parts[0].split()[0]
        seq = parts[1].replace("\n", "").strip() if len(parts) > 1 else ""
        if seq:
            sequences.append((header, seq))
    return sequences
