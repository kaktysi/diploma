def read_sequence_from_file(file_path):
    with open(file_path, 'r') as file:
        sequence = file.read().strip()
    return sequence
