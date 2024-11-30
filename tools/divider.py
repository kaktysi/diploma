import os

if not os.path.exists("../binary_sequences"):
    os.makedirs("../binary_sequences")

binary_file_path = "../sequences/phys.bin"

with open(binary_file_path, "rb") as file:
    chunk_size = 1024 * 1024  # 1 МБ
    chunk_count = 100
    for i in range(chunk_count):
        chunk = file.read(chunk_size)
        # filename = f"../binary_sequences/sequence_phys_{i}.bin"
        #
        # with open(filename, "wb") as chunk_file:
        #     chunk_file.write(chunk)
        binary_sequence = ''.join(format(byte, '08b') for byte in chunk)
        sequence_filename = f"../binary_sequences/sequence_phys_{i}.txt"

        with open(sequence_filename, "w") as sequence_file:
            sequence_file.write(binary_sequence)