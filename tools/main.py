from tools.Test import (MDRNTest, MDRPTest, EmptyBoxesTest, RunsTest, LongRunsTest, ApproximateEntropyTest,
                        ScalarMultiplyTest, SpectralTest)
from tools.tester import test_binary_sequences
import time

start_time = time.time()

Battery_Belt = "../binary_sequences/sequence_BelT"
Battery_lfsr_5 = "../binary_sequences/sequence_lfsr"
Battery_lfsr_30 = "../binary_sequences/sequence_lfsr(x^30+x^29+x^9+1)"
Battery_phys = "../binary_sequences/sequence_phys"

output_template = "../test_results/phys/test_results_phys_"

file_template = Battery_phys

test_binary_sequences(file_template, MDRNTest, f"{output_template}1.txt", 9)
test_binary_sequences(file_template, MDRPTest, f"{output_template}2.txt", 11)
test_binary_sequences(file_template, EmptyBoxesTest, f"{output_template}3.txt", 18)
test_binary_sequences(file_template, RunsTest, f"{output_template}4.txt")
test_binary_sequences(file_template, LongRunsTest, f"{output_template}5.txt", 128)
test_binary_sequences(file_template, ApproximateEntropyTest, f"{output_template}6.txt", 8)
test_binary_sequences(file_template, ScalarMultiplyTest, f"{output_template}7.txt", 4, 10)
test_binary_sequences(file_template, SpectralTest, f"{output_template}8.txt")

current_time = time.time() - start_time
print(f"Прошло {current_time:.2f} секунд.")