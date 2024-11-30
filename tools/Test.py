import math
from abc import ABC, abstractmethod
from collections import Counter

import numpy
import numpy as np
from numpy.fft import fft
from scipy.special import gammaincc
from scipy.stats import chi2, norm


class Test(ABC):
    def __init__(self, sequence):
        self.sequence = sequence

    @abstractmethod
    def run_test(self):
        pass


class MDRNTest(Test):
    def __init__(self, sequence, L):
        super().__init__(sequence)
        self.L = L

    def run_test(self):
        n = len(self.sequence)
        m = n // self.L  # Число фрагментов разбиения

        if m < 50 or m * (1 / 2) * (self.L ** 2) < 10:
            raise ValueError("Условия применения теста не выполняются.")

        fragments = [self.sequence[i * self.L:(i + 1) * self.L] for i in range(m)]
        fragment_counts = Counter(fragments)
        all_possible_fragments = [f"{i:0{self.L}b}" for i in range(2 ** self.L)]
        expected_count = m * (2 ** (-self.L))
        chi_square = 0

        for fragment in all_possible_fragments:
            observed_count = fragment_counts.get(fragment, 0)
            chi_square += ((observed_count - expected_count) ** 2) / expected_count

        p_value = 1 - chi2.cdf(chi_square, 2 ** self.L - 1)

        return p_value


class MDRPTest(Test):
    def __init__(self, sequence, L):
        super().__init__(sequence)
        self.L = L

    def run_test(self):
        n = len(self.sequence)

        if n < 20 * 2 ** self.L:
            raise ValueError("Условия применения теста не выполняются.")

        # Overlapping fragments for L
        fragments_L = []
        for i in range(n):
            if i + self.L <= n:
                fragment = self.sequence[i:i + self.L]
            else:
                tmp = self.L - len(self.sequence[i:n])
                fragment = self.sequence[i:n] + self.sequence[0:tmp]
            fragments_L.append(fragment)

        frequencies_L = Counter(fragments_L)
        gamma_L = sum((freq - n * 2 ** (-self.L)) ** 2 / (n * 2 ** (-self.L)) for v, freq in frequencies_L.items())

        # Overlapping fragments for L-1
        fragments_L_minus_1 = []
        for i in range(n):
            if i + self.L <= n:
                fragment = self.sequence[i:i + self.L - 1]
            else:
                tmp = self.L - len(self.sequence[i:n]) - 1
                fragment = self.sequence[i:n] + self.sequence[0:tmp]
            fragments_L_minus_1.append(fragment)

        frequencies_L_minus_1 = Counter(fragments_L_minus_1)
        gamma_L_minus_1 = sum(
            (freq - n * 2 ** (-(self.L - 1))) ** 2 / (n * 2 ** (-(self.L - 1))) for v, freq in frequencies_L_minus_1.items())

        # Calculate MDPR statistic
        MDPR_statistic = gamma_L - gamma_L_minus_1

        # Calculate p-value using chi-squared distribution
        P = 1 - chi2.cdf(MDPR_statistic, (2 ** self.L) - (2 ** (self.L - 1)))

        return P


class EmptyBoxesTest(Test):
    def __init__(self, sequence, L):
        super().__init__(sequence)
        self.L = L

    def run_test(self):
        n = len(self.sequence)
        m = n // self.L
        fragments = [self.sequence[i * self.L:(i + 1) * self.L] for i in range(m)]
        fragment_counts = Counter(fragments)

        all_possible_fragments = [f"{i:0{self.L}b}" for i in range(2 ** self.L)]
        mu_0 = all_possible_fragments.__len__() - fragment_counts.__len__()
        SEB = m - (2 ** self.L - mu_0)

        lambd = m / (2 ** self.L)
        lambda_pi = m ** 2 / (2 * 2 ** self.L)
        mu = 2 ** self.L * math.exp(-lambd)
        sigma_sq = 2 ** self.L * math.exp(-lambd) * (1 - (1 + lambd) * math.exp(-lambd))
        p_value = 0
        if lambd <= 1 / 32:
            p_value = chi2.cdf(SEB, 2 * lambda_pi)
        elif 1 / 32 < lambd < 5:
            p_value = norm.cdf(-(mu_0 - mu) / math.sqrt(sigma_sq))

        return p_value


class RunsTest(Test):
    def __init__(self, sequence):
        super().__init__(sequence)

    def run_test(self):
        n = len(self.sequence)

        # Find series lengths for 0's and 1's
        series_lengths = {0: [], 1: []}
        i = 0
        while i < n:
            current_bit = self.sequence[i]
            length = 1
            while i + length < n and self.sequence[i + length] == current_bit:
                length += 1
            if i == 0 or i + length == n or self.sequence[i - 1] != current_bit and self.sequence[i + length] != current_bit:
                series_lengths[int(current_bit)].append(length)
            i += length

        # Calculate theoretical frequencies (mu)
        mu = []
        i = 1
        while True:
            mu_i = (n - i + 3) / (2 ** (i + 2))
            if mu_i < 5:
                break
            mu.append(mu_i)
            i += 1
        k = i - 1  # Maximum length for series

        # Calculate observed frequencies (nu_0 and nu_1)
        nu_0 = [0] * k
        nu_1 = [0] * k

        for length in series_lengths[0]:
            if length <= k:
                nu_0[length - 1] += 1

        for length in series_lengths[1]:
            if length <= k:
                nu_1[length - 1] += 1

        # Calculate Sruns statistic
        Sruns = 0
        for i in range(len(mu)):
            Sruns += (nu_0[i] - mu[i]) ** 2 / mu[i] + (nu_1[i] - mu[i]) ** 2 / mu[i]

        # Calculate p-value using chi-squared distribution
        P = 1 - chi2.cdf(Sruns, 2 * k - 2)

        return P


class LongRunsTest(Test):
    def __init__(self, sequence, L):
        super().__init__(sequence)
        self.L = L

    def run_test(self):
        """
        Note that this description is taken from the NIST documentation [1]
        [1] http://csrc.nist.gov/publications/nistpubs/800-22-rev1a/SP800-22rev1a.pdf
        The focus of the tests is the longest run of ones within M-bit blocks. The purpose of this tests is to determine
        whether the length of the longest run of ones within the tested sequences is consistent with the length of the
        longest run of ones that would be expected in a random sequence. Note that an irregularity in the expected
        length of the longest run of ones implies that there is also an irregularity ub tge expected length of the long
        est run of zeroes. Therefore, only one test is necessary for this statistical tests of randomness
        :param sequence: a binary string
        :return: the p-value from the test
        """
        if len(self.sequence) < 128:
            print("\t", "Not enough data to run test!")
            return -1.0
        elif self.L == 8:
            k = 3
            v_values = [1, 2, 3, 4]
            pik_values = [0.21484375, 0.3671875, 0.23046875, 0.1875]
        elif self.L == 128:
            k = 5
            v_values = [4, 5, 6, 7, 8, 9]
            pik_values = [0.1174035788, 0.242955959, 0.249363483, 0.17517706, 0.102701071, 0.112398847]
        elif self.L == 512:
            k = 5
            v_values = [6, 7, 8, 9, 10, 11]
            pik_values = [0.1170, 0.2460, 0.2523, 0.1755, 0.1015, 0.1077]
        elif self.L == 1000:
            k = 5
            v_values = [7, 8, 9, 10, 11, 12]
            pik_values = [0.1307, 0.2437, 0.2452, 0.1714, 0.1002, 0.1088]
        elif self.L == 10000:
            k = 6
            v_values = [10, 11, 12, 13, 14, 15, 16]
            pik_values = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]

        # Work out the number of blocks, discard the remainder
        # pik = [0.2148, 0.3672, 0.2305, 0.1875]
        num_blocks = math.floor(len(self.sequence) / self.L)
        frequencies = numpy.zeros(k + 1)
        block_start, block_end = 0, self.L
        for i in range(num_blocks):
            # Slice the binary string into a block
            block_data = self.sequence[block_start:block_end]
            # Keep track of the number of ones
            max_run_count, run_count = 0, 0
            for j in range(0, self.L):
                if block_data[j] == '1':
                    run_count += 1
                    max_run_count = max(max_run_count, run_count)
                else:
                    max_run_count = max(max_run_count, run_count)
                    run_count = 0
            max_run_count = max(max_run_count, run_count)
            if max_run_count < v_values[0]:
                frequencies[0] += 1
            for j in range(k):
                if max_run_count == v_values[j]:
                    frequencies[j] += 1
            if max_run_count > v_values[k - 1]:
                frequencies[k] += 1
            block_start += self.L
            block_end += self.L
        # print(frequencies)
        chi_squared = 0
        for i in range(len(frequencies)):
            chi_squared += (pow(frequencies[i] - (num_blocks * pik_values[i]), 2.0)) / (num_blocks * pik_values[i])
        p_val = gammaincc(float(k / 2), float(chi_squared / 2))
        return p_val


class ApproximateEntropyTest(Test):
    def __init__(self, sequence, L):
        super().__init__(sequence)
        self.L = L

    def run_test(self):
        """
        Calculates the approximate entropy for a binary sequence.

        As with the Serial test of Section 2.11, the focus of this test is the frequency of all possible overlapping
        m-bit patterns across the entire sequence. The purpose of the test is to compare the frequency of overlapping
        blocks of two consecutive/adjacent lengths (m and m+1) against the expected result for a random sequence.

        :param sequence: a binary string
        :param L: the length of the pattern (m)
        :return: the P-value
        """
        n = len(self.sequence)

        # Extend the sequence with the first m+1 bits to handle wrapping in the frequency counting
        self.sequence += self.sequence[:self.L + 1]

        # Initialize frequency arrays for patterns of length L and L+1
        vobs_one = np.zeros(2 ** self.L)
        vobs_two = np.zeros(2 ** (self.L + 1))

        # Count the occurrences of each pattern in the sequence
        for i in range(n):
            vobs_one[int(self.sequence[i:i + self.L], 2)] += 1
            vobs_two[int(self.sequence[i:i + self.L + 1], 2)] += 1

        # Calculate the test statistics and p-values
        vobs = [vobs_one, vobs_two]
        sums = np.zeros(2)

        for i in range(2):
            for count in vobs[i]:
                if count > 0:
                    sums[i] += count * math.log(count / n)
        sums /= n

        # Approximate entropy
        ape = sums[0] - sums[1]
        chi_squared = 2.0 * n * (math.log(2) - ape)
        p_val = gammaincc(2 ** (self.L - 1), chi_squared / 2.0)

        return p_val


class ScalarMultiplyTest(Test):
    def __init__(self, sequence, m, K):
        super().__init__(sequence)
        self.m = m
        self.K = K

    def run_test(self):
        # Шаг 1: Работаем с последовательностью как со строкой
        n = len(self.sequence)
        L = self.m * self.K
        M = n // L  # Число фрагментов

        # Разбиение строки на фрагменты длины L
        fragments = [self.sequence[i * L:(i + 1) * L] for i in range(M)]

        # Шаг 2: Вычисление скалярных произведений
        Y = []
        for frag in fragments:
            # Разбиение фрагмента на K подфрагментов длины m
            subfrags = [frag[i * self.m:(i + 1) * self.m] for i in range(self.K)]
            y_frag = []
            for j in range(1, self.K):
                # Преобразование подфрагментов в числа (0 или 1) и вычисление скалярного произведения
                subfrag_0 = np.array([int(bit) for bit in subfrags[0]])
                subfrag_j = np.array([int(bit) for bit in subfrags[j]])
                yj = np.dot(subfrag_0, subfrag_j)
                y_frag.append(yj)
            Y.append(y_frag)

        # Шаг 3: Определение экстремальной статистики
        Ymax = [max(y_frag) for y_frag in Y]

        # Шаг 4: Эмпирические частоты
        max_Ymax = max(Ymax)
        f = [Ymax.count(k) for k in range(max_Ymax + 1)]

        # Шаг 5: Теоретические частоты
        q = []
        for k in range(max_Ymax + 1):
            if k == 0:
                q0 = (2 ** -self.m) * (1 + 2 ** (1 - self.K)) ** self.m
                q.append(q0)
            else:
                qk = 2 ** -self.m * sum([math.comb(self.m, l) * (sum([2 ** -l * math.comb(l, y) for y in range(k + 1)]) ** (self.K - 1) - sum(
                    [2 ** -l * math.comb(l, y) for y in range(k)]) ** (self.K - 1)) for l in range(self.m + 1)])
                q.append(qk)

        # Шаг 6: Статистика теста
        Mq = [M * qk for qk in q]
        S_chi2 = sum([(f[k] - Mq[k]) ** 2 / Mq[k] for k in range(len(f))])

        # Шаг 7: P-значение
        p_value = 1 - chi2.cdf(S_chi2, df=self.m)

        return p_value


class SpectralTest(Test):
    def __init__(self, sequence):
        super().__init__(sequence)

    def run_test(self):
        # Шаг 1: Преобразование последовательности X (0 -> -1, 1 -> 1)
        Y = np.array([2 * int(x) - 1 for x in self.sequence])

        # Шаг 2: Применение БПФ
        S = fft(Y)

        # Шаг 3: Вычисление модуля комплексных чисел
        n = len(self.sequence)
        S_abs = np.abs(S[:n // 2])  # Используем половину значений из-за симметрии БПФ

        # Шаг 4: Определение порога h
        # h = np.sqrt(3 * n)
        h = np.sqrt(2.995732274 * n)

        # Шаг 5: Построение последовательности случайных величин и вычисление S_SC
        u = np.array([1 if s < h else 0 for s in S_abs])
        S_SC = np.sum(u)

        # Шаг 6: Вычисление математического ожидания и дисперсии
        p = 0.95
        mu = n * p / 2
        sigma = np.sqrt(n * p * (1 - p) / 4)

        # Шаг 7: Вычисление P-значения
        P_value = 2 * (1 - norm.cdf(np.abs((S_SC - mu) / sigma)))

        return P_value
