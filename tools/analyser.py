from math import comb


#тестирование
nk = [73, 23, 4, 0, 0, 0, 0, 0]


def E_n_k(n, k, p):
    return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))


def P_ear(n, p):
    sum_result = 0
    # Summing from k=0 to k=1
    for k in range(1):
        n_k = nk[k]  # assuming η_k is k for this example
        E_n_k_value = 100 * E_n_k(n, k, p)  # calculate E(η_k)
        if E_n_k_value != 0:  # To avoid division by zero
            term = ((n_k - E_n_k_value) ** 2) / E_n_k_value
            sum_result += term

    # Summing for k = 2 to n
    n_k = 0
    E_n_k_total = 0
    for k in range(1, n + 1):
        n_k += k
        E_n_k_total += 100 * E_n_k(n, k, p)

    if E_n_k_total != 0:
        term = ((n_k - E_n_k_total) ** 2) / E_n_k_total
        sum_result += term

    return sum_result


n = 8
p = 0.05

for k in range(n + 1):
    result = 100 * E_n_k(n, k, p)
    print(f"E(n({k})) = {result}")

print(f"P_ear = {P_ear(n, p)}")
