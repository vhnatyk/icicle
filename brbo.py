def bit_reverse(n, num_bits):
    """Return the bit-reversed value of n given a total of num_bits."""
    result = 0
    for _ in range(num_bits):
        result = (result << 1) | (n & 1)
        n >>= 1
    return result

# def ntt_butterfly_rbo(a, b, w, mod, log_n):
#     """NTT butterfly operation with RBO."""
#     i = bit_reverse(a, log_n)
#     j = bit_reverse(b, log_n)
#     A = (a[i] + a[j]) % mod
#     B = (a[i] - a[j]) * w % mod
#     a[i], a[j] = A, B

def ntt_butterfly_rbo(a, i, j, w, mod, log_n):
    """NTT butterfly operation with RBO."""
    i_rbo = bit_reverse(i, log_n)
    j_rbo = bit_reverse(j, log_n)
    A = (a[i_rbo] + a[j_rbo]) % mod
    B = (a[i_rbo] - a[j_rbo]) * w % mod
    a[i_rbo], a[j_rbo] = A, B


def ntt_radix2_rbo_inplace(a, w, mod):
    """In-place radix-2 NTT using Cooley-Tukey algorithm with integrated RBO in butterfly."""
    n = len(a)
    log_n = (n - 1).bit_length()

    for m in range(1, log_n + 1):
        m_len = 2**m
        half_m_len = m_len // 2
        for k in range(0, n, m_len):
            w_m = 1  # twiddle factor
            for j in range(half_m_len):
                ntt_butterfly_rbo(a, k + j, k + j + half_m_len, w_m, mod, log_n)
                w_m = w_m * w % mod

    return a

# Tests
def test_bit_reverse():
    assert bit_reverse(1, 3) == 4
    assert bit_reverse(2, 3) == 2
    assert bit_reverse(3, 3) == 6
    assert bit_reverse(4, 3) == 1
    print("bit_reverse tests passed!")

# def test_butterfly():
#     a = [2, 3, 1, 0]
#     w = 3
#     mod = 7
#     log_n = 2
#     ntt_butterfly_rbo(a, 0, 2, w, mod, log_n)
#     assert a == [4, 3, 6, 0]
#     print("ntt_butterfly_rbo tests passed!")

# def test_butterfly():
#     a = [2, 3, 1, 0]
#     w = 3
#     mod = 7
#     log_n = 2
#     ntt_butterfly_rbo(a, 0, 2, w, mod, log_n)
#     assert a == [5, 3, 6, 0]
#     print("ntt_butterfly_rbo tests passed!")

def test_butterfly():
    a = [0, 2, 1, 3]
    w = 1
    mod = 7
    log_n = 2
    ntt_butterfly_rbo(a, 0, 2, w, mod, log_n)
    print(f"a after ntt: {a}")
    assert a == [5, 4, 1, 0]
    print("ntt_butterfly_rbo tests passed!")



def naive_ntt(a, w, mod):
    """Naive (slow) implementation of NTT for testing."""
    n = len(a)
    result = [0] * n
    for k in range(n):
        s = 0
        for j in range(n):
            s += a[j] * pow(w, k * j, mod)
        result[k] = s % mod
    return result

def test_ntt():
    #a = [1, -1, 1, -1]
    a = [1, 0, 1, 0]
    w = 1
    mod = 4
    result = ntt_radix2_rbo_inplace(list(a), w, mod)
    naive_result = naive_ntt(a, w, mod)
    
    print(f"naive ntt {naive_result}!")
    print(f"result ntt {result}!")
    assert result == naive_result
    print("ntt_radix2_rbo_inplace tests passed!")

test_bit_reverse()
test_butterfly()
test_ntt()
