import math
import util


def mod_inverse(a, mod):
    return pow(a, -1, mod)


def reverse_bits(number, bit_length):
    # Reverses the bits of `number` up to `bit_length`.
    reversed = 0
    for i in range(0, bit_length):
        if (number >> i) & 1:
            reversed |= 1 << (bit_length - 1 - i)
    return reversed

# def cooley_tukey_ntt(a, w, mod):
#     n = len(a)
#     if n <= 1:
#         return a
#     a0 = a[::2]
#     a1 = a[1::2]
#     w_sq = pow(w, 2, mod)
#     a0 = ntt(a0, w_sq, mod)
#     a1 = ntt(a1, w_sq, mod)
#     f = 1
#     a_out = [0] * n
#     for i in range(n // 2):
#         t = f * a1[i] % mod
#         a_out[i] = (a0[i] + t) % mod
#         a_out[i + n//2] = (a0[i] - t) % mod
#         f = f * w % mod
#     return a_out


# def get_root_of_unity(order: int, g, mod) -> int:
#     """
#     Returns a root of unity of order "order"
#     """
#     print(f"({mod} - 1) % {order} == 0")
#     assert (mod - 1) % order == 0
#     return pow(g, (mod - 1) // order, mod)

def get_root_of_unity_o(order: int, g, modolus_p) -> int:
    assert (modolus_p - 1) % order == 0
    return pow(g, (modolus_p - 1) // order, modolus_p)


def get_root_of_unity(n, g, mod):
    k = int(math.floor(math.log2(n)))
    return get_root_of_unity_o(int(pow(2, k+1)), g, mod)


def cooley_tukey_ntt_nn(inp, P, w):
    """Cooley-Tukey NTT algorithm."""
    ret = inp
    N = len(ret)
    # rev_rbo2(ret, N)

    M = 2
    iters = int(math.log2(N))
    for l in range(iters):
        for i in range(0, N, M):
            g = 0
            for j in range(0, M // 2):
                k = i + j + (M // 2)
                U = ret[i + j]
                V = ret[k] * pow(w, i*j, P)  # omegas[g]
                ret[i + j] = (U + V) % P
                ret[k] = (U - V) % P
                g = g + N // M
        M = M * 2

    return ret


def rev_rbo2(ret, N):
    bit_length = N.bit_length() - 1
    rev_rbo(ret, N, bit_length)


def rev_rbo(ret, N, bit_length):
    for i in range(N):
        rev_i = reverse_bits(i, bit_length)
        if rev_i > i:
            ret[i] ^= ret[rev_i]
            ret[rev_i] ^= ret[i]
            ret[i] ^= ret[rev_i]


def ntt(a, w, mod):
    # rev_rbo2(a, len(a))
    return cooley_tukey_ntt_nn(a, mod, w)


def intt(a, w, mod):
    n = len(a)
    n_inv = mod_inverse(n, mod)
    w_inv = mod_inverse(w, mod)
    a = ntt(a, w_inv, mod)  # Compute INTT via NTT
    for i in range(n):
        a[i] = a[i] * n_inv % mod
    return a


def batch_ntt(a, n1, w, mod):
    # print("#####batch######")
    for i in range(0, len(a), n1):
        tmp = a[i:i+n1]
        # print(f"{tmp}")
        ntt(tmp, w, mod)
        a[i:i+n1] = tmp


def matrix_print(a, n1):
    for i in range(0, len(a), n1):
        tmp = a[i:i+n1]
        print(f"{tmp}")


def batch_rbo(a, n1):
    # print("#####batch######")
    for i in range(0, len(a), n1):
        tmp = a[i:i+n1]
        # print(f"{tmp}")
        rev_rbo2(tmp, n1)
        a[i:i+n1] = tmp


def bailey_ntt(a, n1, g, mod, is_print=False):
    if is_print:
        print("#####bailey######")

    a_len = len(a)
    n2 = a_len // n1
    w_n1 = get_root_of_unity(n1, g, mod)
    w_n2 = get_root_of_unity(n2, g, mod)
    w_full = get_root_of_unity(a_len, g, mod)

    if is_print:
        print(f"input      : {a} n1: {n1} n2: {n2}")

    util.in_place_transpose(a, n1, n2)

    if is_print:
        print(f"transpose  : {a} n1: {n1} n2: {n2}")

    batch_rbo(a, n1)
    if is_print:
        print(f"rbo1       : {a} n1: {n1} n2: {n2}")
    
    batch_ntt(a, n1, w_n1, mod)
    if is_print:
        print(f"batch1     : {a} n1: {n1} n2: {n2}")

    for j in range(n2):
        for i in range(n1):
            a[i * n1 + j] = ((pow(w_full, i * j, mod) % mod)
                             * a[i * n1 + j] % mod) % mod
    if is_print:
        print(f"* w^(i*j)  : {a} n1: {n1} n2: {n2}")

    util.in_place_transpose(a, n2, n1)
    if is_print:
        print(f"transpose 2: {a} n1: {n1} n2: {n2}")

    batch_rbo(a, n2)
    if is_print:
        print(f"rbo2       : {a} n1: {n1} n2: {n2}")

    batch_ntt(a, n2, w_n2, mod)
    if is_print:
        print(f"batch2     : {a} n1: {n1} n2: {n2}")

    util.in_place_transpose(a, n1, n2)
    if is_print:
        print(f"transpose 3: {a} n1: {n1} n2: {n2}")


if __name__ == '__main__':
    mod = 521  # A small prime for example
    mod = 1048589  # A small prime for example
    mod = 7340033  # A small prime for example
    # This should be a 4th root of unity modulo 17 (2^4 % 17 = 1, 3^2 % 17 != 1)
    # g = 3
    g = 5

    a = [1, 1, 16, 16]  # Input array
    a = [1, 1, 1, 1]  # Input array
    # a = [1, 16, 1, 16]  # Input array
    # a = [1, 1, 0, 0]  # Input array

    w = get_root_of_unity(len(a), g, mod)

    # print(f"Original array: {a}")

    # Perform NTT
    a_ntt = ntt(a.copy(), w, mod)
    # print(f"After NTT: {a_ntt}")

    # Perform INTT
    a_intt = intt(a_ntt, w, mod)
    # print(f"After INTT: {a_intt}")

    assert a == a_intt, f"Expected {a}, got {a_intt}"
    # print("Test passed!")

    #################
    # print("\n###Batch ntt###\n")

    a1 = [1, 1, 16, 16]  # Input array
    a2 = [1, 1, 0, 0]  # Input array
    # a = [1, 16, 1, 16]  # Input array
    # a = [1, 1, 0, 0]  # Input array

    a1a2 = a1+a2
    # print(f"Original array: {a1a2}")

    # Perform NTT
    ntt(a1, w, mod)
    ntt(a2, w, mod)

    a1a2_batch = a1a2.copy()
    batch_ntt(a1a2_batch, len(a1), w, mod)

    a1a2 = a1+a2
    assert a1a2 == a1a2_batch, f"Expected {a1a2}, got {a1a2_batch}"
    # print("Test passed!")

    #################
    print("\n###Bailey ntt###\n")

    a1 = [1, 16,]  # Input array
    a2 = [0, 1,]  # Input array
    # a = [1, 16, 1, 16]  # Input array
    # a = [1, 1, 0, 0]  # Input array
    # a = [1, 0, 1, 0]  # Input array
    # a = [1, 0, 0, 0]  # Input array

    # a1a2 = [i % mod for i in range(256)]
    a1a2 = [(i+1) % mod for i in range(16)] * 16

    a1a2_copy = a1a2.copy()
    # print(f"Original array: {a1a2_copy}")

    # rev_rbo2(a1a2_copy, 2)
    w = get_root_of_unity(len(a1a2), g, mod)
    
    rev_rbo2(a1a2_copy, len (a1a2_copy))
    
    ntt(a1a2_copy, w, mod)
    # print(f"ntt array: {a1a2_copy}")

    n1 = 1 << int(math.log2(len(a1a2))//2)

 
    bailey_ntt(a1a2, n1, g, mod, False)
    rev_rbo2(a1a2, len (a1a2))
    
    
    if a1a2 != a1a2_copy:
        print("!!!fail!!!1")
    print("Expected")
    matrix_print(a1a2_copy, n1)
    print("got")
    matrix_print(a1a2, n1)

    if a1a2 != a1a2_copy:

        assert False

    ####################
    print("\n###Bailey ntt 2###\n")

    a1 = [1, 16, 1, 16]  # Input array
    a2 = [1, 1, 0, 0]  # Input array
    a3 = [1, 0, 1, 0]  # Input array
    a4 = [1, 0, 0, 0]  # Input array

    a1a2 = a1+a2+a3+a4
    print(f"Original array: {a1a2}")

    a1a2_copy = a1a2.copy()

    ntt(a1a2_copy, w, mod)
    print(f"ntt array: {a1a2_copy}")
    bailey_ntt(a1a2, int(math.sqrt(len(a1a2))), w, mod)

    print("Expected")
    matrix_print(a1a2, n1)
    print("got")
    if a1a2 != a1a2_copy:
        matrix_print(a1a2_copy, n1)
