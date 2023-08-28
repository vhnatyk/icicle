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


def cooley_tukey_ntt(inp, P, w):
    """Cooley-Tukey NTT algorithm."""
    ret = inp
    N = len(ret)
    bit_length = N.bit_length() - 1

    for i in range(N):
        rev_i = reverse_bits(i, bit_length)
        if rev_i > i:
            ret[i] ^= ret[rev_i]
            ret[rev_i] ^= ret[i]
            ret[i] ^= ret[rev_i]

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


def ntt(a, w, mod):
    return cooley_tukey_ntt(a, mod, w)


def intt(a, w, mod):
    n = len(a)
    n_inv = mod_inverse(n, mod)
    w_inv = mod_inverse(w, mod)
    a = ntt(a, w_inv, mod)  # Compute INTT via NTT
    for i in range(n):
        a[i] = a[i] * n_inv % mod
    return a


def batch_ntt(a, n1, w, mod):
    print("#####batch######")
    for i in range(0, len(a), n1):
        tmp = a[i:i+n1]
        print(f"{tmp}")
        a[i:i+n1] = ntt(tmp, w, mod)


def bailey_ntt(a, n1, w, mod):
    print("#####bailey######")

    n2 = len(a) // n1

    print(f"n1:{n1} n2: {n2}")

    # util.in_place_transpose(a, n1, n2)

    batch_ntt(a, n1, w, mod)

    for j in range(n2):
        for i in range(n1):
            a[i * n1 + j] = ((pow(w, i*j, mod) % mod) * a[i * n1 + j]) % mod

    util.in_place_transpose(a, n2, n1)

    batch_ntt(a, n2, w, mod)

    # util.in_place_transpose(a, n1, n2)


if __name__ == '__main__':
    mod = 17  # A small prime for example
    # This should be a 4th root of unity modulo 17 (3^4 % 17 = 1, 3^2 % 17 != 1)
    w = 3
    a = [1, 1, 16, 16]  # Input array
    a = [1, 1, 1, 1]  # Input array
    # a = [1, 16, 1, 16]  # Input array
    # a = [1, 1, 0, 0]  # Input array

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
    a2 = [1, 1, 1, 1]  # Input array
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
    # print("\n###Bailey ntt###\n")

    a1 = [1, 16,]  # Input array
    a2 = [1, 1,]  # Input array
    # a = [1, 16, 1, 16]  # Input array
    # a = [1, 1, 0, 0]  # Input array

    a1a2 = a1+a2
    # print(f"Original array: {a1a2}")

    a1a2_copy = a1a2.copy()

    ntt(a1a2_copy, w, mod)
    print(f"ntt array: {a1a2_copy}")

    bailey_ntt(a1a2, len(a1a2)//2, w, mod)

    assert a1a2 == a1a2_copy, f"Expected {a1a2_copy}, got {a1a2}"
    print("Test passed!")
