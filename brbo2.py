import math


def mod_inverse(a, mod):
    return pow(a, -1, mod)

def reverse_bits(number, bit_length):
    # Reverses the bits of `number` up to `bit_length`.
    reversed = 0
    for i in range(0, bit_length):
        if (number >> i) & 1: reversed |= 1 << (bit_length - 1 - i)
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
                V = ret[k] * pow(w, i*j, P) # omegas[g]
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

if __name__ == '__main__':
    mod = 17  # A small prime for example
    w = 3  # This should be a 4th root of unity modulo 17 (3^4 % 17 = 1, 3^2 % 17 != 1)
    a = [1, 1, 16, 16]  # Input array
    a = [1, 1, 1, 1]  # Input array
    # a = [1, 16, 1, 16]  # Input array
    # a = [1, 1, 0, 0]  # Input array

    print(f"Original array: {a}")

    # Perform NTT
    a_ntt = ntt(a.copy(), w, mod)
    print(f"After NTT: {a_ntt}")

    # Perform INTT
    a_intt = intt(a_ntt, w, mod)
    print(f"After INTT: {a_intt}")

    assert a == a_intt, f"Expected {a}, got {a_intt}"
    print("Test passed!")
