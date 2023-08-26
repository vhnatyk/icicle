def mod_inverse(a, mod):
    return pow(a, -1, mod)

def cooley_tukey_ntt(a, w, mod):
    n = len(a)
    if n <= 1:
        return a
    a0 = a[::2]
    a1 = a[1::2]
    w_sq = pow(w, 2, mod)
    a0 = cooley_tukey_ntt(a0, w_sq, mod)
    a1 = cooley_tukey_ntt(a1, w_sq, mod)
    f = 1
    a_out = [0] * n
    for i in range(n // 2):
        t = f * a1[i] % mod
        a_out[i] = (a0[i] + t) % mod
        a_out[i + n//2] = (a0[i] - t) % mod
        f = f * w % mod
    return a_out

def cooley_tukey_intt(a, w, mod):
    n = len(a)
    n_inv = mod_inverse(n, mod)
    w_inv = mod_inverse(w, mod)
    a = cooley_tukey_ntt(a, w_inv, mod)  # Compute INTT via NTT
    for i in range(n):
        a[i] = a[i] * n_inv % mod
    return a

if __name__ == '__main__':
    mod = 17  # A small prime for example
    w = 3  # This should be a 4th root of unity modulo 17 (3^4 % 17 = 1, 3^2 % 17 != 1)
    a = [1, 1, 16, 16]  # Input array

    print(f"Original array: {a}")

    # Perform NTT
    a_ntt = cooley_tukey_ntt(a, w, mod)
    print(f"After NTT: {a_ntt}")

    # Perform INTT
    a_intt = cooley_tukey_intt(a_ntt, w, mod)
    print(f"After INTT: {a_intt}")

    assert a == a_intt, f"Expected {a}, got {a_intt}"
    print("Test passed!")
