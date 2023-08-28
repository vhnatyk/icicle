




def in_place_transpose(arr, col_cnt, row_len):
    def A(i, j):
        return i * row_len + j

    n = col_cnt * row_len
    visited = [False] * n

    for start in range(n):
        if visited[start]:
            continue

        cycle = [start]
        i = start
        visited[start] = True

        while True:
            i = (i % row_len) * col_cnt + i // row_len
            if i == start:
                break
            visited[i] = True
            cycle.append(i)

        for k in reversed(range(1, len(cycle))):
            arr[cycle[k]], arr[cycle[k-1]] = arr[cycle[k-1]], arr[cycle[k]]
            
import random

def is_prime(n, k=5):  # number of tests
    if n <= 1 or (n % 2 == 0 and n > 2):
        return False
    # write n as d*2^r + 1
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    # Witness loop
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

# Function to calculate (a ** b) % p
def pow_mod(a, b, p):
    result = 1
    base = a % p
    while b > 0:
        if b % 2 == 1:
            result = (result * base) % p
        b //= 2
        base = (base * base) % p
    return result

# Function to find the prime factors of n
def prime_factors(n):
    factors = []
    # Count the number of 2's that divide n
    while n % 2 == 0:
        factors.append(2)
        n = n // 2
    # n must be odd at this point, skip one element (Note i = i +2)
    for i in range(3, int(n**0.5) + 1, 2):
        # While i divides n, append i and divide n
        while n % i == 0:
            factors.append(i)
            n = n // i
    # If n is a prime greater than 2
    if n > 2:
        factors.append(n)
    return list(set(factors))  # Remove duplicates

# Function to find a generator for a prime p
def find_generator(p):
    factors = prime_factors(p - 1)
    for candidate in range(2, p):
        if all(pow_mod(candidate, (p - 1) // factor, p) != 1 for factor in factors):
            return candidate
    return None  # This should never happen for a prime p


# Test the function
if __name__ == "__main__":
    #####################################
    arr = [1, 2, 3, 4, 5, 6]
    n1, n2 = 2, 3  # Dimensions before transpose
    in_place_transpose(arr, n1, n2)
    print("After in-place transpose:", arr)  # Should print [1, 4, 2, 5, 3, 6]
    assert arr == [1, 4, 2, 5, 3, 6], f"Expected [1, 4, 2, 5, 3, 6], got {arr}"
    print("Transpose test passed!")
    
    arr = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    expected = [1, 4, 7, 2, 5, 8, 3, 6, 9]
    n1, n2 = 3, 3  # Dimensions before transpose
    in_place_transpose(arr, n1, n2)
    print("After in-place transpose:", arr)  # Should print [1, 4, 2, 5, 3, 6]
    assert arr == expected, f"Expected {expected}, got {arr}"
    print("Transpose test passed!")
    
    # Start at 2^20 + 1
    
    rr = 2 ** 20
    # rr = 2 ** 9
    # n = rr + 1
    n = rr + 1

    # Search for a prime greater than 2^20
    while (not is_prime(n)) or (n - 1) % rr != 0:
        n += 1

    print(f"A prime greater than {rr} is {n}.")
    
    # Find a generator for the prime 521
    p = n
    g = find_generator(p)
    assert (p - 1) % 4 == 0
    print("A generator for p =", p, "is g =", g)
    
