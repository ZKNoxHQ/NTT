from test_cases import TEST_CASES


def bit_reverse_order(a):
    '''Reorders the given array in reverse-bit order.'''
    num_bits = len(bin(len(a) - 1)) - 2
    result = [0] * len(a)
    for i in range(len(a)):
        rev_index = int(bin(i)[2:].zfill(num_bits)[::-1], 2)
        result[rev_index] = a[i]
    return result


f = open("ntt_constants.py", "w")
f.write("# File generated with `python generate_constants.py`.\n")
f.write(
    "# Precomputations for NTT (reverse table of ψ and its inverse for different n.\n"
)
f.write("ψ_rev = dict()\n")
f.write("ψ = dict()\n")
f.write("ψ_inv_rev = dict()\n")
f.write("ψ_inv = dict()\n")
f.write("n_inv = dict()\n")


for (q, two_adicity) in TEST_CASES:
    # list of roots of cyclotomic polynomials
    if q == 12*1024 + 1:
        ψ = 1826  # a root of the 2¹¹-th cyclotomic polynomial
    elif q == 3329:
        ψ = 17  # a root of the 2⁸-th cyclotomic polynomial
    elif q == 2013265921:
        ψ = 137  # a root of the 2²⁷-th cyclotomic polynomial
    elif q == 18446744069414584321:
        ψ = 7277203076849721926  # a root of the 2³²-th cyclotomic polynomial
    else:
        print("NOT DEFINED YET")
    n = 1 << (k-1)
    assert pow(ψ, 2*n, q) == 1 and pow(ψ, n, q) != 1

    ψ_inv = pow(ψ, -1, q)
    assert (ψ*ψ_inv) % q == 1

    # Precompute powers of ψ to speedup main NTT process.
    ψ_table = [1] * n
    ψ_inv_table = [1] * n
    for i in range(1, n):
        ψ_table[i] = ((ψ_table[i-1] * ψ) % q)
        ψ_inv_table[i] = ((ψ_inv_table[i-1] * ψ_inv) % q)

    # Change the lists into bit-reverse order.
    ψ_rev = bit_reverse_order(ψ_table)
    ψ_inv_rev = bit_reverse_order(ψ_inv_table)

    f.write("ψ_rev[{}] = {}\n".format(q, ψ_rev))
    f.write("ψ[{}] = {}\n".format(q, ψ_table))
    f.write("ψ_inv_rev[{}] = {}\n".format(q, ψ_inv_rev))
    f.write("ψ_inv[{}] = {}\n".format(q, ψ_inv_table))

    f.write("# inverse of powers of 2 mod q\n")
    f.write("n_inv[{}] = {{\n".format(q))
    for k in range(1, 11):
        f.write("\t {}: {},\n".format(1 << k, pow(1 << k, -1, q)))
    f.write("}\n")
f.close()
