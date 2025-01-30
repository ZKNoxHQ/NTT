def bit_reverse_order(a):
    '''
    Reorders the given array in reverse-bit order.
    [1, 2, 3, 4, 5, 6, 8] -> [1, 5, 3, 7, 2, 6, 4, 8]
    '''
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
f.write("ψ_inv_rev = dict()\n")
f.write("n_inv = dict()\n")


for q in [12 * 1024 + 1, 3329]:  # Falcon, Kyber
    # list of roots of cyclotomic polynomials
    if q == 12*1024 + 1:
        ψ = 1826  # a root of the 2¹¹-th cyclotomic polynomial
        k = 10
    elif q == 3329:
        ψ = 17  # a root of the 2⁸-th cyclotomic polynomial
        k = 7
    else:
        print("NOT DEFINED YET")
    n = 1 << k
    assert pow(ψ, n, q) == q-1
    assert pow(ψ, 2*n, q) == 1

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
    # f.write("ψ_{} = {}\n".format(q, ψ_table))
    f.write("ψ_inv_rev[{}] = {}\n".format(q, ψ_inv_rev))
    # f.write("ψ_{}_inv = {}\n".format(q, ψ_inv_table))

    f.write("# inverse of powers of 2 mod q\n")
    f.write("n_inv[{}] = {{\n".format(q))
    for k in range(1, 11):
        f.write("\t {}: {},\n".format(1 << k, pow(1 << k, -1, q)))
    f.write("}\n")
f.close()
