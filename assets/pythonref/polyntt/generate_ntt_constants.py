from polyntt.test_cases import TEST_CASES
from polyntt.utils import bit_reverse_order, sqrt_mod


#
# Generate constants for the iterative case
#

f = open("polyntt/ntt_constants_iterative.py", "w")
f.write("# File generated with `python polyntt/generate_ntt_constants.py`.\n")
f.write(
    "# Precomputations for NTT.\n\n"
)

ψ_table = dict()
ψ_inv_table = dict()
ψ_rev = dict()
ψ_inv_rev = dict()
n_inv = dict()

for (q, two_adicity) in TEST_CASES:
    # list of roots of cyclotomic polynomials
    if q == 12*1024 + 1:
        ψ = 1826  # a root of the 2¹¹-th cyclotomic polynomial
    elif q == 3329:
        ψ = 3296  # a root of the 2⁷-th cyclotomic polynomial
    elif q == 2013265921:
        # a root of the 2¹⁰ cyclotomic polynomial (but could be larger 2-adicity)
        ψ = 1538055801
        # ψ = 137  # a root of the 2²⁷-th cyclotomic polynomial
    elif q == 18446744069414584321:
        ψ = 7277203076849721926  # a root of the 2³²-th cyclotomic polynomial
    elif q == 8380417:  # Dilithium
        ψ = 2926054
    else:
        print("NOT DEFINED YET")
    n = 1 << (two_adicity-1)
    assert pow(ψ, 2*n, q) == 1 and pow(ψ, n, q) != 1

    ψ_inv = pow(ψ, -1, q)
    assert (ψ*ψ_inv) % q == 1

    # Precompute powers of ψ to speedup main NTT process.
    ψ_table[q] = [1] * n
    ψ_inv_table[q] = [1] * n
    for i in range(1, n):
        ψ_table[q][i] = ((ψ_table[q][i-1] * ψ) % q)
        ψ_inv_table[q][i] = ((ψ_inv_table[q][i-1] * ψ_inv) % q)

    # Change the lists into bit-reverse order.
    ψ_rev[q] = bit_reverse_order(ψ_table[q])
    ψ_inv_rev[q] = bit_reverse_order(ψ_inv_table[q])

# writing ψ
f.write("# Dictionary containing the powers ψ, a 2^n-th root of unity.\n")
f.write("ψ = {\n")
for (q, two_adicity) in TEST_CASES:
    f.write("\t# ψ = {}, ψ has multiplicative order {}.\n".format(
        ψ_table[q][1], 1 << two_adicity))
    f.write("\t{} : {},\n".format(q, ψ_table[q]))
f.write("}\n\n")

# writing ψ_inv
f.write("# Dictionary containing the powers of ψ_inv.\n")
f.write("ψ_inv = {\n")
for (q, two_adicity) in TEST_CASES:
    f.write("\t # ψ_inv = {}, ψ*ψ_inv = 1.\n".format(ψ_inv_table[q][1]))
    f.write("\t{} : {},\n".format(q, ψ_inv_table[q]))
f.write("}\n\n")

# writing ψ_rev
f.write(
    "# The table ψ, but in bit-reversed order, i.e. the i-th element corresponds to ψ^{BitReversed(i)}.\n")
f.write("ψ_rev = {\n")
for (q, two_adicity) in TEST_CASES:
    f.write("\t{} : {},\n".format(q, ψ_rev[q]))
f.write("}\n\n")

# writing ψ_rev_inv
f.write(
    "# The table ψ_inv, but in bit-reversed order, i.e. the i-th element corresponds to ψ^{BitReversed(-i)}.\n")
f.write("ψ_inv_rev = {\n")
for (q, two_adicity) in TEST_CASES:
    f.write("\t{} : {},\n".format(q, ψ_inv_rev[q]))
f.write("}\n\n")

# writing n_inv
f.write("# The inverses of powers of 2 mod q\n")
f.write("n_inv = {\n")
for (q, two_adicity) in TEST_CASES:
    f.write("\t{}: {{\n".format(q))
    # n_inv[{}] = {{\n".format(q))
    for j in range(1, two_adicity+1):
        f.write("\t\t{}: {},\n".format(1 << j, pow(1 << j, -1, q)))
    f.write("\t},\n")
f.write("}")

f.close()

#
# Generate constants for the recursive case
#

file = open("polyntt/ntt_constants_recursive.py", 'w')
file.write(
    "# Roots of the cyclotomic polynomials mod q for the recursive ntt implementation\n")
file.write("# File generated using `generate_contants_recursive.sage`\n")
file.write(
    "# roots_dict_mod[q][n] corresponds to the roots of x^{2n} + 1 mod q\n")
file.write("roots_dict_mod = {\n")

for (q, two_adicity) in TEST_CASES:
    file.write("\t{}: {{\n".format(q))
    phi_roots_Zq = [sqrt_mod(-1, q), q-sqrt_mod(-1, q)]
    for k in range(1, two_adicity):
        file.write("\t\t{} : {},\n".format(1 << k, phi_roots_Zq))
        phi_roots_Zq = sum([[sqrt_mod(elt, q), q - sqrt_mod(elt, q)]
                           for elt in phi_roots_Zq], [])
    file.write("\t},\n")
file.write("}\n")
file.close()
