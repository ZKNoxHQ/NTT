from polyntt.field.m31 import M31ExtensionField
from polyntt.field.prime_field import PrimeField
from polyntt.params import PARAMS
from polyntt.utils import bit_reverse_order  # , sqrt_mod

#
# Generate constants for the iterative case
#

f = open("polyntt/polynomial/ntt_constants_iterative.py", "w")
f.write("# File generated with `python polyntt/generate_ntt_constants.py`.\n")
f.write(
    "# Precomputations for NTT.\n\n"
)

ψ_table = dict()
ψ_inv_table = dict()
ψ_rev = dict()
ψ_inv_rev = dict()
n_inv = dict()

for (q, two_adicity) in PARAMS:
    # list of roots of cyclotomic polynomials
    if q == 12*1024 + 1:
        # Falcon
        # ψ is a root of the 2¹¹-th cyclotomic polynomial
        ψ = 7
    elif q == 3329:
        # Kyber
        # ψ is a root of the 2⁷-th cyclotomic polynomial
        ψ = 3296
    elif q == 8380417 and two_adicity == 9:
        # Dilithium
        # ψ is a root of the 2⁹-th cyclotomic polynomial
        ψ = 1753
    elif q == 2013265921:
        # BabyBear
        # ψ is a root of the 2⁹-th cyclotomic polynomial
        # (larger 2-adicity can be considered)
        ψ = 16303300
    elif q == 2**31-2**24+1:
        # KoalaBear
        # ψ is a root of the 2⁹-th cyclotomic polynomial
        # (larger 2-adicity can be considered)
        ψ = 1506173
    elif q == 2**31 - 1:
        # Mersenne 31
        # ψ is a root of the 2⁹-th cyclotomic polynomial
        # (larger 2-adicity can be considered)
        # defined over Fp² with i = sqrt(-i)
        ψ = [13610297, 1064696601]
    else:
        print("NOT DEFINED YET")

    if q != 2**31-1:
        F = PrimeField(q)
    else:
        F = M31ExtensionField()

    n = 1 << (two_adicity-1)
    ψ = F(ψ)
    toto = ψ
    for i in range(two_adicity-1):
        toto = toto * toto
    assert toto != F(1)
    assert toto * toto == F(1)

    ψ_inv = F(ψ).inverse()
    assert ψ*ψ_inv == F(1)

    # Precompute powers of ψ to speedup main NTT process.
    ψ_table[q] = [F(1)] * n
    ψ_inv_table[q] = [F(1)] * n
    for i in range(1, n):
        ψ_table[q][i] = ψ_table[q][i-1] * ψ
        ψ_inv_table[q][i] = ψ_inv_table[q][i-1] * ψ_inv

    # Change the lists into bit-reverse order.
    ψ_rev[q] = bit_reverse_order(ψ_table[q])
    ψ_inv_rev[q] = bit_reverse_order(ψ_inv_table[q])


# writing ψ
f.write("# Dictionary containing the powers ψ, a 2^n-th root of unity.\n")
f.write("ψ = {\n")
for (q, two_adicity) in PARAMS:
    f.write("\t# ψ = {}, ψ has multiplicative order {}.\n".format(
        ψ_table[q][1], 1 << two_adicity))
    f.write("\t{} : {},\n".format(q, ψ_table[q]))
f.write("}\n\n")

# writing ψ_inv
f.write("# Dictionary containing the powers of ψ_inv.\n")
f.write("ψ_inv = {\n")
for (q, two_adicity) in PARAMS:
    f.write("\t # ψ_inv = {}, ψ*ψ_inv = 1.\n".format(ψ_inv_table[q][1]))
    f.write("\t{} : {},\n".format(q, ψ_inv_table[q]))
f.write("}\n\n")

# writing ψ_rev
f.write(
    "# The table ψ, but in bit-reversed order, i.e. the i-th element corresponds to ψ^{BitReversed(i)}.\n")
f.write("ψ_rev = {\n")
for (q, two_adicity) in PARAMS:
    f.write("\t{} : {},\n".format(q, ψ_rev[q]))
f.write("}\n\n")

# writing ψ_rev_inv
f.write(
    "# The table ψ_inv, but in bit-reversed order, i.e. the i-th element corresponds to ψ^{BitReversed(-i)}.\n")
f.write("ψ_inv_rev = {\n")
for (q, two_adicity) in PARAMS:
    f.write("\t{} : {},\n".format(q, ψ_inv_rev[q]))
f.write("}\n\n")

# writing n_inv
f.write("# The inverses of powers of 2 mod q\n")
f.write("n_inv = {\n")
for (q, two_adicity) in PARAMS:
    f.write("\t{}: {{\n".format(q))
    # n_inv[{}] = {{\n".format(q))
    for j in range(two_adicity+1):
        if q != 2**31-1:
            f.write("\t\t{}: {},\n".format(1 << j, pow(1 << j, -1, q)))
        else:
            f.write("\t\t{}: {},\n".format(1 << j, [pow(1 << j, -1, q), 0]))
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

for (q, two_adicity) in PARAMS:
    file.write("\t{}: {{\n".format(q))
    if q != 2**31-1:
        phi_roots = [sqrt_mod(-1, q), q-sqrt_mod(-1, q)]
    else:
        phi_roots = [[0, 1], [0, p-1]]

    for k in range(1, two_adicity+(q == 2**31-1)):
        file.write("\t\t{} : {},\n".format(1 << k, phi_roots))
        if q != 2**31-1:
            phi_roots = sum([[sqrt_mod(elt, q), q - sqrt_mod(elt, q)]
                            for elt in phi_roots], [])
        else:
            phi_roots = sum([[sqrt_m31_2(elt), opp2(sqrt_m31_2(elt))]
                             for elt in phi_roots], [])

    file.write("\t},\n")
file.write("}\n")
file.close()
