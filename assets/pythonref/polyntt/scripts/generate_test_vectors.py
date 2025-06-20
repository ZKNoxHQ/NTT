import hashlib
from polyntt.field.m31 import M31ExtensionField
from polyntt.field.prime_field import PrimeField
# from polyntt.polynomial.ntt_iterative import NTTIterative
# from polyntt.poly import Poly
from polyntt.params import PARAMS
from polyntt.polynomial.polynomial import PolynomialRing
from polyntt.polynomial.polynomial_ntt import PolynomialRingNTT


def write_test(f, input, name, expected, gas, final=False):
    f.write("\t{\n")
    f.write("\t\t\"Input\": \"{}\",\n".format(input))
    f.write("\t\t\"Name\": \"{}\",\n".format(name))
    f.write("\t\t\"Expected\": \"{}\",\n".format(expected))
    f.write("\t\t\"Gas\": {}\n".format(gas))
    if final:
        f.write("\t}\n")
    else:
        f.write("\t},\n")


def encode(poly):
    # By default, q is a 32-bit integer.
    # NB: this does not apply to q_plonky2.
    q = poly.parent.F.p
    assert all([[x.coeffs < q for x in c] for c in poly.coeffs])
    size_q = (q.bit_length()+7)//8
    byte_string = b''.join(
        b''.join(num.coeffs.to_bytes(size_q, 'big') for num in c.coeffs)
        for c in poly.coeffs
    )
    return byte_string.hex()


def decode(hex_poly, q, field='Prime'):
    # By default, q is a 32-bit integer.
    # NB: this does not apply to q_plonky2.
    size_q = (q.bit_length()+7)//8
    bytes_poly = bytes.fromhex(hex_poly)
    if field == 'Prime':
        F = PrimeField(q)
        return [F(int.from_bytes(bytes_poly[i:i+size_q], 'big')) for i in range(0, len(bytes_poly), size_q)]
    else:  # M31 quadratic extension
        F = M31ExtensionField()
        return [F([
                int.from_bytes(bytes_poly[i:i+size_q], 'big'),
                int.from_bytes(bytes_poly[i+size_q:i+2*size_q], 'big')
                ]) for i in range(0, len(bytes_poly), 2*size_q)]


def deterministic_poly(F, n, seed="fixed_seed"):
    # This function is used for generating polynomials for the tests.
    # No randomness.
    q = F.p
    if F.degree == 1:
        return [F(int(hashlib.sha256(f"{seed}{i}".encode()).hexdigest(), 16) % q) for i in range(n)]
    else:  # F.degree == 2
        return [F([
            int(hashlib.sha256(f"{seed}{i}{0}".encode()).hexdigest(), 16) % q,
            int(hashlib.sha256(f"{seed}{i}{1}".encode()).hexdigest(), 16) % q
        ]) for i in range(n)]


F = M31ExtensionField()
Fx = PolynomialRing(F, 4)
f = Fx(deterministic_poly(F, 4))
assert decode(encode(f), F.p, field='Ext') == f.coeffs

for (q, two_adicity) in PARAMS:
    for n in [1 << (two_adicity-2), 1 << (two_adicity-1)]:  # for two sizes of polynomials
        file = open("../test_vectors/q{}_n{}.json".format(q, n), "w")
        F = M31ExtensionField() if q == 2**31-1 else PrimeField(q)
        Fx = PolynomialRing(F, n)
        FxNTT = PolynomialRingNTT(F, n)
        f = Fx(deterministic_poly(F, n, seed="seed_f"))
        g = Fx(deterministic_poly(F, n, seed="seed_g"))

        f_ntt = FxNTT.ntt(f)
        g_ntt = FxNTT.ntt(g)
        f_ntt_intt = FxNTT.intt(f_ntt)
        g_ntt_intt = FxNTT.intt(g_ntt)
        assert f_ntt_intt == f
        assert g_ntt_intt == g

        f_ntt_mul_g_ntt = f_ntt * g_ntt
        f_mul_g = FxNTT.intt(f_ntt_mul_g_ntt)
        assert f_mul_g == f*g

        f_ntt_add_g_ntt = f_ntt + g_ntt
        f_ntt_sub_g_ntt = f_ntt - g_ntt

        file.write("[\n")
        # 1. ntt
        #   Input: f,g
        #   Output: ntt(f), ntt(g).
        print(f.coeffs[0].coeffs)
        print(q)
        input = encode(f)+encode(g)
        name = "q{}_n{}_{{ntt_of_two_polynomials}}".format(q, n)
        expected = encode(f_ntt)+encode(g_ntt)
        gas = 600
        write_test(file, input, name, expected, gas)

        # 2. intt
        #   Input: f_ntt,g_ntt
        #   Output: f, g
        input = encode(f_ntt)+encode(g_ntt)
        name = "q{}_n{}_{{intt_of_two_polynomials}}".format(q, n)
        expected = encode(f_ntt_intt)+encode(g_ntt_intt)
        gas = 600
        write_test(file, input, name, expected, gas)

        # 3. vec_mul
        #   Input: f_ntt, g_ntt,
        #   Output: f_ntt*g_ntt.
        input = encode(f_ntt)+encode(g_ntt)
        name = "q{}_n{}_{{vec_mul}}".format(q, n)
        expected = encode(f_ntt_mul_g_ntt)
        gas = 300
        write_test(file, input, name, expected, gas)

        # 4. pol_mul
        #   Input: f,g,
        #   Output: f*g.
        input = encode(f)+encode(g)
        name = "q{}_n{}_{{pol_mul}}".format(q, n)
        expected = encode(f_mul_g)
        gas = 600 + 600 + 300 + 600
        write_test(file, input, name, expected, gas)

        # 5. vec_add
        #   Input: f_ntt, g_ntt,
        #   Output: f_ntt + g_ntt.
        input = encode(f_ntt)+encode(g_ntt)
        name = "q{}_n{}_{{vec_add}}".format(q, n)
        expected = encode(f_ntt_add_g_ntt)
        gas = 50  # TODO
        write_test(file, input, name, expected, gas)

        # 6. vec_sub
        #   Input: f_ntt, g_ntt,
        #   Output: f_ntt - g_ntt.
        input = encode(f_ntt)+encode(g_ntt)
        name = "q{}_n{}_{{vec_sub}}".format(q, n)
        expected = encode(f_ntt_sub_g_ntt)
        gas = 100  # TODO
        write_test(file, input, name, expected, gas, final=True)

        file.write("]\n")
    file.close()
