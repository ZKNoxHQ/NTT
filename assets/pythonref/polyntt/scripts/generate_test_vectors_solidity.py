import hashlib
from polyntt.field.m31 import M31ExtensionField
from polyntt.field.prime_field import PrimeField
from polyntt.polynomial.polynomial_ntt import PolynomialRingNTT
from polyntt.params import PARAMS
from polyntt.scripts.generate_test_vectors import encode, deterministic_poly

for (q, two_adicity) in PARAMS:
    F = M31ExtensionField() if q == 2**31-1 else PrimeField(q)
    for n in [1 << (two_adicity-2), 1 << (two_adicity-1)]:  # for two sizes of polynomials
        Fx = PolynomialRingNTT(F, n)
        file = open("../test_vectors/q{}_n{}.sol".format(q, n), "w")

        file.write(
            "// File generated using ../pythonref/scripts/generate_test_vectors_solidity.py\n\n")
        f = Fx(deterministic_poly(F, n, seed="seed_f"))
        g = Fx(deterministic_poly(F, n, seed="seed_g"))

        f_ntt = Fx.ntt(f)
        g_ntt = Fx.ntt(g)
        f_ntt_intt = Fx.intt(f_ntt)
        g_ntt_intt = Fx.intt(g_ntt)
        assert f_ntt_intt == f
        assert g_ntt_intt == g

        f_ntt_mul_g_ntt = f_ntt*g_ntt
        f_mul_g = Fx.intt(f_ntt_mul_g_ntt)
        assert f_mul_g == f*g

        f_ntt_add_g_ntt = f_ntt + g_ntt
        f_ntt_sub_g_ntt = f_ntt - g_ntt

        # 1. ntt
        #   Input: f,g
        #   Output: ntt(f), ntt(g).
        file.write("// ntt of f and g;\n")
        file.write("uint256[] f = {};\n".format(f))
        file.write("uint256[] g = {};\n".format(g))
        file.write("uint256[] f_ntt = {};\n".format(f_ntt))
        file.write("uint256[] g_ntt = {};\n".format(g_ntt))
        file.write("\n")

        # 2. intt
        #   Input: f_ntt,g_ntt
        #   Output: f, g
        file.write("// intt of f_ntt and g_ntt;\n")
        file.write("uint256[] f_ntt = {};\n".format(f_ntt))
        file.write("uint256[] g_ntt = {};\n".format(g_ntt))
        file.write("uint256[] f = {};\n".format(f))
        file.write("uint256[] g = {};\n".format(g))
        file.write("\n")

        # 3. vec_mul
        #   Input: f_ntt, g_ntt,
        #   Output: f_ntt*g_ntt.
        file.write("// vec_mul of f_ntt and g_ntt;\n")
        file.write("uint256[] f_ntt = {};\n".format(f_ntt))
        file.write("uint256[] g_ntt = {};\n".format(g_ntt))
        file.write("uint256[] f_ntt_mul_g_ntt = {};\n".format(f_ntt_mul_g_ntt))
        file.write("\n")

        # 4. pol_mul
        #   Input: f,g,
        #   Output: f*g.
        file.write("// pol_mul of f and g;\n")
        file.write("uint256[] f = {};\n".format(f))
        file.write("uint256[] g = {};\n".format(g))
        file.write("uint256[] f_mul_g = {};\n".format(f_mul_g))
        file.write("\n")

        # 5. vec_add
        #   Input: f_ntt, g_ntt,
        #   Output: f_ntt + g_ntt.
        file.write("// vec_add of f_ntt and g_ntt;\n")
        file.write("uint256[] f_ntt = {};\n".format(f_ntt))
        file.write("uint256[] g_ntt = {};\n".format(g_ntt))
        file.write("uint256[] f_ntt_add_g_ntt = {};\n".format(f_ntt_add_g_ntt))
        file.write("\n")

        # 6. vec_sub
        #   Input: f_ntt, g_ntt,
        #   Output: f_ntt - g_ntt.
        file.write("// vec_sub of f_ntt and g_ntt;\n")
        file.write("uint256[] f_ntt = {};\n".format(f_ntt))
        file.write("uint256[] g_ntt = {};\n".format(g_ntt))
        file.write("uint256[] f_ntt_sub_g_ntt = {};\n".format(f_ntt_sub_g_ntt))
        file.write("\n")
    file.close()
