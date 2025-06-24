import unittest
from time import time
from polyntt.field.extension_field import ExtensionField
from polyntt.field.prime_field import PrimeField
from polyntt.field.m31 import M31ExtensionField, M31Field
from polyntt.params import PARAMS
from polyntt.polynomial.polynomial_ntt import PolynomialNTT, PolynomialRingNTT
from polyntt.polynomial.polynomial import PolynomialRing


class TestPolynomialNTT(unittest.TestCase):
    def setUp(self):
        self.F = M31ExtensionField()
        self.ring = PolynomialRingNTT(self.F, n=4)
        # self.p1 = PolynomialNTT(self.ring, [self.F(1), self.F(2), self.F(3)])
        # self.p2 = PolynomialNTT(self.ring, [self.F(4), self.F(5)])

    def test_mul_with_ntt_fp(self):
        # a test with BabyBear
        F = PrimeField(p=2**31 - 2**24 + 1)
        n = 256
        Fx = PolynomialRing(F, n)
        Fy = PolynomialRingNTT(F, n)
        for i in range(10):
            P = Fx.random()
            Q = Fx.random()
            P_mul_Q = P * Q
            P = Fy(P.coeffs)
            Q = Fy(Q.coeffs)
            P_mul_Q_2 = P*Q
            self.assertEqual(P_mul_Q, P_mul_Q_2)

    def test_mul_with_ntt_fp2(self):
        F = M31ExtensionField()
        for k in range(1, 9):
            n = 1 << k
            Fx = PolynomialRing(F, n)
            P = Fx.random()
            Q = Fx.random()
            P_mul_Q = P*Q
            Fy = PolynomialRingNTT(F, n)
            P_ntt = Fy(P.coeffs)
            Q_ntt = Fy(Q.coeffs)
            PQ_ntt = P_ntt * Q_ntt
            self.assertEqual(P_mul_Q, PQ_ntt)

    def test_compare_using_extension(self):

        nreps = 100

        # Polynomials of degree 256 over F_BB
        F = PrimeField(p=2013265921)
        Fx = PolynomialRingNTT(F, 1 << 2)
        P = Fx.random()
        Q = Fx.random()
        t = time()
        for i in range(nreps):
            R = P*Q
        print(
            "\nMult with Fp[x]/(x²⁵⁶+1): {:.0f}μs".format(10**6 * (time() - t)/nreps))

        # Polynomials of degree 128 over F_BB²
        F2 = ExtensionField(F, 2, 11)
        F2x = PolynomialRingNTT(F2, 1 << 1)
        P2 = F2x.random()
        Q2 = F2x.random()
        t = time()
        for i in range(nreps):
            R2 = P2*Q2
        print(
            "Mult with Fp²[x]/(x¹²⁸+1): {:.0f}μs".format(10**6*(time() - t)/nreps))

    def test_compare_ntt(self):

        nreps = 3

        # Polynomials of degree 256
        k = 8

        # Over BabyBear, degree 256
        F1 = PrimeField(p=2**31-2**24+1)
        F1x = PolynomialRingNTT(F1, 1 << k)
        P = F1x.random()
        t = time()
        for i in range(nreps):
            P_ntt = F1x.ntt(P)
        print(
            "\nNTT with F_BB[x]/(x²⁵⁶+1): {:.0f}μs".format(10**6 * (time() - t)/nreps))

        # Over M31, degree 128
        F2 = M31ExtensionField()
        F2x = PolynomialRingNTT(F2, 1 << (k-1))
        P2 = F2x.random()
        t = time()
        for i in range(nreps):
            P2_ntt = F2x.ntt(P2)
        print(
            "NTT with F_M31²[x]/(x¹²⁸+1): {:.0f}μs".format(10**6*(time() - t)/nreps))

    def test_ntt_intt(self, iterations=100):
        """Test if ntt and intt are indeed inverses of each other."""
        for (q, k) in PARAMS:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                F = M31ExtensionField() if q == 2**31-1 else PrimeField(q)
                Fx = PolynomialRingNTT(F, n)
                for i in range(iterations):
                    f = Fx([F.random() for j in range(n)])
                    self.assertEqual(Fx.intt(Fx.ntt(f)), f)

    def test_ntt_linearity(self, iterations=100):
        """Test the linearity of NTT."""
        for (q, k) in PARAMS:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                F = M31ExtensionField() if q == 2**31-1 else PrimeField(q)
                Fx = PolynomialRingNTT(F, n)
                for i in range(iterations):
                    f = Fx([F.random() for j in range(n)])
                    g = Fx([F.random() for j in range(n)])
                    λ = F.random()
                    μ = F.random()
                    λ_f_plus_μ_g = Fx([(λ*x+μ*y) for (x, y) in zip(f, g)])
                    f_ntt = Fx.ntt(f)
                    g_ntt = Fx.ntt(g)
                    self.assertEqual(
                        Fx.ntt(λ_f_plus_μ_g),
                        Fx([(λ*x+μ*y)
                           for (x, y) in zip(f_ntt, g_ntt)], ntt=True)
                    )

    # NTT_WITHOUT_MOD IS NOT IMPLEMENTED HERE
    # def test_ntt_without_mod(self, iterations=1):
    #     """Test for a more efficient intt."""
    #     # only for FALCON for now
    #     q, k = PARAMS[0]
    #     n = 1 << (k-1)
    #     F = PrimeField(q)
    #     Fx = PolynomialRingNTT(F, n)
    #     for i in range(iterations):
    #         f = [F.random() for j in range(n)]
    #         self.assertEqual(
    #             Fx.ntt(f),
    #             Fx.ntt_without_mod(f)
    #         )

    # INTT_WITHOUT_MOD IS NOT IMPLEMENTED HERE
    # def test_intt_without_mod(self, iterations=100):
    #     """Test for a more efficient intt."""
    #     q, k = PARAMS[0]
    #     n = 1 << (k-1)
    #     F = PrimeField(q)
    #     Fx = PolynomialRingNTT(F, n)
    #     for i in range(iterations):
    #         f = [F.random() for j in range(n)]
    #         self.assertEqual(
    #             Fx.intt(Fx.ntt(f)),
    #             Fx.intt_without_mod(Fx.ntt(f))
    #         )

    def test_parent(self):
        """Test the linearity of NTT."""
        for (q, k) in PARAMS:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                F = M31ExtensionField() if q == 2**31-1 else PrimeField(q)
                Fx = PolynomialRingNTT(F, n)
                for i in range(100):
                    f = Fx([F.random() for j in range(n)])
                    g = Fx([F.random() for j in range(n)])
                    f_plus_g = f+g
                    self.assertEqual(type(f), type(f_plus_g))

    def test_mul_with_integer(self):
        F = PrimeField(2**31-2**24+1)
        R = PolynomialRingNTT(F, 256)
        for i in range(100):
            P = R.random()
            Q = P * R([12] + [0 for i in range(255)])
            self.assertEqual(Q, P*12)
            self.assertEqual(Q, 12*P)


if __name__ == '__main__':
    unittest.main()
