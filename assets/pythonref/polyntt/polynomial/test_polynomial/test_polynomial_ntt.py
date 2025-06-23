import unittest
from time import time
from polyntt.field.extension_field import ExtensionField
from polyntt.field.prime_field import PrimeField
from polyntt.field.m31 import M31ExtensionField
from polyntt.params import PARAMS
from polyntt.polynomial.polynomial_ntt import PolynomialNTT, PolynomialRingNTT
from polyntt.polynomial.polynomial import PolynomialRing


class TestPolynomialNTT(unittest.TestCase):
    def setUp(self):
        self.F = M31ExtensionField()
        self.ring = PolynomialRingNTT(self.F, n=4)
        # self.p1 = PolynomialNTT(self.ring, [self.F(1), self.F(2), self.F(3)])
        # self.p2 = PolynomialNTT(self.ring, [self.F(4), self.F(5)])

    # def test_mul_with_ntt_fp(self):
    #     # a test with BabyBear
    #     F = PrimeField(p=2**31 - 2**24 + 1)
    #     Fx = PolynomialRing(F, n=4)
    #     P = Fx.random()
    #     Q = Fx.random()
    #     P_mul_Q = P.schoolbook_mul(Q)
    #     Fy = PolynomialRingNTT(F, n=4)
    #     P_ntt = Fy.ntt(P.coeffs)
    #     Q_ntt = Fy.ntt(Q.coeffs)
    #     PQ_ntt = [a*b for (a, b) in zip(P_ntt, Q_ntt)]
    #     PQ = Fy.intt(PQ_ntt)
    #     self.assertEqual(P_mul_Q, Fx(PQ))

    def test_mul_with_ntt_fp(self):
        # a test with BabyBear
        F = PrimeField(p=12289)
        # F = PrimeField(p=2**31 - 2**24 + 1)
        for k in range(1, 9):
            n = 1 << k
            Fx = PolynomialRingNTT(F, n)
            P = Fx.random()
            Q = Fx.random()
            P_mul_Q = P*Q
            P_ntt = Fx.ntt(P)
            Q_ntt = Fx.ntt(Q)
            PQ_ntt = P_ntt * Q_ntt
            PQ = Fx.intt(PQ_ntt)
            self.assertEqual(P_mul_Q, PQ)

    def test_mul_with_ntt_fp2(self):
        F = M31ExtensionField()
        for k in range(1, 9):
            n = 1 << k
            Fx = PolynomialRing(F, n)
            P = Fx.random()
            Q = Fx.random()
            P_mul_Q = P*Q
            Fy = PolynomialRingNTT(F, n)
            P_ntt = Fy.ntt(P)
            Q_ntt = Fy.ntt(Q)
            PQ_ntt = P_ntt * Q_ntt
            PQ = Fy.intt(PQ_ntt)
            self.assertEqual(P_mul_Q, PQ)

    def test_compare_using_extension(self):

        nreps = 100

        # Polynomials of degree 256 over F_BB
        F = PrimeField(p=2013265921)
        Fx = PolynomialRingNTT(F, 1 << 8)
        P = Fx.random()
        Q = Fx.random()
        t = time()
        for i in range(nreps):
            R = P*Q
        print(
            "\nMult with Fp[x]/(x²⁵⁶+1): {:.0f}μs".format(10**6 * (time() - t)/nreps))

        # Polynomials of degree 128 over F_BB²
        F2 = ExtensionField(F, 2, 11)
        # F2x = PolynomialRingNTT(F2, 1 << 7)
        F2x = PolynomialRingNTT(F2, 1 << 1)
        print('comp P2')
        P2 = F2x.random()
        print('p2')
        print(P2)
        Q2 = F2x.random()
        t = time()
        for i in range(nreps):
            R2 = P2*Q2
        print(
            "Mult with Fp²[x]/(x¹²⁸+1): {:.0f}μs".format(10**6*(time() - t)/nreps))

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

    def test_edge(self):
        P = PolynomialNTT(self.ring, [self.F(1)])


if __name__ == '__main__':
    unittest.main()
