# -*- coding: utf-8 -*-
from random import randint
from polyntt.ntt_iterative import NTTIterative
from polyntt.poly import Poly
from polyntt.ntt_mert import NTTMert
import unittest
from polyntt.test_cases import TEST_CASES  # We test only Falcon parameters


class TestNTTMert(unittest.TestCase):
    def shortDescription(self):
        return None  # This prevents unittest from printing docstrings

    def test_ntt_intt(self, iterations=100):
        """Test if ntt and intt are indeed inverses of each other."""
        for (q, k) in [TEST_CASES[1]]:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                T = NTTMert(q)
                for i in range(iterations):
                    f = [randint(0, T.q-1) for j in range(n)]
                    self.assertEqual(T.intt(T.ntt(f)), f)

    def test_ntt_linearity(self, iterations=100):
        """Test the linearity of NTT."""
        for (q, k) in [TEST_CASES[1]]:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                T = NTTMert(q)
                for i in range(iterations):
                    f = [randint(0, T.q - 1) for j in range(n)]
                    g = [randint(0, T.q - 1) for j in range(n)]
                    λ = randint(0, T.q-1)
                    μ = randint(0, T.q-1)
                    λ_f_plus_μ_g = [(λ*x+μ*y) % T.q for (x, y) in zip(f, g)]
                    f_ntt = T.ntt(f)
                    g_ntt = T.ntt(g)
                    self.assertEqual(
                        T.ntt(λ_f_plus_μ_g),
                        [(λ*x+μ*y) % T.q for (x, y) in zip(f_ntt, g_ntt)]
                    )

    def test_mul_pwc(self, iterations=100):
        q = 3329
        n = 64
        Mert = NTTMert(q)
        for i in range(iterations):
            f = [randint(0, q-1) for i in range(n)]
            F = Mert.ntt(f)
            g = [randint(0, q-1) for i in range(n)]
            G = Mert.ntt(g)
            # check multiplication mod x^n-1
            F_mul_G = Mert.vec_mul(F, G)
            f_mul_g = Mert.intt(F_mul_G)
            self.assertEqual(f_mul_g, Mert.mul_schoolbook_x_n_minus_one(f, g))

    def test_mul_nwc(self, iterations=100):
        q = 3329
        n = 64
        Mert = NTTMert(q)
        for i in range(iterations):
            f = [randint(0, q-1) for i in range(n)]
            fp = Mert.prec(f)
            Fp = Mert.ntt(fp)
            g = [randint(0, q-1) for i in range(n)]
            gp = Mert.prec(g)
            Gp = Mert.ntt(gp)
            # check multiplication mod x^n+1
            Fp_mul_Gp = Mert.vec_mul(Fp, Gp)
            fp_mul_gp = Mert.intt(Fp_mul_Gp)
            f_mul_g = Mert.post(fp_mul_gp)
            self.assertEqual(f_mul_g, (Poly(f, q) * Poly(g, q)).coeffs)

    def test_ntt_iter_ntt_mert_prec_post_consistency(self, iterations=100):
        q = 3329
        n = 64
        Mert = NTTMert(q)
        T = NTTIterative(q)
        for i in range(iterations):
            f = [randint(0, q-1) for i in range(n)]
            F = Mert.ntt(f)
            g = [randint(0, q-1) for i in range(n)]
            self.assertEqual(Mert.ntt(Mert.prec(f)), T.ntt(f))
            self.assertEqual(Mert.post(Mert.intt(F)), T.intt(F))
            self.assertEqual(Mert.ntt(f), T.ntt(Mert.post(f)))
            self.assertEqual(Mert.intt(F), Mert.prec(T.intt(F)))

            F1 = T.ntt(Mert.post(f))
            G1 = T.ntt(Mert.post(g))
            F1_mul_G1 = T.vec_mul(F1, G1)
            f1_mul_g1 = Mert.prec(T.intt(F1_mul_G1))
            self.assertEqual(
                f1_mul_g1, Mert.mul_schoolbook_x_n_minus_one(f, g))
