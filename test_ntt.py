# -*- coding: utf-8 -*-
from random import randint
from ntt import NTT
import unittest

TEST_CASES = [
    (3329, 8), (12289, 16),  # Add more cases as needed
]


class TestNTT(unittest.TestCase):

    def test_ntt_intt(self, rep=100):
        """Test if ntt and intt are indeed inverses of each other."""
        for (q, n) in TEST_CASES:
            with self.subTest(q=q, n=n):
                T = NTT(q)
                for i in range(rep):
                    f = [randint(0, T.q-1) for j in range(n)]
                    assert T.intt(T.ntt(f)) == f

    def test_ntt_linearity(self, iterations=100):
        """Test the linearity of NTT."""
        for (q, n) in TEST_CASES:
            with self.subTest(q=q, n=n):
                T = NTT(q)
                for i in range(iterations):
                    f = [randint(0, T.q - 1) for j in range(n)]
                    g = [randint(0, T.q - 1) for j in range(n)]
                    λ = randint(0, T.q-1)
                    μ = randint(0, T.q-1)
                    λ_f_plus_μ_g = [(λ*x+μ*y) % T.q for (x, y) in zip(f, g)]
                    f_ntt = T.ntt(f)
                    g_ntt = T.ntt(g)
                    assert T.ntt(λ_f_plus_μ_g) == [(λ*x+μ*y) %
                                                   T.q for (x, y) in zip(f_ntt, g_ntt)]
