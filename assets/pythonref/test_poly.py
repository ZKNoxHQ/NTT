# -*- coding: utf-8 -*-
from random import randint
from poly import Poly
import unittest
from test_cases import TEST_CASES


class TestPoly(unittest.TestCase):

    def shortDescription(self):
        return None  # This prevents unittest from printing docstrings

    # def add_sub(self, q, n):
    #     with self.subTest(msg="Test of add_sub with q={} and n = {}".format(q, n)):
    #         f = Poly([randint(0, q-1) for _ in range(n)], q)
    #         g = Poly([randint(0, q-1) for _ in range(n)], q)
    #         f_plus_g = f+g
    #         self.assertEqual(f_plus_g - g, f)

    # def test_all(self, iterations=100):
    #     for (q, k) in TEST_CASES:
    #         n = 1 << (k-1)
    #         self.add_sub(q, n)

    def test_add_sub(self, iterations=100):
        """Test if ntt and intt are indeed inverses of each other."""
        for (q, k) in TEST_CASES:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                for i in range(iterations):
                    f = Poly([randint(0, q-1) for _ in range(n)], q)
                    g = Poly([randint(0, q-1) for _ in range(n)], q)
                    f_plus_g = f+g
                    self.assertEqual(f_plus_g - g, f)

    def test_mod_q(self, iterations=100):
        """ Test if the reduction mod q works."""
        for (q, k) in TEST_CASES:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                zero = Poly([0 for _ in range(n)], q)
                for i in range(iterations):
                    f = Poly([q*randint(0, q-1) for _ in range(n)], q)
                    self.assertEqual(f, zero)

    def test_mul(self, iterations=100):
        """Compare FFT multiplication with schoolbook multiplication."""
        for (q, k) in TEST_CASES:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                if k > 8:
                    self.skipTest("Schoolbook is too slow for k>8.")
                else:
                    for i in range(iterations):
                        f = Poly([randint(0, q-1) for _ in range(n)], q)
                        g = Poly([randint(0, q-1) for _ in range(n)], q)
                        f_mul_g = f*g
                        self.assertEqual(f_mul_g, f.mul_schoolbook(g))

    def test_div(self, iterations=100):
        """Test the diviison."""
        for (q, k) in TEST_CASES:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                for i in range(iterations):
                    # random f
                    f = Poly([randint(0, q-1) for _ in range(n)], q)
                    # invertible random g
                    g = Poly(f.NTT.intt([randint(1, q-1)
                                         for _ in range(n)]), q)
                    h = f.div(g)
                    self.assertEqual(h * g, f)

    def test_mul_pwc(self, iterations=100):
        """Test the multiplication modulo x^n+1."""
        for (q, k) in TEST_CASES:
            n = 1 << (k-1)
            with self.subTest(q=q, k=k):
                if k > 8:
                    self.skipTest("Schoolbook is too slow for k>8.")
                else:
                    for i in range(1, iterations):
                        # random f,g
                        f = Poly([randint(0, q-1) for _ in range(n)], q)
                        g = Poly([randint(0, q-1) for _ in range(n)], q)
                        self.assertEqual(f.mul_pwc(g), f.mul_schoolbook_pwc(g))
