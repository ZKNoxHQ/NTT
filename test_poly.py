# -*- coding: utf-8 -*-
from random import randint
from poly import Poly
import unittest


class TestPoly(unittest.TestCase):

    def test_add_sub(self, rep=100):
        """Test if ntt and intt are indeed inverses of each other."""
        q = 3329
        n = 8
        for i in range(rep):
            f = Poly([randint(0, q-1) for _ in range(n)], q, 0xdeadbeef)
            g = Poly([randint(0, q-1) for _ in range(n)], q, 0xdeadbeef)
            f_plus_g = f+g
            assert f_plus_g - g == f

    def test_mod_q(self, rep=100):
        """ Test if the reduction mod q works."""
        q = 3329
        n = 8
        zero = Poly([0 for _ in range(n)], q, 0xdeadbeef)
        for i in range(rep):
            f = Poly([q*randint(0, q-1) for _ in range(n)], q, 0xdeadbeef)
            assert f == zero

    def test_mul(self, rep=100):
        """Compare FFT multiplication with schoolbook multiplication."""
        q = 3329
        n = 8
        for i in range(rep):
            f = Poly([randint(0, q-1) for _ in range(n)], q, 0xdeadbeef)
            g = Poly([randint(0, q-1) for _ in range(n)], q, 0xdeadbeef)
            f_mul_g = f*g
            assert f_mul_g == f.mul_schoolbook(g)

    def test_div(n, rep=100):
        """Test the diviison."""
        q = 3329
        n = 8
        for i in range(rep):
            # random f
            f = Poly([randint(0, q-1) for _ in range(n)], q, 0xdeadbeef)
            # invertible random g
            g = Poly(f.NTT.intt([randint(1, q-1)
                     for _ in range(n)]), q, 0xdeadbeef)
            h = f.div(g)
            assert h * g == f
