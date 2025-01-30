# -*- coding: utf-8 -*-
from ntt_constants import n_inv, ψ_rev, ψ_inv_rev, ψ as Ψ, ψ_inv as Ψ_inv
from random import randint
from poly import Poly
from math import log
import unittest
from ntt import NTT

q = 3329


class TestMert(unittest.TestCase):
    def intReverse(self, a, n):
        b = ('{:0'+str(n)+'b}').format(a)
        return int(b[::-1], 2)

    def ntt_mert(self, A):
        # from https://github.com/acmert/ntt-based-polmul/blob/3476c2d9f22eac7696f2785eba036ffae1b4460e/baseline/ntt.py#L206
        N = len(A)
        B = [_ for _ in A]

        for s in range(int(log(N, 2)), 0, -1):
            m = 2**s
            for k in range(int(N/m)):
                # TW = pow(ω, self.intReverse(k, int(log(N, 2))-s)*int(m/2), q)
                TW = ψ_rev[q][k]
                for j in range(int(m/2)):
                    u = B[k*m+j]
                    t = (TW*B[k*m+j+int(m/2)]) % q

                    B[k*m+j] = (u+t) % q
                    B[k*m+j+int(m/2)] = (u-t) % q

        return B

    def intt_mert(self, A):
        # from
        # https://github.com/acmert/ntt-based-polmul/blob/3476c2d9f22eac7696f2785eba036ffae1b4460e/baseline/ntt.py#L353
        # and
        # https://github.com/acmert/ntt-based-polmul/blob/3476c2d9f22eac7696f2785eba036ffae1b4460e/baseline/ntt.py#L502
        N = len(A)
        B = [_ for _ in A]

        m = 1
        v = N
        d = 1

        while v > 1:
            for jf in range(m):
                j = jf
                jt = 0
                while j < (N-1):
                    # bit-reversing jt
                    # TW = pow(ω_inv, self.intReverse(jt, int(log(N >> 1, 2))), q)
                    TW = ψ_inv_rev[q][jt]

                    temp = B[j]

                    B[j] = (temp + B[j+d]) % q
                    B[j+d] = (temp - B[j+d])*TW % q

                    jt = jt+1
                    j = j + 2*m
            m = 2*m
            v = int(v/2)
            d = 2*d

        return [(x*n_inv[q][N]) % q for x in B]

    def mul_schoolbook_x_n_minus_one(self, A, B):
        """Multiplication of two polynomials using the schoolbook algorithm."""
        f = A
        g = B
        n = len(f)
        assert n == len(g)
        C = [0] * (2 * n)
        D = [0] * (n)
        for j, f_j in enumerate(f):
            for k, g_k in enumerate(g):
                C[j+k] = (C[j+k] + f_j * g_k) % q
        # reduction modulo x^n  + 1
        for i in range(n):
            D[i] = (C[i] + C[i+n]) % q
        return D

    def prec(self, A):
        n = len(A)
        # list of roots for the precomputations
        Ψ0 = Ψ[q][::len(Ψ[q])//n]
        return [(x * y) % q for (x, y) in zip(A, Ψ0)]

    def post(self, A):
        n = len(A)
        # list of roots for the precomputations
        Ψ0_inv = Ψ_inv[q][::len(Ψ_inv[q])//n]
        return [(x * y) % q for (x, y) in zip(A, Ψ0_inv)]

    def test_mert(self):

        f = [1, 2, 3, 4]
        # f = [randint(0, q-1) for i in range(4)]
        F = self.ntt_mert(f)
        back_f = self.intt_mert(F)
        for i in range(len(f)):
            assert back_f[i] == f[i]
        g = [5, 6, 7, 8]
        # g = [randint(0, q-1) for i in range(4)]
        G = self.ntt_mert(g)
        back_g = self.intt_mert(G)
        for i in range(len(g)):
            assert back_g[i] == g[i]

        # check multiplication mod x^n-1
        F_mul_G = [(x*y) % q for (x, y) in zip(F, G)]
        f_mul_g = self.intt_mert(F_mul_G)
        assert f_mul_g == self.mul_schoolbook_x_n_minus_one(f, g)
        # check multiplication mod x^n+1
        fp = self.prec(f)
        gp = self.prec(g)
        Fp = self.ntt_mert(fp)
        Gp = self.ntt_mert(gp)
        Fp_mul_Gp = [(x*y) % q for (x, y) in zip(Fp, Gp)]
        fp_mul_gp = self.intt_mert(Fp_mul_Gp)
        f_mul_g_2 = self.post(fp_mul_gp)
        assert f_mul_g_2 == (Poly(f, q) * Poly(g, q)).coeffs

        # now, the other way around, using ntt and intt instead of ntt_mert and intt_mert
        T = NTT(q)
        assert self.ntt_mert(self.prec(f)) == T.ntt(f)
        assert self.post(self.intt_mert(F)) == T.intt(F)
        assert self.ntt_mert(f) == T.ntt(self.post(f))
        assert self.intt_mert(F) == self.prec(T.intt(F))

        # now we compute the multiplication mod x^n-1 using T.(i)ntt:
        F1 = T.ntt(self.post(f))
        G1 = T.ntt(self.post(g))
        F1_mul_G1 = [(x*y) % q for (x, y) in zip(F1, G1)]
        f1_mul_g1 = self.prec(T.intt(F1_mul_G1))
        assert f1_mul_g1 == f_mul_g
