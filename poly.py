"""This file contains the polynomial arithmetic implementation."""
from random import randint
from ntt import NTT
from ntt_constants import *


class Poly:
    def __init__(self, coeffs, q, polmod):
        self.coeffs = coeffs
        self.q = q
        self.NTT = NTT(q)
        self.polmod = polmod

    def __eq__(self, other):
        for (a, b) in zip(self.coeffs, other.coeffs):
            if (a-b) % self.q != 0:
                return False
        return True

    def __add__(self, other):
        f = self.coeffs
        g = other.coeffs
        assert len(f) == len(g)
        deg = len(f)
        return Poly([(f[i] + g[i]) % self.q for i in range(deg)], self.q, self.polmod)

    def __neg__(self):
        """Negation of a polynomials (any representation)."""
        f = self.coeffs
        deg = len(f)
        return Poly([(- f[i]) % self.q for i in range(deg)], self.q, self.polmod)

    def __sub__(self, other):
        """Substraction of two polynomials (any representation)."""
        return self + (-other)

    def __mul__(self, other):
        """Multiplication of two polynomials (coefficient representation)."""
        f = self.coeffs
        g = other.coeffs
        T = self.NTT
        f_ntt = T.ntt(f)
        g_ntt = T.ntt(g)
        return Poly(T.intt(T.mul_ntt(f_ntt, g_ntt)), self.q, self.polmod)

    def mul_schoolbook(self, other):
        """Multiplication of two polynomials using the schoolbook algorithm."""
        f = self.coeffs
        g = other.coeffs
        n = len(f)
        assert n == len(g)
        C = [0] * (2 * n)
        D = [0] * (n)
        for j, f_j in enumerate(f):
            for k, g_k in enumerate(g):
                C[j+k] = (C[j+k] + f_j * g_k) % self.q
        # reduction modulo x^n  + 1
        for i in range(n):
            D[i] = (C[i] - C[i+n]) % self.q  # TODO CHANGE DEPENDING ON POLMOD
        return Poly(D, self.q, self.polmod)

    def div(self, other):
        """Division of two polynomials (coefficient representation)."""
        try:
            f = self.coeffs
            g = other.coeffs
            T = self.NTT
            f_ntt = T.ntt(f)
            g_ntt = T.ntt(g)
            return Poly(T.intt(T.div_ntt(f_ntt, g_ntt)), self.q, self.polmod)
        except ZeroDivisionError:
            raise

    # def adj(f):
    #     """Ajoint of a polynomial (coefficient representation)."""
    #     return intt(adj_ntt(ntt(f)))
