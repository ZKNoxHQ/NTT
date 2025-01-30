"""This file contains the polynomial arithmetic implementation."""
from random import randint
from ntt import NTT
from ntt_constants import ψ, ψ_inv


class Poly:
    def __init__(self, coeffs, q, polmod=0):
        self.coeffs = coeffs
        self.q = q
        self.NTT = NTT(q)
        self.polmod = polmod
        self.ψ = ψ[q]
        self.ψ_inv = ψ_inv[q]

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
            D[i] = (C[i] - C[i+n]) % self.q
        return Poly(D, self.q, self.polmod)

    def mul_PWC(self, other):
        # TODO DEFINE PWC SOMEWHERE
        """Multiplication of `self` by `other` modulo x^n -1 (and not x^n+1)."""
        f = self.coeffs
        g = other.coeffs
        n = len(f)
        # pre-processing
        # list of roots for the precomputations
        ψ0_inv = self.ψ_inv[::len(self.ψ_inv)//n]
        fp = Poly([(x * y) % self.q for (x, y) in zip(f, ψ0_inv)], self.q)
        gp = Poly([(x * y) % self.q for (x, y) in zip(g, ψ0_inv)], self.q)
        fp_mul_gp = fp*gp
        # post processing
        ψ0 = self.ψ[::len(self.ψ)//n]
        f_mul_g = [(x * y) % self.q for (x, y) in zip(fp_mul_gp.coeffs, ψ0)]
        return Poly(f_mul_g, self.q)

    def mul_schoolbook_PWC(self, other):
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
        # reduction modulo x^n  - 1
        for i in range(n):
            D[i] = (C[i] + C[i+n]) % self.q
        return Poly(D, self.q)

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
