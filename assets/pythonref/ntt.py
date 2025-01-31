"""This file contains an iterative implementation of the NTT.

The NTT implemented here is for polynomials in Z_q[x]/(phi), with:
- The integer modulus q = 12 * 1024 + 1 = 12289
- The polynomial modulus phi = x ** n + 1, with n a power of two, n =< 1024
"""
from ntt_constants import *


class NTT:

    def __init__(self, q):
        """Implements Number Theoretic Transform for fast polynomial multiplication."""
        self.q = q
        self.ψ_rev = ψ_rev[q]  # TODO
        self.ψ_inv_rev = ψ_inv_rev[q]  # TODO
        """This value is the ratio between:
        - The degree n
        - The number of complex coefficients of the NTT
        While here this ratio is 1, it is possible to develop a short NTT such that it is 2.
        """
        self.ntt_ratio = 1

    def ntt(self, f):
        # following eprint 2016/504 Algorithm 1
        a = [_ for _ in f]
        n = len(a)
        t = n
        m = 1
        while m < n:
            t //= 2
            for i in range(m):
                j1 = 2*i*t
                j2 = j1+t-1
                S = self.ψ_rev[m+i]
                for j in range(j1, j2+1):
                    U = a[j]
                    V = a[j+t]*S
                    a[j] = (U+V) % self.q
                    a[j+t] = (U-V) % self.q
            m = 2*m
        return a

    def intt(self, f_ntt):
        # following eprint 2016/504 Algorithm 2
        a = [_ for _ in f_ntt]
        n = len(a)
        t = 1
        m = n
        while m > 1:
            j1 = 0
            h = m//2
            for i in range(h):
                j2 = j1+t-1
                S = self.ψ_inv_rev[h+i]
                for j in range(j1, j2+1):
                    U = a[j]
                    V = a[j+t]
                    a[j] = (U+V) % self.q
                    a[j+t] = ((U-V) * S) % self.q
                j1 += 2*t
            t *= 2
            m //= 2
        for j in range(n):
            a[j] = (a[j] * n_inv[self.q][n]) % self.q
        return a

    # def add_ntt(self, f_ntt, g_ntt):
    #     """Addition of two polynomials (NTT representation)."""
    #     return add_zq(f_ntt, g_ntt)

    # def sub_ntt(self, f_ntt, g_ntt):
    #     """Substraction of two polynomials (NTT representation)."""
    #     return sub_zq(f_ntt, g_ntt)

    def mul_ntt(self, f_ntt, g_ntt):
        """Multiplication of two polynomials (coefficient representation)."""
        assert len(f_ntt) == len(g_ntt)
        deg = len(f_ntt)
        return [(f_ntt[i] * g_ntt[i]) % self.q for i in range(deg)]

    def div_ntt(self, f_ntt, g_ntt):
        """Division of two polynomials (NTT representation)."""
        assert len(f_ntt) == len(g_ntt)
        deg = len(f_ntt)
        if any(elt == 0 for elt in g_ntt):
            raise ZeroDivisionError
        return [(f_ntt[i] * pow(g_ntt[i], -1, self.q)) % self.q for i in range(deg)]

    # def adj_ntt(self, f_ntt):
    #     """Ajoint of a polynomial (NTT representation)."""
    #     deg = len(f_ntt)
    #     return [f_ntt[i].conjugate() for i in range(deg)]
