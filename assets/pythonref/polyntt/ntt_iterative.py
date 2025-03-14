"""This file contains an iterative implementation of the NTT.

The NTT implemented here is for polynomials in Z_q[x]/(phi), with:
- The integer modulus q = 12 * 1024 + 1 = 12289
- The polynomial modulus phi = x ** n + 1, with n a power of two, n =< 1024
"""
from polyntt.ntt_constants_iterative import *
from polyntt.ntt import NTT


class NTTIterative(NTT):

    def __init__(self, q):
        """Implements Number Theoretic Transform for fast polynomial multiplication."""
        self.q = q
        # can be removed if run for nodes, increases efficiency of Poly.mul_pwc with a larger storage
        self.ψ = ψ[q]
        # can be removed if run for nodes, increases efficiency of Poly.mul_pwc with a larger storage
        self.ψ_inv = ψ_inv[q]
        # useful for efficiency (even in nodes)
        self.ψ_rev = ψ_rev[q]
        # useful for efficiency (even in nodes)
        self.ψ_inv_rev = ψ_inv_rev[q]
        # ratio between degree n and number of complex coefficients of the NTT
        # while here this ratio is 1, it is possible to develop a short NTT such that it is 2.
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

    def ntt_without_mod(self, f):
        # following eprint 2016/504 Algorithm 1
        a = [_ for _ in f]
        n = len(a)
        t = n
        m = 1
        ell = 0
        while m < n:
            t //= 2
            for i in range(m):
                j1 = 2*i*t
                j2 = j1+t-1
                S = self.ψ_rev[m+i]
                for j in range(j1, j2+1):
                    U = a[j]
                    V = a[j+t]*S
                    V_neg = (self.q << ((ell+1)*15)) - V
                    a[j] = (U+V)
                    a[j+t] = (U+V_neg)
                    assert a[j+t] < (1 << 256)
                    assert a[j+t] >= 0
            m = 2*m
            ell += 1
        return [elt % self.q for elt in a]

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

    def intt_without_mod(self, f_ntt):
        # following eprint 2016/504 Algorithm 2
        # with modular reduction only in the final loop
        # at every step, U,V < 2**ell * q**(ell+1)
        # using powers of 2, U,V < 2**ell * 2**(14*ell) * q
        # so -V = q << (15*ell)
        a = [_ for _ in f_ntt]
        n = len(a)
        t = 1
        m = n
        ell = 0
        while m > 1:
            j1 = 0
            h = m//2
            for i in range(h):
                j2 = j1+t-1
                S = self.ψ_inv_rev[h+i]
                for j in range(j1, j2+1):
                    U = a[j]
                    V = a[j+t]
                    V_neg = (self.q << (ell*15)) - V
                    a[j] = (U+V)
                    a[j+t] = ((U+V_neg) * S)
                    assert a[j+t] >= 0
                    assert a[j+t] < (1 << 256)
                j1 += 2*t
            t *= 2
            m //= 2
            ell += 1
        for j in range(n):
            a[j] = (a[j] * n_inv[self.q][n]) % self.q
        return a
