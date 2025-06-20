"""This file contains an iterative implementation of the NTT.

The NTT implemented here is for polynomials in Z_q[x]/(x^n+1),
with n a power of two, n =< 1024
"""
from polyntt.polynomial.ntt_constants import ψ_rev, ψ_inv_rev, n_inv
from polyntt.polynomial.polynomial import PolynomialRing, Polynomial


class PolynomialRingNTT(PolynomialRing):
    def __init__(self, F, n):
        """Implements Number Theoretic Transform for fast polynomial multiplication."""
        self.F = F
        self.n = n
        self.element = PolynomialNTT
        # useful for efficiency (even in nodes)
        self.ψ_rev = ψ_rev[F.p]
        # useful for efficiency (even in nodes)
        self.ψ_inv_rev = ψ_inv_rev[F.p]
        # inverse of n mod q for intt
        self.n_inv = n_inv[F.p]

    # def __call__(self, coefficients, ntt=False):
    #     if isinstance(coefficients, int):
    #         return self.element(self, [coefficients], ntt=ntt)
    #     if not isinstance(coefficients, list):
    #         raise TypeError(
    #             f"Polynomials should be constructed from a list of integers, of length at most d = {self.n}"
    #         )
    #     return self.element(self, coefficients, ntt=ntt)

    def ntt(self, f):
        # following eprint 2016/504 Algorithm 1
        a = [elt for elt in f.coeffs]
        n = len(a)
        t = n
        m = 1
        while m < n:
            t //= 2
            for i in range(m):
                j1 = 2*i*t
                j2 = j1+t-1
                S = self.F(self.ψ_rev[m+i])
                for j in range(j1, j2+1):
                    U = a[j]
                    V = a[j+t]*S
                    a[j] = (U+V)
                    a[j+t] = (U-V)
            m = 2*m
        return self(a, ntt=True)

    def intt(self, f_ntt):
        # following eprint 2016/504 Algorithm 2
        a = [elt for elt in f_ntt.coeffs]
        n = len(a)
        t = 1
        m = n
        while m > 1:
            j1 = 0
            h = m//2
            for i in range(h):
                j2 = j1+t-1
                S = self.F(self.ψ_inv_rev[h+i])
                for j in range(j1, j2+1):
                    U = a[j]
                    V = a[j+t]
                    a[j] = (U+V)
                    a[j+t] = ((U-V) * S)
                j1 += 2*t
            t *= 2
            m //= 2
        for j in range(n):
            a[j] = a[j] * self.F(self.n_inv[n])
        return self(a, ntt=False)


class PolynomialNTT(Polynomial):
    def __init__(self, parent, coeffs, ntt=False):
        super().__init__(parent, coeffs)
        self.ntt = ntt

    def __repr__(self):
        if self.ntt:
            gen = 'ζ'
        else:
            gen = 'X'

        return " + ".join(
            f"({c})" if i == 0 else (
                f"({c})·{gen}^{i}" if i > 1 else f"({c})·{gen}")
            for i, c in enumerate(self.coeffs) if c != self.parent.F(0)
        ) or "0"

    def __iter__(self):
        return iter(self.coeffs)

    def __mul__(self, other):
        assert self.ntt == other.ntt
        if self.ntt:
            n = self.parent.n
            a = self.coeffs
            while len(a) < n:
                a += [self.parent.F.zero()]
            b = other.coeffs
            while len(b) < n:
                b += [self.parent.F.zero()]
            new_coeffs = [self.parent.F.zero() for _ in range(n)]
            for i in range(n):
                new_coeffs[i] = a[i] * b[i]
            return self.parent(new_coeffs, ntt=True)
        else:
            # multiplication using ntt and intt
            self_ntt = self.parent.ntt(self)
            other_ntt = other.parent.ntt(other)
            return self.parent.intt(self_ntt * other_ntt)

    def div(self, other):
        assert self.ntt == other.ntt

        if self.ntt:
            n = self.parent.n
            a = self.coeffs
            while len(a) < n:
                a += [self.parent.F.zero()]
            b = other.coeffs
            while len(b) < n:
                b += [self.parent.F.zero()]
            new_coeffs = [self.parent.F.zero() for _ in range(n)]
            for i in range(n):
                new_coeffs[i] = a[i] / b[i]
            return self.parent(new_coeffs)
        else:
            return self * other.inverse()

    def inverse(self):
        if self.ntt:
            a = self.coeffs
            try:
                a_invs = self.F.batch_inversion(a)
            except ZeroDivisionError:
                raise
            return self.parent(a_invs)
        else:
            print("Not implemented yet")
            return ZeroDivisionError

    def __call__(self, x):
        return NotImplementedError

    def __eq__(self, other):
        return self.coeffs == other.coeffs  # and self.parent == other.parent

    def __str__(self):
        return self.__repr__()
