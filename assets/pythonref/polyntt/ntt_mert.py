"""This file contains an iterative implementation of the NTT using Mert implementation.
It is optimized for x^n-1, but requires pre- and post-computations for ot x^n+1.
We provide it here for comparison.
"""
from polyntt.ntt_constants import n_inv, ψ_rev, ψ_inv_rev, ψ as Ψ, ψ_inv as Ψ_inv
from polyntt.poly import Poly
from math import log
from polyntt.ntt_iterative import NTTIterative
from polyntt.ntt import NTT


class NTTMert(NTT):

    def __init__(self, q):
        """Implements Number Theoretic Transform for fast polynomial multiplication."""
        self.q = q
        self.ψ = ψ_rev[q]
        self.ψ_inv = ψ_inv_rev[q]
        self.ψ_rev = ψ_rev[q]
        self.ψ_inv_rev = ψ_inv_rev[q]
        # ratio between degree n and number of complex coefficients of the NTT
        # while here this ratio is 1, it is possible to develop a short NTT such that it is 2.
        self.ntt_ratio = 1

    def intReverse(self, a, n):
        b = ('{:0'+str(n)+'b}').format(a)
        return int(b[::-1], 2)

    def ntt(self, A):
        # from https://github.com/acmert/ntt-based-polmul/blob/3476c2d9f22eac7696f2785eba036ffae1b4460e/baseline/ntt.py#L206
        N = len(A)
        B = [_ for _ in A]

        for s in range(int(log(N, 2)), 0, -1):
            m = 2**s
            for k in range(int(N/m)):
                # TW = pow(ω, self.intReverse(k, int(log(N, 2))-s)*int(m/2), self.q)
                TW = ψ_rev[self.q][k]
                for j in range(int(m/2)):
                    u = B[k*m+j]
                    t = (TW*B[k*m+j+int(m/2)]) % self.q

                    B[k*m+j] = (u+t) % self.q
                    B[k*m+j+int(m/2)] = (u-t) % self.q

        return B

    def intt(self, A):
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
                    # TW = pow(ω_inv, self.intReverse(jt, int(log(N >> 1, 2))), self.q)
                    TW = ψ_inv_rev[self.q][jt]

                    temp = B[j]

                    B[j] = (temp + B[j+d]) % self.q
                    B[j+d] = (temp - B[j+d])*TW % self.q

                    jt = jt+1
                    j = j + 2*m
            m = 2*m
            v = int(v/2)
            d = 2*d

        return [(x*n_inv[self.q][N]) % self.q for x in B]

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
                C[j+k] = (C[j+k] + f_j * g_k) % self.q
        # reduction modulo x^n  + 1
        for i in range(n):
            D[i] = (C[i] + C[i+n]) % self.q
        return D

    def prec(self, A):
        n = len(A)
        # list of roots for the precomputations
        Ψ0 = Ψ[self.q][::len(Ψ[self.q])//n]
        return [(x * y) % self.q for (x, y) in zip(A, Ψ0)]

    def post(self, A):
        n = len(A)
        # list of roots for the precomputations
        Ψ0_inv = Ψ_inv[self.q][::len(Ψ_inv[self.q])//n]
        return [(x * y) % self.q for (x, y) in zip(A, Ψ0_inv)]

    def poly_mul_mert(self, A, B):
        A_ntt = self.ntt(self.prec(A))
        B_ntt = self.ntt(self.prec(B))
        A_B_ntt = [(x*y) % self.q for (x, y) in zip(A_ntt, B_ntt)]
        return self.post(self.intt(A_B_ntt))


q = 3329
Mert = NTTMert(q)
f = [1, 2, 3, 4]
# f = [randint(0, q-1) for i in range(4)]
F = Mert.ntt(f)
back_f = Mert.intt(F)
for i in range(len(f)):
    assert back_f[i] == f[i]
g = [5, 6, 7, 8]
# g = [randint(0, q-1) for i in range(4)]
G = Mert.ntt(g)
back_g = Mert.intt(G)
for i in range(len(g)):
    assert back_g[i] == g[i]

# check multiplication mod x^n-1
F_mul_G = [(x*y) % q for (x, y) in zip(F, G)]
f_mul_g = Mert.intt(F_mul_G)
assert f_mul_g == Mert.mul_schoolbook_x_n_minus_one(f, g)
# check multiplication mod x^n+1
fp = Mert.prec(f)
gp = Mert.prec(g)
Fp = Mert.ntt(fp)
Gp = Mert.ntt(gp)
Fp_mul_Gp = [(x*y) % q for (x, y) in zip(Fp, Gp)]
fp_mul_gp = Mert.intt(Fp_mul_Gp)
f_mul_g_2 = Mert.post(fp_mul_gp)
assert f_mul_g_2 == (Poly(f, q) * Poly(g, q)).coeffs

# now, the other way around, using ntt and intt instead of ntt_mert and intt_mert
T = NTTIterative(q)
assert Mert.ntt(Mert.prec(f)) == T.ntt(f)
assert Mert.post(Mert.intt(F)) == T.intt(F)
assert Mert.ntt(f) == T.ntt(Mert.post(f))
assert Mert.intt(F) == Mert.prec(T.intt(F))

# now we compute the multiplication mod x^n-1 using T.(i)ntt:
F1 = T.ntt(Mert.post(f))
G1 = T.ntt(Mert.post(g))
F1_mul_G1 = [(x*y) % q for (x, y) in zip(F1, G1)]
f1_mul_g1 = Mert.prec(T.intt(F1_mul_G1))
assert f1_mul_g1 == f_mul_g
