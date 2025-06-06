"""This file contains a recursive implementation of the NTT.

The NTT implemented here is for polynomials in Z_q[x]/(phi), with:
- The integer modulus q = 12 * 1024 + 1 = 12289
- The polynomial modulus phi = x ** n + 1, with n a power of two, n =< 1024

The code is voluntarily very similar to the code of the FFT.
It is probably possible to use templating to merge both implementations.
"""

from polyntt.ntt import NTT
from polyntt.utils import inv_mod
from polyntt.ntt_constants_recursive import roots_dict_mod
from polyntt.m31_2 import mul2, add2, sub2, inv2, sqr1, sqr1_inv, i2

def merge(f_list):
    """Merge two polynomials into a single polynomial f.

    Args:
        f_list: a list of polynomials

    Format: coefficient

    Function from Thomas Prest repository
    """
    f0, f1 = f_list
    n = 2 * len(f0)
    f = [0] * n
    f[::2] = f0
    f[1::2] = f1
    return f


class NTTRecursive(NTT):

    def __init__(self, q):
        """Implements Number Theoretic Transform for fast polynomial multiplication."""
        self.q = q
        # i2 is the inverse of 2 mod q
        self.i2 = inv_mod(2, self.q)
        # sqr1 is a square root of (-1) mod q (currently, sqr1 = 1479)
        self.sqr1 = roots_dict_mod[q][2][0]
        self.roots_dict_mod = roots_dict_mod[q]
        # ratio between degree n and number of complex coefficients of the NTT
        # while here this ratio is 1, it is possible to develop a short NTT such that it is 2.
        self.ntt_ratio = 1
        if q==2**31-1:
            self.ntt = self.ntt_m31
            self.intt = self.intt_m31
            self.split_ntt = self.split_ntt_m31
            self.merge_ntt = self.merge_ntt_m31

    def split_ntt(self, f_ntt):
        """Split a polynomial f in two or three polynomials.

        Args:
            f_ntt: a polynomial

        Format: NTT
        """
        n = len(f_ntt)
        w = self.roots_dict_mod[n]
        f0_ntt = [0] * (n // 2)
        f1_ntt = [0] * (n // 2)
        for i in range(n // 2):
            f0_ntt[i] = (self.i2 * (f_ntt[2 * i] + f_ntt[2 * i + 1])) % self.q
            f1_ntt[i] = (self.i2 * (f_ntt[2 * i] - f_ntt[2 * i + 1])
                         * inv_mod(w[2 * i], self.q)) % self.q
        return [f0_ntt, f1_ntt]

    def merge_ntt(self, f_list_ntt):
        """Merge two or three polynomials into a single polynomial f.

        Args:
            f_list_ntt: a list of polynomials

        Format: NTT
        """
        f0_ntt, f1_ntt = f_list_ntt
        n = 2 * len(f0_ntt)
        w = self.roots_dict_mod[n]
        f_ntt = [0] * n
        for i in range(n // 2):
            f_ntt[2 * i + 0] = (f0_ntt[i] + w[2 * i] * f1_ntt[i]) % self.q
            f_ntt[2 * i + 1] = (f0_ntt[i] - w[2 * i] * f1_ntt[i]) % self.q
        return f_ntt

    def ntt(self, f):
        """Compute the NTT of a polynomial.

        Args:
            f: a polynomial

        Format: input as coefficients, output as NTT
        """
        n = len(f)
        if (n > 2):
            f0, f1 = f[::2], f[1::2]
            f0_ntt = self.ntt(f0)
            f1_ntt = self.ntt(f1)
            f_ntt = self.merge_ntt([f0_ntt, f1_ntt])
        elif (n == 2):
            f_ntt = [0] * n
            f_ntt[0] = (f[0] + self.sqr1 * f[1]) % self.q
            f_ntt[1] = (f[0] - self.sqr1 * f[1]) % self.q
        return f_ntt

    def intt(self, f_ntt):
        """Compute the inverse NTT of a polynomial.

        Args:
            f_ntt: a NTT of a polynomial

        Format: input as NTT, output as coefficients
        """
        n = len(f_ntt)
        if (n > 2):
            f0_ntt, f1_ntt = self.split_ntt(f_ntt)
            f0 = self.intt(f0_ntt)
            f1 = self.intt(f1_ntt)
            f = merge([f0, f1])
        elif (n == 2):
            f = [0] * n
            f[0] = (self.i2 * (f_ntt[0] + f_ntt[1])) % self.q
            f[1] = (self.i2 * inv_mod(self.sqr1, self.q)
                    * (f_ntt[0] - f_ntt[1])) % self.q
        return f

    ###
    # M31 specific implementation
    ###

    def split_ntt_m31(self, f_ntt):
        """Split a polynomial f in two or three polynomials.
        """
        n = len(f_ntt)
        w = self.roots_dict_mod[n]
        f0_ntt = [0] * (n // 2)
        f1_ntt = [0] * (n // 2)
        for i in range(n // 2):
            f0_ntt[i] = mul2(add2(f_ntt[2*i], f_ntt[2*i+1]), [i2, 0])
            inv_w_2i = inv2(w[2*i])
            f1_ntt[i] = mul2(mul2([i2, 0], inv_w_2i),
                            sub2(f_ntt[2*i], f_ntt[2*i+1]))
        return [f0_ntt, f1_ntt]


    def merge_ntt_m31(self, f_list_ntt):
        """Merge two polynomials into a single polynomial f.
        """
        f0_ntt, f1_ntt = f_list_ntt
        n = len(f0_ntt)
        w = self.roots_dict_mod[2*n]
        f_ntt = [0] * 2*n
        for i in range(n):
            f_ntt[2 * i + 0] = add2(f0_ntt[i],
                                    mul2(w[2 * i], f1_ntt[i]))
            f_ntt[2 * i + 1] = sub2(f0_ntt[i],
                                    mul2(w[2 * i], f1_ntt[i]))
        return f_ntt


    def ntt_m31(self, f):
        """Compute the NTT of a polynomial.
        """
        # first, convert into Fp² elements
        # if needed
        if not(isinstance(f[0], list)):
            # this map is Fp[x]/(x^n+1) -> Fp²[y]/(y^{n/2}+1)
            # by x -> ωy
            a = f
            r = []
            n = len(a)//2
            ω = self.roots_dict_mod[2*n][0]
            ω_i = [1, 0]
            for i in range(n):
                r.append(mul2([a[i], self.q-a[i+n]], ω_i))
                ω_i = mul2(ω_i, ω)
        else:
            r = [_ for _ in f]
        n = len(r)
        if (n > 2):
            f0, f1 = r[::2], r[1::2]
            f0_ntt = self.ntt_m31(f0)
            f1_ntt = self.ntt_m31(f1)
            f_ntt = self.merge_ntt([f0_ntt, f1_ntt])
        elif (n == 2):
            f_ntt = [0] * n
            f_ntt[0] = add2(r[0], mul2(r[1], sqr1))
            f_ntt[1] = sub2(r[0], mul2(r[1], sqr1))
        return f_ntt


    def intt_m31(self, f_ntt):
        """Compute the inverse NTT of a polynomial.
        """
        # first, iNTT in Fp²
        n = len(f_ntt)
        if (n > 2):
            f0_ntt, f1_ntt = self.split_ntt(f_ntt)
            f0 = self.intt_m31(f0_ntt)
            f1 = self.intt_m31(f1_ntt)
            f = merge([f0, f1])
        elif (n == 2):
            f = [0] * n
            f[0] = mul2(add2(f_ntt[0], f_ntt[1]), [i2, 0])
            f[1] = mul2(sub2(f_ntt[0], f_ntt[1]), mul2([i2, 0], sqr1_inv))
        # second, convert back to Fp
        # if needed
        if len(f) == 2**6: # TODO FIX THIS BECAUSE ITS NOT WORKING FROM A GENERAL POINT OF VIEW
            a = [_ for _ in f]
            r = []
            s = []
            n = len(a)
            ω = self.roots_dict_mod[2*n][0]
            ω_inv = inv2(ω)
            ω_inv_i = [1, 0]
            for i in range(n):
                c = mul2(a[i], ω_inv_i)
                ω_inv_i = mul2(ω_inv_i, ω_inv)
                r.append(c[0])
                s.append(self.q-c[1])
            return r+s
        else:
            return f