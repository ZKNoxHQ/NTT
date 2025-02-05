"""This file contains a recursive implementation of the NTT.

The NTT implemented here is for polynomials in Z_q[x]/(phi), with:
- The integer modulus q = 12 * 1024 + 1 = 12289
- The polynomial modulus phi = x ** n + 1, with n a power of two, n =< 1024
"""


"""This file contains an implementation of the NTT.

The NTT implemented here is for polynomials in Z_q[x]/(phi), with:
- The integer modulus q = 12 * 1024 + 1 = 12289
- The polynomial modulus phi = x ** n + 1, with n a power of two, n =< 1024

The code is voluntarily very similar to the code of the FFT.
It is probably possible to use templating to merge both implementations.
"""


"""i2 is the inverse of 2 mod q."""

from ntt_constants_recursive import roots_dict_Zq
i2 = 6145
q = 12 * 1024 + 1

""" sqr1 is a square root of (-1) mod q (currently, sqr1 = 1479)."""
sqr1 = roots_dict_Zq[2][0]


def split(f):
    """Split a polynomial f in two polynomials.

    Args:
        f: a polynomial

    Format: coefficient

    Function from Thomas Prest repository
    """
    n = len(f)
    f0 = [f[2 * i + 0] for i in range(n // 2)]
    f1 = [f[2 * i + 1] for i in range(n // 2)]
    return [f0, f1]


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
    for i in range(n // 2):
        f[2 * i + 0] = f0[i]
        f[2 * i + 1] = f1[i]
    return f


def xgcd(a, b):
    """ Returns gcd(a, b), and x, y such that ax + by = gcd(a, b) """
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def inv_mod_q(elt):
    """
    Thomas Prest stores the inverses mod q, but in the long term, we will consider a larger q,
    and thus we do not store the inverses mod q (it would require a too large storage).
    """
    _, inv_elt, _ = xgcd(elt, q)
    assert (inv_elt * elt) % q == 1
    return inv_elt


class NTTRecursive:

    def __init__(self, q):
        """Implements Number Theoretic Transform for fast polynomial multiplication."""
        self.q = q
        # ratio between degree n and number of complex coefficients of the NTT
        # while here this ratio is 1, it is possible to develop a short NTT such that it is 2.
        self.ntt_ratio = 1

    def split_ntt(self, f_ntt):
        """Split a polynomial f in two or three polynomials.

        Args:
            f_ntt: a polynomial

        Format: NTT
        """
        n = len(f_ntt)
        w = roots_dict_Zq[n]
        f0_ntt = [0] * (n // 2)
        f1_ntt = [0] * (n // 2)
        for i in range(n // 2):
            f0_ntt[i] = (i2 * (f_ntt[2 * i] + f_ntt[2 * i + 1])) % q
            f1_ntt[i] = (i2 * (f_ntt[2 * i] - f_ntt[2 * i + 1])
                         * inv_mod_q(w[2 * i])) % q
        return [f0_ntt, f1_ntt]

    def merge_ntt(self, f_list_ntt):
        """Merge two or three polynomials into a single polynomial f.

        Args:
            f_list_ntt: a list of polynomials

        Format: NTT
        """
        f0_ntt, f1_ntt = f_list_ntt
        n = 2 * len(f0_ntt)
        w = roots_dict_Zq[n]
        f_ntt = [0] * n
        for i in range(n // 2):
            f_ntt[2 * i + 0] = (f0_ntt[i] + w[2 * i] * f1_ntt[i]) % q
            f_ntt[2 * i + 1] = (f0_ntt[i] - w[2 * i] * f1_ntt[i]) % q
        return f_ntt

    def ntt(self, f):
        """Compute the NTT of a polynomial.

        Args:
            f: a polynomial

        Format: input as coefficients, output as NTT
        """
        n = len(f)
        if (n > 2):
            f0, f1 = split(f)
            f0_ntt = self.ntt(f0)
            f1_ntt = self.ntt(f1)
            f_ntt = self.merge_ntt([f0_ntt, f1_ntt])
        elif (n == 2):
            f_ntt = [0] * n
            f_ntt[0] = (f[0] + sqr1 * f[1]) % q
            f_ntt[1] = (f[0] - sqr1 * f[1]) % q
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
            f[0] = (i2 * (f_ntt[0] + f_ntt[1])) % q
            f[1] = (i2 * inv_mod_q(sqr1) * (f_ntt[0] - f_ntt[1])) % q
        return f

    def vec_add(self, f_ntt, g_ntt):
        """Addition of two polynomials(NTT representation)."""
        return [(x+y) % self.q for (x, y) in zip(f_ntt, g_ntt)]

    def vec_sub(self, f_ntt, g_ntt):
        """Substraction of two polynomials(NTT representation)."""
        return self.vec_add(f_ntt, [(-x) % self.q for x in g_ntt])

    def vec_mul(self, f_ntt, g_ntt):
        """Multiplication of two polynomials(coefficient representation)."""
        assert len(f_ntt) == len(g_ntt)
        deg = len(f_ntt)
        return [(f_ntt[i] * g_ntt[i]) % self.q for i in range(deg)]

    def vec_div(self, f_ntt, g_ntt):
        """Division of two polynomials(NTT representation)."""
        assert len(f_ntt) == len(g_ntt)
        deg = len(f_ntt)
        if any(elt == 0 for elt in g_ntt):
            raise ZeroDivisionError
        inv_g_ntt = [pow(g_ntt[i], -1, self.q) % self.q for i in range(deg)]
        return self.vec_mul(f_ntt, inv_g_ntt)
