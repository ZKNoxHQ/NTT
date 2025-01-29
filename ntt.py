"""This file contains an iterative implementation of the NTT.

The NTT implemented here is for polynomials in Z_q[x]/(phi), with:
- The integer modulus q = 12 * 1024 + 1 = 12289
- The polynomial modulus phi = x ** n + 1, with n a power of two, n =< 1024
"""
from copy import copy
from ntt_constants import ψ_12289_rev, ψ_12289_inv_rev, n_12289_inv

q = 12*1024 + 1
ψ_rev, ψ_inv_rev, n_inv = ψ_12289_rev, ψ_12289_inv_rev, n_12289_inv


def ntt(f):
    # following eprint 2016/504 Algorithm 1
    a = copy(f)
    n = len(a)
    t = n
    m = 1
    while m < n:
        t //= 2
        for i in range(m):
            j1 = 2*i*t
            j2 = j1+t-1
            S = ψ_rev[m+i]
            for j in range(j1, j2+1):
                U = a[j]
                V = a[j+t]*S
                a[j] = (U+V) % q
                a[j+t] = (U-V) % q
        m = 2*m
    return a


def intt(f_ntt):
    a = copy(f_ntt)
    # following eprint 2016/504 Algorithm 2
    n = len(a)
    t = 1
    m = n
    while m > 1:
        j1 = 0
        h = m//2
        for i in range(h):
            j2 = j1+t-1
            S = ψ_inv_rev[h+i]
            for j in range(j1, j2+1):
                U = a[j]
                V = a[j+t]
                a[j] = (U+V) % q
                a[j+t] = ((U-V) * S) % q
            j1 += 2*t
        t *= 2
        m //= 2
    for j in range(n):
        a[j] = (a[j] * n_inv[n]) % q
    return a


def add_zq(f, g):
    """Addition of two polynomials (coefficient representation)."""
    assert len(f) == len(g)
    deg = len(f)
    return [(f[i] + g[i]) % q for i in range(deg)]


def neg_zq(f):
    """Negation of a polynomials (any representation)."""
    deg = len(f)
    return [(- f[i]) % q for i in range(deg)]


def sub_zq(f, g):
    """Substraction of two polynomials (any representation)."""
    return add_zq(f, neg_zq(g))


def mul_zq(f, g):
    """Multiplication of two polynomials (coefficient representation)."""
    return intt(mul_ntt(ntt(f), ntt(g)))


def div_zq(f, g):
    """Division of two polynomials (coefficient representation)."""
    try:
        return intt(div_ntt(ntt(f), ntt(g)))
    except ZeroDivisionError:
        raise


# def adj(f):
#     """Ajoint of a polynomial (coefficient representation)."""
#     return intt(adj_ntt(ntt(f)))


def add_ntt(f_ntt, g_ntt):
    """Addition of two polynomials (NTT representation)."""
    return add_zq(f_ntt, g_ntt)


def sub_ntt(f_ntt, g_ntt):
    """Substraction of two polynomials (NTT representation)."""
    return sub_zq(f_ntt, g_ntt)


def mul_ntt(f_ntt, g_ntt):
    """Multiplication of two polynomials (coefficient representation)."""
    assert len(f_ntt) == len(g_ntt)
    deg = len(f_ntt)
    return [(f_ntt[i] * g_ntt[i]) % q for i in range(deg)]


def div_ntt(f_ntt, g_ntt):
    """Division of two polynomials (NTT representation)."""
    assert len(f_ntt) == len(g_ntt)
    deg = len(f_ntt)
    if any(elt == 0 for elt in g_ntt):
        raise ZeroDivisionError
    return [(f_ntt[i] * pow(g_ntt[i], -1, q)) % q for i in range(deg)]


# def adj_ntt(f_ntt):
#     """Ajoint of a polynomial (NTT representation)."""
#     deg = len(f_ntt)
#     return [f_ntt[i].conjugate() for i in range(deg)]

"""This value is the ratio between:
    - The degree n
    - The number of complex coefficients of the NTT
While here this ratio is 1, it is possible to develop a short NTT such that it is 2.
"""
ntt_ratio = 1
