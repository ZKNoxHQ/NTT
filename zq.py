"""This file contains the Zq arithmetic.
"""
from copy import copy
from ntt import ntt, intt, mul_ntt, div_ntt, q
from ntt_constants import *


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
    f_ntt = ntt(f)
    g_ntt = ntt(g)
    return intt(mul_ntt(f_ntt, g_ntt))


def mul_schoolbook_zq(f, g):
    """Multiplication of two polynomials using the schoolbook algorithm."""
    n = len(f)
    assert n == len(g)
    C = [0] * (2 * n)
    D = [0] * (n)
    for j, f_j in enumerate(f):
        for k, g_k in enumerate(g):
            C[j+k] = (C[j+k] + f_j * g_k) % q
    # reduction modulo x^n  + 1
    for i in range(n):
        D[i] = (C[i] - C[i+n]) % q
    return D


def div_zq(f, g):
    """Division of two polynomials (coefficient representation)."""
    try:
        return intt(div_ntt(ntt(f), ntt(g)))
    except ZeroDivisionError:
        raise

# def adj(f):
#     """Ajoint of a polynomial (coefficient representation)."""
#     return intt(adj_ntt(ntt(f)))
