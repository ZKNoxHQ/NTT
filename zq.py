"""This file contains the Zq arithmetic.
"""
from copy import copy
from ntt import ntt, intt, mul_ntt, div_ntt, q


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
