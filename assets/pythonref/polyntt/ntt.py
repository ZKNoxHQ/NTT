class NTT:
    """Base class for Number Theoretic Transform"""

    def ntt(self, poly):
        raise NotImplementedError("Subclasses must implement NTT")

    def intt(self, poly):
        raise NotImplementedError(
            "Subclasses must implement inverse NTT")

    def vec_add(self, f_ntt, g_ntt):
        """Addition of two polynomials (NTT representation)."""
        return [(x+y) % self.q for (x, y) in zip(f_ntt, g_ntt)]

    def vec_sub(self, f_ntt, g_ntt):
        """Substraction of two polynomials (NTT representation)."""
        return self.vec_add(f_ntt, [(-x) % self.q for x in g_ntt])

    def vec_mul(self, f_ntt, g_ntt):
        """Multiplication of two polynomials (NTT representation)."""
        assert len(f_ntt) == len(g_ntt)
        deg = len(f_ntt)
        return [(f_ntt[i] * g_ntt[i]) % self.q for i in range(deg)]

    def vec_div(self, f_ntt, g_ntt):
        """Division of two polynomials (NTT representation)."""
        assert len(f_ntt) == len(g_ntt)
        deg = len(f_ntt)
        if any(elt == 0 for elt in g_ntt):
            raise ZeroDivisionError
        inv_g_ntt = batch_modular_inversion(g_ntt, self.q)
        return self.vec_mul(f_ntt, inv_g_ntt)


def batch_modular_inversion(elements, q):
    """Compute batch inversion of a list of elements mod q."""
    n = len(elements)
    if n == 0:
        return []

    # Prefix products
    prefix = [None] * n
    prefix[0] = elements[0]
    for i in range(1, n):
        prefix[i] = (prefix[i - 1] * elements[i]) % q

    # Iinverse of the total product
    total_inv = pow(prefix[-1], -1, q)

    # Individual inverses using the prefix products
    inverses = [None] * n
    inverses[-1] = total_inv
    for i in range(n - 2, -1, -1):
        inverses[i] = (inverses[i + 1] * elements[i + 1]) % q

    # Final inverses
    for i in range(1, n):
        inverses[i] = (inverses[i] * prefix[i - 1]) % q

    return inverses
