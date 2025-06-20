class FieldElement:
    def __init__(self, field):
        self.field = field
        self.coeffs = []

    def __iter__(self):
        return iter(self.coeffs)

    def __add__(self, other): raise NotImplementedError
    def __sub__(self, other): raise NotImplementedError
    def __mul__(self, other): raise NotImplementedError
    def __truediv__(self, other): raise NotImplementedError
    def inverse(self): raise NotImplementedError
    def __neg__(self): return NotImplementedError
    def __eq__(self, other): return NotImplementedError


class Field:
    def __call__(self, coeffs): raise NotImplementedError
    def extension(self, degree): raise NotImplementedError
    def random(self): raise NotImplementedError
    def zero(self): raise NotImplementedError

    def batch_inversion(self, elements):
        """Compute batch inversion of a list of field elements."""
        n = len(elements)
        if n == 0:
            return []
        # Prefix products
        prefix = [None] * n
        prefix[0] = elements[0]
        for i in range(1, n):
            prefix[i] = prefix[i - 1] * elements[i]
        # Inverse of the total product
        total_inv = prefix[-1].inverse()
        # Individual inverses using the prefix products
        inverses = [None] * n
        inverses[-1] = total_inv
        for i in range(n - 2, -1, -1):
            inverses[i] = inverses[i + 1] * elements[i + 1]
        # Final inverses
        for i in range(1, n):
            inverses[i] = inverses[i] * prefix[i - 1]
        return inverses
