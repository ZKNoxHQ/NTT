from polyntt.field.field import Field, FieldElement
from polyntt.field.extension_field import ExtensionField
from polyntt.utils import xgcd
from random import randint


class PrimeField(Field):
    def __init__(self, p):
        assert self._is_prime(p)
        self.p = p
        self.degree = 1

    def __call__(self, coeffs):
        if isinstance(coeffs, PrimeFieldElement):
            return coeffs
        return PrimeFieldElement(coeffs % self.p, self)

    def extension(self, degree, alpha=None):
        if degree == 1:
            return self
        return ExtensionField(self, degree, alpha)

    def _is_prime(self, n):
        # found on https://stackoverflow.com/questions/15285534/isprime-function-for-python-language
        if n == 2 or n == 3:
            return True
        if n < 2 or n % 2 == 0:
            return False
        if n < 9:
            return True
        if n % 3 == 0:
            return False
        r = int(n**0.5)
        # since all primes > 3 are of the form 6n ± 1
        # start with f=5 (which is prime)
        # and test f, f+2 for being prime
        # then loop by 6.
        f = 5
        while f <= r:
            if n % f == 0:
                return False
            if n % (f+2) == 0:
                return False
            f += 6
        return True

    def random(self):
        return self(randint(0, self.p-1))

    def zero(self):
        return self(0)


class PrimeFieldElement(FieldElement):
    def __init__(self, coeffs, field):
        super().__init__(field)
        self.coeffs = coeffs % field.p

    def __iter__(self):
        return iter([self.coeffs])

    def __add__(self, other):
        return self.field(self.coeffs + other.coeffs)

    def __sub__(self, other):
        return self.field(self.coeffs - other.coeffs)

    def __mul__(self, other):
        return self.field(self.coeffs * other.coeffs)

    def mul_by_small(self, alpha):
        # useful for extension field
        # alpha is an int
        if alpha == 0:
            return self.field(0)
        if alpha > 0:
            ε = 1
            nalpha = alpha
        else:
            ε = -1
            nalpha = -alpha
        res = self
        for i in range(1, nalpha):
            res += self
        return res

    def __truediv__(self, other):
        return self * other.inverse()

    def __neg__(self):
        return self.field(self.field.p - self.coeffs)

    def inverse(self):
        _, inv_elt, _ = xgcd(self.coeffs, self.field.p)
        return self.field(inv_elt)

    def __eq__(self, other):
        if isinstance(other, FieldElement):
            return self.coeffs == other.coeffs and self.field == self.field
        if other == 0:
            return self.coeffs == 0
        return NotImplemented

    def __repr__(self):
        return str(self.coeffs)
