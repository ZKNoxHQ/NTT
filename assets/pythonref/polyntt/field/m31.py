from polyntt.field.field import FieldElement
from polyntt.field.prime_field import PrimeField
from polyntt.field.extension_field import ExtensionField, ExtensionFieldElement


class M31Field(PrimeField):
    def __init__(self):
        self.p = (1 << 31) - 1  # 2^31 - 1

    def extension(self, degree):
        if degree == 1:
            return self
        return ExtensionField(self, degree, -1)

    def __call__(self, coeffs):

        if isinstance(coeffs, M31Element):
            return coeffs
        return M31Element(coeffs, self)

    def _reduce(self, x):
        p = self.p
        x = (x & p) + (x >> 31)
        if x >= p:
            x -= p
        return x


class M31Element(FieldElement):
    def __init__(self, coeffs, field):
        self.field = field
        self.coeffs = field._reduce(coeffs)

    def __add__(self, other):
        val = self.coeffs + other.coeffs
        if val >= self.field.p:
            val -= self.field.p
        return M31Element(val, self.field)

    def __sub__(self, other):
        val = self.coeffs - other.coeffs
        if val < 0:
            val += self.field.p
        return M31Element(val, self.field)

    def __neg__(self):
        return M31Element(0 if self.coeffs == 0 else self.field.p - self.coeffs, self.field)

    def __mul__(self, other):
        return M31Element(self.field._reduce(self.coeffs * other.coeffs), self.field)

    def inverse(self):
        # Use Fermat's little theorem since p is prime
        return M31Element(pow(self.coeffs, self.field.p - 2, self.field.p), self.field)

    def __truediv__(self, other):
        return self * other.inverse()

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __repr__(self):
        return str(self.coeffs)


class M31ExtensionElement(ExtensionFieldElement):
    def __mul__(self, other):
        # (a+bX)(c+dX) = ac - bd + X * (ad + bc)
        a, b = self.coeffs
        c, d = other.coeffs
        # Karatsuba
        ac = a*c
        bd = b*d
        t = (a+b)*(c+d)
        return self.field([(ac - bd).coeffs, (t-ac-bd).coeffs])

    def inverse(self):
        # 1/(a+bX) = (a-bX) / (a²+b²) as X² = -1
        a, b = self.coeffs
        norm_inv = self.field.base(a*a + b*b).inverse()
        return self.field([a*norm_inv, -b*norm_inv])


class M31ExtensionField(ExtensionField):
    def __init__(self):
        super().__init__(M31Field(), 2, -1)

    def __call__(self, coeffs):
        if isinstance(coeffs, M31ExtensionElement):
            return coeffs
        elif isinstance(coeffs, (tuple, list)) and len(coeffs) == 2:
            return M31ExtensionElement([self.base(coeffs[0]), self.base(coeffs[1])], self)
        else:
            return M31ExtensionElement([self.base(coeffs), self.base(0)], self)
