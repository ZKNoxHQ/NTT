from polyntt.field.field import Field, FieldElement

'''
For now, we support only quadratic extension defined with a polynomial of the form x²-α
'''


class ExtensionField(Field):
    def __init__(self, base_field, degree, alpha, check_alpha=False):
        assert degree == 2, "Only quadratic supported for now"
        self.p = base_field.p
        self.base = base_field
        self.degree = degree
        if check_alpha:
            assert pow(alpha, (self.base.p-1)//2, self.base.p) != 1
        self.alpha = alpha  # integer

    def __call__(self, coeffs):
        list_coeffs = []
        if isinstance(coeffs, ExtensionFieldElement):
            return coeffs
        if not (isinstance(coeffs, list)):
            list_coeffs = [coeffs]
        else:
            list_coeffs = coeffs
        while len(list_coeffs) != self.degree:
            list_coeffs += [0]
        return ExtensionFieldElement(list_coeffs, self)

    def random(self):
        return self([self.base.random() for i in range(self.degree)])

    def zero(self):
        return self([self.base.zero() for i in range(self.degree)])


class ExtensionFieldElement(FieldElement):
    def __init__(self, coeffs, field):
        super().__init__(field)
        if isinstance(coeffs, type(field)):
            return coeffs
        if isinstance(coeffs, int):
            coeffs = [coeffs]
        while len(coeffs) != field.degree:
            coeffs.append(0)
        self.coeffs = [field.base(c) for c in coeffs]

    def __add__(self, other):
        # a+bX + c+dX = a+c + (b+c)X
        a, b = self.coeffs
        c, d = other.coeffs
        return self.field([(a+c).coeffs, (b+d).coeffs])

    def __neg__(self):
        a, b = self.coeffs
        return self.field([-a, -b])

    def __sub__(self, other):
        a, b = self.coeffs
        c, d = other.coeffs
        return self.field([a-c, b-d])

    def __mul__(self, other):
        # (a+bX)(c+dX) = ac + alpha * bd + X * (ad + bc)
        a, b = self.coeffs
        c, d = other.coeffs
        # Karatsuba
        ac = a*c
        bd = b*d
        t = (a+b)*(c+d)
        return self.field([(ac + bd.mul_by_small(self.field.alpha)).coeffs, (t-ac-bd).coeffs])

    def inverse(self):
        # 1/(a+bX) = (a-bX) / (a²-αb²) as X² = α
        a, b = self.coeffs
        norm_inv = self.field.base(
            a*a - (b*b).mul_by_small(self.field.alpha)).inverse()
        return self.field([a*norm_inv, -b*norm_inv])

    def __truediv__(self, other):
        return self * other.inverse()

    def __eq__(self, other):
        if isinstance(other, FieldElement):
            return self.coeffs == other.coeffs and self.field == self.field
        if other == 0:
            return self.coeffs == 0
        return NotImplemented

    def __repr__(self):
        # return f"{self.coeffs[0]} + {self.coeffs[1]}·α"
        return f"[{self.coeffs[0]}, {self.coeffs[1]}]"

    def is_zero(self):
        return all([x.is_zero() for x in self.coeffs])
