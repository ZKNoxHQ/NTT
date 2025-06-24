class PolynomialRing:
    """
    Initialise the polynomial ring:

        R = F / (X^n + 1)
    """

    def __init__(self, F, n):
        self.F = F
        self.n = n
        self.element = Polynomial

    def gen(self):
        """
        Return the generator `x` of the polynomial ring
        """
        return self([0, 1])

    def random(self):
        """
        Compute a random element of the polynomial ring with coefficients in the
        canonical range: ``[0, q-1]``
        """
        coefficients = [self.F.random() for i in range(self.n)]
        return self(coefficients)

    def __call__(self, coefficients, ntt=False):
        if isinstance(coefficients, int):
            return self.element(self, [coefficients], ntt=ntt)
        if not isinstance(coefficients, list):
            raise TypeError(
                f"Polynomials should be constructed from a list of integers, of length at most d = {self.n}"
            )
        return self.element(self, coefficients, ntt=ntt)

    def __repr__(self):
        return f"Univariate Polynomial Ring in x over Finite Field of size {self.F.p}^{self.F.degree} with modulus x^{self.n} + 1"

    def uncompact_256(self, lst, m):
        # splits the elements of lst (of 256 bits) into lists of m bits
        if not (1 <= m <= 256) or (256 % m != 0):
            raise ValueError(
                "b must be a divisor of 256 and in the range 1-256.")

        chunk_count = 256 // m
        mask = (1 << m) - 1  # Mask to extract b-bit chunks
        result = []
        for num in lst:
            chunk = []
            for i in range(chunk_count):
                chunk.append(num >> (i*m) & mask)
            result.extend(chunk)

        if self.F.degree == 2:
            return self([[result[2*i], result[2*i+1]] for i in range(len(result)//2)])
        else:
            return self(result)


class Polynomial:
    def __init__(self, parent, coeffs, ntt=False):
        """coeffs is a list of field elements or raw values; lowest degree first"""
        self.parent = parent
        while len(coeffs) != self.parent.n:
            coeffs.append(0)
        self.coeffs = [parent.F(elt) for elt in coeffs]
        # self._trim(coeffs)]
        self.ntt = ntt

    # def _trim(self, coeffs):
    #     while coeffs and coeffs[-1] == self.parent.F(0):
    #         coeffs.pop()
    #     return coeffs or [self.parent.F(0)]

    def is_zero(self):
        """
        Return if polynomial is zero: f = 0
        """
        return all(c.is_zero() for c in self.coeffs)

    def is_constant(self):
        """
        Return if polynomial is constant: f = c
        """
        return all(c == 0 for c in self.coeffs[1:])

    def degree(self):
        return len(self.coeffs) - 1

    def __repr__(self):
        return " + ".join(
            f"({c})" if i == 0 else (f"({c})·X^{i}" if i > 1 else f"({c})·X")
            for i, c in enumerate(self.coeffs)
        ) or "0"

    def __add__(self, other):
        max_len = max(len(self.coeffs), len(other.coeffs))
        result = [
            (self.coeffs[i] if i < len(self.coeffs) else self.parent.F(0)) +
            (other.coeffs[i] if i < len(other.coeffs)
             else self.parent.F(0))
            for i in range(max_len)
        ]
        return self.parent(result)

    def __neg__(self):
        """
        Returns -f, by negating all coefficients
        """
        neg_coeffs = [-x for x in self.coeffs]
        return self.parent(neg_coeffs)

    def __sub__(self, other):
        max_len = max(len(self.coeffs), len(other.coeffs))
        result = [
            (self.coeffs[i] if i < len(self.coeffs) else self.parent.F(0)) -
            (other.coeffs[i] if i < len(other.coeffs)
             else self.parent.F(0))
            for i in range(max_len)
        ]
        return self.parent(result)

    def __mul__(self, other):
        if isinstance(other, int):
            return self * self.parent(other)
        n = self.parent.n
        a = self.coeffs
        while len(a) < n:
            a += [self.parent.F.zero()]
        b = other.coeffs
        while len(b) < n:
            b += [self.parent.F.zero()]
        new_coeffs = [self.parent.F.zero() for _ in range(n)]
        for i in range(n):
            for j in range(0, n - i):
                new_coeffs[i + j] += a[i] * b[j]
        for j in range(1, n):
            for i in range(n - j, n):
                new_coeffs[i + j - n] -= a[i] * b[j]
        return self.parent(new_coeffs)

    def __call__(self, x):
        """Evaluate polynomial at field element x"""
        result = self.parent.F(0)
        power = self.parent.F(1)
        for coeff in self.coeffs:
            result += coeff * power
            power *= x
        return result

    def __eq__(self, other):
        if other == 0:
            return self.is_zero()
        else:
            return self.coeffs == other.coeffs and self.parent == other.parent

    def __divmod__(self, divisor):
        """Return (quotient, remainder) of division with remainder"""
        if divisor.degree() < 0:
            raise ZeroDivisionError()

        dividend = self.coeffs[:]
        divisor_deg = divisor.degree()
        divisor_lead = divisor.coeffs[-1]
        quotient = [self.parent.F(0)] * \
            (len(dividend) - divisor_deg + 1)

        while len(dividend) >= len(divisor.coeffs):
            coeff = dividend[-1] / divisor_lead
            deg = len(dividend) - len(divisor.coeffs)
            quotient[deg] = coeff

            # Subtract coeff * divisor * x^deg
            for i in range(len(divisor.coeffs)):
                dividend[deg + i] -= coeff * divisor.coeffs[i]
            dividend = self._trim(dividend)

        return self(self.parent, quotient, self.parent.F), self(self.parent, dividend)

    def __mod__(self, divisor):
        return divmod(self, divisor)[1]

    def __floordiv__(self, divisor):
        return divmod(self, divisor)[0]

    # def _add_(self, other):
    #     if isinstance(other, type(self)):
    #         new_coeffs = [
    #             self._add_mod_q(x, y) for x, y in zip(self.coeffs, other.coeffs)
    #         ]
    #     elif isinstance(other, int):
    #         new_coeffs = self.coeffs.copy()
    #         new_coeffs[0] = self._add_mod_q(new_coeffs[0], other)
    #     else:
    #         raise NotImplementedError(
    #             "Polynomials can only be added to each other")
    #     return new_coeffs

    # def __add__(self, other):
    #     new_coeffs = self._add_(other)
    #     return self.parent(new_coeffs)

    # def __radd__(self, other):
    #     return self.__add__(other)

    # def __iadd__(self, other):
    #     self = self + other
    #     return self

    # def _sub_(self, other):
    #     if isinstance(other, type(self)):
    #         new_coeffs = [
    #             self._sub_mod_q(x, y) for x, y in zip(self.coeffs, other.coeffs)
    #         ]
    #     elif isinstance(other, int):
    #         new_coeffs = self.coeffs.copy()
    #         new_coeffs[0] = self._sub_mod_q(new_coeffs[0], other)
    #     else:
    #         raise NotImplementedError(
    #             "Polynomials can only be subtracted from each other"
    #         )
    #     return new_coeffs

    # def __sub__(self, other):
    #     new_coeffs = self._sub_(other)
    #     return self.parent(new_coeffs)

    # def __rsub__(self, other):
    #     return -self.__sub__(other)

    # def __isub__(self, other):
    #     self = self - other
    #     return self

    # def __mul__(self, other):
    #     if isinstance(other, type(self)):
    #         new_coeffs = self._schoolbook_multiplication(other)
    #     elif isinstance(other, int):
    #         new_coeffs = [(c * other) % self.parent.q for c in self.coeffs]
    #     else:
    #         raise NotImplementedError(
    #             "Polynomials can only be multiplied by each other, or scaled by integers"
    #         )
    #     return self.parent(new_coeffs)

    # def __rmul__(self, other):
    #     return self.__mul__(other)

    # def __imul__(self, other):
    #     self = self * other
    #     return self

    # def __pow__(self, n):
    #     if not isinstance(n, int):
    #         raise TypeError(
    #             "Exponentiation of a polynomial must be done using an integer."
    #         )

    #     # Deal with negative scalar multiplication
    #     if n < 0:
    #         raise ValueError(
    #             "Negative powers are not supported for elements of a Polynomial Ring"
    #         )
    #     f = self
    #     g = self.parent(1)
    #     while n > 0:
    #         if n % 2 == 1:
    #             g = g * f
    #         f = f * f
    #         n = n // 2
    #     return g

    # def __eq__(self, other):
    #     if isinstance(other, type(self)):
    #         return self.coeffs == other.coeffs
    #     elif isinstance(other, int):
    #         if self.is_constant() and (other % self.parent.q) == self.coeffs[0]:
    #             return True
    #     return False

    # def __getitem__(self, idx):
    #     return self.coeffs[idx]

    # def __repr__(self):
    #     if self.is_zero():
    #         return "0"

    #     info = []
    #     for i, c in enumerate(self.coeffs):
    #         if c != 0:
    #             if i == 0:
    #                 info.append(f"{c}")
    #             elif i == 1:
    #                 if c == 1:
    #                     info.append("x")
    #                 else:
    #                     info.append(f"{c}*x")
    #             else:
    #                 if c == 1:
    #                     info.append(f"x^{i}")
    #                 else:
    #                     info.append(f"{c}*x^{i}")
    #     return " + ".join(info)

    def __str__(self):
        return self.__repr__()

    def compact_256(self, m):
        # compact a list of n small element of m bits into n*m/256 elements of 256 bits
        # (assuming 2^log_m = m is a divisor of n)
        if self.parent.F.degree == 2:
            a = [elt.coeffs for field_element in self.coeffs for elt in field_element]
        else:  # prime field
            a = [x.coeffs for x in self.coeffs]
        assert m < 256
        assert len(a) % m == 0
        for elt in a:
            assert elt < (1 << m)
        b = [0] * (len(a) * m // 256)
        for i in range(len(a)):
            b[(i * m) // 256] |= a[i] << ((i % (256//m)) * m)
        return b
