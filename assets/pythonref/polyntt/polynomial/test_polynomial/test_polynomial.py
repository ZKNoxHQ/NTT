import unittest
from polyntt.polynomial.polynomial import PolynomialRing
from polyntt.field.m31 import M31ExtensionField


class TestPolynomialOverM31_2(unittest.TestCase):
    def setUp(self):
        # Initializes F_{p^2} where p = 2³¹ - 1
        self.field = M31ExtensionField()
        self.zero = self.field([0, 0])
        self.one = self.field([1, 0])
        self.omega = self.field([0, 1])
        self.Fpx = PolynomialRing(self.field, 4)

    def test_repr(self):
        poly = self.Fpx([self.zero, self.one, self.omega])
        s = repr(poly)
        self.assertIn("X^2", s)
        self.assertIn("X", s)

    def test_addition(self):
        f = self.Fpx([self.one, self.omega])
        g = self.Fpx([self.omega, self.one])
        h = f + g
        self.assertEqual(h.coeffs[0], self.one + self.omega)
        self.assertEqual(h.coeffs[1], self.omega + self.one)

    def test_subtraction(self):
        f = self.Fpx([self.one, self.omega])
        g = self.Fpx([self.omega, self.one])
        h = f - g
        self.assertEqual(h.coeffs[0], self.one - self.omega)
        self.assertEqual(h.coeffs[1], self.omega - self.one)

    def test_multiplication(self):
        f = self.Fpx([self.one, self.omega])
        g = self.Fpx([self.one, self.one])
        h = f*g
        # (1 + ωX)(1 + X) = 1 + (ω+1)X + ωX^2
        expected = [
            self.one,
            self.omega + self.one,
            self.omega
        ]
        self.assertEqual(h.coeffs, expected)

    def test_multiplication_mod(self):
        # Fp[x]/(x⁴+1)
        f = self.Fpx([self.one, self.zero, self.zero, self.omega])
        g = self.Fpx([self.one, self.zero, self.one, self.zero])
        h = f*g
        # (1 + ωX³)(1 + X²) = 1 + X² + ωX³ + ωX⁵ = 1 + X² + ωX³ + ω(-X) = 1 - ωX + Χ² + ωX³
        expected = [
            self.one,
            -self.omega,
            self.one,
            self.omega
        ]
        self.assertEqual(h.coeffs, expected)

    def test_evaluation(self):
        poly = self.Fpx([self.one, self.omega, self.one])
        x = self.omega
        val = poly(x)
        # Should be 1 + ω·x + x^2
        expected = self.one + self.omega * x + x * x
        self.assertEqual(val, expected)

    # def test_divmod(self):
    #     f = self.Fpx([self.one, self.zero, self.omega])  # ω·X^2 + 1
    #     g = self.Fpx([self.one, self.one])  # X + 1
    #     q, r = divmod(f, g)
    #     self.assertTrue(isinstance(q, Polynomial))
    #     self.assertTrue(isinstance(r, Polynomial))
    #     self.assertEqual(f, q * g + r)

    def test_debug(self):
        # P = self.Fpx([1, 2, 3])
        Q = self.Fpx([self.field([1, 1]), 2, 3])
        # self.assertEqual(P, Q)


if __name__ == "__main__":
    unittest.main()
