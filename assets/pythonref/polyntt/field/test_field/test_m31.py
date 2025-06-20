import unittest
from time import time
from polyntt.field.m31 import M31Field, M31ExtensionField
from polyntt.field.prime_field import PrimeField


class TestM31Field(unittest.TestCase):
    def setUp(self):
        self.F = M31Field()
        self.p = (1 << 31) - 1
        self.a = self.F(123456789)
        self.b = self.F(987654321)
        self.zero = self.F(0)
        self.one = self.F(1)

    def test_addition(self):
        c = self.a + self.b
        self.assertEqual(c.coeffs, (self.a.coeffs + self.b.coeffs) % self.p)

    def test_subtraction(self):
        c = self.b - self.a
        self.assertEqual(c.coeffs, (self.b.coeffs - self.a.coeffs) % self.p)

    def test_negation(self):
        neg = -self.a
        self.assertEqual((self.a + neg).coeffs, 0)

    def test_multiplication(self):
        c = self.a * self.b
        self.assertEqual(c.coeffs, (self.a.coeffs * self.b.coeffs) % self.p)

    def test_inverse_and_division(self):
        inv = self.a.inverse()
        self.assertEqual((self.a * inv).coeffs, 1)
        div = self.b / self.a
        self.assertEqual((div * self.a).coeffs, self.b.coeffs)

    def test_equality(self):
        self.assertEqual(self.a, self.F(self.a.coeffs))
        self.assertNotEqual(self.a, self.b)


class TestM31ExtensionField(unittest.TestCase):
    def setUp(self):
        self.F = M31ExtensionField()
        self.base = self.F.base
        self.a = self.F([1, 2])  # 1 + 2·ω
        self.b = self.F([3, 4])  # 3 + 4·ω
        self.zero = self.F([0, 0])
        self.one = self.F([1, 0])
        self.omega = self.F([0, 1])

    def test_addition(self):
        c = self.a + self.b
        expected = self.F([
            self.base(1 + 3),
            self.base(2 + 4)
        ])
        self.assertEqual(c, expected)

    def test_subtraction(self):
        c = self.b - self.a
        expected = self.F([
            self.base(3 - 1),
            self.base(4 - 2)
        ])
        self.assertEqual(c, expected)

    def test_negation(self):
        neg = -self.a
        self.assertEqual(self.a + neg, self.zero)

    def test_multiplication(self):
        c = self.a * self.b
        # (1 + 2w)(3 + 4w) = (1*3 + 2*4*(-1)) + (1*4 + 2*3)w = (-5 + 10w)
        expected = self.F([
            self.base(1*3 - 2*4),
            self.base(1*4 + 2*3)
        ])
        self.assertEqual(c, expected)

    def test_inverse_and_division(self):
        inv = self.a.inverse()
        prod = self.a * inv
        self.assertEqual(prod, self.one)

        div = self.b / self.a
        self.assertEqual(div * self.a, self.b)

    def test_equality(self):
        self.assertEqual(self.a, self.F([1, 2]))
        self.assertNotEqual(self.a, self.b)


if __name__ == "__main__":
    unittest.main()
