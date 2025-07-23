import unittest
from polyntt.field.prime_field import PrimeField
from polyntt.field.extension_field import ExtensionField, ExtensionFieldElement


class TestExtensionField(unittest.TestCase):
    def setUp(self):
        self.base = PrimeField(97)
        self.ext = ExtensionField(self.base, 2, alpha=5)
        self.a = self.ext([3, 4])  # 3 + 4·α
        self.b = self.ext([10, 20])  # 10 + 20·α
        self.zero = self.ext([0, 0])
        self.one = self.ext([1, 0])

    def test_addition(self):
        self.assertEqual(self.a + self.b, self.ext([13, 24]))

    def test_subtraction(self):
        self.assertEqual(self.b - self.a, self.ext([7, 16]))
        self.assertEqual(self.a - self.b, self.ext([-7 % 97, -16 % 97]))

    def test_negation(self):
        self.assertEqual(-self.a, self.ext([-3 % 97, -4 % 97]))

    def test_multiplication(self):
        # (3 + 4α)(10 + 20α) = 3*10 + 5*4*20 + α(3*20 + 4*10)
        real = (3*10 + 5*4*20) % 97
        imag = (3*20 + 4*10) % 97
        self.assertEqual(self.a * self.b, self.ext([real, imag]))

    def test_inverse_and_division(self):
        ainv = self.a.inverse()
        self.assertEqual(self.a * ainv, self.one)
        self.assertEqual(self.a / self.b * self.b, self.a)

    def test_equality(self):
        self.assertEqual(self.a, self.ext([3, 4]))
        self.assertNotEqual(self.a, self.b)

    def test_repr(self):
        self.assertEqual(repr(self.ext([1, 2])), "[1, 2]")

    def test_call_with_int_and_element(self):
        self.assertEqual(self.ext(7), self.ext([7, 0]))
        self.assertEqual(self.ext(self.base(9)), self.ext([9, 0]))

    def test_call_with_list_of_field_elements(self):
        self.assertEqual(
            self.ext([self.base(1), self.base(2)]),
            self.ext([1, 2]),
        )
        self.assertEqual(
            self.ext([self.base(1), 2]),
            self.ext([1, 2]),
        )


if __name__ == "__main__":
    unittest.main()
