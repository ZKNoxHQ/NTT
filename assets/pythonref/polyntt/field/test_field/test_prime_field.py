import unittest
from polyntt.field.m31 import M31Field
from polyntt.field.prime_field import PrimeField
from polyntt.field.extension_field import ExtensionField
from time import time


class TestPrimeField(unittest.TestCase):
    def setUp(self):
        self.p = 97  # a small prime
        self.F = PrimeField(self.p)
        self.a = self.F(3)
        self.b = self.F(10)
        self.zero = self.F(0)
        self.one = self.F(1)

    def test_addition(self):
        self.assertEqual(self.a + self.b, self.F(13))
        self.assertEqual(self.F(90) + self.F(10), self.F(3))

    def test_subtraction(self):
        self.assertEqual(self.b - self.a, self.F(7))
        self.assertEqual(self.a - self.b, self.F(90))  # 3 - 10 mod 97

    def test_multiplication(self):
        self.assertEqual(self.a * self.b, self.F(30))
        self.assertEqual(self.a * self.F(0), self.zero)

    def test_negation(self):
        self.assertEqual(-self.a, self.F(94))
        self.assertEqual(-self.F(0), self.zero)

    def test_inverse(self):
        a_inv = self.a.inverse()
        self.assertEqual(self.a * a_inv, self.one)

    def test_division(self):
        div = self.b / self.a
        self.assertEqual(div * self.a, self.b)

    def test_equality(self):
        self.assertEqual(self.a, self.F(3))
        self.assertNotEqual(self.a, self.b)

    def test_repr(self):
        self.assertEqual(repr(self.a), "3")

    def test_coeffs_mod_reduction(self):
        self.assertEqual(self.F(100), self.F(3))

    def test_extension(self):
        ext = self.F.extension(2, 5)
        self.assertIsInstance(ext, ExtensionField)
        self.assertEqual(ext.base.p, self.p)

    def test_is_prime_true(self):
        primes = [2, 3, 5, 7, 11, 97, 127]
        for p in primes:
            self.assertTrue(PrimeField(p)._is_prime(p))

    def test_is_prime_false(self):
        composites = [0, 1, 4, 6, 9, 100]
        for n in composites:
            with self.assertRaises(AssertionError):
                PrimeField(n)

    def test_eq_with_int(self):
        one = self.F(1)
        self.assertEqual(one, 1)

    def test_bench(self):
        nreps = 10

        F1 = PrimeField(p=2**31-2**24+1)
        a = F1.random()
        b = F1.random()

        t = time()
        for i in range(nreps):
            c = a._reduce()
        print(
            "BabyBear prime field multiplication:: {:.0f}ns".format(10**9*(time() - t)/nreps))

        F2 = M31Field()
        a = F2.random()
        b = F2.random()

        t = time()
        for i in range(nreps):
            c = a._reduce()
        print(
            "M31 prime field multiplication:: {:.0f}ns".format(10**9*(time() - t)/nreps))


if __name__ == "__main__":
    unittest.main()
