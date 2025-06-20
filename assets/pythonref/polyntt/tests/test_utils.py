# -*- coding: utf-8 -*-
from random import randint
import unittest
from polyntt.field.prime_field import PrimeField


class TestUtils(unittest.TestCase):
    def shortDescription(self):
        return None  # This prevents unittest from printing docstrings

    def test_batched_modular_multiplication(self, iterations=100):
        """Test the batched modular multiplication."""
        q = 3329
        F = PrimeField(q)
        n = 256
        for i in range(iterations):
            L = [F(randint(1, q-1)) for i in range(n)]  # non-zeros!
            M = F.batch_inversion(L)
            for (l, m) in zip(L, M):
                self.assertEqual(l*m, F(1))
