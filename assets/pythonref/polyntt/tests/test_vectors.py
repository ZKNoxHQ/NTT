# This file import the json test vectors and verify the expected output.

import json
from polyntt.field.m31 import M31ExtensionField
from polyntt.field.prime_field import PrimeField
from polyntt.params import PARAMS
import unittest
from polyntt.polynomial.polynomial_ntt import PolynomialRingNTT
from polyntt.scripts.generate_test_vectors import decode


class TestVectors(unittest.TestCase):
    def test_vectors(self):
        """Run tests on the test vectors."""
        for (q, two_adicity) in PARAMS:
            F = M31ExtensionField() if q == 2**31-1 else PrimeField(q)

            # for two sizes of polynomials
            for n in [1 << (two_adicity-2), 1 << (two_adicity-1)]:
                Fx = PolynomialRingNTT(F, n)

                with open("../test_vectors/q{}_n{}.json".format(q, n), 'r') as file:
                    for test in json.load(file):
                        input = decode(
                            test['Input'], q, field='Ext' if q == 2**31-1 else 'Prime')
                        mid = len(input)//2
                        in1, in2 = input[:mid], input[mid:]

                        name = test['Name'].split('{')[1][:-1]

                        expected = decode(
                            test['Expected'],  q, field='Ext' if q == 2**31-1 else 'Prime')
                        gas = test['Gas']

                        if name == 'ntt':
                            mid = len(expected)//2
                            ex1, ex2 = expected[:mid], expected[mid:]
                            self.assertEqual(Fx.ntt(in1), ex1)
                            self.assertEqual(Fx.ntt(in2), ex2)
                        if name == 'intt':
                            mid = len(expected)//2
                            ex1, ex2 = expected[:mid], expected[mid:]
                            self.assertEqual(Fx.intt(in1), ex1)
                            self.assertEqual(Fx.intt(in2), ex2)
                        if name == 'vec_mul':
                            self.assertEqual(Fx(in1, ntt=True)
                                             * Fx(in2, ntt=True), Fx(expected))
                        if name == 'pol_mul':
                            self.assertEqual(Fx(in1) *
                                             Fx(in2), Fx(expected))
                        if name == 'vec_add':
                            self.assertEqual(Fx(in1) + Fx(in2), Fx(expected))
                        if name == 'vec_sub':
                            self.assertEqual(Fx(in1) - Fx(in2), Fx(expected))
