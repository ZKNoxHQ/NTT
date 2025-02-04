# This file import the json test vectors and verify the expected output.

import json
from ntt import NTT
from poly import Poly
from test_cases import TEST_CASES


def decode(hex_poly, q):
    # By default, q is a 32-bit integer.
    # NB: this does not apply to q_baby_bear nor q_plonky2.
    size_q = (q.bit_length()+7)//8
    bytes_poly = bytes.fromhex(hex_poly)
    return [int.from_bytes(bytes_poly[i:i+size_q], 'big') for i in range(0, len(bytes_poly), size_q)]


for (q, two_adicity) in TEST_CASES:

    for n in [1 << (two_adicity-2), 1 << (two_adicity-1)]:  # for two sizes of polynomials

        with open("../test_vectors/q{}_n{}.json".format(q, n), 'r') as file:
            T = NTT(q)
            for test in json.load(file):
                input = decode(test['Input'], q)
                mid = len(input)//2
                in1, in2 = input[:mid], input[mid:]

                name = test['Name'].split('{')[1][:-1]

                expected = decode(test['Expected'],  q)
                gas = test['Gas']

                if name == 'ntt':
                    mid = len(expected)//2
                    ex1, ex2 = expected[:mid], expected[mid:]
                    assert T.ntt(in1) == ex1 and T.ntt(in2) == ex2
                if name == 'intt':
                    mid = len(expected)//2
                    ex1, ex2 = expected[:mid], expected[mid:]
                    assert T.intt(in1) == ex1 and T.intt(in2) == ex2
                if name == 'vec_mul':
                    assert T.mul_ntt(in1, in2) == expected
                if name == 'pol_mul':
                    assert Poly(in1, q) * \
                        Poly(in2, q) == Poly(expected, q)
                if name == 'vec_add':
                    assert T.add_ntt(in1, in2) == expected
                if name == 'vec_sub':
                    assert T.sub_ntt(in1, in2) == expected
