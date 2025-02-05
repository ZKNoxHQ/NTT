from random import randint
from time import time
from generate_test_vectors import deterministic_poly
from poly import Poly
from test_cases import TEST_CASES

iterations = 100
print("iterations:{}\n".format(iterations))

print("q\tn\tIterative\t\tRecursive")
for (q, two_adicity) in [TEST_CASES[1]]:

    for n in [1 << (two_adicity-2), 1 << (two_adicity-1)]:  # for two sizes of polynomials
        print("{}\t{}".format(q, n), end='\t')
        p1 = Poly(deterministic_poly(q, n), q)
        p2 = Poly(deterministic_poly(q, n), q)
        t1 = time()
        for i in range(iterations):
            p3 = p1*p2
        t2 = time()
        p3_iterative = p3
        print("{:.2f} ms".format((t2-t1) * 10**3), end='\t\t')

        p1 = Poly(deterministic_poly(q, n), q, 'NTTRecursive')
        p2 = Poly(deterministic_poly(q, n), q, 'NTTRecursive')
        t3 = time()
        for i in range(iterations):
            p3 = p1*p2
        t4 = time()
        p3_recursive = p3
        print("{:.2f} ms".format((t4-t3) * 10**3))

        assert p3_iterative == p3_recursive
