from time import time
from polyntt.generate_test_vectors import deterministic_poly
from polyntt.poly import Poly
from polyntt.test_cases import TEST_CASES


class BenchIterativeRecursive:
    def bench_iterative_recursive():
        iterations = 100
        print("iterations:{}\n".format(iterations))

        print("q\tn\tIterative\t\tRecursive")
        for (q, two_adicity) in TEST_CASES:

            # for two sizes of polynomials
            for n in [1 << (two_adicity-2), 1 << (two_adicity-1)]:
                print("{}\t{}".format(q, n), end='\t')
                p1 = Poly(deterministic_poly(q, n), q)
                p2 = Poly(deterministic_poly(q, n), q)
                t1 = time()
                for i in range(iterations):
                    p3 = p1*p2
                t2 = time()
                p3_iterative = p3
                print("{:.2f} μs".format(
                    (t2-t1) * 10**6/iterations), end='\t\t')

                p1 = Poly(deterministic_poly(q, n), q, 'NTTRecursive')
                p2 = Poly(deterministic_poly(q, n), q, 'NTTRecursive')
                t3 = time()
                for i in range(iterations):
                    p3 = p1*p2
                t4 = time()
                p3_recursive = p3
                print("{:.2f} μs".format((t4-t3) * 10**6/iterations))

                assert p3_iterative == p3_recursive


BenchIterativeRecursive.bench_iterative_recursive()
