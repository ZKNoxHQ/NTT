from polyntt.ntt_mert import NTTMert
from time import time
from polyntt.generate_test_vectors import deterministic_poly
from polyntt.poly import Poly
from polyntt.test_cases import TEST_CASES

iterations = 100
print("iterations:{}\n".format(iterations))

print("        q\tn\tPatrick\t\tSchoolbook\tMert")
for (q, two_adicity) in TEST_CASES:

    for n in [1 << (two_adicity-2), 1 << (two_adicity-1)]:  # for two sizes of polynomials
        print("{:10}\t{}".format(q, n), end='\t')
        p1 = Poly(deterministic_poly(q, n), q)
        p2 = Poly(deterministic_poly(q, n), q)
        t1 = time()
        for i in range(iterations):
            p3 = p1*p2
        t2 = time()
        p3_patrick = p3
        print("{:.2f} μs".format((t2-t1) * 10**6/iterations), end='\t')

        p1 = Poly(deterministic_poly(q, n), q)
        p2 = Poly(deterministic_poly(q, n), q)
        t3 = time()
        for i in range(iterations):
            p3 = p1.mul_schoolbook(p2)
        t4 = time()
        p3_schoolbook = p3
        print("{:.2f} μs".format((t4-t3) * 10**6/iterations), end='\t')

        Mert = NTTMert(q)
        p1 = Poly(deterministic_poly(q, n), q)
        p2 = Poly(deterministic_poly(q, n), q)
        t5 = time()
        for i in range(iterations):
            p3 = Mert.poly_mul_mert(p1.coeffs, p2.coeffs)
        t6 = time()
        p3_mert = p3
        print("{:.2f} μs".format((t6-t5) * 10**6/iterations))

        assert p3_mert == p3_patrick.coeffs and p3_mert == p3_schoolbook.coeffs
