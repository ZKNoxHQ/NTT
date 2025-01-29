"""
This file implements tests for various parts of the Falcon.py library.

Test the code with:
> make test
"""
from ntt import intt, mul_zq, div_zq, ntt, q
from random import randint
from timeit import default_timer as timer


def test_ntt_intt(n, iterations=10):
    """Test NTT is inverse of iNTT."""
    for i in range(iterations):
        f = [randint(0, q - 1) for j in range(n)]
        f_ntt = ntt(f)
        assert intt(f_ntt) == f
    return True


def test_ntt_linearity(n, iterations=10):
    """Test the linearity of NTT."""
    for i in range(iterations):
        f = [randint(0, q - 1) for j in range(n)]
        g = [randint(0, q - 1) for j in range(n)]
        λ = randint(0, q-1)
        μ = randint(0, q-1)
        λ_f_plus_μ_g = [(λ*x+μ*y) % q for (x, y) in zip(f, g)]
        f_ntt = ntt(f)
        g_ntt = ntt(g)
        assert ntt(λ_f_plus_μ_g) == [(λ*x+μ*y) %
                                     q for (x, y) in zip(f_ntt, g_ntt)]
    return True


def test_ntt(n, iterations=10):
    """Test the NTT."""
    for i in range(iterations):
        f = [randint(0, q - 1) for j in range(n)]
        g = [randint(0, q - 1) for j in range(n)]
        h = mul_zq(f, g)
        try:
            k = div_zq(h, f)
            if k != g:
                print("(f * g) / f =", k)
                print("g =", g)
                print("mismatch")
                return False
        except ZeroDivisionError:
            continue
    return True


def wrapper_test(my_test, name, n, iterations):
    """
    Common wrapper for tests. Run the test, print whether it is successful,
    and if it is, print the running time of each execution.
    """
    d = {True: "OK    ", False: "Not OK"}
    start = timer()
    rep = my_test(n, iterations)
    end = timer()
    message = "Test {name}".format(name=name)
    message = message.ljust(20) + ": " + d[rep]
    if rep is True:
        diff = end - start
        msec = round(diff * 1000 / iterations, 3)
        message += " ({msec} msec / execution)".format(msec=msec).rjust(30)
    print(message)


def test(n, iterations=500):
    """A battery of tests."""
    # wrapper_test(test_fft, "FFT", n, iterations)
    wrapper_test(test_ntt, "NTT", n, iterations)
    wrapper_test(test_ntt_intt, "NTT(iNTT) = 1", n, iterations)
    wrapper_test(test_ntt_linearity, "NTT linearity", n, iterations)
    print("")


# Run all the tests
if (__name__ == "__main__"):

    for i in range(2, 8):  # 11):
        n = (1 << i)
        it = 100
        print("Test battery for n = {n}".format(n=n))
        test(n, it)
