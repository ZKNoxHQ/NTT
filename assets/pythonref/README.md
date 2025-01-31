# NTT
Generic implementation of the Number Theoretic Transform in the context of cryptography applications.

We provide tests for two examples of polynomial rings:
* Falcon `q = 12*1024+1`, working modulo `x¹⁰²⁴+1`,
* Kyber `q = 3329`, working modulo `x⁸+1`.

The implementation requires the file `ntt_constants.py`, generated using `python generate_constants.py`.

## Install
```
make install
```

## Tests
For running all tests:
```
make test
```
For running a test of `test_ntt.py` (resp. `test_poly.py`):
```
make test TEST=ntt # resp. TEST=poly
```
