# Test vectors

## Generation
In order to generate the test vectors:
```
cd ../pythonref/
python generate test_vectors.py
```
This will create `.json` files.

## Format of the test vectors

The test data is declared in a 'json' file:

```
[
    {
        "Input": bytes of one (or two) polynomials interpreted as byte concatenation,
        "Name": the name of the test,
        "Expected": bytes of one (or two) expected polynomials,
        "Gas": the cost of the gas
    },
    ...
]
```

All byte(s) fields are encoded as strings, in hexadecimal encoding.

### Example
For Dilithium-256, `q = 8380417` and polynomials have degree `255` (`n=256`).
The first test [here](q8380417.json) checks the `ntt` of two polynomials `f` and `g` separately.
We can get `f` and `g` from the byte as follows:
* Input of the first test in hexadecimal:
```
2e000a890124029d2cf50ad3023d1a3d1dfd18..a9161029d413f22a7c1f7a11de2cfc19b1181d
```
* Polynomials $f$ and $g$ in hexadecimal:
```
2e000a890124029d2cf50ad3023d1a3d1dfd18..aa119c071322e42dd825a20b7e0abf01b62810
1ffd081e19680e372a70108f001322b904e30a..a9161029d413f22a7c1f7a11de2cfc19b1181d
```
* Coefficients of $f$ and $g$ modulo $q$:
```
[3014666, 8978724, 171308, 16059091, ... , 2466315, 8260287,  112168,  16]
[2096392, 1972584, 931626, 7344271,  ... , 2062865, 14560508, 1683736, 29]
```
* Polynomials $f(X)$ and $g(X)$:
$$f(x) = 3014666+ 8978724X + ... + 112168 X^{254} + 16X^{255}, \\g(x) = 2096392+1972584X+ ... + 1683736X^{254} + 29X^{255}.$$

## Tests
An example of python test is done in `../pythonref/test_vectors.py`.
