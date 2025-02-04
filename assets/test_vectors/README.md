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

All byte(s) fields are encoded as strings, hexadecimal encoding.

An example of python test is done in `../pythonref/test_vectors.py`.