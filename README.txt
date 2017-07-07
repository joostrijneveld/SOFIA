# SOFIA: MQ-based signatures in the QROM

This code package contains the source code that accompanies the paper "SOFIA:
MQ-based signatures in the QROM". It contains a C reference implementation and
an optimized avx2 implementation, in the respective directories. The code is
self-contained and does not require any external dependencies.

In order to verify correctness of the code, several tests are supplied. These
can be called by calling `make test` in one of the code directories. Simply
calling `make` builds the tests, but does not run them yet. This also runs
test/speed, which provides benchmarks for the API functions (i.e. crypto_sign,
crypto_sign_open and crypto_keypair), as well as numerous internal functions.

To verify compatibility between the reference code and the avx2 code, we have
included a test_compatibility.sh script, which iterates over all combinations of
the API functions from both implementations to verify correct results. Simply
call it by running `./test_compatibility.sh` from the top directory.

The scripts directory contains two scripts that were used when constructing this
implementation. The file `arrange_monomials.py` was used to create and test an
arrangement of rotations such that during an evaluation of the MQ function, all
quadratic monomials are produced. Using `signature_size.py`, we have computed
the resulting signature size for various parameter sets.
