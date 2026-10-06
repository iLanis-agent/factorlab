# FactorLab - integer factorization workbench

Type any integer up to 2^63-1 and see what it is made of. A BigInt engine runs
deterministic Miller-Rabin primality testing, trial division, and Pollard-Brent
rho entirely in your browser - no libraries, no network calls.

**Live app:** https://ilanis-agent.github.io/factorlab/app.html

## What it computes

- Primality: Miller-Rabin with the 12-witness set `[2..37]`, deterministic for
  n < 3,317,044,064,679,887,385,961,981 (Jaeschke / standard result)
- Full factorization via trial division + Pollard-Brent rho
- Arithmetic functions: phi, sigma, tau, omega, Omega, Mobius, divisor list
- Factor tree, prime neighbors, Goldbach pair (even n)
- Classifications: twin/cousin/sexy primes, Mersenne primes and composites,
  semiprime, sphenic, Carmichael (via Korselt's criterion), perfect numbers,
  and strong liars - composites that fool single-base Miller-Rabin tests

## Tests

`tests/oracle.py` factors a 42-number corpus with **GNU factor** and derives
everything else in pure Python (totients, divisors, liar bases via its own
Miller-Rabin, classifications). `tests/run_tests.js` runs the JS engine over
the same corpus - **764 checks**.

Corpus includes Carmichael numbers (561, 41041, 9746347772161...), strong
pseudoprimes (2047, 3215031751 - liar for bases 2,3,5,7), Mersenne primes and
composites, a perfect number, 20!, 2^63-1, and an 18-digit semiprime.

    python3 tests/oracle.py   # needs GNU factor (coreutils)
    node tests/run_tests.js

## Files

- `engine.js` - BigInt Miller-Rabin + Pollard-Brent rho + arithmetic functions
- `app.html` - the lab UI
- `tests/` - corpus, oracle, runner
