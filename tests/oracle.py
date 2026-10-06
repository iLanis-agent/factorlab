#!/usr/bin/env python3
"""FactorLab oracle: GNU `factor` for factorizations, pure python for the
arithmetic functions and Miller-Rabin liar bases. Writes expected.json."""
import json, math, subprocess

CORPUS = [
    '2','3','4','97','561','1105','2047','3277','8128','8191','131071',
    '65537','99991','1000003','1000000007','1000000009','1234567891',
    '3215031751','41041','825265','321197185','9746347772161',
    '2821','6601','410041','62745','987654321','1234567891011',
    '2432902008176640000','9007199254740991','9223372036854775807',
    '999999999999999999','999999999999999989','999999866000004473',
    '1000000000000037','2147483647','2147483648','67280421310721',
    '600851475143','104729','1299709','15485863'
]

def factor(n):
    out = subprocess.run(['factor', str(n)], capture_output=True, text=True, check=True).stdout
    terms = out.split(':')[1].split()
    return [int(t) for t in terms]

def is_prime(n):
    if n < 2: return False
    return len(factor(n)) == 1

def witnesses_pass(n, a):
    if n % 2 == 0: return False
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    x = pow(a % n, d, n)
    if x in (1, n - 1): return True
    for _ in range(r - 1):
        x = x * x % n
        if x == n - 1: return True
    return False

def analyze(n):
    fs = factor(n)
    grouped = {}
    for p in fs: grouped[p] = grouped.get(p, 0) + 1
    items = sorted(grouped.items())
    prime = len(fs) == 1
    phi = 1; sigma = 1; tau = 1; sqfree = True
    for p, e in items:
        phi *= (p - 1) * (p ** (e - 1) if e > 1 else 1)
        sigma *= (p ** (e + 1) - 1) // (p - 1)
        tau *= e + 1
        if e > 1: sqfree = False
    omega = len(items)
    Omega = sum(grouped.values())
    def np_(x):
        c = x + 1
        if c <= 2: return 2
        if c % 2 == 0: c += 1
        while not is_prime(c): c += 2
        return c
    def pp(x):
        if x <= 2: return None
        c = x - 1
        if c == 2: return 2
        if c % 2 == 0: c -= 1
        while c >= 3:
            if is_prime(c): return c
            c -= 2
        return None
    liars = [] if prime else [b for b in (2,3,5,7,11,13) if witnesses_pass(n, b)]
    carm = False
    if not prime and sqfree:
        carm = all((n - 1) % (p - 1) == 0 for p, _ in items)
    np1 = n + 1
    mersenne = None
    if np1 & (np1 - 1) == 0 and np1.bit_length() - 1 >= 2:
        mersenne = np1.bit_length() - 1
    classes = []
    if prime:
        classes.append('prime')
        if mersenne is not None and is_prime(mersenne): classes.append('mersenne-prime')
        if n > 2:
            if is_prime(n - 2) or is_prime(n + 2): classes.append('twin-prime')
            if is_prime(n - 4) or is_prime(n + 4): classes.append('cousin-prime')
            if is_prime(n - 6) or is_prime(n + 6): classes.append('sexy-prime')
    else:
        classes.append('composite')
        if Omega == 2: classes.append('semiprime')
        if Omega == 3 and sqfree: classes.append('sphenic')
        if carm: classes.append('carmichael')
        if sigma == 2 * n: classes.append('perfect')
        if mersenne is not None: classes.append('mersenne-composite')
        if liars: classes.append('strong-liar')
    gb = None
    if n > 2 and n % 2 == 0:
        p = 2
        while p <= n // 2:
            if is_prime(p) and is_prime(n - p): gb = [p, n - p]; break
            p = np_(p)
    divs = None
    if tau <= 5000:
        divs = [1]
        for p, e in items:
            cur = divs[:]; pk = 1
            for _ in range(e):
                pk *= p
                divs += [d * pk for d in cur]
        divs.sort()
    return dict(
        n=str(n), isPrime=prime,
        factors=[[str(p), e] for p, e in items],
        phi=str(phi), sigma=str(sigma), tau=str(tau),
        omega=omega, Omega=Omega, squarefree=sqfree,
        mobius=(0 if not sqfree else (1 if omega % 2 == 0 else -1)),
        divisors=[str(d) for d in divs] if divs else None,
        prevPrime=(str(pp(n)) if pp(n) else None), nextPrime=str(np_(n)),
        liarBases=liars, goldbach=([str(gb[0]), str(gb[1])] if gb else None),
        mersenneExp=mersenne, classes=classes,
        digits=len(str(n)), bits=n.bit_length())

out = [analyze(int(c)) for c in CORPUS]
json.dump(out, open('tests/expected.json', 'w'), indent=1)
print(len(out), 'corpus entries')
for e in out:
    print(e['n'], 'PRIME' if e['isPrime'] else 'comp', ','.join(e['classes']), 'liars:', e['liarBases'])
