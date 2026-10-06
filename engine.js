/* FactorLab engine - integer factorization workbench.
   BigInt arithmetic, deterministic Miller-Rabin (witness set valid for
   n < 3,317,044,064,679,887,385,961,981), Pollard-Brent rho.
   No dependencies. Works in browser and node. */
(function(root){
'use strict';
const B = BigInt;
const WITNESSES = [2,3,5,7,11,13,17,19,23,29,31,37].map(B);

function powmod(b, e, m){
  b %= m; let r = 1n;
  while (e > 0n){ if (e & 1n) r = r*b % m; b = b*b % m; e >>= 1n; }
  return r;
}
const SMALL_PRIMES = (function(){
  const lim = 1000, s = new Uint8Array(lim).fill(1); s[0]=s[1]=0;
  for (let i=2;i*i<lim;i++) if (s[i]) for (let j=i*i;j<lim;j+=i) s[j]=0;
  const out=[]; for (let i=2;i<lim;i++) if (s[i]) out.push(i);
  return out;
})();

function witnessPasses(n, a, d, r){
  let x = powmod(a % n, d, n);
  if (x === 1n || x === n-1n) return true;
  for (let i=1n;i<r;i++){
    x = x*x % n;
    if (x === n-1n) return true;
  }
  return false;
}

function millerDecompose(n){ // n-1 = d * 2^r
  let d = n-1n, r = 0n;
  while ((d & 1n) === 0n){ d >>= 1n; r++; }
  return [d, r];
}

function isPrime(n){
  if (n < 2n) return false;
  for (const p of SMALL_PRIMES){ const pb=B(p); if (n === pb) return true; if (n % pb === 0n) return false; }
  const [d, r] = millerDecompose(n);
  for (const a of WITNESSES){
    if (a % n === 0n) continue;
    if (!witnessPasses(n, a, d, r)) return false;
  }
  return true;
}

/* Strong-pseudoprime test against one base (for the liar display) */
function isStrongPseudoprime(n, a){
  if (n < 2n || n % 2n === 0n) return false;
  const [d, r] = millerDecompose(n);
  return witnessPasses(n, B(a), d, r);
}

function gcd(a, b){ while (b){ const t = a % b; a = b; b = t; } return a < 0n ? -a : a; }

/* Pollard-Brent rho; deterministic c sequence to avoid flaky paths */
function rho(n, c0){
  if (n % 2n === 0n) return 2n;
  const c = B(c0);
  let y = 2n, r = 1n, q = 1n, m = 128n, g = 1n, x = 0n, ys = 0n;
  const f = v => (v*v + c) % n;
  while (g === 1n){
    x = y;
    for (let i=0n;i<r;i++) y = f(y);
    let k = 0n;
    while (k < r && g === 1n){
      ys = y;
      const lim = (m < r-k) ? m : r-k;
      for (let i=0n;i<lim;i++){ y = f(y); q = q * ((x>y)?(x-y):(y-x)) % n; }
      g = gcd(q, n);
      k += m;
    }
    r *= 2n;
  }
  if (g === n){
    do { ys = f(ys); g = gcd((x>ys)?(x-ys):(ys-x), n); } while (g === 1n);
  }
  return g;
}

function factorInternal(n, out, stats){
  if (n === 1n) return;
  for (const p of SMALL_PRIMES){
    const pb = B(p);
    while (n % pb === 0n){ out.push(pb); n /= pb; stats.trialDivisions++; }
    if (n === 1n) return;
    if (pb*pb > n) break;
  }
  if (n === 1n) return;
  if (isPrime(n)){ out.push(n); return; }
  let d = n, c = 1;
  while (d === n){ d = rho(n, c++); stats.rhoCalls++; }
  factorInternal(d, out, stats);
  factorInternal(n/d, out, stats);
}

function factor(n){
  const stats = { trialDivisions: 0, rhoCalls: 0 };
  const fs = [];
  factorInternal(n, fs, stats);
  fs.sort((a,b)=> a<b?-1:a>b?1:0);
  const grouped = [];
  for (const p of fs){
    if (grouped.length && grouped[grouped.length-1][0] === p) grouped[grouped.length-1][1]++;
    else grouped.push([p, 1]);
  }
  return { factors: grouped.map(([p,e])=>[p.toString(), e]), stats };
}

function fromFactors(grouped){
  // grouped: [[BigInt prime, int exp], ...]
  let phi = 1n, sigma = 1n, tau = 1n, k = 0, squarefree = true;
  for (const [p, e] of grouped){
    phi *= (p - 1n) * (e > 1 ? p ** B(e - 1) : 1n);
    sigma *= (p ** B(e + 1) - 1n) / (p - 1n);
    tau *= B(e + 1);
    k++;
    if (e > 1) squarefree = false;
  }
  return { phi, sigma, tau, omega: k, squarefree };
}

function divisorsList(grouped, cap){
  let divs = [1n];
  for (const [p, e] of grouped){
    const cur = divs.slice();
    let pk = 1n;
    for (let i=1;i<=e;i++){
      pk *= p;
      for (const d of cur) divs.push(d * pk);
      if (divs.length > cap) return null;
    }
  }
  divs.sort((a,b)=> a<b?-1:a>b?1:0);
  return divs;
}

function nextPrime(n){
  let c = n + 1n;
  if (c <= 2n) return 2n;
  if (c % 2n === 0n) c += 1n;
  while (!isPrime(c)) c += 2n;
  return c;
}
function prevPrime(n){
  if (n <= 2n) return null;
  let c = n - 1n;
  if (c === 2n) return 2n;
  if (c % 2n === 0n) c -= 1n;
  while (c >= 3n){ if (isPrime(c)) return c; c -= 2n; }
  return null;
}

function goldbach(n){
  if (n <= 2n || n % 2n !== 0n) return null;
  for (let p = 2n; p <= n/2n; p = nextPrime(p)){
    if (isPrime(n - p)) return [p, n - p];
  }
  return null;
}

function analyze(input){
  const t0 = Date.now();
  const n = typeof input === 'bigint' ? input : B(String(input).trim());
  if (n < 2n) throw new Error('Enter an integer of 2 or more');
  if (n > 9223372036854775807n) throw new Error("Range is 2 to 2^63-1 for this lab");
  const prime = isPrime(n);
  const { factors, stats } = factor(n);
  const grouped = factors.map(([p, e]) => [B(p), e]);
  const ar = fromFactors(grouped);
  const cap = 5000;
  const divs = ar.tau <= B(cap) ? divisorsList(grouped, cap) : null;
  const liarBases = prime ? [] : [2,3,5,7,11,13].filter(b => isStrongPseudoprime(n, b));
  // Carmichael via Korselt: composite, squarefree, and (p-1)|(n-1) for all p
  let carmichael = false;
  if (!prime && ar.squarefree){
    carmichael = grouped.every(([p]) => (n - 1n) % (p - 1n) === 0n);
  }
  const Omega = grouped.reduce((s, [,e]) => s + e, 0);
  const mersenneExp = (function(){
    const np = n + 1n;
    if ((np & (np - 1n)) !== 0n) return null;
    const e = np.toString(2).length - 1;
    return e >= 2 ? e : null;
  })();
  const classes = [];
  if (prime){
    classes.push('prime');
    if (mersenneExp !== null && isPrime(B(mersenneExp))) classes.push('mersenne-prime');
    if (n > 2n){
      if (isPrime(n - 2n)) classes.push('twin-prime');
      if (isPrime(n + 2n)) classes.push('twin-prime');
      if (isPrime(n - 4n)) classes.push('cousin-prime');
      if (isPrime(n + 4n)) classes.push('cousin-prime');
      if (isPrime(n - 6n)) classes.push('sexy-prime');
      if (isPrime(n + 6n)) classes.push('sexy-prime');
    }
  } else {
    classes.push('composite');
    if (Omega === 2) classes.push('semiprime');
    if (Omega === 3 && ar.squarefree) classes.push('sphenic');
    if (carmichael) classes.push('carmichael');
    if (ar.sigma === 2n * n) classes.push('perfect');
    if (mersenneExp !== null) classes.push('mersenne-composite');
    if (liarBases.length) classes.push('strong-liar');
  }
  const gb = goldbach(n);
  return {
    n: n.toString(),
    isPrime: prime,
    factors: factors,
    phi: ar.phi.toString(),
    sigma: ar.sigma.toString(),
    tau: ar.tau.toString(),
    omega: ar.omega,
    Omega: Omega,
    squarefree: ar.squarefree,
    mobius: ar.squarefree ? (ar.omega % 2 === 0 ? 1 : -1) : 0,
    divisors: divs ? divs.map(d => d.toString()) : null,
    divisorCountShown: divs ? divs.length : null,
    prevPrime: prevPrime(n) ? prevPrime(n).toString() : null,
    nextPrime: nextPrime(n).toString(),
    liarBases: liarBases,
    goldbach: gb ? [gb[0].toString(), gb[1].toString()] : null,
    mersenneExp: mersenneExp,
    classes: classes,
    digits: n.toString().length,
    bits: n.toString(2).length,
    ms: Date.now() - t0,
    stats: stats
  };
}

const api = { analyze, isPrime: n => isPrime(typeof n === 'bigint' ? n : B(String(n))), factor: n => factor(typeof n === 'bigint' ? n : B(String(n))) };
if (typeof module !== 'undefined' && module.exports) module.exports = api;
root.FactorLab = api;
})(typeof self !== 'undefined' ? self : globalThis);
