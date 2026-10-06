/* FactorLab node runner: engine vs GNU-factor-derived expected.json */
'use strict';
const fs = require('fs');
const path = require('path');
const F = require(path.join(__dirname, '..', 'engine.js'));
const expected = JSON.parse(fs.readFileSync(path.join(__dirname, 'expected.json'), 'utf8'));

let checks = 0, fails = [];
function chk(c, m){ checks++; if (!c) fails.push(m); }

for (const e of expected){
  const a = F.analyze(e.n);
  chk(a.isPrime === e.isPrime, `${e.n}: isPrime ${a.isPrime} != ${e.isPrime}`);
  chk(JSON.stringify(a.factors) === JSON.stringify(e.factors), `${e.n}: factors ${JSON.stringify(a.factors)} != ${JSON.stringify(e.factors)}`);
  chk(a.phi === e.phi, `${e.n}: phi ${a.phi} != ${e.phi}`);
  chk(a.sigma === e.sigma, `${e.n}: sigma ${a.sigma} != ${e.sigma}`);
  chk(a.tau === e.tau, `${e.n}: tau ${a.tau} != ${e.tau}`);
  chk(a.omega === e.omega, `${e.n}: omega ${a.omega} != ${e.omega}`);
  chk(a.Omega === e.Omega, `${e.n}: Omega ${a.Omega} != ${e.Omega}`);
  chk(a.squarefree === e.squarefree, `${e.n}: squarefree mismatch`);
  chk(a.mobius === e.mobius, `${e.n}: mobius ${a.mobius} != ${e.mobius}`);
  chk(a.prevPrime === e.prevPrime, `${e.n}: prevPrime ${a.prevPrime} != ${e.prevPrime}`);
  chk(a.nextPrime === e.nextPrime, `${e.n}: nextPrime ${a.nextPrime} != ${e.nextPrime}`);
  chk(JSON.stringify(a.liarBases) === JSON.stringify(e.liarBases), `${e.n}: liars ${JSON.stringify(a.liarBases)} != ${JSON.stringify(e.liarBases)}`);
  chk(a.mersenneExp === e.mersenneExp, `${e.n}: mersenneExp ${a.mersenneExp} != ${e.mersenneExp}`);
  chk(JSON.stringify(a.classes) === JSON.stringify(e.classes), `${e.n}: classes ${JSON.stringify(a.classes)} != ${JSON.stringify(e.classes)}`);
  chk(a.digits === e.digits, `${e.n}: digits mismatch`);
  chk(a.bits === e.bits, `${e.n}: bits mismatch`);
  if (e.divisors){
    chk(JSON.stringify(a.divisors) === JSON.stringify(e.divisors), `${e.n}: divisors mismatch (len ${a.divisors && a.divisors.length} vs ${e.divisors.length})`);
  } else {
    chk(a.divisors === null, `${e.n}: expected capped divisors`);
  }
  if (e.goldbach){
    chk(!!a.goldbach, `${e.n}: missing goldbach`);
    if (a.goldbach){
      const [p, q] = a.goldbach.map(BigInt);
      chk(p + q === BigInt(e.n), `${e.n}: goldbach sum wrong`);
      chk(F.isPrime(p) && F.isPrime(q), `${e.n}: goldbach terms not prime`);
    }
  } else {
    chk(a.goldbach === null, `${e.n}: unexpected goldbach`);
  }
}
console.log(`${checks} checks, ${fails.length} failures`);
if (fails.length){ fails.forEach(f => console.log('FAIL', f)); process.exit(1); }
console.log('ALL PASS');
