"""L16 emergent-symbol structure probe (lead-run).

For sampled classes and prime-member k in 0..K: compute the FULL emergent
prime list (every unfrozen odd-valuation prime of P(b) with its Hilbert
symbol (x0,d0)_p), including the +1 primes smooth_emergent drops.

Per k row: {class, k, Q, em_full: [[p, sym]], verdict}. Then analyse:
 - distribution of |E_full| per class (rung predictor)
 - for fixed small emergent p: sequence s_p(k) defined where p|P(b(k)):
   is (d0(k)/p) periodic in k mod p? (it must be: b mod p decides it)
 - correlation structure: P(zero) vs product law 2^{-|E|+1} (parity-corrected)
"""
import json
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, '/Users/jinleic/jinleic-workspace/math/h10q')
from h10q import (_L11_CLASSES, _L12_ESCAPE, _is_prime, _sun_h, _l10_P,
                  _l10_supp, hilbert)
from l13_filter import pattern, _small_primes_upto

D = json.load(open('/Users/jinleic/jinleic-workspace/math/h10q/data/l13h_all_closures.json'))
recs = {(r['cell'][0], tuple(r['cell'][1])): r for r in D['records']}

# sample: two outliers, three typical k0 classes, one factorint multi-E
SAMPLES = [
    ((79, (5, 1)), 600),      # outlier k=1694, prime rung, a=1 w=3(4)
    ((83, (-5, 1)), 600),     # outlier k=1511, prime rung, a=1 w=3(4)
    ((3, (3, 1)), 600),       # smallest emergent prime (19 dig) observed
    ((79, (5, 1)), 0),        # placeholder dedup guard
    ((67, (-1, 1)), 600),     # ESC, alt, k=0 closure
    ((59, (7, 3)), 600),      # ESC alt k=51
    ((7, (-3, 1)), 600),      # ESC factorint outlier k=184
    ((11, (7, 3)), 600),      # ESC factorint k=183
]
SAMPLES = [s for i, s in enumerate(SAMPLES) if s[0] not in {x[0] for x in SAMPLES[:i]}]

SMALL = _small_primes_upto(10**6 + 1)
out = open('/Users/jinleic/jinleic-workspace/math/h10q/data/l16_char.jsonl', 'w')
out.write(json.dumps({'meta': {'probe': 'L16 emergent-symbol structure',
                               'date': '2026-08-18',
                               'note': 'full emergent list (signs +1/-1) per prime-member k'}}) + '\n')

def full_emergent(a, z, tau, b):
    A, delta, alpha, Z, D0, s = pattern(a, z, tau)
    c0 = _sun_h(a, b, Z)
    M0 = 16 - delta*c0*c0 - 32*A*b*s*s
    P = _l10_P(a, Z, D0, A, delta, s)
    v = sum(c*b**i for i, c in enumerate(P))
    x0, d0 = alpha*M0, alpha*2*b
    frozen = {2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta) | _l10_supp(b)
    rn, rd = abs(v.numerator), v.denominator
    for p in sorted(frozen):
        while rn % p == 0:
            rn //= p
        while rd % p == 0:
            rd //= p
    emf = []                       # full odd-valuation unfrozen primes, signed
    for p in SMALL:
        if p in frozen:
            continue
        e = e2 = 0
        while rn % p == 0:
            rn //= p; e += 1
        while rd % p == 0:
            rd //= p; e2 += 1
        if (e - e2) % 2:
            emf.append([int(p), int(hilbert(x0, d0, p))])
    return emf, (rn, rd)

tG = time.time()
for cell, K in SAMPLES:
    r = recs[cell]
    a = F(r['a']); A = 1 + 4*a*a; tau = (1 + 2*a*a)/A
    z = F(cell[0])*F(*cell[1])
    a_, delta, alpha, Z, D0, s = pattern(a, z, tau)
    rows, n_k0 = [], 0
    per_p = {}                       # p -> list of (k, sym)
    for k in range(K + 1):
        Q = r['q1'] + k*int(r['N'])
        if not _is_prime(Q):
            continue
        b = F(r['eps']*r['f']*Q)
        emf, cof = full_emergent(a, z, tau, b)
        neg = [p for p, sy in emf if sy == -1]
        verdict = 'zero' if not neg and cof == (1, 1) else (
                  'cofactor-big' if not neg else 'bad')
        rows.append({'cell': list(cell), 'k': k, 'Q': str(Q),
                     'em_full': emf, 'cofactor_digits': [len(str(cof[0])), len(str(cof[1]))],
                     'verdict': verdict})
        for p, sy in emf:
            if p < 10**6:
                per_p.setdefault(p, []).append((k, sy))
        out.write(json.dumps(rows[-1]) + '\n')
    # periodicity check per recurring small p: sym determined by k mod p?
    per_p_summary = {}
    for p, lst in per_p.items():
        if len(lst) >= 4:
            tab = {}
            ok = True
            for k, sy in lst:
                res = k % p
                if res in tab and tab[res] != sy:
                    ok = False
                    break
                tab[res] = sy
            per_p_summary[str(p)] = {'n': len(lst),
                                     'neg': sum(1 for _, sy in lst if sy == -1),
                                     'residue_determined_mod_p': ok}
    n_zero = sum(1 for rw in rows if rw['verdict'] == 'zero')
    n_bad = sum(1 for rw in rows if rw['verdict'] == 'bad')
    lens = [len(rw['em_full']) for rw in rows]
    big = sum(1 for rw in rows if rw['cofactor_digits'][0] > 2)
    out.write(json.dumps({'type': 'class-summary', 'cell': list(cell),
                          'family': r['family'], 'k0_class_member': r['k_zero'],
                          'prime_members_tested': len(rows),
                          'k_zero_found_here': [row['k'] for row in rows if row['verdict'] == 'zero'][:4],
                          'n_zero': n_zero, 'n_bad': n_bad,
                          'em_len_hist': {str(x): lens.count(x) for x in sorted(set(lens))},
                          'recurring_small_emergent': per_p_summary,
                          'elapsed_s': round(time.time()-tG, 1)}) + '\n')
    print(f"{cell}: {len(rows)} prime members ({time.time()-tG:.0f}s) "
          f"zero={n_zero} bad={n_bad} lens={sorted(set(lens))}", flush=True)
out.write(json.dumps({'summary': {'wall_s': round(time.time()-tG, 1)}}) + '\n')
out.close()
print('L16 PROBE COMPLETE', flush=True)
