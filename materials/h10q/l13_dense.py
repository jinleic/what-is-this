#!/usr/bin/env python3
"""Dense unstructured small-|b| scan for the last open grid cell (89,(-2,1)).

Companion to l13_filter.py (READ-ONLY import; that file is not modified).
The structured box in l13_filter.attack only ever tested b of the form
    b = eps * f^e * q   or   b = eps * f^e / q
with f in the odd supp union and q prime.  This scan covers EVERY odd
integer b with |b| in [3, 20001) -- no shape restriction -- thereby
crossing two-prime wild sets and the tau=0 branch ('zero') that the
structured box did not cover densely.

Grid point: w = 89, unit (-2, 1)  =>  z = 89 * (-2/1) = -178.
a in {1, 17, 25, 33}: pool a-values of meta.supp_union_odd_per_a, taken
as vp-clean raw integers (all odd, so vp(a, 2) == 0, asserted at start).
tau in { sq = (1+2a^2)/(1+4a^2),  can = 2a/(1+4a^2),  one = 1, zero = 0 }
-- exact Fractions (plain int `/` would produce floats and crash
_l10_supp).

Per candidate exactly one engine call:
    st, dt = smooth_emergent(a, z, tau, b)        # y default 10**6
status 'zero' => every frozen and emergent Hilbert symbol is +1 and the
stripped remainder is exactly 1 (certified SOLUBLE candidate): print
'*** HIT' and append a full, independently re-verifiable JSON row (all
inputs recorded exactly; re-running smooth_emergent on the recorded
(a, z, tau, b) reproduces the decision bit-for-bit).
Degenerate patterns (status 'none': D==0, delta==0, _sun_h None, M0==0,
or v==0) are counted, not recorded.

Output: data/l13_dense_smallb.jsonl
    - one JSON row per 'zero' hit (written as found, flush on write)
    - data/l13_dense_cofactorbig.jsonl: one row
      {'a','tau'(name),'b','rn','rd'} per 'cofactor-big' candidate
      (stripped remainder pair from smooth_emergent detail) for the
      stage-2 smoothness ladder; row count asserted == counter
    - LAST line {'summary': {...}} with exact counters:
      tested, aligned (rows past the frozen-symbol gate:
      zero + bad + cofactor-big), zero, bad, cofactor-big, align-fail,
      none.  Identities asserted before the summary is written:
      aligned == zero + bad + cofactor-big
      tested == aligned + align-fail + none

CPU policy: single process; time.sleep(0.01) every 200 candidates.
Run from math/h10q/:   python3 l13_dense.py
"""
import json
import time
from fractions import Fraction as F

from l13_filter import smooth_emergent
from h10q import vp

A_POOL = (1, 17, 25, 33)
Z = F(89) * F(-2, 1)            # cell (89, (-2,1)); Z == -178
B_LO, B_HI = 3, 20001           # odd |b| in [3, 20001), both signs
OUT = 'data/l13_dense_smallb.jsonl'
OUT2 = 'data/l13_dense_cofactorbig.jsonl'
Y = 10 ** 6

TAUS = (
    ('sq',   lambda a: (F(1) + 2 * a * a) / (F(1) + 4 * a * a)),
    ('can',  lambda a: F(2 * a, 1 + 4 * a * a)),
    ('one',  lambda a: F(1)),
    ('zero', lambda a: F(0)),
)

SYMBOL_PASS = ('zero', 'bad', 'cofactor-big')
COUNT_KEYS = ('tested', 'aligned', 'zero', 'bad', 'cofactor-big',
              'align-fail', 'none')


def main():
    counts = dict.fromkeys(COUNT_KEYS, 0)
    T0 = time.time()
    n_cof = 0
    with open(OUT, 'w') as out, open(OUT2, 'w') as out2:
        for a_ in A_POOL:
            a = F(a_)
            assert vp(a, 2) == 0, f"a={a_} not vp-clean"
            for tname, tf in TAUS:
                tau = tf(a)
                for bmag in range(B_LO, B_HI, 2):
                    for eps in (1, -1):
                        b = F(eps * bmag)
                        counts['tested'] += 1
                        if counts['tested'] % 200 == 0:
                            time.sleep(0.01)
                        st, dt = smooth_emergent(a, Z, tau, b)
                        if st in SYMBOL_PASS:
                            counts['aligned'] += 1
                        counts[st] += 1
                        if st == 'zero':
                            em, (rn, rd), smooth = dt
                            row = {
                                'kind': 'dense-smallb-hit',
                                'cell': [89, [-2, 1]],
                                'z': str(Z.numerator),
                                'a': a_,
                                'tau_name': tname,
                                'tau': '%d/%d' % (tau.numerator, tau.denominator),
                                'b': eps * bmag,
                                'y': Y,
                                'status': st,
                                'symbols': 'all frozen + emergent Hilbert '
                                           'symbols +1 (smooth_emergent gate)',
                                'emergent_bad_primes': [str(p) for p in em],
                                'stripped_remainder_num': str(rn),
                                'stripped_remainder_den': str(rd),
                                'smooth': smooth,
                                'reverify': ("from math/h10q run: import l13_filter; "
                                             "l13_filter.smooth_emergent(F(%d), F(-178), "
                                             "F(%d, %d), F(%d)) -> ('zero', ([], (1, 1), True))"
                                             % (a_, tau.numerator, tau.denominator,
                                                eps * bmag)),
                            }
                            out.write(json.dumps(row) + '\n')
                            out.flush()
                            print(f"*** HIT a={a_} tau={tname} b={eps * bmag} "
                                  f"tau={tau} (row appended)", flush=True)
                        if st == 'cofactor-big':
                            _, (rn, rd), _ = dt
                            out2.write(json.dumps({'a': a_, 'tau': tname,
                                                   'b': str(b), 'rn': str(rn),
                                                   'rd': str(rd)}) + '\n')
                            n_cof += 1
                n_t, n_z = counts['tested'], counts['zero']
                print(f"  a={a_} tau={tname}: tested={n_t} zero={n_z} "
                      f"aligned={counts['aligned']} bad={counts['bad']} "
                      f"cofactor-big={counts['cofactor-big']} "
                      f"align-fail={counts['align-fail']} none={counts['none']} "
                      f"({time.time() - T0:.0f}s)", flush=True)
        assert counts['aligned'] == counts['zero'] + counts['bad'] + counts['cofactor-big']
        assert counts['tested'] == counts['aligned'] + counts['align-fail'] + counts['none']
        assert n_cof == counts['cofactor-big'], \
            f"cofactor rows written {n_cof} != counter {counts['cofactor-big']}"
        summary = {
            'summary': {
                'cell': [89, [-2, 1]],
                'z': '-178',
                'a_pool': list(A_POOL),
                'taus': ['sq=(1+2a^2)/(1+4a^2)', 'can=2a/(1+4a^2)',
                         'one=1', 'zero=0'],
                'b_box': 'every odd |b| in [3,20001), both signs, '
                         'lex order (a, tau, |b|, sign)',
                'candidates_expected': len(A_POOL) * len(TAUS)
                                       * 2 * ((B_HI - B_LO) // 2),
                'y': Y,
                'engine': 'l13_filter.smooth_emergent (read-only import)',
                'cofactor_rows_file': OUT2 + ' (one row per cofactor-big '
                                      'candidate; row count asserted equal '
                                      'to counts[cofactor-big])',
                'counts': counts,
                'identities': ['aligned == zero + bad + cofactor-big',
                               'tested == aligned + align-fail + none',
                               'counts are EXACT counters over the finite '
                               'box above; a 0-hit pass is evidence only '
                               'about this box'],
                'wall_seconds': round(time.time() - T0, 1),
                'generated': time.strftime('%Y-%m-%d %H:%M:%S'),
                'cpu_policy': 'single process, time.sleep(0.01) every 200 '
                              'candidates',
            },
        }
        out.write(json.dumps(summary) + '\n')
        out.flush()
    print(f"DENSE SCAN DONE: tested={counts['tested']} zero={counts['zero']} "
          f"aligned={counts['aligned']} bad={counts['bad']} "
          f"cofactor-big={counts['cofactor-big']} "
          f"cofactor_rows={n_cof} "
          f"align-fail={counts['align-fail']} none={counts['none']} "
          f"wall={time.time() - T0:.0f}s", flush=True)


if __name__ == '__main__':
    main()
