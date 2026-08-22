"""Lead replay of ALL 293 class closures (data/l13h_all_closures.json).
Per record: rebuild member from authority class data or record fields;
assert ladder verdict 'zero'; recompute ramified(x0,d0)==[] with raised
brent budget; output data/l13h_replay.jsonl one verdict per line.
"""
import json
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, '/Users/jinleic/jinleic-workspace/math/h10q')
import h10q
from h10q import (_L11_CLASSES, _L12_ESCAPE, _is_prime, _sun_h, ramified,
                  FactorBudget, PrimalityBound)
import l12_class
from l13_filter import smooth_emergent, cofactor_decide
from l13h_scan import excluded_primes_L11

_orig = h10q._brent
h10q._brent = lambda n, budget=4000000: _orig(n, budget=budget)
import l13_filter                       # shares h10q module object
import math as _math
from h10q import hilbert as _hil

D = json.load(open('/Users/jinleic/jinleic-workspace/math/h10q/data/l13h_all_closures.json'))
out = open('/Users/jinleic/jinleic-workspace/math/h10q/data/l13h_replay.jsonl', 'w')
n_ok = ram_ok = ram_ref = 0
t0g = time.time()
for r in D['records']:
    w, ut = r['cell'][0], tuple(r['cell'][1])
    key = (w, ut)
    z = F(w)*F(*ut)
    fam = r['family']
    # authority class data for the member base (N etc. match record for
    # waves 1/2; wave-3 alt/deepk records carry their own a/eps/f/q1/N)
    a = F(r['a'])
    eps = r['eps']; f = r['f']; q1 = r['q1']; N = int(r['N'])
    A = 1 + 4*a*a
    tau = (1 + 2*a*a)/A
    # authority cross-check on non-alt records
    if 'alt' not in r['source'] and 'deepk' not in r['source']:
        if fam == 'L11':
            row = _L11_CLASSES[key]
            assert F(*row[0]) == a and row[1] == eps and row[2] == q1 \
                and row[3] == N, (key, 'authority mismatch')
        else:
            f0, eps0, q0 = _L12_ESCAPE[key]
            assert f0 == f and eps0 == eps and q0 == q1, (key, 'esc mismatch')
    Q = q1 + r['k_zero']*N
    assert str(Q) == r['Q'], (key, 'Q mismatch', Q, r['Q'])
    assert _is_prime(Q), (key, 'Q not proved prime')
    b = F(eps*f*Q)
    st, dt = smooth_emergent(a, z, tau, b)
    if st == 'zero':
        verdict, rung = 'zero', 'square'
    elif st == 'cofactor-big':
        verdict, info = cofactor_decide(a, z, tau, b, detail=dt)
        rung = info.get('rung')
    else:
        verdict, rung = st, None
    assert verdict == 'zero', (key, 'LADDER VERDICT', verdict)
    delta = 1 - (1 + 4*a*a)*tau*tau
    alpha = -delta*(1 + 4*a*a)
    c0 = _sun_h(a, b, z**3)
    M0 = 16 - delta*c0*c0 - 32*(1 + 4*a*a)*b*((a-1)/2)**2
    x0, d0 = alpha*M0, alpha*2*b
    try:
        ram = ramified(x0, d0)
        assert ram == [], (key, 'ramified NONEMPTY', ram)
        rs = 'empty'
        ram_ok += 1
    except (FactorBudget, PrimalityBound) as e:
        rs = 'refused:' + type(e).__name__
        ram_ref += 1
    out.write(json.dumps({'cell': r['cell'], 'family': fam,
                          'ladder': verdict, 'rung': rung,
                          'ramified': rs}) + '\n')
    out.flush()
    n_ok += 1
    if n_ok % 25 == 0:
        print(f"{n_ok}/{len(D['records'])} replayed "
              f"({time.time()-t0g:.0f}s)", flush=True)
out.write(json.dumps({'summary': {'n': n_ok, 'ramified_empty': ram_ok,
                                  'ramified_refused': ram_ref,
                                  'wall_s': round(time.time()-t0g, 1)}}) + '\n')
out.close()
print(f"REPLAY-ALL COMPLETE: {n_ok}/293 ladder-zero verified; ramified "
      f"empty {ram_ok}, refused {ram_ref} ({time.time()-t0g:.0f}s)",
      flush=True)
