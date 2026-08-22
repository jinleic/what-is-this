#!/usr/bin/env python3
"""L11 branch-completion probe (2026-08-16).

Search-side tool (peer of l9_steer.py): NOT part of the verified suite.
Findings are frozen into h10q.py only after in-suite re-verification.

The tied block Theta_* currently disjoins exactly two branches,
tau = 0 (alpha = -A) and tau = 2a/A (alpha = -1), where alpha = -delta_tau*A
and delta_tau = 1 - A tau^2.  Two more rational tau give the remaining two
square classes of <-1, A>:

    tau = 1          -> delta = -4a^2,      alpha = 4a^2*A  == A  (mod squares)
    tau = (A+1)/(2A) -> delta = -4a^4/A,    alpha = 4a^4    == 1  (mod squares)

Both clear to polynomials (multiply the second by A), so the disjunction is
still a single product of conics in the same (y, r): the witness count is
untouched.  Soundness is tau-uniform (Sun identities 2.2-2.4 hold for
arbitrary tau; Prop 2.1 never mentions tau), which this file stress-tests
against bad z, i.e. z with numerator +-2^k (no odd prime), where the
definition must NOT fire.
"""
import json
import random
import sys
import time
from fractions import Fraction as F

from h10q import (_sun_h, _sun_tied_solvable, _sun_target_tau, vp, hilbert,
                  FactorBudget, PrimalityBound, factorint, primerange,
                  legendre)

TAUS = {
    'tau0': lambda a, A: F(0),          # alpha = -A   (canonical, w = 3 mod 4)
    'tauC': lambda a, A: 2 * a / A,     # alpha = -1   (canonical, w = 1 mod 4)
    'tau1': lambda a, A: F(1),          # alpha = A
    'tauD': lambda a, A: (A + 1) / (2 * A),   # alpha = square
}


def alpha_of(a, tau):
    A = 1 + 4 * a * a
    return -(1 - A * tau * tau) * A


def bad_z_pool():
    """z NOT in any m_w for odd w  <=>  numerator of z is +-2^k."""
    return [F(sg * 2 ** k, d) for sg in (1, -1) for k in range(0, 5)
            for d in (1, 3, 5, 7, 9, 15, 21, 25)]


def b_pool(nmax=25, mmax=15):
    out = []
    for n in range(-nmax, nmax + 1, 2):
        for m in range(1, mmax + 1, 2):
            q = F(n, m)
            if q not in (0, 1) and q not in out:
                out.append(q)
    out.sort(key=lambda q: (abs(q.numerator) * q.denominator, q < 0))
    return out


def s_pool():
    return [F(n, m) for m in (1, 3, 5, 7, 9) for n in range(-8, 9)
            if F(n, m) != F(-1, 2) and vp(1 + 2 * F(n, m), 2) == 0]


def soundness_probe(trials=4000, seed=20260816, budget=45.0):
    """Bad-z soundness: the new branches must never fire off m_w."""
    rng = random.Random(seed)
    Z, B, S = bad_z_pool(), b_pool(), s_pool()
    stat = {k: {'decided': 0, 'refused': 0, 'soluble': 0} for k in TAUS}
    hits = []
    t0 = time.time()
    for i in range(trials):
        if time.time() - t0 > budget:
            break
        s = rng.choice(S)
        a = 1 + 2 * s
        b = rng.choice(B)
        if vp(b, 2) != 0:
            continue
        z = rng.choice(Z)
        A = 1 + 4 * a * a
        c = _sun_h(a, b, z ** 3)
        if c is None:
            continue
        for nm, f in TAUS.items():
            tau = f(a, A)
            try:
                r = _sun_tied_solvable(a, b, c, tau)
            except (FactorBudget, PrimalityBound):
                stat[nm]['refused'] += 1
                continue
            stat[nm]['decided'] += 1
            if r:
                stat[nm]['soluble'] += 1
                if len(hits) < 8:
                    hits.append((nm, str(s), str(b), str(z)))
    return stat, hits


def target_place_probe(wmax=100, seed=7):
    """W1 at the target place on each branch: is alpha a square unit at w?

    Reports, per branch, the fraction of (w, a) pairs whose target-place
    criterion chi_w(alpha) = 1 holds.  tauD should be unconditional.
    """
    from h10q import legendre, primerange
    rng = random.Random(seed)
    S = s_pool()
    out = {k: [0, 0] for k in TAUS}
    for w in primerange(3, wmax):
        for s in S:
            a = 1 + 2 * s
            A = 1 + 4 * a * a
            for nm, f in TAUS.items():
                al = alpha_of(a, f(a, A))
                num, den = al.numerator, al.denominator
                if num % w == 0 or den % w == 0:
                    continue
                chi = legendre(num, w) * legendre(den, w)
                out[nm][1] += 1
                out[nm][0] += 1 if chi == 1 else 0
    return out


def class_scan(cells=None, qmax=400, out_path='data/l11_classes.jsonl'):
    """Aligned-class certificates for every grid cell on the square branch.

    On tau_d the square class of alpha is trivial, so the frozen set is
    {2} and the L10a structure applies at EVERY cell -- including the
    w = 3 mod 4 half that L10b walls off for the canonical prime-b family.
    """
    from h10q import _l10_class_cert, _l9_grid
    cells = cells or sorted(_l9_grid())
    S, QS = s_pool(), [q for q in primerange(41, qmax)]
    rows, misses = [], []
    t0 = time.time()
    for i, (w, ut) in enumerate(cells):
        u = F(*ut)
        z = F(w) * u
        found = None
        for s in sorted(S, key=lambda q: (abs(q.numerator) * q.denominator,)):
            a = 1 + 2 * s
            if a == 0 or vp(a, 2) != 0:
                continue
            A = 1 + 4 * a * a
            tau = (A + 1) / (2 * A)
            for eps in (1, -1):
                for q1 in QS:
                    r = _l10_class_cert(a, z, tau, eps, q1)
                    if r and r['ok']:
                        found = ((s.numerator, s.denominator), eps, q1, r['N'],
                                 [str(p) for p in r['S']])
                        break
                if found:
                    break
            if found:
                break
        if found:
            rows.append({'w': w, 'u': list(ut), 's': list(found[0]),
                         'eps': found[1], 'q1': found[2], 'N': found[3],
                         'S': found[4]})
        else:
            misses.append((w, list(ut)))
        if (i + 1) % 50 == 0:
            print(f"  [{i+1}/{len(cells)}] found={len(rows)} miss={len(misses)}"
                  f" {time.time()-t0:.0f}s", flush=True)
    with open(out_path, 'w') as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"class scan on the square branch: {len(rows)}/{len(cells)} cells"
          f" certified, {len(misses)} misses -> {out_path}")
    print("  misses:", misses[:20])
    return rows, misses


def assembly_scan(cells, budget=90.0, deep=False):
    """Global solubility of the tied conic on each branch, per cell.

    Reports which branches produce a PROVED-soluble admissible witness.
    """
    from h10q import _l9_steered_solvable
    S, B = s_pool(), b_pool(45 if deep else 25, 21 if deep else 15)
    out = {}
    for (w, ut) in cells:
        u = F(*ut)
        z = F(w) * u
        t0 = time.time()
        hits = {}
        for s in S:
            a = 1 + 2 * s
            if a == 0 or vp(a, 2) != 0:
                continue
            A = 1 + 4 * a * a
            for nm, f in TAUS.items():
                if nm in hits:
                    continue
                tau = f(a, A)
                for b in B:
                    if time.time() - t0 > budget:
                        break
                    if vp(b, 2) != 0:
                        continue
                    try:
                        ok, mech, wilds = _l9_steered_solvable(a, b, z, tau)
                    except (FactorBudget, PrimalityBound):
                        continue
                    if ok:
                        hits[nm] = (str(s), str(b), mech)
                        break
            if time.time() - t0 > budget:
                break
        out[(w, ut)] = hits
        print(f"  cell w={w} u={u}: " +
              (", ".join(f"{k}:{v[0]}|{v[1]}|{v[2]}" for k, v in hits.items())
               or "NONE"), flush=True)
    return out


def d_scan(cells, dpool=None, budget=120.0, deep=True):
    """Assembly search over the SQUARE-BRANCH FAMILY tau_d, not just d = 1.

    A finite set of d is a legal finite disjunction (L11a/L11c), so this
    measures whether the branch parameter is an independent assembly lever.
    """
    from h10q import _l9_steered_solvable
    dpool = dpool or [F(1), F(-1), F(2), F(-2), F(3), F(1, 2), F(-1, 2),
                      F(5), F(1, 3), F(3, 2), F(-3), F(4), F(2, 3)]
    S, B = s_pool(), b_pool(45 if deep else 25, 21 if deep else 15)
    for (w, ut) in cells:
        u = F(*ut)
        z = F(w) * u
        t0, hits = time.time(), []
        for s in S:
            a = 1 + 2 * s
            if a == 0 or vp(a, 2) != 0:
                continue
            A = 1 + 4 * a * a
            for d in dpool:
                lam = (A - d * d) / (2 * d)
                if lam == 0:
                    continue                 # alpha = 0 would force delta = 0
                tau = (A + d * d) / (2 * d * A)
                assert -(1 - A * tau * tau) * A == lam * lam
                for b in B:
                    if time.time() - t0 > budget:
                        break
                    if vp(b, 2) != 0:
                        continue
                    try:
                        ok, mech, wilds = _l9_steered_solvable(a, b, z, tau)
                    except (FactorBudget, PrimalityBound):
                        continue
                    if ok:
                        hits.append((str(s), str(d), str(b), mech))
                        break
                if hits:
                    break
            if hits or time.time() - t0 > budget:
                break
        print(f"  cell w={w} u={u}: " +
              (f"SOLVED s={hits[0][0]} d={hits[0][1]} b={hits[0][2]}"
               f" ({hits[0][3]})" if hits else "none"), flush=True)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if mode in ('all', 'sound'):
        budget = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
        stat, hits = soundness_probe(trials=200000, budget=budget)
        print("== bad-z soundness probe (definition must NOT fire) ==")
        for k, v in stat.items():
            print(f"  {k}: decided {v['decided']:5d}  refused {v['refused']:4d}"
                  f"  SOLUBLE {v['soluble']}   alpha class"
                  f" {'-A' if k=='tau0' else '-1' if k=='tauC' else 'A' if k=='tau1' else 'square'}")
        print("  hits:", hits)
        sys.stdout.flush()
    if mode in ('all', 'target'):
        out = target_place_probe()
        print("== target-place criterion chi_w(alpha) = 1 ==")
        for k, (g, n) in out.items():
            print(f"  {k}: {g}/{n} = {100.0*g/max(n,1):.1f}%")
        sys.stdout.flush()
    if mode == 'classes':
        qmax = int(sys.argv[2]) if len(sys.argv) > 2 else 400
        class_scan(qmax=qmax)
    if mode == 'dscan':
        from h10q import _l9_grid
        which = sys.argv[2] if len(sys.argv) > 2 else 'residual'
        budget = float(sys.argv[3]) if len(sys.argv) > 3 else 120.0
        resid = [(29, (-2, 1)), (31, (-2, 1)), (41, (1, 7)), (61, (2, 1)),
                 (61, (-2, 7)), (67, (2, 1)), (89, (-2, 1))]
        cells = resid if which == 'residual' else sorted(_l9_grid())[:30]
        print(f"== square-branch family scan (tau_d) on {len(cells)} cells ==")
        d_scan(cells, budget=budget)
    if mode == 'assembly':
        from h10q import _l9_grid
        which = sys.argv[2] if len(sys.argv) > 2 else 'residual'
        deep = len(sys.argv) > 4 and sys.argv[4] == 'deep'
        budget = float(sys.argv[3]) if len(sys.argv) > 3 else 90.0
        if which == 'residual':
            cells = [(29, (-2, 1)), (31, (-2, 1)), (41, (1, 7)), (53, (1, 3)),
                     (61, (2, 1)), (61, (-2, 7)), (67, (2, 1)), (89, (-2, 1))]
        elif which == 'w3':
            cells = [k for k in sorted(_l9_grid()) if k[0] % 4 == 3][:24]
        else:
            cells = sorted(_l9_grid())
        print(f"== global assembly scan on {len(cells)} cells"
              f" ({'deep' if deep else 'standard'} pools) ==")
        assembly_scan(cells, budget=budget, deep=deep)


if __name__ == '__main__':
    main()
