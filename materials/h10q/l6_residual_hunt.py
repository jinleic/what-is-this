#!/usr/bin/env python3
"""Deep witness hunt on the 6 grid cells with no known soluble tied conic.

v2: per-cell wall-clock budgets, progress prints, member/candidate caps.
Cells 4 class-aligned (L11h): enumerate certified class members Q = q1 + kN.
Cells 2 walled: wide (branch, eps, e, q <= 1000) scan at a = 1, then a = 3.
Every True is cross-checked with _l9_steered_solvable.  Refusals logged only.
"""
import time
from fractions import Fraction as F
from h10q import (_l10_class_cert, _l7_tied_status, _l9_steered_solvable,
                  vp, factorint, primerange, _is_prime, _L11_CLASSES,
                  FactorBudget, PrimalityBound)

T0 = time.time()
CELLS = [(41, (1, 7)), (61, (2, 1)), (61, (-2, 7)), (89, (-2, 1)),
         (31, (-2, 1)), (67, (2, 1))]


def check(a, b, z, tau, tag):
    try:
        st = _l7_tied_status(a, b, z, tau)
    except (FactorBudget, PrimalityBound):
        return 'budget'
    if st is True:
        try:
            st2 = _l9_steered_solvable(a, b, z, tau)
        except (FactorBudget, PrimalityBound):
            st2 = 'budget'
        print(f"*** WITNESS cell={tag} a={a} b={b} branch={tau} steered={st2}",
              flush=True)
        return True
    return st


def hunt_class(w, ut, budget):
    z, a_, eps, q1, N, p = F(w)*F(*ut), *_L11_CLASSES[(w, ut)]
    a = F(*a_)
    A = 1 + 4*a*a
    tau = (1 + 2*a*a)/A
    r = _l10_class_cert(a, z, tau, eps, q1)
    assert r is not None and r['ok'] and r['N'] == N
    t = q1
    tested = refused = falses = 0
    for k in range(1, 4001):
        t += N
        if time.time() > budget:
            break
        try:
            if not _is_prime(t):
                continue
        except (FactorBudget, PrimalityBound):
            refused += 1
            continue
        tested += 1
        if tested % 200 == 0:
            print(f"  {(w,ut)}: {tested} members, {falses} decided-False,"
                  f" {refused} refusals, k={k}", flush=True)
        ans = check(a, F(eps*t), z, tau, (w, ut))
        if ans is True:
            return ('WITNESS', k, eps*t, tested, refused)
        if ans == 'budget':
            refused += 1
        else:
            falses += 1
    return ('none', None, None, tested, refused, falses)


def hunt_walled(w, ut, budget):
    z = F(w)*F(*ut)
    zs = [f for f in sorted(factorint(abs(z.numerator))) if f != 2]
    tested = refused = falses = 0
    for a_ in (1, 3):
        a = F(a_)
        A = 1 + 4*a*a
        for tau in ((1 + 2*a*a)/A, 2*a/A, F(1)):
            for f in zs:
                for e_ in (1, 2):
                    for eps in (1, -1):
                        g = f**e_ * eps
                        for q in primerange(41, 1001):
                            if time.time() > budget:
                                print(f"  {(w,ut)}: budget out at a={a_}"
                                      f" tau={float(tau):.2f} f={f} e={e_}"
                                      f" eps={eps} q={q} ({tested} tested)",
                                      flush=True)
                                return ('none', None, None, tested, refused,
                                        falses)
                            b = F(g*q)
                            if vp(b, 2) != 0:
                                continue
                            tested += 1
                            if tested % 400 == 0:
                                print(f"  {(w,ut)}: {tested} tested,"
                                      f" {falses} False, {refused} refused",
                                      flush=True)
                            ans = check(a, b, z, tau, (w, ut))
                            if ans is True:
                                return ('WITNESS', (a_, float(tau), g, q),
                                        None, tested, refused)
                            if ans == 'budget':
                                refused += 1
                            else:
                                falses += 1
    return ('done-none', None, None, tested, refused, falses)


if __name__ == '__main__':
    print(f"residual hunt v2 on {len(CELLS)} cells", flush=True)
    per = 60.0*9
    for w, ut in CELLS:
        cell_t0 = time.time()
        if (w, ut) in _L11_CLASSES:
            res = hunt_class(w, ut, cell_t0 + per)
            kind = 'class'
        else:
            res = hunt_walled(w, ut, cell_t0 + per)
            kind = 'walled'
        print(f"== cell {(w,ut)} [{kind}]: {res}", flush=True)
    print("done", flush=True)
