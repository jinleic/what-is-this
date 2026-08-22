#!/usr/bin/env python3
"""L6 witness search: tied (a = 1+2s) bridge certificates for odd primes.

For target input z = w, evaluate the bridge at z^3 = w^3 and find
(s, b, tau) with
  - s 2-integral and w-integral, a = 1 + 2s, A = 1 + 4a^2,
  - b = m*w, m odd, m coprime to w  (so b is a 2-adic unit, v_w(b) = 1),
  - Delta = Delta_fin(A, 2b) is exactly {2,w} and c = h(a,b,z^3) is Delta-integral,
  - the TIED conic  -dA y^2 - 16B r^2 = 16 - d c^2 - 16AB s^2  (d = 1 - A tau^2,
    B = 2b, s = (a-1)/2) is solvable over Q  [Hasse-Minkowski, proven-primality
    factoring engine only].
A hit is a 6-unknown-architecture certificate: the Psi_tau witness s doubles
as the Phi-parameter a = 1+2s, eliminating one existential unknown.  In
--canonical mode the target-place theorem forces the single branch
tau = 0 (w = 3 mod 4) or tau = 2a/A (w = 1 mod 4).

Usage: python3 l6_search.py [--canonical] OUT.jsonl w1 w2 w3 ...
"""
import json
import sys
from fractions import Fraction as Fr

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from h10q import (FactorBudget, PrimalityBound, delta_fin, legendre, represents_q,
                  _sun_h, _sun_integral, _sun_target_tau, vp)


def _v_ok(x, p):
    return x == 0 or vp(x, p) >= 0


S_MENU = sorted(
    {Fr(0)} | {Fr(n, d) for d in (1, 3, 5, 7, 9, 11, 13)
               for n in range(-4, 5) if n and abs(Fr(n, d)) <= 3},
    key=lambda t: (abs(t), t.denominator, t))
TAU_MENU = [Fr(0), Fr(1, 5), Fr(1, 4), Fr(1, 3), Fr(2, 5), Fr(1, 2),
            Fr(3, 5), Fr(2, 3), Fr(3, 4), Fr(1), Fr(5, 4), Fr(3, 2),
            Fr(2), Fr(5, 2), Fr(3), Fr(4)]




def tied_conic_solvable(a, b, c, tau):
    A, B = 1 + 4 * a * a, 2 * b
    d = 1 - A * tau * tau
    if d == 0:
        return None
    s = (a - 1) / 2
    M = 16 - d * c * c - 16 * A * B * s * s
    return True if M == 0 else represents_q([-d * A, -16 * B], M)


def l6_search(w, kmax=15, canonical=False):
    bridge_arg = Fr(w) ** 3
    for s in S_MENU:
        if not (_v_ok(s, 2) and _v_ok(s, w)):
            continue
        a = 1 + 2 * s
        A = 1 + 4 * a * a
        if canonical and (vp(A, w) != 0 or legendre(A, w) != -1):
            continue
        for m in (x for k in range(1, kmax + 1, 2) for x in (k, -k)):
            if m % w == 0:
                continue
            b = Fr(m * w)
            c = _sun_h(a, b, bridge_arg)
            if c is None:
                continue
            try:
                D = delta_fin(A, 2 * b)
                if sorted(D) != [2, w] or not _sun_integral(c, D):
                    continue
            except (FactorBudget, PrimalityBound):
                continue
            taus = (_sun_target_tau(a, w),) if canonical else TAU_MENU
            for tau in taus:
                try:
                    if tied_conic_solvable(a, b, c, tau):
                        return s, b, tau, sorted(D)
                except (FactorBudget, PrimalityBound):
                    continue        # budget artifact must not kill the cell
    return None


def main():
    args = sys.argv[1:]
    canonical = bool(args and args[0] == "--canonical")
    if canonical:
        args = args[1:]
    out_path, ws = args[0], [int(x) for x in args[1:]]
    with open(out_path, "a") as out:
        for w in ws:
            import time
            t0 = time.time()
            hit = l6_search(w, canonical=canonical)
            rec = {"w": w, "sec": round(time.time() - t0, 1)}
            if hit is None:
                rec["found"] = False
            else:
                s, b, tau, D = hit
                rec.update(found=True,
                           s=[s.numerator, s.denominator],
                           b=[b.numerator, b.denominator],
                           tau=[tau.numerator, tau.denominator],
                           delta=D)
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print(f"w={w}: {rec}", flush=True)


if __name__ == "__main__":
    main()
