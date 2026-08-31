#!/usr/bin/env python3
"""Certified per-side cost invariance under the signed-permutation sandwich.

THEOREM (monomial transfer). Let g be a signed permutation of the 9 input
coordinates (equivalently a monomial 9x9 matrix with entries in {-1,0,1},
exactly one nonzero per row/column). For any factor map F: 9 -> 23 given by 23
target vectors, C(F) = C(g·F) where C is the minimum number of x±y gates in the
model (inputs free, sign changes and copies free).

PROOF. A circuit for F is a sequence of gates whose values are integer vectors
in Z^9. Apply g to every wire value. Inputs map to ± inputs (free), each gate
g(x ± y) = g(x) ± g(y) is again a legal single gate, and the 23 targets map
exactly onto the 23 targets of g·F. So the image is a circuit for g·F with the
same gate count, giving C(g·F) <= C(F); applying g^{-1} (also monomial) gives
the converse. QED.

COROLLARY. d(F) = d(g·F) and the floor decision "does a d(F)-gate schedule
exist" has the same answer for F and g·F, so all certified quantities of gate B
transfer verbatim across the whole 48^3 signed-permutation sandwich group.

This script MACHINE-CHECKS the corollary on sampled sandwiches with the
memoized complete floor-DFS (the same procedure that produced the gate-B floor
results: states 33 / 116 / 66 on the paper's three blocks)."""
import sys
import json
import time
from pathlib import Path

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/src")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/scratch")

import gatec_sweep as gs
from gatec_decomps import LOADERS
from gate_b_floor import prep, subset_dfs, canon

OUT = Path("/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/"
           "2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56"
           "/invariance_check.json")


def floor_decide(targets):
    """Complete floor decision at T=d: returns (d, schedule_exists, states)."""
    classes, reps = prep(list(targets))
    d = len(classes)
    ok, stats, _ = subset_dfs(classes, reps)
    return d, ok, stats["states"]


def main():
    samples = [(0, 0, 0), (7, 13, 29), (3, 40, 11), (47, 2, 19), (12, 12, 12)]
    res = {}
    t0 = time.time()
    for name in ["paper55", "sun56", "mws59", "stapleton60"]:
        U0, V0, W0 = LOADERS[name]()
        base = None
        rows = []
        for (xi, yi, zi) in samples:
            X, Y, Z = gs.SP[xi], gs.SP[yi], gs.SP[zi]
            U2, V2, W2 = gs.sandwich(U0, V0, W0, X, Y, Z)
            trip = []
            for side, T in (("U", U2), ("V", V2), ("Wfac", W2)):
                d, ok, st = floor_decide(T)
                trip.append({"side": side, "d": d, "floor_exists": ok, "states": st})
            rows.append({"triple": [xi, yi, zi], "sides": trip})
            key = [(t["d"], t["floor_exists"]) for t in trip]
            if base is None:
                base = key
            elif key != base:
                print(f"INVARIANCE VIOLATION {name} at {(xi,yi,zi)}: {key} != {base}")
                json.dump(res, OUT.open("w"), indent=1)
                sys.exit(1)
            print(f"{name} {(xi,yi,zi)}: " +
                  " ".join(f"{t['side']}(d={t['d']},floor={t['floor_exists']},st={t['states']})"
                           for t in trip), flush=True)
        res[name] = {"invariant_key": [[k[0], k[1]] for k in base], "rows": rows}
    res["_theorem"] = ("monomial transfer: C(F) = C(gF) for signed-permutation g; "
                       "verified on 5 sandwich triples x 4 decompositions x 3 sides")
    res["_seconds"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(res, indent=1))
    print("INVARIANCE CHECK PASSED", OUT)


if __name__ == "__main__":
    main()
