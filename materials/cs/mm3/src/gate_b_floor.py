"""Gate B core — direction-count floor + complete subset-DFS floor decision.

Model (paper Section 2): linear SLP, inputs free (for a factor map: the 9
matrix entries e_0..e_8; gates = x+y / x-y of previous quantities cost 1;
sign changes and copies free).

Lemma (floor): every gate of any valid circuit creates at most one new
output DIRECTION (sign class of its value). Hence C(F) >= d(F), where d(F)
counts needed nonzero direction classes minus +/- input directions. If a
circuit attains exactly d(F) gates, EVERY gate creates a new needed class,
so every gate value lies in {+/- t : t a needed target}. Therefore:
    a floor schedule exists  <=>  C(F) <= d(F).
This module decides exactly that question by exhaustive subset search with
memoization (state = bitmask of created classes; transitions = creatable
classes given available values). Complete and terminating: 2^d states.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

from tensor_data import U_BLOCK_PRINTED, V_BLOCK_PRINTED, W_BLOCK_PRINTED

R = 23
N = 9


def factor_targets():
    U, V, W = [], [], []
    for r in range(R):
        U.append(tuple(U_BLOCK_PRINTED[i][r] for i in range(N)))
        V.append(tuple(V_BLOCK_PRINTED[i][r] for i in range(N)))
        W.append(tuple(W_BLOCK_PRINTED[i][r] for i in range(N)))
    return U, V, W


def canon(v):
    v = tuple(v)
    neg = tuple(-x for x in v)
    return v if v <= neg else neg


def prep(targets):
    """Returns (needed classes list, REP lists)."""
    input_dirs = {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}
    classes = sorted({canon(t) for t in targets if any(t)} - input_dirs)
    base = sorted(set(classes) | input_dirs)
    reps = {}
    for c in classes:
        found = []
        for a in base:
            for b in base:
                ok = False
                for sa in (1, -1):
                    for sb in (1, -1):
                        s = tuple(sa * a[i] + sb * b[i] for i in range(N))
                        if s == c or s == tuple(-x for x in c):
                            ok = True
                            break
                    if ok:
                        break
                if ok:
                    found.append((a, b))
        # dedupe unordered pairs
        seen = set()
        dd = []
        for pr in found:
            key = tuple(sorted(pr))
            if key not in seen:
                seen.add(key)
                dd.append(pr)
        reps[c] = dd
    return classes, reps


def subset_dfs(classes, reps):
    """Complete search: can all classes be created? Returns
    (achievable, trace) where trace records reachable-state statistics."""
    idx = {c: i for i, c in enumerate(classes)}
    n = len(classes)
    inputs_mask = 0  # inputs always available, not tracked
    # REP as bit-pairs over indices (with inputs as sentinel index n+...)
    IN = n  # sentinel used only in data below
    base = sorted(set(classes) | {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)})
    bidx = {c: i for i, c in enumerate(base)}
    rep_pairs_idx = {}
    for c in classes:
        prs = []
        for (a, b) in reps[c]:
            ia, ib = bidx[a], bidx[b]
            prs.append((ia, ib))
        rep_pairs_idx[idx[c]] = prs
    nbase = len(base)

    from functools import lru_cache
    seen_states = set()
    order = []

    def avail(mask_pair, classpos):
        """class with base-index `classpos` available given created-class
        bitmask+input usage? Implemented by checking: is classpos an input
        (always available) or its class created?"""
        pass

    # Instead: represent availability of base classes in a bitmask av over
    # base (nbase bits). Inputs set initially.
    init_av = 0
    for c in base:
        if c in {canon(tuple(1 if j == i else 0 for j in range(N))) for i in range(N)}:
            init_av |= 1 << bidx[c]
    full = (1 << n) - 1

    memo_ok = {}
    stats = {"states": 0}

    def can_complete(av, created):
        key = (av, created)
        if key in memo_ok:
            return memo_ok[key]
        stats["states"] += 1
        if created == full:
            memo_ok[key] = True
            return True
        result = False
        for ci in range(n):
            bit = 1 << ci
            if created & bit:
                continue
            # try to create class ci from available values
            for (ia, ib) in rep_pairs_idx[ci]:
                if (av >> ia) & 1 and (av >> ib) & 1:
                    nav = av | (1 << bidx[classes[ci]])
                    if can_complete(nav, created | bit):
                        result = True
                        order.append(classes[ci])
                        break
            if result:
                break
        memo_ok[key] = result
        return result

    ok = can_complete(init_av, 0)
    return ok, stats, list(reversed([c for c in order]))


def main():
    U_t, V_t, W_t = factor_targets()
    out = {}
    for name, T in (("U", U_t), ("V", V_t), ("Wfactor", W_t)):
        classes, reps = prep(T)
        d = len(classes)
        ok, stats, trace = subset_dfs(classes, reps)
        out[name] = {"d(F)": d, "floor_schedule_exists": ok,
                     "reachable_states": stats["states"],
                     "witness_order_if_any": trace if ok else None}
    print(json.dumps(out, indent=2))
    Path(__file__).parent.joinpath("gate_b_floor.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
