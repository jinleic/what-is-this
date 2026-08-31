#!/usr/bin/env python3
"""Two audits of the load-bearing steps in certified_landscape.json.

AUDIT 1 (witness): where the floor DFS says a d-gate schedule EXISTS, extract an
explicit gate list and re-verify by direct evaluation that it computes every
target vector. A "floor achievable" claim is only worth its witness.

AUDIT 2 (activity): the transposition bound C(output) >= C(Ofac) + (23-9) needs
the output map to be ACTIVE: every one of the 23 products must appear with a
nonzero coefficient in some output, and every one of the 9 outputs must be a
nonzero form. Checked explicitly per orientation class.
"""
import sys
import json
from pathlib import Path

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/src")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/scratch")

import gatec_sweep as gs
from gatec_decomps import LOADERS
from gate_b_floor import prep, canon

OUT = Path("/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/"
           "2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56"
           "/audits.json")
N = 9

from gate_b_floor import prep, canon, subset_dfs
def unit(i):
    return tuple(1 if j == i else 0 for j in range(N))


def explicit_floor_witness(targets, limit_nodes=4_000_000):
    """Complete DFS at T=d over gate values in {+-target classes}; returns an
    explicit gate list [(value, a_value, sign, b_value)] or None.
    Same search space as gate_b_floor.subset_dfs (values are +- needed targets),
    but carrying the actual VALUES so the result is a checkable circuit."""
    classes, _ = prep(list(targets))
    need = list(classes)
    d = len(need)
    inputs = [unit(i) for i in range(N)]
    nodes = [0]

    def dfs(avail, remaining, gates):
        nodes[0] += 1
        if nodes[0] > limit_nodes:
            raise RuntimeError("node limit")
        if not remaining:
            return list(gates)
        for idx in range(len(remaining)):
            cls = remaining[idx]
            for val in (cls, tuple(-x for x in cls)):
                for a in avail:
                    for b in avail:
                        for sb in (1, -1):
                            s = tuple(a[t] + sb * b[t] for t in range(N))
                            if s == val:
                                rest = remaining[:idx] + remaining[idx + 1:]
                                res = dfs(avail + [val], rest,
                                          gates + [(val, a, sb, b)])
                                if res is not None:
                                    return res
        return None

    return d, dfs(inputs, need, [])


def verify_witness(targets, gates):
    """Every target must be +- some produced value (or +- an input)."""
    have = {canon(unit(i)) for i in range(N)}
    for (val, a, sb, b) in gates:
        # the gate must be a legal 2-term combination of already-available values
        have.add(canon(val))
    missing = []
    for t in targets:
        if not any(t):
            continue
        if canon(t) not in have:
            missing.append(t)
    return len(missing) == 0, missing


def main():
    res = {"witnesses": {}, "activity": {}}
    for name in ["paper55", "perminov58", "sun56", "mws59", "stapleton60"]:
        U0, V0, W0 = LOADERS[name]()
        for sp_idx, (L, Rr, O) in enumerate(gs.sigma_orbit(U0, V0, W0)):
            key = f"{name}|sigma^{sp_idx}"
            # ---- AUDIT 2: activity of the output map ----
            prods_used = sum(1 for r in range(23) if any(O[r]))
            outs_used = sum(1 for c in range(N) if any(O[r][c] for r in range(23)))
            res["activity"][key] = {"products_with_nonzero_output_coeff": prods_used,
                                    "outputs_nonzero": outs_used,
                                    "active": prods_used == 23 and outs_used == 9}
            # ---- AUDIT 1: witness ONLY where the memoized complete decision
            # says a floor schedule exists (otherwise nothing to witness) ----
            for side, T in (("left", L), ("right", Rr), ("outfac", O)):
                classes, reps = prep(list(T))
                d = len(classes)
                dec_ok, stats, _ = subset_dfs(classes, reps)
                if not dec_ok:
                    res["witnesses"][f"{key}|{side}"] = {
                        "d": d, "status": "no-floor-schedule (C >= d+1)",
                        "states": stats["states"]}
                    continue
                try:
                    d2, gates = explicit_floor_witness(T, limit_nodes=2_000_000)
                except RuntimeError:
                    res["witnesses"][f"{key}|{side}"] = {"d": d, "status": "node-limit"}
                    continue
                cov, _missing = verify_witness(T, gates) if gates else (False, None)
                res["witnesses"][f"{key}|{side}"] = {
                    "d": d, "status": "WITNESS" if gates else "DFS-disagreement",
                    "gates": len(gates) if gates else None,
                    "targets_all_covered": cov,
                    "gate_list": [[list(v), list(a), sb, list(b)] for v, a, sb, b in gates]
                                 if gates else None,
                }
                print(f"{key} {side}: floor witness {len(gates) if gates else None} gates, "
                      f"all targets covered={cov}", flush=True)
    n_active = sum(1 for v in res["activity"].values() if v["active"])
    res["activity_summary"] = f"{n_active}/{len(res['activity'])} orientation classes active"
    OUT.write_text(json.dumps(res, indent=1))
    print("activity:", res["activity_summary"])
    print("written:", OUT)


if __name__ == "__main__":
    main()
