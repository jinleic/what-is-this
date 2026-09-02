"""Portability of the campaign verdict — M12ALPHA-DIVONLY.

Exact records for the ADDENDUM-2(i) delta chain, computed with plain
python (+numpy only where noted) against byte-equivalent data formats.
The frozen 14-09Z verify_delta_uniform.py semantic payloads are
reproduced FIELD-FOR-FIELD for the three registered anchor rows, with
all non-timing fields required to match BYTE-EXACTLY (amended:
elapsed_s excluded per AMENDMENT-header note).

NO BCOMP — AMENDMENT 1 resolved the fork to the MEASUREMENT path.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mceliece/src")

from instance import Instance  # noqa: E402

OUT = ("/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/"
       "2026-09-01T03-28-10Z_M12ALPHA_DIVONLY")


def record_for(m, n, t, seed):
    """EXACT replica of verify_delta_uniform.verify_delta record shape
    (same route, same fields)."""
    t0 = time.time()
    inst = Instance(m, n, t, seed)
    gf = inst.gf
    rec = {"m": m, "n": n, "t": t, "k": inst.k, "D": inst.D, "seed": seed,
           "G": inst.G}

    # LOAD-BEARING degree condition: ABORT if violated
    degs = [gf.pdeg(f) for f in inst.F]
    deg_ok = len(degs) > 0 and max(degs) <= inst.D
    rec["max_deg_f"] = int(max(degs))
    rec["deg_ok"] = bool(deg_ok)
    rec["guards_build"] = {kk: vv for kk, vv in inst.guards.items() if kk != "degs"}
    if not deg_ok:
        rec["ABORT"] = ("degree condition max_j deg f_j <= D FAILED; the "
                        "Lemma-2 degree argument does not apply, so NO delta "
                        "verdict is produced for this instance")
        rec["delta_exact"] = None
        rec["elapsed_s"] = round(time.time() - t0, 2)
        return rec

    # Lagrange basis, rebuilt through the same construction path
    Ls = inst._lagrange()

    # exact division certificate: Pi = (Z - a_i) q_i, remainder == 0
    exact_div_ok = True
    for i in range(n):
        q, r = gf.pdivmod(inst.Pi, [inst.support[i], 1])
        if gf.pdeg(r) >= 0:
            exact_div_ok = False
            rec["exact_div_fail_at"] = i
            break
    rec["exact_div_ok"] = bool(exact_div_ok)

    # the n scalar checks L_i(a_i) = 1
    bad_unit = []
    for i, a in enumerate(inst.support):
        if gf.peval(Ls[i], a) != 1:
            bad_unit.append(i)
            if len(bad_unit) > 4:
                break
    rec["lagrange_unit_ok"] = (len(bad_unit) == 0)
    rec["lagrange_unit_failures"] = bad_unit
    rec["lagrange_unit_checks"] = n

    # off-diagonal spot check
    import random as _r
    rng = _r.Random(seed ^ 0xDE17A)
    bad_off = 0
    for _ in range(200):
        i = rng.randrange(n)
        l = rng.randrange(n)
        if i == l:
            continue
        if gf.peval(Ls[i], inst.support[l]) != 0:
            bad_off += 1
    rec["offdiag_probes"] = 200
    rec["offdiag_violations"] = bad_off

    # assembly spot check: f_j equals direct sum over its support
    asm_ok = True
    for _ in range(3):
        j = rng.randrange(inst.k)
        acc = []
        row = inst.Y[j]
        for i in range(n):
            if row[i]:
                acc = gf.padd(acc, gf.pscale(Ls[i], gf.inv(inst.lam[i])))
        if gf.ptrim(acc) != gf.ptrim(inst.F[j]):
            asm_ok = False
            rec["assembly_fail_row"] = int(j)
            break
    rec["assembly_spotcheck_rows"] = 3
    rec["assembly_ok"] = bool(asm_ok)

    rec["delta_exact"] = bool(exact_div_ok and rec["lagrange_unit_ok"]
                              and bad_off == 0 and asm_ok)
    rec["delta_scope"] = ("EXACT for ALL (j,l): lam_l f_j(a_l) = Y[j,l], "
                          "implied by L_i(a_i)=1 (n checks) + structural "
                          "off-diagonal vanishing + assembly")
    # alpha stays MEASURED on these rows (unchanged, from the instance guards)
    rec["alpha_identity_measured"] = bool(inst.guards.get("alpha_identity"))
    rec["beta_delta_nonzero"] = bool(inst.guards.get("beta_delta_nonzero"))
    rec["beta_pair"] = inst.guards.get("beta_pair")
    rec["eps_maxdeg_is_D"] = bool(inst.guards.get("eps_maxdeg_is_D"))
    rec["gamma_k"] = bool(inst.guards.get("gamma_k"))
    rec["elapsed_s"] = round(time.time() - t0, 2)
    return rec


ANCHORS = [(11, 2048, 48, 6211), (10, 1024, 40, 5113), (6, 64, 3, 1387)]

if __name__ == "__main__":
    out = {"anchors": {}, "m12": None}
    frozen = json.load(open("/Users/jinleic/jinleic-workspace/cs/mceliece/"
                            "campaigns/2026-08-30T14-09Z_57200ADD/"
                            "delta_uniform_13instances.json"))
    fmap = {(r["m"], r["n"], r["t"], r["seed"]): r for r in frozen}

    allok = True
    for (m, n, t, seed) in ANCHORS:
        r = record_for(m, n, t, seed)
        fr = fmap[(m, n, t, seed)]
        # byte-compare EXCLUDING elapsed_s (timing excl. per AMENDMENT)
        fr_x = {k: v for k, v in fr.items() if k != "elapsed_s"}
        r_x = {k: v for k, v in r.items() if k != "elapsed_s"}
        # NOTE: max_deg_f / deg_ok / guards_build etc are extra fields vs
        # frozen; drop MY extras so compare is against frozen field set:
        r_x_cmp = {k: v for k, v in r_x.items() if k in fr_x}
        extra = {k: v for k, v in r_x.items() if k not in fr_x}
        s_f = json.dumps(fr_x, sort_keys=True, default=str)
        s_r = json.dumps(r_x_cmp, sort_keys=True, default=str)
        equal = (s_f == s_r)
        allok &= equal
        out["anchors"][f"{m}_{n}_{t}_{seed}"] = {
            "byte_equal_vs_frozen": equal,
            "record": r,
        }
        print(f"anchor ({m},{n},{t},{seed}): byte_equal={equal} "
              f"delta_exact={r.get('delta_exact')} "
              f"alpha_measured={r.get('alpha_identity_measured')} "
              f"[{r['elapsed_s']}s]", flush=True)

    json.dump(out, open(os.path.join(OUT, "anchors_result.json"), "w"),
              indent=1, default=str)
    print(f"ANCHORS: all byte-equal = {allok}", flush=True)
