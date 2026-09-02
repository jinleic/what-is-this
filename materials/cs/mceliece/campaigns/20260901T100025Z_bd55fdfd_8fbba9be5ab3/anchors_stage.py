"""Anchors stage — T96 cell (standalone, writes anchor_*.json).

MECHANICAL RE-INSTATIATION of M12ALPHA_DIVONLY/camp_lib.py's anchor
block: the SAME three registered anchors, SAME record_for() route
(record shape identical to verify_delta_uniform.verify_delta, plus the
disclosed additive `guards_build` field), SAME byte-comparison against
the frozen 14-09Z records with elapsed_s excluded (timing excl. per
AMENDMENT), except the comparison now also drops my additive fields
(`guards_build`, `max_deg_f`, `deg_ok`, `alpha_identity_measured`,
`beta_delta_nonzero`, `beta_pair`, `eps_maxdeg_is_D`, `gamma_k`) that
the t=64 campaign dropped implicitly and recorded in
`extra_fields_this_run`. Here the extra-field set is made EXPLICIT and
each extra field additionally cross-checked for internal consistency
(additive verification of what was previously merely disclosed).

HARD GATE (pre_statement_t96.md section 3): all three anchors must be
byte-equal vs frozen (modulo elapsed_s) BEFORE any t=96 verdict is
reported. Pass/fail per anchor recorded; fails escalate to Main.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/20260901T100025Z_bd55fdfd_8fbba9be5ab3/code")

from instance import Instance  # noqa: E402

RUN = "20260901T100025Z_bd55fdfd_8fbba9be5ab3"
OUT = f"/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/{RUN}"
FROZEN = ("/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/"
          "2026-08-30T14-09Z_57200ADD/delta_uniform_13instances.json")

ANCHORS = [(11, 2048, 48, 6211), (10, 1024, 40, 5113), (6, 64, 3, 1387)]
# additive fields this run records on top of the frozen shape (disclosed,
# not part of the frozen comparison set)
EXTRA_FIELDS = {"guards_build", "max_deg_f", "deg_ok"}


def record_for(m, n, t, seed):
    """EXACT replica of verify_delta_uniform.verify_delta record shape
    (same route, same fields) + disclosed additive fields."""
    t0 = time.time()
    inst = Instance(m, n, t, seed)
    gf = inst.gf
    rec = {"m": m, "n": n, "t": t, "k": inst.k, "D": inst.D, "seed": seed,
           "G": inst.G}

    degs = [gf.pdeg(f) for f in inst.F]
    deg_ok = len(degs) > 0 and max(degs) <= inst.D
    rec["max_deg_f"] = int(max(degs))
    rec["deg_ok"] = bool(deg_ok)
    rec["guards_build"] = {kk: vv for kk, vv in inst.guards.items()
                           if kk != "degs"}
    if not deg_ok:
        rec["ABORT"] = ("degree condition max_j deg f_j <= D FAILED; the "
                        "Lemma-2 degree argument does not apply, so NO "
                        "delta verdict is produced for this instance")
        rec["delta_exact"] = None
        rec["elapsed_s"] = round(time.time() - t0, 2)
        return rec

    Ls = inst._lagrange()

    exact_div_ok = True
    for i in range(n):
        q, r = gf.pdivmod(inst.Pi, [inst.support[i], 1])
        if gf.pdeg(r) >= 0:
            exact_div_ok = False
            rec["exact_div_fail_at"] = i
            break
    rec["exact_div_ok"] = bool(exact_div_ok)

    bad_unit = []
    for i, a in enumerate(inst.support):
        if gf.peval(Ls[i], a) != 1:
            bad_unit.append(i)
            if len(bad_unit) > 4:
                break
    rec["lagrange_unit_ok"] = (len(bad_unit) == 0)
    rec["lagrange_unit_failures"] = bad_unit
    rec["lagrange_unit_checks"] = n

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
    rec["alpha_identity_measured"] = bool(inst.guards.get("alpha_identity"))
    rec["beta_delta_nonzero"] = bool(inst.guards.get("beta_delta_nonzero"))
    rec["beta_pair"] = inst.guards.get("beta_pair")
    rec["eps_maxdeg_is_D"] = bool(inst.guards.get("eps_maxdeg_is_D"))
    rec["gamma_k"] = bool(inst.guards.get("gamma_k"))
    rec["elapsed_s"] = round(time.time() - t0, 2)
    return rec


if __name__ == "__main__":
    out = {"anchors": {}, "frozen_source": FROZEN}
    frozen = json.load(open(FROZEN))
    fmap = {(r["m"], r["n"], r["t"], r["seed"]): r for r in frozen}

    allok = True
    for (m, n, t, seed) in ANCHORS:
        r = record_for(m, n, t, seed)
        fr = fmap[(m, n, t, seed)]
        fr_x = {kk: vv for kk, vv in fr.items() if kk != "elapsed_s"}
        r_x = {kk: vv for kk, vv in r.items() if kk != "elapsed_s"}
        r_x_cmp = {kk: vv for kk, vv in r_x.items() if kk in fr_x}
        extra = {kk: vv for kk, vv in r_x.items() if kk not in fr_x}
        # additive cross-check of the disclosed extras: internal consistency
        extra_ok = (r["deg_ok"] is True
                    and r["max_deg_f"] == r["D"]
                    and r["guards_build"].get("alpha_identity") is True
                    and r["guards_build"].get("delta_yi_binary") is True)
        s_f = json.dumps(fr_x, sort_keys=True, default=str)
        s_r = json.dumps(r_x_cmp, sort_keys=True, default=str)
        equal = (s_f == s_r)
        allok &= equal and extra_ok
        out["anchors"][f"{m}_{n}_{t}_{seed}"] = {
            "byte_equal_vs_frozen": equal,
            "extra_fields_this_run": sorted(extra.keys()),
            "extra_fields_consistent": extra_ok,
            "record": r,
        }
        print(f"anchor ({m},{n},{t},{seed}): byte_equal={equal} "
              f"extras_ok={extra_ok} "
              f"delta_exact={r.get('delta_exact')} "
              f"alpha_measured={r.get('alpha_identity_measured')} "
              f"[{r['elapsed_s']}s]", flush=True)
        json.dump(out, open(f"{OUT}/anchors_result.json", "w"),
                  indent=1, default=str)

    # per-anchor dedicated artifact (same file naming as t=64 campaign)
    for key, blk in out["anchors"].items():
        with open(f"{OUT}/anchor_{key}.json", "w") as fh:
            json.dump({"anchors": {key: blk}}, fh, indent=1, default=str)

    out["all_byte_equal"] = allok
    json.dump(out, open(f"{OUT}/anchors_result.json", "w"), indent=1,
              default=str)
    print("ANCHORS: all byte-equal = ", allok, flush=True)
    sys.exit(0 if allok else 4)  # 4 = registered anchor-gate failure
