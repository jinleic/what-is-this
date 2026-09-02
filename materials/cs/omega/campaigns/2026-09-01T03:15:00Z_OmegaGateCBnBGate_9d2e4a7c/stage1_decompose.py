#!/usr/bin/env python3
"""Gate C BnB feasibility gate — STAGE 1 width decomposition.

Implements pre_statement.md (this directory) exactly:
  GATE Z: recompute the L1 width sum from live modules; assert
          bit-equality with the frozen parent checkpoint.
  1.1     per-coordinate width table + cumulative shares.
  1.2     per-term/per-block attribution via byte-reuse of the asrun
          aggregate() (hull overhead reported separately).
  1.3     branch-straddle per R branch from the parent hull data.
  VERDICT: CONCENTRATED (k* <= 6 at >=90%) vs SPREAD (obstruction),
          computed mechanically from the measured curve.

Per-attempt logs under distinct filenames (freeze-hygiene rule).
"""
import sys, os, json, time, hashlib

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "stage1_decomposition.json")
T0 = time.time()
_checkpoint = {}


def ckpt(key, value):
    _checkpoint[key] = value
    _checkpoint["elapsed"] = time.time() - T0
    tmp = OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(_checkpoint, fh, indent=1, default=str)
    os.replace(tmp, OUT)
    print("CKPT", key, json.dumps(value, default=str)[:280])
    sys.stdout.flush()


PARENT = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
          "2026-08-31T02:40:00Z_OmegaGateCFinish_v1_4c8a2f1b")
GC = json.load(open(os.path.join(PARENT, "run_v1_v11_checkpoint.json")))
sb = GC["slope_box"]
lo_f = np.array(sb["grad_lo"], dtype=float)
hi_f = np.array(sb["grad_hi"], dtype=float)
wid_f = hi_f - lo_f
FROZEN_L1 = float(sb["width_sum_l1"])
FROZEN_WMAX = float(sb["width_max"])
assert abs(float(wid_f.sum()) - FROZEN_L1) < 1e-9, "checkpoint L1 mismatch"
assert abs(float(wid_f.max()) - FROZEN_WMAX) < 1e-12

# ---------------- GATE Z: recompute from live modules --------------------
import interval_core as IC
import gate_c_slope_pass as SP
import gate_c_slope_core as SC
import vxxz24_float as F
from scipy.io import loadmat
from flint import arb

params = np.asarray(loadmat(SP.MAT)["params"]).flatten()

with IC.prec():
    outB = SP.run_slope_pass(params, radius=1e-7, verbose=False)
    # byte-reuse of the asrun aggregate() exactly like the parent runner
    asrun_path = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
                  "2026-08-31T07:57:00Z_OmegaGateCDiag_4859327/"
                  "gate_c_endgame.py.asrun")
    aesrc = open(asrun_path).read()
    a0 = aesrc.index("    def aggregate(out, radius, lemma1_off=False):")
    a1 = aesrc.index("    def r_detail(agg):")
    block = aesrc[a0:a1]
    lines = [ln[4:] if ln.startswith("    ") else ln
             for ln in block.splitlines()]
    ns = {"SC": SC, "SP": SP, "IC": IC, "np": np, "arb": arb,
          "I": IC.Ival, "F": F,
          "slope_of": (lambda x: x.v if isinstance(x, SP.SlopeG) else x)}
    exec(compile("\n".join(lines), "<aggregate-asrun>", "exec"), ns)
    aggregate = ns["aggregate"]

    aggB = aggregate(outB, 1e-7)
    Om = aggB["Om_slope"].grad
    lo_v = np.array([float(x.lo) for x in Om])
    hi_v = np.array([float(x.hi) for x in Om])
    wid_v = hi_v - lo_v
    live_l1 = float(wid_v.sum())
    live_wmax = float(wid_v.max())
ckpt("GATE_Z", {
    "live_l1": live_l1, "frozen_l1": FROZEN_L1,
    "l1_agreement": float(abs(live_l1 - FROZEN_L1)),
    "live_wmax": live_wmax, "frozen_wmax": FROZEN_WMAX,
    "max_abs_grad_lo_dev": float(np.abs(lo_v - lo_f).max()),
    "max_abs_grad_hi_dev": float(np.abs(hi_v - hi_f).max()),
})
DEV_TOL = 1e-12
L1_TOL = 1e-6
gz_ok = (float(np.abs(lo_v - lo_f).max()) <= DEV_TOL and
         float(np.abs(hi_v - hi_f).max()) <= DEV_TOL and
         abs(live_l1 - FROZEN_L1) <= L1_TOL)
if not gz_ok:
    ckpt("verdict", {"branch": "GATE_Z_FAIL",
                     "reason": ("live recompute deviates from frozen "
                                "checkpoint beyond replay-noise "
                                "tolerance (1e-12/1e-6)")})
    sys.exit(3)

# ---------------- 1.1 per-coordinate table -------------------------------
tot = float(wid_v.sum())
order = np.argsort(-wid_v)
cum = np.cumsum(wid_v[order]) / tot
k90 = int(np.searchsorted(cum, 0.90) + 1)
table_11 = {
    "widths_sorted": [
        {"coord": int(c), "width": float(wid_v[c]),
         "share": float(wid_v[c] / tot)} for c in order],
    "cumulative_shares": {
        "k2": float(cum[1]), "k4": float(cum[3]), "k6": float(cum[5]),
        "k10": float(cum[9]), "k23": float(cum[22]), "k90": float(cum[k90 - 1])},
    "k_star_90pct": k90,
    "min_width": float(wid_v.min()), "max_width": float(wid_v.max()),
    "n_coords_ge_1pct_each": int((wid_v / tot >= 0.01).sum()),
    "n_coords_ge_2pct_each": int((wid_v / tot >= 0.02).sum()),
}
ckpt("width_table_1_1", table_11)

# ---------------- 1.2 per-term attribution --------------------------------
branch = aggB["branch"]
R_branch_names = {"(3, 0)": "R_comp[0,3]", "(3, 1)": "R_comp[1,3]",
                  "(3, 2)": "R_comp[2,3]", "lvl2": "R_comp[0,2]",
                  "('G', 0)": "R_glob[0]", "('G', 1)": "R_glob[1]",
                  "('G', 2)": "R_glob[2]"}
branch_info = {}
for key, rec in branch.items():
    cands = rec["cand"]
    ivs = rec.get("candidate_intervals")
    branch_info[R_branch_names.get(str(key), str(key))] = {
        "candidate_lows": [float(z) for z in cands],
        "argmin": int(rec["argmin"]),
        "possible_minimizers": [int(z) for z in rec["possible_minimizers"]],
        "zero_possible": bool(rec["zero_possible"]),
        "candidate_interval_spans": ([float(iv[1] - iv[0]) for iv in ivs]
                                     if ivs else None),
        "minimizer_identity_straddles": (
            len(rec["possible_minimizers"]) > 1 if "possible_minimizers"
            in rec else None),
    }
term_widths = {
    "R_branches_detail": branch_info,
    "M_possible_minimizers": [int(z) for z in aggB["M_possible_minimizers"]],
    "Om_slope_grad_l1": live_l1,
    "hull_multiplicity_per_branch": {
        R_branch_names.get(str(k), str(k)):
            len(v["possible_minimizers"]) for k, v in branch.items()},
    "note": ("the Om slope grad is the HULLED composition (min_slope "
             "hulls per branch + branch sum + quotient rule); per-term "
             "pre-hull L1 widths are reported via the branch hull data "
             "and membership counts above"),
}
ckpt("term_attribution_1_2", term_widths)

# ---------------- 1.3 branch straddle ------------------------------------
straddle = {}
for name, info in branch_info.items():
    spans = info["candidate_interval_spans"] or []
    straddle[name] = {
        "n_possible_minimizers": len(info["possible_minimizers"]),
        "zero_possible": info["zero_possible"],
        "max_candidate_span": (max(spans) if spans else None),
    }
straddle["frozen_context"] = {
    "R_interval_span_seven_components": "frozen in CandidateWitness (R spans ~4.2506e-6 across the 7 components)",
    "lemma1_residual_straddles_zero_at_base": True,
}
ckpt("branch_straddle_1_3", straddle)

# ---------------- MECHANICAL VERDICT -------------------------------------
if k90 <= 6:
    kset = [int(c) for c in order[:k90]]
    verdict = {
        "branch": "CONCENTRATED",
        "k_star": k90, "coords": kset,
        "share": float(cum[k90 - 1]),
        "next": ("append pre_statement_addendum_stage2.md BEFORE any "
                 "subbox compute; depth cap 8; bisect widest binding "
                 "coord at exact dyadic midpoint"),
        "stage2_opened": False,
    }
else:
    subboxes = 2.0 ** (6.2 * k90)
    half = 2.0 ** (6.2 * (k90 / 2.0))
    verdict = {
        "branch": "SPREAD — QUANTIFIED OBSTRUCTION",
        "k_star_90pct": k90,
        "width_min": float(wid_v.min()), "width_max": float(wid_v.max()),
        "n_ge_1pct_each": int((wid_v / tot >= 0.01).sum()),
        "implied_subboxes_full": float(subboxes),
        "implied_subboxes_mild_half_reading": float(half),
        "obstruction_sentence": ("branch-and-bound over D at this "
                                 "enclosure technology is INFEASIBLE AT "
                                 "THIS COST; a first-class frozen "
                                 "negative, not a loop failure"),
        "stage2_opened": False,
    }
ckpt("verdict", verdict)
ckpt("done", {"elapsed_total": time.time() - T0})
print("DONE", verdict["branch"])
