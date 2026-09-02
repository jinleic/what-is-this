#!/usr/bin/env python3
"""Gate C V1 re-run under the corrected v11 interval core — ONE pass.

Implements pre_statement.md (this directory) decision tree, fixed before
any compute:
  Step 1  frozen-protocol control on the v11 core (stage_b_rung2, with and
          without Lemma-1) -> C-gate vs 2.3715538358350803 <= 5e-13.
  Step 2  live replay point pass (radius 0) + endgame aggregate ->
          V1-gate vs frozen raw <= 5e-13 and float-equality of the seven
          R branches / R / M against the GateCDiag pre-v11 records.
  Step 3  (V1 pass only) box pass r=1e-7, scaled midpoint LP (u-form),
          ANY-y rigorous width, width/signal ratio, decision line.
  Branch Z: control drift -> STOP, regression finding, nothing else runs.

No source file in src/ is imported-and-edited; everything is reused as it
stands (vxxz24 31657c92..., interval_core v11 ce70f959..., slope pass
c6b8cef2..., endgame baaecb8b...). Checkpoint writer is the only new code
and performs no arithmetic beyond json/float extraction of already computed
values.
"""
import sys, os, json, time, hashlib

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "run_v1_v11_checkpoint.json")
T0 = time.time()
_checkpoint = {}


def ckpt(key, value):
    _checkpoint[key] = value
    _checkpoint["elapsed"] = time.time() - T0
    tmp = OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(_checkpoint, fh, indent=1, default=str)
    os.replace(tmp, OUT)
    print("CKPT", key, json.dumps(value, default=str)[:300])
    sys.stdout.flush()


# ------------------------------------------------------------------ anchors
FROZEN_RAW = 2.3715538358350803
FROZEN_R = 2.8170035674609757
FROZEN_M = 2.0942543887102634
TOL = 5e-13

# GateCDiag pre-v11 authoritative records (endgame_gatecdiag_final.json,
# sha256 f4699bbd..., frozen in the GateCDiag campaign + scratch twin).
GC = json.load(open("/Users/jinleic/jinleic-workspace/cs/omega/scratch/"
                    "endgame_gatecdiag_final.json"))
PRE_V11 = GC["V1_frozen_machine_equality"]["real"]
PRE_V11_RAW = GC["V1_point_value"]["Om_raw"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


SRC = "/Users/jinleic/jinleic-workspace/cs/omega/src"
provenance = {
    "interval_core": sha256(f"{SRC}/interval_core.py"),
    "vxxz24_float": sha256(f"{SRC}/vxxz24_float.py"),
    "gate_c_slope_core": sha256(f"{SRC}/gate_c_slope_core.py"),
    "gate_c_slope_pass": sha256(f"{SRC}/gate_c_slope_pass.py"),
    "gate_c_endgame": sha256(f"{SRC}/gate_c_endgame.py"),
}
ckpt("provenance", provenance)

assert provenance["vxxz24_float"] == \
    "31657c92a1f841805e7355d5b5390f76e0d7ca2d04b1a3f2ce9fc204d9d5d890"
assert provenance["interval_core"] == \
    "ce70f95922f3d5098e892e1b7047570ebb5267826b51f2776c91b56d9934adb5", \
    "live core is not v11 — environment regression before any run"
assert provenance["gate_c_slope_pass"] == \
    "c6b8cef29716c2157660fad4b85b4259ccf940f0664f17f0f2991152964f0dfc"
assert provenance["gate_c_endgame"] == \
    "baaecb8b2743c1a6cacde67a972a269d5593b3f5573fe0cd498c8cf52b369093"

# ================= STEP 1 — frozen-protocol control on the v11 core =======
import interval_core as IC
import vxxz24_float as F
from scipy.io import loadmat
from flint import arb

params = np.asarray(loadmat(
    "/Users/jinleic/jinleic-workspace/cs/omega/scratch/rmmcode/data/"
    "K100_2.37155181.mat")["params"]).flatten()

ctrl = {}
with IC.prec():
    # frozen stage_b_rung2.run(include_lemma1=True) — run in-process,
    # unmodified byte source executed via importlib against live deps.
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "stage_b_rung2_frozen",
        "/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
        "2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/stage_b_rung2.py")
    sb2 = importlib.util.module_from_spec(spec)
    t1 = time.time()
    ctrl["with"] = spec.loader.exec_module(sb2) or sb2.run(include_lemma1=True)
    ctrl["with_elapsed"] = time.time() - t1
    ckpt("ctrl_with_v11", {
        "raw": ctrl["with"]["omega_cert_upper_raw"],
        "R_sum": ctrl["with"]["R_sum_low_certified"],
        "M_low": ctrl["with"]["M_low_certified"],
        "R_detail": ctrl["with"]["R_detail"],
        "elapsed": ctrl["with_elapsed"],
        "evidence": "COMPUTATIONAL-EVIDENCE (fresh run under v11 core)",
    })
    t2 = time.time()
    ctrl["without"] = sb2.run(include_lemma1=False)
    ctrl["without_elapsed"] = time.time() - t2
    ckpt("ctrl_without_v11", {
        "raw": ctrl["without"]["omega_cert_upper_raw"],
        "R_sum": ctrl["without"]["R_sum_low_certified"],
        "R_detail": ctrl["without"]["R_detail"],
        "elapsed": ctrl["without_elapsed"],
        "evidence": "COMPUTATIONAL-EVIDENCE",
    })

ctrl_with_raw = ctrl["with"]["omega_cert_upper_raw"]
ctrl_delta = ctrl_with_raw - FROZEN_RAW
c_gate = abs(ctrl_delta) <= TOL
ckpt("C_gate", {"ctrl_with_raw": ctrl_with_raw, "delta": ctrl_delta,
                "tol": TOL, "ok": bool(c_gate)})
if not c_gate:
    ckpt("decision", {
        "branch": "Z",
        "reason": ("frozen-protocol control drifted outside 5e-13 under "
                   "the v11 core; environment/regression finding, "
                   "pre-statement Step-1 stop rule"),
    })
    print("BRANCH Z — control fail; stopping per pre-statement")
    sys.exit(3)

# ================= STEP 2 — replay V1 on the same process/core ============
import gate_c_slope_pass as SP

with IC.prec():
    out0 = SP.run_slope_pass(params, radius=0.0, verbose=False)
    # Byte-for-byte REUSE of the as-run GateCDiag aggregate(): extract
    # the aggregate+r_detail function block from the asrun source at
    # runtime and exec it in a namespace holding the same live modules.
    asrun_path = ("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
                  "2026-08-31T07:57:00Z_OmegaGateCDiag_4859327/"
                  "gate_c_endgame.py.asrun")
    aesrc = open(asrun_path, "r").read()
    a0 = aesrc.index("    def aggregate(out, radius, lemma1_off=False):")
    a1 = aesrc.index("    def r_detail(agg):")
    block = aesrc[a0:a1]
    # dedent one level: inside the file these defs live in a with-block;
    # strip exactly 4 leading spaces from every line (blank lines pass).
    lines = []
    for ln in block.splitlines():
        lines.append(ln[4:] if ln.startswith("    ") else ln)
    agg_code = "\n".join(lines)
    import gate_c_slope_core as SC
    import vxxz24_float as F

    def _slope_of(x):
        return x.v if isinstance(x, SP.SlopeG) else x

    ns = {"SC": SC, "SP": SP, "IC": IC, "np": np, "arb": arb,
          "I": IC.Ival, "F": F, "slope_of": _slope_of}
    exec(compile(agg_code, "<aggregate-asrun>", "exec"), ns)
    aggregate = ns["aggregate"]

    agg0 = aggregate(out0, 0.0)
    agg0_off = aggregate(out0, 0.0, lemma1_off=True)
    replay_raw = agg0["Om_raw"]

branch = agg0["branch"]
det = {}
for k, v in (("R_comp[0,3]", branch[(3, 0)]), ("R_comp[1,3]", branch[(3, 1)]),
             ("R_comp[2,3]", branch[(3, 2)]), ("R_comp[0,2]", branch["lvl2"]),
             ("R_glob[0]", branch[("G", 0)]), ("R_glob[1]", branch[("G", 1)]),
             ("R_glob[2]", branch[("G", 2)])):
    det[k] = max(0.0, min(v["cand"]))

v1_delta = replay_raw - FROZEN_RAW
v1_ok = abs(v1_delta) <= TOL
branches_equal = det == PRE_V11["R_detail"]
r_equal = agg0["total_R"] == PRE_V11["R_sum"]
m_equal = agg0["M_low"] == PRE_V11["M_low"]
raw_equal_prev11 = replay_raw == PRE_V11_RAW

ckpt("V1_v11", {
    "replay_raw": replay_raw,
    "frozen": FROZEN_RAW,
    "delta": v1_delta,
    "ok": bool(v1_ok),
    "branches_float_equal_pre_v11_records": bool(branches_equal),
    "R_float_equal": bool(r_equal),
    "M_float_equal": bool(m_equal),
    "raw_float_equal_pre_v11": bool(raw_equal_prev11),
    "R_detail": det,
    "R_sum": agg0["total_R"],
    "M_low": agg0["M_low"],
    "lemma1_off_raw": agg0_off["Om_raw"],
    "lemma1_off_R_sum": agg0_off["total_R"],
})
if not (v1_ok and branches_equal):
    # ---- branch B: pre-registered single diagnostic (no source edit) ----
    gm = out0["gm"]
    eps_rows = []
    for r in range(3):
        dm = [x.val for x in gm.dist_max[r].v]
        rows = []
        for i, shp in enumerate(gm.shapes):
            if float(dm[i].hi) <= 0:
                continue
            g_val = gm.lam_sum[r].v[0].val.lo - 1.0
            for d in range(3):
                g_val += gm.lam_margin[(r, d)].v[shp[d]].val.lo
            diff = IC.ln_iv(dm[i]) - IC.const(float(g_val))
            e = diff.pos()
            rows.append({"shape": list(shp), "eps_hi": float(e.hi),
                         "dm_hi": float(dm[i].hi)})
        rows.sort(key=lambda z: -z["eps_hi"])
        eps_rows.append(rows[:3])
    ckpt("V1_fail_diagnostic", {
        "branch": "B",
        "note": ("v11 replay divergence; per-region top-eps shapes vs "
                 "frozen lemma1_eps_by_group (glob) — see manifest"),
        "top_eps_glob": eps_rows,
        "frozen_lemma1_eps_max": 3.076684958606817e-06,
    })
    ckpt("decision", {"branch": "B",
                      "reason": "V1 or branch-equality failed under v11"})
    print("BRANCH B — V1 mismatch; diagnostic recorded; stop per tree")
    sys.exit(4)

# ================= STEP 3 — LP/remainder stages (V1 passed) ==============
RADIUS = 1e-7
outB = SP.run_slope_pass(params, radius=RADIUS, verbose=False)
d0B = outB["pm"].start[outB["gm"].dist_id[0]]
radius_iv = IC.const(RADIUS)
leafB = outB["gm"].dist[0].v
old_inward = {"lower": [], "upper": []}
# Endpoint identity is asserted under the SAME prec envelope the leaves
# were constructed in (arb equality compares balls; a different ambient
# workprec rebuilds the ball radius and breaks object equality without
# any endpoint disagreement — diagnostics-only probe recorded this).
with IC.prec():
    for i, leaf in enumerate(leafB):
        x = float(params[d0B + i])
        center_iv = IC.const(x)
        lo_true = (center_iv - radius_iv).lo
        hi_true = (center_iv + radius_iv).hi
        assert leaf.val.lo == lo_true, (i, float(lo_true), float(leaf.val.lo))
        assert leaf.val.hi == hi_true, (i, float(hi_true), float(leaf.val.hi))
        if IC.const(x - RADIUS).lo > lo_true:
            old_inward["lower"].append(i)
        if IC.const(x + RADIUS).hi < hi_true:
            old_inward["upper"].append(i)
ckpt("box_endpoints", {
    "mode": "exact Arb IC.const(center) +/- IC.const(radius)",
    "old_binary64_counterfactual_inward": old_inward,
    "evidence": "MACHINE-VERIFIED endpoint assertions, all 45 leaves",
})
aggB = aggregate(outB, RADIUS)
Om_gradB = aggB["Om_slope"].grad
lo_b = np.array([float(x.lo) for x in Om_gradB])
hi_b = np.array([float(x.hi) for x in Om_gradB])
wid_b = hi_b - lo_b
ckpt("slope_box", {
    "width_max": float(wid_b.max()),
    "width_sum_l1": float(wid_b.sum()),
    "grad_lo": lo_b.tolist(),
    "grad_hi": hi_b.tolist(),
    "elapsed": time.time() - T0,
    "evidence": "MACHINE-VERIFIED (interval AD, v11 core)",
})

# scaled midpoint LP (candidate only)
from scipy.optimize import linprog
Vex = json.load(open("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
                     "2026-08-30T11:24:41Z_17bdaeb6_099b1224d603/"
                     "kernel_basis_V_exact_rational.json"))
V = np.column_stack([[num / den for num, den in Vex[str(j)]]
                     for j in range(21)])
A = np.load("/Users/jinleic/jinleic-workspace/cs/omega/campaigns/"
            "2026-08-30T11:24:41Z_17bdaeb6_099b1224d603/margin_matrix_M.npy")
assert V.shape == (45, 21) and A.shape == (27, 45)
A_int = A.astype(np.int64)
V_int = V.astype(np.int64)
assert np.array_equal(A, A_int) and np.array_equal(V, V_int)
assert np.array_equal(A_int @ V_int, np.zeros((27, 21), dtype=np.int64))
assert np.array_equal(V_int.sum(axis=0), np.zeros(21, dtype=np.int64))

grad_mid = (lo_b + hi_b) / 2.0
gv_mid = V.T @ grad_mid
Aub = np.vstack([V, -V])
rLPs = linprog(gv_mid, A_ub=Aub, b_ub=np.ones(90),
               bounds=[(None, None)] * 21, method="highs")
assert rLPs.status == 0
lin_mid_scaled = RADIUS * float(rLPs.fun)
ckpt("LP_scaled_midpoint", {
    "min": lin_mid_scaled, "status": int(rLPs.status),
    "evidence": "COMPUTATIONAL-EVIDENCE (float candidate, never certified)",
})

# rigorous ANY-y instrument
def float_up(x):
    return float(np.nextafter(float(x), np.inf))


lo_fit = np.nextafter(lo_b, -np.inf)
hi_fit = np.nextafter(hi_b, np.inf)
fit_obj = np.concatenate([np.zeros(27), np.ones(45)])
fit_A = np.vstack([
    np.hstack([A.T, -np.eye(45)]),
    np.hstack([-A.T, -np.eye(45)]),
])
fit_b = np.concatenate([lo_fit, -hi_fit])
fit = linprog(fit_obj, A_ub=fit_A, b_ub=fit_b,
              bounds=[(None, None)] * 27 + [(0, None)] * 45, method="highs")
assert fit.status == 0
y = np.asarray(fit.x[:27])
I = IC.Ival
z_iv = I(arb(0), arb(0))
residual = z_iv
for i, g_iv in enumerate(Om_gradB):
    aty = z_iv
    for j in np.flatnonzero(A[:, i]):
        aty = aty + IC.const(float(y[j]))
    residual = residual + (g_iv - aty).pos()
residual_zero_y = z_iv
for i, g_iv in enumerate(Om_gradB):
    residual_zero_y = residual_zero_y + g_iv.pos()
swing_iv = IC.const(RADIUS) * residual
swing_zero_iv = IC.const(RADIUS) * residual_zero_y
swing_upper = float_up(swing_iv.up())
swing_zero_upper = float_up(swing_zero_iv.up())
assert swing_upper < swing_zero_upper
lin_min = -swing_upper
lin_max = swing_upper
width_cert = float_up(2.0 * swing_upper)
SIGNAL = 1.5816497000997742e-07
ratio = width_cert / SIGNAL
if lin_min >= 0.0:
    decision = ("PASS-a: no feasible improving direction in D certified "
                "by the ANY-y slope enclosure")
elif lin_max < 0.0:
    decision = ("PASS-b candidate: enclosure strictly below Omega0 over "
                "D — exact witness extraction required before any claim")
else:
    decision = ("FAILURE TO CERTIFY: ANY-y slope enclosure straddles zero "
                "over D; 21-dimensional question remains OPEN")
ckpt("ANY_y_v11", {
    "evidence": "MACHINE-VERIFIED",
    "instrument": "R*sup_G ||G-A^T y||_1; arbitrary y",
    "fit_objective_float": float(fit.fun),
    "residual_l1_hi_buffered": float_up(residual.up()),
    "zero_y_swing_upper": swing_zero_upper,
    "swing_upper": swing_upper,
    "lin_min": lin_min, "lin_max": lin_max,
    "width": width_cert,
})
signal_check = {"signal": SIGNAL, "ratio": ratio,
                "predecessor_ratio_122x": 122.0,
                "gatecdiag_ratio_73.86x": 73.85977508485647}
ckpt("signal_ratio", signal_check)
ckpt("decision", {
    "branch": "A",
    "lin_min": lin_min, "lin_max": lin_max, "width": width_cert,
    "ratio": ratio, "Om0_raw": replay_raw,
    "midpoint_LP_candidate": lin_mid_scaled,
    "decision": decision,
    "scope": ("D caps 12.8x short of published 2.37155181; route can "
              "never produce a record; 21-dim kernel question UNCHANGED "
              "and OPEN; gate C PARTIAL regardless"),
    "elapsed_total": time.time() - T0,
})
print("DONE", decision)
