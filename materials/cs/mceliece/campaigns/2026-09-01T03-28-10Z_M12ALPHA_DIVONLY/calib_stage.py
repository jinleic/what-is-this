"""Calibration stage — M12ALPHA-DIVONLY campaign, pre-verdict.

Measures, IN-PROCESS (time.process_time / perf_counter only; never ps):
  1. full unmodified guard-chain run at (11,2048,48,6211) — the anchor
     re-time (sanity vs frozen 285.91 s; this includes the FULL alpha
     grid at m=11 scale);
  2. build CPU+wall for (12,3488,64,16384);
  3. per-full-check alpha cost at m=12 scale (theInstanceAlphaOrGrid —
     the alpha block of instance.py exactly as in src, timed over a few
     support points and one per-coordinate set) and the full-grid
     projection, via the same code path.

Writes calib.json in the campaign dir. NO verdict computed here.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mceliece/src")

from instance import Instance  # noqa: E402

OUT = ("/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/"
       "2026-09-01T03-28-10Z_M12ALPHA_DIVONLY")

rec = {"tool": "calib", "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def cpu_s():
    return time.process_time()


# ---------------- 1. anchor re-time at m=11 (FULL guard chain, unmodified)
t0w, t0c = time.perf_counter(), cpu_s()
inst11 = Instance(11, 2048, 48, 6211)   # includes full guards incl alpha grid
w11, c11 = time.perf_counter() - t0w, cpu_s() - t0c
rec["anchor_m11"] = {
    "m": 11, "n": 2048, "t": 48, "seed": 6211,
    "guards": {k: v for k, v in inst11.guards.items() if k != "degs"},
    "wall_s": round(w11, 2), "cpu_s": round(c11, 2),
    "frozen_elapsed_s_ref": 285.91,
}
print("m11 anchor:", rec["anchor_m11"], flush=True)

# ---------------- 2. build at m=12
t0w, t0c = time.perf_counter(), cpu_s()
inst12 = Instance(12, 3488, 64, 16384)
w12b, c12b = time.perf_counter() - t0w, cpu_s() - t0c
rec["build_m12"] = {
    "m": 12, "n": 3488, "t": 64, "seed": 16384,
    "k": inst12.k, "D": inst12.D,
    "guards_partial": {k: v for k, v in inst12.guards.items()
                       if k in ("gamma_k", "eps_gcd_const", "eps_maxdeg_is_D",
                                 "beta_delta_nonzero", "beta_pair")},
    "wall_s": round(w12b, 2), "cpu_s": round(c12b, 2),
}
print("m12 build:", rec["build_m12"], flush=True)

# ---------------- 3. alpha-grid cost probe at m=12 scale (UNMODIFIED path)
# Reproduce EXACTLY the alpha block's inner loop from instance.py
# (vectorized Fa/Fpa over all k coordinates at one support point a), timed.
from fastfield import EField  # noqa: E402
import census as _c  # noqa: E402  (lucas_w lives here)

ef = inst12.ef
gf = inst12.gf
D = inst12.D
Dp1 = D + 1
k = inst12.k

# coefficient matrices exactly as in _guards
Dmat = np.zeros((k, Dp1), dtype=np.uint16)
Cmat2 = np.zeros((k, Dp1), dtype=np.uint16)
for j, f in enumerate(inst12.F):
    cf = np.asarray(f + [0] * max(0, Dp1 - len(f)), dtype=np.uint16)[:Dp1]
    Cmat2[j] = cf
    fp = np.zeros(Dp1, dtype=np.uint16)
    for e in range(1, len(cf)):
        if e % 2 == 1:
            fp[e - 1] = cf[e]
    Dmat[j] = fp

Pi_np = np.asarray(inst12.Pi, dtype=np.uint16)
PiD_np = np.asarray(inst12.PiD, dtype=np.uint16)
Gn = np.asarray(inst12.G, dtype=np.uint16)

nz = np.nonzero(Cmat2.any(axis=0))[0]
print(f"m12 coeff matrix nonempty degree slots: {nz.size}", flush=True)


def alpha_full_check_at(i):
    """one support point i: vectorized over all k coords (unmodified inner op)"""
    a = inst12.support[i]
    w0 = _c.lucas_w(ef, Dp1 - 1, 0, a)
    Fa = np.bitwise_xor.reduce(ef.MUL[Cmat2[:, nz], w0[nz][None, :]], axis=1)
    Fpa = np.bitwise_xor.reduce(ef.MUL[Dmat[:, nz], w0[nz][None, :]], axis=1)
    wp = _c.lucas_w(ef, len(Pi_np) - 1, 0, a)
    pia = int(np.bitwise_xor.reduce(ef.MUL[Pi_np, wp]))
    w1 = _c.lucas_w(ef, len(PiD_np) - 1, 0, a)
    pida = int(np.bitwise_xor.reduce(ef.MUL[PiD_np, w1]))
    Ga_ = int(np.bitwise_xor.reduce(ef.MUL[Gn, _c.lucas_w(ef, len(Gn) - 1, 0, a)]))
    lhs = ef.MUL[pia, Fpa] ^ ef.MUL[pida, Fa]
    rhs = ef.MUL[ef.MUL[Ga_, Ga_], ef.MUL[Fa, Fa]]
    return bool(np.array_equal(lhs, rhs))


probe_pts = [0, 1, 2, 3, 4]
t0c = cpu_s()
bools = [alpha_full_check_at(i) for i in probe_pts]
dtc = cpu_s() - t0c
per_pt = dtc / len(probe_pts)
rec["alpha_cost_m12"] = {
    "probe_points": probe_pts,
    "per_point_cpu_s": round(per_pt, 4),
    "all_true_on_probes": bool(all(bools)),
    "projected_full_grid_cpu_s": round(per_pt * 3488, 1),
}
print("m12 alpha cost:", rec["alpha_cost_m12"], flush=True)

rec["budget_projection"] = {
    "budget_cpu_s": 10800,
    "projected_alpha_grid_cpu_s": rec["alpha_cost_m12"]["projected_full_grid_cpu_s"],
    "m11_grid_cpu_s_reference": rec["anchor_m11"]["cpu_s"],
}
rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
json.dump(rec, open(os.path.join(OUT, "calib.json"), "w"), indent=1, default=str)
print("calib written", flush=True)
