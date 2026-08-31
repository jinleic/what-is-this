"""Gate C — stage (a) FLOAT branch-and-bound (foundation for the Arb leg).

Per the mandated ordering (README two-stage protocol + my pre_statement
§4): this screen runs FIRST and must be recorded CLEAN before the Arb leg
starts. It is float64-only — no interval arithmetic anywhere.

Pre-declared probe set (from pre_statement §3, axes selected from the
program's owned group structure, NOT adaptive):
  * p* itself (anchor).
  * ± all-coordinates and ± per-interesting-group (region_prop, omega,
    single_mat_size, glob dist r=0,1,2) at the pre-registered radius
    r = 1e-7.
Each probe records omega_float = (target − R_sum)/M evaluated in float64.
Nothing here is certified; it calibrates the enclosure for stage (b).
"""
import sys, json, time
import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/omega/src")
import vxxz24_float as F
from scipy.io import loadmat

MAT = ("/Users/jinleic/jinleic-workspace/cs/omega/scratch/rmmcode/"
       "data/K100_2.37155181.mat")
RADIUS = 1e-7
TARGET = 4 * np.log(7.0)
CERT_PSTAR = 2.3715538358544617  # frozen rung-2 with-arm endpoint
Q = 5.0
L = 3


def build():
    pm = F.ParamManager()
    ws = F.Workspace(pm, Q, 1.0, 0.0, L)
    params = np.asarray(loadmat(MAT)["params"]).flatten()
    pm.set_value(params)
    return pm, ws, params


def _scalar(z):
    if isinstance(z, float):
        return z
    if hasattr(z, "v"):
        z = z.v
    return float(np.asarray(z).reshape(-1)[0])


def eval_R_M(pm, ws, x):
    pm.set_value(x)
    c, ce, value = ws.evaluate()
    gm = ws.globstage
    total_R = 0.0
    # level-3 component groups
    for l in range(3, L + 1):
        for r in range(3):
            num_block = [0.0] * 3
            pen = 0.0
            pcomp = 0.0
            for t in ws.parts[l - 1]:
                nb = t.num_block_contribution[r]
                num_block = [a + b for a, b in zip(num_block, nb)]
                pen = pen + _scalar(t.hash_penalty_term[r])
                pcomp = pcomp + _scalar(t.p_comp[r])
            for tt in range(3):
                cand = num_block[tt] - (pcomp if tt == r else pen)
                if cand > 0:
                    total_R += cand
    # level-2 symmetric
    num_block = [0.0] * 3
    for t in ws.parts[1]:
        num_block = [a + b for a, b in zip(num_block, t.num_block_contribution)]
    if min(num_block) > 0:
        total_R += min(num_block)
    # global
    for r in range(3):
        num_block = [_scalar(z) for z in gm.num_block[r]]
        pen = _scalar(gm.hash_penalty_term[r])
        pcomp = _scalar(gm.p_comp[r])
        for tt in range(3):
            cand = num_block[tt] - (pcomp if tt == r else pen)
            if cand > 0:
                total_R += cand
    M = min(_scalar(z) for z in gm.mat_size)
    cfe = max(_scalar(z) for z in c)
    cem = max(abs(_scalar(z)) for z in ce)
    return total_R, M, cfe, cem


if __name__ == "__main__":
    t0 = time.time()
    pm, ws, params = build()
    n = len(params)
    R0, M0, cfe0, cem0 = eval_R_M(pm, ws, params)
    om0 = (TARGET - R0) / M0
    print("p*_screen", json.dumps({"omega_float": om0, "R": R0, "M": M0,
                                   "cF_max": cfe0, "ceF_max": cem0}))
    gm = ws.globstage
    interesting = {
        "region_prop": gm.region_prop_id,
        "omega": ws.omega_id,
        "single_mat_size": ws.single_mat_size_id,
        "glob_dist_0": gm.dist_id[0],
        "glob_dist_1": gm.dist_id[1],
        "glob_dist_2": gm.dist_id[2],
        "glob_dmax_0": gm.dist_max_id[0],
        "glob_dmax_1": gm.dist_max_id[1],
        "glob_dmax_2": gm.dist_max_id[2],
    }
    probes = [("all+", np.ones(n)), ("all-", -np.ones(n))]
    for name, gid in interesting.items():
        s, z = pm.start[gid], pm.size[gid]
        d = np.zeros(n)
        d[s:s + z] = 1.0
        probes.append((name + "+", d))
        probes.append((name + "-", -d))
    results = []
    for name, d in probes:
        Rc, Mc, cfe, cem = eval_R_M(pm, ws, params + RADIUS * d)
        om = (TARGET - Rc) / Mc
        results.append({"probe": name, "omega_float": om,
                        "delta_omega": om - om0,
                        "max_c": cfe, "max_ceq": cem})
        print(json.dumps(results[-1]))
        sys.stdout.flush()
    out = {"screen": "stage_a_float", "radius": RADIUS,
           "omega_pstar_float": om0, "target_4ln7": TARGET,
           "R_pstar": R0, "M_pstar": M0,
           "cert_pstar": CERT_PSTAR, "probes": results,
           "wallclock_s": time.time() - t0}
    with open(sys.argv[1] if len(sys.argv) > 1 else "stage_a_screen.json",
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("stage (a) screen complete")
