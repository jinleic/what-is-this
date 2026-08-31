"""gate_b_search2.py — rank ladder with border-rank instrumentation.

Main's directive (2026-08-30): track factor-norm growth alongside residual
at every iteration. A residual -> 0 with diverging factor norms is a
border-rank witness, NOT a rank certificate.

Per (rank, seed) record:
  - best residual (rel F-norm);
  - factor norms at the best point: max_s ||A[:,s]||, ||B[:,s]||, ||C[:,s]||;
  - ratio maxfactor / residual (grows like a power of 1/residual in the
    degenerate case; stays O(1) for a genuine decomposition);
  - iterates at checkpoint iterations: (residual, maxnorm, ratio) triples.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, "..")
sys.path.insert(0, "../scratch/upstream_ref/verify")
from octonion_core import octonion_tensor

TO = octonion_tensor().astype(float)
TF = float(np.linalg.norm(TO))
BASE_SEED = 20260829
STARTS = int(os.environ.get("STARTS", "256"))
RANKS = [int(v) for v in os.environ.get("RANKS", "24").split(",")]
MAX_NFEV = int(os.environ.get("MAX_NFEV", "4000"))
CERT_TOL = 1e-13
WCAP = float(os.environ.get("WCAP_H", "6")) * 3600.0

CHK = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]


def factor_norms(x, r):
    offB = 8 * r
    A = x[:offB].reshape(8, r)
    B = x[offB:2 * offB].reshape(8, r)
    C = x[2 * offB:].reshape(8, r)
    nA = float(np.linalg.norm(A, axis=0).max())
    nB = float(np.linalg.norm(B, axis=0).max())
    nC = float(np.linalg.norm(C, axis=0).max())
    return max(nA, nB, nC), (nA, nB, nC)


def resid_flat(x, Tt, r):
    offB = 8 * r
    A = x[:offB].reshape(8, r)
    B = x[offB:2 * offB].reshape(8, r)
    C = x[2 * offB:].reshape(8, r)
    return (np.einsum('is,js,ks->ijk', A, B, C) - Tt).ravel()


out = {}
for r in RANKS:
    t_r = time.time()
    rng = np.random.default_rng([BASE_SEED, r])
    best = None
    results = []
    for st in range(STARTS):
        scale = (TF / r) ** 0.5
        x0 = rng.standard_normal(3 * 8 * r) * scale

        traj = []

        def cb(xk, rec=traj, it=[0]):
            it[0] += 1
            if it[0] - 1 in CHK or it[0] in CHK:
                f = resid_flat(xk, TO, r)
                v = float(np.linalg.norm(f)) / TF
                mx, _ = factor_norms(xk, r)
                rec.append((it[0], v, mx, mx / max(v, 1e-300)))
            # keep the last iterate always
            rec.append((it[0], None, None, None)) if False else None

        res = least_squares(resid_flat, x0, args=(TO, r), method='trf',
                            xtol=3e-16, ftol=3e-16, gtol=3e-16,
                            max_nfev=MAX_NFEV, callback=cb)
        v = float(np.linalg.norm(res.fun)) / TF
        mx, per = factor_norms(res.x, r)
        rec = {"start": st,
               "best_rel": v,
               "max_factor_norm": mx,
               "factor_norms": list(per),
               "ratio_norm_over_resid": mx / max(v, 1e-300),
               "nfev": int(res.nfev),
               "trajectory": traj}
        results.append(rec)
        if best is None or v < best[0]:
            best = (v, res.x.copy(), mx)
        if v <= CERT_TOL:
            break
        if time.time() - t_r > WCAP:
            results.append({"start_cutoff_wallclock": True})
            break
    v, xbest, mxb = best
    np.savez(f"/tmp/rank{r}_cand.npz", x=xbest, rel=v, maxnorm=mxb,
             r=r)
    rec_rank = {"rank": r, "n_starts": len(results), "best_rel": v,
                "best_max_factor_norm": mxb,
                "ratio": mxb / max(v, 1e-300),
                "seed_base": BASE_SEED, "elapsed_s": time.time() - t_r,
                "results": results}
    with open(f"/tmp/rank{r}_results.json", "w") as f:
        json.dump(rec_rank, f, indent=1, default=float)
    print(f"RANK {r}: best_rel={v:.3e} max_factor_norm={mxb:.3e} "
          f"ratio={mxb/max(v,1e-300):.3e} ({rec_rank['elapsed_s']:.0f}s, "
          f"{len(results)} starts)", flush=True)
    out[str(r)] = {k: rec_rank[k] for k in
                   ("rank", "n_starts", "best_rel", "best_max_factor_norm",
                    "ratio", "seed_base", "elapsed_s")}
print("B_SEARCH2_JSON " + json.dumps(out, default=float))
