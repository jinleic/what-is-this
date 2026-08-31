"""Gate B bounded rank-ladder search with border-rank diagnostics.

This is the actual pre-statement engine: CP-ALS (<=2000 sweeps) followed by
Levenberg-Marquardt/Gauss-Newton (<=300 rounds), PCG64 seeded starts, ranks
24,23,22,21,20.  A nonzero residual is never a lower-bound certificate.
A residual tending to zero with unbounded factors is only border-rank
behaviour.  Every seed record includes residual, factor-column norms, and
maxnorm/residual at checkpoints and at its best point.
"""
import json
import os
import time

import numpy as np

BASE_SEED = 20260829
STARTS = int(os.environ.get("STARTS", "256"))
RANKS = [int(v) for v in os.environ.get("RANKS", "24").split(",")]
ALS_MAX = int(os.environ.get("ALS_MAX", "2000"))
LM_MAX = int(os.environ.get("LM_MAX", "300"))
CERT_TOL = float(os.environ.get("CERT_TOL", "1e-13"))
WCAP_H = float(os.environ.get("WCAP_H", "8"))


def octonion_tensor():
    def t4(i, j):
        if i == 0:
            return j, 1
        if j == 0:
            return i, 1
        if i == j:
            return 0, -1
        return 6 - i - j, (1 if (i, j) in ((1, 2), (2, 3), (3, 1)) else -1)

    def hmul(u, v):
        out = [0] * 4
        for i in range(4):
            for j in range(4):
                if u[i] and v[j]:
                    k, sg = t4(i, j)
                    out[k] += sg * u[i] * v[j]
        return out

    def omul(x, y):
        a, b, c, d = x[:4], x[4:], y[:4], y[4:]
        dc = [d[0], -d[1], -d[2], -d[3]]
        cc = [c[0], -c[1], -c[2], -c[3]]
        return ([hmul(a, c)[k] - hmul(dc, b)[k] for k in range(4)]
                + [hmul(d, a)[k] + hmul(b, cc)[k] for k in range(4)])

    E = np.eye(8, dtype=int)
    return np.asarray([[omul(E[i], E[j]) for j in range(8)] for i in range(8)],
                       dtype=float)


T = octonion_tensor()
TF = float(np.linalg.norm(T))
N1 = N2 = N3 = 8
M = N1 * N2 * N3

# Fixed flattening e=(i,j,k), identical to the certification scripts.
I_IDX = np.repeat(np.arange(N1), N2 * N3)
J_IDX = np.tile(np.repeat(np.arange(N2), N3), N1)
K_IDX = np.tile(np.arange(N3), N1 * N2)


def factor_stats(x, r, residual):
    oB, oC = N1 * r, (N1 + N2) * r
    aa = x[:oB].reshape(N1, r)
    bb = x[oB:oC].reshape(N2, r)
    cc = x[oC:].reshape(N3, r)
    norms = (np.linalg.norm(aa, axis=0),
             np.linalg.norm(bb, axis=0),
             np.linalg.norm(cc, axis=0))
    per = (float(norms[0].max()), float(norms[1].max()),
           float(norms[2].max()))
    mx = max(per)
    prod = float(np.max(norms[0] * norms[1] * norms[2]))
    return mx, per, mx / max(residual, 1e-300), prod


def residual(x, r):
    oB, oC = N1 * r, (N1 + N2) * r
    aa = x[:oB].reshape(N1, r)
    bb = x[oB:oC].reshape(N2, r)
    cc = x[oC:].reshape(N3, r)
    return (np.einsum("is,js,ks->ijk", aa, bb, cc) - T).ravel()


def jacobian(x, r):
    oB, oC = N1 * r, (N1 + N2) * r
    aa = x[:oB].reshape(N1, r)
    bb = x[oB:oC].reshape(N2, r)
    cc = x[oC:].reshape(N3, r)
    j = np.empty((M, 3 * N1 * r), dtype=float)
    bc = bb[J_IDX, :] * cc[K_IDX, :]
    ac = aa[I_IDX, :] * cc[K_IDX, :]
    ab = aa[I_IDX, :] * bb[J_IDX, :]
    for i in range(N1):
        sel = I_IDX == i
        j[np.ix_(sel, np.arange(i * r, (i + 1) * r))] = bc[sel]
    for q in range(N2):
        sel = J_IDX == q
        cols = oB + np.arange(q * r, (q + 1) * r)
        j[np.ix_(sel, cols)] = ac[sel]
    for k in range(N3):
        sel = K_IDX == k
        cols = oC + np.arange(k * r, (k + 1) * r)
        j[np.ix_(sel, cols)] = ab[sel]
    return j


def als(x, r, max_sweeps, rng):
    """CP-ALS with a small ridge for numerical stability."""
    oB, oC = N1 * r, (N1 + N2) * r
    aa = x[:oB].reshape(N1, r).copy()
    bb = x[oB:oC].reshape(N2, r).copy()
    cc = x[oC:].reshape(N3, r).copy()
    eye = np.eye(r)
    traj = []
    last = np.inf
    for it in range(1, max_sweeps + 1):
        g = (bb.T @ bb) * (cc.T @ cc)
        rhs = np.einsum("ijk,js,ks->si", T, bb, cc)
        aa = np.linalg.solve(g + 1e-9 * eye, rhs).T
        g = (aa.T @ aa) * (cc.T @ cc)
        rhs = np.einsum("ijk,is,ks->sj", T, aa, cc)
        bb = np.linalg.solve(g + 1e-9 * eye, rhs).T
        g = (aa.T @ aa) * (bb.T @ bb)
        rhs = np.einsum("ijk,is,js->sk", T, aa, bb)
        cc = np.linalg.solve(g + 1e-9 * eye, rhs).T
        xx = np.concatenate([aa.ravel(), bb.ravel(), cc.ravel()])
        f = residual(xx, r)
        rel = float(np.linalg.norm(f) / TF)
        if it in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2000):
            mx, _, rat, prod = factor_stats(xx, r, rel)
            traj.append({"phase": "ALS", "iter": it, "rel": rel,
                         "max_factor_norm": mx, "max_column_product": prod,
                         "ratio": rat})
        if rel <= CERT_TOL:
            return xx, rel, it, traj
        # If the update stalls for a long time, stop ALS and let LM take over;
        # this is an algorithmic stop, not a mathematical conclusion.
        if it > 50 and abs(last - rel) < 1e-13:
            return xx, rel, it, traj
        last = rel
    return np.concatenate([aa.ravel(), bb.ravel(), cc.ravel()]), rel, it, traj


def lm(x, r, max_rounds):
    """Analytic-Jacobian trust-region LM/Gauss-Newton refinement.

    The residual has 512 equations and 3*8*r unknowns, so scipy's trust
    region solver with ``tr_solver='lsmr'`` handles the underdetermined
    normal equations without a singular dense inverse.  The callback records
    factor norms at every accepted iterate; a near-zero residual with norms
    diverging is explicitly classified as border-rank evidence, never as a
    rank certificate.
    """
    from scipy.optimize import least_squares

    trace = []

    def callback(xk, *args):
        f = residual(xk, r)
        rel = float(np.linalg.norm(f) / TF)
        mx, per, ratio, prod = factor_stats(xk, r, rel)
        trace.append({"phase": "LM", "iter": len(trace) + 1, "rel": rel,
                      "max_factor_norm": mx, "factor_norms": per,
                      "max_column_product": prod, "ratio": ratio})

    res = least_squares(lambda z: residual(z, r), x,
                        jac=lambda z: jacobian(z, r), method="trf",
                        tr_solver="lsmr", x_scale="jac",
                        max_nfev=max_rounds,
                        xtol=3e-16, ftol=3e-16, gtol=3e-16,
                        callback=callback)
    xx = np.asarray(res.x, dtype=float)
    ff = residual(xx, r)
    rel = float(np.linalg.norm(ff) / TF)
    # SciPy can terminate before callback on a final trial; include that
    # endpoint so the saved record always pairs its residual with norms.
    if not trace or trace[-1]["rel"] != rel:
        mx, per, ratio, prod = factor_stats(xx, r, rel)
        trace.append({"phase": "LM-final", "iter": len(trace) + 1,
                      "rel": rel, "max_factor_norm": mx,
                      "factor_norms": per, "max_column_product": prod,
                      "ratio": ratio})
    return xx, rel, int(res.nfev), trace, int(res.status)


def main():
    summary = {}
    for r in RANKS:
        t0 = time.time()
        rng = np.random.default_rng([BASE_SEED, r])
        best = None
        rows = []
        for st in range(STARTS):
            scale = (TF / r) ** (1.0 / 3.0)
            seed = [BASE_SEED, r, st]
            rng_start = np.random.default_rng(seed)
            x0 = rng_start.standard_normal(3 * 8 * r) * scale
            xa, ra, na, ta = als(x0, r, ALS_MAX, rng_start)
            if ra > CERT_TOL:
                xl, rl, nl, tl, lm_status = lm(xa, r, LM_MAX)
            else:
                xl, rl, nl, tl, lm_status = xa, ra, 0, [], 1
            mx, per, ratio, prod = factor_stats(xl, r, rl)
            row = {"start": st, "seed": [BASE_SEED, r, st],
                   "als_rel": ra, "als_sweeps": na,
                   "lm_rel": rl, "lm_rounds": nl,
                   "lm_status": lm_status, "best_rel": rl,
                   "max_factor_norm": mx, "factor_norms": per,
                   "max_column_product": prod,
                   "ratio_norm_over_resid": ratio,
                   "als_trajectory": ta, "lm_trajectory": tl}
            rows.append(row)
            if best is None or rl < best[0]:
                best = (rl, xl.copy(), mx, per, ratio, prod)
            if rl <= CERT_TOL:
                # save this as a candidate; do NOT claim it until Krawczyk
                # certification passes, and stop this rank's declared starts
                # only on the explicit success criterion.
                break
            if time.time() - t0 > WCAP_H * 3600:
                rows.append({"wallclock_cutoff": True, "at_start": st})
                break
        rel, xb, mx, per, ratio, prod = best
        np.savez(f"gate_b_rank{r}_best.npz", x=xb, rel=rel, r=r,
                 max_factor_norm=mx, max_column_product=prod)
        rec = {"rank": r, "starts_attempted": len(rows),
               "best_rel": rel, "best_max_factor_norm": mx,
               "best_factor_norms": per, "best_max_column_product": prod,
               "best_ratio": ratio, "seed": [BASE_SEED, r],
               "elapsed_s": time.time() - t0, "rows": rows}
        with open(f"gate_b_rank{r}_results.json", "w") as f:
            json.dump(rec, f, indent=1, default=float)
        summary[str(r)] = {k: rec[k] for k in
                           ("rank", "starts_attempted", "best_rel",
                            "best_max_factor_norm", "best_factor_norms",
                            "best_max_column_product", "best_ratio",
                            "seed", "elapsed_s")}
        print(f"RANK {r}: best rel={rel:.3e}; maxfactor={mx:.3e}; "
              f"ratio={ratio:.3e}; starts={len(rows)}; "
              f"elapsed={rec['elapsed_s']:.1f}s", flush=True)
    print("B_SEARCH3_JSON " + json.dumps(summary, default=float))


if __name__ == "__main__":
    main()
