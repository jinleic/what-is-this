"""gate_b_search.py — rank-24 (then lower) numerical search for a CP
decomposition of the octonion multiplication tensor T_O, per the
pre-declared budget in pre_statement.md:

* ranks r = 24 first (then 23, 22, 21, 20 if time allows);
* 256 seeded starts per rank (PCG64, seed stream [BASE_SEED + r, start]);
* 2000 ALS sweeps cap per start, then up to 300 Levenberg–Marquardt rounds;
* stop-early on relative residual <= 1e-14;
* wall-clock cap per rank: 8 h;
* SUCCESS-CANDIDATE for certification: rel-F residual <= 1e-13.

A FAILED search yields NO claim about rank. Only Krawczyk-certified
existence counts. Output: per-rank JSON + best factors as .npz + a final
B_SEARCH_JSON line on stdout.
"""
import json
import os
import sys
import time

import numpy as np

BASE_SEED = 20260829
STARTS = int(os.environ.get("STARTS", "256"))
RANKS = [int(v) for v in os.environ.get("RANKS", "24").split(",")]
ALS_MAX = int(os.environ.get("ALS_MAX", "2000"))
GN_MAX = int(os.environ.get("GN_MAX", "300"))
REL_TOL = 1e-14
CERT_TOL = 1e-13
WCAP = float(os.environ.get("WCAP_H", "8")) * 3600.0

np.random.seed(0)


def octonion_tensor_local():
    def t4(i, j):
        if i == 0:
            return j, 1
        if j == 0:
            return i, 1
        if i == j:
            return 0, -1
        k = 6 - i - j
        return k, (1 if (i, j) in ((1, 2), (2, 3), (3, 1)) else -1)

    def hv4(u, v):
        out = [0] * 4
        for i in range(4):
            if u[i] == 0:
                continue
            for j in range(4):
                if v[j] == 0:
                    continue
                k, s = t4(i, j)
                out[k] += s * u[i] * v[j]
        return out

    def h_mul(x, y):
        a, b = x[:4], x[4:]
        c, d = y[:4], y[4:]
        dc = [d[0], -d[1], -d[2], -d[3]]
        cc = [c[0], -c[1], -c[2], -c[3]]
        left = [hv4(a, c)[k] - hv4(dc, b)[k] for k in range(4)]
        right = [hv4(d, a)[k] + hv4(b, cc)[k] for k in range(4)]
        return left + right

    E8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
    return np.array([[[h_mul(E8[i], E8[j])[k] for k in range(8)]
                      for j in range(8)] for i in range(8)], dtype=float)


T = octonion_tensor_local()
TF = float(np.linalg.norm(T))


def cp_als(Tt, r, rng, max_sw, rel_tol, ridge=1e-8):
    n1, n2, n3 = Tt.shape
    scale = (np.linalg.norm(Tt) / r) ** 0.5
    A = rng.standard_normal((n1, r)) * scale
    B = rng.standard_normal((n2, r)) * scale
    Cc = rng.standard_normal((n3, r)) * scale
    I = np.eye(r)

    def resid(Av, Bv, Cv):
        R = Tt - np.einsum('is,js,ks->ijk', Av, Bv, Cv)
        return float(np.linalg.norm(R))

    res = resid(A, B, Cc)
    div = False
    sweeps = 0
    for sw in range(max_sw):
        sweeps = sw + 1
        G = (B.T @ B) * (Cc.T @ Cc)
        M = np.einsum('ijk,js,ks->sk', Tt, B, Cc)   # (r, n1)
        A = np.linalg.solve(G + ridge * I, M)       # (r, n1), rows = A[:,s]
        A = A.T
        G = (A.T @ A) * (Cc.T @ Cc)
        M = np.einsum('ijk,is,ks->sk', Tt, A, Cc)   # (r, n2)
        B = np.linalg.solve(G + ridge * I, M).T     # (n2, r)
        G = (A.T @ A) * (B.T @ B)
        M = np.einsum('ijk,is,js->sk', Tt, A, B)    # (r, n3)
        Cc = np.linalg.solve(G + ridge * I, M).T    # (n3, r)
        new = resid(A, B, Cc)
        if not np.isfinite(new) or new > 1e6 * max(res, 1.0):
            div = True
            break
        res = new
        if res / TF <= rel_tol:
            break
    return A, B, Cc, res / TF, sweeps, div


def build_index(n1, n2, n3):
    i_idx = np.repeat(np.arange(n1), n2 * n3)
    j_idx = np.tile(np.repeat(np.arange(n2), n3), n1)
    k_idx = np.tile(np.arange(n3), n1 * n2)
    return i_idx, j_idx, k_idx


def gn_lm(Tt, A, B, Cc, max_rounds, rel_tol):
    """Gauss-Newton with Levenberg damping on the full residual vector.

    Exact trilinear Jacobian, assembled column-block-wise:
      J[e, i*r + s] = B[j,s] C[k,s],
      J[e, n1*r + j*r + s] = A[i,s] C[k,s],
      J[e, (n1+n2)*r + k*r + s] = A[i,s] B[j,s].
    """
    n1, n2, n3 = Tt.shape
    r = A.shape[1]
    N = (n1 + n2 + n3) * r
    offB, offC = n1 * r, (n1 + n2) * r
    M = n1 * n2 * n3
    i_idx, j_idx, k_idx = build_index(n1, n2, n3)
    x = np.concatenate([A.ravel(), B.ravel(), Cc.ravel()])
    lam = 1e-10

    def unpack(xv):
        return (xv[:offB].reshape(n1, r),
                xv[offB:offC].reshape(n2, r),
                xv[offC:].reshape(n3, r))

    def F(xv):
        Av, Bv, Cv = unpack(xv)
        return (np.einsum('is,js,ks->ijk', Av, Bv, Cv) - Tt).ravel()

    def jac(xv):
        Av, Bv, Cv = unpack(xv)
        J = np.empty((M, N))
        BA = np.empty((M, r))
        BB = np.empty((M, r))
        BC = np.empty((M, r))
        for s in range(r):
            BA[:, s] = Bv[j_idx, s] * Cv[k_idx, s]
            BB[:, s] = Av[i_idx, s] * Cv[k_idx, s]
            BC[:, s] = Av[i_idx, s] * Bv[j_idx, s]
        for i in range(n1):
            sel = i_idx == i
            J[np.ix_(sel, np.arange(i * r, (i + 1) * r))] = BA[sel, :]
        for j in range(n2):
            sel = j_idx == j
            J[np.ix_(sel, offB + np.arange(j * r, (j + 1) * r))] = BB[sel, :]
        for k in range(n3):
            sel = k_idx == k
            J[np.ix_(sel, offC + np.arange(k * r, (k + 1) * r))] = BC[sel, :]
        return J

    f = F(x)
    res = float(np.linalg.norm(f) / TF)
    it = 0
    for it in range(1, max_rounds + 1):
        Jm = jac(x)
        g = Jm.T @ f
        JJ = Jm.T @ Jm
        step_ok = False
        for _ in range(30):
            try:
                dx = np.linalg.solve(JJ + lam * np.eye(N), -g)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            xnew = x + dx
            fnew = F(xnew)
            resnew = float(np.linalg.norm(fnew) / TF)
            if resnew < res:
                x, f, res = xnew, fnew, resnew
                lam = max(lam * 0.3, 1e-16)
                step_ok = True
                break
            lam *= 5
            if lam > 1e8:
                break
        if not step_ok or res <= rel_tol or not np.isfinite(res):
            break
    return unpack(x), res, it


if __name__ == "__main__":
    out = {}
    for r in RANKS:
        t_r = time.time()
        best = (None, np.inf)
        results = []
        for st in range(STARTS):
            rng = np.random.default_rng([BASE_SEED + r, st])
            A, B, Cc, rel, sweeps, div = cp_als(T, r, rng, ALS_MAX, REL_TOL)
            if div or not np.isfinite(rel):
                results.append({"start": st, "diverged": True})
                continue
            rec = {"start": st, "als_rel": rel, "als_sweeps": sweeps}
            if rel < 1e-6:
                (A2, B2, C2), rel2, gnits = gn_lm(T, A, B, Cc, GN_MAX, REL_TOL)
                rec.update({"gn_rel": rel2, "gn_rounds": gnits})
            else:
                rel2 = rel
                rec.update({"gn_rel": None})
            results.append(rec)
            if best[0] is None or rel2 < best[1]:
                best = ((A.copy(), B.copy(), Cc.copy()), rel2)
            if best[1] <= REL_TOL:
                break
            if time.time() - t_r > WCAP:
                results.append({"start_cutoff_wallclock": True,
                                "at_start": st})
                break
        rec_rank = {
            "rank": r,
            "n_starts": len([x for x in results if "diverged" in x or
                             "als_rel" in x]),
            "best_rel": float(best[1]) if best[0] is not None else None,
            "seed_base": BASE_SEED + r,
            "elapsed_s": time.time() - t_r,
            "results": results,
        }
        out[str(r)] = rec_rank
        if best[0] is not None:
            A, B, Cc = best[0]
            np.savez(f"gate_b_rank{r}_best.npz", A=A, B=B, C=Cc,
                     rel=best[1])
        with open(f"gate_b_rank{r}_results.json", "w") as f:
            json.dump(rec_rank, f, indent=1, default=float)
        print(f"RANK {r}: best_rel = {rec_rank['best_rel']:.3e} "
              f"({rec_rank['elapsed_s']:.1f}s, "
              f"{rec_rank['n_starts']} starts)", flush=True)
    print("B_SEARCH_JSON " + json.dumps(
        {k: {"best_rel": v["best_rel"], "starts": v["n_starts"],
             "seed_base": v["seed_base"], "elapsed_s": v["elapsed_s"]}
         for k, v in out.items()}, default=float))
