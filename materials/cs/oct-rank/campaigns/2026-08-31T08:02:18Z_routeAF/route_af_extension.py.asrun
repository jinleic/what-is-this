"""Real-diagonal commuting-extension search for Route-F rank 13.

For slices (I,A,B), a rank-r decomposition with nonzero first coordinates is
P Q = I, P diag(d) Q = A, P diag(e) Q = B.  Write P=E S and
Q=S^{-1}E^T with real invertible S; then a solution yields the real CP factors
([1,d,e], P, Q^T).  This script is numerical discovery only.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from route_af_search import T, TF, P as NSLICES, B as N, RANK, balance, factor_stats

R = RANK
K = R - N
NV = R * R + 2 * R
NSEEDS = int(os.environ.get("NSEEDS", "10"))
MAX_NFEV = int(os.environ.get("MAX_NFEV", "5000"))
BASE_SEED = 20260931
TRIGGER = 1e-11
HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parent / "campaigns" / "2026-08-31T08:02:18Z_routeAF"
OUT = Path(os.environ.get("OUT", CAMPAIGN / "rank13_extension_candidate.npz"))
A_TARGET = T[1]
B_TARGET = T[2]
E = np.zeros((N, R))
E[:, :N] = np.eye(N)


def unpack(z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return z[:R * R].reshape(R, R), z[R * R:R * R + R], z[-R:]


def compression(s: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q = np.linalg.solve(s, E.T)
    p = s[:N, :]
    return p @ (d[:, None] * q), p, q


def residual(z: np.ndarray, targets: tuple[np.ndarray, np.ndarray] = (A_TARGET, B_TARGET)) -> np.ndarray:
    s, d, e = unpack(z)
    try:
        md, _, _ = compression(s, d)
        me, _, _ = compression(s, e)
    except np.linalg.LinAlgError:
        return np.full(2 * N * N, 1e12)
    if not np.all(np.isfinite(md)) or not np.all(np.isfinite(me)):
        return np.full(2 * N * N, 1e12)
    return np.concatenate([(md - targets[0]).ravel(), (me - targets[1]).ravel()])


def jacobian(z: np.ndarray) -> np.ndarray:
    s, d, e = unpack(z)
    q = np.linalg.solve(s, E.T)
    p = s[:N, :]
    sinv = np.linalg.inv(s)
    j = np.zeros((2 * N * N, NV), dtype=float)
    for u in range(R):
        for v in range(R):
            col = u * R + v
            wd = p @ (d[:, None] * sinv[:, u:u + 1])
            we = p @ (e[:, None] * sinv[:, u:u + 1])
            left_d = -wd[:, 0]
            left_e = -we[:, 0]
            if u < N:
                left_d[u] += d[v]
                left_e[u] += e[v]
            j[:N * N, col] = np.outer(left_d, q[v]).ravel()
            j[N * N:, col] = np.outer(left_e, q[v]).ravel()
    od = R * R
    for v in range(R):
        outer = np.outer(p[:, v], q[v]).ravel()
        j[:N * N, od + v] = outer
        j[N * N:, od + R + v] = outer
    return j


def cp_factors(z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    s, d, e = unpack(z)
    q = np.linalg.solve(s, E.T)
    p = s[:N, :]
    af = np.vstack([np.ones(R), d, e])
    bf = p
    cf = q.T
    return af, bf, cf


def cp_vector(z: np.ndarray) -> np.ndarray:
    a, b, c = cp_factors(z)
    return balance(np.concatenate([a.ravel(), b.ravel(), c.ravel()]))


def row_for(z: np.ndarray, result, seed: int, rel0: float) -> dict:
    x = cp_vector(z)
    rr = residual(z)
    rel = float(np.linalg.norm(rr) / TF)
    s, _, _ = unpack(z)
    rank_j = int(np.linalg.matrix_rank(jacobian(z), tol=1e-9))
    return {"seed": seed, "rel_before": rel0, "rel": rel,
            "nfev": int(result.nfev), "status": int(result.status),
            "message": str(result.message), "cond_S": float(np.linalg.cond(s)),
            "jacobian_rank_at_1e-9": rank_j, **factor_stats(x)}


def planted_control() -> dict:
    rng = np.random.default_rng(BASE_SEED - 1)
    s0 = rng.normal(size=(R, R))
    s0 += 3.0 * np.eye(R)
    d0 = rng.normal(size=R)
    e0 = rng.normal(size=R)
    z0 = np.concatenate([s0.ravel(), d0, e0])
    md, _, _ = compression(s0, d0)
    me, _, _ = compression(s0, e0)
    zstart = z0 + 1e-5 * rng.normal(size=NV)
    result = least_squares(lambda z: residual(z, (md, me)), zstart,
                           jac=lambda z: jacobian(z), method="trf",
                           tr_solver="lsmr", x_scale="jac", max_nfev=300,
                           ftol=3e-15, xtol=3e-15, gtol=3e-15)
    rel = float(np.linalg.norm(residual(result.x, (md, me))) /
                np.linalg.norm(np.concatenate([md.ravel(), me.ravel()])))
    if rel >= TRIGGER:
        raise AssertionError(f"commuting-extension plant not recovered: {rel}")
    return {"caught": True, "rel": rel, "nfev": int(result.nfev)}


def start(seed: int) -> np.ndarray:
    rng = np.random.default_rng(BASE_SEED + seed)
    # A nonnormal but controlled real basis is necessary: an orthogonal basis
    # can only compress real diagonals to symmetric matrices, whereas targets
    # A and B are skew-symmetric.
    s = rng.normal(scale=0.45, size=(R, R)) + np.eye(R)
    while np.linalg.cond(s) > 1e4:
        s = rng.normal(scale=0.45, size=(R, R)) + np.eye(R)
    d = rng.normal(scale=1.0, size=R)
    e = rng.normal(scale=1.0, size=R)
    return np.concatenate([s.ravel(), d, e])


def main() -> None:
    # Analytic Jacobian startup anchor against central differences.
    zcheck = start(0)
    direction = np.random.default_rng(BASE_SEED).normal(size=NV)
    h = 1e-6
    fd = (residual(zcheck + h * direction) - residual(zcheck - h * direction)) / (2 * h)
    jd = jacobian(zcheck) @ direction
    jac_rel = float(np.linalg.norm(fd - jd) / max(1.0, np.linalg.norm(fd)))
    if jac_rel >= 1e-7:
        raise AssertionError(f"extension Jacobian anchor failed: {jac_rel}")
    print("AF_EXTENSION_ANCHOR " + json.dumps({
        "n": N, "r": R, "k": K, "jacobian_directional_rel": jac_rel,
        "planted_rank13": planted_control(),
    }))

    best: tuple[float, np.ndarray, dict] | None = None
    for seed in range(1, NSEEDS + 1):
        z0 = start(seed)
        rel0 = float(np.linalg.norm(residual(z0)) / TF)
        result = least_squares(residual, z0, jac=jacobian, method="trf",
                               tr_solver="lsmr", x_scale="jac",
                               max_nfev=MAX_NFEV, ftol=3e-15, xtol=3e-15,
                               gtol=3e-15)
        z = np.asarray(result.x, dtype=float)
        row = row_for(z, result, seed, rel0)
        print("AF_EXTENSION_ROW " + json.dumps(row), flush=True)
        if best is None or row["rel"] < best[0]:
            best = (row["rel"], z, row)
        if row["rel"] < TRIGGER:
            break
    assert best is not None
    z = best[1]
    x = cp_vector(z)
    af, bf, cf = cp_factors(z)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, z=z, x=x, A=af, B=bf, C=cf, rel=best[0],
             row=json.dumps(best[2]), target=T,
             convention="T[p,b,c]=sum_s A[p,s]B[b,s]C[c,s]; conjugated Route-F target tau direct-sum tau")
    print("AF_EXTENSION_BEST " + json.dumps({"path": str(OUT), **best[2]}))
    if best[0] >= TRIGGER:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
