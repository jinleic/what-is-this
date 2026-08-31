"""Square-subsystem refinement of Route-F rank-13 numerical candidates.

At each round, pivoted QR selects 192 independent columns from the full
192x247 CP Jacobian.  The other 55 coordinates are frozen and dense LM solves
the resulting square system.  A candidate is still only numerical until the
same selected subsystem passes route_af_certify.py in exact rational
arithmetic.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.linalg import qr
from scipy.optimize import least_squares

from route_af_search import (T, TF, M, RANK, balance, factor_stats, jacobian,
                             residual_for)

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parent / "campaigns" / "2026-08-31T08:02:18Z_routeAF"
ROUNDS = int(os.environ.get("ROUNDS", "4"))
MAX_NFEV = int(os.environ.get("MAX_NFEV", "1500"))
TRIGGER = 1e-11


def select_columns(x: np.ndarray) -> tuple[np.ndarray, dict]:
    j = jacobian(x)
    _, rmat, piv = qr(j, mode="economic", pivoting=True)
    selected = np.asarray(piv[:M], dtype=np.int64)
    js = j[:, selected]
    sv = np.linalg.svd(js, compute_uv=False)
    info = {"smallest_singular": float(sv[-1]),
            "largest_singular": float(sv[0]),
            "condition": float(sv[0] / sv[-1]),
            "rank_at_1e-10": int(np.linalg.matrix_rank(js, tol=1e-10))}
    if info["rank_at_1e-10"] != M:
        raise RuntimeError(f"selected subsystem singular: {info}")
    return selected, info


def refine(target: np.ndarray, x0: np.ndarray, rounds: int = ROUNDS,
           max_nfev: int = MAX_NFEV) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    x = balance(np.asarray(x0, dtype=float))
    trace: list[dict] = []
    last_selected = None
    for round_no in range(1, rounds + 1):
        selected, sinfo = select_columns(x)
        fixed = x.copy()

        def embed(y: np.ndarray) -> np.ndarray:
            z = fixed.copy()
            z[selected] = y
            return z

        def fun(y: np.ndarray) -> np.ndarray:
            return residual_for(target, embed(y), RANK)

        def jac(y: np.ndarray) -> np.ndarray:
            return jacobian(embed(y))[:, selected]

        rel0 = float(np.linalg.norm(fun(x[selected])) / np.linalg.norm(target))
        result = least_squares(fun, x[selected], jac=jac, method="lm",
                               max_nfev=max_nfev, ftol=3e-15, xtol=3e-15,
                               gtol=3e-15, x_scale="jac")
        x = embed(np.asarray(result.x, dtype=float))
        rel = float(np.linalg.norm(residual_for(target, x, RANK)) /
                    np.linalg.norm(target))
        row = {"round": round_no, "rel_before": rel0, "rel": rel,
               "nfev": int(result.nfev), "status": int(result.status),
               "message": str(result.message), "selection": sinfo,
               **factor_stats(x)}
        trace.append(row)
        print("AF_REFINE_ROUND " + json.dumps(row), flush=True)
        last_selected = selected
        if rel < TRIGGER:
            break
        x = balance(x)
    assert last_selected is not None
    return x, last_selected, trace


def planted_control() -> dict:
    rng = np.random.default_rng(2026083101)
    truth = balance(rng.normal(size=(3 + 8 + 8) * RANK))
    a = truth[:3 * RANK].reshape(3, RANK)
    b = truth[3 * RANK:11 * RANK].reshape(8, RANK)
    c = truth[11 * RANK:].reshape(8, RANK)
    target = np.einsum("pr,br,cr->pbc", a, b, c)
    start = truth + 1e-5 * rng.normal(size=truth.shape)
    got, selected, trace = refine(target, start, rounds=2, max_nfev=300)
    rel = float(np.linalg.norm(residual_for(target, got, RANK)) /
                np.linalg.norm(target))
    if rel >= TRIGGER:
        raise AssertionError(f"square refiner did not catch plant: {rel}")
    return {"caught": True, "rel": rel, "rounds": len(trace),
            "selected_count": int(len(selected))}


def main() -> None:
    print("AF_REFINE_ANCHOR " + json.dumps({"planted_rank13": planted_control()}))
    starts = [
        ("random_best_seed2", CAMPAIGN / "rank13_candidate.npz"),
        ("extension_best_seed2", CAMPAIGN / "rank13_extension_candidate.npz"),
    ]
    best = None
    for label, path in starts:
        if not path.exists():
            continue
        data = np.load(path)
        x0 = np.asarray(data["x"], dtype=float)
        print("AF_REFINE_START " + json.dumps({
            "label": label, "path": str(path),
            "rel": float(np.linalg.norm(residual_for(T, x0, RANK)) / TF),
            **factor_stats(x0),
        }), flush=True)
        x, selected, trace = refine(T, x0)
        rel = float(np.linalg.norm(residual_for(T, x, RANK)) / TF)
        out = CAMPAIGN / f"rank13_refined_{label}.npz"
        np.savez(out, x=x, selected=selected, target=T, rel=rel,
                 trace=json.dumps(trace),
                 convention="T[p,b,c]=sum_s A[p,s]B[b,s]C[c,s]; globally conjugated Route-F target tau direct-sum tau")
        row = {"label": label, "path": str(out), "rel": rel,
               "selected_count": int(len(selected)), **factor_stats(x)}
        print("AF_REFINE_RESULT " + json.dumps(row), flush=True)
        if best is None or rel < best[0]:
            best = (rel, row)
    if best is None:
        raise FileNotFoundError("no numerical candidates")
    print("AF_REFINE_BEST " + json.dumps(best[1]))
    if best[0] >= TRIGGER:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
