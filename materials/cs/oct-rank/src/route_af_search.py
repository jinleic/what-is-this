"""Deterministic rank-13 search for the Route-F tensor (L_1,L_i,L_j).

Discovery only: floating-point outputs are COMPUTATIONAL-EVIDENCE until
route_af_certify.py proves an exact Krawczyk enclosure.  Convention:
Cayley--Dickson basis (1,i,j,k,l,il,jl,kl); tensor entry T[p,b,c] is the
coefficient of e_c in e_p*e_b for p in (0,1,2).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

RANK = 13
P, B, C = 3, 8, 8
M = P * B * C
BASE_SEED = 20260831
MAX_NFEV = int(os.environ.get("MAX_NFEV", "5000"))
N_RANDOM = int(os.environ.get("N_RANDOM", "5"))
CERT_TRIGGER = 1e-11
HERE = Path(__file__).resolve().parent
TARGET = HERE.parent / "campaigns" / "2026-08-31T08:02:18Z_routeAF"
OUT = Path(os.environ.get("OUT", TARGET / "rank13_candidate.npz"))


def qmul(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.array([
        x[0] * y[0] - x[1] * y[1] - x[2] * y[2] - x[3] * y[3],
        x[0] * y[1] + x[1] * y[0] + x[2] * y[3] - x[3] * y[2],
        x[0] * y[2] - x[1] * y[3] + x[2] * y[0] + x[3] * y[1],
        x[0] * y[3] + x[1] * y[2] - x[2] * y[1] + x[3] * y[0],
    ], dtype=float)


def build_target() -> np.ndarray:
    e4 = np.eye(4)
    tau = np.empty((3, 4, 4), dtype=float)
    for p in range(3):
        for b in range(4):
            tau[p, b, :] = qmul(e4[p], e4[b])
    target = np.zeros((3, 8, 8), dtype=float)
    target[:, :4, :4] = tau
    target[:, 4:, 4:] = tau
    assert np.array_equal(target[0], np.eye(8))
    assert np.count_nonzero(target) == 24
    assert all(np.linalg.matrix_rank(target.reshape(P, -1)) == P for _ in [0])
    assert np.linalg.matrix_rank(np.moveaxis(target, 1, 0).reshape(B, -1)) == B
    assert np.linalg.matrix_rank(np.moveaxis(target, 2, 0).reshape(C, -1)) == C
    return target


T = build_target()
TF = float(np.linalg.norm(T))
PI = np.repeat(np.arange(P), B * C)
BI = np.tile(np.repeat(np.arange(B), C), P)
CI = np.tile(np.arange(C), P * B)


def split(x: np.ndarray, rank: int = RANK) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ob = P * rank
    oc = (P + B) * rank
    return x[:ob].reshape(P, rank), x[ob:oc].reshape(B, rank), x[oc:].reshape(C, rank)


def join(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return np.concatenate([a.ravel(), b.ravel(), c.ravel()])


def residual_for(target: np.ndarray, x: np.ndarray, rank: int = RANK) -> np.ndarray:
    a, b, c = split(x, rank)
    return (np.einsum("pr,br,cr->pbc", a, b, c) - target).ravel()


def residual(x: np.ndarray) -> np.ndarray:
    return residual_for(T, x)


def jacobian(x: np.ndarray, rank: int = RANK) -> np.ndarray:
    a, b, c = split(x, rank)
    nv = (P + B + C) * rank
    j = np.zeros((M, nv), dtype=float)
    ob = P * rank
    oc = (P + B) * rank
    bc = b[BI, :] * c[CI, :]
    ac = a[PI, :] * c[CI, :]
    ab = a[PI, :] * b[BI, :]
    for p in range(P):
        rows = PI == p
        j[np.ix_(rows, np.arange(p * rank, (p + 1) * rank))] = bc[rows]
    for q in range(B):
        rows = BI == q
        cols = ob + np.arange(q * rank, (q + 1) * rank)
        j[np.ix_(rows, cols)] = ac[rows]
    for q in range(C):
        rows = CI == q
        cols = oc + np.arange(q * rank, (q + 1) * rank)
        j[np.ix_(rows, cols)] = ab[rows]
    return j


def balance(x: np.ndarray, rank: int = RANK) -> np.ndarray:
    a, b, c = (z.copy() for z in split(x, rank))
    for s in range(rank):
        na, nb, nc = (np.linalg.norm(a[:, s]), np.linalg.norm(b[:, s]),
                      np.linalg.norm(c[:, s]))
        if min(na, nb, nc) == 0:
            continue
        g = (na * nb * nc) ** (1.0 / 3.0)
        a[:, s] *= g / na
        b[:, s] *= g / nb
        c[:, s] *= g / nc
    return join(a, b, c)


def factor_stats(x: np.ndarray, rank: int = RANK) -> dict[str, float]:
    a, b, c = split(x, rank)
    per = [float(np.linalg.norm(z[:, s])) for s in range(rank) for z in (a, b, c)]
    products = [float(np.linalg.norm(a[:, s]) * np.linalg.norm(b[:, s]) *
                      np.linalg.norm(c[:, s])) for s in range(rank)]
    return {"max_factor_norm": max(per), "max_column_product": max(products),
            "min_column_product": min(products)}


def als(target: np.ndarray, x: np.ndarray, rank: int, sweeps: int = 250) -> np.ndarray:
    a, b, c = (z.copy() for z in split(x, rank))
    eye = np.eye(rank)
    for _ in range(sweeps):
        g = (b.T @ b) * (c.T @ c)
        rhs = np.einsum("pbc,br,cr->rp", target, b, c)
        a = np.linalg.solve(g + 1e-10 * eye, rhs).T
        g = (a.T @ a) * (c.T @ c)
        rhs = np.einsum("pbc,pr,cr->rb", target, a, c)
        b = np.linalg.solve(g + 1e-10 * eye, rhs).T
        g = (a.T @ a) * (b.T @ b)
        rhs = np.einsum("pbc,pr,br->rc", target, a, b)
        c = np.linalg.solve(g + 1e-10 * eye, rhs).T
        x = balance(join(a, b, c), rank)
        a, b, c = split(x, rank)
    return join(a, b, c)


def lm(target: np.ndarray, x: np.ndarray, rank: int, max_nfev: int) -> tuple[np.ndarray, dict]:
    if target.shape != (P, B, C):
        raise ValueError(target.shape)
    fun = lambda z: residual_for(target, z, rank)
    jac = (lambda z: jacobian(z, rank)) if rank == RANK else "2-point"
    res = least_squares(fun, balance(x, rank), jac=jac, method="trf",
                        tr_solver="lsmr", x_scale="jac", max_nfev=max_nfev,
                        ftol=3e-15, xtol=3e-15, gtol=3e-15)
    xx = balance(np.asarray(res.x, dtype=float), rank)
    rel = float(np.linalg.norm(fun(xx)) / np.linalg.norm(target))
    row = {"rel": rel, "nfev": int(res.nfev), "status": int(res.status),
           "message": str(res.message), **factor_stats(xx, rank)}
    return xx, row


def planted_control() -> dict:
    rng = np.random.default_rng(BASE_SEED - 1)
    a = rng.normal(size=(P, RANK))
    b = rng.normal(size=(B, RANK))
    c = rng.normal(size=(C, RANK))
    truth = balance(join(a, b, c))
    aa, bb, cc = split(truth)
    plant = np.einsum("pr,br,cr->pbc", aa, bb, cc)
    start = truth + 1e-4 * rng.normal(size=truth.shape)
    got, row = lm(plant, start, RANK, 400)
    row["caught"] = bool(row["rel"] < CERT_TRIGGER)
    if not row["caught"]:
        raise AssertionError(f"rank-13 plant not recovered: {row}")
    # Retain `got` to prevent an accidental no-op control from passing merely
    # because the unperturbed truth was evaluated.
    row["distance_from_start"] = float(np.linalg.norm(got - start))
    return row


def random_start(seed: int) -> np.ndarray:
    rng = np.random.default_rng(BASE_SEED + seed)
    x = rng.normal(scale=0.7, size=(P + B + C) * RANK)
    return als(T, x, RANK)


def main() -> None:
    print("AF_SEARCH_ANCHOR " + json.dumps({
        "tensor_shape": list(T.shape), "target_norm": TF,
        "nonzeros": int(np.count_nonzero(T)),
        "flattening_ranks": [
            int(np.linalg.matrix_rank(T.reshape(P, -1))),
            int(np.linalg.matrix_rank(np.moveaxis(T, 1, 0).reshape(B, -1))),
            int(np.linalg.matrix_rank(np.moveaxis(T, 2, 0).reshape(C, -1))),
        ],
        "planted_rank13": planted_control(),
    }))
    best: tuple[float, np.ndarray, dict] | None = None
    for seed in range(1, N_RANDOM + 1):
        start = random_start(seed)
        rel0 = float(np.linalg.norm(residual(start)) / TF)
        x, row = lm(T, start, RANK, MAX_NFEV)
        row.update({"kind": "random+ALS", "seed": seed, "rel_before_lm": rel0})
        print("AF_SEARCH_ROW " + json.dumps(row), flush=True)
        if best is None or row["rel"] < best[0]:
            best = (row["rel"], x, row)
        if row["rel"] < CERT_TRIGGER:
            break
    assert best is not None
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, x=best[1], rel=best[0], row=json.dumps(best[2]),
             target=T, convention="T[p,b,c]=coefficient e_c in e_p*e_b; p=0,1,2; globally conjugated to tau direct-sum tau")
    print("AF_SEARCH_BEST " + json.dumps({"path": str(OUT), **best[2]}))
    if best[0] >= CERT_TRIGGER:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
