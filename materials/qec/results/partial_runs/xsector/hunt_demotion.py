"""Small-lattice hunt: does d_X(Q) < d_X(P) ever occur for a valid PBB?

Mechanism (Theorem J.4): X-logicals of Q split into
  survivors  x in Xcen(P)\\S_X(P) with [C D]x^T = 0   (weight >= d_X(P)), and
  demoted    x in rowspace(H_X) with [C D]x^T = 0 but x not in S_X(Q).
Hence d_X(Q) < d_X(P) iff some demoted element has weight < d_X(P).
We enumerate small BB parents and small valid perturbations and hunt for
exactly that event.  GF(2)/numpy only, one thread, no SAT.
"""
from __future__ import annotations

import itertools
import json
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import monomial_matrix, poly_matrix  # noqa: E402
from qec_research.gf2.linalg import rank_np, nullspace_np, matmul as gf2matmul  # noqa: E402


class MonoCache:
    def __init__(self, ell, m):
        self.ell, self.m = ell, m
        self.cache = {}

    def mono(self, terms):
        key = tuple(sorted(terms))
        if key not in self.cache:
            self.cache[key] = poly_matrix(self.ell, self.m, list(key))
        return self.cache[key]


def min_wt_coset(basis, exclude, cap_bits=22):
    """Min wt over span(basis) minus span(exclude); None if basis too big."""
    d = basis.shape[0]
    if d > cap_bits:
        return None
    if exclude is not None and exclude.shape[0]:
        exn = nullspace_np(exclude)
    else:
        exn = None
    best = None
    for start in range(1, 1 << d, 1 << 18):
        idx = np.arange(start, min(start + (1 << 18), 1 << d), dtype=np.int64)
        V = np.zeros((idx.size, basis.shape[1]), dtype=np.uint8)
        for j in range(d):
            V ^= (((idx >> j) & 1)[:, None] & basis[j]).astype(np.uint8)
        if exn is not None and exn.shape[0]:
            V = V[~(((V @ exn.T) % 2).sum(axis=1) == 0)]
        if V.shape[0]:
            b = int(V.sum(axis=1).min())
            best = b if best is None else min(best, b)
    return best


def dX_support(HZ, SXP_perp, n, cap):
    """Min wt x with H_Z x^T = 0 and x not in rowspace(H_X); cap-limited."""
    for w in range(1, cap + 1):
        combos = itertools.combinations(range(n), w)
        while True:
            chunk = list(itertools.islice(combos, 100_000))
            if not chunk:
                break
            X = np.zeros((len(chunk), n), dtype=np.uint8)
            for r, sup in enumerate(chunk):
                X[r, list(sup)] = 1
            ok = ~(((X @ HZ.T) % 2).any(axis=1))
            if ok.any():
                Xin = X[ok]
                if SXP_perp.shape[0]:
                    Xin = Xin[~(((Xin @ SXP_perp.T) % 2).sum(axis=1) == 0)]
                if Xin.shape[0]:
                    return w, Xin[0]
    return None, None


def span_outside(kerM, HX, SXQ, N):
    """Min weight of {lam[A B] : lam in span(kerM), lam != 0} minus span(SXQ)."""
    d = kerM.shape[0]
    if d > 16:
        return None, None
    ns = nullspace_np(SXQ) if SXQ.shape[0] else None
    best, bestx = None, None
    for start in range(1, 1 << d, 1 << 14):
        idx = np.arange(start, min(start + (1 << 14), 1 << d), dtype=np.int64)
        V = np.zeros((idx.size, kerM.shape[1]), dtype=np.uint8)
        for j in range(d):
            V ^= (((idx >> j) & 1)[:, None] & kerM[j]).astype(np.uint8)
        X = (V @ HX) % 2
        w = X.sum(axis=1)
        cand = np.argsort(w, kind="stable")
        for r in cand:
            wr = int(w[r])
            if wr == 0 or (best is not None and wr >= best):
                break
            if ns is not None and ns.shape[0]:
                if (((X[r] @ ns.T) % 2).sum()) == 0:
                    continue  # inside S_X(Q)
            best, bestx = wr, X[r]
            break
    return best, (bestx if bestx is None else bestx[None, :])
def perturbations(ell, m, max_mons=3):
    """Anchored (C_terms, D_terms) with 1..max_mons total monomials."""
    grid = [(a, b) for a in range(ell) for b in range(m)]
    seen = set()
    out = []
    for total in range(1, max_mons + 1):
        for split in range(total + 1):
            for Cs in itertools.combinations(grid, split):
                for Ds in itertools.combinations(grid, total - split):
                    allm = sorted(Cs + Ds)
                    t = allm[0]
                    def shift(s):
                        return tuple(sorted(((a - t[0]) % ell, (b - t[1]) % m) for a, b in s))
                    key = (shift(Cs), shift(Ds))
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append(key)
    return out


def parent_menu(ell, m):
    As, Bs = [], []
    for a in range(1, ell):
        for b in range(1, m):
            As.append(((a, 0), (0, b)))
            Bs.append(((0, b), (a, 0)))
            As.append(((0, 0), (a, 0), (0, b)))
            Bs.append(((0, 0), (0, b), (a, 0)))
    dedupA, dedupB = [], []
    for t in As:
        if t not in dedupA:
            dedupA.append(t)
    for t in Bs:
        if t not in dedupB:
            dedupB.append(t)
    return dedupA, dedupB


def main():
    lattices = [(2, 3), (3, 2), (2, 4), (4, 2), (3, 3), (2, 5), (5, 2), (2, 6), (6, 2), (3, 4), (4, 3)]
    stats = {"parents": 0, "parents_used": 0, "valid_perturbations": 0,
             "with_demotion": 0, "counterexamples": []}
    best = None
    for (ell, m) in lattices:
        dim = ell * m
        N = 2 * dim
        mc = MonoCache(ell, m)
        perts = perturbations(ell, m, 2 if dim > 9 else 3)
        As, Bs = parent_menu(ell, m)
        seen = set()
        for At in As:
            for Bt in Bs:
                A, B = mc.mono(At), mc.mono(Bt)
                key = (A.tobytes(), B.tobytes())
                if key in seen:
                    continue
                seen.add(key)
                stats["parents"] += 1
                HX = np.hstack([A, B])
                HZ = np.hstack([B.T, A.T])
                rX, rZ = rank_np(HX), rank_np(HZ)
                kP = N - rX - rZ
                if kP < 2:
                    continue
                Zcen = nullspace_np(HX)
                dXP, witP = dX_support(HZ, Zcen, N, 7)
                if dXP is None or dXP < 2:
                    continue
                stats["parents_used"] += 1
                for (Ct, Dt) in perts:
                    C, D = mc.mono(Ct), mc.mono(Dt)
                    M = (gf2matmul(A, C.T) ^ gf2matmul(B, D.T)).astype(np.uint8)
                    if (M ^ M.T).any():
                        continue
                    stats["valid_perturbations"] += 1
                    CD = np.hstack([C, D])
                    L = nullspace_np(HX.T)
                    Delta = (L @ CD) % 2 if L.shape[0] else np.zeros((0, N), np.uint8)
                    SZD = np.vstack([HZ, Delta]) if Delta.shape[0] else HZ
                    Wperp = nullspace_np(SZD)
                    G = (CD @ Wperp.T) % 2
                    Lam = nullspace_np(G.T)
                    SXQ = (Lam @ HX) % 2 if Lam.shape[0] else np.zeros((0, N), np.uint8)
                    kerM = nullspace_np(M.T)  # rows: {lam : lam M = 0}
                    if kerM.shape[0] == 0:
                        continue
                    mind, minx = span_outside(kerM, HX, SXQ, N)
                    if mind is None:
                        continue
                    stats["with_demotion"] += 1
                    if mind < dXP:
                        rec = {
                            "ell": ell, "m": m,
                            "A": [list(t) for t in At], "B": [list(t) for t in Bt],
                            "C": [list(t) for t in Ct], "D": [list(t) for t in Dt],
                            "k_P": int(kP), "d_X_P": int(dXP), "w_dem": int(mind),
                            "support": sorted(int(i) for i in np.nonzero(minx[0])[0]),
                        }
                        stats["counterexamples"].append(rec)
                        if best is None or mind - dXP < best[0]:
                            best = (mind - dXP, rec)
                            print("COUNTEREXAMPLE", json.dumps(rec), flush=True)
    print(json.dumps({k: v for k, v in stats.items() if k != "counterexamples"}))
    print("num counterexamples:", len(stats["counterexamples"]))
    out = {"stats": {k: (v if k != "counterexamples" else len(v)) for k, v in stats.items()},
           "best": best[1] if best else None,
           "all": stats["counterexamples"][:50]}
    (ROOT / "results/partial_runs/xsector/hunt.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
