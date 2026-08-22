"""Enriched small-lattice hunt for d_X(Q) < d_X(P) (stabilizer demotion).

d_X(Q) < d_X(P) iff some x in rowspace(H_X) has [C D]x^T = 0 (x in Xcen(Q))
but x not in S_X(Q) ("demoted") with wt(x) < d_X(P).  Hunt over small BB
parents and anchored valid perturbations of up to 5 monomials.
GF(2)/numpy only, one thread, no SAT.
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

from qec_research.codes.bicycle import poly_matrix  # noqa: E402
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


def dX_support(HZ, SXP_perp, n, cap):
    """Min wt x with H_Z x^T = 0 and x not in rowspace(H_X); cap-limited."""
    for w in range(1, cap + 1):
        combos = itertools.combinations(range(n), w)
        while True:
            chunk = list(itertools.islice(combos, 200_000))
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


def span_outside(kerM, HX, SXQ):
    """Min wt of {lam[A B] : lam in span(kerM)}\\span(SXQ); vectorized."""
    d = kerM.shape[0]
    if d == 0 or d > 18:
        return None, None
    ns = nullspace_np(SXQ) if SXQ.shape[0] else None
    best, bestx = None, None
    for start in range(1, 1 << d, 1 << 16):
        idx = np.arange(start, min(start + (1 << 16), 1 << d), dtype=np.int64)
        V = np.zeros((idx.size, kerM.shape[1]), dtype=np.uint8)
        for j in range(d):
            V ^= (((idx >> j) & 1)[:, None] & kerM[j]).astype(np.uint8)
        X = (V @ HX) % 2
        w = X.sum(axis=1)
        keep = w > 0
        if ns is not None and ns.shape[0]:
            keep &= ((X @ ns.T) % 2).sum(axis=1) != 0
        if keep.any():
            ww = w[keep]
            r = int(np.argmin(ww))
            if best is None or int(ww[r]) < best:
                best, bestx = int(ww[r]), X[keep][r]
    return best, (None if bestx is None else bestx[None, :])


def perturbations(ell, m, max_mons):
    grid = [(a, b) for a in range(ell) for b in range(m)]
    seen, out = set(), []
    for total in range(1, max_mons + 1):
        for split in range(total + 1):
            for Cs in itertools.combinations(grid, split):
                for Ds in itertools.combinations(grid, total - split):
                    allm = sorted(Cs + Ds)
                    t = allm[0]
                    def shift(s):
                        return tuple(sorted(((a - t[0]) % ell, (b - t[1]) % m) for a, b in s))
                    key = (shift(Cs), shift(Ds))
                    if key not in seen:
                        seen.add(key)
                        out.append(key)
    return out


def menu(ell, m):
    xs = list(range(1, ell))
    ys = list(range(1, m))
    S = []
    for a in xs:
        for b in ys:
            S.append(((a, 0), (0, b)))
            S.append(((0, 0), (a, 0), (0, b)))
            S.append(((a, 0), (0, b), (0, m - 1) if b != m - 1 else (0, 1)))
    for a in xs:
        S.append(((0, 0), (a, 0)))
        for a2 in xs:
            if a2 > a:
                S.append(((0, 0), (a, 0), (a2, 0)))
                S.append(((a, 0), (a2, 0)))
    for b in ys:
        S.append(((0, 0), (0, b)))
        for b2 in ys:
            if b2 > b:
                S.append(((0, 0), (0, b), (0, b2)))
                S.append(((0, b), (0, b2)))
    for a in xs:
        for b in ys:
            for b2 in ys:
                if b2 > b:
                    S.append(((a, 0), (0, b), (0, b2)))
                    pass
    for b in ys:
        for a in xs:
            for a2 in xs:
                if a2 > a:
                    S.append(((0, b), (a, 0), (a2, 0)))
    dedup, out = set(), []
    for t in S:
        k = tuple(sorted(t))
        if k not in dedup:
            dedup.add(k)
            out.append(k)
    return out


def main():
    spec = [
        ((2, 5), 3), ((5, 2), 3), ((2, 6), 3), ((3, 4), 3), ((4, 3), 3), ((6, 2), 3),
        ((2, 2), 5), ((2, 3), 5), ((3, 2), 5), ((2, 4), 4), ((4, 2), 4), ((3, 3), 4),
    ]
    stats = {"parents": 0, "parents_used": 0, "cand_perts": 0, "valid_perturbations": 0,
             "with_demotion": 0, "counterexamples": []}
    best = None
    rng = np.random.default_rng(7)
    for (ell, m), max_mons in spec:
        dim, N = ell * m, 2 * ell * m
        mc = MonoCache(ell, m)
        perts = perturbations(ell, m, max_mons)
        if len(perts) > 4000:
            ii = rng.choice(len(perts), 4000, replace=False)
            perts = [perts[i] for i in sorted(ii)]
        As = menu(ell, m)
        seen = set()
        for At in As:
            for Bt in As:
                A, B = mc.mono(At), mc.mono(Bt)
                key = (A.tobytes(), B.tobytes())
                if key in seen:
                    continue
                seen.add(key)
                stats["parents"] += 1
                HX = np.hstack([A, B])
                HZ = np.hstack([B.T, A.T])
                kP = N - rank_np(HX) - rank_np(HZ)
                if kP < 2:
                    continue
                Zcen = nullspace_np(HX)
                dXP, _ = dX_support(HZ, Zcen, N, 8)
                if dXP is None or dXP < 2:
                    continue
                stats["parents_used"] += 1
                for (Ct, Dt) in perts:
                    stats["cand_perts"] += 1
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
                    kerM = nullspace_np(M.T)
                    mind, minx = span_outside(kerM, HX, SXQ)
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
        print(f"lattice ({ell},{m}) done: {json.dumps({k: v for k, v in stats.items() if k != 'counterexamples'})}", flush=True)
    print(json.dumps({k: v for k, v in stats.items() if k != "counterexamples"}))
    print("num counterexamples:", len(stats["counterexamples"]))
    out = {"stats": {k: (v if k != "counterexamples" else len(v)) for k, v in stats.items()},
           "best": best[1] if best else None, "all": stats["counterexamples"][:100]}
    (ROOT / "results/partial_runs/xsector/hunt_enriched.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
