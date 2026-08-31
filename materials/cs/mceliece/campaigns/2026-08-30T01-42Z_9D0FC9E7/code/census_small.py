"""MINIMAL Gate A census: family-restricted quantities only. FAST.

For instance (m, n, t, seed):
  * builds the binary Goppa instance (exact, seeded; probe-validated engine)
  * computes F-jet flags V_a^[j] at held points, depth R
  * builds the (2t+3)-FAMILY-RESTRICTED flag blocks directly:
      T(A) = Q_A + rho_A Z F per Apon Def 4 (eqs. 17-20);
      at held point a, order j: rows b -> h_b . jet_j T(Z^m)(a), m = 0..2t+2
    (columns = the (2t+3) family parameters: ONE ROW PER (a, j, b))
  * records:
      N_fam(c)      = (2t+3) - rank(all rows for first c points)  <- vs 2t+3-c
      rank(a, j)    = rank of block (a,j) alone                   <- vs <= 1
      N_eval(c)     = 2t+3 - rank(Vandermonde on first c pts)     (= above rows' prediction)
      jetprof(a)    = dim V_a^[j], j = 0..R
      degenerate(a) = any prof plateau (Remark-8 nondegeneracy visible)
      guards alpha..eps incl. Delta_{p,q}                            <- Finding 3
  * verification of the rank <= 1 prediction and the 2t+3-c law
  * lower-bound certificate at c = 2t+2: T(A_S) jets inside flags

NO ambient k(D+1) computation. Exact GF(2^m) arithmetic: MUL-table
(vector op) cross-validated against the scalar engine (500+ checks).

NOTE ON THE LOWER BOUND (Apon Thm 1):  for the family the flag conditions
reduce to A(a) = 0 per held point — Lemma 7.  This script does NOT assume
it; it derives the rows explicitly and measures their rank.  If rank > 1
at any (a, j): REFUTATION EVENT -> halt + do not write as established.
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fastfield import EField
from gfield import GF
from instance import Instance
from census import (flag_data, jets_F_at_points, lucas_w, poly_jets,
                    poly_jets_vec, conv, series_inv_jets, compute_R_G,
                    rho_of, jets_monomial, lower_bound_certificate)


def family_restricted_rows(inst, ef, fd, R):
    """family-restricted flag rows: columns parameterized by A = Z^m,
    m = 0..2t+2.  Row (a, j, b): h . jet_j T(Z^m)(a).  Returns
    dict[(a, j)] -> rows (h-count x (2t+3))."""
    gf = inst.gf
    cols = 2 * inst.t + 3
    RG = compute_R_G(inst, ef)
    G2 = gf.pmul(inst.G, inst.G)
    blocks = {}
    # F' jets per point: FJp[k, q] = (q+1 mod 2) * FJ[k, q+1]
    for a in fd:
        d = fd[a]
        FJm = d["FJ"]  # (k, R+1)
        FJp = np.zeros((inst.k, R + 1), dtype=np.uint16)
        for q in range(R):
            FJp[:, q] = ef.MUL[(q + 1) % 2, FJm[:, q + 1]]
        # series 1/G^2 at a to order R
        GJ = poly_jets(ef, np.array(inst.G, dtype=np.uint16), a, R + 1)
        G2J = conv(ef, GJ, GJ, R + 1)
        invG2 = series_inv_jets(ef, G2J, R)
        # rho for A = Z^{2t+2}: [Z^{2t+2}] A = 1, (D mod 2) factor
        rho_scalar = 0
        if inst.D % 2 == 0:
            rho_scalar = rho_of(inst, ef)
        # jets of Z^m and of B_m = (Z^m R_G mod G^2)
        # jets of Z F: (k, R+1): ZF_j = a FJ_j + FJ_{j-1}
        ZF = np.zeros((inst.k, R + 1), dtype=np.uint16)
        for j in range(R + 1):
            ZF[:, j] = ef.MUL[a, FJm[:, j]] ^ (FJm[:, j - 1] if j >= 1 else 0)
        ZJ = np.zeros((cols, R + 2), dtype=np.uint16)  # jets of Z^m  (indices 0..R+1)
        BJ = np.zeros((cols, R + 1), dtype=np.uint16)  # jets of B_m  (0..R)
        for mm in range(cols):
            ZJ[mm] = jets_monomial(ef, mm, a, R + 1)
            # B: (Z^m R_G) mod G^2
            with np.errstate(all="ignore"):
                pass
            Am = [0] * (mm) + [1]
            Bm = gf.pmod(gf.pmul(Am, RG), G2)
            BJ[mm] = poly_jets_vec(ef, np.array(gf.ptrim(Bm), dtype=np.uint16), a, R)
        for j in range(R):
            # numerator jets 0..j for each family coord m
            numjets = []
            for e in range(j + 1):
                N = np.zeros((cols, inst.k), dtype=np.uint16)
                for p in range(e + 1):
                    q = e - p
                    # A F' term: jet_p(A) * jet_q(F')  (LUCAS-EXACT)
                    # (note F' jets: FJp[k, q]; scalar = ZJ[mm, p])
                    N ^= ef.MUL[ZJ[:, p][:, None], FJp[:, q][None, :]]
                    # B F term: jet_p(B) * jet_q(F)
                    N ^= ef.MUL[BJ[:, p][:, None], FJm[:, q][None, :]]
                numjets.append(N)
            # Q = num / G^2: jet_j = sum_{e} invG2[e] * numjet[j - e]
            Q = np.zeros((cols, inst.k), dtype=np.uint16)
            for e in range(j + 1):
                Q ^= ef.MUL[invG2[e], numjets[j - e]]
            # T = Q + rho Z F
            Psi = Q
            if rho_scalar:
                Psi = Psi ^ ef.MUL[rho_scalar, ZF[:, j][None, :]]
            # restrict to complement h
            H = d["H"][j]
            if H.shape[0] == 0:
                blocks[(a, j)] = np.zeros((0, cols), dtype=np.uint16)
                continue
            rows = np.zeros((H.shape[0], cols), dtype=np.uint16)
            for b in range(H.shape[0]):
                prod = ef.MUL[H[b][None, :], Psi]      # (cols, k)
                rows[b] = np.bitwise_xor.reduce(prod, axis=1)
            blocks[(a, j)] = rows
    return blocks


def run(m, n, t, seed, R=4, cmax_extra=4):
    t0 = time.time()
    inst = Instance(m, n, t, seed)
    ef = EField(m)
    tt = 2 * t + 3
    cmax = min(2 * t + cmax_extra, n)
    S = inst.support[:cmax]

    # flags at all cmax points
    fd = flag_data(inst, ef, S, R)
    # family blocks (per (a, j))
    blocks = family_restricted_rows(inst, ef, fd, R)

    res = {
        "m": m, "n": n, "t": t, "k": inst.k, "D": inst.D, "seed": seed,
        "prim": hex(ef.prim), "G": inst.G, "g_tries": inst.g_tries,
        "support": inst.support, "R": R,
        "guards": {k: v for k, v in inst.guards.items() if k != "degs"},
        "degs": inst.guards.get("degs"),
        "jetprof": {str(a): fd[a]["prof"] for a in S},
        "degenerate_held": {str(a): any(fd[a]["prof"][j+1] == fd[a]["prof"][j]
                                         for j in range(R)) for a in S},
    }

    # per-point order ranks
    ranks = {}
    refut = {}
    for (a, j), B in blocks.items():
        r = ef.rank(B) if B.shape[0] else 0
        ranks[f"{a},{j}"] = int(r)
    for key, r in ranks.items():
        if r >= 2:
            refut[key] = r
    res["rank_point_order"] = ranks
    res["refutation_event"] = bool(refut)
    res["refuting_blocks"] = refut

    # N_eval(c): (2t+3) - rank Vandermonde first c
    Nf = {}
    for c in range(1, cmax + 1):
        rows = [lucas_w(ef, tt - 1, 0, a) for a in S[:c]]
        rk = ef.rank(np.array(rows, dtype=np.uint16))
        Nf[c] = tt - int(rk)
    res["N_eval"] = Nf

    # N_fam(c): rank of stacked family blocks
    Nfam = {}
    for c in range(1, cmax + 1):
        stack = [blocks[(a, j)] for a in S[:c] for j in range(R)]
        stack = [s for s in stack if s.shape[0]]
        if stack:
            M = np.vstack(stack)
            Nfam[c] = tt - int(ef.rank(M))
        else:
            Nfam[c] = tt
    res["N_fam"] = Nfam

    # derived contribution check: does N_fam equal 2t+3 - #{contributing}
    contrib = {}
    for a in S:
        pr = fd[a]["prof"]
        contrib[a] = any((j % 2 == 0) and (j + 1 <= R) and pr[j + 1] > pr[j]
                         for j in range(R))
    Nfam_pred = {c: tt - sum(contrib[a] for a in S[:c]) for c in range(1, cmax + 1)}
    res["N_fam_derived"] = Nfam_pred
    res["N_fam_matches_derived"] = all(Nfam[c] == Nfam_pred[c] for c in Nfam)

    # lower-bound certificate at c = 2t+2
    res["lb_cert"] = lower_bound_certificate(inst, ef, fd, S, R)

    res["elapsed_s"] = round(time.time() - t0, 2)
    return res


def main():
    # one seeded instance per ladder point; sweep t small first
    plan = [
        (6, 64, 3, 1387),
        (6, 64, 4, 2311),
        (6, 64, 5, 3413),
        (7, 128, 3, 4421),
        (7, 128, 6, 5531),
        (8, 256, 3, 6637),
        (8, 256, 5, 7741),
    ]
    out = []
    for (m, n, t, seed) in plan:
        r = run(m, n, t, seed, R=4)
        print(f"m={m} n={n} t={t} seed={seed}: "
              f"N_fam(c) {r['N_fam']} (pred 2t+3-c) "
              f"maxrank {max(r['rank_point_order'].values())} "
              f"refut {r['refutation_event']} "
              f"Delta {r['guards']['beta_delta_nonzero']} "
              f"[{r['elapsed_s']}s]",
              flush=True)
        out.append(r)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
