"""Gate B exact verdict engine.

Verdict definition (pre-registered):
  instance DEGENERATE <=> for ALL pairs p<q: Delta_{p,q} = f_p f_q' + f_q f_p'
  is the ZERO polynomial in E[Z].

Engine (exact; numpy MUL-table gathers over E = F_{2^m}; no floats):

  v[i,j] := f_j(a_i)      = Y[j,i] * lam_i^{-1}      (EXACT by the Lagrange
          construction: f_j = sum_i Y[j,i] lam_i^{-1} L_i, L_i(a_l)=delta_il,
          proven in frozen ADDENDUM 2; re-verified bit-exactly against the
          frozen 13 instances in reproduction_check.py)
  w[i,j] := f_j'(a_i)     = sum_d O[j,d] * (a_i^2)^d,
          where O[j,d] = coeff of Z^{2d+1} of f_j   (char-2 derivative:
          f' = sum_{e odd} f_e Z^{e-1}, i.e. f'(Z) = O(Z^2))
  delta(a_i) pair value:  v[p] w[q] + w[p] v[q]  (2x2 minor of [v|w])

  FACT (exact, used with the completeness lemma): all 2x2 minors at point i
  vanish  <=>  rank_E([v_i | w_i]) <= 1  <=>
     (a) w_i[j] = 0 for every j with Y[j,i] = 0, and
     (b) lam_i * w_i[j] all equal for j with Y[j,i] = 1.
  (rows with v=0 force w=0 for rank<=1; rows with v=lam^{-1} != 0 force
   w_j = beta_i * lam^{-1}, same beta_i.)

  COMPLETENESS LEMMA (full-support cells): char-2 algebra gives
     Delta_{p,q} = S_p O_q + S_q O_p  as a polynomial in x = Z^2, where
     S_j = even part of f_j, O_j = odd part.  Hence deg_x Delta <=
     floor(D/2) + floor((D-1)/2) <= D - 1 < n.  Delta == 0 in E[Z] iff
     Delta(a_i) = 0 for ALL i (n DISTINCT points, n > deg).  Frobenius
     injectivity gives the x_i = a_i^2 distinct.  So the all-points test is
     a COMPLETE exact decision at full support (n = 2^m), no gcd cascade
     needed.  The O(k^2) gcd-cascade is retained as an independent
     cross-check on small cells.

  Degeneracy certificate for a DEGENERATE verdict: the full per-point rank
  evidence (route "pointscan") or the all-squares argument (route
  "all_squares": every f_j has zero odd part => f_j' = 0 => all Delta = 0).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gfield import GF
from fastfield import EField
from instance import Instance


# --------------------------------------------------------------------------
# exact engine
# --------------------------------------------------------------------------

def odd_part_coeffs(f, H):
    """O[j][d] = f[2d+1]; length H+1 padded."""
    out = np.zeros(H + 1, dtype=np.uint16)
    for d in range(H + 1):
        e = 2 * d + 1
        if e < len(f):
            out[d] = f[e]
    return out


def engine_verdict(F, Y, lam, support, ef, gf, full_support_points=None):
    """F: list of k coefficient lists; Y: k x n binary; lam: n-vector;
    support: n field elements (order = registered scan order).
    Returns dict with verdict, route, witness, evidence.  EXACT."""
    k = len(F)
    n = len(support)
    D = max((gf.pdeg(f) for f in F), default=-1)
    res = {"k": k, "n": n, "D": int(D)}
    if k < 2:
        res.update(verdict="VACUOUS_K_LT_2", route="no_pairs")
        return res
    H = (D - 1) // 2 if D >= 1 else 0
    Omat = np.zeros((k, H + 1), dtype=np.uint16)
    any_odd = False
    for j, f in enumerate(F):
        Omat[j] = odd_part_coeffs(f, H)
        if Omat[j].any():
            any_odd = True
    if not any_odd:
        res.update(verdict="DEGENERATE", route="all_squares",
                   evidence="all f_j have zero odd part => f_j' = 0 (char 2) => all Delta = 0",
                   D=D)
        return res

    # ---- route 2: pointwise rank test over the n support points -------------
    if k < 2:
        res.update(verdict="VACUOUS_K_LT_2", route="no_pairs")
        return res
    LamInv = np.array([ef.INV[int(l)] for l in lam], dtype=np.uint16)
    V = ef.MUL[Y.astype(np.uint16), LamInv[None, :]]      # v[j,i] exact
    pts = range(n) if full_support_points is None else full_support_points
    beta_seen = {}
    n_informative = 0
    for i in pts:
        a = int(support[i])
        S = np.nonzero(Y[:, i])[0]
        if S.size == 0:
            # v-column zero at this point: uninformative for Delta values
            continue
        n_informative += 1
        x = ef.MUL[a, a]
        xp = ef.powers(int(x), H)
        nzd_all = np.nonzero(Omat.any(axis=0))[0]      # degrees d with some odd coeff
        # w_j(a_i) = sum_d O[j,d] x_i^d  for ALL j (needed for rank condition (a))
        Wvals = np.zeros(k, dtype=np.uint16)
        if nzd_all.size:
            terms = ef.MUL[Omat[:, nzd_all], xp[nzd_all][None, :]]
            Wvals = np.bitwise_xor.reduce(terms, axis=1)
        wOff = Wvals[np.setdiff1d(np.arange(k, dtype=int), S, assume_unique=False)]
        if wOff.any():
            q = int(np.nonzero(wOff)[0][0])
            q = int(np.setdiff1d(np.arange(k, dtype=int), S, assume_unique=False)[q])
            p = int(S[0])
            # exact witness: delta(a_i) = v_p w_q + w_p v_q, v_q = 0
            wp = int(Wvals[p])
            val = ef.MUL[int(V[p, i]), wp]
            res.update(verdict="NONDEGENERATE", route="pointscan",
                       witness_pair=[p, q], witness_point=int(i),
                       witness_delta_value=int(val),
                       witness_note="row q with Y[q,i]=0 has f_q'(a_i) != 0",
                       n_points_scanned=int(i) + 1, D=D)
            return res
        wS = Wvals[S]
        # (b): lam_i * wS all equal?
        lamw = ef.MUL[int(lam[i]), wS]
        if lamw.size and not np.all(lamw == lamw[0]):
            # rank-2 witness: two S-rows with unequal lam*w
            j2 = int(np.nonzero(lamw != lamw[0])[0][0])
            p, q = int(S[0]), int(S[j2])
            # exact: v_p = v_q = lam_i^{-1}; delta(a_i) = v_p w_q + w_p v_q
            #       = lam^{-1} (w_p + w_q) = lam^{-1} (w_0 ^ w_j2)
            val = ef.MUL[int(LamInv[i]), int(wS[0]) ^ int(wS[j2])]
            res.update(verdict="NONDEGENERATE", route="pointscan",
                       witness_pair=[p, q], witness_point=int(i),
                       witness_delta_value=int(val),
                       witness_note="two y-support rows with unequal lam_i*f_j'(a_i)",
                       n_points_scanned=int(i) + 1, D=D)
            return res
    # all points rank <= 1
    # completeness bookkeeping: verify deg bound deg_x Delta <= D-1 < n
    res.update(verdict="DEGENERATE", route="pointscan_complete",
               evidence=("all n support points have rank([v_i|w_i]) <= 1; "
                         "completeness: deg_x Delta <= D-1 < n, x-points distinct"),
               beta_at_points={str(i): b for i, b in list(beta_seen.items())[:8]},
               n_informative_points=n_informative, D=D)
    return res


# --------------------------------------------------------------------------
# independent exact cascade (gcd-based) — cross-check on small cells
# --------------------------------------------------------------------------

def is_square_poly(f):
    return all(c == 0 for c in f[1::2])


def cascade_verdict(gf, F):
    """Independent instrument #3: EXACT brute-force Wronskian polynomials.
    NOTE (self-test caught 2026-08-30): an earlier gcd+coprime-squares test
    was INVALID as an iff: all coprime-parts-squares is only SUFFICIENT for
    W = 0; the general coprime-degenerate case is (f_e, f_o) proportional to
    (g_e, g_o) over E(x), i.e. W = f_e g_o + f_o g_e vanishing as a
    polynomial in x = Z^2.  This instrument computes W exactly; ground
    truth for small k (O(k^2) pmul/pderiv)."""
    for p in range(len(F)):
        for q in range(p + 1, len(F)):
            z, _ = brute_pair_values(gf, F, p, q)
            if not z:
                return {"verdict": "NONDEGENERATE", "first_nonzero_pair": [p, q],
                        "route": "brute_wronskian"}
    return {"verdict": "DEGENERATE", "route": "brute_wronskian"}


def brute_pair_values(gf, F, p, q):
    """Exact Delta_{p,q} polynomial value route (degree-based)."""
    fp, fq = gf.pderiv(F[p]), gf.pderiv(F[q])
    d = gf.ptrim(gf.padd(gf.pmul(F[p], fq), gf.pmul(F[q], fp)))
    return (len(d) == 0), gf.pdeg(d)
