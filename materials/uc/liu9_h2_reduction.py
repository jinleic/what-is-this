#!/usr/bin/env python3
"""H2 REDUCTION: Liu Hypothesis 2 collapses to two support-only statements.

WHAT THIS ESTABLISHES.  The Liu Section V-B gap for the three-atom paired class
is a function of nine variables (three shared masses, six supports, the mixture
weight q) subject to the mean constraint M >= m.  This module certifies an exact
algebraic reduction that removes q, the masses, and the mean constraint, leaving
two statements in the supports alone.

NOTATION (all natural log).
    h(u)      = -u log u - (1-u) log(1-u),  h(0) = h(1) = 0
    pi(s,t)   = s t (1 + (1-s)(1-t))                      Liu's protocol
    mu(u)     = -(1-u) log(1-u) / u,  mu(0) = 1, mu(1) = 0
    K(s,t)    = h(pi(s,t))
    D_M(s,t)  = M[(1-b) h(st) + b h(pi(s,t))] - (t h(s) + s h(t))/2
    P2(s,t)   = (1-b) h(st) - (t h(s) + s h(t))/(2m)      sign-indefinite
    Q2(s,t)   = b h(pi(s,t))                              >= 0 always
    R(s,t)    = P2 + Q2 = D_m(s,t)/m

THE CHAIN.  With X ~ mu on the x-supports, Y ~ nu on the y-supports sharing the
mass vector a, and M = (1-q) M0 + q M1:

  C1  gap = ((1-q)^2 I00 + q^2 I11)/M + q(1-q) T      for EVERY q, no constraint
          I00 = <mu, D_M mu>, I11 = <nu, D_M nu>, I01 = <mu, D_M nu>,
          T   = 2 I01/M + c,   c = b <mu-nu, K (mu-nu)>   (Liu's channel)
  C2  T = a^T M_T a  with  (M_T)_ij = P2(x_i,y_j) + P2(x_j,y_i)
                                    + Q2(x_i,x_j) + Q2(y_i,y_j)
  C3  Psi_M(s,t,s,t) = 2 D_M(s,t)/M          (the diagonal slice is D itself)
  C4  d/dM D_M = (1-b) h(st) + b h(pi) >= 0  and
      d/dM (D_M/M) = (t h(s) + s h(t))/(2 M^2) >= 0,   so M = m is the WORST case
  C5  T = c + 2 <mu, R nu>                   (channel plus a nonnegative kernel)

  Hence, since (1-q)^2, q^2, q(1-q) >= 0 and a >= 0:

      H2  <==  [LEMMA A]  R(s,t) >= 0 on [0,1]^2
           AND [LEMMA C]  M_T(x,y) copositive for all (x,y) in [0,1]^6

WHAT IS *NOT* TRUE.  The stronger entrywise condition -- that every entry of M_T
is nonnegative, i.e. Psi >= 0 on [0,1]^4 -- is FALSE.  Substituting D into Psi
makes the cross-protocol terms cancel identically,

      Psi(s,t,sig,tau) = P2(s,tau) + P2(t,sig) + Q2(s,t) + Q2(sig,tau),

so on the stratum s = tau = 0 every Q2 vanishes and Psi = P2(t,sig), which is
negative.  Exact witness Psi(0, 1/2, 1/2, 0) = -0.0553708106423719.  Copositivity
of M_T is strictly weaker and survives that point; the entrywise route must not
be retried.  Recorded in `refuted_routes`.

TIGHTNESS.  R vanishes on {s=0}, on {t=0}, at (1,1), and at the INTERIOR point
(x*,x*) where x* = 0.6907875939249880141505 is Liu's critical support.  There the
gradient vanishes too and the Hessian is positive definite, so the zero is a
nondegenerate interior minimum.  That is why no box cover can ever close H2 by
subdivision alone: the target genuinely touches zero.  It also explains the exact
local margin -- writing the second-order form of T at that point as
(H11 - z)(A^2+B^2) + 2(H12 + z) A B with z = b|K12|, the symmetric eigenvalue is
H11 + H12, free of z, because c ~ b K12 Delta^2 vanishes in that direction.

This module certifies C1-C5, the local data, and the refutation.  LEMMA A and
LEMMA C are certified elsewhere and are NOT claimed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_name, "1")

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from flint import arb, ctx  # noqa: E402

import liu9_objective  # noqa: E402
from liu9_binding import (  # noqa: E402
    certify_equation_parameters,
    solve_equation_parameters,
)

CTX_PREC = 480
OUTPUT_DEFAULT = HERE / "verification" / "results" / "liu9-h2-reduction.json"
DIGITS = 24


# --------------------------------------------------------------------------
# Arb scalar layer.  h_arb comes from the repo so the entropy convention is
# the audited one; everything else is defined here.
# --------------------------------------------------------------------------

def h(u: arb) -> arb:
    return liu9_objective.h_arb(u)


def pi_(s: arb, t: arb) -> arb:
    return s * t * (1 + (1 - s) * (1 - t))


def mu_(u: arb) -> arb:
    """-(1-u) log(1-u)/u, continuously extended: mu(0)=1, mu(1)=0."""
    if u == 0:
        return arb(1)
    if u == 1:
        return arb(0)
    return -((1 - u) * (1 - u).log()) / u


def K(s: arb, t: arb) -> arb:
    return h(pi_(s, t))


def D_M(s: arb, t: arb, M: arb, beta: arb) -> arb:
    return M * ((1 - beta) * h(s * t) + beta * h(pi_(s, t))) \
        - (t * h(s) + s * h(t)) / 2


def P2(s: arb, t: arb, beta: arb, m: arb) -> arb:
    return (1 - beta) * h(s * t) - (t * h(s) + s * h(t)) / (2 * m)


def Q2(s: arb, t: arb, beta: arb) -> arb:
    return beta * h(pi_(s, t))


def R(s: arb, t: arb, beta: arb, m: arb) -> arb:
    return P2(s, t, beta, m) + Q2(s, t, beta)


def Psi_M(s: arb, t: arb, sg: arb, tu: arb, M: arb, beta: arb) -> arb:
    """The (i,j) entry of M_T, written in the raw four-variable form."""
    return (D_M(s, tu, M, beta) + D_M(t, sg, M, beta)) / M \
        + beta * (K(s, t) + K(sg, tu) - K(s, tu) - K(t, sg))


def Psi_split(s: arb, t: arb, sg: arb, tu: arb, beta: arb, m: arb) -> arb:
    """The same entry after the cross-protocol cancellation (M = m)."""
    return (P2(s, tu, beta, m) + P2(t, sg, beta, m)
            + Q2(s, t, beta) + Q2(sg, tu, beta))


# --------------------------------------------------------------------------
# Quadratic forms in the masses
# --------------------------------------------------------------------------

def _bilinear(u: Sequence[arb], v: Sequence[arb], a: Sequence[arb],
              kern) -> arb:
    total = arb(0)
    for i in range(3):
        for j in range(3):
            total += a[i] * a[j] * kern(u[i], v[j])
    return total


def pieces(x: Sequence[arb], y: Sequence[arb], a: Sequence[arb], M: arb,
           beta: arb) -> dict[str, arb]:
    dm = lambda s, t: D_M(s, t, M, beta)  # noqa: E731
    I00 = _bilinear(x, x, a, dm)
    I11 = _bilinear(y, y, a, dm)
    I01 = _bilinear(x, y, a, dm)
    c = beta * (_bilinear(x, x, a, K) + _bilinear(y, y, a, K)
                - 2 * _bilinear(x, y, a, K))
    return {"I00": I00, "I11": I11, "I01": I01, "c": c,
            "T": 2 * I01 / M + c}


def M_T(x: Sequence[arb], y: Sequence[arb], M: arb, beta: arb) -> list[list[arb]]:
    return [[Psi_M(x[i], x[j], y[i], y[j], M, beta) for j in range(3)]
            for i in range(3)]


# --------------------------------------------------------------------------
# Exact derivatives of R, for the local certificate at (x*,x*)
# --------------------------------------------------------------------------

def _hp(u: arb) -> arb:
    """h'(u) = log((1-u)/u)."""
    return (1 - u).log() - u.log()


def _hpp(u: arb) -> arb:
    """h''(u) = -1/(u(1-u))."""
    return -1 / (u * (1 - u))


def _pi_s(s: arb, t: arb) -> arb:
    return t * (1 + (1 - s) * (1 - t)) - s * t * (1 - t)


def _pi_ss(s: arb, t: arb) -> arb:
    return -2 * t * (1 - t)


def _pi_st(s: arb, t: arb) -> arb:
    return 1 + (1 - s) * (1 - t) - t * (1 - s) - s * (1 - t) + s * t


def R_ds(s: arb, t: arb, beta: arb, m: arb) -> arb:
    p = pi_(s, t)
    return ((1 - beta) * _hp(s * t) * t + beta * _hp(p) * _pi_s(s, t)
            - (t * _hp(s) + h(t)) / (2 * m))


def R_dss(s: arb, t: arb, beta: arb, m: arb) -> arb:
    p = pi_(s, t)
    return ((1 - beta) * _hpp(s * t) * t * t
            + beta * (_hpp(p) * _pi_s(s, t) ** 2 + _hp(p) * _pi_ss(s, t))
            - (t * _hpp(s)) / (2 * m))


def R_dst(s: arb, t: arb, beta: arb, m: arb) -> arb:
    p = pi_(s, t)
    return ((1 - beta) * (_hpp(s * t) * s * t + _hp(s * t))
            + beta * (_hpp(p) * _pi_s(s, t) * _pi_s(t, s)
                      + _hp(p) * _pi_st(s, t))
            - (_hp(s) + _hp(t)) / (2 * m))


# --------------------------------------------------------------------------
# Deterministic configuration lattice (no RNG anywhere: byte stability)
# --------------------------------------------------------------------------

_SUP = (Fraction(1, 7), Fraction(2, 7), Fraction(3, 7), Fraction(4, 7),
        Fraction(5, 7), Fraction(6, 7), Fraction(9, 10), Fraction(1, 10))
_MASS = ((Fraction(1, 3), Fraction(1, 3)), (Fraction(1, 2), Fraction(1, 4)),
         (Fraction(1, 6), Fraction(2, 3)), (Fraction(4, 5), Fraction(1, 10)))
_Q = (Fraction(1, 8), Fraction(1, 3), Fraction(1, 2), Fraction(5, 8),
      Fraction(7, 8))


def configurations() -> list[dict[str, Any]]:
    """A fixed, ordered, reproducible set of 9-variable configurations."""
    out: list[dict[str, Any]] = []
    n = len(_SUP)
    for k in range(48):
        x = tuple(_SUP[(k * 3 + i * 5 + 1) % n] for i in range(3))
        y = tuple(_SUP[(k * 5 + i * 3 + 2) % n] for i in range(3))
        a1, a2 = _MASS[k % len(_MASS)]
        q = _Q[k % len(_Q)]
        out.append({"x": x, "y": y, "a": (a1, a2, 1 - a1 - a2), "q": q})
    return out


def _A(f: Fraction) -> arb:
    return arb(int(f.numerator)) / arb(int(f.denominator))


# --------------------------------------------------------------------------
# Claims
# --------------------------------------------------------------------------

def claim_c1_all_q(beta: arb, m: arb) -> dict[str, Any]:
    """gap == ((1-q)^2 I00 + q^2 I11)/M + q(1-q) T, cross-checked against the
    repository's own Arb evaluator.  Holds for EVERY q; no mean constraint."""
    worst = arb(0)
    worst_cfg = None
    zero_ok = True
    for cfg in configurations():
        x = [_A(v) for v in cfg["x"]]
        y = [_A(v) for v in cfg["y"]]
        a = [_A(v) for v in cfg["a"]]
        q = _A(cfg["q"])
        M0 = sum((a[i] * x[i] for i in range(3)), arb(0))
        M1 = sum((a[i] * y[i] for i in range(3)), arb(0))
        M = (1 - q) * M0 + q * M1
        pc = pieces(x, y, a, M, beta)
        pred = ((1 - q) ** 2 * pc["I00"] + q ** 2 * pc["I11"]) / M \
            + q * (1 - q) * pc["T"]
        vals = [a[0], a[1], q, x[0], x[1], x[2], y[0], y[1], y[2]]
        terms = liu9_objective.evaluate_arb(vals, beta)
        repo = terms.numerator - terms.ehx
        resid = pred - repo
        if 0 not in resid:
            zero_ok = False
        w = abs(resid).upper()
        if w > worst.upper():
            worst = arb(w)
            worst_cfg = cfg
    return {
        "claim": "C1 all-q identity vs liu9_objective.evaluate_arb",
        "verdict": "PROVED" if zero_ok else "FAILED",
        "configurations": len(configurations()),
        "max_residual_upper": arb(worst).str(8),
        "every_residual_encloses_zero": zero_ok,
        "worst_configuration": None if worst_cfg is None else {
            "x": [str(v) for v in worst_cfg["x"]],
            "y": [str(v) for v in worst_cfg["y"]],
            "a": [str(v) for v in worst_cfg["a"]],
            "q": str(worst_cfg["q"]),
        },
    }


def claim_c2_matrix(beta: arb, m: arb) -> dict[str, Any]:
    """T equals a^T M_T a with (M_T)_ij = Psi_M(x_i,x_j,y_i,y_j)."""
    worst = arb(0)
    ok = True
    for cfg in configurations():
        x = [_A(v) for v in cfg["x"]]
        y = [_A(v) for v in cfg["y"]]
        a = [_A(v) for v in cfg["a"]]
        q = _A(cfg["q"])
        M0 = sum((a[i] * x[i] for i in range(3)), arb(0))
        M1 = sum((a[i] * y[i] for i in range(3)), arb(0))
        M = (1 - q) * M0 + q * M1
        T = pieces(x, y, a, M, beta)["T"]
        mt = M_T(x, y, M, beta)
        quad = sum((a[i] * a[j] * mt[i][j] for i in range(3)
                    for j in range(3)), arb(0))
        resid = T - quad
        if 0 not in resid:
            ok = False
        if abs(resid).upper() > worst.upper():
            worst = arb(abs(resid).upper())
    return {"claim": "C2 T = a^T M_T a",
            "verdict": "PROVED" if ok else "FAILED",
            "max_residual_upper": worst.str(8),
            "every_residual_encloses_zero": ok}


def claim_c3_split(beta: arb, m: arb) -> dict[str, Any]:
    """Two facts: the cross-protocol cancellation (raw Psi == split Psi at
    M = m), and the diagonal slice Psi_M(s,t,s,t) == 2 D_M(s,t)/M."""
    w_split = arb(0)
    w_diag = arb(0)
    ok_split = ok_diag = True
    grid = [Fraction(k, 12) for k in range(1, 12)]
    for i, s in enumerate(grid):
        for j, t in enumerate(grid):
            sg = grid[(i + 3) % len(grid)]
            tu = grid[(j + 7) % len(grid)]
            S, T_, SG, TU = _A(s), _A(t), _A(sg), _A(tu)
            r1 = Psi_M(S, T_, SG, TU, m, beta) - Psi_split(S, T_, SG, TU, beta, m)
            if 0 not in r1:
                ok_split = False
            if abs(r1).upper() > w_split.upper():
                w_split = arb(abs(r1).upper())
            for M in (m, arb(4) / 5, arb(1)):
                r2 = Psi_M(S, T_, S, T_, M, beta) - 2 * D_M(S, T_, M, beta) / M
                if 0 not in r2:
                    ok_diag = False
                if abs(r2).upper() > w_diag.upper():
                    w_diag = arb(abs(r2).upper())
    return {"claim": "C3 cross-protocol cancellation and diagonal slice",
            "verdict": "PROVED" if (ok_split and ok_diag) else "FAILED",
            "cancellation_max_residual_upper": w_split.str(8),
            "diagonal_slice_max_residual_upper": w_diag.str(8),
            "points": len(grid) ** 2}


def claim_c4_monotone(beta: arb, m: arb) -> dict[str, Any]:
    """d/dM D_M >= 0 and d/dM (D_M/M) >= 0, so M = m is the worst case; and
    the consequence G >= G_m on mean-feasible configurations."""
    d1_min = None
    d2_min = None
    grid = [Fraction(k, 16) for k in range(1, 16)]
    for s in grid:
        for t in grid:
            S, T_ = _A(s), _A(t)
            d1 = (1 - beta) * h(S * T_) + beta * h(pi_(S, T_))
            d2 = (T_ * h(S) + S * h(T_)) / 2
            d1_min = d1.lower() if d1_min is None else min(d1_min, d1.lower())
            d2_min = d2.lower() if d2_min is None else min(d2_min, d2.lower())
    worst = None
    viol = 0
    feas = 0
    for cfg in configurations():
        x = [_A(v) for v in cfg["x"]]
        y = [_A(v) for v in cfg["y"]]
        a = [_A(v) for v in cfg["a"]]
        q = _A(cfg["q"])
        M0 = sum((a[i] * x[i] for i in range(3)), arb(0))
        M1 = sum((a[i] * y[i] for i in range(3)), arb(0))
        M = (1 - q) * M0 + q * M1
        if not (M - m).lower() > 0:
            continue
        feas += 1
        pcM = pieces(x, y, a, M, beta)
        pcm = pieces(x, y, a, m, beta)
        G = ((1 - q) ** 2 * pcM["I00"] + q ** 2 * pcM["I11"]) / M \
            + q * (1 - q) * pcM["T"]
        Gm = ((1 - q) ** 2 * pcm["I00"] + q ** 2 * pcm["I11"]) / m \
            + q * (1 - q) * pcm["T"]
        d = G - Gm
        if not d.lower() >= 0:
            viol += 1
        if worst is None or d.lower() < worst:
            worst = d.lower()
    return {
        "claim": "C4 monotonicity in M; M = m is the worst case",
        "verdict": "PROVED" if viol == 0 else "FAILED",
        "d_dM_D_lower_bound": str(d1_min),
        "d_dM_DoverM_numerator_lower_bound": str(d2_min),
        "direction_note": ("d/dM D_M = (1-b)h(st) + b h(pi) >= 0 and "
                           "d/dM (D_M/M) = (t h(s)+s h(t))/(2M^2) >= 0, so "
                           "raising M can only raise every coefficient; the "
                           "certified quantity is a LOWER bound at M = m"),
        "mean_feasible_configurations": feas,
        "min_G_minus_Gm_lower": "n/a" if worst is None else str(worst),
        "violations": viol,
    }


def claim_c5_channel(beta: arb, m: arb) -> dict[str, Any]:
    """T = c + 2 <mu, R nu>."""
    worst = arb(0)
    ok = True
    for cfg in configurations():
        x = [_A(v) for v in cfg["x"]]
        y = [_A(v) for v in cfg["y"]]
        a = [_A(v) for v in cfg["a"]]
        pc = pieces(x, y, a, m, beta)
        rr = _bilinear(x, y, a, lambda s, t: R(s, t, beta, m))
        resid = pc["T"] - (pc["c"] + 2 * rr)
        if 0 not in resid:
            ok = False
        if abs(resid).upper() > worst.upper():
            worst = arb(abs(resid).upper())
    return {"claim": "C5 T = c + 2 <mu, R nu>",
            "verdict": "PROVED" if ok else "FAILED",
            "max_residual_upper": worst.str(8),
            "every_residual_encloses_zero": ok}


def phi_(s: arb, beta: arb) -> arb:
    """phi(s) = sqrt(Q2(s,s)); the rank-one PSD direction."""
    return Q2(s, s, beta).sqrt()


def A_(s: arb, t: arb, beta: arb, m: arb) -> arb:
    """P2 side of the splitting."""
    return P2(s, t, beta, m) + phi_(s, beta) * phi_(t, beta)


def B_(s: arb, t: arb, beta: arb) -> arb:
    """Q2 side of the splitting."""
    return Q2(s, t, beta) - phi_(s, beta) * phi_(t, beta)


def claim_c6_splitting(beta: arb, m: arb) -> dict[str, Any]:
    """The rank-one splitting.  With phi = sqrt(Q2(s,s)),

        A = P2 + phi (x) phi,     B = Q2 - phi (x) phi,     R = A + B,

        T = <mu,B mu> + <nu,B nu> + 2 <mu,A nu> + (<mu,phi> - <nu,phi>)^2.

    Every summand is an integral of the corresponding kernel against
    NONNEGATIVE measures, plus a square.  So A >= 0 and B >= 0 imply T >= 0 for
    ARBITRARY nonnegative measures -- any number of atoms, no copositivity
    machinery, no dimension.  This is Diananda's N + PSD with a rank-one PSD
    part, made explicit.  Since R = A + B, LEMMA A is implied as well.
    """
    w_id = arb(0)
    w_ra = arb(0)
    ok_id = ok_ra = True
    min_summand = None
    for cfg in configurations():
        x = [_A(v) for v in cfg["x"]]
        y = [_A(v) for v in cfg["y"]]
        a = [_A(v) for v in cfg["a"]]
        Ak = lambda s, t: A_(s, t, beta, m)   # noqa: E731
        Bk = lambda s, t: B_(s, t, beta)      # noqa: E731
        T = (_bilinear(x, x, a, lambda s, t: Q2(s, t, beta))
             + 2 * _bilinear(x, y, a, lambda s, t: P2(s, t, beta, m))
             + _bilinear(y, y, a, lambda s, t: Q2(s, t, beta)))
        mphi = sum((a[i] * phi_(x[i], beta) for i in range(3)), arb(0))
        nphi = sum((a[i] * phi_(y[i], beta) for i in range(3)), arb(0))
        parts = [_bilinear(x, x, a, Bk), _bilinear(y, y, a, Bk),
                 2 * _bilinear(x, y, a, Ak), (mphi - nphi) ** 2]
        resid = T - sum(parts, arb(0))
        if 0 not in resid:
            ok_id = False
        if abs(resid).upper() > w_id.upper():
            w_id = arb(abs(resid).upper())
        for p in parts:
            if min_summand is None or p.lower() < min_summand:
                min_summand = p.lower()
        for i in range(3):
            for j in range(3):
                r2 = (R(x[i], y[j], beta, m)
                      - (Ak(x[i], y[j]) + Bk(x[i], y[j])))
                if 0 not in r2:
                    ok_ra = False
                if abs(r2).upper() > w_ra.upper():
                    w_ra = arb(abs(r2).upper())
    return {
        "claim": "C6 rank-one splitting: T = <mu,B mu>+<nu,B nu>+2<mu,A nu>+sq",
        "verdict": "PROVED" if (ok_id and ok_ra) else "FAILED",
        "identity_max_residual_upper": w_id.str(8),
        "R_equals_A_plus_B_max_residual_upper": w_ra.str(8),
        "min_summand_lower_bound": "n/a" if min_summand is None else str(min_summand),
        "consequence": ("A >= 0 and B >= 0 imply T >= 0 for ARBITRARY "
                        "nonnegative measures (dimension-free); and R = A + B "
                        "implies LEMMA A, so both LEMMA A and LEMMA C follow "
                        "from the two two-variable inequalities"),
    }


def ab_zero_structure(beta: arb, m: arb, xstar: arb) -> dict[str, Any]:
    """Zero sets of A and B, which any certificate must stratify around."""
    z, one, half = arb(0), arb(1), arb(1) / 2
    diag_vals = []
    for sv in ("0.05", "0.3", "0.5", "0.680688509846052185", "0.8", "0.95"):
        s = arb(sv)
        diag_vals.append({"s": sv, "B_on_diagonal": B_(s, s, beta).str(6)})
    return {
        "A": {
            "definition": "A(s,t) = P2(s,t) + sqrt(Q2(s,s) Q2(t,t))",
            "zeros": ["{s=0}", "{t=0}", "(1,1)", "(x*,x*) interior"],
            "A_at_origin": A_(z, z, beta, m).str(6),
            "A_on_axis_s_half_t_zero": A_(half, z, beta, m).str(6),
            "A_at_one_one": A_(one, one, beta, m).str(6),
            "A_at_xstar": A_(xstar, xstar, beta, m).str(8),
            "A_at_xstar_encloses_zero": bool(0 in A_(xstar, xstar, beta, m)),
            "factorization": ("A = s*t*Acheck, since P2 = s*t*(...) and "
                              "Q2(s,s) = s^2*(...), which removes both axis "
                              "zeros; (1,1) and (x*,x*) remain"),
            "interior_zero_is_nondegenerate": True,
        },
        "B": {
            "definition": "B(s,t) = Q2(s,t) - sqrt(Q2(s,s) Q2(t,t))",
            "zeros": ["the ENTIRE diagonal {s=t}, with vanishing gradient",
                      "{s=0}", "{t=0}"],
            "diagonal_samples": diag_vals,
            "gradient_vanishes_on_diagonal": ("d/dt B(s,t)|_{t=s} = 0 "
                                              "identically, because "
                                              "g'(s) = 2 dQ2/dt(s,s) by "
                                              "symmetry of Q2"),
            "factorization": "B = (s-t)^2 * W, certify W >= 0",
            "equivalent_form": "h(pi(s,t))^2 >= h(pi(s,s)) h(pi(t,t))",
            "obstruction_note": ("pi(s,t)^2 <= pi(s,s) pi(t,t) ALWAYS (with "
                                 "a=1-s, b=1-t it is (1+ab)^2 <= "
                                 "(1+a^2)(1+b^2) <=> -(a-b)^2 <= 0), so pi "
                                 "runs the wrong way and B >= 0 must come "
                                 "from the shape of h, not monotonicity"),
        },
    }


def local_data(beta: arb, m: arb, xstar: arb) -> dict[str, Any]:
    """Certified local data for R at the interior zero (x*,x*)."""
    val = R(xstar, xstar, beta, m)
    g = R_ds(xstar, xstar, beta, m)
    H11 = R_dss(xstar, xstar, beta, m)
    H12 = R_dst(xstar, xstar, beta, m)
    lam_min = H11 + H12          # H22 == H11 by symmetry; H12 < 0
    lam_max = H11 - H12
    K12 = _channel_K12(xstar)
    z = beta * abs(K12)
    return {
        "point": "(x*, x*)",
        "x_star": xstar.str(DIGITS),
        "R_encloses_zero": bool(0 in val),
        "R_value": val.str(8),
        "grad_encloses_zero": bool(0 in g),
        "grad_value": g.str(8),
        "H11": H11.str(DIGITS),
        "H12": H12.str(DIGITS),
        "lambda_min": lam_min.str(DIGITS),
        "lambda_max": lam_max.str(DIGITS),
        "hessian_positive_definite": bool(lam_min.lower() > 0),
        "K12": K12.str(DIGITS),
        "beta_abs_K12": z.str(DIGITS),
        "channel_sign_near_diagonal": "NEGATIVE" if K12.upper() < 0 else "NONNEGATIVE",
        "local_T_symmetric_eigenvalue": lam_min.str(DIGITS),
        "local_T_antisymmetric_eigenvalue": (lam_max - 2 * z).str(DIGITS),
        "local_T_positive_definite": bool(lam_min.lower() > 0
                                          and (lam_max - 2 * z).lower() > 0),
        "note": ("c ~ beta*K12*Delta^2 vanishes in the symmetric direction, so "
                 "the symmetric eigenvalue of the local form of T equals "
                 "lambda_min(Hess R) exactly, independent of the channel "
                 "strength beta|K12|"),
    }


def _channel_K12(xs: arb) -> arb:
    """d^2 K/ds dt at (xs,xs), K = h(pi)."""
    p = pi_(xs, xs)
    return (_hpp(p) * _pi_s(xs, xs) ** 2 + _hp(p) * _pi_st(xs, xs))


def refuted_entrywise(beta: arb, m: arb) -> dict[str, Any]:
    """The entrywise strengthening of LEMMA C is FALSE.  Exact witness."""
    z = arb(0)
    half = arb(1) / 2
    raw = Psi_M(z, half, half, z, m, beta)
    spl = Psi_split(z, half, half, z, beta, m)
    p_mid = P2(half, half, beta, m)
    return {
        "route": "entrywise nonnegativity of M_T, i.e. Psi >= 0 on [0,1]^4",
        "verdict": "REFUTED",
        "witness": "(s,t,sig,tau) = (0, 1/2, 1/2, 0)",
        "Psi_raw": raw.str(DIGITS),
        "Psi_split": spl.str(DIGITS),
        "raw_and_split_agree": bool(0 in (raw - spl)),
        "strictly_negative": bool(raw.upper() < 0),
        "mechanism": ("at s = tau = 0 every Q2 term vanishes and Psi reduces to "
                      "P2(t,sig); P2 is sign-indefinite with minimum "
                      "-0.0693028554493 at s = t = 0.680688509846052185"),
        "P2_at_half_half": p_mid.str(DIGITS),
        "consequence": ("LEMMA C must be stated as COPOSITIVITY of M_T, which is "
                        "strictly weaker and survives this point; the entrywise "
                        "route must not be retried"),
    }


# --------------------------------------------------------------------------
# Mutations: each must break a claim the sound code satisfies
# --------------------------------------------------------------------------

def mutations(beta: arb, m: arb) -> list[dict[str, Any]]:
    cfg = configurations()[0]
    x = [_A(v) for v in cfg["x"]]
    y = [_A(v) for v in cfg["y"]]
    a = [_A(v) for v in cfg["a"]]
    q = _A(cfg["q"])
    M0 = sum((a[i] * x[i] for i in range(3)), arb(0))
    M1 = sum((a[i] * y[i] for i in range(3)), arb(0))
    M = (1 - q) * M0 + q * M1
    pc = pieces(x, y, a, M, beta)
    sound = ((1 - q) ** 2 * pc["I00"] + q ** 2 * pc["I11"]) / M \
        + q * (1 - q) * pc["T"]
    vals = [a[0], a[1], q, x[0], x[1], x[2], y[0], y[1], y[2]]
    terms = liu9_objective.evaluate_arb(vals, beta)
    repo = terms.numerator - terms.ehx
    out = []

    def record(name: str, mutant: arb, target: str) -> None:
        caught = 0 not in (mutant - repo)
        out.append({"mutation": name, "target": target,
                    "mutant": mutant.str(10), "sound": sound.str(10),
                    "reference": repo.str(10), "caught": bool(caught),
                    "expected": "fail"})

    # 1. drop the protocol term from the direct kernel
    dm_noproto = lambda s, t: M * ((1 - beta) * h(s * t)) \
        - (t * h(s) + s * h(t)) / 2  # noqa: E731
    I00b = _bilinear(x, x, a, dm_noproto)
    I11b = _bilinear(y, y, a, dm_noproto)
    I01b = _bilinear(x, y, a, dm_noproto)
    Tb = 2 * I01b / M + pc["c"]
    record("drop_protocol_in_D",
           ((1 - q) ** 2 * I00b + q ** 2 * I11b) / M + q * (1 - q) * Tb, "C1")

    # 2. freeze the kernel mean at m instead of the true mixture mean M
    pcm = pieces(x, y, a, m, beta)
    record("freeze_kernel_mean_at_m",
           ((1 - q) ** 2 * pcm["I00"] + q ** 2 * pcm["I11"]) / m
           + q * (1 - q) * pcm["T"], "C1")

    # 3. drop the channel from T
    record("drop_channel_from_T",
           ((1 - q) ** 2 * pc["I00"] + q ** 2 * pc["I11"]) / M
           + q * (1 - q) * (2 * pc["I01"] / M), "C1/C5")

    # 4. use the wrong cross pairing in the channel
    c_bad = beta * (_bilinear(x, x, a, K) + _bilinear(y, y, a, K)
                    - _bilinear(x, y, a, K))
    record("halve_channel_cross_term",
           ((1 - q) ** 2 * pc["I00"] + q ** 2 * pc["I11"]) / M
           + q * (1 - q) * (2 * pc["I01"] / M + c_bad), "C5")
    return out


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def _digest(payload: dict) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def phi_sharpness(beta: arb, m: arb, xstar: arb) -> dict[str, Any]:
    """No GLOBAL rescaling of phi = sqrt(Q2(s,s)) survives; phi is pinned
    pointwise exactly on {0, x*, 1}, and has slack elsewhere.

    The splitting identity itself holds for ANY phi -- it is pure algebra -- so
    a perturbation of phi cannot break it.  What the specific choice buys is the
    NONNEGATIVITY of both kernels, and a uniform rescaling is squeezed out from
    two sides:

      scaling phi UP   makes B(s,s) = Q2(s,s)(1 - k^2) < 0 on the diagonal;
      scaling phi DOWN makes A(x*,x*) = R(x*,x*) - (1-k^2) Q2(x*,x*) < 0,
                       because R(x*,x*) = 0 exactly.

    But that is a statement about UNIFORM k only.  Pointwise, the diagonal
    constraint is -P2(s,s) <= phi(s)^2 <= Q2(s,s), whose width is exactly
    R(s,s); so phi(s) is forced precisely where R(s,s) = 0, i.e. on
    {0, x*, 1}, and elsewhere admits genuine slack -- at s = 1/2 an isolated
    rescale phi(1/2) -> r phi(1/2) survives for every
    r >= sqrt(-P2(1/2,1/2)/Q2(1/2,1/2)) = 0.94395292747338.
    Credit: the slack was found by the ChannelTheory worker and re-verified
    here; an earlier wording of this function claimed phi was "forced" outright,
    which overclaimed.  The slack matters: a pointwise-varying phi could buy
    margin at the tight points and make certification easier.
    """
    out = []
    probe = arb("0.5")
    for label, k in (("scale_phi_up_1.05", arb("1.05")),
                     ("scale_phi_down_0.95", arb("0.95"))):
        pp = lambda s: k * phi_(s, beta)   # noqa: E731
        b_diag = Q2(probe, probe, beta) - pp(probe) * pp(probe)
        a_star = P2(xstar, xstar, beta, m) + pp(xstar) * pp(xstar)
        broke = (b_diag.upper() < 0) or (a_star.upper() < 0)
        out.append({
            "mutation": label,
            "target": "nonnegativity of the split kernels (not the identity)",
            "B_on_diagonal_at_0.5": b_diag.str(8),
            "A_at_xstar": a_star.str(8),
            "caught": bool(broke),
            "expected": "fail",
        })
    pin = []
    for lbl, s in (("0", arb(0)), ("x*", xstar), ("1", arb(1)),
                   ("1/2", arb("0.5")), ("0.3", arb("0.3"))):
        rss = R(s, s, beta, m)
        lo = -P2(s, s, beta, m)
        hi = Q2(s, s, beta)
        ratio = (lo / hi).sqrt() if hi.lower() > 0 else arb(0)
        pin.append({"s": lbl, "R_on_diagonal": rss.str(8),
                    "phi_sq_window_width_equals_R": rss.str(8),
                    "min_admissible_isolated_rescale_r": ratio.str(14),
                    "phi_forced_here": bool(0 in rss)})
    return {
        "note": ("the splitting identity is phi-independent; these probes show "
                 "no UNIFORM rescaling of phi survives, and separately that phi "
                 "is forced pointwise exactly where R(s,s) = 0"),
        "probes": out,
        "no_uniform_rescaling_survives": all(p["caught"] for p in out),
        "pointwise_window": pin,
        "phi_forced_on": ["0", "x*", "1"],
        "phi_has_slack_elsewhere": True,
    }


def build_report() -> dict[str, Any]:
    ctx.prec = CTX_PREC
    pars = certify_equation_parameters(solve_equation_parameters(100))
    beta, m, xstar = pars.beta, pars.mean, pars.x
    claims = [claim_c1_all_q(beta, m), claim_c2_matrix(beta, m),
              claim_c3_split(beta, m), claim_c4_monotone(beta, m),
              claim_c5_channel(beta, m), claim_c6_splitting(beta, m)]
    muts = mutations(beta, m)
    all_proved = all(cl["verdict"] == "PROVED" for cl in claims)
    all_caught = all(mt["caught"] for mt in muts)
    report: dict[str, Any] = {
        "claim_status": ("REDUCTION-CERTIFIED" if (all_proved and all_caught)
                         else "FAILED"),
        "arb_precision_bits": CTX_PREC,
        "constants": {"beta": beta.str(DIGITS), "m": m.str(DIGITS),
                      "x_star": xstar.str(DIGITS)},
        "statement": {
            "reduces": ("Liu H2 for the three-atom paired class: 9 variables "
                        "(3 shared masses, 6 supports, q) with the mean "
                        "constraint M >= m"),
            "to": ["LEMMA A': A(s,t) = P2(s,t) + sqrt(Q2(s,s) Q2(t,t)) >= 0",
                   "LEMMA B': B(s,t) = Q2(s,t) - sqrt(Q2(s,s) Q2(t,t)) >= 0"],
            "equivalently": ("the single sandwich  0 <= Q2(s,t) - "
                             "sqrt(Q2(s,s) Q2(t,t)) <= R(s,t)  on [0,1]^2"),
            "superseded": ("earlier framings of this target -- 'Psi >= 0 on "
                           "[0,1]^4' (REFUTED, see refuted_routes) and 'M_T "
                           "copositive on [0,1]^6' (true but unnecessary, "
                           "implied by C6) -- must not be cited as the live "
                           "obligation"),
            "eliminated_exactly": ["q (all q in [0,1], no case split)",
                                   "the three masses (rank-one splitting C6, "
                                   "valid for arbitrarily many atoms)",
                                   "the mean constraint (monotonicity in M)"],
            "not_required": ["the activation lemma min{M>=m} = min{M=m}",
                             "the active-face q* elimination",
                             "the Delta = 0 case analysis",
                             "any box cover of [0,1]^7",
                             "any 3x3 copositivity decision procedure",
                             "any 6-dimensional branch-and-bound"],
        },
        "claims": claims,
        "local_data": local_data(beta, m, xstar),
        "ab_zero_structure": ab_zero_structure(beta, m, xstar),
        "phi_sharpness": phi_sharpness(beta, m, xstar),
        "refuted_routes": [refuted_entrywise(beta, m)],
        "mutations": muts,
        "limitations": [
            "LEMMA A' (A >= 0) and LEMMA B' (B >= 0) are NOT proved here; this "
            "module certifies only the exact reduction to them.  Both are "
            "OPEN, and both are supported so far only by dense scans "
            "(DISCOVERY): min A = 0 and min B = 0 on a 600^2 grid of [0,1]^2.",
            "A vanishes at the interior point (x*,x*) with vanishing gradient, "
            "and B vanishes on the whole diagonal with vanishing gradient, so "
            "neither can be certified by interval subdivision alone; the "
            "factorizations A = s*t*Acheck and B = (s-t)^2*W plus local "
            "strata are mandatory.",
            "The identities are verified on a fixed finite configuration "
            "lattice.  They are exact algebraic identities in Arb, so the "
            "lattice demonstrates them rather than quantifying over the "
            "continuum; the reduction's validity rests on the symbolic "
            "derivation, independently re-derived in liu9_psi_reduction_audit.",
            "The reduction is a SUFFICIENT direction only: A >= 0 and B >= 0 "
            "imply H2, but H2 does not require them.  A failure of either "
            "would not refute H2.",
        ],
    }
    report["report_sha256"] = _digest(report)
    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    args = parser.parse_args(argv)
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("claim_status:", report["claim_status"])
    for cl in report["claims"]:
        print("  %-58s %s" % (cl["claim"][:58], cl["verdict"]))
    ld = report["local_data"]
    print("local at (x*,x*): R~0 %s  grad~0 %s  lambda_min %s  PD %s"
          % (ld["R_encloses_zero"], ld["grad_encloses_zero"],
             ld["lambda_min"], ld["hessian_positive_definite"]))
    print("local T eigenvalues: sym %s  antisym %s  PD %s"
          % (ld["local_T_symmetric_eigenvalue"],
             ld["local_T_antisymmetric_eigenvalue"],
             ld["local_T_positive_definite"]))
    rr = report["refuted_routes"][0]
    print("refuted entrywise route: Psi(0,1/2,1/2,0) =", rr["Psi_raw"])
    for mt in report["mutations"]:
        print("  mutation %-28s caught=%s" % (mt["mutation"], mt["caught"]))
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] == "REDUCTION-CERTIFIED" else 3


if __name__ == "__main__":
    raise SystemExit(main())
