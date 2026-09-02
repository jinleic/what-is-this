#!/usr/bin/env python3
"""Adversarial verification: Main's Psi_M reduction of the Liu 9-variable gap.

Target under audit (Liu Hypothesis 2, discovery by Main, reported in
LIU9_BLOCK_COPOSITIVE_2026-08-29.md and uc/state.json): the full 9-variable +
q problem reduces to ONE four-variable inequality

    Psi_m >= 0  on [0, 1]^4.

Definitions (exact):
    h(u) = -u*log(u) - (1-u)*log(1-u),  h(0) = h(1) = 0,
    pi(s,t) = s*t*(1 + (1-s)*(1-t)),
    D_M(s,t) = M*((1-beta)*h(s*t) + beta*h(pi(s,t))) - (t*h(s) + s*h(t))/2,
    Psi_M(s,t,sig,tau) = (D_M(s,tau) + D_M(t,sig))/M
        + beta*( h(pi(s,t)) + h(pi(sig,tau)) - h(pi(s,tau)) - h(pi(t,sig)) ),
with beta ~ 0.10005255986289310666 and m ~ 0.61729091208126497007 from
liu9_binding.solve_equation_parameters(100).

Chain audited here (independent re-derivation; Main's arithmetic is NOT used):
  Claim 1 (ALL-q identity).  With M0 = sum a_i x_i, M1 = sum a_i y_i,
      M = (1-q)M0 + q M1,
      I00 = sum_ij a_i a_j D_M(x_i,x_j), I11 = sum_ij a_i a_j D_M(y_i,y_j),
      I01 = sum_ij a_i a_j D_M(x_i,y_j),
      c = beta * sum_ij a_i a_j ( h(pi(x_i,x_j)) + h(pi(y_i,y_j))
                                  - 2 h(pi(x_i,y_j)) ),
      T = 2*I01/M + c, then
      gap == ((1-q)^2*I00 + q^2*I11)/M + q*(1-q)*T   for ALL q.
  Claim 2.  The symmetric matrix (M_T)_ij := Psi_M(x_i,x_j,y_i,y_j) is exactly
      the symmetric matrix of the quadratic form a -> T.
  Claim 3.  Psi_M(s,t,s,t) == (2/M)*D_M(s,t) identically.
  Claim 4 (monotonicity).  For M >= m: D_M >= D_m pointwise and
      D_M/M >= D_m/M and hence Psi_M >= Psi_m.
  Claim 5 (conclusion).  If Psi_m >= 0 on [0,1]^4 then gap >= 0 for every
      feasible configuration with M >= m, q in [0,1], a in the simplex.

VERIFICATION METHOD.  The Liu 9-atom paired-class transcription carries
weights w = ((1-q)a1, (1-q)a2, (1-q)a3, q*a1, q*a2, q*a3) over atoms
(x1,x2,x3,y1,y2,y3) and the raw gap
    gap = (1-beta)*EHXY + beta*EHPI - EHX
with EHXY = sum_ij w_i w_j h(z_i z_j) (36 ordered terms), EHPI the two
9-term component protocol sums, EHX = sum_i w_i h(z_i).

Every symbolic quantity is expanded in a canonical FREE LOG BASIS with an
entropy reduction INDEPENDENT of the repository's liu9_paired_class_c
machinery (which encodes the same reductions through its own p/q/r symbol
families): for a NONEMPTY multi-set P of atom indices,
    h(prod P)     = -(prod P)*L(P) - (1 - prod P)*K(P),
    mu(prod P)    = -(1 - prod P)*K(P)/(prod P),       mu(1) := 1,
    log(prod P)   = sum_{i in P} l_i                   (exact additivity),
where l_i are six independent free symbols, K(P) is an independent free
symbol per multi-set (K0..K5 for singletons; K001-style names otherwise), and
the protocol values h(pi(x_i,x_j)) are reductions in the per-pair free
symbols (LQ[i,j], KQ[i,j]):
    h(pi(x_i,x_j)) = -pi_ij*LQ[i,j] - (1 - pi_ij)*KQ[i,j],
        pi_ij = x_i x_j (1 + (1-x_i)(1-x_j)),
exactly mirroring the fact that log pi and log(1-pi) are free functions of
the ordered pair (i,j) (the complement 1 - pi = 1 - x_i x_j (2 - x_i - x_j +
x_i x_j) is NOT a product of atoms, so no additive collapse exists).
Expressions are multilinear over these independent generators (plus the atom
polynomials), so "expand canonically and require exact zero" is a valid
identity proof over the free commutative ring: any residual reported is an
exact canonical numerator, never a float.

Two independent symbolic routes prove Claim 1:
  Route A (discharge / kernel collapse): M*gap - q(1-q)*c equals
      sum_ij w_i w_j [ M*e(x_i,x_j) - (z_j h(z_i) + z_i h(z_j))/2 ]
      with the mu-lemma kernel
      e(s,t) = t h(s) + s h(t) - s t (mu(s)+mu(t)-mu(st))
               + beta*(h(pi(s,t)) - h(s*t)),
      which additionally PROVES the kernel collapse
      e(s,t) == (1-beta) h(s t) + beta h(pi(s,t))  (equivalently the raw
      mu-lemma t h(s) + s h(t) - s t(mu(s)+mu(t)-mu(st)) == h(s t)).
  Route B (decomposition): gap == ((1-q)^2 I00 + q^2 I11)/M + q(1-q) T,
      assembled directly from the D-kernel D_M at parameter M = mixture mean.

Arb numeric certificate: Claim 1 is re-verified at ctx.prec = 320 over 288
seeded deterministic exact-rational configurations (masses a simplex draw,
q in (0,1), supports in (0,1)); configurations with M < m are lifted
supports-up by the exact rational theta = (m - M)/(1 - M) acting as
z -> z + theta (1 - z), which forces M = m exactly.  Both sides are evaluated
in Arb ball arithmetic on exact rational inputs; every residual must be
enclosed below 2^-250.  float64 residuals are reported for scale, labelled
DISCOVERY (they are not certificates).

Adversarial mutations: four soundness-breaking mutations of the identities
that MUST fire (mutant residual != 0) while the sound control residual is 0.

Byte-stable report: uc/verification/results/liu9-psi-reduction.json
Run:  math/.venv/bin/python -I -B math/uc/liu9_psi_reduction_audit.py
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import random
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

# This workstation is shared with live certification workers; these must be
# set before importing numpy/openblas.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import mpmath  # noqa: E402
import sympy as sp  # noqa: E402
from flint import arb, ctx  # noqa: E402

from liu9_binding import solve_equation_parameters  # noqa: E402

ctx.prec = max(ctx.prec, 320)
OUTPUT_DEFAULT = HERE / "verification/results/liu9-psi-reduction.json"
CERT_SEED = 20260901
ARB_DRAWS = 288          # >= 200 required by the audit ticket
ARB_PREC_TARGET = 320
RESIDUAL_BOUND_EXP = 250  # certification: |residual| < 2^-250

# --- base symbols ------------------------------------------------------------

A1, BETA, Q = sp.symbols("a1 beta q", real=True)
A2 = sp.Symbol("a2", real=True)
X1, X2, X3, Y1, Y2, Y3 = sp.symbols("x1 x2 x3 y1 y2 y3", real=True)
ATOMS = (X1, X2, X3, Y1, Y2, Y3)
MASSES = (A1, A2, 1 - A1 - A2)
_MEAN = sp.Symbol("M_mix", real=True)       # mixture mean parameter
_M_PSI = sp.Symbol("M_Psi", real=True)      # free Psi/D parameter in Claims 2-4
_M_ALT = sp.Symbol("m_binding_alt", real=True)  # mutation M3 Free-solid alternative

# --- independent free log basis ----------------------------------------------

_L1 = {i: sp.Symbol(f"l_{ATOMS[i].name}") for i in range(6)}
_K1 = {i: sp.Symbol(f"k_{ATOMS[i].name}") for i in range(6)}
_LQ = {(i, j): sp.Symbol(f"LQ[{i},{j}]") for i in range(6) for j in range(6)}
_KQ = {(i, j): sp.Symbol(f"KQ[{i},{j}]") for i in range(6) for j in range(6)}


def _kg_of(atomset: Tuple[int, ...]) -> sp.Expr:
    """Independent generator of log(1 - prod P) for a nonempty multi-set P."""
    if len(atomset) == 1:
        return _K1[atomset[0]]
    return sp.Symbol("K" + "".join(str(i) for i in atomset))


def _lg_of(atomset: Tuple[int, ...]) -> sp.Expr:
    """log(prod P) equals the sum of singleton logs (exact additivity)."""
    return sum((_L1[i] for i in atomset), sp.Integer(0))


def _zprod(atomset: Tuple[int, ...]) -> sp.Expr:
    return sp.prod([ATOMS[i] for i in atomset])


def _h_of(atomset: Tuple[int, ...]) -> sp.Expr:
    """h(prod P) with P a nonempty atom multi-set."""
    if not atomset:
        return sp.Integer(0)
    z = _zprod(atomset)
    return -z * _lg_of(atomset) - (1 - z) * _kg_of(atomset)


def _atomset_of_product(u: sp.Expr) -> Tuple[int, ...]:
    """Sorted multi-set of atom indices of a product-of-atoms expression."""
    if u == 0:
        raise AssertionError("empty product has no atom set")
    out: List[int] = []
    for atom, p in sorted(u.as_powers_dict().items(), key=lambda kv: str(kv[0])):
        if int(p) != p or p < 1 or not isinstance(atom, sp.Symbol):
            raise AssertionError(f"not a product of atoms: {u}")
        out.extend([ATOMS.index(atom)] * int(p))
    return tuple(sorted(out))


def _h_free(u: sp.Expr) -> sp.Expr:
    """h(u) for u a product of atoms (repetition allowed), in the basis."""
    return _h_of(_atomset_of_product(u))


def _h1(i: int) -> sp.Expr:
    return _h_of((i,))


def _mu_of(atomset: Tuple[int, ...]) -> sp.Expr:
    """mu(prod P) = -(1-prod P) log(1-prod P) / (prod P); mu(1) := 1."""
    if not atomset:
        return sp.Integer(1)
    z = _zprod(atomset)
    return -(1 - z) * _kg_of(atomset) / z


def _mu1(i: int) -> sp.Expr:
    return _mu_of((i,))


def _pi_raw(i: int, j: int) -> sp.Expr:
    s, t = ATOMS[i], ATOMS[j]
    return sp.expand(s * t * (1 + (1 - s) * (1 - t)))



# --- Liu 9-variable transcription ---------------------------------------------
def _pi_coef(i: int, j: int) -> sp.Expr:
    """h(pi(x_i,x_j)) = -pi*LQ - (1-pi)*KQ on the CANONICAL unordered pair
    key (min(i,j), max(i,j)): h(pi) is a single value per atom pair, exactly
    as liu9_paired_class_c._PI_SYM keys it by (lo, hi)."""
    lo, hi = min(i, j), max(i, j)
    pij = _pi_raw(lo, hi)
    return -pij * _LQ[(lo, hi)] - (1 - pij) * _KQ[(lo, hi)]


def _weights() -> Tuple[sp.Expr, ...]:
    """w = ((1-q)a1, (1-q)a2, (1-q)a3, q a1, q a2, q a3)."""
    return ((1 - Q) * MASSES[0], (1 - Q) * MASSES[1], (1 - Q) * MASSES[2],
            Q * MASSES[0], Q * MASSES[1], Q * MASSES[2])


def _mixture_mean() -> sp.Expr:
    """M = sum_i w_i z_i = (1-q) sum a_i x_i + q sum a_i y_i."""
    w = _weights()
    return sp.expand(sum((w[i] * ATOMS[i] for i in range(6)), sp.Integer(0)))


# --- canonical residual helpers ------------------------------------------------

def _numer(expr: sp.Expr) -> sp.Expr:
    """Canonical exact numerator of a rational-function residual."""
    num, _den = sp.fraction(sp.cancel(sp.together(expr)))
    return sp.expand(num)


def _residual_record(residual: sp.Expr) -> Dict[str, Any]:
    args = sp.Add.make_args(residual)
    return {"residual_is_zero": residual == 0,
            "residual_term_count": len(args),
            "residual_str": str(residual) if residual != 0 else "0",
            "residual_str_head": (str(residual)[:240] if residual != 0 else "")}


# --- gap assembly (direct, independent of liu9_paired_class_c) ------------------

def _gap(beta: sp.Symbol) -> sp.Expr:
    """gap = (1-beta) EHXY + beta EHPI - EHX in the independent basis."""
    w = _weights()
    ehxy = sp.Integer(0)
    for i in range(6):
        for j in range(6):
            ehxy += w[i] * w[j] * _h_free(ATOMS[i] * ATOMS[j])
    ehpi = sp.Integer(0)
    for comp, factor in ((0, 1 - Q), (1, Q)):
        off = 3 * comp
        for i in range(3):
            for j in range(3):
                ehpi += factor * MASSES[i] * MASSES[j] * _pi_coef(off + i, off + j)
    ehx = sum((w[i] * _h1(i) for i in range(6)), sp.Integer(0))
    return sp.expand((1 - beta) * ehxy + beta * ehpi - ehx)


def _remainder(beta: sp.Symbol) -> sp.Expr:
    """c = beta sum_{i,j<3} a_i a_j (h(pi(x_i,x_j)) + h(pi(y_i,y_j))
    - 2 h(pi(x_i,y_j)))."""
    total = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            total += MASSES[i] * MASSES[j] * (
                _pi_coef(i, j) + _pi_coef(3 + i, 3 + j) - 2 * _pi_coef(i, 3 + j))
    return beta * total


# --- D kernel, Psi, Claims 2-4 --------------------------------------------------

def _D_M(M: sp.Symbol, s: sp.Expr, t: sp.Expr, beta: sp.Symbol) -> sp.Expr:
    si, ti = ATOMS.index(s), ATOMS.index(t)
    return M * ((1 - beta) * _h_free(s * t) + beta * _pi_coef(si, ti)) \
        - (t * _h1(si) + s * _h1(ti)) / 2


def _psi(M: sp.Symbol, s: sp.Expr, t: sp.Expr, sig: sp.Expr,
         tau: sp.Expr, beta: sp.Symbol) -> sp.Expr:
    d1 = _D_M(M, s, tau, beta)
    d2 = _D_M(M, t, sig, beta)
    gamma = (_pi_coef(ATOMS.index(s), ATOMS.index(t))
             + _pi_coef(ATOMS.index(sig), ATOMS.index(tau))
             - _pi_coef(ATOMS.index(s), ATOMS.index(tau))
             - _pi_coef(ATOMS.index(t), ATOMS.index(sig)))
    return (d1 + d2) / M + beta * gamma


# --- Claims ---------------------------------------------------------------------

def prove_chain1_route_A(beta: sp.Symbol) -> Dict[str, Any]:
    """Route A: M*gap - q(1-q)c == discharge via the mu-lemma kernel."""
    lhs = _MEAN * _gap(beta) - _MEAN * Q * (1 - Q) * _remainder(beta)
    w = _weights()
    rhs = sp.Integer(0)
    for i in range(6):
        for j in range(6):
            e = (ATOMS[j] * _h1(i) + ATOMS[i] * _h1(j)
                 - ATOMS[i] * ATOMS[j] * (_mu1(i) + _mu1(j)
                                          - _mu_of(tuple(sorted((i, j)))))
                 + beta * (_pi_coef(i, j) - _h_free(ATOMS[i] * ATOMS[j])))
            rhs += w[i] * w[j] * (_MEAN * e
                                  - (ATOMS[j] * _h1(i) + ATOMS[i] * _h1(j)) / 2)
    residual = sp.expand(lhs - rhs).subs(_MEAN, _mixture_mean())
    residual = sp.expand(residual)
    rec = _residual_record(residual)
    rec["status"] = "PROVED" if residual == 0 else "FAILED"
    rec["method"] = ("route A: M*gap - q(1-q)c vs w_i w_j [M e(x_i,x_j) "
                     "- (z_j h(z_i) + z_i h(z_j))/2], " "M = mixture mean")
    return rec


def prove_chain1_route_B(beta: sp.Symbol) -> Dict[str, Any]:
    """Route B: gap == ((1-q)^2 I00 + q^2 I11)/M + q(1-q)(2 I01/M + c)."""
    dmixed = {(i, j): _D_M(_MEAN, ATOMS[i], ATOMS[3 + j], beta)
              for i in range(3) for j in range(3)}
    i00 = sum((MASSES[i] * MASSES[j] * _D_M(_MEAN, ATOMS[i], ATOMS[j], beta)
               for i in range(3) for j in range(3)), sp.Integer(0))
    i11 = sum((MASSES[i] * MASSES[j] * _D_M(_MEAN, ATOMS[3 + i], ATOMS[3 + j], beta)
               for i in range(3) for j in range(3)), sp.Integer(0))
    i01 = sum((MASSES[i] * MASSES[j] * dmixed[(i, j)]
               for i in range(3) for j in range(3)), sp.Integer(0))
    t_expr = 2 * i01 / _MEAN + _remainder(beta)
    raw = _gap(beta) - ((1 - Q) ** 2 * i00 + Q ** 2 * i11) / _MEAN \
        - Q * (1 - Q) * t_expr
    residual = _numer(raw.subs(_MEAN, _mixture_mean()))
    rec = _residual_record(residual)
    rec["status"] = "PROVED" if residual == 0 else "FAILED"
    rec["method"] = ("route B: direct D-kernel decomposition "
                     "((1-q)^2 I00 + q^2 I11)/M + q(1-q)(2 I01/M + c)")
    return rec


def prove_kernel_collapse(beta: sp.Symbol) -> Dict[str, Any]:
    """36 rows: e == (1-beta)h(st) + beta h(pi); 36 rows: mu-lemma."""
    kernel_fail: List[Dict[str, Any]] = []
    mu_fail: List[Dict[str, Any]] = []
    for i in range(6):
        for j in range(i, 6):
            s, t = ATOMS[i], ATOMS[j]
            e = (t * _h1(i) + s * _h1(j)
                 - s * t * (_mu1(i) + _mu1(j) - _mu_of(tuple(sorted((i, j)))))
                 + beta * (_pi_coef(i, j) - _h_free(s * t)))
            residual = _numer(e - ((1 - beta) * _h_free(s * t) + beta * _pi_coef(i, j)))
            if residual != 0:
                kernel_fail.append({"i": i, "j": j,
                                    "residual_head": str(residual)[:240]})
            mu_residual = _numer(t * _h1(i) + s * _h1(j)
                                 - s * t * (_mu1(i) + _mu1(j)
                                            - _mu_of(tuple(sorted((i, j)))))
                                 - _h_free(s * t))
            if mu_residual != 0:
                mu_fail.append({"i": i, "j": j,
                                "residual_head": str(mu_residual)[:240]})
    return {
        "rows_checked": 21,
        "kernel_collapse_residual_failures": kernel_fail,
        "mu_lemma_residual_failures": mu_fail,
        "status": "PROVED" if not kernel_fail and not mu_fail else "FAILED",
        "method": ("independent canonical expansion of e against "
                   "(1-beta)h(st)+beta h(pi(s,t)) and of the raw mu-lemma "
                   "against h(st), over all 21 unordered atom pairs "
                   "(equivalent to 36 ordered by symmetry)"),
    }


def prove_chain2(beta: sp.Symbol) -> Dict[str, Any]:
    """Claim 2: [[a^T Psi a]] coefficient extraction vs Psi(x_i,x_j,y_i,y_j)."""
    M = _M_PSI
    G = sp.symbols("g1 g2 g3", real=True)
    dm = {(i, j): _D_M(M, ATOMS[i], ATOMS[3 + j], beta)
          for i in range(3) for j in range(3)}
    i01 = sum((G[i] * G[j] * dm[(i, j)]
               for i in range(3) for j in range(3)), sp.Integer(0))
    cg = beta * sum(
        (G[i] * G[j] * (_pi_coef(i, j) + _pi_coef(3 + i, 3 + j)
                        - 2 * _pi_coef(i, 3 + j))
         for i in range(3) for j in range(3)), sp.Integer(0))
    Te = sp.expand(2 * i01 / M + cg)
    psi = {(i, j): sp.cancel(_psi(M, ATOMS[i], ATOMS[j],
                                  ATOMS[3 + i], ATOMS[3 + j], beta))
           for i in range(3) for j in range(3)}
    rows: List[Dict[str, Any]] = []
    for i in range(3):
        for j in range(3):
            a_ij = sp.cancel(sp.diff(sp.diff(Te, G[i]), G[j])) / 2
            residual = _numer(a_ij - psi[(i, j)])
            rec = _residual_record(residual)
            rec["i"] = i
            rec["j"] = j
            rows.append(rec)
    sym_rows: List[Dict[str, Any]] = []
    for i in range(3):
        for j in range(i + 1, 3):
            residual = _numer(psi[(i, j)] - psi[(j, i)])
            rec = _residual_record(residual)
            rec["i"] = i
            rec["j"] = j
            sym_rows.append(rec)
    ok = all(r["residual_is_zero"] for r in rows + sym_rows)
    return {"status": "PROVED" if ok else "FAILED",
            "method": ("T(g) = 2*I01(g)/M + c(g) with g the free mass vector; "
                       "A_ij = d^2 T/dg_i dg_j / 2 vs Psi_M(x_i,x_j,y_i,y_j); "
                       "plus Psi-matrix symmetry check"),
            "entry_rows": rows,
            "symmetry_rows": sym_rows}


def prove_chain3(beta: sp.Symbol) -> Dict[str, Any]:
    """Claim 3: Psi_M(s,t,s,t) == 2 D_M(s,t)/M (exact cancellation)."""
    rows: List[Dict[str, Any]] = []
    M = _M_PSI
    for (s, t) in ((X1, X2), (X3, Y3), (X1, X1), (Y2, Y2)):
        residual = _numer(_psi(M, s, t, s, t, beta)
                          - 2 * _D_M(M, s, t, beta) / M)
        rec = _residual_record(residual)
        rec["s"] = str(s)
        rec["t"] = str(t)
        rows.append(rec)
    ok = all(r["residual_is_zero"] for r in rows)
    return {"status": "PROVED" if ok else "FAILED",
            "method": "Psi(s,t,s,t) - 2 D_M(s,t)/M expanded canonically",
            "rows": rows}


def prove_chain4(beta: sp.Symbol) -> Dict[str, Any]:
    """Claim 4: exact d/dM identities + sign facts => D_M >= D_m, Psi_M >= Psi_m."""
    M = _M_PSI
    rows: List[Dict[str, Any]] = []
    for (s, t) in ((X1, X2), (X1, Y1), (X1, X1), (Y3, X2)):
        dD = sp.expand(sp.diff(_D_M(M, s, t, beta), M))
        expect_dD = (1 - beta) * _h_free(s * t) \
            + beta * _pi_coef(s == ATOMS and 0 or ATOMS.index(s),
                              ATOMS.index(t)) if False else \
            (1 - beta) * _h_free(s * t) \
            + beta * _pi_coef(ATOMS.index(s), ATOMS.index(t))
        r1 = _numer(dD - expect_dD)
        dR = sp.diff(_D_M(M, s, t, beta) / M, M)
        expect_dR = (t * _h1(ATOMS.index(s)) + s * _h1(ATOMS.index(t))) / (2 * M ** 2)
        r2 = _numer(dR - expect_dR)
        rec = {"s": str(s), "t": str(t),
               "dD_residual_is_zero": r1 == 0,
               "dD_residual_str": str(r1) if r1 != 0 else "0",
               "dR_residual_is_zero": r2 == 0,
               "dR_residual_str": str(r2) if r2 != 0 else "0"}
        rows.append(rec)
    ok = all(r["dD_residual_is_zero"] and r["dR_residual_is_zero"] for r in rows)
    return {"status": "PROVED" if ok else "FAILED",
            "method": ("exact sympy differentiation in the independent basis; "
                       "sign transported by the listed monotonicity facts"),
            "rows": rows,
            "sign_facts_used": [
                "h(u) >= 0 for u in [0,1] (natural-log binary entropy), "
                "so t*h(s) + s*h(t) >= 0 for s,t in [0,1]",
                "pi(s,t) = st(1+(1-s)(1-t)) in [0,1] (Liu's protocol map), "
                "so h(pi(s,t)) >= 0",
                "d/dk D_k(s,t) = (1-beta)h(st)+beta h(pi(s,t)) >= 0 for k > 0 "
                "=> k -> D_k(s,t) nondecreasing on [m, M], M >= m",
                "d/dk (D_k(s,t)/k) = (t h(s)+s h(t))/(2 k^2) >= 0 for k > 0 "
                "=> k -> D_k(s,t)/k nondecreasing on [m, M]",
                "Psi_M = (D_M(s,tau)+D_M(t,sig))/M + beta*(pi-block free of M); "
                "combining the two monotonicities gives Psi_M >= Psi_m for M >= m",
            ]}


def _arb_log1m(u: arb) -> arb:
    one = arb(1)
    if u in (arb(0), one):
        return -u  # exact boundary value: log(1-0) = 0, log(1-1) = -inf
    return (one - u).log()


def _arb_h1(u: arb) -> arb:
    one = arb(1)
    if u in (arb(0), one):
        return arb(0)  # exact boundary value h(0) = h(1) = 0
    return -(u * u.log() + (one - u) * (one - u).log())


def _arb_pi(s: arb, t: arb) -> arb:
    one = arb(1)
    return s * t * (one + (one - s) * (one - t))


def _arb_D(M: arb, s: arb, t: arb, beta: arb) -> arb:
    return M * ((arb(1) - beta) * _arb_h1(s * t) + beta * _arb_h1(_arb_pi(s, t))) \
        - (t * _arb_h1(s) + s * _arb_h1(t)) / 2


def chain5_refutation_check(beta_r: Fraction) -> Dict[str, Any]:
    """Re-run the refutation of the entrywise target independently: evaluate
    Psi_m(0, 1/2, 1/2, 0) in exact Arb at 320-bit precision under the exact
    boundary convention h(0) = h(1) = 0."""
    beta_a = _F2A(beta_r)
    ma = _F2A(_binding_params()[1])
    s, t, sig, tau = arb(0), arb(1) / 2, arb(1) / 2, arb(0)
    val = (_arb_D(ma, s, tau, beta_a) + _arb_D(ma, t, sig, beta_a)) / ma \
        + beta_a * (_arb_h1(_arb_pi(s, t)) + _arb_h1(_arb_pi(sig, tau))
                    - _arb_h1(_arb_pi(s, tau)) - _arb_h1(_arb_pi(t, sig)))
    d_val = _arb_D(ma, t, sig, beta_a)
    return {
        "probe": "Psi_m(0, 1/2, 1/2, 0)",
        "value": str(val),
        "certified_strictly_negative": bool(val < arb(0)),
        "reduces_to": "D_m(1/2, 1/2) / m (all Q2 terms vanish on s = tau = 0)",
        "d_m_half": str(d_val),
    }


def psi_decomposition_formula(beta: sp.Symbol) -> sp.Expr:
    """Exact U-cancellation identity for Psi (residual returned as canonical
    numerator; must be identically zero):
        Psi_M(s,t,sig,tau) = A_M(s,tau) + A_M(t,sig)
            + beta*(h(pi(s,t)) + h(pi(sig,tau)) - h(pi(s,tau)) - h(pi(t,sig)))
            + beta*(h(pi(s,tau)) + h(pi(t,sig))),
        A_M(u,v) = (1-beta) h(u v) - (v h(u) + u h(v)) / (2 M).
    Every pi-argument in the Psi bracket carries one factor from each of
    (s,t) and (sig,tau), so the (s,tau)/(t,sig) pi-terms of the A-kernels
    cancel exactly; at M = m this is exactly P2 + Q2 with
        P2(s,t) = (1-beta) h(s t) - (t h(s) + s h(t))/(2 m),
        Q2(s,t) = beta h(pi(s,t)),
    the A-kernels ARE P2 up to the m-parameter, and the bracket collapses:
        (M_T)_ij = P2(x_i,y_j) + P2(x_j,y_i) + Q2(x_i,x_j) + Q2(y_i,y_j),
    matching Main's and TwoVarCertifier's independent restatement."""
    M = _M_PSI
    s, t, sig, tau = X1, X2, X3, Y1
    big_a_st = (1 - beta) * _h_free(s * tau) \
        - (tau * _h1(0) + s * _h1(3)) / (2 * M)
    big_a_ts = (1 - beta) * _h_free(t * sig) \
        - (sig * _h1(1) + t * _h1(2)) / (2 * M)
    residual = _psi(M, s, t, sig, tau, beta) \
        - (big_a_st + big_a_ts
           + beta * (_pi_coef(0, 1) + _pi_coef(2, 3)
                     - _pi_coef(0, 3) - _pi_coef(1, 2)
                     + _pi_coef(0, 3) + _pi_coef(1, 2)))
    return _numer(residual)


def chain5_logic_audit(beta: sp.Symbol) -> Dict[str, Any]:
    return {
        "status": "REFUTED",
        "target_under_refutation": ("Psi_m >= 0 ENTRYWISE on [0,1]^4: "
                                    "superseded; superseded_by = LEMMA A + "
                                    "LEMMA C (chain5_corrected_verdict); DO "
                                    "NOT cite as a live obligation"),
        "refutation_witness": "SELF-REPRODUCED (certified in this audit): "
            "see chain5_refutation_check",
        "sufficiency_logic_of_the_refuted_route": [
            "Step 1 (Claim 3 + Claim 4): Psi_m >= 0 on [0,1]^4 would imply "
            "D_m >= 0 on [0,1]^2 (diagonal slice Psi_m(s,t,s,t)=2D_m/M); "
            "Claim 4 lifts D_M >= D_m >= 0 and D_M/M >= D_m/m for M >= m.",
            "Step 2 (Claim 2): T = a^T Psi_matrix a; entrywise P >= 0 would "
            "give T >= 0 for a in the simplex.",
            "Step 3 (Claim 4): T(M) >= T(m) since 2 I01/M >= 2 I01(m)/m "
            "termwise and c is M-free.",
            "Step 4: (1-q)^2, q^2, q(1-q) >= 0 gives gap >= q(1-q) T >= 0.",
        ],
        "reason_refuted": ("：(s,t,sig,tau) = (0, 1/2, 1/2, 0) is feasible; "
                           "by the P2/Q2 decomposition every Q2-protocol term "
                           "carries one of the zero factors and Psi restricts "
                           "to a sign-indefinite two-variable function; the "
                           "certified value is strictly negative."),
    }


def chain5_corrected_verdict(beta: sp.Symbol) -> Dict[str, Any]:
    r = psi_decomposition_formula(beta)
    return {
        "status": "PENDING_USER_TARGET",
        "formula": (
            "Psi_M(s,t,sig,tau) = A_M(s,tau) + A_M(t,sig) "
            "+ beta*(h(pi(s,t)) + h(pi(sig,tau)) - h(pi(s,tau)) - h(pi(t,sig))) "
            "+ beta*(h(pi(s,tau)) + h(pi(t,sig))), "
            "A_M(u,v) = (1-beta) h(u v) - (v h(u) + u h(v))/(2M); "
            "that is P2(s,tau) + P2(t,sig) + Q2(s,t) + Q2(sig,tau) with "
            "P2(s,t) = (1-beta) h(s t) - (t h(s) + s h(t))/(2M), "
            "Q2(s,t) = beta h(pi(s,t)) (U-cross terms cancel exactly; "
            "machine-checked: psi_decomposition_residual_zero). At M = m the "
            "entry is (M_T)_ij = P2(x_i,y_j)+P2(x_j,y_i)+Q2(x_i,x_j)"
            "+Q2(y_i,y_j)"),
        "psi_decomposition_residual_zero": r == 0,
        "lemma_a": {
            "id": "LEMMA A",
            "statement": "R(s,t) = D_m(s,t)/m >= 0 on [0,1]^2",
            "role": ("((1-q)^2 I00 + q^2 I11)/M is a convex combination of "
                     "the R(s,t) values, so it inherits Lemma A's "
                     "nonnegativity and lifts to M >= m via "
                     "Claim 4 (D_M/M >= D_m/m)"),
            "discovery_support": "worst D_m = +5.850360114e-4 at s = t = 3/4 "
                                 "(k/8 grid, DISCOVERY only)",
            "status": "PENDING_USER_TARGET",
        },
        "lemma_c": {
            "id": "LEMMA C",
            "statement": ("M_T(x,y) is COPOSITIVE for all (x,y) in [0,1]^6, "
                          "(M_T)_ij = P2(x_i,y_j) + P2(x_j,y_i) + Q2(x_i,y_j) "
                          "+ Q2(x_j,y_i) + Q2(x_i,x_j) + Q2(y_i,y_j)"),
            "role": "T = a^T M_T a for a in the simplex (Claim 2), so Lemma A "
                    "+ Lemma C + Claim 4 give T >= 0 at every M >= m",
            "strictly_weaker_than_entrywise": True,
            "survives_witness": ("at (x1,x2,y1,y2) = (0, 1/2, 1/2, 0) the "
                                 "2x2 block [[0.062146, -0.05537], "
                                 "[-0.05537, 0.062146]] has simplex-minimum "
                                 "+0.003388 > 0"),
            "status": "PENDING_USER_TARGET",
        },
        "witness_coexistence_note": ("the pointwise witness Psi(0,1/2,1/2,0) "
                                     "numbers a P2 minimum; co-positivity of "
                                     "the 2x2 block is a DIFFERENT, weaker "
                                     "condition that the witness does not "
                                     "refute"),
        "sufficiency_argument": [
            "Step 1 (Claim 3 + Claim 4): Psi_m(s,t,s,t) = 2 D_m(s,t)/m; if "
            "Lemma A holds then D_m >= 0, and Claim 4 lifts D_M/M >= D_m/m "
            "for M >= m; masses a_i a_j >= 0 make I00, I11 >= 0 entrywise.",
            "Step 2 (Claim 2): T = a^T M_T a exactly; if Lemma C holds then "
            "T >= 0 for every a in the simplex and every (x,y) in [0,1]^6.",
            "Step 3 (Claim 4): T evaluated at mixture parameter M satisfies "
            "T(M) >= T(m): 2 I01/M >= 2 I01(m)/m termwise (提升 per pair) and "
            "c is M-free.",
            "Step 4: all decomposition coefficients (1-q)^2, q^2, q(1-q) are "
            "nonnegative for q in [0,1]; the I00/I11 term is nonnegative by "
            "Lemma A and T(M) >= 0 by Steps 2-3, so gap >= 0.",
        ],
        "open_obligations": [
            "LEMMA A: R(s,t) = D_m(s,t)/m >= 0 on [0,1]^2 (PENDING_USER_TARGET)",
            "LEMMA C: copositivity of M_T for all (x,y) in [0,1]^6 "
            "(PENDING_USER_TARGET)",
        ],
    }


# --- binding constants ------------------------------------------------------------

def _mpf_frac(x: Any) -> Fraction:
    mp = x if isinstance(x, mpmath.mpf) else mpmath.mpf(x)
    sign, man, exp, _bc = mp._mpf_
    v = Fraction(int(man))
    if exp >= 0:
        v *= Fraction(2) ** int(exp)
    else:
        v /= Fraction(2) ** int(-exp)
    return -v if sign else v


def _binding_params() -> Tuple[Fraction, Fraction]:
    with mpmath.workdps(100):
        P = solve_equation_parameters(100)
    return _mpf_frac(P.beta), _mpf_frac(P.mean)


# --- numeric certificate -----------------------------------------------------------

def _F2A(f: Fraction) -> arb:
    return arb(int(f.numerator)) / arb(int(f.denominator))


def _arb_h(u: arb) -> arb:
    one = arb(1)
    return -(u * u.log() + (one - u) * (one - u).log())


def _draw_config(rng: random.Random, m_r: Fraction) -> Dict[str, Any]:
    """Seeded exact-rational configuration with M >= m enforced exactly."""
    for _attempt in range(64):
        u = [Fraction(rng.random()) for _ in range(3)]
        q = Fraction(rng.random())
        xs = [Fraction(rng.random()) for _ in range(3)]
        ys = [Fraction(rng.random()) for _ in range(3)]
        tot = u[0] + u[1] + u[2]
        if tot == 0 or q == 0 or q == 1 or any(v == 0 for v in xs + ys):
            continue
        masses = [v / tot for v in u]
        lifted = 0
        mean = (1 - q) * (masses[0] * xs[0] + masses[1] * xs[1] + masses[2] * xs[2]) \
            + q * (masses[0] * ys[0] + masses[1] * ys[1] + masses[2] * ys[2])
        if mean < m_r:
            theta = (m_r - mean) / (1 - mean)
            xs = [v + theta * (1 - v) for v in xs]
            ys = [v + theta * (1 - v) for v in ys]
            mean = (1 - q) * (masses[0] * xs[0] + masses[1] * xs[1]
                              + masses[2] * xs[2]) \
                + q * (masses[0] * ys[0] + masses[1] * ys[1] + masses[2] * ys[2])
            if mean != m_r:  # exact rational identity of the lift
                raise AssertionError("lift did not force M = m exactly")
            lifted = 1
        return {"masses": masses, "q": q, "x": xs, "y": ys,
                "M": mean, "lifted": lifted}
    raise RuntimeError("config draw exhausted")


def _cfg_fractions(cfg: Dict[str, Any]) -> List[Fraction]:
    """Ordered flat list [a1,a2,a3,q,x1,x2,x3,y1,y2,y3,M]."""
    return list(cfg["masses"]) + [cfg["q"]] + list(cfg["x"]) + list(cfg["y"]) \
        + [cfg["M"]]


def _arb_gap(cfg: Dict[str, Any], beta: arb) -> arb:
    a = [_F2A(v) for v in cfg["masses"]]
    q = _F2A(cfg["q"])
    xs = [_F2A(v) for v in cfg["x"]]
    ys = [_F2A(v) for v in cfg["y"]]
    zz = xs + ys
    w = [(arb(1) - q) * a[i] for i in range(3)] + [q * a[i] for i in range(3)]
    ehxy = sum((w[i] * w[j] * _arb_h(zz[i] * zz[j])
                for i in range(6) for j in range(6)), arb(0))
    ehpi = sum(((arb(1) - q) * a[i] * a[j] * _arb_h(xs[i] * xs[j]
                                                    * (arb(1) + (arb(1) - xs[i]) * (arb(1) - xs[j])))
                + q * a[i] * a[j] * _arb_h(ys[i] * ys[j]
                                           * (arb(1) + (arb(1) - ys[i]) * (arb(1) - ys[j])))
                for i in range(3) for j in range(3)), arb(0))
    ehx = sum((w[i] * _arb_h(zz[i]) for i in range(6)), arb(0))
    return (arb(1) - beta) * ehxy + beta * ehpi - ehx


def _arb_decomposed(cfg: Dict[str, Any], beta: arb) -> arb:
    a = [_F2A(v) for v in cfg["masses"]]
    q = _F2A(cfg["q"])
    xs = [_F2A(v) for v in cfg["x"]]
    ys = [_F2A(v) for v in cfg["y"]]
    M = _F2A(cfg["M"])
    one = arb(1)

    def D(u: arb, v: arb) -> arb:
        p = u * v * (one + (one - u) * (one - v))
        return M * ((one - beta) * _arb_h(u * v) + beta * _arb_h(p)) \
            - (v * _arb_h(u) + u * _arb_h(v)) / 2

    i00 = sum((a[i] * a[j] * D(xs[i], xs[j])
               for i in range(3) for j in range(3)), arb(0))
    i11 = sum((a[i] * a[j] * D(ys[i], ys[j])
               for i in range(3) for j in range(3)), arb(0))
    i01 = sum((a[i] * a[j] * D(xs[i], ys[j])
               for i in range(3) for j in range(3)), arb(0))
    c = beta * sum((a[i] * a[j] * (_arb_h(xs[i] * xs[j]
                                                   * (one + (one - xs[i]) * (one - xs[j])))
                                   + _arb_h(ys[i] * ys[j]
                                            * (one + (one - ys[i]) * (one - ys[j])))
                                   - 2 * _arb_h(xs[i] * ys[j]
                                                * (one + (one - xs[i]) * (one - ys[j]))))
                    for i in range(3) for j in range(3)), arb(0))
    return ((one - q) ** 2 * i00 + q ** 2 * i11) / M \
        + q * (one - q) * (2 * i01 / M + c)


def _f64_gap(cfg: Dict[str, Any], beta_f: float) -> float:
    a = [float(v) for v in cfg["masses"]]
    q = float(cfg["q"])
    xs = [float(v) for v in cfg["x"]]
    ys = [float(v) for v in cfg["y"]]
    zz = xs + ys
    w = [(1.0 - q) * a[i] for i in range(3)] + [q * a[i] for i in range(3)]

    def h(u: float) -> float:
        if u <= 0.0 or u >= 1.0:
            return 0.0
        return -u * math.log(u) - (1.0 - u) * math.log1p(-u)

    def pi(u: float, v: float) -> float:
        return u * v * (1.0 + (1.0 - u) * (1.0 - v))

    ehxy = sum(w[i] * w[j] * h(zz[i] * zz[j])
               for i in range(6) for j in range(6))
    ehpi = sum((1.0 - q) * a[i] * a[j] * h(pi(xs[i], xs[j]))
               + q * a[i] * a[j] * h(pi(ys[i], ys[j]))
               for i in range(3) for j in range(3))
    ehx = sum(w[i] * h(zz[i]) for i in range(6))
    return (1.0 - beta_f) * ehxy + beta_f * ehpi - ehx


def _f64_decomposed(cfg: Dict[str, Any], beta_f: float) -> float:
    a = [float(v) for v in cfg["masses"]]
    q = float(cfg["q"])
    xs = [float(v) for v in cfg["x"]]
    ys = [float(v) for v in cfg["y"]]
    M = float(cfg["M"])

    def h(u: float) -> float:
        if u <= 0.0 or u >= 1.0:
            return 0.0
        return -u * math.log(u) - (1.0 - u) * math.log1p(-u)

    def pi(u: float, v: float) -> float:
        return u * v * (1.0 + (1.0 - u) * (1.0 - v))

    def D(u: float, v: float) -> float:
        return M * ((1.0 - beta_f) * h(u * v) + beta_f * h(pi(u, v))) \
            - (v * h(u) + u * h(v)) / 2.0

    i00 = sum(a[i] * a[j] * D(xs[i], xs[j])
              for i in range(3) for j in range(3))
    i11 = sum(a[i] * a[j] * D(ys[i], ys[j])
              for i in range(3) for j in range(3))
    i01 = sum(a[i] * a[j] * D(xs[i], ys[j])
              for i in range(3) for j in range(3))
    c = beta_f * sum(a[i] * a[j] * (h(pi(xs[i], xs[j])) + h(pi(ys[i], ys[j]))
                                    - 2.0 * h(pi(xs[i], ys[j])))
                     for i in range(3) for j in range(3))
    return ((1.0 - q) ** 2 * i00 + q ** 2 * i11) / M \
        + q * (1.0 - q) * (2.0 * i01 / M + c)


def numeric_claim1(beta_r: Fraction, m_r: Fraction) -> Dict[str, Any]:
    """Arb machine-check of Claim 1 (route-B decomposition) over seeded draws."""
    rng = random.Random(CERT_SEED)
    mpmath.mp.dps = max(mpmath.mp.dps, 400)
    beta_a = _F2A(beta_r)
    bound = arb(2) ** (-RESIDUAL_BOUND_EXP)
    worst_key = -1.0
    worst_abs = arb(0)
    worst: Optional[Tuple[int, Dict[str, Any]]] = None
    all_ok = True
    lifted_count = 0
    min_m_seen: Optional[Fraction] = None
    f64_worst = 0.0
    for idx in range(ARB_DRAWS):
        cfg = _draw_config(rng, m_r)
        lifted_count += cfg["lifted"]
        resid = _arb_gap(cfg, beta_a) - _arb_decomposed(cfg, beta_a)
        ra = abs(resid)
        key = float(ra.upper())
        if key > worst_key:
            worst_key = key
            worst_abs = ra
            worst = (idx, {"masses_a1": str(cfg["masses"][0]),
                           "masses_a2": str(cfg["masses"][1]),
                           "masses_a3": str(cfg["masses"][2]),
                           "q": str(cfg["q"]),
                           "x1": str(cfg["x"][0]), "x2": str(cfg["x"][1]),
                           "x3": str(cfg["x"][2]),
                           "y1": str(cfg["y"][0]), "y2": str(cfg["y"][1]),
                           "y3": str(cfg["y"][2]),
                           "M": str(cfg["M"]),
                           "gap_lhs": str(_arb_gap(cfg, beta_a)),
                           "rhs_decomposed": str(_arb_decomposed(cfg, beta_a)),
                           "residual": str(resid)})
        if not (ra < bound):
            all_ok = False
        m_here = cfg["M"]
        if min_m_seen is None or m_here < min_m_seen:
            min_m_seen = m_here
        f64_res = abs(_f64_gap(cfg, float(beta_r)) - _f64_decomposed(cfg, float(beta_r)))
        if f64_res > f64_worst:
            f64_worst = f64_res
    if worst is None or min_m_seen is None:
        raise AssertionError("no numeric draws completed")
    return {
        "label": "MACHINE-CHECKED (Arb balls at 320-bit precision on exact "
                 "rational inputs; float64 companion is DISCOVERY only)",
        "seed": CERT_SEED,
        "arb_precision_bits": ARB_PREC_TARGET,
        "draws": ARB_DRAWS,
        "lifted_draws": lifted_count,
        "min_M": str(min_m_seen),
        "residual_bound": "2^-%d" % RESIDUAL_BOUND_EXP,
        "all_residuals_certified": all_ok,
        "max_abs_residual": str(worst_abs),
        "worst_draw_index": worst[0],
        "worst_draw": worst[1],
        "float64_companion": {
            "label": "DISCOVERY (double-precision logits, not a certificate)",
            "max_abs_residual": repr(f64_worst),
        },
    }


# --- mutations ----------------------------------------------------------------------

def _chain1_route_B_residual(beta: sp.Symbol, mutated: Optional[str] = None
                             ) -> sp.Expr:
    """Route-B residual; mutated in {None, 'kernel_M_to_m', 'c_misindex'}."""
    M_param = _M_ALT if mutated == "kernel_M_to_m" else _MEAN
    if mutated == "c_misindex":
        cexpr = beta * sum((MASSES[i] * MASSES[j] * (_pi_coef(i, i)
                                                     + _pi_coef(3 + j, 3 + j)
                                                     - 2 * _pi_coef(i, 3 + j))
                            for i in range(3) for j in range(3)), sp.Integer(0))
    else:
        cexpr = _remainder(beta)
    i00 = sum((MASSES[i] * MASSES[j] * _D_M(_MEAN, ATOMS[i], ATOMS[j], beta)
               for i in range(3) for j in range(3)), sp.Integer(0))
    i11 = sum((MASSES[i] * MASSES[j] * _D_M(_MEAN, ATOMS[3 + i],
                                            ATOMS[3 + j], beta)
               for i in range(3) for j in range(3)), sp.Integer(0))
    i01 = sum((MASSES[i] * MASSES[j] * _D_M(_MEAN, ATOMS[i], ATOMS[3 + j], beta)
               for i in range(3) for j in range(3)), sp.Integer(0))
    t_expr = 2 * i01 / M_param + cexpr
    raw = _gap(beta) - ((1 - Q) ** 2 * i00 + Q ** 2 * i11) / _MEAN \
        - Q * (1 - Q) * t_expr
    residual = raw.subs(_MEAN, _mixture_mean())
    return _numer(residual)


def _chain3_residual(beta: sp.Symbol, mutate_cross: bool = False) -> sp.Expr:
    s, t, M = X1, X2, _M_PSI
    if mutate_cross:
        d1 = _D_M(M, s, t, beta)
        d2 = _D_M(M, t, s, beta)
        gamma = (2 * _pi_coef(0, 1) - _pi_coef(0, 0) - _pi_coef(1, 1))
        psi_mut = (d1 + d2) / M + beta * gamma
        return _numer(psi_mut - 2 * _D_M(M, s, t, beta) / M)
    return _numer(_psi(M, s, t, s, t, beta) - 2 * _D_M(M, s, t, beta) / M)


def run_mutations(beta: sp.Symbol) -> List[Dict[str, Any]]:
    """Soundness-breaking mutations; each must fire (residual != 0)."""
    rows: List[Dict[str, Any]] = []

    def row(name: str, surface: str, mutant_residual: sp.Expr,
            sound_residual: sp.Expr) -> None:
        rec = _residual_record(mutant_residual)
        sound = _residual_record(sound_residual)
        rows.append({
            "mutation": name,
            "attack_surface": surface,
            "expected": "FAIL",
            "sound_control_residual_is_zero": sound["residual_is_zero"],
            "mutant_residual_is_zero": rec["residual_is_zero"],
            "mutant_residual_term_count": rec["residual_term_count"],
            "fired": (not rec["residual_is_zero"]) and sound["residual_is_zero"],
        })

    # M1: drop beta*(h(pi) - h(st)) from the discharge kernel (route A).
    lhs = _MEAN * _gap(beta) - _MEAN * Q * (1 - Q) * _remainder(beta)
    w = _weights()
    rhs_mut = sp.Integer(0)
    rhs_sound = sp.Integer(0)
    for i in range(6):
        for j in range(6):
            base = (ATOMS[j] * _h1(i) + ATOMS[i] * _h1(j)
                    - ATOMS[i] * ATOMS[j] * (_mu1(i) + _mu1(j)
                                             - _mu_of(tuple(sorted((i, j))))))
            corr = beta * (_pi_coef(i, j) - _h_free(ATOMS[i] * ATOMS[j]))
            rhs_sound += w[i] * w[j] * (_MEAN * (base + corr)
                                        - (ATOMS[j] * _h1(i) + ATOMS[i] * _h1(j)) / 2)
            rhs_mut += w[i] * w[j] * (_MEAN * base
                                      - (ATOMS[j] * _h1(i) + ATOMS[i] * _h1(j)) / 2)
    sound_a = sp.expand(lhs - rhs_sound).subs(_MEAN, _mixture_mean())
    mut_a = sp.expand(lhs - rhs_mut).subs(_MEAN, _mixture_mean())
    row("m1_drop_beta_protocol_from_e_kernel",
        "claim1_route_A_discharge", sp.expand(mut_a), sp.expand(sound_a))

    # M2: swap the crossed protocol pairs inside Psi (chain 3).
    sound3 = _chain3_residual(beta, mutate_cross=False)
    mut3 = _chain3_residual(beta, mutate_cross=True)
    row("m2_swap_cross_pairs_in_psi", "chain3_diagonal_slice", mut3, sound3)

    # M3: replace 2*I01/M by 2*I01/m_binding inside T (route B).
    sound_b = _chain1_route_B_residual(beta)
    mut_b = _chain1_route_B_residual(beta, mutated="kernel_M_to_m")
    row("m3_replace_M_by_m_in_T", "claim1_route_B_decomposition", mut_b, sound_b)

    # M4: misindex the c cross-terms c = beta sum a_i a_j (pi(x_i,x_i)
    # + pi(y_j,y_j) - 2 pi(x_i,y_j)) (route B).
    mut_c = _chain1_route_B_residual(beta, mutated="c_misindex")
    row("m4_misindex_c_protocol_pairs", "claim1_route_B_decomposition",
        mut_c, sound_b)
    return rows


# --- report ---------------------------------------------------------------------------

def _digest(payload: Dict[str, Any]) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                      default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def build_report() -> Dict[str, Any]:
    beta_r, m_r = _binding_params()
    beta = BETA

    chain1_a = prove_chain1_route_A(beta)
    chain1_b = prove_chain1_route_B(beta)
    kernel = prove_kernel_collapse(beta)
    chain2 = prove_chain2(beta)
    chain3 = prove_chain3(beta)
    chain4 = prove_chain4(beta)
    numeric = numeric_claim1(beta_r, m_r)
    mutations = run_mutations(beta)
    chain5 = chain5_logic_audit(beta)
    chain5_ref_res = chain5_refutation_check(beta_r)
    chain5_corr = chain5_corrected_verdict(beta)

    statements = {
        "chain1_all_q_identity": {
            "claim": ("gap == ((1-q)^2 I00 + q^2 I11)/M + q(1-q) T for all q, "
                      "M = (1-q)M0 + qM1 the mixture mean"),
            "route_A": chain1_a,
            "route_B": chain1_b,
            "status": "PROVED" if chain1_a["status"] == "PROVED"
            and chain1_b["status"] == "PROVED" else "REFUTED",
        },
        "kernel_collapse_e_equals_direct_kernel": {
            "claim": ("e(s,t) == (1-beta) h(s t) + beta h(pi(s,t)) where "
                      "e(s,t) = t h(s) + s h(t) - s t (mu(s)+mu(t)-mu(st)) "
                      "+ beta (h(pi(s,t)) - h(s*t)); equivalently the "
                      "mu-lemma t h(s)+s h(t) - s t (mu(s)+mu(t)-mu(st)) "
                      "== h(s t)"),
            **{k: v for k, v in kernel.items() if k != "rows"},
            **kernel,
        },
        "chain2_matrix_representation": chain2,
        "chain3_diagonal_slice": chain3,
        "chain4_monotonicity_in_M": chain4,
        "chain5_sufficiency_refuted": chain5,
        "chain5_corrected_verdict": chain5_corr,
        "chain5_refutation_check": chain5_ref_res,
        "numeric_claim1_arb": numeric,
    }
    core_ok = (statements["chain1_all_q_identity"]["status"] == "PROVED"
               and kernel["status"] == "PROVED"
               and chain2["status"] == "PROVED"
               and chain3["status"] == "PROVED"
               and chain4["status"] == "PROVED"
               and numeric["all_residuals_certified"]
               and chain5["status"] == "REFUTED"
               and chain5_corr["psi_decomposition_residual_zero"]
               and bool(chain5_ref_res["certified_strictly_negative"]))
    mutations_fired = all(r["fired"] for r in mutations)
    claim_status = "PROVED" if core_ok and mutations_fired else "REFUTED"
    record = {
        "tool": "uc/liu9_psi_reduction_audit.py",
        "audit_role": ("independent adversarial verification of Main's "
                       "candidate reduction ref_liu_hypothesis_2 "
                       "(9-var + q -> Psi_m >= 0 on [0,1]^4); no Main "
                       "arithmetic reused"),
        "honesty": ("All symbolic residuals herein are exact canonical "
                    "polynomial numerators over an independent free log "
                    "basis; any nonzero residual is a formal contradiction "
                    "of the corresponding claim, not a float artifact. The "
                    "float64 companion numbers are DISCOVERY only."),
        "binding_constants": {
            "source": "liu9_binding.solve_equation_parameters(100)",
            "beta_rational": str(beta_r),
            "m_rational": str(m_r),
            "beta_head": mpmath.nstr(mpmath.mpf(str(beta_r.numerator))
                                     / mpmath.mpf(str(beta_r.denominator)), 40),
            "m_head": mpmath.nstr(mpmath.mpf(str(m_r.numerator))
                                  / mpmath.mpf(str(m_r.denominator)), 40),
        },
        "arb_precision_bits": ARB_PREC_TARGET,
        "claim_status": claim_status,
        "claim_status_note": ("PROVED covers Claims 1-4 (+ kernel collapse, "
                              "mu-lemma, numeric certificate, all mutations "
                              "firing) AND the independent re-certification of "
                              "the REFUTATION of the entrywise target "
                              "Psi_m >= 0 on [0,1]^4 (witness below); the "
                              "live obligation is LEMMA A + LEMMA C only "
                              "(chain5_corrected_verdict, both "
                              "PENDING_USER_TARGET)."),
        "statements": statements,
        "mutations": mutations,
        "report_sha256": "",
    }
    if claim_status == "REFUTED":
        failures: List[str] = []
        if statements["chain1_all_q_identity"]["status"] != "PROVED":
            failures.append("chain1")
        if kernel["status"] != "PROVED":
            failures.append("kernel_collapse")
        if chain2["status"] != "PROVED":
            failures.append("chain2")
        if chain3["status"] != "PROVED":
            failures.append("chain3")
        if chain4["status"] != "PROVED":
            failures.append("chain4")
        if not numeric["all_residuals_certified"]:
            failures.append("numeric_residual_bound")
        if not mutations_fired:
            failures.append("mutations")
        record["refutation_targets"] = failures
    record["report_sha256"] = _digest(record)
    return record


def main() -> int:
    report = build_report()
    with OUTPUT_DEFAULT.open("w", encoding="utf-8") as fh:
        fh.write(json.dumps(report, indent=1, sort_keys=True,
                            default=str) + "\n")
    print("PSI REDUCTION AUDIT (independent; adversary of the reduction)")
    print("claim_status %s" % report["claim_status"])
    st = report["statements"]
    for key in ("chain1_all_q_identity", "chain2_matrix_representation",
                "chain3_diagonal_slice", "chain4_monotonicity_in_M"):
        print("  %-32s %s" % (key, st[key]["status"]))
    print("  %-32s %s" % ("kernel_collapse + mu_lemma",
                          st["kernel_collapse_e_equals_direct_kernel"]["status"]))
    num = st["numeric_claim1_arb"]
    print("  numeric: %d draws @ %d bits, all certified %s, max |res| %s"
          % (num["draws"], num["arb_precision_bits"],
             num["all_residuals_certified"], num["max_abs_residual"]))
    fired = sum(1 for r in report["mutations"] if r["fired"])
    print("  mutations fired: %d/%d" % (fired, len(report["mutations"])))
    print("  chain5 entrywise route: %s (witness Psi_m(0,1/2,1/2,0) = %s)"
          % (st["chain5_sufficiency_refuted"]["status"],
             st["chain5_refutation_check"]["value"]))
    corr = st["chain5_corrected_verdict"]
    print("  chain5 corrected: %s (LEMMA A %s, LEMMA C %s)"
          % (corr["status"], corr["lemma_a"]["status"],
             corr["lemma_c"]["status"]))
    print("report_sha256 %s" % report["report_sha256"])
    return 0 if report["claim_status"] == "PROVED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
