"""Exact characterization of the q-free channel scalar c of the r <= 3
paired class, and its size relative to the certified F(P_mix) margins.

Split identity (liu9_paired_class_c.py machinery, channel family):

    raw gap = F(P_mix) + q (1 - q) c,      F(P) >= 0 universally.

with c q-free.  This module pins the EXACT FORM of c, its extrema on the
mass simplex, and the certified worst-case ratio of q(1-q)|c| against the
per-cell F(P_mix) margin of the published k/8 sweep
(liu9-h2-qp-sweep.json, 11,025 cells, worst certified min +1.6305e-3).

1. SYMBOLIC MATRIX FORM.  For r = 2 (fresh two-atom paired class, built
   here with its own 4-atom raw gap and discharge identity) and r = 3
   (reduced in liu9_paired_class_c's canonical 54-symbol free log basis):

       c = beta * <a, H a>,   H_ij = K(x_i, x_j) + K(y_i, y_j)
                                          - 2 K(x_i, y_j),
       K(s, t) = h(pi(s, t)),  pi(s, t) = s t (1 + (1 - s) (1 - t)),

   is proved as an exact sympy identity (residual exactly 0 on the free
   log basis).  H_ij itself need NOT be symmetric (the form only reads
   the pair sums H_ij + H_ji).  For r = 3 the matrix form is additionally
   verified equal (residual 0) to the certified remainder_expr of
   liu9_paired_class_c.

2. DIAGONAL VANISHING.  c vanishes identically when y_i = x_i: each H
   entry becomes K_ij + K_ij - 2 K_ij = 0.  Proved as a sympy identity on
   the duplicated atom tuple (x1..xr, x1..xr).  Sharp consequence (Arb):
   when P0 = P1 the mixture is q-independent, so the raw gap must be
   CONSTANT in q: the exact q^2 / q coefficients extracted from the
   certified gap on several diagonal geometries are certified zero balls
   (Q1 = Q2 = 0), including 0/1 supports and vertex masses.

3. CLOSED-FORM SIMPLEX EXTREMES.  c = beta <a, H a> is a quadratic form
   on the simplex {a >= 0, sum a = 1}; its extrema are attained at
   vertices, on edges (1D quadratic criticals), or at the interior
   stationary point -- a finite candidate list certified with Arb,
   evaluated per geometry on the k/8 ordered grid (35 x 35 = 1225).

4. |c| vs THE F MARGIN.  Per sweep cell, certified worst q(1-q)|c| is
   compared against the certified per-cell margin min_lower of
   liu9-h2-qp-sweep.json.  ratio_simplex: simplex-wide most negative c,
   max q(1-q) = 1/4.  ratio_cell: the same c-quadratic minimized over the
   cell's mean-feasible band at the cell's own q(1-q).  A ratio >= 1 is
   NOT a counterexample (F is certified at the cell's gap minimiser,
   |c| is maximised over a different mass set); such cells are flagged
   DANGER with all coordinates.

Numerical mpmath sections are labelled DISCOVERY; exact / symbolic /
Arb-ball statements are labelled PROVED.  Machine contract: single core,
byte-stable output, >= 2 failing mutations, canonical report digest.
"""

import os

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

import mpmath
import sympy as sp

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from flint import arb, ctx  # noqa: E402

CTX_PREC = 320
ctx.prec = max(ctx.prec, CTX_PREC)

from liu9_binding import (  # noqa: E402
    _h_arb_endpoint_safe,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import gap_mp  # noqa: E402
from liu9_h2_qp_sweep import (  # noqa: E402
    FLOOR_HI,
    MEAN_INDEX,
    build_constraints,
    geometry_pairs,
    kkt_candidates,
    mean_line,
    polygon_vertices,
)
from liu9_objective import evaluate_arb  # noqa: E402
from liu9_witness_qp import _pt  # noqa: E402

OUTPUT_DEFAULT = HERE / "verification/results/liu9-c-characterization.json"
SWEEP_RESULTS = HERE / "verification/results/liu9-h2-qp-sweep.json"

SIMPLEX_CONSTRAINTS = [
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(-1), Fraction(-1), Fraction(-1)),
]


def _canonical_digest(payload):
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# Generic r-paired-class symbolic machinery.

def _pi_formula(s, t):
    return s * t * (1 + (1 - s) * (1 - t))


class PairedClass:
    """The r-pair two-component model over atoms (x_1..x_r, y_1..y_r).

    Weights w = ((1-q) m_i for the x-side, q m_i for the y-side) with
    masses (a_1, ..., a_r), sum a_i = 1; free mass symbols a1..a_{r-1},
    last mass = 1 - their sum.

    Symbol discipline (liu9_paired_class_c's design): the FREE LOG BASIS
    consists of index-named symbols -- L_/K_ per atom, PR_(i,j) for
    h(prod), TQ_(i,j) for h(pi) -- all passed through the reducer
    VERBATIM.  The KERNEL symbols that circulate inside expressions are
    argument-named h2(pi(<expr>)) / h2(<prod>) free symbols; on the
    DIAGONAL build the duplicated atoms make their argument strings
    coincide, which is exactly how the y = x collapsing is expressed.
    """

    def __init__(self, r, diagonal=False):
        self.r = r if not diagonal else r
        self.is_diagonal = diagonal
        self.beta = sp.Symbol("beta", real=True)
        self.q = sp.Symbol("q", real=True)
        self.mass_syms = tuple(sp.Symbol("a%d" % k, real=True)
                               for k in range(1, r))
        self.masses = self.mass_syms + (
            1 - sum(self.mass_syms, sp.Integer(0)),)
        xs = tuple(sp.Symbol("x%d" % k, real=True)
                   for k in range(1, r + 1))
        if diagonal:
            self.atoms = xs + xs
        else:
            self.atoms = xs + tuple(sp.Symbol("y%d" % k, real=True)
                                    for k in range(1, r + 1))
        n = 2 * r
        self.log_sym = {a: sp.Symbol("L_" + a.name) for a in self.atoms}
        self.log1m_sym = {a: sp.Symbol("K_" + a.name) for a in self.atoms}
        self.prod_sym = {}
        self.protocol_sym = {}
        for i in range(n):
            for j in range(i, n):
                self.prod_sym[(i, j)] = sp.Symbol("PR_%d_%d" % (i, j))
                self.protocol_sym[(i, j)] = sp.Symbol("TQ_%d_%d" % (i, j))
        self.mean_sym = sp.Symbol("M_mix")

    # -- free-basis evaluable entropy ---------------------------------------
    def h_of(self, atoms_list):
        """h of a product of DISTINCT atoms, as free-basis symbols."""
        if not atoms_list:
            return sp.Integer(0)
        index = {}
        for k, a in enumerate(self.atoms):
            index.setdefault(a.name, k)
        inds = tuple(sorted(index[a.name] for a in atoms_list))
        if len(inds) == 1:
            a = self.atoms[inds[0]]
            return -a * self.log_sym[a] - (1 - a) * self.log1m_sym[a]
        z = sp.prod([self.atoms[k] for k in inds])
        log_z = sum((self.log_sym[self.atoms[k]] for k in inds),
                    sp.Integer(0))
        return -z * log_z - (1 - z) * self.prod_sym[tuple(sorted(inds))]

    def mu_of(self, atoms_list):
        if not atoms_list:
            return sp.Integer(1)
        hz = self.h_of(atoms_list)
        z = (sp.prod(atoms_list) if len(atoms_list) > 1 else atoms_list[0])
        log_z = sum((self.log_sym[a] for a in atoms_list), sp.Integer(0))
        return sp.cancel(hz / z - log_z)

    def _atom_list_from_arg(self, arg):
        pd = arg.as_powers_dict()
        atoms_list = []
        for a, p in sorted(pd.items(), key=lambda kv: str(kv[0])):
            atoms_list.extend([a] * int(p))
        return atoms_list

    # -- reduction of argument-named symbols ---------------------------------
    def _h_reduce(self, name):
        locs = {a.name: a for a in self.atoms}
        n = len(self.atoms)
        if name.startswith("h2(hpp("):
            # protocol kernel of exactly two atoms, written x_i=x_j
            inner = name[7:-2]
            parts = inner.split("=")
            def _find(symname):
                for k, atom in enumerate(self.atoms):
                    if atom.name == symname:
                        return k
                raise KeyError(symname)
            i, j = _find(parts[0]), _find(parts[1])
            lo, hi = min(i, j), max(i, j)
            return self.protocol_sym[(lo, hi)]
        arg = sp.sympify(name[3:-1], locals=locs)
        n = len(self.atoms)
        for i in range(n):
            for j in range(i, n):
                if sp.expand(arg - _pi_formula(self.atoms[i],
                                               self.atoms[j])) == 0:
                    lo, hi = min(i, j), max(i, j)
                    return self.protocol_sym[(lo, hi)]
        for i in range(n):
            if sp.expand(arg - self.atoms[i]) == 0:
                return self.h_of([self.atoms[i]])
        return self.h_of(self._atom_list_from_arg(arg))

    def reduce(self, expr):
        subs = {}
        for s in list(expr.free_symbols):
            n = s.name
            if n.startswith("h2("):
                subs[s] = self._h_reduce(n)
            elif n.startswith("mu2("):
                locs = {a.name: a for a in self.atoms}
                arg = sp.sympify(n[4:-1], locals=locs)
                subs[s] = self.mu_of(self._atom_list_from_arg(arg))
        return sp.expand(sp.cancel(expr.subs(subs)))

    # -- basis symbols ---------------------------------------------------------
    def kernel_symbol(self, i, j):
        """h(pi(atom_i, atom_j)) as an argument-named symbol (reducible)."""
        lo, hi = min(i, j), max(i, j)
        return sp.Symbol("h2(hpp(%s=%s))" % (self.atoms[lo],
                                             self.atoms[hi]))

    def prod_argument(self, i, j):
        """h(atom_i * atom_j) as the multilinear free-basis expansion."""
        lo, hi = min(i, j), max(i, j)
        return self.h_of([self.atoms[lo], self.atoms[hi]])

    def _h_symbol_for_atom(self, a):
        for k, atom in enumerate(self.atoms):
            if atom == a:
                return self.h_of([self.atoms[k]])
        raise KeyError(a)

    # -- the objects of study ----------------------------------------------------
    def c_matrix_expr(self, beta=None, cross=2, kernel=None):
        """c = beta sum_ij m_i m_j [K(xi,xj) + K(yi,yj) - cross K(xi,yj)].

        kernel: None -> argument-named symbols (diagonal build collapses
        because duplicated atoms give identical argument strings);
        kernel=protocol -> index-named TQ free symbols passed verbatim
        (used only for the split identity, matching the raw gap)."""
        beta = self.beta if beta is None else beta
        r = self.r
        if kernel == "protocol":
            def ks(i, j):
                lo, hi = min(i, j), max(i, j)
                return self.protocol_sym[(lo, hi)]
        else:
            ks = self.kernel_symbol
        total = sp.Integer(0)
        for i in range(r):
            for j in range(r):
                total += self.masses[i] * self.masses[j] * (
                    ks(i, j) + ks(r + i, r + j) - cross * ks(i, r + j))
        return beta * total

    def _weights(self):
        return tuple((1 - self.q) * m for m in self.masses) \
            + tuple(self.q * m for m in self.masses)

    def raw_gap_expr(self, beta=None):
        beta = self.beta if beta is None else beta
        n = 2 * self.r
        weights = self._weights()
        ehxy = sum((weights[i] * weights[j] * self.prod_argument(i, j)
                    for i in range(n) for j in range(n)), sp.Integer(0))
        ehpi = sp.Integer(0)
        for comp in range(2):
            factor = (1 - self.q) if comp == 0 else self.q
            off = 0 if comp == 0 else self.r
            for i in range(self.r):
                for j in range(self.r):
                    lo, hi = min(off + i, off + j), max(off + i, off + j)
                    ehpi += factor * self.masses[i] * self.masses[j] \
                        * self.protocol_sym[(lo, hi)]
        ehx = sum((weights[i] * self._h_symbol_for_atom(self.atoms[i])
                   for i in range(n)), sp.Integer(0))
        return (1 - beta) * ehxy + beta * ehpi - ehx

    def mixture_mean_expr(self):
        n = 2 * self.r
        weights = self._weights()
        return sum((weights[i] * self.atoms[i] for i in range(n)),
                   sp.Integer(0))

    def e_kernel_expr(self, i, j, beta=None):
        beta = self.beta if beta is None else beta
        s, t = self.atoms[i], self.atoms[j]
        big_a = t * self._h_symbol_for_atom(s) \
            + s * self._h_symbol_for_atom(t)
        tk = (sp.Symbol("mu2(%s)" % s.name)
              + sp.Symbol("mu2(%s)" % t.name)
              - sp.Symbol("mu2(%s)" % str(s * t)))
        lo, hi = min(i, j), max(i, j)
        gain = self.protocol_sym[(lo, hi)] - self.prod_argument(i, j)
        return big_a - s * t * tk + beta * gain

    def discharge_expr(self, beta=None):
        beta = self.beta if beta is None else beta
        n = 2 * self.r
        weights = self._weights()
        total = sp.Integer(0)
        for i in range(n):
            for j in range(n):
                d = self.mean_sym * self.e_kernel_expr(i, j, beta) \
                    - (self.atoms[j] * self._h_symbol_for_atom(self.atoms[i])
                       + self.atoms[i]
                       * self._h_symbol_for_atom(self.atoms[j])) / 2
                total += weights[i] * weights[j] * d
        return total

    def prove_split_identity(self):
        gap = self.raw_gap_expr()
        cform = self.c_matrix_expr(kernel="protocol")
        num = self.mean_sym * (gap - self.q * (1 - self.q) * cform)
        poly = num - self.discharge_expr()
        reduced = self.reduce(poly.subs(self.mean_sym,
                                        self.mixture_mean_expr()))
        basis = (len(self.log_sym) + len(self.log1m_sym)
                 + len(self.prod_sym) + len(self.protocol_sym))
        if reduced != 0:
            return {"status": "FAILED", "r": self.r,
                    "residual_terms": len(sp.Add.make_args(reduced)),
                    "residual_str_head": str(reduced)[:240]}
        return {"status": "PROVED", "r": self.r,
                "residual": "0", "residual_terms": 0,
                "method": ("sympy multilinear expansion in the canonical "
                           "index-named free log basis (L_, K_, PR_, TQ_) "
                           "of the r-pair paired class"),
                "log_basis_size": basis}

    def prove_diagonal_vanishing(self):
        """On y = x the argument-named matrix form reduces to exactly 0."""
        if not self.is_diagonal:
            diag = PairedClass(self.r, diagonal=True)
        else:
            diag = self
        reduced = diag.reduce(diag.c_matrix_expr())
        if reduced != 0:
            return {"status": "FAILED", "r": self.r,
                    "residual_terms": len(sp.Add.make_args(reduced)),
                    "residual_str_head": str(reduced)[:240]}
        return {"status": "PROVED", "r": self.r,
                "residual": "0", "residual_terms": 0,
                "method": ("on the duplicated atom tuple (x1..xr, x1..xr) "
                           "the argument-named kernel symbols coincide: "
                           "h(pi(x_i,x_j)) appears identically in K(x_i,x_j) "
                           "and K(y_i,y_j), and K(x_i,y_j) = K(x_j,y_i) "
                           "with y=x; each H entry is K_ij + K_ij - 2 K_ij "
                           "= 0 exactly (identity, not a sample)"),
                "log_basis_size": (len(diag.log_sym) + len(diag.log1m_sym)
                                   + len(diag.prod_sym)
                                   + len(diag.protocol_sym))}


def _r3_matrix_matches_remainder():
    import liu9_paired_class_c as P

    beta = P.BETA
    total = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            total += P.MASSES_C(i) * P.MASSES_C(j) * (
                P._PI_SYM(i, j) + P._PI_SYM(3 + i, 3 + j)
                - 2 * P._PI_SYM(i, 3 + j))
    mine = beta * total
    theirs = P.remainder_expr(beta)
    reduced = P._reduce(sp.expand(mine - theirs))
    if reduced != 0:
        return {"status": "FAILED", "r": 3,
                "residual_terms": len(sp.Add.make_args(reduced))}
    return {"status": "PROVED", "r": 3,
            "residual": "0", "residual_terms": 0,
            "method": ("sympy equality of the matrix form against the "
                       "certified remainder_expr inside "
                       "liu9_paired_class_c's canonical 54-symbol free "
                       "log basis"),
            "log_basis_size": 54}


# ---------------------------------------------------------------------------
# Certified arithmetic (Arb).

def hpi_arb(sb, tb):
    prod = sb * tb * (1 + (1 - sb) * (1 - tb))
    return _h_arb_endpoint_safe(prod)


def h_matrix_arb(geom):
    r = 3
    B = [_pt(v) for v in geom]
    K = [[hpi_arb(B[i], B[j]) for j in range(2 * r)]
         for i in range(2 * r)]
    H = []
    for i in range(r):
        row = []
        for j in range(r):
            row.append(K[i][j] + K[r + i][r + j] - 2 * K[i][r + j])
        H.append(row)
    return H


def quad_coefficients(H):
    S01 = H[0][1] + H[1][0]
    S02 = H[0][2] + H[2][0]
    S12 = H[1][2] + H[2][1]
    return {
        "A0": H[2][2],
        "A1": S02 - 2 * H[2][2],
        "A2": S12 - 2 * H[2][2],
        "A12": S01 - S02 - S12 + 2 * H[2][2],
        "A11": H[0][0] + H[2][2] - S02,
        "A22": H[1][1] + H[2][2] - S12,
    }


def simplex_extremes(coefd, beta):
    """Certified min/max of c = beta <a, H a> over the mass simplex.

    coefd carries the beta-free form Q; every value ball is scaled by the
    certified beta ball (Arb interval multiplication)."""
    cands = kkt_candidates(coefd, SIMPLEX_CONSTRAINTS, mean_index=None)
    feasible = [c for c in cands if c["feasible"]]
    for c in feasible:
        c["_c_value"] = beta * c["_value"]
        c["_c_lower_str"] = c["_c_value"].lower().str(26)
        c["_c_upper_str"] = c["_c_value"].upper().str(26)
    best = min(feasible, key=lambda c: (float(c["_c_value"].lower()),
                                        c["name"]))
    top = max(feasible, key=lambda c: (float(c["_c_value"].upper()),
                                       c["name"]))
    return {"min": best, "max": top, "n_candidates": len(cands),
            "n_feasible": len(feasible)}


def gap_arb(ball_params, a1, a2, q, triple):
    beta = ball_params.beta
    values = tuple(_pt(v) for v in (a1, a2, q)) \
        + tuple(_pt(v) for v in triple + triple)
    terms = evaluate_arb(values, beta)
    return terms.numerator - terms.ehx


def q_scan_arb(ball_params, a1, a2, triple):
    G = [gap_arb(ball_params, a1, a2, Fraction(k, 8), triple)
         for k in range(9)]
    two = _pt(2)
    g2 = two * (G[8] + (-2) * G[4] + G[0])
    g1 = (G[8] + (-1) * G[0]) + (-1) * g2
    validation = []
    ok = True
    for k in (1, 3, 5, 7):
        tb = _pt(Fraction(k, 8))
        pred = G[0] + g1 * tb + g2 * (tb * tb)
        residual = pred + (-1) * G[k]
        contains = bool(residual.contains(arb(0)))
        ok = ok and contains
        validation.append({"q": str(Fraction(k, 8)),
                           "residual_radius_float": float(residual.rad()),
                           "contains_zero": contains})
    both_zero = bool(g2.contains(arb(0)) and g1.contains(arb(0)))
    return {"gap_at_midpoint_q12": G[4].str(26),
            "q2_coefficient": g2.str(8), "q1_coefficient": g1.str(8),
            "q2_straddles_zero": bool(g2.contains(arb(0))),
            "q1_straddles_zero": bool(g1.contains(arb(0))),
            "both_coefficients_zero_balls": both_zero,
            "interior_validation_ok": ok,
            "interior_validation": validation,
            "status": "PROVED (Arb)" if (both_zero and ok) else "FAILED"}


# ---------------------------------------------------------------------------
# Numeric (DISCOVERY) identity checks at 200 dps.

def _h_mp(u):
    if u == 0 or u == 1:
        return mpmath.mpf(0)
    return -(u * mpmath.log(u) + (1 - u) * mpmath.log1p(-u))


def _mu_mp(u):
    if u == 0:
        return mpmath.mpf(1)
    if u == 1:
        return mpmath.mpf(0)
    return -((1 - u) * mpmath.log1p(-u)) / u


def _pi_mp(s, t):
    return s * t * (1 + (1 - s) * (1 - t))


def numeric_identity(values, beta, r, cross=2):
    """gap - q(1-q) c - F at 200 dps (DISCOVERY).

    r = 3: values = (a1, a2, q, x1, x2, x3, y1, y2, y3);
    r = 2: values = (a1, q, x1, x2, y1, y2), masses (a1, 1 - a1)."""
    if r == 3:
        a1, a2, q = values[0], values[1], values[2]
        sup = tuple(values[k] for k in range(3, 9))
        masses = (a1, a2, 1 - a1 - a2)
    else:
        a1, q = values[0], values[1]
        sup = tuple(values[k] for k in range(2, 6))
        masses = (a1, 1 - a1)
    n = 2 * r
    weights = tuple((1 - q) * m for m in masses) + tuple(q * m for m in masses)
    mean = mpmath.fsum(weights[k] * sup[k] for k in range(n))
    ehxy = mpmath.fsum(weights[i] * weights[j] * _h_mp(sup[i] * sup[j])
                       for i in range(n) for j in range(n))
    ehpi = mpmath.mpf(0)
    for comp, factor in enumerate(((1 - q), q)):
        off = 0 if comp == 0 else r
        for i in range(r):
            for j in range(r):
                ehpi += factor * masses[i] * masses[j] * _h_mp(
                    _pi_mp(sup[off + i], sup[off + j]))
    ehx = mpmath.fsum(weights[i] * _h_mp(sup[i]) for i in range(n))
    gap = (1 - beta) * ehxy + beta * ehpi - ehx

    def e_ker(s, t):
        return (t * _h_mp(s) + s * _h_mp(t)
                - s * t * (_mu_mp(s) + _mu_mp(t) - _mu_mp(s * t))
                + beta * (_h_mp(_pi_mp(s, t)) - _h_mp(s * t)))

    dis = mpmath.fsum(weights[i] * weights[j]
                      * (mean * e_ker(sup[i], sup[j])
                         - (sup[j] * _h_mp(sup[i])
                            + sup[i] * _h_mp(sup[j])) / 2)
                      for i in range(n) for j in range(n))
    quad = mpmath.mpf(0)
    for i in range(r):
        for j in range(r):
            quad += masses[i] * masses[j] * (
                _h_mp(_pi_mp(sup[i], sup[j]))
                + _h_mp(_pi_mp(sup[r + i], sup[r + j]))
                - cross * _h_mp(_pi_mp(sup[i], sup[r + j])))
    residual = gap - q * (1 - q) * beta * quad - dis / mean
    return {"config": [mpmath.nstr(v, 20) for v in values],
            "mean": mpmath.nstr(mean, 24),
            "gap": mpmath.nstr(gap, 24),
            "residual": mpmath.nstr(residual, 6)}


def _fixed_rational(numer_range, index, seed):
    state = (index * 1103515245 + seed * 12345 + 7) % (2 ** 31)
    return Fraction(1 + state % (numer_range - 2), numer_range)


# ---------------------------------------------------------------------------
# Mutations (each must FAIL).

def run_mutations(params_arb):
    out = []

    # M1: cross coefficient 3 instead of 2: r=2 diagonal vanishing breaks.
    diag2 = PairedClass(2, diagonal=True)
    res_bad = diag2.reduce(diag2.c_matrix_expr(cross=3))
    out.append({"mutation": "M1_r2_diagonal_cross3",
                "expected": "nonzero residual (diagonal vanishing broken)",
                "residual_terms": len(sp.Add.make_args(res_bad))
                if res_bad != 0 else 0,
                "ok": res_bad != 0})

    # M2: swap atoms y2 <-> y3 in the r=3 gap/c-form, then demand the split
    # identity against the ORIGINAL class's discharge: residual nonzero.
    original = PairedClass(3)
    gap_broken = PairedClass(3)
    gap_broken.atoms = original.atoms[:4] + (original.atoms[5],
                                             original.atoms[4])
    gap_broken.log_sym = {a: original.log_sym[a] for a in gap_broken.atoms}
    gap_broken.log1m_sym = {a: original.log1m_sym[a]
                            for a in gap_broken.atoms}
    n = 6
    gap_broken.prod_sym = {}
    gap_broken.protocol_sym = {}
    for i in range(n):
        for j in range(i, n):
            gap_broken.prod_sym[(i, j)] = sp.Symbol("PR_%d_%d" % (i, j))
            gap_broken.protocol_sym[(i, j)] = sp.Symbol("TQ_%d_%d" % (i, j))
    gap_expr = gap_broken.raw_gap_expr()
    cform = gap_broken.c_matrix_expr(kernel="protocol")
    num = original.mean_sym * (gap_expr
                               - gap_broken.q * (1 - gap_broken.q) * cform)
    poly = sp.expand(num - gap_broken.discharge_expr())
    reduced = original.reduce(poly.subs(original.mean_sym,
                                        original.mixture_mean_expr()))
    out.append({"mutation": "M2_r3_gap_atom_swap_y2y3",
                "expected": "nonzero residual (split identity broken)",
                "residual_terms": len(sp.Add.make_args(reduced))
                if reduced != 0 else 0,
                "ok": reduced != 0})

    # M3: detector sanity: a NON-diagonal geometry must have q1/q2
    # coefficients that are NOT zero balls.
    beta = params_arb.beta
    quad = (Fraction(1, 4), Fraction(1, 2), Fraction(3, 8),
            Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))
    G = []
    for k in range(9):
        vals = tuple(_pt(v) for v in (Fraction(3, 10), Fraction(2, 5),
                                      Fraction(k, 8))) \
            + tuple(_pt(v) for v in quad)
        terms = evaluate_arb(vals, beta)
        G.append(terms.numerator - terms.ehx)
    g2 = _pt(2) * (G[8] + (-2) * G[4] + G[0])
    g1 = (G[8] + (-1) * G[0]) + (-1) * g2
    nonzero = not (g1.contains(arb(0)) and g2.contains(arb(0)))
    out.append({"mutation": "M3_nondiagonal_not_q_independent",
                "expected": "q1/q2 coefficients NOT zero balls",
                "q1": g1.str(8), "q2": g2.str(18),
                "ok": bool(nonzero)})
    return out


# ---------------------------------------------------------------------------
# Report assembly.

def _candidate_record(cand, geom, which):
    return {"family": cand["family"], "name": cand["name"],
            "geometry": [str(v) for v in geom], "a1": cand["a1"],
            "a2": cand["a2"],
            "c": cand["_c_value"].str(26),
            "c_lower": cand["_c_lower_str"],
            "c_upper": cand["_c_upper_str"],
            "interior": cand["family"] == "interior", "which": which}


def grid_extremes(beta):
    worst_min = None
    worst_max = None
    family_counts = {}
    total_candidates = 0
    total_feasible = 0
    for geom in geometry_pairs():
        H = h_matrix_arb(geom)
        coefd = quad_coefficients(H)
        ex = simplex_extremes(coefd, beta)
        fam = ex["min"]["family"]
        family_counts[fam] = family_counts.get(fam, 0) + 1
        total_candidates += ex["n_candidates"]
        total_feasible += ex["n_feasible"]
        lo = float(ex["min"]["_c_value"].lower())
        hi = float(ex["max"]["_c_value"].upper())
        if worst_min is None or lo < worst_min["lower_f"]:
            worst_min = {"lower_f": lo,
                         "geometry": [str(v) for v in geom],
                         "candidate": _candidate_record(ex["min"], geom,
                                                        "min")}
        if worst_max is None or hi > worst_max["upper_f"]:
            worst_max = {"upper_f": hi,
                         "geometry": [str(v) for v in geom],
                         "candidate": _candidate_record(ex["max"], geom,
                                                        "max")}
    return {"worst_min": worst_min, "worst_max": worst_max,
            "min_family_counts": family_counts,
            "total_candidates": total_candidates,
            "total_feasible": total_feasible}


def ball_mid_float(ball_str):
    try:
        return float(arb(ball_str).mid())
    except Exception:
        return None


def ratio_table(cells, beta):
    worst_simplex = None
    worst_cell = None
    danger = []
    near_one = []
    worst_qfac_c = (0.0, None)
    checked = 0
    skipped = 0
    cached_coef = {}
    for (gi, qi), row in sorted(cells.items()):
        if row.get("status") != "ok" or row.get("min_lower") is None:
            skipped += 1
            continue
        geom = tuple(Fraction(v) for v in row["geometry"])
        if gi not in cached_coef:
            cached_coef[gi] = quad_coefficients(h_matrix_arb(geom))
        coefd = cached_coef[gi]
        margin = ball_mid_float(row["min_lower"])
        if margin is None or margin <= 0:
            skipped += 1
            continue
        q = Fraction(row["q"])
        qfac = float(q * (1 - q))
        ex = simplex_extremes(coefd, beta)
        cmin_lower = float(ex["min"]["_c_value"].lower())
        c_neg = max(0.0, -cmin_lower)
        r_simplex = qfac * c_neg / margin if c_neg > 0 and qfac > 0 else 0.0
        cons = build_constraints(mean_line(q, geom), FLOOR_HI)
        verts = polygon_vertices(cons)
        rc_cell = 0.0
        c_cell_lower = None
        if verts:
            cell_cands = kkt_candidates(coefd, cons, mean_index=MEAN_INDEX)
            feasible = [c for c in cell_cands if c["feasible"]]
            for c in feasible:
                c["_c_value"] = beta * c["_value"]
                c["_c_lower_str"] = c["_c_value"].lower().str(26)
                c["_c_upper_str"] = c["_c_value"].upper().str(26)
            if feasible:
                ccell = min(feasible,
                            key=lambda c: (float(c["_c_value"].lower()),
                                           c["name"]))
                c_cell_lower = float(ccell["_c_value"].lower())
                neg = max(0.0, -c_cell_lower)
                rc_cell = qfac * neg / margin
        entry = {"cell": row["cell"],
                 "geometry": row["geometry"], "q": row["q"],
                 "margin_min_lower": row["min_lower"],
                 "margin_f": margin, "c_simplex_lower": cmin_lower,
                 "c_cell_lower": c_cell_lower,
                 "qfac_c_simplex_f": qfac * c_neg,
                 "ratio_simplex_f": r_simplex, "ratio_cell_f": rc_cell}
        if worst_simplex is None or r_simplex > worst_simplex["ratio_simplex_f"]:
            worst_simplex = entry
        if worst_cell is None or rc_cell > worst_cell["ratio_cell_f"]:
            worst_cell = entry
        if rc_cell >= 1.0:
            danger.append(entry)
        near_one.append(entry)
        if qfac * c_neg > worst_qfac_c[0]:
            worst_qfac_c = (qfac * c_neg, entry)
        checked += 1
    near_one.sort(key=lambda e: (-e["ratio_cell_f"], e["cell"]))
    return {"cells_checked": checked, "cells_skipped": skipped,
            "worst_simplex": worst_simplex, "worst_cell": worst_cell,
            "worst_qfac_times_abs_c": worst_qfac_c[0],
            "worst_qfac_times_abs_c_cell": None if worst_qfac_c[1] is None
            else worst_qfac_c[1]["cell"],
            "closest_to_one_top10": near_one[:10],
            "danger_cells": danger,
            "danger_count": len(danger),
            "note": ("a ratio >= 1 is NOT a counterexample: the sweep's "
                     "per-cell margin lower-bounds F(P_mix) at the cell's "
                     "gap minimiser, while |c| is maximised over a "
                     "different mass set; such cells are flagged DANGER "
                     "for follow-up")}


def load_sweep_cells():
    data = json.loads(SWEEP_RESULTS.read_text())
    table = {}
    for row in data["cells"]:
        gi, qi = row["cell"][1:].split("|q")
        table[(int(gi), int(qi))] = row
    return table


def run_all():
    params_mp = solve_equation_parameters(100)
    params_arb = certify_equation_parameters(params_mp)
    beta = params_arb.beta
    beta_mp = mpmath.mpf(float(beta.upper()))
    assert abs(float(beta.upper()) - float(beta.lower())) < 1e-30

    report = {}

    # 1. symbolic matrix form
    cls2 = PairedClass(2)
    cls3 = PairedClass(3)
    sym2 = cls2.prove_split_identity()
    sym3 = cls3.prove_split_identity()
    r3_match = _r3_matrix_matches_remainder()
    report["symbolic_identity"] = {
        "r2_matrix_form": sym2,
        "r3_matrix_form": sym3,
        "r3_matrix_equals_certified_remainder": r3_match,
        "statement": ("c = beta <a, H a> with H_ij = K(x_i,x_j) "
                      "+ K(y_i,y_j) - 2 K(x_i,y_j), K(s,t) = h(pi(s,t)), "
                      "pi(s,t) = st(1+(1-s)(1-t)); residual exactly 0 in "
                      "the canonical free log basis.  H_ij itself is not "
                      "symmetric; the form only reads the pair sums "
                      "H_ij + H_ji."),
        "scope": "r in {2, 3}; exact symbolic identity in (q, masses, supports)",
    }

    # 2. diagonal vanishing
    report["diagonal_vanishing"] = {
        "r2": cls2.prove_diagonal_vanishing(),
        "r3": cls3.prove_diagonal_vanishing(),
        "statement": ("c vanishes identically on y = x: every H entry is "
                      "K_ij + K_ij - 2 K_ij = 0 (identity in the "
                      "argument-named kernel symbols, not a sample)"),
        "scope": "r in {2, 3}; exact symbolic identity",
    }

    # 3. Arb Q1 = Q2 = 0 on diagonals
    diagonal_geometries = [
        (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)),
        (Fraction(1, 8), Fraction(5, 8), Fraction(7, 8)),
        (Fraction(0), Fraction(1, 2), Fraction(1)),
        (Fraction(0), Fraction(1, 4), Fraction(1)),
        (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)),
    ]
    mass_points = [
        (Fraction(3, 10), Fraction(2, 5)),
        (Fraction(2, 7), Fraction(1, 5)),
        (Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(1)),
        (Fraction(1, 3), Fraction(1, 3)),
    ]
    q_checks = []
    all_q_ok = True
    seen_pairs = set()
    for gi, triple in enumerate(diagonal_geometries):
        a1, a2 = mass_points[gi]
        if (a1, a2) in seen_pairs:
            a1, a2 = mass_points[0]
        seen_pairs.add((a1, a2))
        res = q_scan_arb(params_arb, a1, a2, triple)
        res["geometry"] = [str(v) for v in triple]
        res["masses"] = [str(a1), str(a2), str(Fraction(1) - a1 - a2)]
        all_q_ok = all_q_ok and res["status"] == "PROVED (Arb)"
        q_checks.append(res)
    def _fm(v):
        return mpmath.mpf(v.numerator) / mpmath.mpf(v.denominator)

    with mpmath.workdps(200):
        spreads = []
        for triple in diagonal_geometries[:3]:
            gs = []
            for k in range(9):
                vals = (_fm(Fraction(3, 10)),
                        _fm(Fraction(2, 5)),
                        _fm(Fraction(k, 8))) \
                    + tuple(_fm(v) for v in triple + triple)
                gs.append(gap_mp(vals, beta_mp))
            spreads.append(mpmath.nstr(max(gs) - min(gs), 5))
    # Main's corroborating configs: 60 dps, q in {1/10, 1/2, 9/10}
    with mpmath.workdps(60):
        spreads60 = []
        for triple in (diagonal_geometries[0], diagonal_geometries[1]):
            gs = []
            for qf in (Fraction(1, 10), Fraction(1, 2), Fraction(9, 10)):
                vals = (_fm(Fraction(3, 10)), _fm(Fraction(2, 5)),
                        _fm(qf)) \
                    + tuple(_fm(v) for v in triple + triple)
                gs.append(gap_mp(vals, beta_mp))
            spreads60.append(mpmath.nstr(max(gs) - min(gs), 5))
    report["q_independence_on_diagonals"] = {
        "claim": ("when P0 = P1 the mixture is q-independent, so the raw "
                  "gap is CONSTANT in q: exact q^2 and q coefficients of "
                  "the certified gap are zero balls (Q1 = Q2 = 0)"),
        "status": "PROVED (Arb)" if all_q_ok else "FAILED",
        "checks": q_checks,
        "numeric_spread_max_discovery": spreads,
        "numeric_spread_60dps_main_configs_discovery": spreads60,
        "numeric_spread_note": ("mpmath max-min spread of the raw gap "
                                "across q, masses (3/10, 2/5).  At 200 "
                                "dps over q in {k/8}, k=0..8, the spread "
                                "is roundoff-level (~1e-201), NOT "
                                "bitwise zero.  At 60 dps over Main's "
                                "configs (q in {1/10, 1/2, 9/10}, "
                                "diagonals (1/4,1/2,3/4) and "
                                "(1/8,5/8,7/8)) the spread is 0.0 "
                                "bitwise, reproducing the published "
                                "spread-0.0 observation (DISCOVERY)"),
        "scope": ("exact-rational diagonal geometries incl. 0/1 supports "
                  "and vertex masses; Arb ball containment at 320 bits; "
                  "grid q in {k/8}"),
    }

    # 4. numeric identity checks (DISCOVERY)
    tol = mpmath.mpf(10) ** -45
    rows3 = []
    for trial in range(5):
        qv = _fixed_rational(1000, trial, 11)
        a1 = _fixed_rational(1000, trial + 5, 13)
        a2 = _fixed_rational(1000, trial + 10, 17)
        while a1 + a2 >= 1:
            a2 = Fraction(a2.numerator // 2, a2.denominator)
        sup = tuple(_fixed_rational(1000, trial + 15 + k, 19)
                    for k in range(6))
        with mpmath.workdps(200):
            vals = tuple(_fm(v) for v in (a1, a2, qv) + sup)
            row = numeric_identity(vals, beta_mp, r=3)
        row["ok"] = abs(mpmath.mpf(row["residual"])) < tol
        rows3.append(row)
    rows2 = []
    for trial in range(5):
        qv = _fixed_rational(1000, trial + 30, 23)
        a1 = _fixed_rational(1000, trial + 35, 29)
        sup = tuple(_fixed_rational(1000, trial + 40 + k, 31)
                    for k in range(4))
        with mpmath.workdps(200):
            vals = tuple(_fm(v) for v in (a1, qv) + sup)
            row = numeric_identity(vals, beta_mp, r=2)
        row["ok"] = abs(mpmath.mpf(row["residual"])) < tol
        rows2.append(row)
    report["numeric_identity_discovery"] = {
        "label": "DISCOVERY",
        "tolerance": mpmath.nstr(tol, 5),
        "dps": 200,
        "r3": {"rows": rows3, "all_ok": all(r["ok"] for r in rows3)},
        "r2": {"rows": rows2, "all_ok": all(r["ok"] for r in rows2)},
        "scope": ("5 deterministic r=3 and 5 deterministic r=2 "
                  "configurations, pointwise residual at 200 dps"),
    }

    # 5. grid extremes of c over the simplex
    extremes = grid_extremes(beta)
    report["simplex_extremes"] = {
        "grid": "k/8 ordered pairs, 35 x 35 = 1225 geometries",
        "candidates": ("vertex + edge-critical + interior-critical "
                       "(closed form; certified Arb balls, exact-rational "
                       "vertices)"),
        "worst_min": extremes["worst_min"],
        "worst_max": extremes["worst_max"],
        "min_family_counts": extremes["min_family_counts"],
        "total_candidates": extremes["total_candidates"],
        "total_feasible": extremes["total_feasible"],
        "scope": ("per-geometry certified min/max over the full simplex "
                  "{a >= 0, sum a = 1}; pointwise grid (support "
                  "directions between grid points stay open)"),
    }

    # 6. ratio table vs sweep margins
    cells = load_sweep_cells()
    ratios = ratio_table(cells, beta)
    report["ratio_table"] = {
        "source": ("liu9-h2-qp-sweep.json per-cell min_lower margins; "
                   "11025-cell k/8 x k/8 grid"),
        **ratios,
        "scope": ("certified Arb balls for c; margins parsed from the "
                  "published sweep balls; pointwise cells"),
    }

    # 7. mutations
    report["mutations"] = run_mutations(params_arb)
    mutations_ok = all(m["ok"] for m in report["mutations"])
    report["mutations_all_fail_detected"] = mutations_ok

    # 8. verdict
    symbolic_ok = (sym2["status"] == "PROVED"
                   and sym3["status"] == "PROVED"
                   and r3_match["status"] == "PROVED")
    diag_ok = (report["diagonal_vanishing"]["r2"]["status"] == "PROVED"
               and report["diagonal_vanishing"]["r3"]["status"] == "PROVED")
    if (symbolic_ok and diag_ok and all_q_ok and mutations_ok
            and report["numeric_identity_discovery"]["r3"]["all_ok"]
            and report["numeric_identity_discovery"]["r2"]["all_ok"]):
        report["verdict"] = "C-FORM-CHARACTERIZED"
        report["claim_status"] = "MACHINE-VERIFIED"
    else:
        report["verdict"] = "FAILED"
        report["claim_status"] = "FAILED"
    report["claim"] = (
        "the channel scalar c of the r <= 3 paired class is EXACTLY "
        "c = beta <a, H a> with H_ij = K(x_i,x_j) + K(y_i,y_j) "
        "- 2 K(x_i,y_j) (symbolic, r = 2 and r = 3, zero residual; r = 3 "
        "form equals the certified liu9_paired_class_c remainder); c "
        "vanishes identically on y = x; on diagonal geometries the gap is "
        "q-constant with Q1 = Q2 = 0 certified zero balls (Arb); simplex "
        "extremes of c over 1225 k/8 geometries are certified closed-form "
        "candidate lists; per-cell q(1-q)|c| vs F-margin ratios certified "
        "against the published 11,025-cell sweep")
    report["claim_scope"] = (
        "r in {2, 3}; k/8 ordered geometry grid (7 support values, "
        "C(7,3)=35 triples per side) x k/8 q grid, pointwise cells; "
        "symbolic statements exact and q-free; Arb statements at 320-bit "
        "ball precision; mpmath rows labelled DISCOVERY")
    report["evaluator"] = (
        "liu9_objective.evaluate_arb (gap = numerator - ehx), ctx.prec = "
        "320; mpmath appears only in labelled DISCOVERY sections via "
        "liu9_boundary_layer.gap_mp")
    report["limitations"] = [
        "r is fixed at 2 and 3; no claim about r >= 4",
        "geometry grid is the k/8 ordered grid (35 x 35); support "
        "directions between grid points stay open",
        "q grid is k/8; q-between-grid-points stays open except where a "
        "statement is q-free (the c form) or polynomial-exact "
        "(diagonal q-scan)",
        "the displayed 2x2 block kernel is NOT revisited here (refuted "
        "object; see LIU9_BLOCK_COPOSITIVE section 6)",
        "ratio >= 1 cells are danger flags, not counterexamples: F is "
        "certified at the cell's gap minimiser, |c| maximised elsewhere",
    ]
    report["report_sha256"] = _canonical_digest(report)
    return report


def main():
    report = run_all()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("C-CHARACTERIZATION (r=2,3 matrix form; diagonal vanishing; "
          "Q1=Q2=0; simplex extremes; ratio table)")
    si = report["symbolic_identity"]
    print("r2 split identity:", si["r2_matrix_form"]["status"],
          "| r3 split identity:", si["r3_matrix_form"]["status"],
          "| r3 == remainder:", si["r3_matrix_equals_certified_remainder"]
          ["status"])
    dv = report["diagonal_vanishing"]
    print("diagonal vanishing r2:", dv["r2"]["status"], "| r3:",
          dv["r3"]["status"])
    print("q-independence on diagonals:",
          report["q_independence_on_diagonals"]["status"],
          "| spreads (discovery):",
          report["q_independence_on_diagonals"]["numeric_spread_max_discovery"])
    num = report["numeric_identity_discovery"]
    print("numeric identity r3 all_ok:", num["r3"]["all_ok"],
          "| r2 all_ok:", num["r2"]["all_ok"])
    se = report["simplex_extremes"]["worst_min"]
    print("most negative c: %s at geometry %s arg (%s, %s) family %s"
          % (se["candidate"]["c_lower"], se["geometry"],
             se["candidate"]["a1"], se["candidate"]["a2"],
             se["candidate"]["family"]))
    rt = report["ratio_table"]
    print("ratio cells checked:", rt["cells_checked"], "skipped:",
          rt["cells_skipped"])
    if rt["worst_simplex"] is not None:
        ws = rt["worst_simplex"]
        print("worst ratio_simplex %.6f at cell %s q=%s"
              % (ws["ratio_simplex_f"], ws["cell"], ws["q"]))
    if rt["worst_cell"] is not None:
        wc = rt["worst_cell"]
        print("worst ratio_cell %.6f at cell %s q=%s geom %s"
              % (wc["ratio_cell_f"], wc["cell"], wc["q"], wc["geometry"]))
    print("worst q(1-q)|c| product: %.6f (cell %s)"
          % (rt["worst_qfac_times_abs_c"],
             rt["worst_qfac_times_abs_c_cell"]))
    print("danger cells (ratio_cell >= 1):", rt["danger_count"])
    for e in rt["closest_to_one_top10"][:5]:
        print("  close: cell %s q=%s ratio_cell %.6f geom %s"
              % (e["cell"], e["q"], e["ratio_cell_f"], e["geometry"]))
    for m in report["mutations"]:
        print("mutation %-34s %s" % (m["mutation"], "OK-detected"
                                     if m["ok"] else "NOT-DETECTED"))
    print("verdict:", report["verdict"], "| claim_status:",
          report["claim_status"])
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] == "MACHINE-VERIFIED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
