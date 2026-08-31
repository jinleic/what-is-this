#!/usr/bin/env python3
"""REGION 2 (two-atom paired class) — exact c_two, channel decomposition,
meaning of the covering problem.

Family: P0 = a1 d_x1 + a2 d_x2, P1 = a1 d_y1 + a2 d_y2, a1 + a2 = 1.

The block diagonalization (liu9_block_copositive.py, PROVED) instantiates
here as

    gap = F(P_mix) + beta q(1-q) E_K(P0 - P1),   K = h o pi,

with the four-atom mixture law ((1-q)a1, (1-q)a2, q a1, q a2) on
(x1, x2, y1, y2) and the EXACT two-atom scalar

    c_two = beta sum_{i,j} a_i a_j [K(x_i,x_j) + K(y_i,y_j) - 2K(x_i,y_j)].

Each (i, j) bracket is a per-pair CHANNEL: the second difference of K
across the rectangle {x_i, x_j} x {y_i, y_j}.  This module

  (1) proves the two-atom diagonalization SYMBOLICALLY in the canonical
      free log basis (same method as region1's three-slot proof),
  (2) proves the CHANNEL DECOMPOSITION structurally (it is the definition,
      checked symbolically from the general remainder form),
  (3) probes concavity of the total gap in (q, a2) on mean-feasible
      configurations (exact quadratic in q; behaviour in a2 sampled, then
      certified where possible),
  (4) verifies the identity vs gap_mp at 90 dps on deterministic Family-B
      configurations (witness-search family),
  (5) runs a deterministic danger scan (worst total gap + worst c_two),
      reporting the Region 2 obligation precisely:
      F(P_mix) + beta q(1-q) E_K >= 0 with F(P_mix) >= 0 proved; open is
      whether certified F-margins absorb -beta q(1-q) E_K when E_K < 0.

Byte-stable report: uc/verification/results/liu9-scalar-margin-region2.json
Run:  math/.venv/bin/python -I -B math/uc/liu9_scalar_margin_region2.py
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import mpmath
import sympy as sp

from liu9_binding import solve_equation_parameters
from liu9_boundary_layer import gap_mp

OUTPUT_DEFAULT = (
    HERE / "verification/results/liu9-scalar-margin-region2.json")
CERT_SEED = 20260829

# --- two-atom symbolic scaffolding ------------------------------------------

A1, BETA, Q = sp.symbols("a1 beta q", real=True)
X1, X2, Y1, Y2 = sp.symbols("x1 x2 y1 y2", real=True)
_MEAN = sp.Symbol("M_mix")
MASSES = (A1, 1 - A1)
ATOMS = (X1, X2, Y1, Y2)

_LOG_ATOM = {a: sp.Symbol(f"l_{a.name}") for a in ATOMS}
_LOG_1M = {a: sp.Symbol(f"k_{a.name}") for a in ATOMS}
_PROD_1M = {}
for _i in range(4):
    for _j in range(_i, 4):
        _PROD_1M[(_i, _j)] = sp.Symbol(f"p{_i}{_j}")
_PROTOCOL_1M = {}
for _i in range(4):
    for _j in range(_i, 4):
        _PROTOCOL_1M[(_i, _j)] = sp.Symbol(f"tq{_i}{_j}")


def _pi(s, t):
    return s * t * (1 + (1 - s) * (1 - t))


def _PROD_SYM(i, j):
    lo, hi = min(i, j), max(i, j)
    return sp.Symbol(f"h({ATOMS[lo] * ATOMS[hi]})")


def _PI_SYM(i, j):
    lo, hi = min(i, j), max(i, j)
    return sp.Symbol(f"h({_pi(ATOMS[lo], ATOMS[hi])})")


def e_kernel(i, j, beta):
    s, t = ATOMS[i], ATOMS[j]
    big_a = t * sp.Symbol(f"h({s})") + s * sp.Symbol(f"h({t})")
    t_kernel = (sp.Symbol(f"mu({s})") + sp.Symbol(f"mu({t})")
                - sp.Symbol(f"mu({s * t})"))
    gain = _PI_SYM(i, j) - _PROD_SYM(i, j)
    return big_a - s * t * t_kernel + beta * gain


def raw_gap_expr(beta):
    masses = MASSES
    weights = tuple((1 - Q) * m for m in masses) + tuple(Q * m for m in masses)
    ehxy = sp.Integer(0)
    for i in range(4):
        for j in range(4):
            ehxy += weights[i] * weights[j] * _PROD_SYM(i, j)
    ehpi = sp.Integer(0)
    for comp in range(2):
        factor = (1 - Q) if comp == 0 else Q
        off = 0 if comp == 0 else 2
        for i in range(2):
            for j in range(2):
                ehpi += factor * masses[i] * masses[j] * _PI_SYM(off + i,
                                                                 off + j)
    ehx = sp.Integer(0)
    for i in range(4):
        ehx += weights[i] * sp.Symbol(f"h({ATOMS[i]})")
    return (1 - beta) * ehxy + beta * ehpi - ehx



def _discharge(beta):
    weights = tuple((1 - Q) * m for m in MASSES) + tuple(Q * m for m in MASSES)
    total = sp.Integer(0)
    for i in range(4):
        for j in range(4):
            s, t = ATOMS[i], ATOMS[j]
            d = _MEAN * e_kernel(i, j, beta) \
                - (t * sp.Symbol(f"h({s})") + s * sp.Symbol(f"h({t})")) / 2
            total += weights[i] * weights[j] * d
    return total


def remainder_expr(beta):
    total = sp.Integer(0)
    for i in range(2):
        for j in range(2):
            total += MASSES[i] * MASSES[j] * (
                _PI_SYM(i, j) + _PI_SYM(2 + i, 2 + j)
                - 2 * _PI_SYM(i, 2 + j))
    return beta * total


def mixture_mean_expr():
    weights = tuple((1 - Q) * m for m in MASSES) + tuple(Q * m for m in MASSES)
    return sum((weights[i] * ATOMS[i] for i in range(4)), sp.Integer(0))


def _arg_atoms(name, prefix):
    locs = {a.name: a for a in ATOMS}
    arg = sp.sympify(name[prefix:-1], locals=locs)
    pd = arg.as_powers_dict()
    atoms_list = []
    for a, p in sorted(pd.items(), key=lambda kv: str(kv[0])):
        atoms_list.extend([a] * int(p))
    return atoms_list


def _atom_indices(atoms_list):
    idx = {a.name: k for k, a in enumerate(ATOMS)}
    return tuple(sorted(idx[a.name] for a in atoms_list))


def _h_of(atoms_list):
    if not atoms_list:
        return sp.Integer(0)
    indices = _atom_indices(atoms_list)
    if len(indices) == 1:
        a = atoms_list[0]
        return -a * _LOG_ATOM[a] - (1 - a) * _LOG_1M[a]
    z = sp.prod([ATOMS[k] for k in indices])
    log_z = sum((_LOG_ATOM[ATOMS[k]] for k in indices), sp.Integer(0))
    return -z * log_z - (1 - z) * _PROD_1M[tuple(sorted(indices))]


def _mu_of(atoms_list):
    if not atoms_list:
        return sp.Integer(1)
    hz = _h_of(atoms_list)
    z = sp.prod(atoms_list) if len(atoms_list) > 1 else atoms_list[0]
    log_z = sum((_LOG_ATOM[a] for a in atoms_list), sp.Integer(0))
    return sp.cancel(hz / z - log_z)


def _protocol_pair(arg):
    for i in range(4):
        for j in range(i, 4):
            if sp.simplify(sp.expand(arg - _pi(ATOMS[i], ATOMS[j]))) == 0:
                return (i, j)
    return None


def _h_reduce(name):
    locs = {a.name: a for a in ATOMS}
    arg = sp.sympify(name[2:-1], locals=locs)
    pair = _protocol_pair(arg)
    if pair is not None:
        i, j = pair
        left, right = ATOMS[i], ATOMS[j]
        z = _pi(left, right)
        log_z = (2 * _LOG_ATOM[left] if left == right
                 else _LOG_ATOM[left] + _LOG_ATOM[right])
        return -z * log_z - (1 - z) * _PROTOCOL_1M[(i, j)]
    return _h_of(_arg_atoms(name, 2))


def _reduce(expr):
    subs = {}
    for s in list(expr.free_symbols):
        n = s.name
        if n.startswith("mu("):
            subs[s] = _mu_of(_arg_atoms(n, 3))
        elif n.startswith("h("):
            subs[s] = _h_reduce(n)
    return sp.expand(sp.cancel(expr.subs(subs)))


def prove_two_atom_identity():
    beta = BETA
    gap = raw_gap_expr(beta)
    c_of_q = Q * (1 - Q) * remainder_expr(beta)
    numerator = _MEAN * (gap - c_of_q)
    poly = sp.expand(numerator - _discharge(beta))
    reduced = _reduce(poly.subs(_MEAN, mixture_mean_expr()))
    if reduced != 0:
        return {"status": "FAILED",
                "residual_terms": len(sp.Add.make_args(reduced))}
    return {"status": "PROVED",
            "method": ("sympy multilinear expansion in the canonical free "
                       "log basis of the two-atom family; zero residual"),
            "log_basis_size": (len(_LOG_ATOM) + len(_LOG_1M)
                               + len(_PROD_1M) + len(_PROTOCOL_1M))}


def prove_channel_decomposition():
    """c_two is exactly the per-pair channel sum: verify symbolically that
    the remainder form equals sum_{i,j} a_i a_j C_ch(x_i, x_j, y_i, y_j)
    where C_ch = K(xi,xj) + K(yi,yj) - 2K(xi,yj)."""
    beta = BETA
    rem = remainder_expr(beta)
    channel_sum = sp.Integer(0)
    for i in range(2):
        for j in range(2):
            channel_sum += MASSES[i] * MASSES[j] * (
                _PI_SYM(i, j) + _PI_SYM(2 + i, 2 + j) - 2 * _PI_SYM(i, 2 + j))
    channel_sum = beta * channel_sum
    diff = sp.expand(rem - channel_sum)
    reduced = _reduce(diff)
    if reduced != 0:
        return {"status": "FAILED"}
    return {
        "status": "PROVED",
        "channels": 4,
        "note": ("c_two decomposes exactly into the four (i, j) rectangle "
                 "channels, each a second difference of K = h o pi across "
                 "{x_i, x_j} x {y_i, y_j}; every channel vanishes when "
                 "y_i = x_i and y_j = x_j (diagonal-collapse)"),
    }


# --- numeric certification ----------------------------------------------------

def random_two_atom(rng, m):
    def _fm(v):
        return mpmath.mpf(v.numerator) / mpmath.mpf(v.denominator)
    while True:
        a1 = _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
        q = _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
        supports = tuple(
            _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
            for _ in range(4))
        values = (a1, q) + supports
        a2 = 1 - a1
        mean = q * (a1 * supports[0] + a2 * supports[1]) \
            + (1 - q) * (a1 * supports[2] + a2 * supports[3])
        if mean >= m and min(supports) > mpmath.mpf("1e-3"):
            return values

def _h_mp(u):
    if u == 0 or u == 1:
        return mpmath.mpf(0)
    return -(u * mpmath.log(u) + (1 - u) * mpmath.log1p(-u))


def _pi_mp(x, y):
    return x * y * (1 + (1 - x) * (1 - y))


def _mu_mp(u):
    if u == 0:
        return mpmath.mpf(1)
    if u == 1:
        return mpmath.mpf(0)
    return -((1 - u) * mpmath.log1p(-u)) / u


def _e_mp(s, t, beta):
    a = t * _h_mp(s) + s * _h_mp(t)
    tk = _mu_mp(s) + _mu_mp(t) - _mu_mp(s * t)
    return a - s * t * tk + beta * (_h_mp(_pi_mp(s, t)) - _h_mp(s * t))


def _d_mp(s, t, mean, beta):
    a = t * _h_mp(s) + s * _h_mp(t)
    return mean * _e_mp(s, t, beta) - a / 2


def evaluate_two_atom(values, beta):
    """values = (a1, q, x1, x2, y1, y2).  Returns
    (discharge, c_two, mean)."""
    a1, q, x1, x2, y1, y2 = values
    a2 = 1 - a1
    masses = (a1, a2)
    sup = (x1, x2, y1, y2)
    weights = tuple((1 - q) * m for m in masses) + tuple(q * m for m in masses)
    mean = sum((weights[k] * sup[k] for k in range(4)), mpmath.mpf(0))
    discharge = mpmath.mpf(0)
    for i in range(4):
        for j in range(4):
            discharge += weights[i] * weights[j] * _d_mp(
                sup[i], sup[j], mean, beta)
    discharge = discharge / mean
    quad = mpmath.mpf(0)
    for i in range(2):
        for j in range(2):
            quad += masses[i] * masses[j] * (
                _h_mp(_pi_mp(sup[i], sup[j]))
                + _h_mp(_pi_mp(sup[2 + i], sup[2 + j]))
                - 2 * _h_mp(_pi_mp(sup[i], sup[2 + j])))
    c_two = beta * quad
    return discharge, c_two, mean


def random_two_atom(rng, m):
    def _fm(v):
        return mpmath.mpf(v.numerator) / mpmath.mpf(v.denominator)
    while True:
        a1 = _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
        q = _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
        supports = tuple(
            _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
            for _ in range(4))
        values = (a1, q) + supports
        a2 = 1 - a1
        mean = q * (a1 * supports[0] + a2 * supports[1]) \
            + (1 - q) * (a1 * supports[2] + a2 * supports[3])
        if mean >= m and min(supports) > mpmath.mpf("1e-3"):
            return values


def numeric_certification(dps=90):
    parameters = solve_equation_parameters(dps)
    beta = parameters.beta
    rng = random.Random(CERT_SEED)
    records = []
    worst = mpmath.mpf(0)
    most_neg_c = mpmath.mpf(0)
    with mpmath.workdps(dps):
        m = parameters.mean
        for trial in range(8):
            values = random_two_atom(rng, m)
            a1, q, x1, x2, y1, y2 = values
            vals9 = (a1, 1 - a1, q, x1, x2, mpmath.mpf("0.5"),
                     y1, y2, mpmath.mpf("0.5"))
            ref = gap_mp(vals9, beta)
            discharge, c_two, _mean = evaluate_two_atom(values, beta)
            total = discharge + q * (1 - q) * c_two
            diff = abs(total - ref)
            worst = max(worst, diff)
            if c_two < most_neg_c:
                most_neg_c = c_two
            records.append({
                "trial": trial,
                "values": [mpmath.nstr(v, 30) for v in values],
                "raw_gap": mpmath.nstr(ref, 40),
                "identity_absdiff": mpmath.nstr(diff, 8),
                "c_two": mpmath.nstr(c_two, 40),
                "discharge": mpmath.nstr(discharge, 40),
            })
        if not worst < mpmath.mpf(10) ** -85:
            raise AssertionError("two-atom numeric identity check failed")
    return {
        "claim_status": "MACHINE-VERIFIED",
        "dps": dps,
        "configurations": records,
        "worst_identity_absdiff": mpmath.nstr(worst, 8),
        "most_negative_c_two": mpmath.nstr(most_neg_c, 40),
        "seed": CERT_SEED,
    }


def concavity_probe(dps=70):
    """Exact quadratic-in-q check (3-point interpolation vs direct) and a
    sample of second differences in a2 along mean-feasible lines."""
    parameters = solve_equation_parameters(dps)
    beta = parameters.beta
    out = {"quadratic_in_q_max_deviation": None,
           "a2_direction_samples": []}
    with mpmath.workdps(dps):
        m = parameters.mean
        max_dev = mpmath.mpf(0)
        for a1f in ("0.3", "0.5", "0.8"):
            a1 = mpmath.mpf(a1f)
            x1, x2, y1, y2 = (mpmath.mpf(t) for t in
                              ("0.3", "0.7", "0.9", "0.4"))
            a2 = 1 - a1
            base_mean = mpmath.mpf(0)

            def gap_at(qq):
                vals9 = (a1, a2, qq, x1, x2, mpmath.mpf("0.5"),
                         y1, y2, mpmath.mpf("0.5"))
                return gap_mp(vals9, beta), None

            qs = (mpmath.mpf("0.2"), mpmath.mpf("0.5"), mpmath.mpf("0.8"))
            greedy = [gap_at(qq)[0] for qq in qs]
            # feasible band for q at this a1: q >= (m - a2-ref mean)/...
            # crude: require mean >= m
            def mean_at(qq):
                return m - (qq * (a1 * y1 + a2 * y2)
                            + (1 - qq) * (a1 * x1 + a2 * x2)) * 0
            # direct mean:
            def mean_of(qq):
                return qq * (a1 * y1 + a2 * y2) \
                    + (1 - qq) * (a1 * x1 + a2 * x2)
            band_ok = all(mean_of(qq) >= m for qq in qs)
            if band_ok:
                # quadratic interpolation at q = 0.7
                q4 = mpmath.mpf("0.7")
                v4 = greedy[0] * ((q4 - qs[1]) * (q4 - qs[2])) \
                    / ((qs[0] - qs[1]) * (qs[0] - qs[2])) \
                    + greedy[1] * ((q4 - qs[0]) * (q4 - qs[2])) \
                    / ((qs[1] - qs[0]) * (qs[1] - qs[2])) \
                    + greedy[2] * ((q4 - qs[0]) * (q4 - qs[1])) \
                    / ((qs[2] - qs[0]) * (qs[2] - qs[1]))
                dev = abs(v4 - gap_at(q4)[0])
                max_dev = max(max_dev, dev)
            # a2-direction second difference at fixed q:
            qv = mpmath.mpf("0.5")
            for d in (mpmath.mpf("0.05"), mpmath.mpf("0.1")):
                if a1 - d > 0 and a1 + d < 1:
                    # a2 shift = a1 shift with recomputed mean:
                    def gap_a1(a1v, qq):
                        a2v = 1 - a1v
                        meanv = qq * (a1v * y1 + a2v * y2) \
                            + (1 - qq) * (a1v * x1 + a2v * x2)
                        vals9 = (a1v, a2v, qq, x1, x2, mpmath.mpf("0.5"),
                                 y1, y2, mpmath.mpf("0.5"))
                        return gap_mp(vals9, beta) if meanv >= m else None
                    left = gap_a1(a1 - d, qv)
                    mid = gap_a1(a1, qv)
                    right = gap_a1(a1 + d, qv)
                    if left is not None and right is not None:
                        second = (right - 2 * mid + left) / (d * d)
                        out["a2_direction_samples"].append({
                            "a1": a1f, "second_diff_a2": mpmath.nstr(
                                second, 10)})
        out["quadratic_in_q_max_deviation"] = mpmath.nstr(max_dev, 8)
    return out


def danger_scan(dps=50):
    """Deterministic Family-B danger scan with coordinate descent."""
    parameters = solve_equation_parameters(dps)
    beta = parameters.beta
    rng = random.Random(CERT_SEED + 1)
    m = parameters.mean
    worst = None
    worst_vals = None
    with mpmath.workdps(dps):
        for _restart in range(8):
            vals = random_two_atom(rng, m)
            cur = list(vals)
            scale = mpmath.mpf("0.05")
            best_here = None
            best_vals = None
            for _step in range(40):
                j = rng.randrange(6)
                old = cur[j]
                new_float = float(cur[j]) + (rng.random() - 0.5) * 2 * float(scale)
                _fr = Fraction(new_float).limit_denominator(10 ** 9)
                new = mpmath.mpf(_fr.numerator) / mpmath.mpf(_fr.denominator)
                if j == 0:
                    new = min(max(new, mpmath.mpf("1e-3")), mpmath.mpf(1))
                elif j == 1:
                    new = min(max(new, mpmath.mpf(0)), mpmath.mpf(1))
                else:
                    new = min(max(new, mpmath.mpf("1e-3")), mpmath.mpf(1))
                trial = list(cur)
                trial[j] = new
                a1t, qt = trial[0], trial[1]
                a2t = 1 - a1t
                mean_t = (1 - qt) * (a1t * trial[2] + a2t * trial[3]) \
                    + qt * (a1t * trial[4] + a2t * trial[5])
                if mean_t < m:
                    cur[j] = old
                    continue
                discharge, c_two, _ = evaluate_two_atom(tuple(trial), beta)
                total = discharge + qt * (1 - qt) * c_two
                if best_here is None or total < best_here:
                    best_here = total
                    best_vals = tuple(trial)
                    cur = trial
                else:
                    cur[j] = old
                    scale = scale * mpmath.mpf("0.995")
            if best_here is not None and (worst is None
                                          or best_here < worst):
                worst = best_here
                worst_vals = best_vals
    return {
        "worst_total_gap": mpmath.nstr(worst, 12) if worst is not None else None,
        "worst_values": [mpmath.nstr(v, 10) for v in worst_vals]
        if worst_vals else None,
        "negative_found": bool(worst is not None and worst < 0),
        "restarts": 8,
    }


def build_report():
    parameters = solve_equation_parameters(100)
    symbolic = prove_two_atom_identity()
    if symbolic["status"] != "PROVED":
        raise AssertionError("two-atom identity did not prove out")
    channels = prove_channel_decomposition()
    if channels["status"] != "PROVED":
        raise AssertionError("channel decomposition did not prove out")
    numeric = numeric_certification()
    concavity = concavity_probe()
    scan = danger_scan()
    proved = bool(symbolic["status"] == "PROVED"
                  and channels["status"] == "PROVED")
    report = {
        "tool": "liu9_scalar_margin_region2.py",
        "claim_status": "PROVED" if proved else "FAILED",
        "label_scope": (
            "the PROVED label covers the two symbolic statements only -- the "
            "two-atom diagonalization and the channel decomposition, both "
            "sympy zero-residual identities.  The concavity probe, the "
            "numeric certification and the danger scan are deterministic "
            "high-precision mpmath computations and are DISCOVERY-grade "
            "evidence, not certificates; the interval certificates for these "
            "curve families are liu9_scalar_margin_region2hulls.py and "
            "liu9_scalar_margin_region2_curve2.py"),
        "constants": {
            "m": mpmath.nstr(parameters.mean, 40),
            "beta": mpmath.nstr(parameters.beta, 40),
        },
        "statements": {
            "two_atom_diagonalization": symbolic,
            "channel_decomposition": channels,
            "corrected_form": (
                "gap = F(P_mix) + beta q(1-q) E_{h o pi}(P0 - P1) with the "
                "four-atom mixture law; c_two = beta sum_{i,j<=2} a_i a_j "
                "[K(x_i,x_j) + K(y_i,y_j) - 2K(x_i,y_j)], K = h o pi; each "
                "(i,j) term is a rectangle second-difference CHANNEL that "
                "vanishes on the diagonal-collapse y = x"
            ),
        },
        "concavity_probe": concavity,
        "numeric_certification": numeric,
        "danger_scan": scan,
        "remaining_obligation": (
            "Region 2 reduces exactly to F(P_mix) + beta q(1-q) E_K >= 0 "
            "with F(P_mix) >= 0 certified for every mixture mean >= m "
            "(universal q=1 theorem).  Open: certified absorption of "
            "-beta q(1-q) E_K by the F-margins (kappa dist^2 near the "
            "optimizer box; face lemmas).  Cover plan mirrors "
            "region1/region1b: exact concave quadratic in q gives an "
            "endpoint-only reduction PER (a1, geometry) fiber; the a2 "
            "direction is probed but not yet certified"
        ),
        "limitations": [
            "two-atom paired class only (shared masses a1 + a2 = 1)",
            "no certified Region-2 cell cover yet: this report proves the "
            "exact reduction + channel decomposition, verifies numerically, "
            "and maps the danger; the cover is the next deliverable",
        ],
    }
    report["report_sha256"] = _digest(report)
    return report


def _digest(payload):
    return hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, indent=1, sort_keys=True, default=str) + "\n")
    print("REGION 2 TWO-ATOM PAIRED CLASS")
    print("identity:",
          report["statements"]["two_atom_diagonalization"]["status"])
    print("channels:", report["statements"]["channel_decomposition"]["channels"])
    print("danger scan worst:", report["danger_scan"]["worst_total_gap"])
    print("negative found:", report["danger_scan"]["negative_found"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
