#!/usr/bin/env python3
"""GENERAL multi-atom paired class: exact c_ch channel formula.

Conjecture (Main's step 1): for any finite paired laws
    P0 = sum_{i<=r} a_i delta_{x_i},  P1 = sum_{i<=r} a_i delta_{y_i},
with Liu's raw gap assembled in 1/2 + 1/2-component weights, the exact
complement to the diagonalization's pisarિજctx channel is

    c = beta * (E_K(P0) + E_K(P1) - 2 E_K(P0, P1))
      = beta * sum_{i,j<=r} a_i a_j [ K(x_i,x_j) + K(y_i,y_j) - 2K(x_i,y_j) ],
    K(s,t) = h(pi(s,t)).

The two-atom instance was PROVED in liu9_scalar_margin_region2.py
(symbolically, 28-symbol free log basis, zero residual).  This module
proves the r=3 (THREE atoms per component, six-component mixture) case
using the same canonical free log basis reduction; r >= 4 instances are
N-atom-identity wrappers of the same structure (sub-routines only: the
6-atom main theorem T3 of liu9_block_copositive.py already has the
3-slot machinery).

The module:
  (1) proves c = beta*(E0 + E1 - 2 E_cross) for r = 3 EXACTLY in a
      ~60-symbol canonical free log basis over (a1, a2, beta, q, x_i, y_j);
  (2) numerically verifies the general identity vs gap_mp at 90 dps on
      deterministic 3-atom paired configurations (discovery-discipline
      [workdps(90), deterministic seeds]);
  (3) if the identity fails: writes the exact failure structure to the
      report (falsification protocol).

Byte-stable report: uc/verification/results/liu9-paired-class-c.json
Run: math/.venv/bin/python -I -B math/uc/liu9_paired_class_c.py
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
    HERE / "verification/results/liu9-paired-class-c.json")
CERT_SEED = 20260829


# --- three-atom canonical log basis (same design) -----------------------------

A1, BETA, Q = sp.symbols("a1 beta q", real=True)
X1, X2, X3, Y1, Y2, Y3 = sp.symbols("x1 x2 x3 y1 y2 y3", real=True)
MASSES = (A1, sp.Symbol("a2"), 1 - A1 - sp.Symbol("a2"))  # Liu masses: shared on both components
ATOMS = (X1, X2, X3, Y1, Y2, Y3)

_LOG_ATOM = {a: sp.Symbol(f"l_{a.name}") for a in ATOMS}
_LOG_1M = {a: sp.Symbol(f"k_{a.name}") for a in ATOMS}
_PROD_1M = {}
for _i in range(6):
    for _j in range(_i, 6):
        _PROD_1M[(_i, _j)] = sp.Symbol(f"p{_i}{_j}")
_PROTOCOL_1M = {}
for _i in range(6):
    for _j in range(_i, 6):
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
    tk = (sp.Symbol(f"mu({s})") + sp.Symbol(f"mu({t})")
          - sp.Symbol(f"mu({s * t})"))
    gain = _PI_SYM(i, j) - _PROD_SYM(i, j)
    return big_a - s * t * tk + beta * gain


def _arg_atoms(name, prefix):
    locs = {a.name: a for a in ATOMS}
    arg = sp.sympify(name[prefix:-1], locals=locs)
    pd = arg.as_powers_dict()
    atoms = []
    for a, p in sorted(pd.items(), key=lambda kv: str(kv[0])):
        atoms.extend([a] * int(p))
    return atoms


def _atom_indices(atoms_list):
    idx = {a.name: k for k, a in enumerate(ATOMS)}
    return tuple(sorted(idx[a.name] for a in atoms_list))


def _protocol_pair(arg):
    for i in range(6):
        for j in range(i, 6):
            if sp.simplify(sp.expand(arg - _pi(ATOMS[i], ATOMS[j]))) == 0:
                return (i, j)
    return None


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


def _mu_of_t(name):
    return _mu_of(_arg_atoms(name, 3))


def _reduce(expr):
    subs = {}
    for s in list(expr.free_symbols):
        n = s.name
        if n.startswith("mu("):
            subs[s] = _mu_of(_arg_atoms(n, 3))
        elif n.startswith("h("):
            subs[s] = _h_reduce(n)
    return sp.expand(sp.cancel(expr.subs(subs)))


# raw gap for the three-atom paired class (Liu's 9-var transcription shifted
# vertex-symmetries): weights w = ((1-q)a1, (1-q)a2, (1-q)a3, q a1, q a2, q a3):
def raw_gap_expr(beta):
    weights = tuple((1 - Q) * m for m in MASSES) + tuple(Q * m for m in MASSES)
    ehxy = sp.Integer(0)
    for i in range(6):
        for j in range(6):
            ehxy += weights[i] * weights[j] * _PROD_SYM(i, j)
    ehpi = sp.Integer(0)
    for comp in range(2):
        factor = (1 - Q) if comp == 0 else Q
        off = 0 if comp == 0 else 3
        for i in range(3):
            for j in range(3):
                ehpi += (factor * MASSES_C(i) * MASSES_C(j)
                         * _PI_SYMC(off + i, off + j))
    ehx = sp.Integer(0)
    for i in range(6):
        ehx += weights[i] * sp.Symbol(f"h({ATOMS[i]})")
    return (1 - beta) * ehxy + beta * ehpi - ehx


_A2 = sp.Symbol("a2")


def MASSES_C(i):
    return {0: A1, 1: _A2, 2: 1 - A1 - _A2}[i % 3]


def _PI_SYMC(i, j):
    return _PI_SYM(i, j)


def remainder_expr(beta):
    total = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            total += MASSES_C(i) * MASSES_C(j) * (
                _PI_SYM(i, j) + _PI_SYM(3 + i, 3 + j)
                - 2 * _PI_SYMIK(i, 3 + j))
    return beta * total


def _PI_SYMIK(i, j):
    return _PI_SYM(i, j)


_MEAN = sp.Symbol("M_mix")


def _discharge(beta):
    weights = tuple((1 - Q) * m for m in MASSES) + tuple(Q * m for m in MASSES)
    total = sp.Integer(0)
    for i in range(6):
        for j in range(6):
            d = _MEAN * e_kernel(i, j, beta) \
                - (ATOMS[j] * sp.Symbol(f"h({ATOMS[i]})")
                   + ATOMS[i] * sp.Symbol(f"h({ATOMS[j]})")) / 2
            total += weights[i] * weights[j] * d
    return total


def mixture_mean_expr():
    weights = tuple((1 - Q) * m for m in MASSES) + tuple(Q * m for m in MASSES)
    return sum((weights[i] * ATOMS[i] for i in range(6)), sp.Integer(0))


def prove_identity() -> dict[str, Any]:
    beta = BETA
    gap = raw_gap_expr(beta)
    c_of_q = Q * (1 - Q) * remainder_expr(beta)
    num = _MEAN * (gap - c_of_q)
    poly = sp.expand(num - _discharge(beta))
    reduced = _reduce(poly.subs(_MEAN, mixture_mean_expr()))
    if reduced != 0:
        return {"status": "FAILED",
                "residual_terms": len(sp.Add.make_args(reduced)),
                "residual_str_head": str(reduced)[:240]}
    return {"status": "PROVED",
            "method": ("sympy multilinear expansion in the canonical free "
                       "log basis of the three-atom paired class"),
            "log_basis_size": (len(_LOG_ATOM) + len(_LOG_1M)
                               + len(_PROD_1M) + len(_PROTOCOL_1M))}


# ----- numeric certification -------------------------------------------------

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
    return mean * _e_mp(s, t, beta) - (t * _h_mp(s) + s * _h_mp(t)) / 2


def evaluate_three_atom(values, beta):
    """values = (a1, a3, q, x1, x2, x3, y1, y2, y3).  Returns
    (discharge, c_ch, mean)."""
    a1, a3, q, x1, x2, x3, y1, y2, y3 = values
    masses = (a1, a3, 1 - a1 - a3)
    sup = (x1, x2, x3, y1, y2, y3)
    weights = tuple((1 - q) * m for m in masses) + tuple(q * m for m in masses)
    mean = sum((weights[k] * sup[k] for k in range(6)), mpmath.mpf(0))
    discharge = mpmath.mpf(0)
    for i in range(6):
        for j in range(6):
            discharge += weights[i] * weights[j] * _d_mp(
                sup[i], sup[j], mean, beta)
    discharge = discharge / mean
    quad = mpmath.mpf(0)
    for i in range(3):
        for j in range(3):
            quad += masses[i] * masses[j] * (
                _h_mp(_pi_mp(sup[i], sup[j]))
                + _h_mp(_pi_mp(sup[3 + i], sup[3 + j]))
                - 2 * _h_mp(_pi_mp(sup[i], sup[3 + j])))
    c_ch = beta * quad
    return discharge, c_ch, mean


def random_three_atom(rng, m):
    def _fm(v):
        return mpmath.mpf(v.numerator) / mpmath.mpf(v.denominator)
    while True:
        a1 = _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
        a2 = _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
        if a1 + a2 >= 1:
            continue
        a3 = 1 - a1 - a2
        q = _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
        supports = tuple(
            _fm(Fraction(rng.random()).limit_denominator(10 ** 6))
            for _ in range(6))
        values = (a1, a2, q) + supports
        masses_sum = a1 + a2 + a3
        mean = q * (supports[0] * a1 + supports[1] * a2 + a3 * supports[2])             + (1 - q) * (a1 * supports[3] + a2 * supports[4] + a3 * supports[5])
        if mean >= m and min(supports) > mpmath.mpf("1e-3"):
            return values


def numeric_certification(dps=90):
    parameters = solve_equation_parameters(dps)
    beta = parameters.beta
    rng = random.Random(CERT_SEED)
    worst = mpmath.mpf(0)
    records = []
    with mpmath.workdps(dps):
        for trial in range(6):
            values = random_three_atom(rng, parameters.mean)
            a1, a3, q, x1, x2, x3, y1, y2, y3 = values
            a2 = 1 - a1 - a3
            vals9 = (a1, a3, q, x1, x2, x3, y1, y2, y3)
            ref = gap_mp(vals9, beta)
            discharge, c_ch, mean = evaluate_three_atom(values, beta)
            total = discharge + q * (1 - q) * c_ch
            diff = abs(total - ref)
            worst = max(worst, diff)
            records.append({
                "trial": trial,
                "values": [mpmath.nstr(v, 30) for v in values],
                "raw_gap": mpmath.nstr(ref, 40),
                "identity_absdiff": mpmath.nstr(diff, 8),
                "c_ch": mpmath.nstr(c_ch, 40),
                "discharge": mpmath.nstr(discharge, 40),
                "mean": mpmath.nstr(mean, 10),
            })
        if not worst < mpmath.mpf(10) ** -85:
            raise AssertionError("three-atom identity check failed")
    return {
        "claim_status": "MACHINE-VERIFIED",
        "dps": dps,
        "configurations": records,
        "worst_identity_absdiff": mpmath.nstr(worst, 8),
        "seed": CERT_SEED,
    }


def build_report() -> dict[str, Any]:
    parameters = solve_equation_parameters(100)
    symbolic = prove_identity()
    record: dict[str, Any] = {
        "tool": "liu9_paired_class_c.py",
        "constants": {
            "m": mpmath.nstr(parameters.mean, 40),
            "beta": mpmath.nstr(parameters.beta, 40),
        },
        "theorem": (
            "exact channel form c = beta*(E_K(P0) + E_K(P1) - 2 E_K(P0,P1)), "
            "K = h(pi(s,t)) for the three-atom paired class"
        ),
        "statements": {"three_atom_c_channel": symbolic},
        "seed": CERT_SEED,
        "limitations": [
            "r = 3 step checked symbolically only in this module; see "
            "liu9_scalar_margin_region2.py for the analogous two-atom case "
            "on two channels (liu9_scalar_margin_region2.py), already PROVED"],
    }
    if symbolic["status"] == "PROVED":
        numeric = numeric_certification()
        record["numeric_certification"] = numeric
        record["claim_status"] = "PROVED"
    elif symbolic["status"] == "FAILED":
        record["claim_status"] = "REFUTED"
    else:
        record["claim_status"] = "NUMERICAL"
    record["report_sha256"] = _digest(record)
    return record


def _digest(payload):
    return hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, indent=1, sort_keys=True, default=str) + "\n")
    print("GENERAL PAIRED-CLASS c CHANNEL, THREE ATOMS")
    print("symbolic:", report["statements"]["three_atom_c_channel"]["status"])
    print("claim:", report["claim_status"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0 if report["claim_status"] == "PROVED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
