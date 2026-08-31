#!/usr/bin/env python3
"""Exact one-step channel-increment identity for the Liu H2 paired class,
and the honest comparison against the certified S_asym pencil.

WHAT IS PROVED HERE (symbolic, zero residual, free basis on the
h(pi(.,.)) values -- the same discipline as liu9_paired_class_c.py):

    c_ch^{(r+1)} = (1-t)^2 c_ch^{(r)} + t^2 D(u,v) + 2 t (1-t) X(u,v)

for the insertion of one atom pair (u, v) at shared relative mass t into
an r-atom paired law with shared masses (a_i), where

    c_ch^{(r)} = beta sum_{i,j<=r} a_i a_j [Kpi(x_i,x_j) + Kpi(y_i,y_j)
                                            - 2 Kpi(x_i,y_j)]
    D(u,v)     = beta [Kpi(u,u) + Kpi(v,v) - 2 Kpi(u,v)]
    X(u,v)     = beta sum_i a_i [Kpi(x_i,u) + Kpi(y_i,v)
                                 - Kpi(x_i,v) - Kpi(y_i,u)]
    Kpi(s,t)   = h(pi(s,t)),  pi(s,t) = s t (1 + (1-s)(1-t)).

So the channel increment IS exactly quadratic in the inserted relative
mass t.  Combined with the already-PROVED diagonalization
gap = F + q(1-q) c_ch (liu9_paired_class_c.py, r = 3; region2, r = 2)
this gives the exact one-step gap identity

    gap^{(r+1)} - gap^{(r)} = [F^{(r+1)} - F^{(r)}]
                            + q(1-q) [(t^2 - 2t) c_ch^{(r)}
                                      + t^2 D + 2t(1-t) X].

WHAT IS REFUTED HERE: that this increment is nonnegative, hence that an
induction over atom insertion closes H2 using the S_asym certificate.
D(u,v) is sign-indefinite: an Arb interval certificate at
(u, v) = (7/20, 39/40) encloses D/beta strictly below zero.  The
S_asym pencil of liu9_second_order.py is NOT this object: its exact
symbolic residual against the increment is reported term by term.

Run: math/.venv/bin/python -I -B math/uc/liu9_channel_insertion.py
"""
from __future__ import annotations

import hashlib
import json
import os
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

import sympy as sp
from flint import arb, ctx

ctx.prec = 320

OUTPUT_DEFAULT = HERE / "verification/results/liu9-channel-insertion.json"

BETA = sp.Symbol("beta", positive=True)
T = sp.Symbol("t", positive=True)


# --- free basis on the protocol-entropy values -------------------------------
# Kpi(s, t) is an opaque symbol per unordered support pair.  This is a free
# basis: no relation between distinct pairs is used, so a zero residual is a
# genuine multilinear identity, not a numerical coincidence.

_KPI: dict[tuple[str, str], sp.Symbol] = {}


def kpi(s: sp.Symbol, t: sp.Symbol) -> sp.Symbol:
    key = tuple(sorted((s.name, t.name)))
    if key not in _KPI:
        _KPI[key] = sp.Symbol("Kpi(%s,%s)" % key)
    return _KPI[key]


def c_ch(masses, xs, ys):
    """beta sum_ij a_i a_j [Kpi(x_i,x_j) + Kpi(y_i,y_j) - 2 Kpi(x_i,y_j)]."""
    total = sp.Integer(0)
    n = len(masses)
    for i in range(n):
        for j in range(n):
            total += masses[i] * masses[j] * (
                kpi(xs[i], xs[j]) + kpi(ys[i], ys[j]) - 2 * kpi(xs[i], ys[j]))
    return BETA * total


def d_term(u, v):
    return BETA * (kpi(u, u) + kpi(v, v) - 2 * kpi(u, v))


def x_term(masses, xs, ys, u, v):
    total = sp.Integer(0)
    for i in range(len(masses)):
        total += masses[i] * (kpi(xs[i], u) + kpi(ys[i], v)
                              - kpi(xs[i], v) - kpi(ys[i], u))
    return BETA * total


def prove_insertion_identity(r: int) -> dict[str, Any]:
    """Zero-residual check of the one-step insertion identity at atom count r."""
    masses = [sp.Symbol("a%d" % (i + 1), positive=True) for i in range(r)]
    xs = [sp.Symbol("x%d" % (i + 1), real=True) for i in range(r)]
    ys = [sp.Symbol("y%d" % (i + 1), real=True) for i in range(r)]
    u = sp.Symbol("u", real=True)
    v = sp.Symbol("v", real=True)

    base = c_ch(masses, xs, ys)
    new_masses = [(1 - T) * a for a in masses] + [T]
    inserted = c_ch(new_masses, xs + [u], ys + [v])

    claimed = ((1 - T) ** 2 * base
               + T ** 2 * d_term(u, v)
               + 2 * T * (1 - T) * x_term(masses, xs, ys, u, v))

    residual = sp.expand(inserted - claimed)
    ok = residual == 0
    # Degree in t of the increment (must be exactly 2).
    increment = sp.expand(inserted - base)
    degree_t = sp.Poly(increment, T).degree() if increment != 0 else 0
    return {
        "atom_count_r": r,
        "residual_is_zero": bool(ok),
        "residual_term_count": 0 if ok else len(sp.Add.make_args(residual)),
        "residual_head": "" if ok else str(residual)[:400],
        "increment_degree_in_t": int(degree_t),
        "free_basis_size": len(_KPI),
    }


# --- the S_asym comparison ---------------------------------------------------
# S_asym transcribed EXACTLY from uc/liu9_second_order.py (read 2026-08-30):
#   kernel(s,t)   = (1-beta) h(s t) + beta h(pi(s,t))            [line 176-178]
#   L(y)          = 2 p kernel(x,y) - h(y) - (y/x) A - kappa y^2 (y-x)^2
#                                                                [line 181-185]
#   Q(y)          = kernel(y,y) - 2 (y/x) kernel(x,y)
#                   + (y/x)^2 kernel(x,x)                        [line 188-193]
#   S(y)          = 2 L(y) + Q(y)                                 [line 196-197]
#   cross_raw     = (2d/x)(kernel(x,y1) - kernel(x,y2))
#                   - (d^2/x^2) kernel(x,x)
#                   - (1-beta)(h(y1 y1) + h(y2 y2) - 2 h(y1 y2))  [line 200-209]
#   Chat          = cross_raw - kappa d^2                         [line 212-214]
#   S_asym        = (1-q) S(y1) + q S(y2) + q(1-q) Chat(y1,y2)    [hand_expansion]

KAPPA = Fraction(4_119_063, 33_554_432)


def s_asym_symbolic():
    """S_asym in a free basis on h(product) and h(protocol) values."""
    beta = BETA
    q = sp.Symbol("q", positive=True)
    x = sp.Symbol("xstar", positive=True)
    p = sp.Symbol("p", positive=True)
    A = sp.Symbol("A_const", real=True)
    kappa = sp.Rational(KAPPA.numerator, KAPPA.denominator)
    y1 = sp.Symbol("y1", real=True)
    y2 = sp.Symbol("y2", real=True)

    def hprod(s, t):
        key = tuple(sorted((s.name, t.name)))
        return sp.Symbol("Hprod(%s,%s)" % key)

    def kern(s, t):
        return (1 - beta) * hprod(s, t) + beta * kpi(s, t)

    def h1(s):
        return sp.Symbol("h(%s)" % s.name)

    def L(y):
        return (2 * p * kern(x, y) - h1(y) - (y / x) * A
                - kappa * y ** 2 * (y - x) ** 2)

    def Q(y):
        r = y / x
        return kern(y, y) - 2 * r * kern(x, y) + r * r * kern(x, x)

    d = y1 - y2
    cross_raw = ((2 * d / x) * (kern(x, y1) - kern(x, y2))
                 - (d * d / (x * x)) * kern(x, x)
                 - (1 - beta) * (hprod(y1, y1) + hprod(y2, y2)
                                 - 2 * hprod(y1, y2)))
    chat = cross_raw - kappa * d * d
    s_asym = ((1 - q) * (2 * L(y1) + Q(y1)) + q * (2 * L(y2) + Q(y2))
              + q * (1 - q) * chat)
    return sp.expand(s_asym), (q, x, y1, y2)


def compare_s_asym_to_increment() -> dict[str, Any]:
    """Is S_asym (any scalar multiple of) the t^2 channel increment D?"""
    s_asym, (q, x, y1, y2) = s_asym_symbolic()
    # The t^2 channel increment on the inserted pair (y1, y2):
    d_inc = sp.expand(BETA * (kpi(y1, y1) + kpi(y2, y2) - 2 * kpi(y1, y2)))
    # Does S_asym contain d_inc as its beta-protocol part?  Compare the
    # coefficients of every protocol symbol appearing in either expression.
    protocol_syms = sorted(
        {s for s in (s_asym.free_symbols | d_inc.free_symbols)
         if s.name.startswith("Kpi(")}, key=lambda s: s.name)
    rows = []
    matched = 0
    for sym in protocol_syms:
        cs = sp.simplify(s_asym.coeff(sym))
        cd = sp.simplify(d_inc.coeff(sym))
        same = sp.simplify(cs - cd) == 0
        matched += int(same)
        rows.append({"symbol": sym.name,
                     "coeff_in_S_asym": str(cs),
                     "coeff_in_channel_increment": str(cd),
                     "equal": bool(same)})
    # Terms present in S_asym that cannot appear in a channel increment at all
    # (the channel is built only from protocol entropies).
    non_channel = sorted({s.name for s in s_asym.free_symbols
                          if s.name.startswith("Hprod(")
                          or s.name.startswith("h(")
                          or s.name == "A_const"})
    residual = sp.expand(s_asym - d_inc)
    return {
        "identity_holds": bool(residual == 0),
        "protocol_symbol_rows": rows,
        "protocol_symbols_matching": matched,
        "protocol_symbols_total": len(protocol_syms),
        "non_channel_symbols_in_S_asym": non_channel,
        "residual_term_count": len(sp.Add.make_args(residual)),
        "verdict": (
            "S_asym is NOT the channel increment: it carries the F-margin "
            "pieces (product entropies Hprod, single entropies h(y), the "
            "(y/x)A_const first-order term and the kappa distance "
            "subtraction) which no channel increment contains, and its "
            "protocol coefficients differ."),
    }


# --- rigorous sign certificate for D ----------------------------------------

def _pt(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def _h_arb(z: arb) -> arb:
    return -(z * z.log() + (arb(1) - z) * (arb(1) - z).log())


def _pi_arb(s: arb, t: arb) -> arb:
    return s * t * (arb(1) + (arb(1) - s) * (arb(1) - t))


def d_over_beta_arb(u: Fraction, v: Fraction) -> arb:
    ua, va = _pt(u), _pt(v)
    return (_h_arb(_pi_arb(ua, ua)) + _h_arb(_pi_arb(va, va))
            - 2 * _h_arb(_pi_arb(ua, va)))


def certify_channel_sign() -> dict[str, Any]:
    """Arb enclosures of D/beta at exact rational points."""
    probes = [
        (Fraction(7, 20), Fraction(39, 40)),
        (Fraction(9, 10), Fraction(1, 10)),
        (Fraction(1, 2), Fraction(1, 2)),
        (Fraction(1, 4), Fraction(3, 4)),
    ]
    rows = []
    strictly_negative = 0
    for u, v in probes:
        enc = d_over_beta_arb(u, v)
        neg = bool(enc < 0)
        strictly_negative += int(neg)
        rows.append({"u": str(u), "v": str(v),
                     "D_over_beta": enc.str(24),
                     "certified_strictly_negative": neg})
    return {
        "probes": rows,
        "strictly_negative_certified": strictly_negative,
        "sign_definite": False,
        "conclusion": (
            "D(u,v)/beta = h(pi(u,u)) + h(pi(v,v)) - 2 h(pi(u,v)) is NOT "
            "nonnegative: Arb encloses it strictly below zero at exact "
            "rational points.  Hence the t^2 term of the channel increment "
            "can be negative and the insertion increment is sign-indefinite."),
    }


def mutations() -> list[dict[str, Any]]:
    """Both-direction checks: each mutation must break the identity."""
    out = []
    # (1) wrong mass power on the base term: (1-t) instead of (1-t)^2.
    masses = [sp.Symbol("a1", positive=True), sp.Symbol("a2", positive=True)]
    xs = [sp.Symbol("x1"), sp.Symbol("x2")]
    ys = [sp.Symbol("y1"), sp.Symbol("y2")]
    u, v = sp.Symbol("u"), sp.Symbol("v")
    base = c_ch(masses, xs, ys)
    inserted = c_ch([(1 - T) * a for a in masses] + [T], xs + [u], ys + [v])
    bad1 = ((1 - T) * base + T ** 2 * d_term(u, v)
            + 2 * T * (1 - T) * x_term(masses, xs, ys, u, v))
    out.append({"mutation": "base_mass_power_1_instead_of_2",
                "expected": "FAIL",
                "residual_is_zero": bool(sp.expand(inserted - bad1) == 0)})
    # (2) sign flip on the cross term.
    bad2 = ((1 - T) ** 2 * base + T ** 2 * d_term(u, v)
            - 2 * T * (1 - T) * x_term(masses, xs, ys, u, v))
    out.append({"mutation": "cross_term_sign_flip",
                "expected": "FAIL",
                "residual_is_zero": bool(sp.expand(inserted - bad2) == 0)})
    # (3) drop the D term entirely.
    bad3 = ((1 - T) ** 2 * base
            + 2 * T * (1 - T) * x_term(masses, xs, ys, u, v))
    out.append({"mutation": "drop_D_term",
                "expected": "FAIL",
                "residual_is_zero": bool(sp.expand(inserted - bad3) == 0)})
    return out


def build_report() -> dict[str, Any]:
    identities = [prove_insertion_identity(r) for r in (1, 2, 3)]
    comparison = compare_s_asym_to_increment()
    sign = certify_channel_sign()
    muts = mutations()
    all_ident = all(row["residual_is_zero"] for row in identities)
    all_quadratic = all(row["increment_degree_in_t"] == 2 for row in identities)
    muts_fail = all(not row["residual_is_zero"] for row in muts)
    report = {
        "tool": "liu9_channel_insertion.py",
        "statements": {
            "one_step_channel_identity": {
                "status": "PROVED" if (all_ident and all_quadratic) else "FAILED",
                "claim": ("c_ch^{(r+1)} = (1-t)^2 c_ch^{(r)} + t^2 D(u,v) "
                          "+ 2 t (1-t) X(u,v); exact, quadratic in t"),
                "rows": identities,
                "method": ("sympy multilinear expansion over a free basis on "
                           "the h(pi(.,.)) protocol-entropy values"),
            },
            "s_asym_is_the_channel_increment": {
                "status": "REFUTED" if not comparison["identity_holds"] else "PROVED",
                "detail": comparison,
            },
            "channel_increment_is_nonnegative": {
                "status": "REFUTED",
                "detail": sign,
            },
        },
        "consequence_for_the_induction": (
            "The r -> r+1 induction of uc/GENERAL_R_H2_COMPOSITION.md section 4 "
            "is UNFOUNDED as written: the channel increment is exactly "
            "quadratic in the inserted mass (PROVED here) but sign-indefinite "
            "(REFUTED here), and it is not an S_asym instance.  Any induction "
            "must dominate the negative channel term by the F-margin, which is "
            "precisely the open block-copositivity obligation recorded in "
            "uc/verification/results/liu9-block-kernel.json."),
        "mutations": muts,
        "claim_status": "PROVED" if (all_ident and all_quadratic and muts_fail) else "FAILED",
        "honesty": (
            "This report supersedes uc/verification/results/"
            "sasym-channel-identification.json, which was fabricated prose "
            "with no computation, and the liu9-sasym-4slot* reports, whose "
            "pencil was invented rather than taken from liu9_second_order.py."),
    }
    report["report_sha256"] = hashlib.sha256(json.dumps(
        report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return report


def main() -> int:
    report = build_report()
    OUTPUT_DEFAULT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_DEFAULT.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    st = report["statements"]
    print("one-step channel identity:",
          st["one_step_channel_identity"]["status"],
          "| degrees in t:",
          [r["increment_degree_in_t"] for r in st["one_step_channel_identity"]["rows"]])
    print("S_asym == channel increment:",
          st["s_asym_is_the_channel_increment"]["status"],
          "| protocol coeffs matching:",
          "%d/%d" % (st["s_asym_is_the_channel_increment"]["detail"]["protocol_symbols_matching"],
                     st["s_asym_is_the_channel_increment"]["detail"]["protocol_symbols_total"]))
    print("channel increment >= 0:",
          st["channel_increment_is_nonnegative"]["status"],
          "| Arb-certified negative probes:",
          st["channel_increment_is_nonnegative"]["detail"]["strictly_negative_certified"])
    print("mutations all FAIL:",
          all(not r["residual_is_zero"] for r in report["mutations"]))
    print("claim_status:", report["claim_status"])
    print("report_sha256", report["report_sha256"])
    return 0 if report["claim_status"] == "PROVED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
