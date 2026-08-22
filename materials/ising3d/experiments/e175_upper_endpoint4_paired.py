"""Exact minimal paired-momentum/four-point SDP for the upper endpoint.

The finite certificate class keeps the aggregate nonzero Fourier power

    R = N^-2 sum_{k != 0} |sigmahat(k)|^2 = 1 - (M/N)^2

and its genuinely four-point moment E[R^2].  A positive exact SDP order floor
would imply a strict upper-endpoint improvement.  The exact optimum is zero.
This module is a producer component; e177 writes the combined JSON artifact.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
UPPER_INFRARED = ROOT / "results" / "bounds" / "upper_infrared.json"
SCRIPT = "experiments/e175_upper_endpoint4_paired.py"
V_TARGET = Fraction(6, 25)
ATANH_TERMS = 20
ETA_CHALLENGE = Fraction(1, 1024)


def fstr(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_decimal(value: str) -> Fraction:
    return Fraction(value)


def decimal_floor(value: Fraction, digits: int = 24) -> str:
    """Return a directed-down decimal for a nonnegative Fraction."""
    if value < 0:
        return "-" + decimal_ceil(-value, digits)
    scale = 10**digits
    integer = value.numerator * scale // value.denominator
    whole, frac = divmod(integer, scale)
    return f"{whole}.{frac:0{digits}d}"


def decimal_ceil(value: Fraction, digits: int = 24) -> str:
    """Return a directed-up decimal for a nonnegative Fraction."""
    if value < 0:
        return "-" + decimal_floor(-value, digits)
    scale = 10**digits
    integer = -((-value.numerator * scale) // value.denominator)
    whole, frac = divmod(integer, scale)
    return f"{whole}.{frac:0{digits}d}"


def atanh_interval(v: Fraction, terms: int = ATANH_TERMS) -> tuple[Fraction, Fraction]:
    """Positive-series enclosure of atanh(v), 0 < v < 1.

    The lower endpoint contains terms n=0,...,terms-1.  In the positive tail,
    1/(2n+1) <= 1/(2*terms+1), giving an exact geometric upper remainder.
    """
    if not Fraction(0) < v < Fraction(1):
        raise ValueError("atanh_interval requires 0 < v < 1")
    low = sum((v ** (2 * n + 1) / (2 * n + 1) for n in range(terms)), Fraction(0))
    first_power = v ** (2 * terms + 1)
    tail = first_power / ((2 * terms + 1) * (1 - v * v))
    return low, low + tail


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(checks: list[dict[str, Any]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def _psd_2x2(a: Fraction, b: Fraction, c: Fraction) -> bool:
    return a >= 0 and c >= 0 and a * c - b * b >= 0


def build_paired_data() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Build the exact paired-momentum certificate and all rational checks."""
    checks: list[dict[str, Any]] = []
    source = json.loads(UPPER_INFRARED.read_text(encoding="utf-8"))
    i3_text = source["data"]["certified_constants"]["I3"]
    i3_lo, i3_hi = map(parse_decimal, i3_text)
    k_lo, k_hi = atanh_interval(V_TARGET)

    _record(checks, "target_atanh_interval_ordered", 0 < k_lo < k_hi, f"K in [{fstr(k_lo)}, {fstr(k_hi)}]")
    _record(
        checks,
        "target_strictly_below_incumbent",
        2 * k_hi < i3_lo,
        f"I3_lo-2*K_hi={fstr(i3_lo - 2 * k_hi)} > 0",
    )

    incumbent_gap_lo = i3_lo / 2 - k_hi
    incumbent_gap_hi = i3_hi / 2 - k_lo
    delta_lo = i3_lo / (2 * k_hi) - 1
    delta_hi = i3_hi / (2 * k_lo) - 1
    _record(
        checks,
        "positive_endpoint_gap",
        0 < incumbent_gap_lo <= incumbent_gap_hi,
        f"I3/2-K >= {fstr(incumbent_gap_lo)}",
    )
    _record(
        checks,
        "summed_ir_cap_redundant",
        0 < delta_lo <= delta_hi,
        f"I3/(2*K)-1 in [{fstr(delta_lo)}, {fstr(delta_hi)}]",
    )

    # Primal SDP: maximize s subject to [[1,s],[s,q]] >= 0 and 0 <= q <= s <= 1.
    # Exact lower witness (s,q)=(1,1); the scalar row 1-s>=0 is the exact dual upper bound.
    s_witness = Fraction(1)
    q_witness = Fraction(1)
    primal_rows = {
        "q_nonnegative": q_witness,
        "s_minus_q": s_witness - q_witness,
        "one_minus_s": 1 - s_witness,
        "moment_determinant": q_witness - s_witness * s_witness,
    }
    primal_feasible = (
        _psd_2x2(Fraction(1), s_witness, q_witness)
        and q_witness >= 0
        and s_witness - q_witness >= 0
        and 1 - s_witness >= 0
    )
    _record(checks, "paired_primal_witness", primal_feasible, f"(s,q)=({s_witness},{q_witness}), rows={primal_rows}")

    optimum_s = Fraction(1)
    optimum_eta = Fraction(0)
    dual_upper = Fraction(1)  # the exact identity 1-s = 1*(1-s)
    # Exact rational dual-certificate feasibility system.  For a requested
    # floor eta, set U=1-eta and seek a,b,c>=0 and Z=[[z00,z01],[z01,z11]]>=0
    # such that
    #   U-s = a*q + b*(s-q) + c*(1-s) + <Z, [[1,s],[s,q]]>.
    # Coefficient matching gives the three affine equalities below.
    dual_eta = Fraction(0)
    dual_a = dual_b = Fraction(0)
    dual_c = Fraction(1)
    dual_z00 = dual_z01 = dual_z11 = Fraction(0)
    dual_coefficients_match = (
        dual_c + dual_z00 == 1 - dual_eta
        and dual_b - dual_c + 2 * dual_z01 == -1
        and dual_a - dual_b + dual_z11 == 0
        and _psd_2x2(dual_z00, dual_z01, dual_z11)
    )
    _record(
        checks,
        "paired_dual_feasibility_identity",
        dual_coefficients_match,
        "eta=0, a=b=0, c=1, Z=0 matches 1-s exactly",
    )
    _record(
        checks,
        "paired_exact_primal_dual_optimum",
        primal_feasible and optimum_s == dual_upper and optimum_eta == 1 - optimum_s,
        "primal s=1 and dual row s<=1 agree exactly",
    )

    challenge_upper = 1 - ETA_CHALLENGE
    challenge_gap = s_witness - challenge_upper
    _record(
        checks,
        "paired_positive_floor_challenge_infeasible",
        challenge_gap == ETA_CHALLENGE > 0,
        f"witness violates s<={fstr(challenge_upper)} by {fstr(challenge_gap)}",
    )

    # A concrete bounded-spin lift: any zero-magnetisation configuration has R=1.
    # For N=8 the checkerboard has zero total spin, so its R and R^2 are both one.
    sites = [(x, y, z) for x in range(2) for y in range(2) for z in range(2)]
    spins = [(-1) ** (x + y + z) for x, y, z in sites]
    magnetisation = sum(spins)
    fourier_power = []
    for jx in range(2):
        for jy in range(2):
            for jz in range(2):
                amplitude = sum(
                    sigma * ((-1) ** (jx * x + jy * y + jz * z))
                    for sigma, (x, y, z) in zip(spins, sites)
                )
                fourier_power.append(amplitude * amplitude)
    n_sites = len(sites)
    r_checkerboard = Fraction(sum(fourier_power[1:]), n_sites * n_sites)
    # The zero mode is index 0 under the loop ordering above.
    _record(
        checks,
        "paired_bounded_spin_lift",
        magnetisation == 0 and fourier_power[0] == 0 and sum(fourier_power) == n_sites**2 and r_checkerboard == 1,
        f"M={magnetisation}, Parseval={sum(fourier_power)}, R={r_checkerboard}",
    )

    data = {
        "class_name": "PM4-aggregate-one-moment",
        "scope": (
            "thermodynamic-limit two-point summed infrared cap plus the exact bounded-spin paired-momentum moments "
            "s=E[R], q=E[R^2] for R=N^-2 sum_{k!=0}|sigmahat(k)|^2; no mode-resolved "
            "four-point, DLR, or random-current constraints"
        ),
        "whole_incumbent_interval_obstruction": {
            "reason": (
                "for every K<=I3/2 the limiting summed infrared cap s<=I3/(2K) has right side >=1, "
                "so it is redundant against s<=1 and the same abstract aggregate (s,q)=(1,1) witness remains feasible"
            ),
            "exact_optimum_s": "1/1",
            "exact_optimum_order_floor_eta": "0/1",
        },
        "target": {
            "v": fstr(V_TARGET),
            "K_definition": "atanh(6/25)",
            "K_interval": [fstr(k_lo), fstr(k_hi)],
            "K_decimal_interval": [decimal_floor(k_lo), decimal_ceil(k_hi)],
            "I3_interval_input": list(i3_text),
            "incumbent_gap_interval": [fstr(incumbent_gap_lo), fstr(incumbent_gap_hi)],
            "incumbent_gap_decimal_lower": decimal_floor(incumbent_gap_lo),
            "delta_star_interval": [fstr(delta_lo), fstr(delta_hi)],
            "input_sha256": {str(UPPER_INFRARED.relative_to(ROOT)): _sha256(UPPER_INFRARED)},
        },
        "primal_sdp": {
            "variables": ["s=E[R]", "q=E[R^2]"],
            "objective": "maximize s (equivalently maximize the obstruction to M_L^2=1-s)",
            "constraints": [
                "[[1,s],[s,q]] is positive semidefinite",
                "q>=0",
                "s-q>=0 (R^2<=R because 0<=R<=1)",
                "1-s>=0",
                "limiting s<=I3/(2*K_target), redundant because certified I3/(2*K_target)>1",
            ],
            "exact_optimum_s": fstr(optimum_s),
            "exact_optimum_order_floor_eta": fstr(optimum_eta),
            "primal_witness": {"s": fstr(s_witness), "q": fstr(q_witness), "moment_matrix": [["1", "1"], ["1", "1"]]},
            "dual_certificate": {"upper_bound": "1", "identity": "1-s = 1*(1-s)"},
        },
        "dual_feasibility_sdp": {
            "objective": "maximize eta",
            "variables": ["eta", "a>=0", "b>=0", "c>=0", "Z=[[z00,z01],[z01,z11]] positive semidefinite"],
            "certificate_identity": (
                "1-eta-s = a*q + b*(s-q) + c*(1-s) "
                "+ <Z,[[1,s],[s,q]]>"
            ),
            "coefficient_equalities": [
                "c+z00=1-eta",
                "b-c+2*z01=-1",
                "a-b+z11=0",
            ],
            "exact_optimal_point": {
                "eta": fstr(dual_eta),
                "a": fstr(dual_a),
                "b": fstr(dual_b),
                "c": fstr(dual_c),
                "Z": [["0", "0"], ["0", "0"]],
            },
            "rational_certificate_rule": (
                "any rational feasible point with eta>0 proves s<=1-eta "
                "by a sum of nonnegative primal rows and a PSD inner product"
            ),
        },
        "certificate_feasibility": {
            "statement": (
                "a rational SDP dual certificate of s<=1-eta with eta>0 applies to every even torus, "
                "hence gives M_L^2>=eta uniformly and certifies K_c<=atanh(6/25)<I3/2"
            ),
            "optimal_eta": fstr(optimum_eta),
            "challenge_eta": fstr(ETA_CHALLENGE),
            "challenge_upper_s": fstr(challenge_upper),
            "exact_infeasibility_gap": fstr(challenge_gap),
        },
        "bounded_spin_obstruction": {
            "configuration": "2x2x2 checkerboard sigma(x)=(-1)^(x1+x2+x3)",
            "magnetisation": str(magnetisation),
            "fourier_powers": [str(value) for value in fourier_power],
            "R": fstr(r_checkerboard),
            "R_squared": fstr(r_checkerboard * r_checkerboard),
        },
    }
    return data, checks


def main() -> int:
    _, checks = build_paired_data()
    for check in checks:
        print(f"{'PASS' if check['passed'] else 'FAIL'} {check['name']} -- {check['detail']}")
    passed = all(check["passed"] for check in checks)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
