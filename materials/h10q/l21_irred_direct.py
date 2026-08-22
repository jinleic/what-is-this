#!/usr/bin/env python3
"""Direct reciprocal and p-adic attack on the constructed-branch octic P.

The script proves symbolic identities, checks their exact specializations on the
353 canonical (cell, a) pairs, and records the remaining gap without promoting
finite modular certificates to a uniform theorem.
"""
from __future__ import annotations

import json
import math
import sys
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import h10q  # noqa: E402  (the stipulated arithmetic authority)


OUT = HERE / "data" / "l21_irred_direct.jsonl"
REPORT = Path("/tmp/l21_irred_direct.md")
A_POOL = (1, 3, 5, 7)
IRRED_PRIME_LIMIT = 1000
PACE_EVERY = 200
PACE_SECONDS = 0.005


class Pacer:
    """Apply the requested low-duty-cycle pause to hot certificate loops."""

    def __init__(self) -> None:
        self.candidates = 0
        self.sleeps = 0

    def tick(self) -> None:
        self.candidates += 1
        if self.candidates % PACE_EVERY == 0:
            time.sleep(PACE_SECONDS)
            self.sleeps += 1


PACER = Pacer()


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def quadratic_character(value: int, prime: int) -> int:
    residue = value % prime
    if residue == 0:
        return 0
    return 1 if pow(residue, (prime - 1) // 2, prime) == 1 else -1


def choose_a(prime: int) -> int:
    """Replay L19/L20's first admissible member of the frozen odd a-pool."""
    for a in A_POOL:
        if quadratic_character(1 + 4 * a * a, prime) == -1:
            return a
    raise AssertionError((prime, "the frozen a-pool has no nonresidue A"))


def evaluate(coefficients: list[F], value: F | int) -> F:
    answer = F(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def poly_mul(left: list[F], right: list[F]) -> list[F]:
    product = [F(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            product[i + j] += left_coefficient * right_coefficient
    return product


def trace_polynomial(P: list[F]) -> list[F]:
    """Return T with P(b)=b^4 T(b+b^-1), for a reciprocal octic."""
    assert len(P) == 9 and all(P[index] == P[8 - index] for index in range(9))
    trace = [
        P[4] + 2 * P[0] - 2 * P[2],
        P[3] - 3 * P[1],
        P[2] - 4 * P[0],
        P[1],
        P[0],
    ]
    lifted = [
        trace[4],
        trace[3],
        trace[2] + 4 * trace[4],
        trace[1] + 3 * trace[3],
        trace[0] + 2 * trace[2] + 6 * trace[4],
        trace[1] + 3 * trace[3],
        trace[2] + 4 * trace[4],
        trace[3],
        trace[4],
    ]
    assert lifted == P
    return trace


def primitive_integer_poly(coefficients: list[F]) -> list[int]:
    denominator = 1
    for coefficient in coefficients:
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in coefficients]
    content = 0
    for coefficient in integers:
        content = math.gcd(content, abs(coefficient))
    assert content > 0
    primitive = [coefficient // content for coefficient in integers]
    if primitive[-1] < 0:
        primitive = [-coefficient for coefficient in primitive]
    return primitive


def exact_irred_prime(P: list[F]) -> int:
    primitive = primitive_integer_poly(P)
    assert len(primitive) == 9 and primitive[-1] != 0
    for prime in h10q.primerange(2, IRRED_PRIME_LIMIT + 1):
        PACER.tick()
        if h10q._l13_irred8(primitive, prime):
            return prime
    raise AssertionError("no exact degree-8 certificate below 1000")


def mod_fraction(value: F, prime: int) -> int:
    assert value.denominator % prime != 0
    return value.numerator * pow(value.denominator, prime - 2, prime) % prime


def lower_newton_polygon(P: list[F], prime: int) -> tuple[list[tuple[int, int]], list[tuple[int, F]]]:
    points = [(degree, h10q.vp(coefficient, prime))
              for degree, coefficient in enumerate(P) if coefficient]
    hull: list[tuple[int, int]] = []
    for point in points:
        while len(hull) >= 2:
            x1, y1 = hull[-2]
            x2, y2 = hull[-1]
            x3, y3 = point
            old_slope = F(y2 - y1, x2 - x1)
            new_slope = F(y3 - y2, x3 - x2)
            if old_slope >= new_slope:
                hull.pop()
            else:
                break
        hull.append(point)
    segments = [
        (x2 - x1, F(y2 - y1, x2 - x1))
        for (x1, y1), (x2, y2) in zip(hull, hull[1:])
    ]
    return hull, segments


def segment_records(segments: list[tuple[int, F]]) -> list[dict[str, Any]]:
    return [
        {
            "length": length,
            "slope": frac_text(slope),
            "slope_denominator": slope.denominator,
        }
        for length, slope in segments
    ]


def segment_key(segments: list[tuple[int, F]]) -> str:
    return ", ".join(f"{length}@{frac_text(slope)}" for length, slope in segments)


def fp2_mul(left: tuple[int, int], right: tuple[int, int], prime: int) -> tuple[int, int]:
    """Multiply x+y*t pairs in F_p[t]/(t^2+5)."""
    x1, y1 = left
    x2, y2 = right
    return ((x1 * x2 - 5 * y1 * y2) % prime,
            (x1 * y2 + y1 * x2) % prime)


def fp2_pow(value: tuple[int, int], exponent: int, prime: int) -> tuple[int, int]:
    answer = (1, 0)
    while exponent:
        if exponent & 1:
            answer = fp2_mul(answer, value, prime)
        value = fp2_mul(value, value, prime)
        exponent >>= 1
    return answer


def trace_residue_obstruction(prime: int) -> dict[str, Any]:
    """Mechanically replay the Q(sqrt(-5)) residue proof for a=1."""
    assert quadratic_character(5, prime) == -1 and prime != 5
    minus_five_character = quadratic_character(-5, prime)
    if minus_five_character == -1:
        plus = fp2_pow((0, 2), (prime * prime - 1) // 2, prime)
        minus = fp2_pow((0, -2 % prime), (prime * prime - 1) // 2, prime)
        assert plus == minus == (prime - 1, 0)
        norm_character = quadratic_character(20, prime)
        assert norm_character == -1
        return {
            "splitting_in_Q_sqrt_minus_5": "inert",
            "minus_5_character": -1,
            "characters_of_plus_minus_2t": [-1, -1],
            "norm_20_character": norm_character,
            "both_quadratic_factors_have_a_local_nonsquare_discriminant": True,
        }

    assert minus_five_character == 1
    root = next(residue for residue in range(1, prime)
                if residue * residue % prime == -5 % prime)
    at_root = [quadratic_character(2 * root, prime),
               quadratic_character(-2 * root, prime)]
    at_conjugate = [at_root[1], at_root[0]]
    assert sorted(at_root) == [-1, 1]
    assert sorted(at_conjugate) == [-1, 1]
    assert all(-1 in pair for pair in zip(at_root, at_conjugate))
    return {
        "splitting_in_Q_sqrt_minus_5": "split",
        "minus_5_character": 1,
        "chosen_t_residue": root,
        "characters_of_plus_minus_2t_at_prime": at_root,
        "characters_at_conjugate_prime": at_conjugate,
        "each_quadratic_factor_has_a_local_nonsquare_discriminant": True,
    }


def is_rational_square(value: F) -> bool:
    if value < 0:
        return False
    numerator_root = math.isqrt(value.numerator)
    denominator_root = math.isqrt(value.denominator)
    return (numerator_root * numerator_root == value.numerator
            and denominator_root * denominator_root == value.denominator)


def pair_record(w: int, unit: tuple[int, int]) -> dict[str, Any]:
    a_int = choose_a(w)
    a = F(a_int)
    z = F(w) * F(*unit)
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    delta = 1 - A * tau * tau
    Z = z**3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    assert a and z and D and A.denominator == 1
    assert delta == -4 * a**4 / A < 0
    assert h10q.vp(z, w) >= 1
    assert quadratic_character(A.numerator, w) == -1

    P = h10q._l10_P(a, Z, D, A, delta, s)
    assert len(P) == 9 and P[-1] != 0
    endpoint = 4 * a**8 * A * Z**4
    assert P[0] == P[8] == endpoint

    perturbation = 32 * A**3 * s**2 * D**2
    core = list(P)
    core[5] += perturbation
    assert all(core[index] == core[8 - index] for index in range(9))
    reciprocal_defect = [P[index] - P[8 - index] for index in range(9)]
    expected_defect = [F(0)] * 9
    expected_defect[3] = perturbation
    expected_defect[5] = -perturbation
    assert reciprocal_defect == expected_defect
    reciprocal = perturbation == 0
    assert reciprocal == (a_int == 1)
    # Equal endpoints force a rational reciprocal twist to have scale 1 and
    # twist k=+/-1; the nonzero b^1 coefficient then forces k=1.
    assert P[1] != 0
    rational_twisted_reciprocal = reciprocal

    mod_coefficients = [mod_fraction(coefficient, w) for coefficient in P]
    expected_mod = [0] * 9
    expected_mod[4] = mod_fraction(16 * D**2 * A**2, w)
    expected_mod[5] = mod_fraction(-32 * A**3 * s**2 * D**2, w)
    assert mod_coefficients == expected_mod
    assert expected_mod[4] != 0
    if expected_mod[5] == 0:
        factorization_type = "1^4"
        factorization_factors = [{"degree": 1, "multiplicity": 4}]
    else:
        factorization_type = "1^4 1"
        factorization_factors = [
            {"degree": 1, "multiplicity": 4},
            {"degree": 1, "multiplicity": 1},
        ]
    mod_degree = max(index for index, coefficient in enumerate(mod_coefficients)
                     if coefficient)
    assert mod_degree in (4, 5)
    primitive = primitive_integer_poly(P)
    assert not h10q._l13_irred8(primitive, w)

    w_vertices, w_segments = lower_newton_polygon(P, w)
    e = h10q.vp(z, w)
    if s == 0:
        expected_segments = [(4, F(-3 * e)), (4, F(3 * e))]
    else:
        assert h10q.vp(s, w) == 0  # true for the frozen bounded a-pool
        expected_segments = [
            (4, F(-3 * e)),
            (1, F(0)),
            (3, F(4 * e)),
        ]
    assert w_segments == expected_segments
    assert not (len(w_segments) == 1 and w_segments[0][1].denominator == 8)

    A_newton: list[dict[str, Any]] = []
    for prime, exponent in sorted(h10q.factorint(A.numerator).items()):
        vertices, segments = lower_newton_polygon(P, prime)
        A_newton.append(
            {
                "prime": prime,
                "v_p_A": exponent,
                "v_p_z": h10q.vp(z, prime),
                "v_p_D": h10q.vp(D, prime),
                "vertices": [[degree, valuation] for degree, valuation in vertices],
                "segments": segment_records(segments),
                "single_slope_denominator_8": (
                    len(segments) == 1 and segments[0][1].denominator == 8
                ),
                "any_slope_denominator_8": any(
                    slope.denominator == 8 for _length, slope in segments
                ),
            }
        )
    assert A_newton and not any(
        item["single_slope_denominator_8"] or item["any_slope_denominator_8"]
        for item in A_newton
    )

    trace_data: dict[str, Any] | None = None
    if reciprocal:
        assert a_int == 1 and A == 5 and quadratic_character(5, w) == -1
        trace = trace_polynomial(P)
        H = [F(-4), F(20), F(-5)]  # 16 - 5*(u-2)^2
        expected_trace = poly_mul(H, H)
        expected_trace = [F(4, 5) * Z**4 * coefficient
                          for coefficient in expected_trace]
        expected_trace[0] += 400 * D**2
        assert trace == expected_trace

        residue_obstruction = trace_residue_obstruction(w)
        factor_64 = 125 * D**2 + 64 * Z**4
        factor_1024 = 125 * D**2 + 1024 * Z**4
        norm_obstruction = factor_64 * factor_1024
        trace_norm = evaluate(trace, 2) * evaluate(trace, -2) / trace[-1]**2
        expected_norm = (F(4, 25) / Z**4) ** 2 * norm_obstruction
        assert trace_norm == expected_norm
        assert not is_rational_square(norm_obstruction)
        assert not is_rational_square(trace_norm)
        trace_data = {
            "trace_polynomial_coefficients_ascending": [
                frac_text(coefficient) for coefficient in trace
            ],
            "trace_formula": "400*D^2+(4/5)*Z^4*(16-5*(u-2)^2)^2",
            "trace_quartic_irreducible_over_Q": True,
            "trace_irreducibility_lemma": "PROVED_UNIFORMLY_FOR_v_w(z)>0_AND_(5|w)=-1",
            "Q_sqrt_minus_5_residue_obstruction": residue_obstruction,
            "lift_norm_square_class": frac_text(norm_obstruction),
            "lift_norm_is_rational_square": False,
            "degree_8_irreducibility_by_trace_norm": True,
            "degree_8_scope": "PROVED for this exact row by the uniform trace lemma plus the checked norm condition",
        }

    certificate_prime = exact_irred_prime(P)
    assert h10q._l13_irred8(primitive, certificate_prime)

    return {
        "type": "constructed-grid-pair",
        "label": "PROVED_EXACT_PER_ROW",
        "uniform_scope": "mechanical evidence unless a symbolic lemma is named",
        "cell": [w, list(unit)],
        "z": frac_text(z),
        "a": a_int,
        "A": frac_text(A),
        "tau_dagger": frac_text(tau),
        "delta_dagger": frac_text(delta),
        "A_character_mod_w": -1,
        "endpoint_common": frac_text(endpoint),
        "endpoint_formula_verified": True,
        "reciprocal": reciprocal,
        "palindromic_core": True,
        "reciprocal_defect": "32*A^3*s^2*D^2*(b^3-b^5)",
        "rational_twisted_reciprocal": rational_twisted_reciprocal,
        "mod_w": {
            "degree_after_reduction": mod_degree,
            "nonzero_coefficients_ascending": {
                str(index): coefficient
                for index, coefficient in enumerate(mod_coefficients) if coefficient
            },
            "identity": "16*A^2*b^4*(1-2*A*s^2*b)",
            "factorization_type": factorization_type,
            "factors": factorization_factors,
            "irreducible_degree_8": False,
        },
        "newton_w": {
            "v_w_z": e,
            "vertices": [[degree, valuation] for degree, valuation in w_vertices],
            "segments": segment_records(w_segments),
            "single_slope_denominator_8": False,
            "local_irreducibility_route": "REFUTED",
        },
        "newton_primes_dividing_A": A_newton,
        "trace_route": trace_data,
        "auxiliary_exact_irreducibility_certificate_prime": certificate_prime,
        "tau_dagger_reducible_instance": False,
    }


def count_table(rows: list[dict[str, Any]]) -> dict[str, Any]:
    a_counts = Counter(str(row["a"]) for row in rows)
    mod_counts = Counter(row["mod_w"]["factorization_type"] for row in rows)
    w_np_counts = Counter(
        ", ".join(
            f'{segment["length"]}@{segment["slope"]}'
            for segment in row["newton_w"]["segments"]
        )
        for row in rows
    )
    A_np_counts: Counter[str] = Counter()
    for row in rows:
        for item in row["newton_primes_dividing_A"]:
            slopes = ", ".join(
                f'{segment["length"]}@{segment["slope"]}'
                for segment in item["segments"]
            )
            A_np_counts[f'a={row["a"]}, p={item["prime"]}: {slopes}'] += 1
    trace_rows = [row for row in rows if row["trace_route"] is not None]
    trace_local_counts = Counter(
        row["trace_route"]["Q_sqrt_minus_5_residue_obstruction"]
        ["splitting_in_Q_sqrt_minus_5"]
        for row in trace_rows
    )
    certificate_counts = Counter(
        str(row["auxiliary_exact_irreducibility_certificate_prime"])
        for row in rows
    )
    return {
        "a_counts": dict(sorted(a_counts.items(), key=lambda item: int(item[0]))),
        "mod_w_factorization_counts": dict(sorted(mod_counts.items())),
        "newton_w_segment_counts": dict(sorted(w_np_counts.items())),
        "newton_A_segment_counts": dict(sorted(A_np_counts.items())),
        "trace_local_case_counts": dict(sorted(trace_local_counts.items())),
        "auxiliary_certificate_prime_counts": dict(
            sorted(certificate_counts.items(), key=lambda item: int(item[0]))
        ),
    }


def write_report(summary: dict[str, Any]) -> None:
    tables = summary["mechanical_tables"]
    lines = [
        "# L21 direct irreducibility attack",
        "",
        "## Verdict",
        "",
        "- **Uniform irreducibility on every admissible constructed parameter:** **OPEN**.",
        "- **Reciprocal structure:** **PROVED exactly.** On the constructed branch, `P` is a reciprocal octic exactly when `a=1`; otherwise it is a reciprocal core with one `b^5` perturbation, and no rational reciprocal twist absorbs it.",
        "- **The proposed place `w`: REFUTED as an irreducibility place.** Reduction modulo `w` is always `16 A^2 b^4(1-2 A s^2 b)`, and the `w`-Newton polygon always has at least two slopes, never one slope of denominator 8.",
        "- **Reciprocal subfamily (`a=1`, `(5|w)=-1`):** the trace quartic is **PROVED uniformly irreducible**. The octic is **PROVED under the explicit norm condition** `Xi(Z)` nonsquare. That condition holds exactly on all 189 canonical `a=1` rows.",
        f"- **Constructed-branch scan:** {summary['pairs']} exact `(cell,a)` pairs, {summary['tau_dagger_exact_irred_certificates']} exact auxiliary mod-prime certificates, and **0 reducible instances** on `tau_dagger`. This is EVIDENCE, not a uniform theorem.",
        "",
        "The lead L21 square-delta reducible-locus lemma is separate. Here `delta_dagger=-4a^4/A<0`, so this report stays entirely on the branch used by the construction.",
        "",
        "## Lemma L21-D1: endpoints and the unique reciprocity defect (PROVED)",
        "",
        "On `tau_dagger=(1+2a^2)/A`,",
        "",
        "```text",
        "delta_dagger = -4a^4/A,",
        "-delta_dagger*a^4*Z^4 = (4a^8/A) Z^4.",
        "```",
        "",
        "`N_g(b)=16a^4b^2-A(b-1)^4` is reciprocal of degree 4. Hence its square is reciprocal of degree 8, as is the central `b^4` term. Its constant and leading coefficients are both `A^2`, so",
        "",
        "```text",
        "P_0 = P_8 = 4 a^8 A Z^4.",
        "```",
        "",
        "Writing `lambda=32 A^3 s^2 D^2`,",
        "",
        "```text",
        "P(b) = R(b) - lambda*b^5,                    R reciprocal,",
        "P(b) - b^8 P(1/b) = lambda*(b^3-b^5).",
        "```",
        "",
        "Thus admissibility (`D != 0`) gives reciprocity iff `s=0`, equivalently `a=1`. For a twisted identity `P(b)=c b^8 P(k/b)` with rational nonzero `c,k`, the equal nonzero endpoints force `c=1` and `k^8=1`; the nonzero `b^1` coefficient forces `k=1`. Therefore no rational twist repairs the nonreciprocal cases.",
        "",
        "## Lemma L21-D2: exact reduction and Newton polygon at `w` (PROVED)",
        "",
        "Let `e=v_w(z)>0`, so `v_w(Z)=3e`. The character condition `(A|w)=-1` makes `aA` a `w`-adic unit, while `D=1-Z-a^2Z^2` is congruent to 1. Therefore",
        "",
        "```text",
        "P(b) = 16 A^2 b^4(1-2 A s^2 b)  (mod w).",
        "```",
        "",
        "The factorization type is `1^4` when `s=0 mod w`, otherwise `1^4 1`; in particular the degree drops to 4 or 5 and `P mod w` is never an irreducible octic. If `s=0`, the lower Newton polygon has `(length,slope)=(4,-3e),(4,3e)`. If `v_w(s)=0`, it has `(4,-3e),(1,0),(3,4e)`. More generally, for `r=v_w(s)>0`, the right side is `(1,2r),(3,(12e-2r)/3)` when `2r<3e`, and `(4,3e)` when `2r>=3e`. None is a single slope with denominator 8. The distinct negative and nonnegative slopes also force local slope factors, so `w` cannot certify local irreducibility.",
        "",
        "The nonresidue condition therefore correlates with a completely reducible low-degree reduction, not an irreducible or large-degree factor pattern.",
        "",
        "## Lemma L21-D3: uniform trace-quartic irreducibility when `a=1` (PROVED)",
        "",
        "Now `A=5`, `s=0`, `(5|w)=-1`, and",
        "",
        "```text",
        "P(b)=b^4 T(u),  u=b+b^-1,",
        "T(u)=400D^2+(4/5)Z^4(16-5(u-2)^2)^2.",
        "```",
        "",
        "Put `E=Q(t)`, `t^2=-5`. Up to the rational factor `5/4`,",
        "",
        "```text",
        "T(u) = (Z^2(16-5(u-2)^2)-10Dt)",
        "       (Z^2(16-5(u-2)^2)+10Dt).",
        "```",
        "",
        "For either quadratic, substitute `v=Z(u-2)`. At a prime of `E` above `w`, its square test reduces to whether `+2t` or `-2t` is a square in the residue field. If `w` is inert in `E`, the norm of either element is 20, whose character is `(5|w)=-1`, so both are nonsquares. If `w` splits, the characters of `2t` and `-2t` multiply to `(20|w)=-1`; conjugating the prime exchanges the signs, so each quadratic has a prime where its residue is nonsquare. Each quadratic is therefore irreducible over `E`. They are conjugate, so their product `T` is irreducible over `Q`.",
        "",
        "This proof is uniform in rational `z` with `v_w(z)>0` and `D!=0`; it is not a sample inference.",
        "",
        "## Lemma L21-D4: an explicit sufficient lift condition (PROVED)",
        "",
        "For a root `u` of the irreducible trace quartic, the reciprocal lift is irreducible iff `u^2-4` is nonsquare in `Q(u)`. Its field norm is",
        "",
        "```text",
        "N(u^2-4) = T(2)T(-2)/lc(T)^2",
        "           = (4/(25 Z^4))^2 * Xi(Z),",
        "Xi(Z) = (125D^2+64Z^4)(125D^2+1024Z^4).",
        "```",
        "",
        "Consequently `Xi(Z)` nonsquare in `Q` is a sufficient condition for degree-8 irreducibility. Exact integer-square tests find it nonsquare on all 189 canonical `a=1` rows. Proving this (or another lift obstruction) uniformly is still open.",
        "",
        "## Mechanical evidence tables",
        "",
        f"All {summary['pairs']} rows are in `math/h10q/data/l21_irred_direct.jsonl`; every identity is checked with exact `Fraction` arithmetic against `h10q._l10_P`.",
        "",
        "### Constructed `a` choices",
        "",
        "| a | pairs |",
        "|---:|---:|",
    ]
    for a, count in tables["a_counts"].items():
        lines.append(f"| {a} | {count} |")
    lines.extend([
        "",
        "### Factorization of `P mod w`",
        "",
        "| factorization type | pairs |",
        "|---|---:|",
    ])
    for factor_type, count in tables["mod_w_factorization_counts"].items():
        lines.append(f"| `{factor_type}` | {count} |")
    lines.extend([
        "",
        "### Lower Newton polygon at `w`",
        "",
        "| `length@slope` sequence | pairs |",
        "|---|---:|",
    ])
    for slopes, count in tables["newton_w_segment_counts"].items():
        lines.append(f"| `{slopes}` | {count} |")
    lines.extend([
        "",
        "### Newton polygons at primes dividing `A`",
        "",
        "| case | pairs |",
        "|---|---:|",
    ])
    for case, count in tables["newton_A_segment_counts"].items():
        lines.append(f"| `{case}` | {count} |")
    lines.extend([
        "",
        "No prime dividing `A` had a slope denominator 8 in this scan; these rows show two length-4 slopes. This is mechanical evidence only outside the generic unit hypotheses.",
        "",
        "### Trace-quartic local cases (`a=1`)",
        "",
        "| behavior of `w` in `Q(sqrt(-5))` | pairs | anomalies |",
        "|---|---:|---:|",
    ])
    for case, count in tables["trace_local_case_counts"].items():
        lines.append(f"| {case} | {count} | 0 |")
    lines.extend([
        "",
        "### Exact per-row octic certificates on `tau_dagger`",
        "",
        f"`h10q._l13_irred8` certified {summary['tau_dagger_exact_irred_certificates']}/{summary['pairs']} rows at auxiliary primes <= {summary['maximum_auxiliary_certificate_prime']}; reducible instances observed: **0**. These certificates corroborate the direct lemmas but do not close their remaining uniform scope.",
        "",
        "## Final status",
        "",
        "**OPEN** for all admissible constructed parameters. The direct route proves a uniform quartic theorem and a conditional octic theorem on the 189-row reciprocal subfamily, but the nonsquare lift condition has not been proved uniformly and the 164 nonreciprocal rows retain only exact per-row certificates. The places `w` and `p|A` do not supply the hoped-for single denominator-8 Newton slope.",
        "",
        f"Replay wall-clock: **{summary['wall_seconds']:.3f} s**. Pacing: {summary['pacing']['candidates']} candidate checks, {summary['pacing']['sleeps']} sleeps of {summary['pacing']['sleep_seconds']} s every {summary['pacing']['every']} candidates.",
    ])
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    grid = sorted(h10q._l9_grid())
    assert len(grid) == 353
    rows = [pair_record(w, unit) for w, unit in grid]
    assert len(rows) >= 200
    assert all(row["A_character_mod_w"] == -1 for row in rows)
    assert all(not row["mod_w"]["irreducible_degree_8"] for row in rows)
    assert all(not row["newton_w"]["single_slope_denominator_8"] for row in rows)
    assert all(not row["tau_dagger_reducible_instance"] for row in rows)

    reciprocal_rows = [row for row in rows if row["reciprocal"]]
    nonreciprocal_rows = [row for row in rows if not row["reciprocal"]]
    assert len(reciprocal_rows) == 189 and len(nonreciprocal_rows) == 164
    assert all(row["trace_route"] is not None for row in reciprocal_rows)
    assert all(row["trace_route"] is None for row in nonreciprocal_rows)
    assert all(row["trace_route"]["degree_8_irreducibility_by_trace_norm"]
               for row in reciprocal_rows)

    elapsed = time.perf_counter() - started
    tables = count_table(rows)
    certificate_primes = [
        row["auxiliary_exact_irreducibility_certificate_prime"] for row in rows
    ]
    summary = {
        "type": "summary",
        "label": "OPEN",
        "uniform_irreducibility_verdict": "OPEN",
        "proved_structural_lemmas": [
            "L21-D1 exact reciprocal core and unique b^5 perturbation",
            "L21-D2 exact mod-w factorization and Newton-polygon law",
            "L21-D3 uniform trace-quartic irreducibility for a=1 and (5|w)=-1",
            "L21-D4 degree-8 irreducibility under explicit trace-norm nonsquare condition",
        ],
        "refuted_routes": [
            "P irreducible modulo w",
            "a single w-Newton slope of denominator 8",
            "a rational twisted reciprocity absorbing the b^5 perturbation",
        ],
        "pairs": len(rows),
        "reciprocal_pairs": len(reciprocal_rows),
        "nonreciprocal_pairs": len(nonreciprocal_rows),
        "trace_norm_condition_checked_nonsquare": len(reciprocal_rows),
        "trace_norm_condition_anomalies": 0,
        "mod_w_factorization_anomalies": 0,
        "newton_w_anomalies": 0,
        "tau_dagger_exact_irred_certificates": len(certificate_primes),
        "tau_dagger_reducible_instances": 0,
        "maximum_auxiliary_certificate_prime": max(certificate_primes),
        "mechanical_tables": tables,
        "wall_seconds": elapsed,
        "pacing": {
            "candidates": PACER.candidates,
            "every": PACE_EVERY,
            "sleep_seconds": PACE_SECONDS,
            "sleeps": PACER.sleeps,
        },
        "kernel_authority": "math/h10q/h10q.py",
        "refusals": 0,
        "refusals_are_evidence": False,
    }
    payload = [
        {
            "type": "meta",
            "label": "OPEN",
            "schema": "l21-irred-direct-v1",
            "source_script": "math/h10q/l21_irred_direct.py",
            "kernel_authority": "math/h10q/h10q.py",
            "branch": "tau_dagger=(1+2a^2)/A, delta_dagger=-4a^4/A",
            "uniform_verdict": "OPEN",
            "refusals_are_evidence": False,
        },
        *rows,
        summary,
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "\n".join(json.dumps(row, sort_keys=True, separators=(",", ":"))
                    for row in payload) + "\n",
        encoding="utf-8",
    )
    write_report(summary)
    print(
        f"L21 direct: {len(rows)} pairs; reciprocal={len(reciprocal_rows)}, "
        f"mod-w irreducible=0, tau_dagger reducible=0; "
        f"uniform verdict OPEN; {elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
