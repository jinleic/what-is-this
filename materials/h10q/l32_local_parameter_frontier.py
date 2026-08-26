#!/usr/bin/env python3
"""L32 automatic-Phi target maps and two global-section obstructions.

This producer proves three exact statements that sharpen the L30/L31 algebraic
frontier without claiming a global rational point:

* two rational maps land in the dyadic Phi conditions for every rational
  parameter and supply smooth constant-two target residues at every odd prime;
* every finite fixed menu of proportional-trace simplifications misses all but
  finitely many target primes;
* the two linear bridge simplifications b=1+/-2a cover only the target residue
  classes p=5,19 (mod 24).

All finite loops replay symbolic identities or exhaust the explicitly bounded
small-prime tail of a uniform character-sum proof.  No bounded search is used as
an existence or nonexistence theorem.
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction as F
from pathlib import Path

from l25_scaled_coupling import bridge_c, rational_sqrt, vp
from l30_quartic_frontier import legendre, primes_below

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "l32_local_parameter_frontier.jsonl"
REPORT = Path("/tmp/l32_local_parameter_frontier.md")


def frac_text(value: F | int) -> str:
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def bareiss_determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    assert size and all(len(row) == size for row in matrix)
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size) if work[row][pivot_index]),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                assert numerator % previous == 0
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def resultant(left: list[int], right: list[int]) -> int:
    left_degree = len(left) - 1
    right_degree = len(right) - 1
    assert left_degree >= 1 and right_degree >= 1
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))
    size = left_degree + right_degree
    matrix: list[list[int]] = []
    for shift in range(right_degree):
        matrix.append(
            [0] * shift
            + left_descending
            + [0] * (size - shift - len(left_descending))
        )
    for shift in range(left_degree):
        matrix.append(
            [0] * shift
            + right_descending
            + [0] * (size - shift - len(right_descending))
        )
    return bareiss_determinant(matrix)


def discriminant(polynomial: list[int]) -> int:
    degree = len(polynomial) - 1
    derivative = [(index + 1) * polynomial[index + 1] for index in range(degree)]
    signed = (-1) ** (degree * (degree - 1) // 2)
    value = signed * resultant(polynomial, derivative)
    assert value % polynomial[-1] == 0
    return value // polynomial[-1]


def safe_s(parameter: F) -> F:
    return parameter / (1 + 2 * parameter * parameter)


def safe_a(parameter: F) -> F:
    return 1 + 2 * safe_s(parameter)


def safe_b(parameter: F, constant: int) -> F:
    square = parameter * parameter
    return (square + constant) / (5 * square + constant)


def automatic_phi_maps() -> dict:
    parameters = {
        F(numerator, denominator)
        for denominator in range(1, 18)
        for numerator in range(-35, 36)
    }
    rows = 0
    valuation_cases: set[tuple[str, int]] = set()
    for parameter in parameters:
        s = safe_s(parameter)
        a = safe_a(parameter)
        assert vp(s, 2) >= 0 if s else True
        assert vp(a, 2) == 0
        valuation_cases.add(("s", 99 if parameter == 0 else max(-2, min(2, vp(parameter, 2)))))
        for constant in (1, 2, -2):
            b = safe_b(parameter, constant)
            assert b
            assert vp(b, 2) == 0
            valuation_cases.add(
                (f"b_{constant}", 99 if parameter == 0 else max(-2, min(2, vp(parameter, 2))))
            )
            rows += 1

    assert rational_sqrt(F(2)) is None
    assert rational_sqrt(F(2, 5)) is None
    return {
        "type": "automatic-Phi-rational-maps",
        "label": "PROVED for every rational parameter by complete 2-adic case analysis",
        "maps": {
            "s(t)": "t/(1+2*t^2)",
            "a(t)": "1+2*s(t)",
            "b_c(u)": "(u^2+c)/(5*u^2+c), c in {1,2,-2}",
        },
        "conclusions": [
            "v2(s(t))>=0 and v2(a(t))=0 for every t in Q",
            "v2(b_c(u))=0 for every u in Q and c in {1,2,-2}",
            "all denominators and numerators used by the b-maps are nonzero over Q",
        ],
        "proof_cases": {
            "s": [
                "v2(t)>=0: 1+2*t^2 is a unit",
                "v2(t)<0: v2(s)=-1-v2(t)>=0",
            ],
            "b_c": [
                "v2(u)>0: numerator and denominator have the same valuation",
                "v2(u)=0: for c=1 both valuations are 1; for c=+/-2 both are 0",
                "v2(u)<0: the u^2 terms dominate with equal valuation",
            ],
        },
        "finite_identity_rows": rows,
        "valuation_case_buckets": sorted([list(row) for row in valuation_cases]),
    }


def target_a_rows(prime: int) -> list[tuple[int, int, int, int]]:
    rows = []
    for parameter in range(prime):
        denominator = (1 + 2 * parameter * parameter) % prime
        if denominator == 0:
            continue
        numerator = (2 * parameter * parameter + 2 * parameter + 1) % prime
        a = numerator * pow(denominator, -1, prime) % prime
        A = (1 + 4 * a * a) % prime
        shifted = (A + 4) % prime
        if A and shifted and legendre(A, prime) == -1 and legendre(shifted, prime) == 1:
            X = next(value for value in range(prime) if value * value % prime == shifted)
            rows.append((parameter, a, A, X))
    return rows


def target_a_map() -> dict:
    # If d=2t^2+1 and n=2t^2+2t+1, then A and A+4 have
    # squareclasses F5=d^2+4n^2 and F9=5d^2+4n^2.
    def f5(t: int) -> int:
        return 20 * t**4 + 32 * t**3 + 36 * t * t + 16 * t + 5

    def f9(t: int) -> int:
        return 36 * t**4 + 32 * t**3 + 52 * t * t + 16 * t + 9

    identity_rows = 0
    for t in map(F, range(-20, 21)):
        d = 2 * t * t + 1
        n = 2 * t * t + 2 * t + 1
        assert d * d + 4 * n * n == f5(int(t))
        assert 5 * d * d + 4 * n * n == f9(int(t))
        a = n / d
        assert 1 + 4 * a * a == F(f5(int(t)), int(d * d))
        assert 5 + 4 * a * a == F(f9(int(t)), int(d * d))
        identity_rows += 1

    f5_coefficients = [5, 16, 36, 32, 20]
    f9_coefficients = [9, 16, 52, 32, 36]
    disc_f5 = discriminant(f5_coefficients)
    disc_f9 = discriminant(f9_coefficients)
    common_resultant = resultant(f5_coefficients, f9_coefficients)
    assert disc_f5 == 17825792
    assert disc_f9 == 2333081600
    assert common_resultant == 268435456

    small_rows: dict[str, dict[str, int]] = {}
    for prime in primes_below(211):
        if prime < 5:
            continue
        rows = target_a_rows(prime)
        assert rows
        parameter, a, A, X = rows[0]
        small_rows[str(prime)] = {
            "count": len(rows),
            "t": parameter,
            "a": a,
            "A": A,
            "X": X,
        }

    checked_primes = 0
    for prime in primes_below(5000):
        if prime < 5:
            continue
        assert target_a_rows(prime)
        checked_primes += 1

    threshold = 211
    lower_bound = (threshold - 13 * math.sqrt(threshold)) / 4 - 4
    assert lower_bound > 0
    return {
        "type": "automatic-Phi-all-target-a-map",
        "label": "PROVED for every odd prime w>=5",
        "map": "a(t)=(2*t^2+2*t+1)/(2*t^2+1)",
        "character_polynomials": {
            "F5": f5_coefficients,
            "F9": f9_coefficients,
            "meaning": "chi(A)=chi(F5), chi(A+4)=chi(F9)",
        },
        "exact_invariants": {
            "disc_F5": disc_f5,
            "disc_F9": disc_f9,
            "resultant": common_resultant,
            "odd_exceptional_primes": [5, 17, 89],
        },
        "uniform_bound": "N_w >= (w-13*sqrt(w))/4-4 > 0 for w>=211",
        "bound_explanation": [
            "Weil bounds 3*sqrt(w), 3*sqrt(w), 7*sqrt(w) for F5, F9, F5*F9",
            "the resultant is 2^28, so the two quartics have no common odd-prime root",
            "zero-value indicator weights contribute at most 4 in total",
            "a denominator root contributes zero to the indicator",
        ],
        "small_prime_exhaustion": small_rows,
        "crosschecked_odd_primes_below_5000": checked_primes,
        "identity_rows": identity_rows,
        "w3_scope": "not covered by this map; L30's separate a=5 strong-Hensel fibre remains the w=3 branch",
    }


def target_b_map() -> dict:
    constants = (1, 2, -2)
    samples: dict[str, dict[str, int | str]] = {}
    checked = 0
    for prime in primes_below(1000):
        if prime < 3:
            continue
        symbols = {constant: legendre(-constant, prime) for constant in constants}
        assert math.prod(symbols.values()) == 1
        choices = [constant for constant, symbol in symbols.items() if symbol == 1]
        assert choices
        constant = choices[0]
        root = next(value for value in range(1, prime) if (value * value + constant) % prime == 0)
        assert (5 * root * root + constant) % prime == (-4 * constant) % prime != 0

        lift = None
        for step in range(prime):
            candidate = root + prime * step
            numerator = candidate * candidate + constant
            if numerator and vp(F(numerator), prime) == 1:
                lift = candidate
                break
        assert lift is not None
        value = safe_b(F(lift), constant)
        assert vp(value, prime) == 1
        assert vp(value, 2) == 0
        if prime < 50:
            samples[str(prime)] = {
                "c": constant,
                "u_mod_w": root,
                "integer_lift": lift,
                "b": frac_text(value),
            }
        checked += 1

    return {
        "type": "automatic-Phi-all-target-b-map",
        "label": "PROVED for every odd prime w",
        "menu": ["(u^2+1)/(5*u^2+1)", "(u^2+2)/(5*u^2+2)", "(u^2-2)/(5*u^2-2)"],
        "character_identity": "(-1|w)*(-2|w)*(2|w)=(4|w)=1",
        "selection": [
            "at least one of -1,-2,2 is a square modulo every odd w",
            "choose u0^2=-c; then the denominator is -4c and the numerator derivative is 2u0",
            "a rational lift has numerator valuation exactly 1, hence v_w(b_c)=1",
        ],
        "samples": samples,
        "crosschecked_odd_primes_below_1000": checked,
    }


def proportional_trace_no_go() -> dict:
    checks = 0
    values = [F(-7, 3), F(-3, 2), F(-1), F(1, 3), F(2), F(5, 2), F(7)]
    for r in values:
        for t in values:
            if not r or not t:
                continue
            a = (t * t - r) / (2 * r * t)
            h = (t * t + r) / (2 * t)
            A = 1 + 4 * a * a
            b_plus = t * t / r
            b_minus = r / (t * t)
            assert h * h - r * r * a * a == r
            for b in (b_plus, b_minus):
                trace = (b - 1) ** 2 / b
                assert trace == 4 * r * a * a
                Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
                assert Ng / (b * b) == 16 * a**4 * (1 - A * r * r)
                for Z in (F(-2), F(1, 3), F(3)):
                    D = 1 - Z - a * a * Z * Z
                    if not a or not D:
                        continue
                    expected_c = 16 * a**6 * Z * Z * (1 - A * r * r) / (A * D)
                    assert bridge_c(a, Z, b) == expected_c
                    checks += 1

    return {
        "type": "fixed-proportional-trace-target-no-go",
        "label": "PROVED for every fixed rational r and every finite fixed menu",
        "ansatz": "(b-1)^2/b=4*r*a^2",
        "parameterization": {
            "h": "(t^2+r)/(2*t)",
            "a": "(t^2-r)/(2*r*t)",
            "b_branches": ["t^2/r", "r/t^2"],
            "bridge": "c=16*a^6*Z^2*(1-A*r^2)/(A*D)",
        },
        "valuation_proof": [
            "for b=t^2/r, v_w(b)=1 gives v_w(r)=2v_w(t)-1 and v_w(a)=-v_w(t)",
            "for b=r/t^2, v_w(b)=1 gives v_w(r)=2v_w(t)+1 and v_w(a)=-v_w(t)-1",
            "in either branch v_w(a)=0 forces v_w(r)=-1",
        ],
        "consequence": "a fixed r can hit unit-a target strata only at denominator primes of r; a finite menu misses every target outside finite denominator support",
        "identity_rows": checks,
    }


def linear_bridge_scope() -> dict:
    checks = 0
    for a in map(F, range(-25, 26)):
        if not a:
            continue
        A = 1 + 4 * a * a
        for sign in (1, -1):
            b = 1 + sign * 2 * a
            if not b:
                continue
            Ng = 16 * a**4 * b * b - A * (b - 1) ** 4
            assert Ng == sign * 64 * a**5
            checks += 1

    accepted = []
    rejected_samples = []
    for prime in primes_below(500):
        if prime < 5:
            continue
        good = legendre(2, prime) == -1 and legendre(6, prime) == 1
        if good:
            assert prime % 24 in (5, 19)
            accepted.append(prime)
        elif len(rejected_samples) < 12:
            rejected_samples.append(prime)

    return {
        "type": "linear-bridge-simplification-scope",
        "label": "PROVED exact identities and target-character limitation",
        "sections": {
            "b=1+2a": "N_g=64*a^5",
            "b=1-2a": "N_g=-64*a^5",
        },
        "target_reduction": "w|b implies A=1+4a^2=2 (mod w); constant-two target regularity needs (2|w)=-1 and (6|w)=+1",
        "prime_classes": "for w>=5 this is equivalent to w=5 or 19 (mod 24)",
        "accepted_primes_below_500": accepted,
        "rejected_sample_primes": rejected_samples,
        "consequence": "the two linear simplifications miss infinitely many target primes and cannot be a uniform L30 section menu",
        "identity_rows": checks,
    }


def frontier_record() -> dict:
    return {
        "type": "L32-strict-frontier",
        "label": "PROVED local parameterization and two infinite ansatz obstructions; global closures remain open",
        "closed": [
            "automatic rational Phi maps with a smooth L30 target base at every odd w>=5",
            "automatic rational b-unit maps with v_w(b)=1 at every odd target",
            "every finite fixed proportional-trace menu is non-uniform in w",
            "the b=1+/-2a simplifications cover only w=5,19 mod 24",
        ],
        "open": [
            "a rational point on the resulting bridge-specialized H2 surface meeting every controlled place",
            "the fixed-family Hilbert-detector asymptotic beyond Bombieri-Vinogradov",
            "AP1 and removal of classical Schinzel H from the six-variable chain",
        ],
        "scope": "the local maps do not control emergent ramification from their numerators and denominators, so they do not prove global completeness",
    }


def build_report(records: list[dict], elapsed: float) -> str:
    a_map = next(row for row in records if row["type"] == "automatic-Phi-all-target-a-map")
    return "\n".join(
        [
            "# L32 local parameter frontier",
            "",
            "## Positive local result",
            "",
            "Two explicit rational maps land in Phi at 2 for every rational parameter.",
            "A degree-eight character-sum calculation plus a finite small-prime tail gives",
            "a smooth constant-two target base at every odd prime w>=5.",
            f"The uniform lower bound is `{a_map['uniform_bound']}`.",
            "A three-map b-menu gives v_w(b)=1 at every odd target by the identity",
            "(-1|w)(-2|w)(2|w)=1.",
            "",
            "## Closed ansatzes",
            "",
            "A fixed proportional-trace parameter r can meet a unit-a target only when",
            "v_w(r)=-1, so every finite fixed menu misses almost every target.",
            "The simplifications b=1+/-2a force A=2 mod w and cover only",
            "w=5 or 19 mod 24.",
            "",
            "## Remaining global step",
            "",
            "The maps do not control emergent ramification at their numerator/denominator primes.",
            "They reduce the base-selection problem but do not provide a rational H2 root,",
            "a fixed-family detector asymptotic, or AP1.",
            "",
            f"Wall time: {elapsed:.3f} s.",
        ]
    ) + "\n"


def main() -> int:
    started = time.perf_counter()
    records = [
        automatic_phi_maps(),
        target_a_map(),
        target_b_map(),
        proportional_trace_no_go(),
        linear_bridge_scope(),
        frontier_record(),
    ]
    elapsed = time.perf_counter() - started
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in records))
    REPORT.write_text(build_report(records, elapsed))
    print(
        "l32_local_parameter_frontier: automatic Phi/all-target maps proved; "
        "two section menus closed; global H2/detector/AP1 remain open; "
        f"rows={len(records)} elapsed={elapsed:.3f}s"
    )
    print(f"wrote {OUT}")
    print(f"wrote {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
