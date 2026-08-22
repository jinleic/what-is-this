"""Standalone verifier for the e88 exact Lee--Yang box certificates.

Run from the repository root:

    timeout 1800 .venv/bin/python tests/test_lee_yang_zeros.py

This file deliberately does not import experiments/e88_lee_yang_zeros.py.
It independently rebuilds the small 3x3x3 polynomial by brute force, checks
an independent zero-field transfer identity, reconstructs the circle
restriction, re-runs exact Sturm counts for both headline fallback boxes at
x=2/3, and re-derives the four-size divided-difference edge calculation.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import sympy as sp

from ising.lattices import hyperrect
from ising.lee_yang import _CRT_PRIMES, field_dos_from_enumeration, rational_field_polynomial
from ising.transfer_matrix import box_broken_bond_poly

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "lee_yang" / "zeros_l5.json"
T = sp.Symbol("t")
COSINE_EVEN_ORDER = 80
EDGE = Fraction(1, 50)


def text_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def parse(value: str) -> Fraction:
    return Fraction(value)


def divide_by_z_plus_one(coefficients: list[int]) -> list[int]:
    degree = len(coefficients) - 1
    assert degree % 2 == 1 and coefficients == coefficients[::-1]
    quotient = [coefficients[0]]
    for coefficient in coefficients[1:-1]:
        quotient.append(coefficient - quotient[-1])
    assert quotient == quotient[::-1]
    assert [quotient[0]] + [quotient[i - 1] + quotient[i] for i in range(1, degree)] + [quotient[-1]] == coefficients
    return quotient


def reduced_power(coefficients: list[int]) -> list[int]:
    values = divide_by_z_plus_one(coefficients) if (len(coefficients) - 1) % 2 else coefficients
    half = (len(values) - 1) // 2
    result = [values[half]]
    previous = [1]
    current = [0, 1]
    for order in range(1, half + 1):
        if order == 1:
            cheb = current
        else:
            cheb = [0] * (len(current) + 1)
            for i, value in enumerate(current):
                cheb[i + 1] += 2 * value
            for i, value in enumerate(previous):
                cheb[i] -= value
            previous, current = current, cheb
        if len(result) < len(cheb):
            result.extend([0] * (len(cheb) - len(result)))
        factor = 2 * values[half - order]
        for i, value in enumerate(cheb):
            result[i] += factor * value
    return result


def sign_at_fraction(power: list[int], value: Fraction) -> int:
    numerator, denominator = value.numerator, value.denominator
    degree = len(power) - 1
    num_power = 1
    den_power = denominator**degree
    total = 0
    for i, coefficient in enumerate(power):
        total += coefficient * num_power * den_power
        if i < degree:
            num_power *= numerator
            den_power //= denominator
    return (total > 0) - (total < 0)


def cosine_bounds(value: Fraction) -> tuple[Fraction, Fraction]:
    assert 0 <= value < 2
    term = Fraction(1)
    partial = term
    upper = None
    for order in range(1, COSINE_EVEN_ORDER + 2):
        term *= -value * value
        term /= (2 * order - 1) * (2 * order)
        partial += term
        if order == COSINE_EVEN_ORDER:
            upper = partial
    assert upper is not None
    lower = partial
    assert lower <= upper
    return lower, upper


def root_certificate_check(record: dict) -> None:
    coefficients = [int(value) for value in record["A_k"]]
    assert coefficients == coefficients[::-1]
    assert min(coefficients) >= 0
    power = reduced_power(coefficients)
    poly = sp.Poly.from_list(list(reversed(power)), T, domain=sp.ZZ)
    expected = (len(coefficients) - 1) // 2
    assert int(poly.count_roots(-1, 1)) == expected
    cert = record["first_zero_certificate"]
    lo, hi = map(parse, cert["t_interval"])
    assert lo < hi
    assert sign_at_fraction(power, lo) * sign_at_fraction(power, hi) < 0
    coarse = poly.intervals()
    matching = []
    for bounds, multiplicity in coarse:
        a, b = Fraction(bounds[0]), Fraction(bounds[1])
        if a <= lo and hi <= b:
            matching.append((a, b, multiplicity))
    assert len(matching) == 1 and matching[0][2] == 1
    theta_lo, theta_hi = map(parse, cert["theta_interval"])
    cos_lo_lower, _ = cosine_bounds(theta_lo)
    _, cos_hi_upper = cosine_bounds(theta_hi)
    assert cos_hi_upper < lo < hi < cos_lo_lower


def independent_small_zero_field_sum(record: dict, numerator: int, denominator: int) -> None:
    coefficients = [int(value) for value in record["A_k"]]
    broken = box_broken_bond_poly((4, 4, 2), periodic=(False, False, False))
    n_bonds = record["n_bonds"]
    expected = sum(
        int(count) * numerator**q * denominator ** (n_bonds - q)
        for q, count in enumerate(broken)
    )
    assert sum(coefficients) == expected
    assert sum(broken) == 1 << record["n_sites"]


def weights_for_sizes(sizes: tuple[int, ...]) -> dict[int, Fraction]:
    nodes = {size: Fraction(1, size) for size in sizes}
    answer = {}
    for size in sizes:
        product = Fraction(1)
        for other in sizes:
            if other != size:
                product *= nodes[size] - nodes[other]
        answer[size] = 1 / product
    return answer


def dd3_edge_interval(angle_intervals: dict[int, tuple[Fraction, Fraction]], power: int):
    sizes = tuple(sorted(angle_intervals))
    weights = weights_for_sizes(sizes)
    a_lo = Fraction(0)
    a_hi = Fraction(0)
    for size in sizes:
        weight = weights[size]
        factor = size**power * weight
        if factor > 0:
            a_lo += factor * angle_intervals[size][0]
            a_hi += factor * angle_intervals[size][1]
        else:
            a_lo += factor * angle_intervals[size][1]
            a_hi += factor * angle_intervals[size][0]
    b = sum((size**power * weights[size] for size in sizes), Fraction(0))
    if b > 0:
        edge = (a_lo / b, a_hi / b)
    else:
        edge = (a_hi / b, a_lo / b)
    frozen = (a_lo - EDGE * b, a_hi - EDGE * b)
    return edge, frozen


def verify_four_size_dd3(artifact: dict) -> None:
    family = artifact["data"]["exact_zero_data"]["x_2_over_3"]["ladder16"]
    intervals = {}
    for key, record in family.items():
        length = int(key.split("x")[-1])
        theta = record["first_zero_certificate"]["theta_interval"]
        intervals[length] = (parse(theta[0]), parse(theta[1]))
    stored = artifact["data"]["four_size_separation"]["x_2_over_3"]["ladder16"]["models"]
    for model in stored:
        edge, frozen = dd3_edge_interval(intervals, int(model["p"]))
        assert [text_fraction(edge[0]), text_fraction(edge[1])] == model["edge_interval_e_star"]
        assert [text_fraction(frozen[0]), text_fraction(frozen[1])] == model["frozen_edge_1_over_50_dd3_interval"]


def verify_resource_wall(artifact: dict) -> None:
    wall = artifact["data"]["resource_wall_5x5x5"]
    assert wall["final_grade_buffer_bytes_full_width"] == (1 << 25) * 126 * 8
    assert wall["final_grade_buffer_bytes_truncated"] == (1 << 25) * 63 * 8
    n_bonds = 300
    for name, numerator, denominator in (
        ("x_2_over_3", 2, 3),
        ("x_7_over_10", 7, 10),
    ):
        bound = (1 << 125) * max(numerator, denominator) ** n_bonds
        product = 1
        count = 0
        sufficient = False
        for prime in _CRT_PRIMES:
            product *= prime
            count += 1
            if product > bound:
                sufficient = True
                break
        detail = wall["per_coupling"][name]
        assert detail["coefficient_bound_bits"] == bound.bit_length()
        assert detail["stored_crt_prime_table_length"] == len(_CRT_PRIMES)
        assert detail["stored_crt_prime_table_sufficient"] == sufficient
        assert detail["minimum_30bit_moduli_by_bit_count"] == (
            bound.bit_length() + 29
        ) // 30
        assert detail["crt_primes_needed_from_stored_table_exact_bound"] == (
            count if sufficient else None
        )


def main() -> None:
    artifact = json.loads(ARTIFACT.read_text())
    assert artifact["provenance"]["script"] == "experiments/e88_lee_yang_zeros.py"
    assert all(item["passed"] for item in artifact["checks"])

    # Fully independent brute-force reconstruction at the largest enumerated cube.
    dos = field_dos_from_enumeration(hyperrect((3, 3, 3), periodic=False))
    for grid_id, numerator, denominator in (("x_2_over_3", 2, 3), ("x_7_over_10", 7, 10)):
        expected = [int(value) for value in rational_field_polynomial(dos, numerator, denominator)]
        stored = artifact["data"]["exact_zero_data"][grid_id]["cube3"]["3x3x3"]["A_k"]
        assert expected == [int(value) for value in stored]

    # Independent zero-field transfer-matrix identity on a 4x4x2 box.
    small = artifact["data"]["exact_zero_data"]["x_2_over_3"]["ladder16"]["4x4x2"]
    independent_small_zero_field_sum(small, 2, 3)

    # Exact Sturm root counts and exact angle containment for both requested fallback boxes.
    headline = artifact["data"]["exact_zero_data"]["x_2_over_3"]
    root_certificate_check(headline["ladder16"]["4x4x5"])  # 5x4x4
    root_certificate_check(headline["ladder20"]["4x5x5"])  # 5x5x4

    verify_four_size_dd3(artifact)
    verify_resource_wall(artifact)
    print("PASS")


if __name__ == "__main__":
    main()
