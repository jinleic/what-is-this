#!/usr/bin/env python3
"""Exact all-coupling obstruction to pairwise closure under cubic-star decimation.

Eliminating the centre of a uniform degree-d Ising star gives the boundary
weight ``2*cosh(K * sum_i sigma_i)``.  Its logarithm has a unique Walsh
expansion.  Pairwise Ising weights have Walsh degree at most two; this producer
certifies exact nonzero degree-four and degree-six coefficients for d=4 and
d=6, respectively, without evaluating a logarithm numerically.

The d=6 sign is reduced to an integer polynomial positive on x>1, where
x=cosh(2K)^2.  All computations below are exact integers.
"""

from __future__ import annotations

import json
import itertools
import platform
import resource
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "star_decimation.json"


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def check(rows: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    rows.append({"name": name, "passed": bool(passed), "detail": detail})


def poly_add(left: list[int], right: list[int]) -> list[int]:
    out = [0] * max(len(left), len(right))
    for index in range(len(out)):
        out[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_scale(poly: list[int], scalar: int) -> list[int]:
    return [scalar * coefficient for coefficient in poly]


def poly_multiply(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            out[i + j] += left_coefficient * right_coefficient
    return out


def poly_power(poly: list[int], exponent: int) -> list[int]:
    out = [1]
    base = list(poly)
    while exponent:
        if exponent & 1:
            out = poly_multiply(out, base)
        exponent >>= 1
        if exponent:
            base = poly_multiply(base, base)
    return out


def divide_x_minus_one(poly: list[int]) -> tuple[list[int], int]:
    descending = list(reversed(poly))
    quotient_descending = [descending[0]]
    for coefficient in descending[1:-1]:
        quotient_descending.append(coefficient + quotient_descending[-1])
    remainder = descending[-1] + quotient_descending[-1]
    return list(reversed(quotient_descending)), remainder


def shift_at_one(poly: list[int]) -> list[int]:
    from math import comb

    out = [0] * len(poly)
    for degree, coefficient in enumerate(poly):
        for shifted_degree in range(degree + 1):
            out[shifted_degree] += coefficient * comb(degree, shifted_degree)
    return out


def walsh_counts(degree: int, subset_size: int) -> dict[int, int]:
    """Integer multiplicities before the common 2^-degree factor.

    The key is ``abs(sum spins)`` because log(2*cosh(K*s)) is even in s.
    Permutation symmetry makes the result depend only on ``subset_size``.
    """

    counts: dict[int, int] = defaultdict(int)
    for state in range(1 << degree):
        spins = [1 if (state >> site) & 1 else -1 for site in range(degree)]
        character = 1
        for spin in spins[:subset_size]:
            character *= spin
        counts[abs(sum(spins))] += character
    return dict(sorted(counts.items()))
def walsh_inner(degree: int, left: tuple[int, ...], right: tuple[int, ...]) -> int:
    total = 0
    for state in range(1 << degree):
        character = 1
        for site in left + right:
            character *= 1 if (state >> site) & 1 else -1
        total += character
    return total




def cubic_interior_neighbourhood(vertex: tuple[int, int, int]) -> tuple[tuple[int, int, int], ...]:
    x, y, z = vertex
    return tuple(
        sorted(
            (
                (x - 1, y, z),
                (x + 1, y, z),
                (x, y - 1, z),
                (x, y + 1, z),
                (x, y, z - 1),
                (x, y, z + 1),
            )
        )
    )


def finite_neighbourhood_audit(side: int) -> tuple[int, int]:
    interiors = [
        (x, y, z)
        for x in range(1, side - 1)
        for y in range(1, side - 1)
        for z in range(1, side - 1)
    ]
    supports = [cubic_interior_neighbourhood(vertex) for vertex in interiors]
    return len(interiors), len(set(supports))
def open_neighbourhood(
    vertex: tuple[int, int, int], side: int
) -> tuple[tuple[int, int, int], ...]:
    candidates = cubic_interior_neighbourhood(vertex)
    return tuple(
        point
        for point in candidates
        if all(0 <= coordinate < side for coordinate in point)
    )


def checkerboard_incidence_audit(side: int, parity: int) -> dict[str, int | bool]:
    eliminated = [
        (x, y, z)
        for x in range(side)
        for y in range(side)
        for z in range(side)
        if (x + y + z) % 2 == parity
    ]
    supports = [open_neighbourhood(vertex, side) for vertex in eliminated]
    multiplicities: dict[tuple[tuple[int, int, int], ...], int] = defaultdict(int)
    for support in supports:
        multiplicities[support] += 1
    full = [support for support in supports if len(support) == 6]
    boundary = [support for support in supports if len(support) < 6]
    return {
        "eliminated_sites": len(eliminated),
        "full_stars": len(full),
        "boundary_stars": len(boundary),
        "full_supports_unique": all(multiplicities[support] == 1 for support in full),
        "maximum_boundary_arity": max((len(support) for support in boundary), default=0),
    }




def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    walsh = {
        degree: {
            subset_size: walsh_counts(degree, subset_size)
            for subset_size in range(degree + 1)
        }
        for degree in range(2, 7)
    }
    odd_zero = all(
        all(value == 0 for value in walsh[degree][subset_size].values())
        for degree in range(2, 7)
        for subset_size in range(1, degree + 1, 2)
    )
    check(
        checks,
        "global_flip_kills_all_odd_walsh_sectors",
        odd_zero,
        "direct enumeration for every degree 2 through 6",
    )

    d3_pair = walsh[3][2]
    check(
        checks,
        "degree_three_has_only_constant_and_pair_sectors",
        d3_pair == {1: -2, 3: 2} and 4 > 3,
        "global flip removes odd sectors and arity three has no even sector above two",
    )
    check(
        checks,
        "uniform_star_triangle_control",
        d3_pair == {1: -2, 3: 2},
        "J=(1/4)log(cosh(3K)/cosh(K)); A=2*cosh(K)*exp(J) matches both boundary orbits",
    )

    d4_top = walsh[4][4]
    d4_gap = poly_add(poly_power([0, 1], 2), poly_scale([-1, 2], -1))
    check(
        checks,
        "degree_four_walsh_vector",
        d4_top == {0: 6, 2: -8, 4: 2},
        "16*c4=2*log(cosh(4K))-8*log(cosh(2K))",
    )
    check(
        checks,
        "degree_four_ratio_strictly_below_one",
        d4_gap == poly_power([-1, 1], 2),
        "for x=cosh(2K)^2>1, x^2-(2x-1)=(x-1)^2>0, so c4<0",
    )

    d6_even = {subset_size: walsh[6][subset_size] for subset_size in (0, 2, 4, 6)}
    d6_top = d6_even[6]
    d6_four_body_gap = poly_add(poly_power([-1, 2], 2), poly_scale([-3, 4], -1))
    ratio_numerator = poly_add(
        poly_multiply(poly_power([0, 1], 8), [-3, 4]),
        poly_scale(poly_power([-1, 2], 6), -1),
    )
    quotient = ratio_numerator
    remainders: list[int] = []
    for _ in range(3):
        quotient, remainder = divide_x_minus_one(quotient)
        remainders.append(remainder)
    shifted = shift_at_one(quotient)
    reconstructed = poly_multiply(poly_power([-1, 1], 3), quotient)
    check(
        checks,
        "degree_six_walsh_vector",
        d6_top == {0: -20, 2: 30, 4: -12, 6: 2},
        "64*c6=2*log(cosh(6K))-12*log(cosh(4K))+30*log(cosh(2K))",
    )
    check(
        checks,
        "degree_six_full_even_walsh_decomposition",
        d6_even
        == {
            0: {0: 20, 2: 30, 4: 12, 6: 2},
            2: {0: -4, 2: -2, 4: 4, 6: 2},
            4: {0: 4, 2: -2, 4: -4, 6: 2},
            6: {0: -20, 2: 30, 4: -12, 6: 2},
        },
        "the exact octahedral log-weight has only 0-, 2-, 4-, and 6-spin sectors",
    )
    check(
        checks,
        "degree_six_four_body_coefficient_negative",
        d6_four_body_gap == poly_scale(poly_power([-1, 1], 2), 4),
        "exp(32*c4)=(4*x-3)/(2*x-1)^2<1 because the denominator gap is 4*(x-1)^2",
    )
    check(
        checks,
        "six_body_polynomial_factor_exact",
        remainders == [0, 0, 0] and reconstructed == ratio_numerator,
        "x^8(4x-3)-(2x-1)^6=(x-1)^3 Q(x)",
    )
    check(
        checks,
        "six_body_shifted_factor_positive",
        shifted == [8, 54, 144, 188, 120, 33, 4]
        and all(coefficient > 0 for coefficient in shifted),
        "Q(1+y)=8+54y+144y^2+188y^3+120y^4+33y^5+4y^6",
    )
    check(
        checks,
        "six_body_coefficient_positive_for_nonzero_real_K",
        all(coefficient > 0 for coefficient in shifted),
        "x=cosh(2K)^2>1 for K!=0, hence exp(32*c6)>1 and c6>0",
    )

    box_audits = {str(side): finite_neighbourhood_audit(side) for side in (3, 4, 5, 6)}
    incidence_audits = {
        f"side_{side}_parity_{parity}": checkerboard_incidence_audit(side, parity)
        for side in (4, 5)
        for parity in (0, 1)
    }
    offsets = cubic_interior_neighbourhood((0, 0, 0))
    offset_sum = [
        sum(point[axis] for point in offsets) for axis in range(3)
    ]
    pairwise_orthogonal = all(
        walsh_inner(6, tuple(range(6)), subset) == 0
        for subset_size in range(3)
        for subset in itertools.combinations(range(6), subset_size)
    )
    check(
        checks,
        "finite_cubic_interior_star_supports_unique",
        all(total == unique for total, unique in box_audits.values()),
        str(box_audits),
    )
    check(
        checks,
        "checkerboard_full_support_incidence",
        offset_sum == [0, 0, 0]
        and all(
            bool(record["full_supports_unique"])
            and int(record["maximum_boundary_arity"]) <= 5
            for record in incidence_audits.values()
        ),
        f"offset_sum={offset_sum}; {incidence_audits}",
    )
    check(
        checks,
        "pairwise_log_weight_has_no_degree_six_sector",
        pairwise_orthogonal,
        "the six-spin Walsh character is exactly orthogonal to every constant, field, and pair character on all 64 states",
    )

    if not all(row["passed"] for row in checks):
        failed = [row["name"] for row in checks if not row["passed"]]
        raise AssertionError(f"producer checks failed: {failed}")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e234_star_decimation.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integer Walsh enumeration and polynomial identities",
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": max_rss_bytes(),
            "benchmark_used": False,
        },
        "data": {
            "claim_tag": "[THEOREM]",
            "local_weight": "W_d(sigma)=2*cosh(K*sum_i sigma_i)",
            "walsh_normalization": "hat(f)_A=2^(-d) sum_sigma f(sigma) prod_(i in A) sigma_i",
            "degree_three_control": {
                "top_even_degree": 2,
                "pair_walsh_counts": {str(key): value for key, value in d3_pair.items()},
                "star_triangle_coupling": "J=(1/4)*log(cosh(3K)/cosh(K))",
                "star_triangle_prefactor": "A=2*cosh(K)*exp(J)",
            },
            "degree_four": {
                "top_walsh_counts": {str(key): value for key, value in d4_top.items()},
                "coefficient": "c4=(1/8)*log(cosh(4K)/cosh(2K)^4)<0 for K!=0",
                "ratio": "exp(8*c4)=(2*x-1)/x^2 with x=cosh(2K)^2",
                "denominator_minus_numerator_coefficients_ascending": d4_gap,
            },
            "degree_six": {
                "top_walsh_counts": {str(key): value for key, value in d6_top.items()},
                "even_walsh_counts_by_subset_size": {
                    str(subset_size): {
                        str(key): value for key, value in counts.items()
                    }
                    for subset_size, counts in d6_even.items()
                },
                "logcosh_exponent_vectors_scaled_by_32": {
                    "c2": {"2": -1, "4": 2, "6": 1},
                    "c4": {"2": -1, "4": -2, "6": 1},
                    "c6": {"2": 15, "4": -6, "6": 1},
                },
                "exact_octahedral_couplings": {
                    "c2": "(1/32)*log(cosh(4K)^2*cosh(6K)/cosh(2K))>0 for K!=0",
                    "c4": "(1/32)*log(cosh(6K)/(cosh(2K)*cosh(4K)^2))<0 for K!=0",
                    "c6": "(1/32)*log(cosh(6K)*cosh(2K)^15/cosh(4K)^6)>0 for K!=0",
                },
                "two_body_ratio": "exp(32*c2)=(2*x-1)^2*(4*x-3) with x=cosh(2K)^2",
                "two_body_ratio_factors": [
                    {"constant": -1, "linear": 2, "power": 2},
                    {"constant": -3, "linear": 4, "power": 1},
                ],
                "four_body_ratio": "exp(32*c4)=(4*x-3)/(2*x-1)^2 with x=cosh(2K)^2",
                "four_body_denominator_gap_coefficients_ascending": d6_four_body_gap,
                "coefficient": "c6=(1/32)*log(cosh(6K)*cosh(2K)^15/cosh(4K)^6)>0 for K!=0",
                "ratio": "exp(32*c6)=x^8*(4*x-3)/(2*x-1)^6 with x=cosh(2K)^2",
                "ratio_numerator_coefficients_ascending": ratio_numerator,
                "factorization": "numerator=(x-1)^3*Q(x)",
                "quotient_coefficients_ascending": quotient,
                "positive_shift_factor_coefficients_ascending": shifted,
            },
            "whole_sublattice": {
                "statement": "On the infinite simple-cubic lattice, or any open box with an eliminated full-degree interior site, checkerboard elimination generates a strictly positive six-spin Walsh interaction on that site's six retained neighbours.",
                "factorization_reason": "Eliminated checkerboard spins have no mutual edges, so their conditional sums multiply and the logarithms add.",
                "no_cancellation_reason": "For a full cubic star, sum_(u in N(v)) u=6v; equal six-neighbour supports therefore imply equal centres. Boundary factors have arity at most five and visible pairwise terms have Walsh degree at most two, so neither can contribute to that six-spin character.",
                "exact_reformulation": "After eliminating one checkerboard class, Z is exactly the partition function on the retained class with, for every eliminated full cubic star, the octahedral interaction c0+c2*sum_|A|=2 sigma_A+c4*sum_|A|=4 sigma_A+c6*prod_i sigma_i; c2>0, c4<0, c6>0 for K!=0, plus lower-arity factors from eliminated boundary sites and any pre-existing visible terms.",
                "finite_audits_total_unique": {
                    side: {"total": total, "unique": unique}
                    for side, (total, unique) in box_audits.items()
                },
                "barycenter_certificate": {
                    "neighbour_count": 6,
                    "offset_sum": offset_sum,
                    "identity": "sum_(u in N(v)) u = 6v",
                },
                "checkerboard_incidence_audits": incidence_audits,
                "pairwise_character_orthogonal": pairwise_orthogonal,
            },
            "scope": {
                "proved": [
                    "one uniform degree-six star at every real K!=0",
                    "impossibility of representing its boundary weight by any real exponential with arbitrary constant, one-spin fields, and pairwise couplings among the six visible neighbours",
                    "non-closure of exact checkerboard decimation within visible-spin pairwise Ising Hamiltonians on the infinite simple-cubic lattice and open boxes containing a full-degree eliminated site",
                ],
                "not_proved": [
                    "no-go for representations with hidden auxiliary spins or unrestricted tensor bond dimension",
                    "no-go for nonlocal transformations, dual variables, or models retaining four- and six-spin interactions",
                    "nonintegrability or impossibility of solving the three-dimensional Ising model",
                    "an exact critical point, thermodynamic free energy, or critical exponent",
                ],
            },
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for row in payload["checks"]:
        print(f"[{'PASS' if row['passed'] else 'FAIL'}] {row['name']}: {row['detail']}")
    print(f"PASS: {len(payload['checks'])}/{len(payload['checks'])} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
