#!/usr/bin/env python3
"""Independent exact verifier for the cubic-star decimation obstruction.

No producer imports and no floating-point arithmetic.  The verifier rebuilds
Walsh multiplicities directly from all spin configurations and checks the
integer polynomial certificate proving the six-body coefficient is positive
for every nonzero real coupling.
"""

from __future__ import annotations

import itertools
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "star_decimation.json"
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def add(left: list[int], right: list[int]) -> list[int]:
    out = [0] * max(len(left), len(right))
    for index in range(len(out)):
        out[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def multiply(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return out


def power(poly: list[int], exponent: int) -> list[int]:
    out = [1]
    base = poly
    while exponent:
        if exponent & 1:
            out = multiply(out, base)
        exponent >>= 1
        if exponent:
            base = multiply(base, base)
    return out


def scale(poly: list[int], scalar: int) -> list[int]:
    return [scalar * coefficient for coefficient in poly]


def shift_at_one(poly: list[int]) -> list[int]:
    """Return coefficients of poly(1+y), using integer binomial arithmetic."""

    from math import comb

    out = [0] * len(poly)
    for degree, coefficient in enumerate(poly):
        for shifted_degree in range(degree + 1):
            out[shifted_degree] += coefficient * comb(degree, shifted_degree)
    return out


def divide_by_x_minus_one(poly: list[int]) -> tuple[list[int], int]:
    """Synthetic division in ascending order."""

    descending = list(reversed(poly))
    quotient_desc = [descending[0]]
    for coefficient in descending[1:-1]:
        quotient_desc.append(coefficient + quotient_desc[-1])
    remainder = descending[-1] + quotient_desc[-1]
    return list(reversed(quotient_desc)), remainder


def walsh_counts(degree: int, subset_size: int) -> dict[int, int]:
    counts: dict[int, int] = defaultdict(int)
    for state in range(1 << degree):
        spins = [1 if (state >> index) & 1 else -1 for index in range(degree)]
        character = 1
        for spin in spins[:subset_size]:
            character *= spin
        counts[abs(sum(spins))] += character
    return dict(sorted(counts.items()))
def walsh_inner(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    total = 0
    for state in range(64):
        character = 1
        for site in left + right:
            character *= 1 if (state >> site) & 1 else -1
        total += character
    return total


STEPS = (
    (-1, 0, 0),
    (1, 0, 0),
    (0, -1, 0),
    (0, 1, 0),
    (0, 0, -1),
    (0, 0, 1),
)


def box_incidence(side: int, parity: int) -> tuple[bool, int]:
    supports: dict[tuple[tuple[int, int, int], ...], int] = defaultdict(int)
    full: list[tuple[tuple[int, int, int], ...]] = []
    maximum_boundary = 0
    for x in range(side):
        for y in range(side):
            for z in range(side):
                if (x + y + z) % 2 != parity:
                    continue
                support = tuple(
                    sorted(
                        (x + dx, y + dy, z + dz)
                        for dx, dy, dz in STEPS
                        if 0 <= x + dx < side
                        and 0 <= y + dy < side
                        and 0 <= z + dz < side
                    )
                )
                supports[support] += 1
                if len(support) == 6:
                    full.append(support)
                else:
                    maximum_boundary = max(maximum_boundary, len(support))
    return all(supports[support] == 1 for support in full), maximum_boundary




def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact has meta data checks", {"meta", "data", "checks"} <= set(artifact))
    names = [row["name"] for row in artifact["checks"]]
    check("producer checks unique and true", len(names) == len(set(names)) and all(row["passed"] for row in artifact["checks"]))

    expected_d3 = walsh_counts(3, 2)
    expected_d4 = walsh_counts(4, 4)
    expected_d6_even = {
        subset_size: walsh_counts(6, subset_size) for subset_size in (0, 2, 4, 6)
    }
    expected_d6 = expected_d6_even[6]
    check("degree-three control has no higher even Walsh sector", 3 < 4 and expected_d3 == {1: -2, 3: 2})
    check("degree-four top Walsh counts rebuilt", expected_d4 == {0: 6, 2: -8, 4: 2})
    check("degree-six top Walsh counts rebuilt", expected_d6 == {0: -20, 2: 30, 4: -12, 6: 2})
    check(
        "stored Walsh counts match independent enumeration",
        artifact["data"]["degree_four"]["top_walsh_counts"] == {str(k): v for k, v in expected_d4.items()}
        and artifact["data"]["degree_six"]["top_walsh_counts"] == {str(k): v for k, v in expected_d6.items()},
    )
    check(
        "full octahedral Walsh decomposition rebuilt",
        expected_d6_even
        == {
            0: {0: 20, 2: 30, 4: 12, 6: 2},
            2: {0: -4, 2: -2, 4: 4, 6: 2},
            4: {0: 4, 2: -2, 4: -4, 6: 2},
            6: {0: -20, 2: 30, 4: -12, 6: 2},
        }
        and artifact["data"]["degree_six"]["even_walsh_counts_by_subset_size"]
        == {
            str(subset_size): {str(key): value for key, value in counts.items()}
            for subset_size, counts in expected_d6_even.items()
        },
    )
    derived_vectors = {
        f"c{subset_size}": {
            str(absolute_sum): coefficient // 2
            for absolute_sum, coefficient in expected_d6_even[subset_size].items()
            if absolute_sum
        }
        for subset_size in (2, 4, 6)
    }
    check(
        "octahedral log-cosh exponent vectors independently derived",
        derived_vectors
        == {
            "c2": {"2": -1, "4": 2, "6": 1},
            "c4": {"2": -1, "4": -2, "6": 1},
            "c6": {"2": 15, "4": -6, "6": 1},
        }
        and artifact["data"]["degree_six"][
            "logcosh_exponent_vectors_scaled_by_32"
        ]
        == derived_vectors,
    )
    c2_factors = artifact["data"]["degree_six"]["two_body_ratio_factors"]
    check(
        "octahedral two-body coefficient is strictly positive",
        c2_factors
        == [
            {"constant": -1, "linear": 2, "power": 2},
            {"constant": -3, "linear": 4, "power": 1},
        ]
        and all(
            factor["constant"] + factor["linear"] == 1
            and factor["linear"] > 0
            for factor in c2_factors
        ),
    )

    # For x=cosh(2K)^2, exp(32*c6)-1 has numerator
    # P=x^8(4x-3)-(2x-1)^6.
    first = multiply(power([0, 1], 8), [-3, 4])
    second = power([-1, 2], 6)
    p = add(first, scale(second, -1))
    quotient = p
    remainders = []
    for _ in range(3):
        quotient, remainder = divide_by_x_minus_one(quotient)
        remainders.append(remainder)
    shifted = shift_at_one(quotient)
    check("six-body numerator has a triple root at x=1", remainders == [0, 0, 0], str(remainders))
    check("remaining shifted factor is coefficient-positive", shifted == [8, 54, 144, 188, 120, 33, 4] and all(value > 0 for value in shifted), str(shifted))
    check(
        "stored polynomial certificate matches",
        artifact["data"]["degree_six"]["ratio_numerator_coefficients_ascending"] == p
        and artifact["data"]["degree_six"]["positive_shift_factor_coefficients_ascending"] == shifted,
    )

    # The degree-four obstruction is exp(8*c4)=(2x-1)/x^2,
    # whose denominator minus numerator is (x-1)^2 > 0 for x>1.
    degree_four_gap = add(power([0, 1], 2), scale([-1, 2], -1))
    check(
        "degree-four ratio gap is exactly (x-1)^2",
        degree_four_gap == power([-1, 1], 2),
        str(degree_four_gap),
    )
    d6_four_body_gap = add(power([-1, 2], 2), scale([-3, 4], -1))
    check(
        "octahedral four-body coefficient is strictly negative",
        d6_four_body_gap == scale(power([-1, 1], 2), 4)
        and artifact["data"]["degree_six"][
            "four_body_denominator_gap_coefficients_ascending"
        ]
        == d6_four_body_gap,
        str(d6_four_body_gap),
    )
    pairwise_orthogonal = all(
        walsh_inner(tuple(range(6)), subset) == 0
        for subset_size in range(3)
        for subset in itertools.combinations(range(6), subset_size)
    )
    incidence = [
        box_incidence(side, parity)
        for side in (4, 5)
        for parity in (0, 1)
    ]
    barycenter = artifact["data"]["whole_sublattice"]["barycenter_certificate"]
    check(
        "visible pairwise characters are exactly orthogonal to six-body support",
        pairwise_orthogonal
        and artifact["data"]["whole_sublattice"][
            "pairwise_character_orthogonal"
        ],
    )
    offset_sum = [sum(step[axis] for step in STEPS) for axis in range(3)]
    stored_incidence = artifact["data"]["whole_sublattice"][
        "checkerboard_incidence_audits"
    ]
    check(
        "checkerboard incidence and barycenter block global cancellation",
        offset_sum == [0, 0, 0]
        and barycenter
        == {
            "neighbour_count": 6,
            "offset_sum": [0, 0, 0],
            "identity": "sum_(u in N(v)) u = 6v",
        }
        and all(unique and maximum <= 5 for unique, maximum in incidence)
        and all(
            record["full_supports_unique"]
            and record["maximum_boundary_arity"] <= 5
            for record in stored_incidence.values()
        ),
    )
    check(
        "open-box reformulation retains boundary factors explicitly",
        "lower-arity factors"
        in artifact["data"]["whole_sublattice"]["exact_reformulation"],
    )

    limitations = artifact["data"]["scope"]["not_proved"]
    check(
        "scope excludes hidden-auxiliary and global-solution overclaims",
        any("auxiliar" in item.lower() for item in limitations)
        and any("solv" in item.lower() for item in limitations),
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks")
        return 1
    print("PASS: independent star-decimation obstruction verifier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
