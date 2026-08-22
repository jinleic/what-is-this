"""Standalone exact checks for e48's parametric Gaussianity certificates.

This file intentionally does not import the experiment.  It reconstructs the polynomial
representative independently, uses full 64x64 Fraction LDL rather than symmetry sectors, and
re-derives two propagated subinterval certificates from scratch.
"""

from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def independent_polynomial_matrix(
    n: int, bonds: list[tuple[int, int]]
) -> list[list[list[int]]]:
    """Expand (2t)^d q^s P diag(q^e) P term by term."""
    dim = 1 << n
    energies = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[i] * spins[j] for i, j in bonds))
    parity = energies[0] % 2
    exponents = [(energy - parity) // 2 for energy in energies]
    shift = -min(exponents)
    degree = max(exponent + shift for exponent in exponents)
    coefficients = [
        [[0] * (2 * n + 2 * degree + 1) for _ in range(dim)] for _ in range(dim)
    ]
    for row in range(dim):
        for col in range(dim):
            entry = coefficients[row][col]
            for state, exponent in enumerate(exponents):
                q_power = exponent + shift
                base = (
                    degree
                    - q_power
                    + (row ^ state).bit_count()
                    + (state ^ col).bit_count()
                )
                for choice in range(q_power + 1):
                    entry[base + 2 * choice] += (
                        2 ** (degree - q_power) * math.comb(q_power, choice)
                    )
    return coefficients


def evaluate(
    coefficients: list[list[list[int]]], t: Fraction
) -> list[list[Fraction]]:
    powers = [t**power for power in range(len(coefficients[0][0]))]
    return [
        [
            sum(
                (Fraction(coefficient) * powers[power]
                 for power, coefficient in enumerate(entry) if coefficient),
                Fraction(0),
            )
            for entry in row
        ]
        for row in coefficients
    ]


def fraction_inertia(matrix: list[list[Fraction]], shift: Fraction) -> int | None:
    """Independent dense symmetric Fraction LDL with no sector decomposition."""
    size = len(matrix)
    work = [
        [matrix[row][col] - (shift if row == col else 0) for col in range(size)]
        for row in range(size)
    ]
    negatives = 0
    for pivot_index in range(size):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot < 0:
            negatives += 1
        pivot_row = work[pivot_index]
        for row_index in range(pivot_index + 1, size):
            factor = work[row_index][pivot_index] / pivot
            if factor:
                row = work[row_index]
                for col_index in range(pivot_index, size):
                    row[col_index] -= factor * pivot_row[col_index]
    return negatives


def derivative_square_upper(
    coefficients: list[list[list[int]]], upper_t: Fraction
) -> Fraction:
    total = Fraction(0)
    for row in coefficients:
        for entry in row:
            derivative = sum(
                (Fraction(power * coefficient) * upper_t ** (power - 1)
                 for power, coefficient in enumerate(entry) if power and coefficient),
                Fraction(0),
            )
            total += derivative * derivative
    return total


def verify_row_from_scratch(
    row: dict[str, object], coefficients: list[list[list[int]]]
) -> dict[str, object]:
    center = Fraction(row["center_t"])
    half_width = Fraction(row["half_width_t"])
    matrix = evaluate(coefficients, center)

    stored_lambdas = row["lambda_center_enclosures"]
    enclosures = [
        (
            Fraction(stored_lambdas[f"lambda{index}"]["lower"]),
            Fraction(stored_lambdas[f"lambda{index}"]["upper"]),
        )
        for index in range(3)
    ]
    endpoint_counts = [
        (fraction_inertia(matrix, lower), fraction_inertia(matrix, upper))
        for lower, upper in enclosures
    ]

    stored_derivative_square = Fraction(row["frobenius_derivative_square_upper"])
    rebuilt_derivative_square = derivative_square_upper(
        coefficients, center + half_width
    )
    derivative_bound = Fraction(row["frobenius_derivative_upper"]["exact"])
    rho = Fraction(row["weyl_radius"]["exact"])
    (l0, u0), (l1, u1), (l2, u2) = enclosures
    global_lower = (l1 - rho) * (l2 - rho) / (u0 + rho)
    global_upper = (u1 + rho) * (u2 + rho) / (l0 - rho)
    tested = row["tested_center_absence_window"]
    tested_lower, tested_upper = Fraction(tested["lower"]), Fraction(tested["upper"])
    decisive_counts = (
        fraction_inertia(matrix, tested_lower),
        fraction_inertia(matrix, tested_upper),
    )

    return {
        "lambda_counts": endpoint_counts,
        "derivative_square_matches": rebuilt_derivative_square == stored_derivative_square,
        "sqrt_outward": derivative_bound * derivative_bound >= rebuilt_derivative_square,
        "rho_matches": rho == derivative_bound * half_width,
        "positive": l0 - rho > 0,
        "distinct": u0 + rho < l1 - rho and u1 + rho < l2 - rho,
        "forced_matches": (
            global_lower
            == Fraction(row["forced_value_whole_subinterval"]["lower"])
            and global_upper
            == Fraction(row["forced_value_whole_subinterval"]["upper"])
        ),
        "tested_contains_weyl_expansion": (
            tested_lower <= global_lower - rho
            and tested_upper >= global_upper + rho
        ),
        "decisive_counts": decisive_counts,
    }


def build_chain_R(n: int, t: Fraction) -> list[list[Fraction]]:
    """Direct original R construction, independent of polynomial scaling."""
    bonds = [(site, site + 1) for site in range(n - 1)]
    q = (1 + t * t) / (2 * t)
    energies = []
    for state in range(1 << n):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[i] * spins[j] for i, j in bonds))
    parity = energies[0] % 2
    diagonal = [q ** ((energy - parity) // 2) for energy in energies]
    powers = [t**power for power in range(n + 1)]
    p_matrix = [
        [powers[(row ^ col).bit_count()] for col in range(1 << n)]
        for row in range(1 << n)
    ]
    return [
        [
            sum(
                (p_matrix[row][middle] * diagonal[middle] * p_matrix[middle][col]
                 for middle in range(1 << n)),
                Fraction(0),
            )
            for col in range(1 << n)
        ]
        for row in range(1 << n)
    ]


print("1. stored artifact schema and embedded checks")
result_path = ROOT / "results" / "spectral" / "parametric_gaussianity.json"
with result_path.open() as handle:
    stored = json.load(handle)
check("top-level envelope", set(stored) == {"provenance", "data", "checks"})
provenance = stored["provenance"]
check(
    "required provenance",
    {"script", "generated_utc", "interpreter", "method"} <= set(provenance),
)
check("interpreter recorded", provenance["interpreter"] == ".venv/bin/python")
check("status is honestly unresolved", stored["data"]["status"] == "UNRESOLVED")
check("all embedded checks pass", all(item["passed"] for item in stored["checks"]))
rows = stored["data"]["certified_subintervals"]
check("four partial intervals stored", len(rows) == 4)
check(
    "each interval has nonzero exact width",
    all(Fraction(row["interval_t"]["lower"]) < Fraction(row["interval_t"]["upper"]) for row in rows),
)
check(
    "partial intervals do not masquerade as a chain",
    all(
        Fraction(rows[index]["interval_t"]["upper"])
        < Fraction(rows[index + 1]["interval_t"]["lower"])
        for index in range(3)
    ),
)

print("\n2. independent dense exact reconstruction of two subintervals")
bonds = list(layer_bonds((2, 3), (False, False)))
coefficients = independent_polynomial_matrix(6, bonds)
check(
    "independent coefficient positivity",
    all(coefficient >= 0 for row in coefficients for entry in row for coefficient in entry),
)
for row in (rows[0], rows[1]):
    rebuilt = verify_row_from_scratch(row, coefficients)
    name = row["name"]
    check(
        f"{name}: lambda enclosure inertia",
        rebuilt["lambda_counts"] == [(0, 1), (1, 2), (2, 3)],
        str(rebuilt["lambda_counts"]),
    )
    check(f"{name}: derivative polynomial bound", rebuilt["derivative_square_matches"])
    check(f"{name}: outward rational square root", rebuilt["sqrt_outward"])
    check(f"{name}: exact Weyl radius", rebuilt["rho_matches"])
    check(f"{name}: lambda0 positive throughout", rebuilt["positive"])
    check(f"{name}: three lowest distinct throughout", rebuilt["distinct"])
    check(f"{name}: forced interval recomputed", rebuilt["forced_matches"])
    check(
        f"{name}: tested window contains Weyl expansion",
        rebuilt["tested_contains_weyl_expansion"],
    )
    expected_count = tuple(row["inertia_counts_at_tested_window_ends"])
    check(
        f"{name}: dense Fraction LDL decisive counts",
        rebuilt["decisive_counts"] == expected_count,
        str(rebuilt["decisive_counts"]),
    )

print("\n3. independent chain control")
control = stored["data"]["control"]
chain_matrix = build_chain_R(4, Fraction(control["t"]))
forced_window = control["forced_window"]
chain_t = Fraction(control["t"])
chain_q = (1 + chain_t * chain_t) / (2 * chain_t)
chain_polynomial_scale = (2 * chain_t) ** 3 * chain_q**2
chain_counts = (
    fraction_inertia(
        chain_matrix, Fraction(forced_window["lower"]) / chain_polynomial_scale
    ),
    fraction_inertia(
        chain_matrix, Fraction(forced_window["upper"]) / chain_polynomial_scale
    ),
)
check(
    "chain control has an eigenvalue in the forced window",
    chain_counts == tuple(control["inertia_counts"]) == (5, 6),
    str(chain_counts),
)
check("control does not certify absence", chain_counts[0] != chain_counts[1])

print()
if FAILURES:
    print(f"FAIL: {FAILURES}")
    raise SystemExit(1)
print("PASS")
