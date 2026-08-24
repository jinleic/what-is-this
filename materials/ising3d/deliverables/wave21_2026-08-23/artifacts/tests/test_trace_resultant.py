#!/usr/bin/env python3
"""Independent exact verifier for the trace-resultant obstruction.

No producer imports.  This verifier reconstructs the rational transfer matrix
from the spin definition, computes traces with SymPy integer matrices rather
than NumPy object matrices, evaluates cyclotomic norms as quotient-ring
multiplication determinants rather than resultants, and uses lexicographic
Groebner elimination rather than the producer's grevlex order.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "trace_resultant.json"
N = 6
T = Fraction(1, 3)
BONDS_2X3 = ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5))
BONDS_CHAIN = tuple((index, index + 1) for index in range(N - 1))
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/")
    return Fraction(int(numerator), int(denominator))


def build_integer_matrix(bonds: tuple[tuple[int, int], ...]) -> tuple[int, sp.Matrix, Fraction, int]:
    q = (1 + T * T) / (2 * T)
    dimension = 1 << N
    energies: list[int] = []
    for state in range(dimension):
        spins = [1 - 2 * ((state >> (N - 1 - site)) & 1) for site in range(N)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    epsilon = energies[0] % 2
    if any((energy - epsilon) % 2 for energy in energies):
        raise AssertionError("energy parity mismatch")
    diagonal = [q ** ((energy - epsilon) // 2) for energy in energies]
    powers = [T**exponent for exponent in range(N + 1)]
    p_matrix = [
        [powers[(row ^ column).bit_count()] for column in range(dimension)]
        for row in range(dimension)
    ]
    rational = [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]
    for row in range(dimension):
        for column in range(dimension):
            rational[row][column] = sum(
                p_matrix[row][middle] * diagonal[middle] * p_matrix[middle][column]
                for middle in range(dimension)
            )
    scale = 1
    for row in rational:
        for value in row:
            scale = lcm(scale, value.denominator)
    integer = sp.Matrix([[int(value * scale) for value in row] for row in rational])
    return scale, integer, q, epsilon


def exact_traces(integer: sp.Matrix) -> list[int]:
    power = sp.eye(integer.rows)
    traces: list[int] = []
    for _ in range(7):
        power *= integer
        traces.append(int(sp.trace(power)))
    return traces


def normalized_targets(
    scale: int, traces_integer: list[int], q: Fraction, epsilon: int
) -> dict[str, Fraction]:
    traces = [
        Fraction(traces_integer[power - 1], scale**power) for power in range(1, 8)
    ]
    centre_square = (1 - T * T) ** (2 * N) / q**epsilon
    return {
        "r1_squared": traces[0] ** 2 / centre_square,
        "r2": traces[1] / centre_square,
        "r3_over_r1": traces[2] / (traces[0] * centre_square),
        "r4": traces[3] / centre_square**2,
        "r5_over_r1": traces[4] / (traces[0] * centre_square**2),
        "r6_over_r2": traces[5] / (traces[1] * centre_square**2),
        "r7_over_r1": traces[6] / (traces[0] * centre_square**3),
    }


def quotient_norm(polynomial: sp.Expr, minimal: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    """Norm of polynomial modulo a monic minimal polynomial via multiplication."""

    degree = sp.Poly(minimal, variable).degree()
    columns: list[list[sp.Expr]] = []
    for basis_power in range(degree):
        remainder = sp.Poly(
            sp.rem(sp.expand(polynomial * variable**basis_power), minimal, variable),
            variable,
        )
        columns.append([remainder.coeff_monomial(variable**row) for row in range(degree)])
    multiplication = sp.Matrix(degree, degree, lambda row, column: columns[column][row])
    return sp.factor(multiplication.det())


def equations(targets: dict[str, Fraction]) -> list[sp.Poly]:
    a5, a4, a3 = sp.symbols("a5 a4 a3")
    u = sp.symbols("u")
    rat = {name: sp.Rational(value.numerator, value.denominator) for name, value in targets.items()}
    c2 = rat["r2"] - rat["r1_squared"] - 64 - 32 * a5 - 16 * a4 - 8 * a3
    c3 = rat["r3_over_r1"] - rat["r1_squared"] - 729 - 243 * a5 - 81 * a4 - 27 * a3
    a2 = c3 / 3 - c2 / 2
    a1 = 3 * c2 / 2 - 2 * c3 / 3
    qpoly = sp.expand(
        u**6
        + a5 * u**5
        + a4 * u**4
        + a3 * u**3
        + a2 * u**2
        + a1 * u
        + rat["r1_squared"]
    )
    minimal_polynomials = (
        ("r4", u**2 - 4 * u + 2),
        ("r5_over_r1", u**2 - 5 * u + 5),
        ("r6_over_r2", u**2 - 4 * u + 1),
        ("r7_over_r1", u**3 - 7 * u**2 + 14 * u - 7),
    )
    return [
        sp.Poly(quotient_norm(qpoly, minimal, u) - rat[name], a5, a4, a3, domain=sp.QQ)
        for name, minimal in minimal_polynomials
    ]


def prefix_units(targets: dict[str, Fraction]) -> list[bool]:
    a5, a4, a3 = sp.symbols("a5 a4 a3")
    polys = equations(targets)
    results: list[bool] = []
    for count in range(1, 5):
        basis = sp.groebner(
            [poly.as_expr() for poly in polys[:count]],
            a5,
            a4,
            a3,
            order="lex",
            domain=sp.QQ,
        )
        results.append([poly.as_expr() for poly in basis.polys] == [sp.Integer(1)])
    return results


def lucas_audit() -> bool:
    tau, u = sp.symbols("tau u")
    lucas = [sp.Integer(2), tau]
    for _ in range(2, 8):
        lucas.append(sp.expand(tau * lucas[-1] - lucas[-2]))
    expected = {
        2: u - 2,
        3: tau * (u - 3),
        4: u**2 - 4 * u + 2,
        5: tau * (u**2 - 5 * u + 5),
        6: (u - 2) * (u**2 - 4 * u + 1),
        7: tau * (u**3 - 7 * u**2 + 14 * u - 7),
    }
    return all(
        sp.expand(lucas[index] - expression.subs(u, tau**2)) == 0
        for index, expression in expected.items()
    )


def synthetic_audit() -> bool:
    u = sp.symbols("u")
    mode_u = [sp.Rational((prime + 1) ** 2, prime) for prime in (2, 3, 5, 7, 11, 13)]
    qpoly = sp.expand(sp.prod(u - value for value in mode_u))
    targets = {
        "r1_squared": Fraction(sp.prod(mode_u)),
        "r2": Fraction(qpoly.subs(u, 2)),
        "r3_over_r1": Fraction(qpoly.subs(u, 3)),
    }
    for name, minimal in (
        ("r4", u**2 - 4 * u + 2),
        ("r5_over_r1", u**2 - 5 * u + 5),
        ("r6_over_r2", u**2 - 4 * u + 1),
        ("r7_over_r1", u**3 - 7 * u**2 + 14 * u - 7),
    ):
        targets[name] = Fraction(quotient_norm(qpoly, minimal, u))
    coefficients = sp.Poly(qpoly, u).all_coeffs()
    substitution = dict(zip(sp.symbols("a5 a4 a3"), coefficients[1:4], strict=True))
    return all(poly.as_expr().subs(substitution) == 0 for poly in equations(targets))


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact has meta/data/checks", set(("meta", "data", "checks")) <= set(artifact))
    names = [record["name"] for record in artifact["checks"]]
    check("artifact checks are nonempty and uniquely named", bool(names) and len(names) == len(set(names)))
    check("artifact producer checks all passed", all(record["passed"] for record in artifact["checks"]))
    check("Lucas polynomials through order seven have declared factors", lucas_audit())
    check("known rational six-mode control satisfies all four norm equations", synthetic_audit())

    scale, integer, q, epsilon = build_integer_matrix(BONDS_2X3)
    traces = exact_traces(integer)
    targets = normalized_targets(scale, traces, q, epsilon)
    stored = artifact["data"]["trace_record"]
    check("fresh open-2x3 matrix scale matches", scale == int(stored["matrix_scale"]), str(scale))
    check("fresh open-2x3 traces one through seven match", traces == stored["integer_traces"])
    check(
        "fresh normalized trace targets match",
        targets == {name: parse_fraction(value) for name, value in stored["targets"].items()},
    )
    units = prefix_units(targets)
    check("trace-four through trace-six ideals remain proper", units[:3] == [False, False, False], str(units))
    check("adding trace seven gives exact unit ideal", units == [False, False, False, True], str(units))

    chain_scale, chain_integer, chain_q, chain_epsilon = build_integer_matrix(BONDS_CHAIN)
    chain_targets = normalized_targets(
        chain_scale, exact_traces(chain_integer), chain_q, chain_epsilon
    )
    check("open-chain control ideal remains proper through trace seven", not prefix_units(chain_targets)[-1])
    check("theorem scope remains one graph and one coupling", artifact["data"]["theorem"]["scope"].startswith("one finite bipartite layer"))
    check("benchmark was not used", artifact["data"]["theorem"]["benchmark_used_for_selection"] is False)

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} check(s): {', '.join(FAILURES)}")
        return 1
    print("PASS: independent trace-resultant verification complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
