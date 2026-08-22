"""Standalone exact and high-precision checks for the Kac--Ward investigation."""

from __future__ import annotations

from fractions import Fraction

import mpmath as mp
import sympy as sp

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    BULK_SC_LOG_P,
    bulk_cubic_trace_monomials,
    closed_nonbacktracking_walk_count,
    exact_kac_ward_polynomial,
    formal_log_coefficients,
    fit_strict_cubic_candidate,
    kac_ward_determinant,
    kac_ward_trace_basis,
    numeric_determinant_prefix,
    polynomial_square,
    strict_cubic_determinant_prefix,
    strict_cubic_equations,
)
from ising.lattices import cubic, square


DPS = 90
TOL = mp.mpf("1e-70")
SQUARE_CASES = ((3, 3), (3, 4), (4, 4), (4, 5))


def _polyval(coefficients: list[int] | tuple[int, ...], value: mp.mpf) -> mp.mpf:
    result = mp.mpf(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def _to_mpc(value: sp.Expr, dps: int) -> mp.mpc:
    evaluated = sp.N(value, dps)
    with mp.workdps(dps):
        return +mp.mpc(str(sp.re(evaluated)), str(sp.im(evaluated)))


def check_2d_control() -> None:
    with mp.workdps(DPS):
        value = mp.mpf("0.173205080756887729352744634150587236694280525381")
        for shape in SQUARE_CASES:
            lattice = square(*shape, periodic=False)
            target = polynomial_square(even_subgraph_polynomial(lattice))
            observed, certificate = exact_kac_ward_polynomial(lattice)
            discrepancies = [a - b for a, b in zip(observed, target, strict=True)]
            assert max(map(abs, discrepancies), default=0) == 0
            assert certificate.max_nonrational_component == 0
            assert certificate.modulus > 2 * certificate.coefficient_bound

            numerical = kac_ward_determinant(lattice, value, dps=DPS)
            expected = _polyval(target, value)
            error = abs(numerical - expected)
            assert error < TOL, (shape, mp.nstr(error, 12))

            for order in (4, 6, 8):
                basis = kac_ward_trace_basis(lattice, order)
                log_p = formal_log_coefficients(even_subgraph_polynomial(lattice), order)
                required = -2 * order * log_p[order]
                assert basis == (required, 0, 0, 0), (shape, order, basis, required)

            print(
                f"2D {shape[0]}x{shape[1]} free: exact max coefficient discrepancy 0; "
                f"{DPS}-dps determinant error {mp.nstr(error, 6)}"
            )


def check_3d_obstruction() -> None:
    a, b, equations, required = strict_cubic_equations(BULK_SC_LOG_P)
    assert bulk_cubic_trace_monomials(4) == {0: 24}
    assert bulk_cubic_trace_monomials(6) == {0: 192, 2: 72}
    assert bulk_cubic_trace_monomials(8) == {0: 1080, 1: 2304, 2: 1056, 4: 144}
    assert required == {4: Fraction(-24), 6: Fraction(-264), 8: Fraction(-3000)}
    basis = sp.groebner([equations[k] for k in (4, 6, 8)], a, b, order="lex")
    assert list(basis) == [1]
    assert list(sp.groebner([equations[k] for k in (4, 6)], a, b, order="lex")) != [1]

    lattice = cubic(3, 3, 2, periodic=False)
    target = polynomial_square(even_subgraph_polynomial(lattice))
    straight, turn = fit_strict_cubic_candidate(lattice)
    symbolic = strict_cubic_determinant_prefix(lattice, straight, turn, 8)
    symbolic_failures = [
        k for k in range(9) if sp.simplify(symbolic[k] - target[k]) != 0
    ]
    assert symbolic_failures[0] == 8

    numeric = numeric_determinant_prefix(
        lattice,
        _to_mpc(straight, DPS),
        _to_mpc(turn, DPS),
        8,
        dps=DPS,
    )
    numeric_failures = [
        k for k in range(9) if abs(numeric[k] - target[k]) > TOL
    ]
    assert numeric_failures[0] == 8
    assert abs(numeric[8] - _to_mpc(symbolic[8], DPS)) < TOL

    assert closed_nonbacktracking_walk_count(lattice, 4) == 160
    assert closed_nonbacktracking_walk_count(lattice, 6) == 936
    assert closed_nonbacktracking_walk_count(lattice, 8) == 6688
    print(
        "3D strict cubic scalar weights: k=4,6 solvable; k=4,6,8 Groebner basis [1]; "
        "independent finite-matrix first coefficient failure v^8"
    )


def main() -> None:
    check_2d_control()
    check_3d_obstruction()
    print("PASS")


if __name__ == "__main__":
    main()
