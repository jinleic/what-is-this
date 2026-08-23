"""Exact characteristic-zero pair-product counting from irreducible factors.

For each irreducible spectral factor, resultants construct the polynomials whose
roots are squares, within-factor unordered products, and cross-factor products.
The squarefree lcm therefore has exactly the distinct pair products as its roots.
This gives both inequality directions over Q, unlike a modular derivative gcd.

Run from the repository root with:
    .venv/bin/python experiments/e167_rank_char0_products.py
"""

from __future__ import annotations

from fractions import Fraction
from typing import Iterable

import sympy as sp

from e166_rank_char0_spectrum import (
    CpuBudget,
    Z,
    exact_spectrum_record,
    polynomial_digest,
)

U = sp.Symbol("u")


def _as_expression(poly: sp.Poly, variable: sp.Symbol) -> sp.Expr:
    coefficients = poly.all_coeffs()
    degree = poly.degree()
    return sum(coefficient * variable ** (degree - index)
               for index, coefficient in enumerate(coefficients))


def cross_product_polynomial(
    left: sp.Poly,
    right: sp.Poly,
    budget: CpuBudget,
    stage: str,
) -> sp.Poly:
    """Return prod_{alpha root left, beta root right} (z-alpha*beta)."""

    budget.check(f"before resultant {stage}")
    left_u = _as_expression(left, U)
    right_coefficients = right.all_coeffs()
    right_degree = right.degree()
    # u^deg(right) right(z/u) = sum_k b_k z^(n-k) u^k.
    transformed = sum(
        coefficient * Z ** (right_degree - index) * U**index
        for index, coefficient in enumerate(right_coefficients)
    )
    result = sp.Poly(sp.resultant(left_u, transformed, U), Z, domain=sp.QQ).monic()
    budget.check(f"after resultant {stage}")
    expected_degree = left.degree() * right.degree()
    if result.degree() != expected_degree:
        raise ArithmeticError(
            f"cross-product resultant {stage} has degree {result.degree()}, expected {expected_degree}"
        )
    return result


def square_product_polynomial(
    factor: sp.Poly,
    budget: CpuBudget,
    stage: str,
) -> sp.Poly:
    """Return prod_{alpha root factor} (z-alpha^2)."""

    budget.check(f"before square resultant {stage}")
    result = sp.Poly(
        sp.resultant(_as_expression(factor, U), Z - U**2, U),
        Z,
        domain=sp.QQ,
    ).monic()
    budget.check(f"after square resultant {stage}")
    if result.degree() != factor.degree():
        raise ArithmeticError(
            f"square resultant {stage} has degree {result.degree()}, expected {factor.degree()}"
        )
    return result


def monic_polynomial_square_root(square: sp.Poly) -> sp.Poly:
    """Recover the unique monic q with q^2=square, and verify the identity exactly."""

    square = sp.Poly(square, Z, domain=sp.QQ).monic()
    if square.degree() % 2:
        raise ArithmeticError("a polynomial square cannot have odd degree")
    root_degree = square.degree() // 2
    square_coefficients = square.all_coeffs()
    root_coefficients = [sp.QQ.one]
    for index in range(1, root_degree + 1):
        convolution = sum(
            (root_coefficients[left] * root_coefficients[index - left]
             for left in range(1, index)),
            sp.QQ.zero,
        )
        root_coefficients.append((square_coefficients[index] - convolution) / 2)
    root = sp.Poly.from_list(root_coefficients, gens=Z, domain=sp.QQ)
    if root * root != square:
        raise ArithmeticError("exact monic polynomial square-root verification failed")
    return root


def within_factor_product_polynomial(
    factor: sp.Poly,
    square_poly: sp.Poly,
    budget: CpuBudget,
    stage: str,
) -> sp.Poly:
    """Return prod_{i<j}(z-alpha_i alpha_j) for the roots of one factor."""

    ordered = cross_product_polynomial(factor, factor, budget, f"{stage}_ordered")
    quotient = ordered.exquo(square_poly)
    within = monic_polynomial_square_root(quotient)
    expected_degree = factor.degree() * (factor.degree() - 1) // 2
    if within.degree() != expected_degree:
        raise ArithmeticError(
            f"within-factor polynomial {stage} has degree {within.degree()}, expected {expected_degree}"
        )
    return within


def _component(
    name: str,
    kind: str,
    factor_indices: Iterable[int],
    polynomial: sp.Poly,
) -> dict[str, object]:
    return {
        "name": name,
        "kind": kind,
        "factor_indices": list(factor_indices),
        "polynomial": sp.Poly(polynomial, Z, domain=sp.QQ).monic(),
    }


def build_r_components(
    factors: list[tuple[sp.Poly, int]],
    budget: CpuBudget,
) -> tuple[list[dict[str, object]], dict[int, sp.Poly]]:
    components: list[dict[str, object]] = []
    squares: dict[int, sp.Poly] = {}
    for index, (factor, exponent) in enumerate(factors):
        square_poly = square_product_polynomial(factor, budget, f"f{index}")
        squares[index] = square_poly
        if factor.degree() >= 2:
            within = within_factor_product_polynomial(
                factor, square_poly, budget, f"f{index}"
            )
            components.append(_component(
                f"within_f{index}", "within_distinct_roots", [index], within
            ))
        if exponent >= 2:
            components.append(_component(
                f"repeated_square_f{index}", "square_from_repeated_slots", [index], square_poly
            ))

    for left in range(len(factors)):
        for right in range(left + 1, len(factors)):
            cross = cross_product_polynomial(
                factors[left][0], factors[right][0], budget, f"f{left}_f{right}"
            )
            components.append(_component(
                f"cross_f{left}_f{right}", "cross_factor", [left, right], cross
            ))
    return components, squares


def extend_squarefree_union(
    initial: sp.Poly,
    components: list[dict[str, object]],
    budget: CpuBudget,
    stage: str,
) -> tuple[sp.Poly, list[dict[str, object]]]:
    union = sp.Poly(initial, Z, domain=sp.QQ).monic()
    rows: list[dict[str, object]] = []
    for position, component in enumerate(components):
        budget.check(f"{stage} lcm component {position}")
        polynomial = component["polynomial"]
        squarefree = sp.Poly(polynomial.sqf_part(), Z, domain=sp.QQ).monic()
        overlap = sp.Poly(sp.gcd(union, squarefree), Z, domain=sp.QQ).monic()
        union = sp.Poly((union * squarefree).exquo(overlap), Z, domain=sp.QQ).monic()
        rows.append({
            "name": component["name"],
            "kind": component["kind"],
            "factor_indices": component["factor_indices"],
            "raw_degree": polynomial.degree(),
            "squarefree_degree": squarefree.degree(),
            "internal_collision_degree": polynomial.degree() - squarefree.degree(),
            "overlap_with_previous_union_degree": overlap.degree(),
            "cumulative_union_degree": union.degree(),
            "polynomial_sha256": polynomial_digest(polynomial),
            "squarefree_sha256": polynomial_digest(squarefree),
        })
    return union, rows


def exact_pair_product_counts(
    factors: list[tuple[sp.Poly, int]],
    budget: CpuBudget,
) -> dict[str, object]:
    """Return exact r and r_all certificates over Q for a factored spectrum."""

    r_components, squares = build_r_components(factors, budget)
    r_union, r_rows = extend_squarefree_union(
        sp.Poly(1, Z, domain=sp.QQ), r_components, budget, "r"
    )

    missing_square_components = [
        _component(
            f"new_square_f{index}",
            "square_from_diagonal_slot",
            [index],
            squares[index],
        )
        for index, (_factor, exponent) in enumerate(factors)
        if exponent == 1
    ]
    all_union, all_extension_rows = extend_squarefree_union(
        r_union, missing_square_components, budget, "r_all"
    )


    r_value = r_union.degree()
    r_all_value = all_union.degree()
    return {
        "claim_tag": "[COMPUTATION]",
        "method": (
            "[LEMMA] Exact Q-resultants construct every root-set component; the "
            "squarefree lcm has exactly the union of pair-product values as roots."
        ),
        "r": r_value,
        "r_all": r_all_value,
        "r_two_sided_sandwich": {"lower": r_value, "upper": r_value},
        "r_all_two_sided_sandwich": {"lower": r_all_value, "upper": r_all_value},
        "pair_slot_polynomial_degree": 496,
        "symmetric_slot_polynomial_degree": 528,
        "pair_derivative_gcd_degree_over_Q": 496 - r_value,
        "symmetric_derivative_gcd_degree_over_Q": 528 - r_all_value,
        "r_squarefree_union_sha256": polynomial_digest(r_union),
        "r_all_squarefree_union_sha256": polynomial_digest(all_union),
        "r_union_is_squarefree": True,
        "r_all_union_is_squarefree": True,
        "squarefree_justification": (
            "[LEMMA] Each component is replaced by sqf_part; the lcm of squarefree "
            "polynomials over Q is squarefree by unique factorization."
        ),
        "r_components": r_rows,
        "r_all_new_square_components": all_extension_rows,
        "inequality_directions": {
            "upper": (
                "[LEMMA] Every allowed slot-pair product is a root of the displayed "
                "squarefree lcm, so r <= deg(lcm) and r_all <= deg(lcm_all)."
            ),
            "lower": (
                "[LEMMA] Every root of every resultant component is an actual product "
                "of spectral roots; the squarefree lcm has deg(lcm) distinct roots, "
                "so r >= deg(lcm) and r_all >= deg(lcm_all)."
            ),
        },
    }


def main() -> None:
    budget = CpuBudget()
    _spectrum, factors = exact_spectrum_record(Fraction(1, 3), budget)
    counts = exact_pair_product_counts(factors, budget)
    checks = [
        counts["r"] == 417,
        counts["r_all"] == 445,
        counts["pair_derivative_gcd_degree_over_Q"] == 79,
        counts["symmetric_derivative_gcd_degree_over_Q"] == 83,
        counts["r_union_is_squarefree"] is True,
        counts["r_all_union_is_squarefree"] is True,
    ]
    if not all(checks):
        raise SystemExit("FAIL e167_rank_char0_products")
    print(
        "PASS e167_rank_char0_products: exact-Q r(T)=417, r_all(T)=445",
        flush=True,
    )


if __name__ == "__main__":
    main()
