"""Exact simultaneous Dolan--Grady and tridiagonal deformation search.

The search uses the smallest common translation/D4-orbit basis containing the
wave-4 F3 cancellation and the Ising bond generator.  Both generators range
independently over that basis in the projective chart with coefficient(X in A)
= coefficient(ZZ_edge in B) = 1.  All arithmetic is exact over Q.

Run from the repository root with
``.venv/bin/python experiments/e36_dg_simultaneous.py``.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import sympy as sp

EXPERIMENT_DIR = Path(__file__).resolve().parent
REPO_ROOT = EXPERIMENT_DIR.parent
RESULT_PATH = REPO_ROOT / "results" / "deformation" / "dg_simultaneous.json"
sys.path.insert(0, str(EXPERIMENT_DIR))

from e17_dg_deformation import (  # noqa: E402
    add,
    base_operators,
    body_weight,
    comm,
    compress_invariant,
    control_cases,
    dg2_residual,
    family_f3,
    format_pauli,
    mul,
    normalized_polynomials,
    orbit_operator,
    pack_word,
    simplify_symbolic_operator,
    smul,
    spatial_group,
    word_orbits,
)

Operator = dict[int, object]


def _equal(left: object, right: object) -> bool:
    return sp.expand(sp.sympify(left) - sp.sympify(right)) == 0


def _coefficient(operator: Operator, n: int, word: tuple[tuple[int, str], ...]) -> object:
    return sp.expand(operator.get(pack_word(word, n), 0))


def wave4_controls() -> dict[str, object]:
    """Reproduce all four frozen controls, raising before ansatz construction."""

    controls = control_cases()
    cycle = controls["cycle_C6"]
    obstruction = controls["path_P6_plus_bond_1_4"]
    assert cycle["dg2_floating_lambda"] == {
        "pauli_equations": 12,
        "unknowns": 1,
        "rank": 1,
        "nullity": 0,
        "solution_exists": True,
        "lambda": "1",
        "quartic_obstruction_terms": 0,
    }
    assert cycle["residual_at_solution_terms"] == 0
    assert cycle["tridiagonal_relation"]["pauli_equations"] == 96
    assert cycle["tridiagonal_relation"]["solution"] == {
        "beta": "2",
        "gamma": "0",
        "rho": "16",
    }
    assert obstruction["dg2_floating_lambda"] == {
        "pauli_equations": 14,
        "unknowns": 1,
        "rank": 1,
        "nullity": 0,
        "solution_exists": False,
        "lambda": None,
        "quartic_obstruction_terms": 2,
    }
    assert obstruction["residual_at_solution_terms"] is None
    assert obstruction["tridiagonal_relation"]["pauli_equations"] == 84
    assert obstruction["tridiagonal_relation"]["solution_exists"] is False

    group = spatial_group(3, 3, True)
    f3_orbits = word_orbits(family_f3(3, 3, True), group)
    f3_orbit = f3_orbits[69]
    assert len(f3_orbit) == 18
    assert format_pauli(pack_word(f3_orbit[0], 9), 9) == "X0 Z1 Z2"
    f3_direction = orbit_operator(f3_orbit, 9)
    n, A, B = base_operators(3, 3, True)
    A_f3 = add(A, smul(Fraction(-1, 2), f3_direction))
    f3_dg2 = dg2_residual(A_f3, B, n, Fraction(1))
    assert not f3_dg2

    first = comm(A_f3, B, n)
    third = comm(A_f3, comm(A_f3, first, n), n)
    witness_word = ((0, "Y"), (1, "Z"), (2, "X"), (3, "X"))
    witness = pack_word(witness_word, n)
    assert third[witness] == 3
    assert first.get(witness, 0) == 0

    return {
        "claim_tag": "[COMPUTATION]",
        "normalization_note": (
            "The frozen control writes ad_B^3(A)=16*lambda_frozen*ad_B(A); "
            "lambda_frozen=1 is lambda_direct=16 in this experiment's relations."
        ),
        "cycle_C6": cycle,
        "path_P6_plus_bond_1_4": obstruction,
        "torus_3x3_F3_DG2_cancellation": {
            "formula": "sum_v X_v [1 - (Z_{v+x}Z_{v-x} + Z_{v+y}Z_{v-y})/2]",
            "orbit_index_in_frozen_F3": 69,
            "orbit_representative": "X0 Z1 Z2",
            "orbit_size": len(f3_orbit),
            "lambda_frozen": "1",
            "lambda_direct": "16",
            "dg2_residual_terms": len(f3_dg2),
        },
        "torus_3x3_F3_DG1_witness": {
            "pauli_word": format_pauli(witness, n),
            "coefficient_in_ad_A_cubed_B": str(third[witness]),
            "coefficient_in_ad_A_B": str(first.get(witness, 0)),
        },
    }


def transverse_f3_orbit(a: int, b: int) -> Operator:
    """O_C=sum_v X_v(Z_{v+x}Z_{v-x}+Z_{v+y}Z_{v-y}) on a torus."""

    if a < 3 or b < 3:
        raise ValueError("the transverse-neighbour orbit requires both periods >= 3")
    n = a * b
    out: Operator = {}
    for x in range(a):
        for y in range(b):
            v = x * b + y
            opposite_pairs = (
                (((x + 1) % a) * b + y, ((x - 1) % a) * b + y),
                (x * b + (y + 1) % b, x * b + (y - 1) % b),
            )
            for u, w in opposite_pairs:
                vector = pack_word(((v, "X"), (u, "Z"), (w, "Z")), n)
                out[vector] = out.get(vector, 0) + 1
    return out


def minimal_basis(a: int, b: int) -> tuple[int, tuple[Operator, Operator, Operator]]:
    n, O_X, O_E = base_operators(a, b, True)
    O_C = transverse_f3_orbit(a, b)
    assert not (set(O_X) & set(O_E) or set(O_X) & set(O_C) or set(O_E) & set(O_C))
    return n, (O_X, O_E, O_C)


def basis_enumeration(a: int, b: int) -> dict[str, object]:
    n, basis = minimal_basis(a, b)
    group = spatial_group(a, b, True)
    names = ("O_X", "O_E", "O_C")
    formulas = (
        "sum_v X_v",
        "sum_{<uv>} Z_u Z_v over nearest-neighbour torus edges",
        "sum_v X_v (Z_{v+x}Z_{v-x}+Z_{v+y}Z_{v-y})",
    )
    rows = []
    cache: dict[int, int] = {}
    for name, formula, operator in zip(names, formulas, basis):
        compressed = compress_invariant(operator, group, n, cache)
        assert len(compressed) == 1
        representative, coefficient = next(iter(compressed.items()))
        assert coefficient == 1
        rows.append(
            {
                "name": name,
                "formula": formula,
                "representative": format_pauli(representative, n),
                "literal_pauli_terms": len(operator),
                "symmetry_orbits": len(compressed),
            }
        )
    return {
        "layer": f"{a}x{b}_torus",
        "sites": n,
        "group": "translations semidirect D4",
        "group_order": len(group),
        "orbits": rows,
    }


def chart_pair(
    a: int,
    b: int,
    alpha: object,
    chi: object,
    delta: object,
    eta: object,
) -> tuple[int, Operator, Operator]:
    n, (O_X, O_E, O_C) = minimal_basis(a, b)
    A = add(O_X, smul(alpha, O_E), smul(chi, O_C))
    B = add(smul(delta, O_X), O_E, smul(eta, O_C))
    return n, A, B


def dg_symbolic_system(a: int = 3, b: int = 3) -> dict[str, object]:
    alpha, chi, delta, eta, mu, lam = sp.symbols("alpha chi delta eta mu lambda")
    n, A, B = chart_pair(a, b, alpha, chi, delta, eta)
    first_a = comm(A, B, n)
    first_b = comm(B, A, n)
    third_a = comm(A, comm(A, first_a, n), n)
    third_b = comm(B, comm(B, first_b, n), n)
    r1 = simplify_symbolic_operator(add(third_a, smul(-mu, first_a)))
    r2 = simplify_symbolic_operator(add(third_b, smul(-lam, first_b)))
    group = spatial_group(a, b, True)
    cache: dict[int, int] = {}
    r1_reduced = compress_invariant(r1, group, n, cache)
    r2_reduced = compress_invariant(r2, group, n, cache)
    variables = (alpha, chi, delta, eta, mu, lam)
    return {
        "n": n,
        "variables": variables,
        "A": A,
        "B": B,
        "first_a": first_a,
        "first_b": first_b,
        "third_a": third_a,
        "third_b": third_b,
        "r1": r1,
        "r2": r2,
        "r1_reduced": r1_reduced,
        "r2_reduced": r2_reduced,
        "equations": normalized_polynomials(
            [*r1_reduced.values(), *r2_reduced.values()], variables
        ),
    }


def eliminate_scalar_multiple(
    third: Operator,
    first: Operator,
    deformation_variables: tuple[sp.Symbol, ...],
) -> list[sp.Expr]:
    """Eliminate k from third=k*first by all exact 2-by-2 minors."""

    keys = sorted(set(third) | set(first))
    rows = [(sp.expand(third.get(key, 0)), sp.expand(first.get(key, 0))) for key in keys]
    conditions = []
    for index, (left, right) in enumerate(rows):
        for other_left, other_right in rows[index + 1 :]:
            minor = sp.expand(left * other_right - other_left * right)
            if minor:
                conditions.append(minor)
    return normalized_polynomials(conditions, deformation_variables)


def support_counts(operator: Operator, n: int) -> dict[str, int]:
    return {
        str(weight): count
        for weight, count in sorted(Counter(body_weight(vector, n) for vector in operator).items())
    }


def dg_analysis() -> dict[str, object]:
    system = dg_symbolic_system()
    n = system["n"]
    alpha, chi, delta, eta, mu, lam = system["variables"]
    deformation_variables = (alpha, chi, delta, eta)
    p = alpha * delta - 1
    q = chi * delta - eta
    r = alpha * eta - chi

    witness_q_word = ((0, "Y"), (1, "Z"), (2, "Z"), (3, "Z"), (6, "X"))
    witness_q = _coefficient(system["r2_reduced"], n, witness_q_word)
    assert _equal(witness_q, 24 * q)
    witness_chi_p_word = ((0, "X"), (1, "Y"), (2, "X"), (3, "Z"))
    witness_chi_p = _coefficient(system["r1_reduced"], n, witness_chi_p_word)
    witness_chi_p_reduced = sp.expand(witness_chi_p.subs(eta, chi * delta))
    assert _equal(witness_chi_p_reduced, 8 * chi * p)
    witness_p_word = ((0, "Y"), (1, "Z"), (2, "Z"), (3, "Z"))
    witness_p = _coefficient(system["r2_reduced"], n, witness_p_word)
    witness_p_reduced = sp.expand(witness_p.subs({chi: 0, eta: 0}))
    assert _equal(witness_p_reduced, 48 * p)

    group = spatial_group(3, 3, True)
    cache: dict[int, int] = {}
    first_a = compress_invariant(system["first_a"], group, n, cache)
    first_b = compress_invariant(system["first_b"], group, n, cache)
    third_a = compress_invariant(system["third_a"], group, n, cache)
    third_b = compress_invariant(system["third_b"], group, n, cache)
    eliminated_dg1 = eliminate_scalar_multiple(third_a, first_a, deformation_variables)
    eliminated_dg2 = eliminate_scalar_multiple(third_b, first_b, deformation_variables)
    eliminated = normalized_polynomials(
        [*eliminated_dg1, *eliminated_dg2], deformation_variables
    )
    groebner = sp.groebner(eliminated, *deformation_variables, order="grevlex")
    groebner_expressions = [sp.factor(poly.as_expr()) for poly in groebner.polys]
    expected = [p**2, p * q, q**2, p * r, q * r, r**2]
    assert {sp.expand(poly.monic().as_expr()) for poly in groebner.polys} == {
        sp.Poly(expr, *deformation_variables, domain=sp.QQ).monic().as_expr()
        for expr in expected
    }

    saturation = {}
    z = sp.Symbol("z")
    for name, cross_product in (("p_nonzero", p), ("q_nonzero", q), ("r_nonzero", r)):
        saturated = sp.groebner(
            [*eliminated, z * cross_product - 1],
            z,
            *deformation_variables,
            order="lex",
        )
        basis = [sp.factor(poly.as_expr()) for poly in saturated.polys]
        assert basis == [sp.Integer(1)]
        saturation[name] = [str(value) for value in basis]

    r1_distinct = normalized_polynomials(system["r1_reduced"].values(), system["variables"])
    r2_distinct = normalized_polynomials(system["r2_reduced"].values(), system["variables"])
    return {
        "claim_tag": "[COMPUTATION]",
        "relations": {
            "DG1": "ad_A^3(B)=mu*ad_A(B)",
            "DG2": "ad_B^3(A)=lambda*ad_B(A)",
            "constant_normalization": "mu and lambda are direct constants (Dolan--Grady point 16)",
        },
        "unknowns": {
            "deformation": [str(value) for value in deformation_variables],
            "linear_relation_constants_eliminated_first": [str(mu), str(lam)],
            "total": len(system["variables"]),
        },
        "orbit_reduction": {
            "DG1": {
                "full_pauli_equations": len(system["r1"]),
                "symmetry_orbit_equations": len(system["r1_reduced"]),
                "distinct_up_to_nonzero_rational_scale": len(r1_distinct),
                "full_support_counts": support_counts(system["r1"], n),
                "orbit_support_counts": support_counts(system["r1_reduced"], n),
            },
            "DG2": {
                "full_pauli_equations": len(system["r2"]),
                "symmetry_orbit_equations": len(system["r2_reduced"]),
                "distinct_up_to_nonzero_rational_scale": len(r2_distinct),
                "full_support_counts": support_counts(system["r2"], n),
                "orbit_support_counts": support_counts(system["r2_reduced"], n),
            },
            "row_slots_both_relations": len(system["r1"]) + len(system["r2"]),
            "distinct_full_pauli_words_union": len(set(system["r1"]) | set(system["r2"])),
            "orbit_row_slots_both_relations": len(system["r1_reduced"]) + len(system["r2_reduced"]),
            "distinct_orbit_representatives_union": len(set(system["r1_reduced"]) | set(system["r2_reduced"])),
            "distinct_polynomials_both_relations": len(system["equations"]),
        },
        "cross_products": {
            "p": "alpha*delta - 1",
            "q": "chi*delta - eta",
            "r": "alpha*eta - chi",
            "nonproportional_condition": "at least one of p,q,r is nonzero",
        },
        "low_support_certificate": [
            {
                "stage": "unsubstituted DG2",
                "pauli_word": format_pauli(pack_word(witness_q_word, n), n),
                "support": 5,
                "coefficient": str(sp.factor(witness_q)),
                "consequence": "q=0",
            },
            {
                "stage": "DG1 after q=0 (eta=chi*delta)",
                "pauli_word": format_pauli(pack_word(witness_chi_p_word, n), n),
                "support": 4,
                "coefficient": str(sp.factor(witness_chi_p_reduced)),
                "consequence": "chi*p=0",
            },
            {
                "stage": "DG2 branch chi=eta=0",
                "pauli_word": format_pauli(pack_word(witness_p_word, n), n),
                "support": 4,
                "coefficient": str(sp.factor(witness_p_reduced)),
                "consequence": "p=0",
            },
        ],
        "linear_elimination": {
            "method": "all 2x2 minors of (third_commutator, first_commutator) after orbit reduction",
            "DG1_polynomials": len(eliminated_dg1),
            "DG2_polynomials": len(eliminated_dg2),
            "combined_distinct_polynomials": len(eliminated),
        },
        "groebner_after_linear_elimination": {
            "domain": "QQ[alpha,chi,delta,eta]",
            "order": "grevlex",
            "basis": [str(value) for value in groebner_expressions],
            "nonproportional_saturation_charts": saturation,
        },
        "exact_variety": {
            "equations": ["alpha*delta - 1 = 0", "chi*delta - eta = 0"],
            "equivalent_operator_statement": "B_prime=delta*A_prime",
            "mu_lambda": "arbitrary because [A_prime,B_prime]=0",
        },
        "decision": (
            "No noncommuting simultaneous Dolan--Grady pair exists in the enumerated chart; "
            "the complete solution variety consists only of proportional commuting generators."
        ),
    }


def tridiagonal_components(A: Operator, B: Operator, n: int) -> tuple[Operator, ...]:
    first = comm(A, B, n)
    third = comm(A, comm(A, first, n), n)
    middle = mul(mul(A, first, n), A, n)
    square_commutator = add(mul(A, first, n), mul(first, A, n))
    return add(third, smul(2, middle)), middle, square_commutator, first


def tridiagonal_symbolic_system(a: int = 3, b: int = 3) -> dict[str, object]:
    alpha, chi, delta, eta = sp.symbols("alpha chi delta eta")
    beta, gamma, rho, gamma_star, rho_star = sp.symbols("beta gamma rho gamma_star rho_star")
    n, A, B = chart_pair(a, b, alpha, chi, delta, eta)
    left = tridiagonal_components(A, B, n)
    right = tridiagonal_components(B, A, n)
    td1 = simplify_symbolic_operator(
        add(left[0], smul(-beta, left[1]), smul(-gamma, left[2]), smul(-rho, left[3]))
    )
    td2 = simplify_symbolic_operator(
        add(right[0], smul(-beta, right[1]), smul(-gamma_star, right[2]), smul(-rho_star, right[3]))
    )
    group = spatial_group(a, b, True)
    cache: dict[int, int] = {}
    td1_reduced = compress_invariant(td1, group, n, cache)
    td2_reduced = compress_invariant(td2, group, n, cache)
    variables = (alpha, chi, delta, eta, beta, gamma, rho, gamma_star, rho_star)
    return {
        "n": n,
        "variables": variables,
        "A": A,
        "B": B,
        "td1": td1,
        "td2": td2,
        "td1_reduced": td1_reduced,
        "td2_reduced": td2_reduced,
        "equations": normalized_polynomials([*td1_reduced.values(), *td2_reduced.values()], variables),
    }


def tridiagonal_analysis() -> dict[str, object]:
    system = tridiagonal_symbolic_system()
    n = system["n"]
    alpha, chi, delta, eta, beta, gamma, rho, gamma_star, rho_star = system["variables"]
    p = alpha * delta - 1
    q = chi * delta - eta
    r = alpha * eta - chi

    beta_q_word = ((0, "Y"), (1, "Z"), (2, "X"), (3, "X"), (4, "X"))
    beta_q = _coefficient(system["td1_reduced"], n, beta_q_word)
    assert _equal(beta_q, 4 * (beta - 2) * q)
    beta_r_word = ((0, "Y"), (1, "Z"), (2, "Z"), (3, "Z"), (4, "X"), (5, "X"))
    beta_r = _coefficient(system["td1_reduced"], n, beta_r_word)
    assert _equal(beta_r, 4 * (beta - 2) * r)
    beta_p_word = ((0, "Y"), (1, "Z"), (2, "X"), (3, "X"))
    beta_p = _coefficient(system["td1_reduced"], n, beta_p_word)
    beta_p_reduced = sp.expand(beta_p.subs({chi: 0, eta: 0}))
    assert _equal(beta_p_reduced, 4 * (beta - 2) * p)

    q_word = ((0, "Y"), (1, "X"), (2, "Z"), (4, "Z"), (6, "Z"))
    q_coefficient = _coefficient(system["td2_reduced"], n, q_word)
    q_at_beta_two = sp.expand(q_coefficient.subs(beta, 2))
    assert _equal(q_at_beta_two, 24 * q)
    chi_p_word = ((0, "X"), (1, "Y"), (2, "X"), (3, "Z"))
    chi_p = _coefficient(system["td1_reduced"], n, chi_p_word)
    chi_p_reduced = sp.expand(chi_p.subs({beta: 2, eta: chi * delta}))
    assert _equal(chi_p_reduced, 8 * chi * p)

    gamma_word = ((1, "Y"), (3, "Z"))
    gamma_coefficient = _coefficient(system["td2_reduced"], n, gamma_word)
    terminal_substitution = {beta: 2, chi: 0, eta: 0}
    gamma_reduced = sp.expand(gamma_coefficient.subs(terminal_substitution))
    assert _equal(gamma_reduced, -8 * gamma_star * p)
    contradiction_word = ((0, "Y"), (1, "Z"), (2, "Z"), (3, "Z"))
    contradiction = _coefficient(system["td2_reduced"], n, contradiction_word)
    contradiction_reduced = sp.expand(contradiction.subs(terminal_substitution))
    assert _equal(contradiction_reduced, -4 * (gamma_star - 12) * p)

    terminal_constant = sp.expand((-4 * (gamma_star - 12)).subs(gamma_star, 0))
    assert terminal_constant == 48
    unit = sp.groebner([terminal_constant], gamma_star, order="lex")
    unit_basis = [poly.as_expr() for poly in unit.polys]
    assert unit_basis == [sp.Integer(1)]

    td1_distinct = normalized_polynomials(system["td1_reduced"].values(), system["variables"])
    td2_distinct = normalized_polynomials(system["td2_reduced"].values(), system["variables"])
    return {
        "claim_tag": "[COMPUTATION]",
        "relations": {
            "TD1": "[A,A^2 B-beta A B A+B A^2-gamma(AB+BA)-rho B]=0",
            "TD2": "[B,B^2 A-beta B A B+A B^2-gamma_star(BA+AB)-rho_star A]=0",
            "Dolan_Grady_point": {"beta": "2", "gamma": "0", "gamma_star": "0", "rho": "16", "rho_star": "16"},
        },
        "unknowns": {
            "deformation": [str(value) for value in (alpha, chi, delta, eta)],
            "linear_relation_parameters": [str(value) for value in (beta, gamma, rho, gamma_star, rho_star)],
            "total": len(system["variables"]),
        },
        "orbit_reduction": {
            "TD1": {
                "full_pauli_equations": len(system["td1"]),
                "symmetry_orbit_equations": len(system["td1_reduced"]),
                "distinct_up_to_nonzero_rational_scale": len(td1_distinct),
                "full_support_counts": support_counts(system["td1"], n),
                "orbit_support_counts": support_counts(system["td1_reduced"], n),
            },
            "TD2": {
                "full_pauli_equations": len(system["td2"]),
                "symmetry_orbit_equations": len(system["td2_reduced"]),
                "distinct_up_to_nonzero_rational_scale": len(td2_distinct),
                "full_support_counts": support_counts(system["td2"], n),
                "orbit_support_counts": support_counts(system["td2_reduced"], n),
            },
            "row_slots_both_relations": len(system["td1"]) + len(system["td2"]),
            "distinct_full_pauli_words_union": len(set(system["td1"]) | set(system["td2"])),
            "orbit_row_slots_both_relations": len(system["td1_reduced"]) + len(system["td2_reduced"]),
            "distinct_orbit_representatives_union": len(set(system["td1_reduced"]) | set(system["td2_reduced"])),
            "distinct_polynomials_both_relations": len(system["equations"]),
        },
        "low_support_certificate": {
            "stage_1_nonproportional_forces_beta_2": [
                {
                    "condition": "q != 0",
                    "pauli_word": format_pauli(pack_word(beta_q_word, n), n),
                    "support": 5,
                    "coefficient": str(sp.factor(beta_q)),
                },
                {
                    "condition": "q=0 and r != 0",
                    "pauli_word": format_pauli(pack_word(beta_r_word, n), n),
                    "support": 6,
                    "coefficient": str(sp.factor(beta_r)),
                },
                {
                    "condition": "q=r=0 and p != 0, hence chi=eta=0",
                    "pauli_word": format_pauli(pack_word(beta_p_word, n), n),
                    "support": 4,
                    "coefficient_after_case_substitution": str(sp.factor(beta_p_reduced)),
                },
            ],
            "stage_2_at_beta_2": [
                {
                    "pauli_word": format_pauli(pack_word(q_word, n), n),
                    "support": 5,
                    "coefficient": str(sp.factor(q_at_beta_two)),
                    "consequence": "q=0",
                },
                {
                    "pauli_word": format_pauli(pack_word(chi_p_word, n), n),
                    "support": 4,
                    "coefficient_after_q_0": str(sp.factor(chi_p_reduced)),
                    "consequence": "on the nonproportional p!=0 chart, chi=eta=0",
                },
            ],
            "stage_3_terminal_inconsistency": [
                {
                    "pauli_word": format_pauli(pack_word(gamma_word, n), n),
                    "support": 2,
                    "coefficient": str(sp.factor(gamma_reduced)),
                    "consequence_on_p_nonzero": "gamma_star=0",
                },
                {
                    "pauli_word": format_pauli(pack_word(contradiction_word, n), n),
                    "support": 4,
                    "coefficient": str(sp.factor(contradiction_reduced)),
                    "consequence_on_p_nonzero": "gamma_star=12",
                },
            ],
        },
        "linear_elimination_then_groebner": {
            "eliminated_variable": "gamma_star=0 from the support-2 terminal row",
            "remaining_exact_constant": str(terminal_constant),
            "groebner_basis": [str(value) for value in unit_basis],
        },
        "exact_variety": {
            "equations": ["alpha*delta - 1 = 0", "chi*delta - eta = 0"],
            "equivalent_operator_statement": "B_prime=delta*A_prime",
            "five_tridiagonal_parameters": "arbitrary because [A_prime,B_prime]=0",
        },
        "decision": (
            "No noncommuting simultaneous tridiagonal/Askey--Wilson pair exists in the enumerated chart; "
            "the complete solution variety is the same proportional commuting component as for Dolan--Grady."
        ),
    }


def dg_residual(A: Operator, B: Operator, n: int, constant: object) -> Operator:
    first = comm(A, B, n)
    return simplify_symbolic_operator(add(comm(A, comm(A, first, n), n), smul(-constant, first)))


def tridiagonal_residual(
    A: Operator,
    B: Operator,
    n: int,
    beta: object,
    gamma: object,
    rho: object,
) -> Operator:
    T0, T1, T2, T3 = tridiagonal_components(A, B, n)
    return simplify_symbolic_operator(add(T0, smul(-beta, T1), smul(-gamma, T2), smul(-rho, T3)))


def multiple_size_checks() -> dict[str, object]:
    coefficients = {
        "alpha": Fraction(1, 2),
        "chi": Fraction(-1, 2),
        "delta": Fraction(2),
        "eta": Fraction(-1),
    }
    rows = []
    for a, b in ((3, 3), (4, 4)):
        n, A, B = chart_pair(a, b, **coefficients)
        commutator = comm(A, B, n)
        dg1 = dg_residual(A, B, n, Fraction(7))
        dg2 = dg_residual(B, A, n, Fraction(11))
        td1 = tridiagonal_residual(A, B, n, Fraction(3), Fraction(5), Fraction(7))
        td2 = tridiagonal_residual(B, A, n, Fraction(3), Fraction(11), Fraction(13))
        assert not commutator and not dg1 and not dg2 and not td1 and not td2

        _, (O_X, O_E, O_C) = minimal_basis(a, b)
        A_f3 = add(O_X, smul(Fraction(-1, 2), O_C))
        f3_dg2 = dg_residual(O_E, A_f3, n, Fraction(16))
        assert not f3_dg2
        f3_first = comm(A_f3, O_E, n)
        f3_third = comm(A_f3, comm(A_f3, f3_first, n), n)
        ratios = {
            Fraction(f3_third.get(vector, 0), coefficient)
            for vector, coefficient in f3_first.items()
            if coefficient
        }
        f3_dg1_solvable = any(
            all(
                f3_third.get(vector, 0) == ratio * f3_first.get(vector, 0)
                for vector in set(f3_first) | set(f3_third)
            )
            for ratio in ratios
        )
        assert f3_dg1_solvable is False
        rows.append(
            {
                "layer": f"{a}x{b}_torus",
                "sites": n,
                "fixed_coefficients": {key: str(value) for key, value in coefficients.items()},
                "identity": "B_prime=2*A_prime",
                "commutator_terms": len(commutator),
                "DG1_residual_terms_at_mu_7": len(dg1),
                "DG2_residual_terms_at_lambda_11": len(dg2),
                "TD1_residual_terms_at_beta_3_gamma_5_rho_7": len(td1),
                "TD2_residual_terms_at_beta_3_gamma_star_11_rho_star_13": len(td2),
                "F3_DG2_residual_terms_at_lambda_16": len(f3_dg2),
                "F3_DG1_floating_constant_solution_exists": f3_dg1_solvable,
            }
        )
    return {
        "claim_tag": "[COMPUTATION]",
        "sizes": rows,
        "interpretation": (
            "The only exact solution component is checked with one fixed rational coefficient tuple on two sizes, "
            "but it is proportional and commuting, hence not an integrable tridiagonal pair.  The noncommuting "
            "F3 point cancels DG2 on both sizes and still fails DG1."
        ),
    }


def run_all(write_results: bool = True) -> dict[str, object]:
    controls = wave4_controls()

    basis_3 = basis_enumeration(3, 3)
    basis_4 = basis_enumeration(4, 4)
    assert [row["literal_pauli_terms"] for row in basis_3["orbits"]] == [9, 18, 18]
    assert basis_3["group_order"] == 72
    assert [row["literal_pauli_terms"] for row in basis_4["orbits"]] == [16, 32, 32]
    assert basis_4["group_order"] == 128

    dg = dg_analysis()
    td = tridiagonal_analysis()
    sizes = multiple_size_checks()

    ansatz = {
        "claim_tag": "[LEMMA]",
        "common_orbit_basis": ["O_X", "O_E", "O_C"],
        "common_orbit_basis_dimension": 3,
        "raw_pair_coefficients": 6,
        "projective_chart": {
            "A_prime": "O_X + alpha*O_E + chi*O_C",
            "B_prime": "delta*O_X + O_E + eta*O_C",
            "free_deformation_coefficients": 4,
            "fixed_nonzero_coefficients": [
                "coefficient(O_X in A_prime)=1",
                "coefficient(O_E in B_prime)=1",
            ],
        },
        "successful_F3_point": {"alpha": "0", "chi": "-1/2", "delta": "0", "eta": "0"},
        "minimality": (
            "O_X, O_E, and O_C have disjoint Pauli supports and are linearly independent. "
            "The F3 pair uses all three, so no smaller common Pauli-orbit basis contains it."
        ),
        "strict_extension": (
            "Both generators independently access all three orbit sums; the F3 point used only "
            "O_X,O_C in A and O_E in B."
        ),
        "enumeration_by_size": [basis_3, basis_4],
        "covered": (
            "Every complex/rational coefficient tuple in the stated four-dimensional chart, "
            "with arbitrary direct DG constants or arbitrary five tridiagonal parameters."
        ),
        "not_covered": (
            "The projective hyperplanes coefficient(O_X in A_prime)=0 or coefficient(O_E in B_prime)=0; "
            "non-invariant coefficients; orbit supports beyond O_X,O_E,O_C; and pairs outside this basis."
        ),
    }

    checks = [
        {
            "name": "wave-4 controls reproduced before search",
            "passed": True,
            "detail": "C6 zero, P6-plus-bond obstruction, F3 DG2 zero, DG1 witness (3,0)",
        },
        {
            "name": "minimal orbit basis enumeration",
            "passed": True,
            "detail": "3 independent single-orbit sums; 6 raw pair coefficients, 4 in chart",
        },
        {
            "name": "DG low-support certificate",
            "passed": "No noncommuting" in dg["decision"],
            "detail": dg["decision"],
        },
        {
            "name": "DG Groebner after linear elimination",
            "passed": all(
                value == ["1"]
                for value in dg["groebner_after_linear_elimination"]["nonproportional_saturation_charts"].values()
            ),
            "detail": "all three nonproportional cross-product charts have unit ideal",
        },
        {
            "name": "tridiagonal low-support and terminal unit certificate",
            "passed": td["linear_elimination_then_groebner"]["groebner_basis"] == ["1"],
            "detail": td["decision"],
        },
        {
            "name": "fixed coefficients checked on two layer sizes",
            "passed": len(sizes["sizes"]) == 2 and all(row["commutator_terms"] == 0 for row in sizes["sizes"]),
            "detail": "same rational tuple checked exactly on 3x3 and 4x4 torus layers",
        },
    ]
    assert all(check["passed"] for check in checks)

    report = {
        "provenance": {
            "script": "experiments/e36_dg_simultaneous.py",
            "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "precision": "exact; no floating-point arithmetic",
            "arithmetic": "Python integers, fractions.Fraction, and SymPy QQ polynomials",
            "sympy_version": sp.__version__,
            "pauli_basis": "phase-free Q_(a|b)=X^a Z^b",
            "frozen_control_source": "experiments/e17_dg_deformation.py",
        },
        "data": {
            "scope": ansatz,
            "controls": controls,
            "dolan_grady": dg,
            "tridiagonal": td,
            "multiple_size_checks": sizes,
            "decision": {
                "claim_tag": "[THEOREM]",
                "statement": (
                    "Within the precisely enumerated four-dimensional projective chart on the minimal "
                    "three-orbit common basis, both the simultaneous Dolan--Grady and simultaneous "
                    "tridiagonal noncommuting varieties are empty.  Their full solution varieties "
                    "consist only of B_prime=delta*A_prime."
                ),
                "integrability": False,
                "reason": "proportional commuting generators are degenerate, not a tridiagonal pair",
                "size_scope": (
                    "The empty noncommuting variety is certified already on 3x3.  Therefore no "
                    "size-independent coefficient tuple in this chart can work on both 3x3 and 4x4. "
                    "The only degenerate component was nevertheless checked on both sizes."
                ),
            },
        },
        "checks": checks,
    }
    if write_results:
        RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
        RESULT_PATH.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    try:
        report = run_all(write_results=True)
        for check in report["checks"]:
            print(f"{check['name']}: {'PASS' if check['passed'] else 'FAIL'}")
        print("PASS" if all(check["passed"] for check in report["checks"]) else "FAIL")
    except Exception as error:
        print(f"failure: {type(error).__name__}: {error}")
        print("FAIL")
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
