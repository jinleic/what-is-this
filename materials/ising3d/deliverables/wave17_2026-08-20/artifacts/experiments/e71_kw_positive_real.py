"""Exact Positivstellensatz no-go for nonnegative scalar Kac--Ward weights.

The full 30-weight family remains unresolved over Q/C.  This experiment proves a
strictly scoped real theorem: already the order-four equation on the free 2x2x2
box has no point in the nonnegative orthant.  It also specializes and verifies
that certificate on every branch of the existing five-edge gauge-tree map.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import sympy as sp

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    ALLOWED_DIRECTION_PAIRS,
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    OPPOSITE_DIRECTION,
    formal_log_coefficients,
    full_weight_finite_system,
)
from ising.lattices import cubic

ROOT = Path(__file__).resolve().parents[1]
FULL_FAMILY_SOURCE = ROOT / "results/kac_ward/full_family.json"
BRANCH_SOURCE = ROOT / "results/kac_ward/branches.json"
OUTPUT = ROOT / "results/kac_ward/positive_real.json"
CERTIFICATE_SHAPE = (2, 2, 2)
CERTIFICATE_ORDER = 4
EXPECTED_EVEN_SUBGRAPH_POLYNOMIAL = [1, 0, 0, 0, 6, 0, 16, 0, 9, 0, 0, 0, 0]


def pair_name(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def rational_record(value: int | Fraction | sp.Rational) -> dict[str, int]:
    rational = sp.Rational(value)
    return {"numerator": int(rational.p), "denominator": int(rational.q)}


def rational_value(record: dict[str, int]) -> sp.Rational:
    return sp.Rational(record["numerator"], record["denominator"])


def rational_squares(value: sp.Rational) -> list[sp.Rational]:
    """Give an exact rational four-squares decomposition of a positive rational."""

    value = sp.Rational(value)
    if value <= 0:
        raise ValueError("the SOS coefficient must be positive")
    numerator = int(value.p)
    denominator = int(value.q)
    target = numerator * denominator
    limit = math.isqrt(target)
    for a in range(limit, -1, -1):
        remainder_a = target - a * a
        for b in range(math.isqrt(remainder_a), -1, -1):
            remainder_b = remainder_a - b * b
            for c in range(math.isqrt(remainder_b), -1, -1):
                remainder_c = remainder_b - c * c
                d = math.isqrt(remainder_c)
                if d * d == remainder_c:
                    squares = [
                        sp.Rational(entry, denominator)
                        for entry in (a, b, c, d)
                        if entry
                    ]
                    if sum(entry * entry for entry in squares) != value:
                        raise AssertionError("internal four-squares reconstruction failed")
                    return squares
    raise AssertionError(f"no rational four-squares decomposition found for {value}")


def term_record(coefficient: sp.Rational, exponents: tuple[int, ...], variables) -> dict:
    return {
        "coefficient": rational_record(coefficient),
        "powers": [
            {"variable": str(variable), "exponent": int(exponent)}
            for variable, exponent in zip(variables, exponents, strict=True)
            if exponent
        ],
    }


def sparse_polynomial_record(expression: sp.Expr, variables) -> dict:
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    terms = [
        term_record(coefficient, exponents, variables)
        for exponents, coefficient in polynomial.terms()
        if any(exponents)
    ]
    constant = polynomial.coeff_monomial((0,) * len(variables))
    return {
        "variables": [str(variable) for variable in variables],
        "constant": rational_record(constant),
        "nonconstant_terms": terms,
        "expanded_expression": str(sp.expand(expression)),
        "sha256_sympy_srepr": hashlib.sha256(
            sp.srepr(sp.expand(expression)).encode("utf-8")
        ).hexdigest(),
    }


def expression_from_term(record: dict, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    expression = rational_value(record["coefficient"])
    for power in record["powers"]:
        expression *= symbols[power["variable"]] ** int(power["exponent"])
    return expression


def preordering_certificate(expression: sp.Expr, variables) -> dict:
    """Return -1 = ideal multiplier * F + an exact rational preordering term."""

    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    zero_exponents = (0,) * len(variables)
    constant = polynomial.coeff_monomial(zero_exponents)
    if constant <= 0:
        raise ValueError("certificate requires a positive constant term")

    positive_terms = []
    for exponents, coefficient in polynomial.terms():
        if not any(exponents):
            continue
        if coefficient <= 0:
            raise ValueError("certificate requires positive nonconstant coefficients")
        normalized = sp.Rational(coefficient, constant)
        scalar_squares = rational_squares(normalized)
        positive_terms.append(
            {
                "coefficient": rational_record(normalized),
                "preordering_generators": [
                    str(variable)
                    for variable, exponent in zip(variables, exponents, strict=True)
                    if exponent % 2
                ],
                "square_monomial": [
                    {"variable": str(variable), "exponent": int(exponent // 2)}
                    for variable, exponent in zip(variables, exponents, strict=True)
                    if exponent >= 2
                ],
                "sos_scalar_squares": [
                    rational_record(entry) for entry in scalar_squares
                ],
                "source_term": term_record(coefficient, exponents, variables),
            }
        )

    certificate = {
        "type": "rational_preordering_Positivstellensatz_identity",
        "equation_ideal_multiplier": rational_record(-sp.Rational(1, constant)),
        "target": rational_record(-1),
        "domain_generators": [f"{variable} >= 0" for variable in variables],
        "positive_terms": positive_terms,
        "certified_lower_bound_for_equation": rational_record(constant),
        "identity_template": (
            "-1 = (-1/constant)*F + sum_alpha "
            "[(sum_j scalar_square_j^2)*(square_monomial)^2*product(odd generators)]"
        ),
    }
    verify_preordering_certificate(expression, variables, certificate)
    return certificate


def verify_preordering_certificate(
    expression: sp.Expr, variables, certificate: dict
) -> None:
    symbols = {str(variable): variable for variable in variables}
    positive_part = sp.Integer(0)
    for term in certificate["positive_terms"]:
        scalar_squares = [
            rational_value(record) for record in term["sos_scalar_squares"]
        ]
        coefficient = sum(value * value for value in scalar_squares)
        if coefficient != rational_value(term["coefficient"]):
            raise AssertionError("stored rational SOS coefficient is incorrect")
        square_monomial = sp.Integer(1)
        for power in term["square_monomial"]:
            square_monomial *= symbols[power["variable"]] ** int(power["exponent"])
        generator_product = sp.Integer(1)
        for name in term["preordering_generators"]:
            generator_product *= symbols[name]
        represented = sp.expand(coefficient * square_monomial**2 * generator_product)
        source_term = expression_from_term(term["source_term"], symbols)
        normalized_source = sp.expand(
            source_term
            / rational_value(certificate["certified_lower_bound_for_equation"])
        )
        if represented != normalized_source:
            raise AssertionError("preordering term does not reconstruct its source monomial")
        positive_part += represented

    identity = sp.expand(
        rational_value(certificate["equation_ideal_multiplier"]) * expression
        + positive_part
    )
    if identity != rational_value(certificate["target"]):
        raise AssertionError(f"Positivstellensatz identity failed: {identity}")


def branch_expression(system, branch: str) -> tuple[tuple[sp.Symbol, ...], sp.Expr]:
    if len(branch) != len(DIRECTION_GAUGE_TREE) or set(branch) - {"0", "1"}:
        raise ValueError(f"invalid branch {branch!r}")
    symbol_map = dict(system.all_symbols)
    substitutions = {
        symbol_map[pair]: int(bit)
        for pair, bit in zip(DIRECTION_GAUGE_TREE, branch, strict=True)
    }
    remaining = tuple(
        symbol for pair, symbol in system.all_symbols if pair not in DIRECTION_GAUGE_TREE
    )
    return remaining, sp.expand(system.equations[0].subs(substitutions))


def monomial_pair_records(system, expression: sp.Expr) -> list[dict]:
    pair_by_symbol = {symbol: pair for pair, symbol in system.all_symbols}
    polynomial = sp.Poly(expression, *system.variables, domain=sp.QQ)
    monomials: list[tuple[frozenset[tuple[int, int]], tuple[str, ...]]] = []
    for exponents, coefficient in polynomial.terms():
        if not any(exponents):
            continue
        if coefficient != 8 or any(exponent not in (0, 1) for exponent in exponents):
            raise AssertionError("unexpected order-four plaquette term")
        pairs = frozenset(
            pair_by_symbol[variable]
            for variable, exponent in zip(system.variables, exponents, strict=True)
            if exponent
        )
        names = tuple(
            str(variable)
            for variable, exponent in zip(system.variables, exponents, strict=True)
            if exponent
        )
        monomials.append((pairs, names))

    by_pairs = {pairs: names for pairs, names in monomials}
    records = []
    consumed: set[frozenset[tuple[int, int]]] = set()
    axis_names = ("x", "x", "y", "y", "z", "z")
    for pairs, names in monomials:
        if pairs in consumed:
            continue
        reversed_pairs = frozenset(
            (OPPOSITE_DIRECTION[next_direction], OPPOSITE_DIRECTION[direction])
            for direction, next_direction in pairs
        )
        if reversed_pairs not in by_pairs or reversed_pairs == pairs:
            raise AssertionError("plaquette monomials did not pair under path reversal")
        axes = sorted({axis_names[index] for pair in pairs for index in pair})
        if len(axes) != 2:
            raise AssertionError("plaquette monomial does not use exactly two axes")
        consumed.update((pairs, reversed_pairs))
        records.append(
            {
                "plane": "".join(axes),
                "representative_product": list(names),
                "reversed_product": list(by_pairs[reversed_pairs]),
                "reversal_map": [
                    {
                        "from": pair_name(pair),
                        "to": pair_name(
                            (
                                OPPOSITE_DIRECTION[pair[1]],
                                OPPOSITE_DIRECTION[pair[0]],
                            )
                        ),
                    }
                    for pair in sorted(pairs)
                ],
            }
        )
    return sorted(records, key=lambda record: record["plane"])


def main() -> None:
    full_source = json.loads(FULL_FAMILY_SOURCE.read_text(encoding="utf-8"))
    branch_source = json.loads(BRANCH_SOURCE.read_text(encoding="utf-8"))
    system = full_weight_finite_system(
        (CERTIFICATE_SHAPE,), (CERTIFICATE_ORDER,), gauge_fix=False
    )
    expression = sp.expand(system.equations[0])
    polynomial_record = sparse_polynomial_record(expression, system.variables)
    global_certificate = preordering_certificate(expression, system.variables)

    even_polynomial = even_subgraph_polynomial(
        cubic(*CERTIFICATE_SHAPE, periodic=False)
    )
    log_coefficients = formal_log_coefficients(even_polynomial, CERTIFICATE_ORDER)
    log_coefficient = log_coefficients[CERTIFICATE_ORDER]

    expected_labels = [
        "".join(str(bit) for bit in bits)
        for bits in itertools.product((0, 1), repeat=len(DIRECTION_GAUGE_TREE))
    ]
    source_records = branch_source["data"]["branches"]
    source_by_label = {record["branch"]: record for record in source_records}
    gauge_names = [pair_name(pair) for pair in DIRECTION_GAUGE_TREE]

    branch_records = []
    for branch in expected_labels:
        remaining, specialized = branch_expression(system, branch)
        certificate = preordering_certificate(specialized, remaining)
        zeros = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "0"
        ]
        ones = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "1"
        ]
        branch_records.append(
            {
                "branch": branch,
                "bits_in_gauge_tree_order": [int(bit) for bit in branch],
                "zero_tree_weights": zeros,
                "normalized_positive_tree_weights": ones,
                "source_complex_Q_status": source_by_label[branch]["status"],
                "positive_real_status": "EMPTY_OVER_R_NONNEGATIVE",
                "claim_tag": "[THEOREM]",
                "physical_scope": (
                    "strict_positive_physical_chart"
                    if branch == "11111"
                    else "nonnegative_boundary_extension"
                ),
                "specialized_equation": sparse_polynomial_record(
                    specialized, remaining
                ),
                "certificate": certificate,
            }
        )

    reversal_pairs = monomial_pair_records(system, expression)
    full_minimal = full_source["data"]["minimal_unresolved_system"]
    source_gauge_order = branch_source["data"]["gauge_tree_order"]
    strict_physical_records = [
        record
        for record in branch_records
        if record["physical_scope"] == "strict_positive_physical_chart"
    ]

    checks = [
        {
            "name": "full_30_weight_parameterization_reconstructed",
            "passed": (
                len(ALLOWED_DIRECTION_PAIRS) == 30
                and len(system.variables) == 30
                and full_minimal["variable_count"] == 25
                and len(DIRECTION_GAUGE_TREE) == 5
            ),
            "detail": {
                "ungauged_variables": len(system.variables),
                "stored_generic_chart_variables": full_minimal["variable_count"],
                "gauge_tree_edges": len(DIRECTION_GAUGE_TREE),
            },
        },
        {
            "name": "source_branch_map_reproduced_exactly",
            "passed": (
                source_gauge_order == gauge_names
                and [record["branch"] for record in source_records] == expected_labels
                and len(source_by_label) == 32
            ),
            "detail": {
                "branch_count": len(source_by_label),
                "gauge_tree_order": gauge_names,
            },
        },
        {
            "name": "cube_even_subgraph_and_log_coefficient_exact",
            "passed": (
                even_polynomial == EXPECTED_EVEN_SUBGRAPH_POLYNOMIAL
                and log_coefficient == 6
            ),
            "detail": {
                "even_subgraph_polynomial": even_polynomial,
                "log_P_v4": rational_record(log_coefficient),
            },
        },
        {
            "name": "cube_order_four_equation_has_positive_margin_48",
            "passed": (
                rational_value(polynomial_record["constant"]) == 48
                and len(polynomial_record["nonconstant_terms"]) == 6
                and all(
                    rational_value(term["coefficient"]) == 8
                    and len(term["powers"]) == 4
                    and all(power["exponent"] == 1 for power in term["powers"])
                    for term in polynomial_record["nonconstant_terms"]
                )
            ),
            "detail": {
                "constant": polynomial_record["constant"],
                "nonconstant_term_count": len(
                    polynomial_record["nonconstant_terms"]
                ),
            },
        },
        {
            "name": "global_rational_Positivstellensatz_identity_rechecked",
            "passed": True,
            "detail": {
                "target": global_certificate["target"],
                "positive_term_count": len(global_certificate["positive_terms"]),
            },
        },
        {
            "name": "all_32_branch_Positivstellensatz_identities_rechecked",
            "passed": (
                len(branch_records) == 32
                and all(
                    record["positive_real_status"] == "EMPTY_OVER_R_NONNEGATIVE"
                    and rational_value(
                        record["certificate"]["certified_lower_bound_for_equation"]
                    )
                    == 48
                    for record in branch_records
                )
            ),
            "detail": {"excluded_nonnegative_branches": len(branch_records)},
        },
        {
            "name": "strict_positive_domain_is_generic_11111_chart",
            "passed": (
                len(strict_physical_records) == 1
                and strict_physical_records[0]["branch"] == "11111"
                and not strict_physical_records[0]["zero_tree_weights"]
            ),
            "detail": {
                "strict_positive_branches": [
                    record["branch"] for record in strict_physical_records
                ]
            },
        },
        {
            "name": "plaquette_products_pair_under_conditional_reversal_conjugation",
            "passed": (
                [record["plane"] for record in reversal_pairs] == ["xy", "xz", "yz"]
                and all(
                    len(record["representative_product"])
                    == len(record["reversed_product"])
                    == 4
                    for record in reversal_pairs
                )
            ),
            "detail": {"paired_planes": [record["plane"] for record in reversal_pairs]},
        },
        {
            "name": "no_floating_point_or_benchmark_decision",
            "passed": True,
            "detail": "All claims use Python integers, fractions, and SymPy QQ; no numerical optimization and no K_c input.",
        },
    ]

    output = {
        "provenance": {
            "script": "experiments/e71_kw_positive_real.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "arithmetic": "exact Python integers/Fraction and SymPy QQ only",
            "numerical_precision": None,
            "numerical_candidate_search": "not run; exact order-four sign certificate is decisive",
            "source_artifacts": [
                "results/kac_ward/full_family.json",
                "results/kac_ward/branches.json",
            ],
            "source_system_sha256": full_minimal["sha256_sympy_srepr"],
        },
        "data": {
            "headline": (
                "[THEOREM] The translation-invariant scalar 30-weight Kac--Ward "
                "equations have no nonnegative-real point; the free 2x2x2 order-four "
                "equation alone has exact lower bound 48."
            ),
            "parameterization": {
                "directions": list(DIRECTION_LABELS),
                "allowed_ordered_pairs": [
                    pair_name(pair) for pair in ALLOWED_DIRECTION_PAIRS
                ],
                "forbidden_transition_rule": "d'=-d (immediate reversal) is absent",
                "weight_count": 30,
                "orbit_description": "6 straight direction pairs plus 24 orthogonal direction pairs",
                "gauge_action": "U(d,d') -> g(d)^(-1) U(d,d') g(d')",
                "gauge_tree_order": gauge_names,
                "generic_chart_variable_count": 25,
            },
            "physical_domain": {
                "claim_tag": "[DEFINITION]",
                "edge_weight": {
                    "symbol": "v=tanh(K)",
                    "finite_ferromagnetic_domain": "K>0, hence 0<v<1",
                    "role": "one common scalar factor in Lambda(v)=v*U; no separate anisotropic edge variables occur in this ansatz",
                },
                "strict_positive_corner_transition_domain": (
                    "U(d,d') is a finite real number >0 for every one of the 30 allowed pairs"
                ),
                "nonnegative_extension_certified": (
                    "U(d,d')>=0 for every allowed pair; zeros are included only as a closed-boundary strengthening"
                ),
                "reality_constraint": "U(d,d')=conjugate(U(d,d')) because every U is real",
                "reversal_conjugation_constraint": (
                    "not imposed by the full-family parameterization and not used in the no-go; "
                    "if a convention additionally imposes U(-d',-d)=conjugate(U(d,d')), "
                    "that is a subset of the excluded positive-real domain"
                ),
                "nonzero_and_branch_scope": (
                    "strict positivity makes all five gauge-tree entries nonzero; a positive "
                    "direction-state gauge normalizes them to branch 11111"
                ),
            },
            "small_box_certificate": {
                "shape": list(CERTIFICATE_SHAPE),
                "order": CERTIFICATE_ORDER,
                "normalization": "F_(B,4)=Tr_B(U^4)+8*[v^4]log(P_B(v))",
                "even_subgraph_polynomial": even_polynomial,
                "log_P_v4": rational_record(log_coefficient),
                "equation": polynomial_record,
                "global_nonnegative_orthant_certificate": global_certificate,
                "exact_conclusion": (
                    "[THEOREM] F_(2x2x2,4)>=48 on U>=0, so F=0 has no "
                    "nonnegative-real solution; on U>0 the six plaquette terms are positive."
                ),
            },
            "sign_parity_and_conjugation_consequences": {
                "oriented_plaquette_pairing": reversal_pairs,
                "real_signed_necessary_condition": (
                    "[LEMMA] Any real signed solution of the cube k=4 equation must have "
                    "the sum of the six oriented plaquette products equal to -6; hence at "
                    "least one product is negative and has an odd number of negative factors."
                ),
                "conditional_reversal_conjugation": (
                    "[LEMMA] If U(-d',-d)=conjugate(U(d,d')), the three paired "
                    "representative holonomies H_xy,H_xz,H_yz obey "
                    "Re(H_xy)+Re(H_xz)+Re(H_yz)=-3."
                ),
                "conditional_unit_modulus_corollary": (
                    "[LEMMA] If the same reversal-conjugation condition holds and every "
                    "corner weight has modulus one, each representative plaquette holonomy "
                    "must equal -1. Unit modulus is not assumed in the no-go theorem."
                ),
            },
            "branch_map": {
                "bit_convention": branch_source["data"]["bit_convention"],
                "branch_count": len(branch_records),
                "strict_positive_physical_branch": "11111",
                "positive_real_status_counts": {
                    "EMPTY_OVER_R_NONNEGATIVE": len(branch_records)
                },
                "source_complex_Q_status_counts": branch_source["data"]["status_counts"],
                "branches": branch_records,
            },
            "candidate_search": {
                "positive_real_candidate_count": 0,
                "high_precision_holdout_status": "not_applicable",
                "reason": "exact rational Positivstellensatz certificate excludes the domain before numerical scouting",
            },
            "scope": (
                "[SCOPE] This proves a positive/nonnegative-real scalar-weight no-go for the "
                "coefficientwise all-box determinant identity. It does not exclude signed-real "
                "or complex corner weights (including planar half-angle phases), does not decide "
                "the 15 source branches unresolved over Q/C, and does not address position-dependent, "
                "matrix-valued, spin-structure-summed, or fixed-v-only constructions."
            ),
        },
        "checks": checks,
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if not all(check["passed"] for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        raise AssertionError(f"failed checks: {failed}")
    print("PASS")


if __name__ == "__main__":
    main()
