"""Independent exact checks for the positive-real Kac--Ward certificates."""

from __future__ import annotations

import hashlib
import itertools
import json
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
RESULT = ROOT / "results/kac_ward/positive_real.json"
SHAPE = (2, 2, 2)
ORDER = 4
EXPECTED_EVEN_POLYNOMIAL = [1, 0, 0, 0, 6, 0, 16, 0, 9, 0, 0, 0, 0]


def pair_name(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def q(record: dict[str, int]) -> sp.Rational:
    assert set(record) == {"numerator", "denominator"}
    assert isinstance(record["numerator"], int)
    assert isinstance(record["denominator"], int)
    assert record["denominator"] > 0
    return sp.Rational(record["numerator"], record["denominator"])


def term_expression(record: dict, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    expression = q(record["coefficient"])
    for power in record["powers"]:
        assert power["variable"] in symbols
        assert isinstance(power["exponent"], int) and power["exponent"] > 0
        expression *= symbols[power["variable"]] ** power["exponent"]
    return sp.expand(expression)


def check_polynomial_record(
    expression: sp.Expr, variables: tuple[sp.Symbol, ...], record: dict
) -> None:
    assert record["variables"] == [str(variable) for variable in variables]
    symbols = {str(variable): variable for variable in variables}
    reconstructed = q(record["constant"]) + sum(
        (term_expression(term, symbols) for term in record["nonconstant_terms"]),
        sp.Integer(0),
    )
    assert sp.expand(reconstructed - expression) == 0
    assert record["expanded_expression"] == str(sp.expand(expression))
    assert record["sha256_sympy_srepr"] == hashlib.sha256(
        sp.srepr(sp.expand(expression)).encode("utf-8")
    ).hexdigest()


def check_certificate(
    expression: sp.Expr, variables: tuple[sp.Symbol, ...], certificate: dict
) -> None:
    assert certificate["type"] == "rational_preordering_Positivstellensatz_identity"
    assert q(certificate["target"]) == -1
    constant = sp.Poly(expression, *variables, domain=sp.QQ).coeff_monomial(
        (0,) * len(variables)
    )
    assert constant > 0
    assert q(certificate["certified_lower_bound_for_equation"]) == constant
    assert q(certificate["equation_ideal_multiplier"]) == -sp.Rational(1, constant)
    assert certificate["domain_generators"] == [
        f"{variable} >= 0" for variable in variables
    ]

    symbols = {str(variable): variable for variable in variables}
    positive_part = sp.Integer(0)
    for term in certificate["positive_terms"]:
        scalar_squares = [q(record) for record in term["sos_scalar_squares"]]
        assert scalar_squares
        coefficient = sum(value * value for value in scalar_squares)
        assert coefficient == q(term["coefficient"]) > 0

        square_monomial = sp.Integer(1)
        for power in term["square_monomial"]:
            assert power["variable"] in symbols
            assert power["exponent"] > 0
            square_monomial *= symbols[power["variable"]] ** power["exponent"]

        generator_product = sp.Integer(1)
        assert len(term["preordering_generators"]) == len(
            set(term["preordering_generators"])
        )
        for name in term["preordering_generators"]:
            assert name in symbols
            generator_product *= symbols[name]

        represented = sp.expand(coefficient * square_monomial**2 * generator_product)
        source = term_expression(term["source_term"], symbols)
        assert represented == sp.expand(source / constant)
        positive_part += represented

    identity = sp.expand(
        q(certificate["equation_ideal_multiplier"]) * expression + positive_part
    )
    assert identity == -1


def load_payloads() -> tuple[dict, dict, dict]:
    return (
        json.loads(FULL_FAMILY_SOURCE.read_text(encoding="utf-8")),
        json.loads(BRANCH_SOURCE.read_text(encoding="utf-8")),
        json.loads(RESULT.read_text(encoding="utf-8")),
    )


def fresh_system():
    return full_weight_finite_system((SHAPE,), (ORDER,), gauge_fix=False)


def check_envelope_parameterization_and_scope() -> None:
    full_source, branch_source, payload = load_payloads()
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e71_kw_positive_real.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert payload["provenance"]["arithmetic"] == (
        "exact Python integers/Fraction and SymPy QQ only"
    )
    assert payload["provenance"]["numerical_precision"] is None
    assert all(check["passed"] for check in payload["checks"])
    assert payload["data"]["headline"].startswith("[THEOREM]")
    assert payload["data"]["scope"].startswith("[SCOPE]")
    assert "does not exclude signed-real or complex" in payload["data"]["scope"]

    parameterization = payload["data"]["parameterization"]
    assert parameterization["directions"] == list(DIRECTION_LABELS)
    assert parameterization["allowed_ordered_pairs"] == [
        pair_name(pair) for pair in ALLOWED_DIRECTION_PAIRS
    ]
    assert parameterization["weight_count"] == len(ALLOWED_DIRECTION_PAIRS) == 30
    assert parameterization["generic_chart_variable_count"] == 25
    assert full_source["data"]["minimal_unresolved_system"]["variable_count"] == 25
    assert payload["provenance"]["source_system_sha256"] == full_source["data"][
        "minimal_unresolved_system"
    ]["sha256_sympy_srepr"]

    physical = payload["data"]["physical_domain"]
    assert physical["claim_tag"] == "[DEFINITION]"
    assert physical["edge_weight"]["finite_ferromagnetic_domain"] == "K>0, hence 0<v<1"
    assert ">0 for every one of the 30 allowed pairs" in physical[
        "strict_positive_corner_transition_domain"
    ]
    assert "not imposed" in physical["reversal_conjugation_constraint"]
    assert "branch 11111" in physical["nonzero_and_branch_scope"]

    gauge_names = [pair_name(pair) for pair in DIRECTION_GAUGE_TREE]
    assert parameterization["gauge_tree_order"] == gauge_names
    assert branch_source["data"]["gauge_tree_order"] == gauge_names
    undirected_edges = {frozenset(pair) for pair in DIRECTION_GAUGE_TREE}
    assert len(undirected_edges) == 5
    reached = {0}
    while True:
        enlarged = reached | {
            vertex
            for edge in undirected_edges
            if edge & reached
            for vertex in edge
        }
        if enlarged == reached:
            break
        reached = enlarged
    assert reached == set(range(6))


def check_cube_polynomial_and_global_certificate() -> None:
    _full_source, _branch_source, payload = load_payloads()
    system = fresh_system()
    assert len(system.variables) == 30
    assert system.labels == ((SHAPE, ORDER),)
    expression = sp.expand(system.equations[0])

    even_polynomial = even_subgraph_polynomial(cubic(*SHAPE, periodic=False))
    assert even_polynomial == EXPECTED_EVEN_POLYNOMIAL
    log_coefficients = formal_log_coefficients(even_polynomial, ORDER)
    assert log_coefficients[ORDER] == 6

    symbol_by_name = {str(symbol): symbol for symbol in system.variables}
    expected_products = (
        ("u_px_py", "u_mx_my", "u_py_mx", "u_my_px"),
        ("u_px_my", "u_mx_py", "u_py_px", "u_my_mx"),
        ("u_px_pz", "u_mx_mz", "u_pz_mx", "u_mz_px"),
        ("u_px_mz", "u_mx_pz", "u_pz_px", "u_mz_mx"),
        ("u_py_pz", "u_my_mz", "u_pz_my", "u_mz_py"),
        ("u_py_mz", "u_my_pz", "u_pz_py", "u_mz_my"),
    )
    expected_expression = sp.Integer(48)
    for product in expected_products:
        expected_expression += 8 * sp.prod(symbol_by_name[name] for name in product)
    assert sp.expand(expression - expected_expression) == 0

    small_box = payload["data"]["small_box_certificate"]
    assert small_box["shape"] == list(SHAPE)
    assert small_box["order"] == ORDER
    assert small_box["even_subgraph_polynomial"] == even_polynomial
    assert q(small_box["log_P_v4"]) == 6
    check_polynomial_record(expression, system.variables, small_box["equation"])
    check_certificate(
        expression,
        system.variables,
        small_box["global_nonnegative_orthant_certificate"],
    )
    certificate = small_box["global_nonnegative_orthant_certificate"]
    assert len(certificate["positive_terms"]) == 6
    for term in certificate["positive_terms"]:
        assert q(term["coefficient"]) == sp.Rational(1, 6)
        assert [q(value) for value in term["sos_scalar_squares"]] == [
            sp.Rational(1, 3),
            sp.Rational(1, 6),
            sp.Rational(1, 6),
        ]
        assert not term["square_monomial"]
        assert len(term["preordering_generators"]) == 4


def check_all_branch_certificates() -> None:
    _full_source, branch_source, payload = load_payloads()
    system = fresh_system()
    symbol_map = dict(system.all_symbols)
    remaining = tuple(
        symbol for pair, symbol in system.all_symbols if pair not in DIRECTION_GAUGE_TREE
    )
    expected_labels = [format(index, "05b") for index in range(32)]
    source_records = branch_source["data"]["branches"]
    records = payload["data"]["branch_map"]["branches"]
    assert [record["branch"] for record in source_records] == expected_labels
    assert [record["branch"] for record in records] == expected_labels
    assert payload["data"]["branch_map"]["branch_count"] == 32
    assert payload["data"]["branch_map"]["positive_real_status_counts"] == {
        "EMPTY_OVER_R_NONNEGATIVE": 32
    }
    assert payload["data"]["branch_map"]["source_complex_Q_status_counts"] == {
        "EMPTY_OVER_Q": 17,
        "RATIONAL_POINT_FOUND": 0,
        "UNRESOLVED": 15,
    }

    source_by_label = {record["branch"]: record for record in source_records}
    gauge_names = [pair_name(pair) for pair in DIRECTION_GAUGE_TREE]
    strict_records = []
    for record in records:
        branch = record["branch"]
        substitutions = {
            symbol_map[pair]: int(bit)
            for pair, bit in zip(DIRECTION_GAUGE_TREE, branch, strict=True)
        }
        specialized = sp.expand(system.equations[0].subs(substitutions))
        zeros = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "0"
        ]
        ones = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "1"
        ]
        assert record["bits_in_gauge_tree_order"] == [int(bit) for bit in branch]
        assert record["zero_tree_weights"] == zeros
        assert record["normalized_positive_tree_weights"] == ones
        assert record["source_complex_Q_status"] == source_by_label[branch]["status"]
        assert record["positive_real_status"] == "EMPTY_OVER_R_NONNEGATIVE"
        check_polynomial_record(specialized, remaining, record["specialized_equation"])
        check_certificate(specialized, remaining, record["certificate"])
        assert q(record["certificate"]["certified_lower_bound_for_equation"]) == 48
        if record["physical_scope"] == "strict_positive_physical_chart":
            strict_records.append(record)
        else:
            assert record["physical_scope"] == "nonnegative_boundary_extension"

    assert len(strict_records) == 1
    assert strict_records[0]["branch"] == "11111"
    assert strict_records[0]["zero_tree_weights"] == []


def check_sign_parity_and_conditional_conjugation() -> None:
    _full_source, _branch_source, payload = load_payloads()
    system = fresh_system()
    expression = sp.Poly(system.equations[0], *system.variables, domain=sp.QQ)
    pair_by_symbol = {symbol: pair for pair, symbol in system.all_symbols}
    names_by_pairs: dict[frozenset[tuple[int, int]], frozenset[str]] = {}
    monomial_sum = sp.Integer(0)
    for exponents, coefficient in expression.terms():
        if not any(exponents):
            assert coefficient == 48
            continue
        assert coefficient == 8
        assert sum(exponents) == 4
        assert all(exponent in (0, 1) for exponent in exponents)
        symbols = [
            variable
            for variable, exponent in zip(system.variables, exponents, strict=True)
            if exponent
        ]
        pairs = frozenset(pair_by_symbol[symbol] for symbol in symbols)
        names_by_pairs[pairs] = frozenset(str(symbol) for symbol in symbols)
        monomial_sum += sp.prod(symbols)
    assert len(names_by_pairs) == 6
    assert sp.expand(system.equations[0] - (8 * monomial_sum + 48)) == 0

    stored_pairs = payload["data"]["sign_parity_and_conjugation_consequences"][
        "oriented_plaquette_pairing"
    ]
    assert [record["plane"] for record in stored_pairs] == ["xy", "xz", "yz"]
    covered: set[frozenset[tuple[int, int]]] = set()
    symbol_pair_by_name = {str(symbol): pair for pair, symbol in system.all_symbols}
    for record in stored_pairs:
        representative = frozenset(
            symbol_pair_by_name[name] for name in record["representative_product"]
        )
        reversed_pairs = frozenset(
            (OPPOSITE_DIRECTION[next_direction], OPPOSITE_DIRECTION[direction])
            for direction, next_direction in representative
        )
        assert representative in names_by_pairs
        assert reversed_pairs in names_by_pairs
        assert frozenset(record["representative_product"]) == names_by_pairs[
            representative
        ]
        assert frozenset(record["reversed_product"]) == names_by_pairs[reversed_pairs]
        covered.update((representative, reversed_pairs))
    assert covered == set(names_by_pairs)

    consequences = payload["data"]["sign_parity_and_conjugation_consequences"]
    assert "equal to -6" in consequences["real_signed_necessary_condition"]
    assert "=-3" in consequences["conditional_reversal_conjugation"]
    assert "must equal -1" in consequences["conditional_unit_modulus_corollary"]
    assert "not assumed" in consequences["conditional_unit_modulus_corollary"]


def run_check(name: str, function) -> None:
    try:
        function()
    except Exception as error:
        print(f"FAIL {name}: {type(error).__name__}: {error}")
        raise
    print(f"PASS {name}")


def main() -> None:
    run_check("envelope, parameterization, and exact scope", check_envelope_parameterization_and_scope)
    run_check("cube polynomial and global Positivstellensatz", check_cube_polynomial_and_global_certificate)
    run_check("all 32 specialized branch certificates", check_all_branch_certificates)
    run_check("sign, parity, and conditional conjugation", check_sign_parity_and_conditional_conjugation)
    print("PASS")


if __name__ == "__main__":
    main()
