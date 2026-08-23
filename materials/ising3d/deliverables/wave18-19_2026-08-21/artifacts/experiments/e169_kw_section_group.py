#!/usr/bin/env python3
"""Exact symmetry audit for the branch-11111 thin-section fibration.

The previous rank-three torus lemma acts before diagonal normalization.  This
module determines which exact symmetries survive on the 13-coordinate chart and
which of those induce maps on the seven held coordinates used by e156.  It is
also imported by e171; running it directly performs the exact audit and prints
PASS only when every group check succeeds.
"""

from __future__ import annotations

import importlib.util
import itertools
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
E139_PATH = ROOT / "experiments/e139_kw_components.py"
UNIT_VALUES = (1, 2, 3, 4)
FREE_NAMES = [
    "u_px_px", "u_py_px", "u_py_py", "u_py_pz", "u_py_mz", "u_my_px", "u_my_mx",
    "u_pz_px", "u_pz_mx", "u_pz_py", "u_pz_my", "u_pz_pz", "u_mz_px", "u_mz_mx",
    "u_mz_py", "u_mz_my",
]
DIRECTION_NAMES = ("px", "mx", "py", "my", "pz", "mz")


def load_e139():
    spec = importlib.util.spec_from_file_location("e139_kw_components_for_e169", E139_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {E139_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pair_name(pair: tuple[int, int]) -> str:
    return f"u_{DIRECTION_NAMES[pair[0]]}_{DIRECTION_NAMES[pair[1]]}"


def build_construction(e139) -> dict[str, Any]:
    """Build the raw 42 equations and both thin-torus polynomial charts."""
    system = e139.full_weight_finite_system(e139.CONSTRUCTION_SHAPES, e139.BASE_ORDERS, gauge_fix=False)
    variables, equations = e139.branch_11111_equations(system)
    symbol = {str(value): value for value in variables}
    free = [symbol[name] for name in FREE_NAMES]
    non_diagonal = [symbol[name] for name in e139.NONDIAGONAL_NAMES]
    thin = e139.thin_torus_substitutions(symbol)
    diagonal = {symbol[name]: 1 for name in e139.DIAGONAL_NAMES}
    primitives = e139.primitive_numerators(equations, variables, thin, free, diagonal)
    return {
        "system": system,
        "variables": variables,
        "equations": equations,
        "labels": list(system.labels),
        "symbol": symbol,
        "free": free,
        "non_diagonal": non_diagonal,
        "thin": thin,
        "diagonal": diagonal,
        "primitives": primitives,
        "sparse42": [e139.sparse_integer_poly(expression, variables) for expression in equations],
    }


def normalized_thin_chart(e139) -> tuple[list[sp.Symbol], dict[tuple[int, int], sp.Expr]]:
    """Return all 30 raw weights on the gauge- and diagonal-normalized thin chart."""
    variables = list(sp.symbols(" ".join(e139.NONDIAGONAL_NAMES)))
    value = dict(zip(e139.NONDIAGONAL_NAMES, variables, strict=True))
    raw: dict[tuple[int, int], sp.Expr] = {}
    gauge = set(e139.DIRECTION_GAUGE_TREE)
    for pair in e139.ALLOWED_DIRECTION_PAIRS:
        name = pair_name(pair)
        if pair in gauge:
            raw[pair] = sp.Integer(1)
        elif name in value:
            raw[pair] = value[name]
        elif name in e139.DIAGONAL_NAMES:
            raw[pair] = sp.Integer(1)
    reflected = {
        "u_mx_my": -1 / value["u_my_px"],
        "u_mx_py": -1 / (value["u_my_mx"] * value["u_py_px"]),
        "u_mx_mx": sp.Integer(1),
        "u_mx_mz": -1 / (value["u_mz_px"] * value["u_pz_mx"]),
        "u_mx_pz": -1 / (value["u_mz_mx"] * value["u_pz_px"]),
        "u_mz_mz": sp.Integer(1),
        "u_my_mz": -1 / (value["u_mz_py"] * value["u_py_pz"] * value["u_pz_my"]),
        "u_my_pz": -1 / (value["u_mz_my"] * value["u_py_mz"] * value["u_pz_py"]),
        "u_my_my": sp.Integer(1),
    }
    for pair in e139.ALLOWED_DIRECTION_PAIRS:
        if pair_name(pair) in reflected:
            raw[pair] = reflected[pair_name(pair)]
    if set(raw) != set(e139.ALLOWED_DIRECTION_PAIRS):
        raise AssertionError("thin chart did not reconstruct all 30 allowed weights")
    return variables, raw


def direction_permutation(
    directions: tuple[tuple[int, int, int], ...],
    axis_permutation: tuple[int, int, int],
    signs: tuple[int, int, int],
) -> tuple[int, ...]:
    index = {tuple(vector): position for position, vector in enumerate(directions)}
    return tuple(
        index[tuple(signs[axis] * vector[axis_permutation[axis]] for axis in range(3))]
        for vector in directions
    )


def restore_tree_gauge(e139, transformed: dict[tuple[int, int], sp.Expr]) -> dict[tuple[int, int], sp.Expr]:
    """Apply the unique direction-state similarity that restores the five tree entries."""
    potential: dict[int, sp.Expr] = {0: sp.Integer(1)}
    pending = list(e139.DIRECTION_GAUGE_TREE)
    while pending:
        progress = False
        next_pending: list[tuple[int, int]] = []
        for source, target in pending:
            weight = transformed[(source, target)]
            if source in potential and target not in potential:
                potential[target] = sp.cancel(potential[source] / weight)
                progress = True
            elif target in potential and source not in potential:
                potential[source] = sp.cancel(potential[target] * weight)
                progress = True
            else:
                next_pending.append((source, target))
        if not progress and next_pending:
            raise AssertionError("gauge tree is not connected")
        pending = next_pending
    if set(potential) != set(range(6)):
        raise AssertionError("gauge tree did not determine all direction potentials")
    normalized = {
        pair: sp.factor(sp.cancel(weight * potential[pair[1]] / potential[pair[0]]))
        for pair, weight in transformed.items()
    }
    if any(sp.cancel(normalized[pair] - 1) != 0 for pair in e139.DIRECTION_GAUGE_TREE):
        raise AssertionError("tree gauge restoration failed")
    return normalized


def point_action(
    e139,
    raw: dict[tuple[int, int], sp.Expr],
    permutation: tuple[int, ...],
) -> dict[tuple[int, int], sp.Expr]:
    inverse = [0] * 6
    for old, new in enumerate(permutation):
        inverse[new] = old
    relabelled = {
        (source, target): raw[(inverse[source], inverse[target])]
        for source, target in e139.ALLOWED_DIRECTION_PAIRS
    }
    return restore_tree_gauge(e139, relabelled)


def monomial_record(expression: sp.Expr, variables: list[sp.Symbol]) -> dict[str, Any]:
    expression = sp.factor(sp.cancel(expression))
    coefficient, factors = expression.as_coeff_Mul()
    powers = factors.as_powers_dict()
    exponents = [int(powers.get(variable, 0)) for variable in variables]
    reconstructed = coefficient * sp.prod(
        variable**exponent for variable, exponent in zip(variables, exponents, strict=True)
    )
    if sp.cancel(expression - reconstructed) != 0:
        raise AssertionError(f"point action is not Laurent monomial: {expression}")
    return {"coefficient": int(coefficient), "exponents": exponents, "expression": str(expression)}


def validate_thin_action(e139, action: dict[tuple[int, int], sp.Expr]) -> None:
    out = {pair_name(pair): expression for pair, expression in action.items()}
    identities = [
        out["u_px_px"] - 1,
        out["u_py_py"] - 1,
        out["u_pz_pz"] - 1,
        out["u_mx_my"] + 1 / out["u_my_px"],
        out["u_mx_py"] + 1 / (out["u_my_mx"] * out["u_py_px"]),
        out["u_mx_mx"] - 1,
        out["u_mx_mz"] + 1 / (out["u_mz_px"] * out["u_pz_mx"]),
        out["u_mx_pz"] + 1 / (out["u_mz_mx"] * out["u_pz_px"]),
        out["u_mz_mz"] - 1,
        out["u_my_mz"] + 1 / (out["u_mz_py"] * out["u_py_pz"] * out["u_pz_my"]),
        out["u_my_pz"] + 1 / (out["u_mz_my"] * out["u_py_mz"] * out["u_pz_py"]),
        out["u_my_my"] - 1,
    ]
    if any(sp.cancel(identity) != 0 for identity in identities):
        raise AssertionError("point action does not preserve the thin chart")


def torus_record(e139, construction: dict[str, Any]) -> dict[str, Any]:
    num16: list[sp.Poly] = []
    for index in e139.SELECTED_NUMERATOR_INDICES:
        reduced = sp.cancel(construction["equations"][index].subs(construction["thin"]))
        numerator, _ = sp.fraction(reduced)
        num16.append(sp.Poly(sp.expand(numerator), *construction["free"], domain=sp.ZZ))
    rows: list[list[int]] = []
    for polynomial in num16:
        terms = polynomial.terms()
        origin = terms[0][0]
        rows.extend(
            [left - right for left, right in zip(monomial, origin, strict=True)]
            for monomial, _coefficient in terms[1:]
        )
    matrix = sp.Matrix(rows)
    kernel = [e139.integer_vector(vector) for vector in matrix.nullspace()]
    diagonal_positions = [FREE_NAMES.index(name) for name in e139.DIAGONAL_NAMES]
    diagonal_matrix = sp.Matrix(
        [[kernel[row][column] for column in diagonal_positions] for row in range(len(kernel))]
    )
    for polynomial in num16:
        for vector in kernel:
            e139.character_dots(polynomial, [vector])
    determinant = int(diagonal_matrix.det())
    return {
        "support_difference_rank": int(matrix.rank()),
        "kernel_vectors": kernel,
        "diagonal_names": list(e139.DIAGONAL_NAMES),
        "diagonal_character_matrix": [[int(entry) for entry in row] for row in diagonal_matrix.tolist()],
        "diagonal_minor_determinant": determinant,
        "diagonal_minor_abs_determinant": abs(determinant),
        "normalized_chart_stabilizer_order": 1,
        "reason": (
            "The diagonal character matrix is unimodular.  Hence fixing all three diagonal "
            "coordinates to one forces all three torus parameters to one over every unit group."
        ),
    }


def analyze_group(e139, construction: dict[str, Any]) -> dict[str, Any]:
    directions = tuple(tuple(vector) for vector in e139.DIRECTIONS_3D)
    variables13, raw = normalized_thin_chart(e139)
    shapes = set(tuple(shape) for shape in e139.CONSTRUCTION_SHAPES)
    labels = set(construction["labels"])
    full_elements: list[dict[str, Any]] = []
    catalog_elements: list[dict[str, Any]] = []
    solve_symbols = set(variables13[:6])
    held_names = list(e139.NONDIAGONAL_NAMES[6:])

    for axis_permutation in itertools.permutations(range(3)):
        mapped_shapes = {
            tuple(shape[axis_permutation[axis]] for axis in range(3)) for shape in shapes
        }
        preserves_catalog = mapped_shapes == shapes
        for signs in itertools.product((1, -1), repeat=3):
            permutation = direction_permutation(directions, axis_permutation, signs)
            element = {
                "axis_permutation": list(axis_permutation),
                "axis_signs": list(signs),
                "direction_permutation": list(permutation),
                "preserves_shape_catalog": preserves_catalog,
            }
            full_elements.append(element)
            if not preserves_catalog:
                continue
            for shape, order in labels:
                image_label = (
                    tuple(shape[axis_permutation[axis]] for axis in range(3)),
                    order,
                )
                if image_label not in labels:
                    raise AssertionError("catalog-preserving point map did not permute equation labels")
            normalized = point_action(e139, raw, permutation)
            validate_thin_action(e139, normalized)
            coordinate_action = {
                name: monomial_record(normalized[next(pair for pair in e139.ALLOWED_DIRECTION_PAIRS if pair_name(pair) == name)], variables13)
                for name in e139.NONDIAGONAL_NAMES
            }
            witnesses = []
            for held_name in held_names:
                record = coordinate_action[held_name]
                dependencies = [
                    e139.NONDIAGONAL_NAMES[index]
                    for index, exponent in enumerate(record["exponents"][:6])
                    if exponent
                ]
                if dependencies:
                    witnesses.append({
                        "held_output": held_name,
                        "solve_dependencies": dependencies,
                        "expression": record["expression"],
                    })
            element = dict(element)
            element["coordinate_action"] = coordinate_action
            element["induces_section_map"] = not witnesses
            element["fibration_obstruction_witness"] = None if not witnesses else witnesses[0]
            catalog_elements.append(element)

    inducing = [element for element in catalog_elements if element["induces_section_map"]]
    identity = [0, 1, 2, 3, 4, 5]
    if len(inducing) != 1 or inducing[0]["direction_permutation"] != identity:
        raise AssertionError("the section-fibration stabilizer was not the identity")
    torus = torus_record(e139, construction)
    checks = [
        {"name": "full_cubic_point_group_order", "passed": len(full_elements) == 48, "detail": len(full_elements)},
        {"name": "construction_catalog_point_group_order", "passed": len(catalog_elements) == 16, "detail": len(catalog_elements)},
        {
            "name": "catalog_axis_permutations",
            "passed": sorted({tuple(row["axis_permutation"]) for row in catalog_elements}) == [(0, 1, 2), (1, 0, 2)],
            "detail": sorted({tuple(row["axis_permutation"]) for row in catalog_elements}),
        },
        {"name": "torus_diagonal_minor_unimodular", "passed": torus["diagonal_minor_abs_determinant"] == 1, "detail": torus["diagonal_character_matrix"]},
        {"name": "tree_has_six_vertices_five_edges", "passed": len(e139.DIRECTION_GAUGE_TREE) == 5, "detail": list(e139.DIRECTION_GAUGE_TREE)},
        {"name": "section_action_order", "passed": len(inducing) == 1, "detail": len(inducing)},
        {"name": "section_orbit_count", "passed": len(UNIT_VALUES) ** 7 == 16_384, "detail": len(UNIT_VALUES) ** 7},
    ]
    return {
        "raw_direction_state_gauge": {
            "action": "U(i,j) -> a_i^(-1) U(i,j) a_j",
            "effective_group": "(G_m)^6 / diagonal G_m",
            "tree": [list(pair) for pair in e139.DIRECTION_GAUGE_TREE],
            "tree_connected": True,
            "tree_chart_residual_order": 1,
            "sign_redundancy_after_tree_normalization": "none",
        },
        "diagonal_torus": torus,
        "cubic_point_group": {
            "full_lattice_group_order": len(full_elements),
            "construction_catalog_group_order": len(catalog_elements),
            "construction_catalog_group_description": (
                "all eight independent axis sign flips, semidirect the identity or x-y axis swap"
            ),
            "excluded_axis_permutations": [
                list(permutation)
                for permutation in itertools.permutations(range(3))
                if permutation not in ((0, 1, 2), (1, 0, 2))
            ],
            "exclusion_reason": (
                "the 14-shape construction catalog is not invariant under the other four axis permutations"
            ),
            "elements": catalog_elements,
        },
        "section_action": {
            "held_names": held_names,
            "solve_names": list(e139.SLICE6_NAMES),
            "criterion": "every transformed held Laurent monomial is independent of all six solve coordinates",
            "inducing_element_count": len(inducing),
            "group_order": 1,
            "group_description": "identity only",
            "orbit_count": 16_384,
            "orbit_sizes": {"1": 16_384},
            "quotient_reduction_factor": 1,
        },
        "checks": checks,
    }


def main() -> None:
    e139 = load_e139()
    construction = build_construction(e139)
    record = analyze_group(e139, construction)
    if not all(check["passed"] for check in record["checks"]):
        raise AssertionError(record["checks"])
    print(
        "catalog point group order=16; section stabilizer order=1; "
        "orbit count=16384; quotient reduction factor=1"
    )
    print("PASS")


if __name__ == "__main__":
    main()
