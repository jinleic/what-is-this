#!/usr/bin/env python3
"""Exact characteristic-zero structure certificate for the open 3x3 Ising layer.

For A=sum_v X_v and B=sum_{uv in E(3x3)} Z_u Z_v this script combines an
exact-Q, all-sector linked upper container with literal nested Pauli-word
closures at one good prime.  The modular computation is used only for lower
bounds; every containment, projector, linkage, and scalar witness is exact
rational arithmetic.
"""
from __future__ import annotations

import argparse
import gc
import importlib.util
import itertools
import json
import math
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e87_algebra_3x3.py"
RESULT = ROOT / "results" / "algebra_structure" / "char0_3x3.json"
E20_PATH = ROOT / "experiments" / "e20_algebra_growth.py"
ROWS = COLS = 3
N = ROWS * COLS
GOOD_PRIME = 2_147_483_647
CHARACTERS = tuple(itertools.product((0, 1), repeat=3))
BONDS = tuple(
    (r * COLS + c, (r + 1) * COLS + c)
    for r in range(ROWS - 1)
    for c in range(COLS)
) + tuple(
    (r * COLS + c, r * COLS + c + 1)
    for r in range(ROWS)
    for c in range(COLS - 1)
)
DIMENSION = 8034
DERIVED_DIMENSION = 8033
FACTOR_SPECS = (
    ("A32", "000", 33, 1088, 32),
    ("A47", "000", 48, 2303, 47),
    ("A58", "010", 59, 3480, 58),
    ("A2", "110", 3, 8, 2),
    ("A15", "110", 16, 255, 15),
    ("A29", "110", 30, 899, 29),
)
EXPECTED_SECTOR_DIMS = {
    "000": 84,
    "001": 84,
    "010": 60,
    "011": 60,
    "100": 60,
    "101": 60,
    "110": 52,
    "111": 52,
}
EXPECTED_BLOCK_DIMS = {
    "000": [1, 2, 33, 48],
    "001": [1, 2, 33, 48],
    "010": [1, 59],
    "011": [1, 59],
    "100": [1, 59],
    "101": [1, 59],
    "110": [1, 2, 3, 16, 30],
    "111": [1, 2, 3, 16, 30],
}
EXPECTED_COMMUTANTS = {
    "000": (7, 4),
    "001": (7, 4),
    "010": (2, 2),
    "011": (2, 2),
    "100": (2, 2),
    "101": (2, 2),
    "110": (8, 5),
    "111": (8, 5),
}
PER_SECTOR_IMAGE_DIMS = {
    "000": 3392,
    "001": 3392,
    "010": 3481,
    "011": 3481,
    "100": 3481,
    "101": 3481,
    "110": 1163,
    "111": 1163,
}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def timed(label: str, timings: dict, memory: dict, function, *args):
    started = time.monotonic()
    before = rss_bytes()
    value = function(*args)
    timings[label] = round(time.monotonic() - started, 6)
    memory[label] = {"rss_before_bytes": before, "rss_peak_bytes": rss_bytes()}
    return value


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(N))


def d4_permutations() -> tuple[tuple[int, ...], ...]:
    identity = tuple(range(N))
    row_flip = tuple((ROWS - 1 - r) * COLS + c for r in range(ROWS) for c in range(COLS))
    col_flip = tuple(r * COLS + (COLS - 1 - c) for r in range(ROWS) for c in range(COLS))
    transpose = tuple(c * COLS + r for r in range(ROWS) for c in range(COLS))
    seen = {identity}
    queue = [identity]
    for current in queue:
        for generator in (row_flip, col_flip, transpose):
            candidate = compose(generator, current)
            if candidate not in seen:
                seen.add(candidate)
                queue.append(candidate)
    return tuple(sorted(seen))


def site_symmetries() -> tuple[tuple[int, int, int, tuple[int, ...]], ...]:
    """The diagonalizable C2(row)xC2(column)xC2(layer-parity subgroup."""
    return tuple(
        (
            flip_rows,
            flip_cols,
            spin_flip,
            tuple(
                (ROWS - 1 - r if flip_rows else r) * COLS
                + (COLS - 1 - c if flip_cols else c)
                for r in range(ROWS)
                for c in range(COLS)
            ),
        )
        for flip_rows, flip_cols, spin_flip in CHARACTERS
    )


def sector_with_gram(character: tuple[int, int, int]) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Exact C2^3 character sector and inherited physical Gram form."""
    size = 1 << N
    symmetries = site_symmetries()

    def act(configuration: int, symmetry: tuple[int, int, int, tuple[int, ...]]) -> int:
        output = 0
        for old, new in enumerate(symmetry[3]):
            if (configuration >> old) & 1:
                output |= 1 << new
        return output ^ ((size - 1) if symmetry[2] else 0)

    seen: set[int] = set()
    representatives: list[int] = []
    orbit_vectors: list[dict[int, int]] = []
    for configuration in range(size):
        if configuration in seen:
            continue
        seen.update(act(configuration, symmetry) for symmetry in symmetries)
        values: dict[int, int] = {}
        for symmetry in symmetries:
            sign = -1 if sum(x * y for x, y in zip(character, symmetry[:3])) & 1 else 1
            image = act(configuration, symmetry)
            values[image] = values.get(image, 0) + sign
        if not any(values.values()):
            continue
        representative = min(configuration for configuration, coefficient in values.items() if coefficient)
        scale = values[representative]
        representatives.append(representative)
        orbit_vectors.append({configuration: coefficient // scale for configuration, coefficient in values.items() if coefficient})

    dimension = len(representatives)
    a = sp.zeros(dimension)
    b = sp.zeros(dimension)
    gram = sp.zeros(dimension)
    for j, orbit in enumerate(orbit_vectors):
        gram[j, j] = sum(coefficient * coefficient for coefficient in orbit.values())
        image: dict[int, int] = {}
        for configuration, coefficient in orbit.items():
            for site in range(N):
                target = configuration ^ (1 << site)
                image[target] = image.get(target, 0) + coefficient
        for i, representative in enumerate(representatives):
            a[i, j] = image.get(representative, 0)
        b[j, j] = sum(
            1 if ((representatives[j] >> u) & 1) == ((representatives[j] >> v) & 1) else -1
            for u, v in BONDS
        )
    assert a.T * gram == gram * a
    assert b.T * gram == gram * b
    return a, b, gram


def primitive_integer_matrix(matrix: sp.MatrixBase) -> sp.Matrix:
    entries = [sp.Rational(value) for value in matrix]
    denominator = sp.ilcm(*(value.q for value in entries)) if entries else 1
    integers = [int(value * denominator) for value in entries]
    content = math.gcd(*(abs(value) for value in integers)) if integers else 1
    if content:
        integers = [value // content for value in integers]
    first = next((value for value in integers if value), 1)
    if first < 0:
        integers = [-value for value in integers]
    return sp.Matrix(matrix.rows, matrix.cols, integers)


def sparse_matrix_certificate(matrix: sp.MatrixBase, store_entries: bool = True) -> dict:
    integer = primitive_integer_matrix(matrix)
    entries = (
        [
            [i, j, int(integer[i, j])]
            for i in range(integer.rows)
            for j in range(integer.cols)
            if integer[i, j]
        ]
        if store_entries
        else None
    )
    return {
        "shape": [integer.rows, integer.cols],
        "entries": entries,
        "nonzero_entries": sum(value != 0 for value in integer),
        "sum_abs_entries": sum(abs(int(value)) for value in integer),
        "weighted_mod_1000000007": sum(
            (i + 1) * (j + 1) * int(integer[i, j])
            for i in range(integer.rows)
            for j in range(integer.cols)
        )
        % 1_000_000_007,
    }


def sparse_nullspace(equations: dict[int, dict[int, object]], rows: int, columns: int) -> list[list[object]]:
    if columns == 0:
        return []
    if rows == 0:
        return [[QQ.one if i == j else QQ.zero for j in range(columns)] for i in range(columns)]
    matrix = DomainMatrix(equations, (rows, columns), QQ)
    return matrix.nullspace().to_Matrix().tolist()


def intertwiner_basis(
    target: tuple[sp.Matrix, sp.Matrix],
    source: tuple[sp.Matrix, sp.Matrix],
    signs: tuple[int, int] = (1, 1),
) -> list[sp.Matrix]:
    """Exact maps T satisfying target_i*T=sign_i*T*source_i."""
    target_a, target_b = target
    source_a, source_b = source
    sign_a, sign_b = signs
    target_dimension, source_dimension = target_a.rows, source_a.rows
    diagonal = (
        target_b == sp.diag(*(target_b[i, i] for i in range(target_dimension)))
        and source_b == sp.diag(*(source_b[i, i] for i in range(source_dimension)))
    )
    pairs = [
        (i, j)
        for i in range(target_dimension)
        for j in range(source_dimension)
        if not diagonal or target_b[i, i] == sign_b * source_b[j, j]
    ]
    variable = {pair: column for column, pair in enumerate(pairs)}
    equations: dict[int, dict[int, object]] = {}
    row_index = 0
    generators = ((target_a, source_a, sign_a),) if diagonal else (
        (target_a, source_a, sign_a),
        (target_b, source_b, sign_b),
    )
    for target_generator, source_generator, sign in generators:
        for i in range(target_dimension):
            for j in range(source_dimension):
                row: dict[int, object] = {}
                for k in range(target_dimension):
                    column = variable.get((k, j))
                    coefficient = target_generator[i, k]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) + coefficient
                for k in range(source_dimension):
                    column = variable.get((i, k))
                    coefficient = source_generator[k, j]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) - sign * coefficient
                row = {column: coefficient for column, coefficient in row.items() if coefficient}
                if row:
                    equations[row_index] = {column: QQ.convert(coefficient) for column, coefficient in row.items()}
                    row_index += 1
    answer = []
    for vector in sparse_nullspace(equations, row_index, len(pairs)):
        matrix = sp.zeros(target_dimension, source_dimension)
        for coefficient, (i, j) in zip(vector, pairs):
            matrix[i, j] = coefficient
        assert target_a * matrix == sign_a * matrix * source_a
        assert target_b * matrix == sign_b * matrix * source_b
        answer.append(primitive_integer_matrix(matrix))
    return answer


def commutant_center(matrices: list[sp.Matrix]) -> list[sp.Matrix]:
    dimension = matrices[0].rows
    count = len(matrices)
    commutators = [
        [matrices[i] * matrices[j] - matrices[j] * matrices[i] for j in range(count)]
        for i in range(count)
    ]
    equations: dict[int, dict[int, object]] = {}
    row_index = 0
    for right in range(count):
        for i in range(dimension):
            for j in range(dimension):
                row = {left: commutators[left][right][i, j] for left in range(count) if commutators[left][right][i, j]}
                if row:
                    equations[row_index] = {column: QQ.convert(value) for column, value in row.items()}
                    row_index += 1
    answer = []
    for vector in sparse_nullspace(equations, row_index, count):
        matrix = sp.zeros(dimension)
        for coefficient, item in zip(vector, matrices):
            matrix += coefficient * item
        answer.append(primitive_integer_matrix(matrix))
    return answer


def central_projectors(actions: tuple[sp.Matrix, sp.Matrix]) -> tuple[list[sp.Matrix], list[sp.Matrix], tuple[int, ...], tuple[tuple[int, sp.Matrix, object], ...]]:
    a, b = actions
    commutant = intertwiner_basis(actions, actions)
    commutative = all(
        commutant[i] * commutant[j] == commutant[j] * commutant[i]
        for i in range(len(commutant))
        for j in range(i)
    )
    center = commutant if commutative else commutant_center(commutant)
    for weights in itertools.product((1, 2, 3), repeat=len(center)):
        if all(weight == 1 for weight in weights):
            continue
        operator = sp.zeros(a.rows)
        for weight, matrix in zip(weights, center):
            operator += weight * matrix
        roots = sp.roots(sp.factor(operator.charpoly().as_expr()))
        if len(roots) != len(center) or not all(root.is_Rational for root in roots):
            continue
        identity = sp.eye(a.rows)
        projectors = []
        for value in roots:
            projector = identity
            for other in roots:
                if other != value:
                    projector = projector * (operator - other * identity) / (value - other)
            assert projector * projector == projector
            assert projector * a == a * projector
            assert projector * b == b * projector
            projectors.append((projector.rank(), projector, value))
        projectors = tuple(sorted(projectors, key=lambda item: (item[0], str(item[2]))))
        assert sum((projector for _rank, projector, _value in projectors), sp.zeros(a.rows)) == identity
        return commutant, center, weights, projectors
    raise AssertionError("no rational central separator")


def restrict_to_projector(actions: tuple[sp.Matrix, sp.Matrix], projector: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix]:
    inclusion = sp.Matrix.hstack(*projector.columnspace())
    left_inverse = (inclusion.T * inclusion).inv() * inclusion.T
    assert inclusion * left_inverse * projector == projector
    restricted = tuple(sp.simplify(left_inverse * action * inclusion) for action in actions)
    b = restricted[1]
    assert b == sp.diag(*(b[i, i] for i in range(b.rows)))
    return restricted


def traceless_actions(actions: tuple[sp.Matrix, sp.Matrix]) -> tuple[sp.Matrix, sp.Matrix]:
    dimension = actions[0].rows
    identity = sp.eye(dimension)
    return tuple(action - action.trace() * identity / dimension for action in actions)


def invariant_form_nullity(actions: tuple[sp.Matrix, sp.Matrix], alternating: bool) -> tuple[int, int]:
    """Exact nullity of X^T F+F X=0 with F symmetric/alternating."""
    a, b = actions
    dimension = a.rows
    assert b == sp.diag(*(b[i, i] for i in range(dimension)))
    sign = -1 if alternating else 1
    pairs = [
        (i, j)
        for i in range(dimension)
        for j in range(i + (1 if alternating else 0), dimension)
        if b[i, i] + b[j, j] == 0
    ]
    variable = {pair: column for column, pair in enumerate(pairs)}

    def form_variable(i: int, j: int) -> tuple[int | None, int]:
        if i <= j:
            return variable.get((i, j)), 1
        return variable.get((j, i)), sign

    equations: dict[int, dict[int, object]] = {}
    row_index = 0
    for i in range(dimension):
        for j in range(dimension):
            row: dict[int, object] = {}
            for k in range(dimension):
                column, orientation = form_variable(k, j)
                coefficient = orientation * a[k, i]
                if column is not None and coefficient:
                    row[column] = row.get(column, 0) + coefficient
                column, orientation = form_variable(i, k)
                coefficient = orientation * a[k, j]
                if column is not None and coefficient:
                    row[column] = row.get(column, 0) + coefficient
            row = {column: coefficient for column, coefficient in row.items() if coefficient}
            if row:
                equations[row_index] = {column: QQ.convert(coefficient) for column, coefficient in row.items()}
                row_index += 1
    return len(pairs), len(sparse_nullspace(equations, row_index, len(pairs)))


def full_rank_witness(matrices: list[sp.Matrix]) -> sp.Matrix | None:
    if not matrices:
        return None
    target = min(matrices[0].shape)
    for weights in itertools.product((1, 2, 3), repeat=len(matrices)):
        candidate = sum((weight * matrix for weight, matrix in zip(weights, matrices)), sp.zeros(*matrices[0].shape))
        if candidate.rank() == target:
            return primitive_integer_matrix(candidate)
    return None


def component_by_dimension(components: dict[str, list[dict]], label: str, dimension: int) -> dict:
    matches = [part for part in components[label] if part["dimension"] == dimension]
    if len(matches) != 1:
        raise AssertionError((label, dimension, len(matches)))
    return matches[0]


def intertwiner_record(source_label: str, target_label: str, dimension: int, source: tuple[sp.Matrix, sp.Matrix], target: tuple[sp.Matrix, sp.Matrix], signs: tuple[int, int]) -> tuple[dict, sp.Matrix]:
    basis = intertwiner_basis(target, source, signs)
    witness = full_rank_witness(basis)
    if witness is None:
        raise AssertionError((source_label, target_label, dimension, signs))
    return {
        "source_sector": source_label,
        "target_sector": target_label,
        "component_dimension_Q": dimension,
        "signs": list(signs),
        "space_dimension_Q": len(basis),
        "rank_Q": witness.rank(),
        "witness": sparse_matrix_certificate(witness),
    }, witness


def exact_q_certificate() -> dict:
    sectors = {key: sector_with_gram(key) for key in CHARACTERS}
    components: dict[str, list[dict]] = {}
    sector_records: dict[str, dict] = {}

    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, gram = sectors[key]
        commutant, center, weights, projectors = central_projectors((a, b))
        blocks = []
        for dimension, projector, eigenvalue in projectors:
            actions = restrict_to_projector((a, b), projector)
            trace_a, trace_b = actions[0].trace(), actions[1].trace()
            traceless = traceless_actions(actions)
            blocks.append(
                {
                    "dimension": dimension,
                    "projector": projector,
                    "eigenvalue": eigenvalue,
                    "actions": actions,
                    "traceless": traceless,
                    "trace_a": trace_a,
                    "trace_b": trace_b,
                    "scalar_a": trace_a / dimension,
                    "traceless_zero": traceless[0] == sp.zeros(dimension) and traceless[1] == sp.zeros(dimension),
                }
            )
        assert a.rows == EXPECTED_SECTOR_DIMS[label]
        assert [block["dimension"] for block in blocks] == EXPECTED_BLOCK_DIMS[label]
        assert (len(commutant), len(center)) == EXPECTED_COMMUTANTS[label]
        components[label] = blocks
        sector_records[label] = {
            "claim_tag": "[COMPUTATION]",
            "module_dimension_Q": a.rows,
            "physical_gram": sparse_matrix_certificate(gram),
            "commutant_dimension_Q": len(commutant),
            "commutant_center_dimension_Q": len(center),
            "central_separator_weights": list(weights),
            "isotypic_blocks": [
                {
                    "dimension_Q": block["dimension"],
                    "central_eigenvalue_Q": str(block["eigenvalue"]),
                    "projector": sparse_matrix_certificate(block["projector"]),
                    "trace_A_Q": str(block["trace_a"]),
                    "trace_B_Q": str(block["trace_b"]),
                    "scalar_A_Q": str(block["scalar_a"]),
                    "traceless_actions_zero": block["traceless_zero"],
                }
                for block in blocks
            ],
            "sector_image_container_dimension_Q": PER_SECTOR_IMAGE_DIMS[label],
            "certification_scope": "exact-Q projector/trace/linkage data; equality of the displayed image container follows from the faithful global 8034-dimensional squeeze",
        }

    factor_controls: dict[str, dict] = {}
    for factor, label, dimension, factor_dimension, rank in FACTOR_SPECS:
        block = component_by_dimension(components, label, dimension)
        commutant_dimension = len(intertwiner_basis(block["traceless"], block["traceless"]))
        symmetric_variables, symmetric_nullity = invariant_form_nullity(block["traceless"], False)
        alternating_variables, alternating_nullity = invariant_form_nullity(block["traceless"], True)
        assert commutant_dimension == 1
        assert symmetric_nullity == alternating_nullity == 0
        factor_controls[factor] = {
            "module_dimension_Q": dimension,
            "factor_dimension_Q": factor_dimension,
            "rank": rank,
            "traceless_commutant_dimension_Q": 1,
            "invariant_symmetric_form_variables_Q": symmetric_variables,
            "invariant_symmetric_form_nullity_Q": 0,
            "invariant_alternating_form_variables_Q": alternating_variables,
            "invariant_alternating_form_nullity_Q": 0,
            "classification": f"sl({dimension},Q), absolute type {factor}",
        }

    # The scalar part of A, assembled before primitive normalization, is one
    # global exact central witness.  B has zero trace on every isotypic block.
    scalar_pieces = []
    residual_checks = []
    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, _gram = sectors[key]
        scalar = sp.zeros(a.rows)
        for block in components[label]:
            scalar += block["scalar_a"] * block["projector"]
            residual_a = block["actions"][0] - block["scalar_a"] * sp.eye(block["dimension"])
            residual_checks.append(residual_a.trace() == 0 and block["actions"][1].trace() == 0)
            if block["traceless_zero"]:
                assert residual_a == sp.zeros(block["dimension"])
                assert block["actions"][1] == sp.zeros(block["dimension"])
        assert scalar * a == a * scalar
        assert scalar * b == b * scalar
        scalar_pieces.append(scalar)
    central = sp.diag(*scalar_pieces)
    central_integer = primitive_integer_matrix(central)
    central_pivot = next(index for index, value in enumerate(central_integer) if value)
    central_scale = sp.Rational(list(central)[central_pivot], list(central_integer)[central_pivot])
    assert central == central_scale * central_integer
    full_a = sp.diag(*(sectors[key][0] for key in CHARACTERS))
    full_b = sp.diag(*(sectors[key][1] for key in CHARACTERS))
    assert central != sp.zeros(1 << N)
    assert central * full_a == full_a * central
    assert central * full_b == full_b * central
    assert all(residual_checks)

    p_signed_records = []
    for source_label, target_label, dimensions in (
        ("000", "001", (33, 48)),
        ("010", "011", (59,)),
        ("100", "101", (59,)),
        ("110", "111", (3, 16, 30)),
    ):
        for dimension in dimensions:
            source = component_by_dimension(components, source_label, dimension)["traceless"]
            target = component_by_dimension(components, target_label, dimension)["traceless"]
            record, p_witness = intertwiner_record(source_label, target_label, dimension, source, target, (-1, 1))
            outer_basis = intertwiner_basis(source, (source[0].T, source[1].T), (1, -1))
            outer_witness = full_rank_witness(outer_basis)
            if outer_witness is None:
                raise AssertionError((source_label, dimension, "outer"))
            assert source[0] * outer_witness == outer_witness * source[0].T
            assert source[1] * outer_witness == -outer_witness * source[1].T
            dual_witness = primitive_integer_matrix(p_witness * outer_witness)
            assert target[0] * dual_witness == -dual_witness * source[0].T
            assert target[1] * dual_witness == -dual_witness * source[1].T
            record["outer_sign_witness"] = {
                "space_dimension_Q": len(outer_basis),
                "rank_Q": outer_witness.rank(),
                "witness": sparse_matrix_certificate(outer_witness),
            }
            record["dual_witness"] = sparse_matrix_certificate(dual_witness)
            p_signed_records.append(record)

    d4_records = []
    for source_label, target_label in (("010", "100"), ("011", "101")):
        source = component_by_dimension(components, source_label, 59)["traceless"]
        target = component_by_dimension(components, target_label, 59)["traceless"]
        record, _witness = intertwiner_record(source_label, target_label, 59, source, target, (1, 1))
        d4_records.append(record)

    factor_dimensions = [item[3] for item in FACTOR_SPECS]
    factor_ranks = [item[4] for item in FACTOR_SPECS]
    assert sum(factor_dimensions) == DERIVED_DIMENSION
    assert sum(factor_ranks) == 183
    linked_container = {
        "type_with_center": "centre(1) + A32 + A47 + A58 + A2 + A15 + A29",
        "dimension_Q": DIMENSION,
        "derived_dimension_Q": DERIVED_DIMENSION,
        "factor_dimensions": factor_dimensions,
        "factor_ranks": factor_ranks,
        "levi_rank": sum(factor_ranks),
        "per_sector_image_dimensions_Q": PER_SECTOR_IMAGE_DIMS,
        "dimension_identity": "1 + (33^2-1) + (48^2-1) + (59^2-1) + (3^2-1) + (16^2-1) + (30^2-1) = 8034",
        "linkage_reason": "The D4 transpose witnesses identify 010 with 100 and 011 with 101. The P-sign witnesses, together with the exact outer witnesses A0*C=C*A0^T and B*C=-C*B^T, identify every P-odd nontrivial block with the contragredient of its P-even partner. A-z and B are block-traceless, while z is the one common scalar line. Hence the full eight-sector faithful image lies in one linked split container, not a sum of eight independent sector algebras.",
    }
    return {
        "claim_tag": "[THEOREM]",
        "sector_certificates_Q": sector_records,
        "factor_controls_Q": factor_controls,
        "intertwiners_Q": {
            "P_signed": p_signed_records,
            "D4_direct": d4_records,
        },
        "central_witness_global_integer": sparse_matrix_certificate(central_integer),
        "central_projection_equals_integer_times": str(central_scale),
        "central_projection_generator": "A",
        "centrality_verified_exactly": True,
        "central_membership_proof": "The exact central matrix z is the blockwise scalar projection of A. Thus A-z and B are in the block-traceless linked semisimple container S. The derived raw-word lower certificate has dimension dim(S)=8033, so g'=S. Hence A-z belongs to g', and z=A-(A-z) belongs to g. The faithful factor projection has kernel Qz.",
        "containing_algebra_Q": linked_container,
        "direct_sector_matrix_word_minors": {
            "claim_tag": "[UNRESOLVED]",
            "status": "not used",
            "reason": "A direct p=65521 dense literal matrix-word minor for the 84-dimensional sector was stopped at a 3600-second wall before a certificate completed. Per-sector equality is instead a consequence of the faithful all-sector exact-Q squeeze and is not represented as a direct single-sector selected minor.",
        },
    }


def modular_word_certificate() -> dict:
    """Actual nested Pauli words modulo a good prime, through saturation."""
    e20 = load_module("_e20_for_e87", E20_PATH)
    model = e20.OrbitModel(ROWS, COLS, False)
    support = model.enumerate_support(100_000)
    assert len(model.group) == 8
    assert len(support) == 8739

    def close(derived: bool) -> dict:
        worker = e20.DenseOrbitEngine(model, GOOD_PRIME, support, None, None)

        def apply_raw(action_index: int, word: tuple[int, np.ndarray]) -> tuple[int, np.ndarray]:
            source_block, row = word
            target_block = source_block ^ (3 if action_index == 0 else 1)
            sources, signs = worker.actions[action_index][source_block]
            result = np.zeros(worker.block_sizes[target_block], dtype=np.int64)
            for source, sign in zip(sources, signs):
                selected = source >= 0
                result[selected] += sign[selected].astype(np.int64) * row[source[selected]].astype(np.int64)
            return target_block, result % GOOD_PRIME

        if derived:
            seed = apply_raw(0, worker.dense_seed(model.zz_terms))
            seeds = [(seed, [-1, 0, 2])]
        else:
            seeds = [
                (worker.dense_seed(model.x_terms), [-1, 0, 1]),
                (worker.dense_seed(model.zz_terms), [-1, 1, 1]),
            ]
        words: list[tuple[int, np.ndarray]] = []
        tree: list[list[int]] = []
        frontier: list[int] = []
        for word, recipe in seeds:
            assert worker.reduce_add(word) is not None
            words.append((word[0], np.asarray(word[1], dtype=np.int64) % GOOD_PRIME))
            tree.append(recipe)
            frontier.append(len(words) - 1)
        while frontier:
            next_frontier: list[int] = []
            for parent in frontier:
                for action_index in (0, 1):
                    candidate = apply_raw(action_index, words[parent])
                    if worker.reduce_add(candidate) is not None:
                        words.append(candidate)
                        tree.append([parent, action_index, tree[parent][2] + 1])
                        next_frontier.append(len(words) - 1)
            frontier = next_frontier
        return {
            "rank_Fp": len(words),
            "maximum_depth": max(recipe[2] for recipe in tree),
            "grading_block_ranks": list(worker.row_counts),
            "tree_checksum_mod_1000000007": sum(
                (index + 1) * (recipe[0] + 2) * (recipe[1] + 2) * (recipe[2] + 1)
                for index, recipe in enumerate(tree)
            )
            % 1_000_000_007,
            "tree": tree,
            "basis_used_bytes": worker.used_basis_bytes,
            "basis_allocated_bytes": worker.allocated_basis_bytes,
        }

    full = close(False)
    gc.collect()
    derived = close(True)
    assert full["rank_Fp"] == DIMENSION
    assert full["grading_block_ranks"] == [1962, 1963, 2147, 1962]
    assert derived["rank_Fp"] == DERIVED_DIMENSION
    assert derived["grading_block_ranks"] == [1962, 1963, 2146, 1962]
    return {
        "claim_tag": "[COMPUTATION]",
        "prime": GOOD_PRIME,
        "arithmetic": "exact F_p; each accepted row is a literal nested [A,w]/2 or [B,w]/2 Pauli word",
        "support_orbits": len(support),
        "full_closure": full,
        "derived_ideal": derived,
        "rational_implication": "a nonzero modular pivot minor of literal integer Pauli words is nonzero over Q, so the displayed ranks are rigorous characteristic-zero lower bounds",
    }


def build_artifact() -> dict:
    timings: dict[str, float] = {}
    memory: dict[str, dict] = {}
    started = time.monotonic()
    exact = timed("exact_Q_sector_container", timings, memory, exact_q_certificate)
    modular = timed("good_prime_raw_word_closures", timings, memory, modular_word_certificate)
    timings["total"] = round(time.monotonic() - started, 6)

    full_lower = modular["full_closure"]["rank_Fp"]
    derived_lower = modular["derived_ideal"]["rank_Fp"]
    container = exact["containing_algebra_Q"]
    full_upper = container["dimension_Q"]
    derived_upper = container["derived_dimension_Q"]
    assert full_lower == full_upper == DIMENSION
    assert derived_lower == derived_upper == DERIVED_DIMENSION

    d4 = d4_permutations()
    bond_set = {tuple(sorted(bond)) for bond in BONDS}
    assert len(d4) == 8
    assert all({tuple(sorted((permutation[u], permutation[v]))) for u, v in BONDS} == bond_set for permutation in d4)

    structural_killing = {
        "claim_tag": "[THEOREM]",
        "rank_Q": DERIVED_DIMENSION,
        "nullity_Q": 1,
        "reconstruction": {
            "zero_center_block_dimension": 1,
            "simple_blocks": [
                {"absolute_type": factor, "dimension": dimension, "rational_form": f"sl({module_dimension},Q)"}
                for factor, _sector, module_dimension, dimension, _rank in FACTOR_SPECS
            ],
            "reason": "the exact linked decomposition is Qz plus split simple sl_d(Q) blocks; the Killing form vanishes on Qz and is nondegenerate on each simple block",
        },
        "direct_adjoint_matrix": {
            "claim_tag": "[UNRESOLVED]",
            "materialized": False,
            "distinction": "the rank is structurally reconstructed from the exact direct sum; no dense 8034x8034 adjoint-trace array or minor was materialized",
        },
    }
    data = {
        "scope": "the single finite open 3x3 layer only; no all-sizes or thermodynamic assertion",
        "field": "Q, with one exact good-prime computation used only for lower bounds",
        "symmetry": {
            "claim_tag": "[COMPUTATION]",
            "d4_order": len(d4),
            "d4_permutations": [list(permutation) for permutation in d4],
            "layer_parity_order": 2,
            "abelian_sector_subgroup": "C2(row reflection) x C2(column reflection) x C2(global spin flip)",
            "sector_dimensions_Q": EXPECTED_SECTOR_DIMS,
            "fixed_generators": True,
            "D4_sector_map_under_transpose": {
                "000": "000", "001": "001", "010": "100", "011": "101",
                "100": "010", "101": "011", "110": "110", "111": "111",
            },
        },
        "dimension_Q": {
            "claim_tag": "[THEOREM]",
            "value": DIMENSION,
            "upper_bound_Q": full_upper,
            "lower_bound_from_Fp_word_minor": full_lower,
        },
        "derived_algebra": {
            "claim_tag": "[THEOREM]",
            "dimension_Q": DERIVED_DIMENSION,
            "semisimple": True,
            "type": "A32 + A47 + A58 + A2 + A15 + A29",
            "upper_bound_Q": derived_upper,
            "lower_bound_from_Fp_word_minor": derived_lower,
        },
        "factor_map": {
            "claim_tag": "[THEOREM]",
            "domain_dimension_Q": DIMENSION,
            "joint_rank_Q": DERIVED_DIMENSION,
            "kernel_dimension_Q": 1,
            "kernel_equals_center": True,
            "central_witness_global_integer": exact["central_witness_global_integer"],
            "central_projection_equals_integer_times": exact["central_projection_equals_integer_times"],
            "proof": exact["central_membership_proof"],
        },
        "solvable_radical": {
            "claim_tag": "[THEOREM]",
            "dimension_Q": 1,
            "equals_center": True,
            "proof": "the quotient by the exact central kernel is the displayed split semisimple direct sum; every solvable ideal maps trivially to it",
        },
        "Levi_type": {
            "claim_tag": "[THEOREM]",
            "type": "A32 + A47 + A58 + A2 + A15 + A29",
            "dimension_Q": DERIVED_DIMENSION,
            "rank": 183,
            "factor_dimensions": [item[3] for item in FACTOR_SPECS],
            "factor_ranks": [item[4] for item in FACTOR_SPECS],
            "classification_method": "exact block-traceless sl_d(Q) containment plus the faithful dimension squeeze; symmetric and alternating invariant-form nullities are exactly zero on every representative natural module",
        },
        "Killing": structural_killing,
        "exact_Q_certificate": exact,
        "modular_word_certificate": modular,
        "resource_observations": {
            "claim_tag": "[COMPUTATION]",
            "direct_adjoint_collection_int64_bytes": DIMENSION ** 3 * 8,
            "direct_adjoint_collection_gib": DIMENSION ** 3 * 8 / 1024 ** 3,
            "direct_sector_matrix_word_minor": exact["direct_sector_matrix_word_minors"],
            "resolution": "the exact all-sector projector/linkage container plus raw Pauli lower certificate avoids both dense adjoint arrays and the stalled direct sector-matrix minor route",
        },
        "timings_seconds": timings,
        "memory_observations": memory,
    }
    checks = [
        {"name": "D4_layer_parity_sectorization", "passed": data["symmetry"]["d4_order"] == 8 and sum(EXPECTED_SECTOR_DIMS.values()) == 512, "detail": "D4 x C2 symmetry with exact C2^3 character sectors"},
        {"name": "same_word_dimension_lower", "passed": full_lower == DIMENSION, "detail": "8034 literal nested Pauli words are independent at the good prime"},
        {"name": "exact_Q_linked_upper", "passed": full_upper == DIMENSION, "detail": "faithful all-sector linked centre plus six split-A factor container"},
        {"name": "derived_dimension_squeeze", "passed": derived_lower == derived_upper == DERIVED_DIMENSION, "detail": "8033 derived raw words meet the exact block-traceless upper bound"},
        {"name": "factor_dimensions", "passed": sum(item[3] for item in FACTOR_SPECS) == DERIVED_DIMENSION, "detail": "A32 + A47 + A58 + A2 + A15 + A29"},
        {"name": "factor_rank", "passed": sum(item[4] for item in FACTOR_SPECS) == 183, "detail": "32+47+58+2+15+29"},
        {"name": "central_factor_kernel", "passed": data["factor_map"]["kernel_dimension_Q"] == 1 and data["factor_map"]["kernel_equals_center"], "detail": "A-scalar witness z and faithful semisimple quotient"},
        {"name": "radical_equals_center", "passed": data["solvable_radical"]["equals_center"], "detail": "one-dimensional central kernel of a split semisimple quotient"},
        {"name": "structural_Killing_rank", "passed": structural_killing["rank_Q"] == DERIVED_DIMENSION and structural_killing["nullity_Q"] == 1, "detail": "structural direct-sum reconstruction; no dense adjoint matrix claimed"},
    ]
    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": "faithful D4/layer-parity sector upper container, exact commutant-center projectors and signed/dual linkages, plus literal raw Pauli-word good-prime lower closures",
            "python": platform.python_version(),
        },
        "data": data,
        "checks": checks,
    }


def run(write: bool = True) -> dict:
    artifact = build_artifact()
    if write:
        RESULT.parent.mkdir(parents=True, exist_ok=True)
        RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return artifact


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args(argv)
    try:
        artifact = run(not args.no_write)
        for check in artifact["checks"]:
            print(("PASS" if check["passed"] else "FAIL") + ": " + check["name"])
        passed = all(check["passed"] for check in artifact["checks"])
        print("PASS" if passed else "FAIL")
        return 0 if passed else 1
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
