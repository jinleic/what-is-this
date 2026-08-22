#!/usr/bin/env python3
"""Clean-room verifier for the open 3x3 two-sum Lie-algebra certificate.

This module deliberately does not import ``experiments/e87_algebra_3x3.py``.
It rebuilds the exact symmetry sectors, commutant-center projectors, scalar
projection, factor linkages, and the two raw D4-orbit Pauli closures.  The
JSON artifact is only a certificate envelope and source of stored witnesses.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_structure" / "char0_3x3.json"
SCRIPT = "experiments/e87_algebra_3x3.py"
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


def orbit_bonds() -> tuple[tuple[int, int], ...]:
    """The independently specified site-major bond traversal used by the raw tree."""
    answer = []
    for r in range(ROWS):
        for c in range(COLS):
            if r + 1 < ROWS:
                answer.append((r * COLS + c, (r + 1) * COLS + c))
            if c + 1 < COLS:
                answer.append((r * COLS + c, r * COLS + c + 1))
    return tuple(answer)


ORBIT_BONDS = orbit_bonds()
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
FACTOR_SPECS = (
    ("A32", "000", 33, 1088, 32),
    ("A47", "000", 48, 2303, 47),
    ("A58", "010", 59, 3480, 58),
    ("A2", "110", 3, 8, 2),
    ("A15", "110", 16, 255, 15),
    ("A29", "110", 30, 899, 29),
)


def check(name: str, condition: bool, detail: str = "") -> None:
    print(("PASS" if condition else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not condition:
        raise AssertionError(name)


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


def integer_certificate(matrix: sp.MatrixBase, store_entries: bool = True) -> dict:
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


def decode_integer_certificate(certificate: dict) -> sp.Matrix:
    if certificate["entries"] is None:
        raise AssertionError("clean-room verification requires matrix entries")
    matrix = sp.zeros(*certificate["shape"])
    for i, j, value in certificate["entries"]:
        matrix[i, j] = value
    assert integer_certificate(matrix) == certificate
    return matrix


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Permutation obtained by applying right, then left, to a site index."""
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


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    output = 0
    for old, new in enumerate(permutation):
        if (mask >> old) & 1:
            output |= 1 << new
    return output


def site_symmetries() -> tuple[tuple[int, int, int, tuple[int, ...]], ...]:
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


def sector(character: tuple[int, int, int]) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Exact C2^3 character block, physical Gram form, and inclusion."""
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
    vectors: list[dict[int, int]] = []
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
        vectors.append({configuration: coefficient // scale for configuration, coefficient in values.items() if coefficient})

    dimension = len(vectors)
    a = sp.zeros(dimension)
    b = sp.zeros(dimension)
    gram = sp.zeros(dimension)
    inclusion = sp.zeros(size, dimension)
    for j, orbit in enumerate(vectors):
        gram[j, j] = sum(coefficient * coefficient for coefficient in orbit.values())
        for configuration, coefficient in orbit.items():
            inclusion[configuration, j] = coefficient
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
    assert inclusion.T * inclusion == gram
    return a, b, gram, inclusion


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
    """Exact T with target_i T = signs_i T source_i, using diagonal B first."""
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


def commutant_dimension(actions: tuple[sp.Matrix, sp.Matrix]) -> int:
    return len(intertwiner_basis(actions, actions))


def invariant_form_nullity(actions: tuple[sp.Matrix, sp.Matrix], alternating: bool) -> tuple[int, int]:
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


def build_components(sectors: dict[tuple[int, int, int], tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]]) -> tuple[dict[str, list[dict]], dict[str, tuple[int, int, tuple[int, ...]]]]:
    components: dict[str, list[dict]] = {}
    metadata: dict[str, tuple[int, int, tuple[int, ...]]] = {}
    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, _gram, _inclusion = sectors[key]
        commutant, center, weights, projectors = central_projectors((a, b))
        metadata[label] = (len(commutant), len(center), weights)
        parts = []
        for rank, projector, eigenvalue in projectors:
            actions = restrict_to_projector((a, b), projector)
            trace_a, trace_b = actions[0].trace(), actions[1].trace()
            traceless = traceless_actions(actions)
            parts.append(
                {
                    "dimension": rank,
                    "projector": projector,
                    "eigenvalue": eigenvalue,
                    "actions": actions,
                    "traceless": traceless,
                    "trace_a": trace_a,
                    "trace_b": trace_b,
                    "scalar_a": trace_a / rank,
                    "traceless_zero": traceless[0] == sp.zeros(rank) and traceless[1] == sp.zeros(rank),
                }
            )
        components[label] = parts
    return components, metadata


def component_by_dimension(components: dict[str, list[dict]], label: str, dimension: int) -> dict:
    matches = [part for part in components[label] if part["dimension"] == dimension]
    if len(matches) != 1:
        raise AssertionError((label, dimension, len(matches)))
    return matches[0]


def verify_intertwiner_record(record: dict, target: tuple[sp.Matrix, sp.Matrix], source: tuple[sp.Matrix, sp.Matrix]) -> bool:
    signs = tuple(record["signs"])
    basis = intertwiner_basis(target, source, signs)
    witness = full_rank_witness(basis)
    stored = decode_integer_certificate(record["witness"])
    return (
        len(basis) == record["space_dimension_Q"]
        and witness is not None
        and witness.rank() == record["rank_Q"]
        and integer_certificate(witness) == record["witness"]
        and stored.rank() == record["rank_Q"]
        and all(left * stored == sign * stored * right for left, right, sign in zip(target, source, signs))
    )


def pauli_permute(value: int, permutation: tuple[int, ...]) -> int:
    mask = (1 << N) - 1
    a, b = value & mask, (value >> N) & mask
    return permute_mask(a, permutation) | (permute_mask(b, permutation) << N)


def commutator_sign_half(generator: int, value: int) -> int:
    mask = (1 << N) - 1
    first = ((((generator >> N) & mask) & (value & mask)).bit_count()) & 1
    second = ((((value >> N) & mask) & (generator & mask)).bit_count()) & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1


def pauli_grade(value: int) -> int:
    mask = (1 << N) - 1
    a, b = value & mask, (value >> N) & mask
    return ((a.bit_count() & 1) << 1) | ((a & b).bit_count() & 1)


class D4OrbitClosure:
    """Independent raw D4-orbit closure, matching only the mathematical rules."""

    def __init__(self, prime: int) -> None:
        self.prime = prime
        self.group = d4_permutations()
        self.x_terms = tuple(1 << i for i in range(N))
        self.zz_terms = tuple(((1 << u) | (1 << v)) << N for u, v in ORBIT_BONDS)
        self.single_terms = self.x_terms + self.zz_terms
        self.cache: dict[int, int] = {}
        self.support = self.enumerate_support()
        self.blocks = tuple(tuple(value for value in self.support if pauli_grade(value) == block) for block in range(4))
        self.block_sizes = tuple(len(block) for block in self.blocks)
        self.index = tuple({value: i for i, value in enumerate(block)} for block in self.blocks)
        self.actions = tuple(
            tuple(self.action_table(terms, source_block) for source_block in range(4))
            for terms in (self.x_terms, self.zz_terms)
        )

    def canonical(self, value: int) -> int:
        cached = self.cache.get(value)
        if cached is None:
            cached = min(pauli_permute(value, permutation) for permutation in self.group)
            self.cache[value] = cached
        return cached

    def enumerate_support(self) -> tuple[int, ...]:
        seen = {self.canonical(term) for term in self.single_terms}
        queue = sorted(seen)
        head = 0
        while head < len(queue):
            value = queue[head]
            head += 1
            for generator in self.single_terms:
                if commutator_sign_half(generator, value):
                    target = self.canonical(value ^ generator)
                    if target not in seen:
                        seen.add(target)
                        queue.append(target)
        mask = (1 << N) - 1
        assert not any(((value >> N) & mask).bit_count() & 1 for value in seen)
        return tuple(reversed(queue))

    def action_table(self, terms: tuple[int, ...], source_block: int) -> tuple[np.ndarray, np.ndarray]:
        target_block = source_block ^ (3 if len(terms) == N else 1)
        sources = np.full((len(terms), self.block_sizes[target_block]), -1, dtype=np.int32)
        signs = np.zeros((len(terms), self.block_sizes[target_block]), dtype=np.int8)
        for term_index, generator in enumerate(terms):
            for target_index, target in enumerate(self.blocks[target_block]):
                source = target ^ generator
                sign = commutator_sign_half(generator, source)
                if sign:
                    source_index = self.index[source_block].get(self.canonical(source))
                    if source_index is None:
                        raise AssertionError("support not closed")
                    sources[term_index, target_index] = source_index
                    signs[term_index, target_index] = sign
        return sources, signs

    def seed(self, terms: tuple[int, ...]) -> tuple[int, np.ndarray]:
        term_multiplicities: dict[int, int] = {}
        for term in terms:
            term_multiplicities[term] = term_multiplicities.get(term, 0) + 1
        coefficients: dict[int, int] = {}
        for term, multiplicity in term_multiplicities.items():
            orbit = self.canonical(term)
            previous = coefficients.setdefault(orbit, multiplicity)
            if previous != multiplicity:
                raise AssertionError("nonuniform coefficient within point-group orbit")
        grades = {pauli_grade(orbit) for orbit in coefficients}
        assert len(grades) == 1
        block = grades.pop()
        row = np.zeros(self.block_sizes[block], dtype=np.int64)
        for orbit, coefficient in coefficients.items():
            row[self.index[block][orbit]] = coefficient
        return block, row

    def closure(self, derived: bool) -> dict:
        rows = tuple(np.empty((size, size), dtype=np.int32) for size in self.block_sizes)
        pivot_row = tuple(np.full(size, -1, dtype=np.int32) for size in self.block_sizes)
        pivot_inverse = tuple(np.zeros(size, dtype=np.int32) for size in self.block_sizes)
        row_counts = [0, 0, 0, 0]
        words: list[tuple[int, np.ndarray]] = []
        recipes: list[list[int]] = []

        def reduce_add(vector: tuple[int, np.ndarray]) -> tuple[int, int] | None:
            block, raw = vector
            value = np.asarray(raw, dtype=np.int64) % self.prime
            while True:
                nonzero = np.flatnonzero(value)
                if not nonzero.size:
                    return None
                lead = int(nonzero[0])
                existing = int(pivot_row[block][lead])
                if existing < 0:
                    row_index = row_counts[block]
                    insertion = int(value[lead])
                    rows[block][row_index] = value.astype(np.int32)
                    pivot_row[block][lead] = row_index
                    pivot_inverse[block][lead] = pow(insertion, self.prime - 2, self.prime)
                    row_counts[block] += 1
                    return block, row_index
                factor = int(value[lead]) * int(pivot_inverse[block][lead]) % self.prime
                basis = rows[block][existing, lead:].astype(np.int64)
                product = basis * factor
                product = (product & self.prime) + (product >> 31)
                product = (product & self.prime) + (product >> 31)
                product[product >= self.prime] -= self.prime
                tail = value[lead:]
                tail -= product
                tail[tail < 0] += self.prime

        def apply_raw(action_index: int, word: tuple[int, np.ndarray]) -> tuple[int, np.ndarray]:
            source_block, row = word
            target_block = source_block ^ (3 if action_index == 0 else 1)
            sources, signs = self.actions[action_index][source_block]
            output = np.zeros(self.block_sizes[target_block], dtype=np.int64)
            for source, sign in zip(sources, signs):
                selected = source >= 0
                output[selected] += sign[selected].astype(np.int64) * row[source[selected]].astype(np.int64)
            return target_block, output % self.prime

        frontier: list[int] = []
        if derived:
            seed = apply_raw(0, self.seed(self.zz_terms))
            seeds = [(seed, [-1, 0, 2])]
        else:
            seeds = [
                (self.seed(self.x_terms), [-1, 0, 1]),
                (self.seed(self.zz_terms), [-1, 1, 1]),
            ]
        for word, recipe in seeds:
            assert reduce_add(word) is not None
            words.append((word[0], np.asarray(word[1], dtype=np.int64) % self.prime))
            recipes.append(recipe)
            frontier.append(len(words) - 1)
        while frontier:
            next_frontier: list[int] = []
            for parent in frontier:
                for action_index in (0, 1):
                    candidate = apply_raw(action_index, words[parent])
                    if reduce_add(candidate) is not None:
                        words.append(candidate)
                        recipes.append([parent, action_index, recipes[parent][2] + 1])
                        next_frontier.append(len(words) - 1)
            frontier = next_frontier
        return {
            "rank_Fp": len(words),
            "maximum_depth": max(recipe[2] for recipe in recipes),
            "grading_block_ranks": row_counts,
            "tree_checksum_mod_1000000007": sum(
                (index + 1) * (recipe[0] + 2) * (recipe[1] + 2) * (recipe[2] + 1)
                for index, recipe in enumerate(recipes)
            )
            % 1_000_000_007,
            "tree": recipes,
        }


def main() -> int:
    started = time.monotonic()
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(artifact) == {"provenance", "data", "checks"}
        and artifact["provenance"]["script"] == SCRIPT
        and all(set(item) == {"name", "passed", "detail"} for item in artifact["checks"]),
    )
    data = artifact["data"]
    exact = data["exact_Q_certificate"]

    d4 = d4_permutations()
    bond_set = {tuple(sorted(bond)) for bond in BONDS}
    check(
        "D4 exact symmetry group",
        len(d4) == 8
        and all({tuple(sorted((permutation[u], permutation[v]))) for u, v in BONDS} == bond_set for permutation in d4)
        and data["symmetry"]["d4_order"] == 8
        and data["symmetry"]["d4_permutations"] == [list(permutation) for permutation in d4],
    )

    sectors = {key: sector(key) for key in CHARACTERS}
    dimensions = {"".join(map(str, key)): sectors[key][0].rows for key in CHARACTERS}
    inclusions = [sectors[key][3] for key in CHARACTERS]
    complete = sp.Matrix.hstack(*inclusions)
    complete_gram = complete.T * complete
    check(
        "exact C2^3 sector completeness",
        dimensions == EXPECTED_SECTOR_DIMS
        and sum(dimensions.values()) == 512
        and all(complete_gram[i, j] == 0 for i in range(512) for j in range(512) if i != j)
        and all(complete_gram[i, i] != 0 for i in range(512)),
    )

    components, metadata = build_components(sectors)
    for label in EXPECTED_SECTOR_DIMS:
        a, b, gram, _inclusion = sectors[tuple(map(int, label))]
        record = exact["sector_certificates_Q"][label]
        commutant_dim, center_dim, weights = metadata[label]
        observed_parts = components[label]
        check(
            f"{label} exact isotypic sector decomposition",
            a.rows == EXPECTED_SECTOR_DIMS[label]
            and [part["dimension"] for part in observed_parts] == EXPECTED_BLOCK_DIMS[label]
            and (commutant_dim, center_dim) == EXPECTED_COMMUTANTS[label]
            and record["module_dimension_Q"] == a.rows
            and record["physical_gram"] == integer_certificate(gram)
            and record["commutant_dimension_Q"] == commutant_dim
            and record["commutant_center_dimension_Q"] == center_dim
            and record["central_separator_weights"] == list(weights)
            and [part["dimension_Q"] for part in record["isotypic_blocks"]] == EXPECTED_BLOCK_DIMS[label]
            and all(
                block["projector"] == integer_certificate(part["projector"])
                and block["trace_A_Q"] == str(part["trace_a"])
                and block["trace_B_Q"] == str(part["trace_b"])
                and block["scalar_A_Q"] == str(part["scalar_a"])
                and block["traceless_actions_zero"] == part["traceless_zero"]
                for block, part in zip(record["isotypic_blocks"], observed_parts)
            ),
        )

    factor_controls = exact["factor_controls_Q"]
    for factor, label, dimension, factor_dimension, rank in FACTOR_SPECS:
        part = component_by_dimension(components, label, dimension)
        commutant = commutant_dimension(part["traceless"])
        symmetric_variables, symmetric_nullity = invariant_form_nullity(part["traceless"], False)
        alternating_variables, alternating_nullity = invariant_form_nullity(part["traceless"], True)
        record = factor_controls[factor]
        check(
            f"{factor} exact split-A control",
            commutant == 1
            and symmetric_nullity == 0
            and alternating_nullity == 0
            and record == {
                "module_dimension_Q": dimension,
                "factor_dimension_Q": factor_dimension,
                "rank": rank,
                "traceless_commutant_dimension_Q": 1,
                "invariant_symmetric_form_variables_Q": symmetric_variables,
                "invariant_symmetric_form_nullity_Q": 0,
                "invariant_alternating_form_variables_Q": alternating_variables,
                "invariant_alternating_form_nullity_Q": 0,
                "classification": f"sl({dimension},Q), absolute type {factor}",
            },
        )

    scalar_pieces = []
    residual_checks = []
    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, _gram, _inclusion = sectors[key]
        scalar = sp.zeros(a.rows)
        for part in components[label]:
            scalar += part["scalar_a"] * part["projector"]
            residual_a = part["actions"][0] - part["scalar_a"] * sp.eye(part["dimension"])
            residual_checks.append(residual_a.trace() == 0 and part["actions"][1].trace() == 0)
        assert scalar * a == a * scalar and scalar * b == b * scalar
        scalar_pieces.append(scalar)
    z = sp.diag(*scalar_pieces)
    full_a = sp.diag(*(sectors[key][0] for key in CHARACTERS))
    full_b = sp.diag(*(sectors[key][1] for key in CHARACTERS))
    z_integer = primitive_integer_matrix(z)
    z_record = exact["central_witness_global_integer"]
    z_stored = decode_integer_certificate(z_record)
    z_scale = sp.Rational(exact["central_projection_equals_integer_times"])
    check(
        "global A-scalar central witness",
        z != sp.zeros(512)
        and z_integer == z_stored
        and z == z_scale * z_stored
        and z * full_a == full_a * z
        and z * full_b == full_b * z
        and all(residual_checks)
        and exact["central_projection_generator"] == "A",
    )

    intertwiners = exact["intertwiners_Q"]
    for record in intertwiners["P_signed"]:
        source = component_by_dimension(components, record["source_sector"], record["component_dimension_Q"])["traceless"]
        target = component_by_dimension(components, record["target_sector"], record["component_dimension_Q"])["traceless"]
        check(
            f"P signed linkage {record['source_sector']}->{record['target_sector']}:{record['component_dimension_Q']}",
            verify_intertwiner_record(record, target, source),
        )
        outer = record["outer_sign_witness"]
        outer_basis = intertwiner_basis(source, (source[0].T, source[1].T), (1, -1))
        outer_witness = full_rank_witness(outer_basis)
        stored_outer = decode_integer_certificate(outer["witness"])
        check(
            f"outer A-sign automorphism {record['source_sector']}:{record['component_dimension_Q']}",
            len(outer_basis) == outer["space_dimension_Q"]
            and outer_witness is not None
            and integer_certificate(outer_witness) == outer["witness"]
            and stored_outer.rank() == record["component_dimension_Q"]
            and source[0] * stored_outer == stored_outer * source[0].T
            and source[1] * stored_outer == -stored_outer * source[1].T,
        )
        dual = decode_integer_certificate(record["dual_witness"])
        check(
            f"dual linkage {record['source_sector']}->{record['target_sector']}:{record['component_dimension_Q']}",
            dual.rank() == record["component_dimension_Q"]
            and target[0] * dual == -dual * source[0].T
            and target[1] * dual == -dual * source[1].T,
        )
    for record in intertwiners["D4_direct"]:
        source = component_by_dimension(components, record["source_sector"], record["component_dimension_Q"])["traceless"]
        target = component_by_dimension(components, record["target_sector"], record["component_dimension_Q"])["traceless"]
        check(
            f"D4 direct linkage {record['source_sector']}->{record['target_sector']}",
            verify_intertwiner_record(record, target, source),
        )

    container = exact["containing_algebra_Q"]
    check(
        "faithful linked exact-Q container",
        container["dimension_Q"] == 8034
        and container["derived_dimension_Q"] == 8033
        and container["type_with_center"] == "centre(1) + A32 + A47 + A58 + A2 + A15 + A29"
        and container["factor_dimensions"] == [1088, 2303, 3480, 8, 255, 899]
        and container["factor_ranks"] == [32, 47, 58, 2, 15, 29]
        and container["levi_rank"] == 183
        and container["per_sector_image_dimensions_Q"] == {
            "000": 3392,
            "001": 3392,
            "010": 3481,
            "011": 3481,
            "100": 3481,
            "101": 3481,
            "110": 1163,
            "111": 1163,
        },
    )

    pauli = D4OrbitClosure(GOOD_PRIME)
    full = pauli.closure(derived=False)
    derived = pauli.closure(derived=True)
    modular = data["modular_word_certificate"]
    for label, observed, expected_rank in (
        ("full_closure", full, 8034),
        ("derived_ideal", derived, 8033),
    ):
        stored = modular[label]
        check(
            f"independent raw D4 Pauli closure {label}",
            observed["rank_Fp"] == stored["rank_Fp"] == expected_rank
            and observed["tree"] == stored["tree"]
            and observed["grading_block_ranks"] == stored["grading_block_ranks"]
            and observed["maximum_depth"] == stored["maximum_depth"]
            and observed["tree_checksum_mod_1000000007"] == stored["tree_checksum_mod_1000000007"],
        )
    check(
        "raw closure support and lower-bound scope",
        modular["prime"] == GOOD_PRIME
        and modular["support_orbits"] == 8739
        and full["grading_block_ranks"] == [1962, 1963, 2147, 1962]
        and derived["grading_block_ranks"] == [1962, 1963, 2146, 1962],
    )

    check(
        "dimension and derived squeeze",
        data["dimension_Q"] == {
            "claim_tag": "[THEOREM]",
            "value": 8034,
            "upper_bound_Q": 8034,
            "lower_bound_from_Fp_word_minor": 8034,
        }
        and data["derived_algebra"]["dimension_Q"] == 8033
        and data["derived_algebra"]["upper_bound_Q"] == 8033
        and data["derived_algebra"]["lower_bound_from_Fp_word_minor"] == 8033,
    )
    check(
        "central quotient radical and structural Killing form",
        data["factor_map"]["kernel_dimension_Q"] == 1
        and data["factor_map"]["kernel_equals_center"] is True
        and data["solvable_radical"] == {
            "claim_tag": "[THEOREM]",
            "dimension_Q": 1,
            "equals_center": True,
            "proof": "the quotient by the exact central kernel is the displayed split semisimple direct sum; every solvable ideal maps trivially to it",
        }
        and data["Killing"]["rank_Q"] == 8033
        and data["Killing"]["nullity_Q"] == 1
        and data["Levi_type"]["type"] == "A32 + A47 + A58 + A2 + A15 + A29"
        and data["Levi_type"]["rank"] == 183,
    )
    check("all producer checks pass", all(item["passed"] for item in artifact["checks"]))
    print(f"WALL_SECONDS={time.monotonic() - started:.3f}")
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
