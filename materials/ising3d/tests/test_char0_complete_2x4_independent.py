#!/usr/bin/env python3
"""Clean-room verifier for the finite open-2x4 characteristic-zero theorem.

This file deliberately does not import any experiment module.  It reconstructs
sector matrices, exact-Q decompositions, literal matrix-word minors, and Pauli
closures from stable library primitives and raw formulas.  The result artifact
is input only for certificate recipes, pivots, witnesses, and claimed fields.
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

from ising.clifford.fast_lie import support_closure

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "algebra_structure" / "char0_complete_2x4.json"
ROWS, COLS, N = 2, 4, 8
GOOD_PRIME = 2_147_483_647
MINOR_PRIME = 65_521
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
X_TERMS = tuple(1 << i for i in range(N))
ZZ_TERMS = tuple(((1 << u) | (1 << v)) << N for u, v in BONDS)
SINGLE_TERMS = X_TERMS + ZZ_TERMS
GROUPS = (tuple(range(len(X_TERMS))), tuple(range(len(X_TERMS), len(SINGLE_TERMS))))


def check(name, condition, detail=""):
    print(("PASS" if condition else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not condition:
        raise AssertionError(name)


def primitive_integer_matrix(matrix):
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


def integer_certificate(matrix, store_entries=True):
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
        ) % 1_000_000_007,
    }


def decode_integer_certificate(certificate):
    if certificate["entries"] is None:
        raise AssertionError("clean-room verifier requires stored global entries")
    matrix = sp.zeros(*certificate["shape"])
    for i, j, value in certificate["entries"]:
        matrix[i, j] = value
    assert integer_certificate(matrix, store_entries=True) == certificate
    return matrix


def site_symmetries():
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


def sector(character):
    """Exact sector generators, Gram form, and physical inclusion columns."""
    size = 1 << N
    symmetries = site_symmetries()

    def act(configuration, symmetry):
        output = 0
        for old, new in enumerate(symmetry[3]):
            if (configuration >> old) & 1:
                output |= 1 << new
        return output ^ ((size - 1) if symmetry[2] else 0)

    seen = set()
    representatives = []
    vectors = []
    for configuration in range(size):
        if configuration in seen:
            continue
        seen.update(act(configuration, symmetry) for symmetry in symmetries)
        values = {}
        for symmetry in symmetries:
            sign = -1 if sum(
                x * y for x, y in zip(character, symmetry[:3])
            ) & 1 else 1
            image = act(configuration, symmetry)
            values[image] = values.get(image, 0) + sign
        if not any(values.values()):
            continue
        representative = min(q for q, coefficient in values.items() if coefficient)
        scale = values[representative]
        representatives.append(representative)
        vectors.append(
            {q: coefficient // scale for q, coefficient in values.items() if coefficient}
        )

    dimension = len(vectors)
    a = sp.zeros(dimension)
    b = sp.zeros(dimension)
    inclusion = sp.zeros(size, dimension)
    gram = sp.zeros(dimension)
    for j, orbit in enumerate(vectors):
        for configuration, coefficient in orbit.items():
            inclusion[configuration, j] = coefficient
        gram[j, j] = sum(coefficient * coefficient for coefficient in orbit.values())
        image = {}
        for configuration, coefficient in orbit.items():
            for site_index in range(N):
                target = configuration ^ (1 << site_index)
                image[target] = image.get(target, 0) + coefficient
        for i, representative in enumerate(representatives):
            a[i, j] = image.get(representative, 0)
        b[j, j] = sum(
            1
            if ((representatives[j] >> u) & 1)
            == ((representatives[j] >> v) & 1)
            else -1
            for u, v in BONDS
        )
    assert inclusion.T * inclusion == gram
    assert inclusion * a == physical_a() * inclusion
    assert inclusion * b == physical_b() * inclusion
    return a, b, gram, inclusion


def physical_a():
    answer = sp.zeros(1 << N)
    for configuration in range(1 << N):
        for site_index in range(N):
            answer[configuration ^ (1 << site_index), configuration] += 1
    return answer


def physical_b():
    answer = sp.zeros(1 << N)
    for configuration in range(1 << N):
        answer[configuration, configuration] = sum(
            1
            if ((configuration >> u) & 1) == ((configuration >> v) & 1)
            else -1
            for u, v in BONDS
        )
    return answer


def gram_split(actions, gram):
    dimension = actions[0].rows
    kernel = sp.Matrix.vstack(*actions).nullspace()
    kernel_matrix = sp.Matrix.hstack(*kernel) if kernel else sp.zeros(dimension, 0)
    complement = (
        (kernel_matrix.T * gram).nullspace()
        if kernel
        else [sp.eye(dimension)[:, j] for j in range(dimension)]
    )
    complement_matrix = sp.Matrix.hstack(*complement)
    change = sp.Matrix.hstack(kernel_matrix, complement_matrix)
    inverse = change.inv()
    conjugated = tuple(inverse * generator * change for generator in actions)
    k = kernel_matrix.cols
    assert all(
        generator[:k, :] == sp.zeros(k, dimension)
        and generator[:, :k] == sp.zeros(dimension, k)
        for generator in conjugated
    )
    return kernel_matrix, complement_matrix, tuple(
        generator[k:, k:] for generator in conjugated
    )


def sparse_nullspace(equations, rows, columns):
    matrix = DomainMatrix(equations, (rows, columns), QQ)
    return matrix.nullspace().to_Matrix().tolist()


def intertwiner_basis(target, source, signs=(1, 1)):
    target_a, target_b = target
    source_a, source_b = source
    sign_a, sign_b = signs
    diagonal = (
        target_b == sp.diag(*(target_b[i, i] for i in range(target_b.rows)))
        and source_b == sp.diag(*(source_b[i, i] for i in range(source_b.rows)))
    )
    pairs = [
        (i, j)
        for i in range(target_a.rows)
        for j in range(source_a.rows)
        if not diagonal or target_b[i, i] == sign_b * source_b[j, j]
    ]
    variable = {pair: column for column, pair in enumerate(pairs)}
    equations = {}
    row_index = 0
    generators = ((target_a, source_a, sign_a),) if diagonal else (
        (target_a, source_a, sign_a),
        (target_b, source_b, sign_b),
    )
    for target_generator, source_generator, sign in generators:
        for i in range(target_a.rows):
            for j in range(source_a.rows):
                row = {}
                for k in range(target_a.rows):
                    column = variable.get((k, j))
                    coefficient = target_generator[i, k]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) + coefficient
                for k in range(source_a.rows):
                    column = variable.get((i, k))
                    coefficient = source_generator[k, j]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) - sign * coefficient
                if row:
                    equations[row_index] = {
                        column: QQ.convert(value)
                        for column, value in row.items()
                        if value
                    }
                    row_index += 1
    answer = []
    for vector in sparse_nullspace(equations, row_index, len(pairs)):
        matrix = sp.zeros(target_a.rows, source_a.rows)
        for coefficient, (i, j) in zip(vector, pairs):
            matrix[i, j] = coefficient
        answer.append(primitive_integer_matrix(matrix))
    return answer


def commutant_basis(actions):
    return intertwiner_basis(actions, actions)


def multiplicity_free_projectors(actions):
    a, b = actions
    commutant = commutant_basis(actions)
    for weights in itertools.product(range(1, 5), repeat=len(commutant)):
        operator = sum(
            (weight * matrix for weight, matrix in zip(weights, commutant)),
            sp.zeros(a.rows),
        )
        eigenvalues = operator.eigenvals()
        if len(eigenvalues) != len(commutant) or not all(
            value.is_Rational for value in eigenvalues
        ):
            continue
        identity = sp.eye(a.rows)
        projectors = []
        for value, multiplicity in eigenvalues.items():
            projector = identity
            for other in eigenvalues:
                if other != value:
                    projector = projector * (operator - other * identity) / (value - other)
            assert projector * projector == projector
            assert projector * a == a * projector and projector * b == b * projector
            projectors.append((multiplicity, projector))
        return tuple(sorted(projectors, key=lambda item: (item[0], list(item[1]))))
    raise AssertionError("no rational separating commutant element")


def restrict_to_projector(actions, projector):
    inclusion = sp.Matrix.hstack(*projector.columnspace())
    left_inverse = (inclusion.T * inclusion).inv() * inclusion.T
    assert inclusion * left_inverse * projector == projector
    return tuple(left_inverse * action * inclusion for action in actions)


def invariant_symmetric_forms(actions):
    dimension = actions[0].rows
    b = actions[1]
    diagonal = b == sp.diag(*(b[i, i] for i in range(dimension)))
    pairs = [
        (i, j)
        for i in range(dimension)
        for j in range(i, dimension)
        if not diagonal or b[i, i] + b[j, j] == 0
    ]
    variable = {pair: column for column, pair in enumerate(pairs)}
    equations = {}
    row_index = 0
    generators = actions[:1] if diagonal else actions
    for generator in generators:
        for i in range(dimension):
            for j in range(dimension):
                row = {}
                for k in range(dimension):
                    column = variable.get(tuple(sorted((k, j))))
                    coefficient = generator[k, i]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) + coefficient
                    column = variable.get(tuple(sorted((i, k))))
                    coefficient = generator[k, j]
                    if column is not None and coefficient:
                        row[column] = row.get(column, 0) + coefficient
                if row:
                    equations[row_index] = {
                        column: QQ.convert(value)
                        for column, value in row.items()
                        if value
                    }
                    row_index += 1
    answer = []
    for vector in sparse_nullspace(equations, row_index, len(pairs)):
        form = sp.zeros(dimension)
        for coefficient, (i, j) in zip(vector, pairs):
            form[i, j] = coefficient
            form[j, i] = coefficient
        form = primitive_integer_matrix(form)
        assert all(
            generator.T * form + form * generator == sp.zeros(dimension)
            for generator in actions
        )
        answer.append(form)
    return answer


def traceless_actions(actions):
    dimension = actions[0].rows
    identity = sp.eye(dimension)
    return tuple(
        action - action.trace() * identity / dimension for action in actions
    )


def full_rank_witness(matrices):
    if not matrices:
        return None
    target = min(matrices[0].shape)
    for weights in itertools.product(range(1, 5), repeat=len(matrices)):
        candidate = sum(
            (weight * matrix for weight, matrix in zip(weights, matrices)),
            sp.zeros(*matrices[0].shape),
        )
        if candidate.rank() == target:
            return primitive_integer_matrix(candidate)
    return None


def rational_matrix_mod(matrix, prime):
    answer = np.zeros(matrix.shape, dtype=np.int64)
    denominators = []
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            value = sp.Rational(matrix[i, j])
            if value.q % prime == 0:
                raise AssertionError("bad reduction prime")
            denominators.append(int(value.q))
            answer[i, j] = (
                (int(value.p) % prime) * pow(int(value.q) % prime, prime - 2, prime)
            ) % prime
    return answer, math.lcm(*denominators)


def raw_word_minor_clean(blocks, stored):
    """Independently generate and eliminate every literal matrix word."""
    prime = stored["prime"]
    converted = []
    denominator_lcms = []
    for a, b in blocks:
        a_mod, a_lcm = rational_matrix_mod(a, prime)
        b_mod, b_lcm = rational_matrix_mod(b, prime)
        converted.append((a_mod, b_mod))
        denominator_lcms.append([a_lcm, b_lcm])
    generators = [
        tuple(block[index] for block in converted) for index in (0, 1)
    ]
    raw_words = []
    recipes = []
    reduced = []
    pivots = {}
    pivot_columns = []
    insertion_pivots = []
    elimination_traces = []
    frontier = []

    def flatten(word):
        return np.concatenate([matrix.ravel() for matrix in word])

    def accept(raw_word, recipe):
        value = flatten(raw_word).copy() % prime
        trace = []
        while True:
            nonzero = np.flatnonzero(value)
            if not nonzero.size:
                return None
            lead = int(nonzero[0])
            previous = pivots.get(lead)
            if previous is None:
                insertion = int(value[lead])
                normalized = value * pow(insertion, prime - 2, prime) % prime
                pivots[lead] = len(reduced)
                reduced.append(normalized)
                raw_words.append(tuple(matrix.copy() for matrix in raw_word))
                recipes.append(recipe)
                pivot_columns.append(lead)
                insertion_pivots.append(insertion)
                elimination_traces.append(trace)
                return len(raw_words) - 1
            coefficient = int(value[lead])
            trace.append([previous, coefficient])
            value = (value - coefficient * reduced[previous]) % prime

    for generator_index, generator in enumerate(generators):
        accepted = accept(generator, [-1, generator_index, 1])
        if accepted is not None:
            frontier.append(accepted)
    head = 0
    while head < len(frontier):
        parent = frontier[head]
        head += 1
        for generator_index, generator in enumerate(generators):
            raw_parent = raw_words[parent]
            candidate = tuple(
                (left @ right - right @ left) % prime
                for left, right in zip(generator, raw_parent)
            )
            accepted = accept(
                candidate,
                [parent, generator_index, recipes[parent][2] + 1],
            )
            if accepted is not None:
                frontier.append(accepted)

    # Explicitly select the stored/computed pivot columns in the actual raw
    # word matrix.  Its determinant is recomputed by an independent dense
    # elimination, rather than inferred from the producer's pivot product.
    selected_columns = np.asarray(pivot_columns, dtype=np.int64)
    selected_minor = np.asarray(
        [flatten(raw_word)[selected_columns] % prime for raw_word in raw_words],
        dtype=np.int64,
    )
    independent_minor = selected_minor.copy()
    independent_determinant = 1
    for column in range(len(raw_words)):
        choices = np.flatnonzero(independent_minor[column:, column])
        assert choices.size
        pivot_row = column + int(choices[0])
        if pivot_row != column:
            independent_minor[[column, pivot_row]] = independent_minor[[pivot_row, column]]
            independent_determinant = -independent_determinant % prime
        pivot = int(independent_minor[column, column])
        independent_determinant = independent_determinant * pivot % prime
        inverse = pow(pivot, prime - 2, prime)
        independent_minor[column] = independent_minor[column] * inverse % prime
        for row_start in range(column + 1, len(raw_words), 64):
            row_indices = np.arange(
                row_start, min(row_start + 64, len(raw_words)), dtype=np.int64
            )
            independent_minor[row_indices] = (
                independent_minor[row_indices]
                - independent_minor[row_indices, column, None]
                * independent_minor[column]
            ) % prime

    determinant = math.prod(insertion_pivots) % prime
    assert independent_determinant == determinant != 0
    observed = {
        "rank_Fp": len(raw_words),
        "recipes": recipes,
        "pivot_columns": pivot_columns,
        "insertion_pivots": insertion_pivots,
        "determinant": determinant,
        "denominator_lcms": denominator_lcms,
        "maximum_depth": max(recipe[2] for recipe in recipes),
        "checksum": sum(
            (index + 1)
            * (recipe[0] + 2)
            * (recipe[1] + 2)
            * (recipe[2] + 1)
            for index, recipe in enumerate(recipes)
        ) % 1_000_000_007,
    }
    assert observed["rank_Fp"] == stored["rank_Fp"]
    assert observed["recipes"] == stored["raw_parent_recipes"]
    assert observed["pivot_columns"] == stored["pivot_columns"]
    assert observed["insertion_pivots"] == stored["unnormalised_insertion_pivots"]
    assert observed["determinant"] == stored["selected_minor_determinant_mod_p"] != 0
    assert observed["denominator_lcms"] == stored["generator_denominator_lcms"]
    assert observed["maximum_depth"] == stored["maximum_word_depth"]
    assert observed["checksum"] == stored["recipe_checksum_mod_1000000007"]
    return observed


def pauli_grade(value):
    mask = (1 << N) - 1
    a = value & mask
    b = (value >> N) & mask
    return ((a.bit_count() & 1) << 1) | ((a & b).bit_count() & 1)


def commutator_sign_half(generator, value):
    mask = (1 << N) - 1
    first = ((((generator >> N) & mask) & (value & mask)).bit_count()) & 1
    second = ((((value >> N) & mask) & (generator & mask)).bit_count()) & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1

def point_group_orbit_count(values):
    permutations = (
        tuple(range(N)),
        tuple((ROWS - 1 - r) * COLS + c for r in range(ROWS) for c in range(COLS)),
        tuple(r * COLS + (COLS - 1 - c) for r in range(ROWS) for c in range(COLS)),
        tuple(
            (ROWS - 1 - r) * COLS + (COLS - 1 - c)
            for r in range(ROWS)
            for c in range(COLS)
        ),
    )
    mask = (1 << N) - 1

    def permute(value, permutation):
        a, b = value & mask, (value >> N) & mask
        new_a = new_b = 0
        for old, new in enumerate(permutation):
            if (a >> old) & 1:
                new_a |= 1 << new
            if (b >> old) & 1:
                new_b |= 1 << new
        return new_a | (new_b << N)

    return len(
        {
            min(permute(value, permutation) for permutation in permutations)
            for value in values
        }
    )


class PauliRecipeClosure:
    """Independent four-grading raw-word closure modulo the good prime."""

    def __init__(self, prime):
        self.prime = prime
        support = support_closure(SINGLE_TERMS, SINGLE_TERMS, N, cap=100_000)
        self.support_orbits = point_group_orbit_count(support)
        self.blocks = tuple(
            tuple(value for value in support if pauli_grade(value) == grade)
            for grade in range(4)
        )
        self.index = tuple(
            {value: index for index, value in enumerate(block)}
            for block in self.blocks
        )
        self.transitions = {}
        for action_index, group in enumerate(GROUPS):
            for source_block, source_values in enumerate(self.blocks):
                target_block = source_block ^ (3 if action_index == 0 else 1)
                pieces = []
                for term_index in group:
                    generator = SINGLE_TERMS[term_index]
                    sources, targets, signs = [], [], []
                    for source_index, value in enumerate(source_values):
                        sign = commutator_sign_half(generator, value)
                        if sign:
                            sources.append(source_index)
                            targets.append(self.index[target_block][generator ^ value])
                            signs.append(sign)
                    pieces.append(
                        (
                            np.asarray(sources, dtype=np.int64),
                            np.asarray(targets, dtype=np.int64),
                            np.asarray(signs, dtype=np.int64),
                        )
                    )
                self.transitions[action_index, source_block] = (target_block, pieces)

    def seed(self, terms):
        block = pauli_grade(terms[0])
        assert all(pauli_grade(term) == block for term in terms)
        row = np.zeros(len(self.blocks[block]), dtype=np.int64)
        for term in terms:
            row[self.index[block][term]] += 1
        return block, row % self.prime

    def apply(self, action_index, word):
        source_block, row = word
        target_block, pieces = self.transitions[action_index, source_block]
        result = np.zeros(len(self.blocks[target_block]), dtype=np.int64)
        for sources, targets, signs in pieces:
            result[targets] += signs * row[sources]
        return target_block, result % self.prime

    def canonical(self, derived=False):
        rows = [[] for _ in range(4)]
        pivots = [{} for _ in range(4)]
        words = []
        recipes = []
        frontier = []

        def accept(word, recipe):
            block, raw = word
            value = raw.copy() % self.prime
            while True:
                nonzero = np.flatnonzero(value)
                if not nonzero.size:
                    return None
                lead = int(nonzero[0])
                previous = pivots[block].get(lead)
                if previous is None:
                    value = value * pow(int(value[lead]), self.prime - 2, self.prime) % self.prime
                    pivots[block][lead] = len(rows[block])
                    rows[block].append(value)
                    words.append((block, raw.copy() % self.prime))
                    recipes.append(recipe)
                    return len(words) - 1
                value = (
                    value - int(value[lead]) * rows[block][previous]
                ) % self.prime

        if derived:
            seed = self.apply(0, self.seed(ZZ_TERMS))
            accepted = accept(seed, [-1, 0, 2])
            assert accepted == 0
            frontier = [accepted]
        else:
            for generator_index, terms in enumerate((X_TERMS, ZZ_TERMS)):
                accepted = accept(self.seed(terms), [-1, generator_index, 1])
                assert accepted is not None
                frontier.append(accepted)
        while frontier:
            next_frontier = []
            for parent in frontier:
                for action_index in (0, 1):
                    candidate = self.apply(action_index, words[parent])
                    accepted = accept(
                        candidate,
                        [parent, action_index, recipes[parent][2] + 1],
                    )
                    if accepted is not None:
                        next_frontier.append(accepted)
            frontier = next_frontier
        return {
            "rank": len(words),
            "tree": recipes,
            "grading_block_ranks": [len(block_rows) for block_rows in rows],
            "support_orbits": sum(len(block) for block in self.blocks),
            "maximum_depth": max(recipe[2] for recipe in recipes),
            "checksum": sum(
                (index + 1)
                * (recipe[0] + 2)
                * (recipe[1] + 2)
                * (recipe[2] + 1)
                for index, recipe in enumerate(recipes)
            ) % 1_000_000_007,
            "support_orbits": self.support_orbits,
        }


def verify_intertwiner_record(record, target, source):
    source = traceless_actions(source)
    target = traceless_actions(target)
    if record["dual"]:
        source = tuple(-matrix.T for matrix in source)
    basis = intertwiner_basis(target, source)
    witness = full_rank_witness(basis)
    stored_witness = decode_integer_certificate(record["witness"])
    return (
        len(basis) == record["space_dimension_Q"]
        and witness is not None
        and witness.rank() == record["rank_Q"]
        and stored_witness.rank() == record["rank_Q"]
        and all(
            target_generator * stored_witness
            == stored_witness * source_generator
            for target_generator, source_generator in zip(target, source)
        )
    )


def main():
    started = time.monotonic()
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(artifact) == {"provenance", "data", "checks"}
        and artifact["provenance"]["script"] == "experiments/e76_char0_complete_2x4.py"
        and all(set(item) == {"name", "passed", "detail"} for item in artifact["checks"]),
    )
    data = artifact["data"]
    exact = data["exact_Q_certificate"]

    sectors = {key: sector(key) for key in CHARACTERS}
    dimensions = [sectors[key][0].rows for key in CHARACTERS]
    inclusions = [sectors[key][3] for key in CHARACTERS]
    complete_change = sp.Matrix.hstack(*inclusions)
    complete_gram = complete_change.T * complete_change
    check(
        "exact sector completeness",
        dimensions == [44, 32, 28, 32, 28, 32, 28, 32]
        and sum(dimensions) == 256
        and complete_gram.det() != 0,
    )

    active = {}
    components = {}
    expected_kernels = {"000": 2, "010": 1, "100": 1, "110": 1}
    expected_forms = {"000": (1, 42), "010": (1, 27), "100": (1, 27), "110": (1, 27)}
    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, gram, _ = sectors[key]
        kernel, complement, active_actions = gram_split((a, b), gram)
        active[label] = active_actions
        record = exact["sector_certificates_Q"][label]
        check(
            f"{label} exact common-kernel split",
            kernel.cols == record["common_kernel_dimension_Q"]
            and integer_certificate(kernel) == record["kernel_basis"]
            and integer_certificate(complement) == record["gram_orthogonal_complement"],
        )
        if label in expected_kernels:
            forms = invariant_symmetric_forms(active_actions)
            check(
                f"{label} rational orthogonal form",
                kernel.cols == expected_kernels[label]
                and len(forms) == expected_forms[label][0]
                and forms[0].rank() == expected_forms[label][1]
                and integer_certificate(forms[0]) == record["invariant_symmetric_form"],
            )
        else:
            projectors = multiplicity_free_projectors((a, b))
            components[label] = [
                restrict_to_projector((a, b), projector)
                for _, projector in projectors
            ]
            rebuilt_dimensions = [dimension for dimension, _ in projectors]
            check(
                f"{label} exact rational projectors",
                rebuilt_dimensions == record["projector_dimensions_Q"]
                and [integer_certificate(projector) for _, projector in projectors]
                == record["projectors"],
            )

    expected_component_dimensions = {
        "001": [4, 4, 24], "111": [4, 4, 24],
        "011": [4, 28], "101": [4, 28],
    }
    check(
        "all odd constituent dimensions",
        {label: [actions[0].rows for actions in blocks] for label, blocks in components.items()}
        == expected_component_dimensions,
    )

    # Exact-Q upper containers and exact intertwiners are rebuilt here; no
    # producer helper or stored dimension inference is used.
    inter = exact["intertwiners_Q"]
    paired = [
        ("001_to_111_componentwise_dual", "001", "111"),
        ("011_to_101_componentwise_dual", "011", "101"),
    ]
    for family, source_label, target_label in paired:
        records = inter[family]
        check(
            f"{family} rebuilt",
            len(records) == len(components[source_label])
            and all(
                verify_intertwiner_record(record, components[target_label][index], components[source_label][index])
                for index, record in enumerate(records)
            ),
        )
    check(
        "shared A3 rebuilt",
        verify_intertwiner_record(
            inter["shared_A3_001_to_011"], components["011"][0], components["001"][0]
        ),
    )
    check(
        "repeated B13 rebuilt",
        verify_intertwiner_record(
            inter["repeated_B13_100_to_110"], active["110"], active["100"]
        ),
    )
    second_direct = intertwiner_basis(
        traceless_actions(components["011"][0]),
        traceless_actions(components["001"][1]),
    )
    second_dual_source = tuple(-matrix.T for matrix in traceless_actions(components["001"][1]))
    second_dual = intertwiner_basis(
        traceless_actions(components["011"][0]), second_dual_source
    )
    check("second A3 has zero direct/dual Hom", not second_direct and not second_dual)

    # Literal raw-word image and joint minors in exactly the rebuilt Q bases.
    image_blocks = {
        label: [active[label]] if label in expected_kernels else [sectors[tuple(map(int, label))][:2]]
        for label in ("000", "001", "010", "011", "100", "101", "110", "111")
    }
    joint_blocks = {
        "010_plus_100": [active["010"], active["100"]],
        "100_plus_110": [active["100"], active["110"]],
        "001_plus_111": [sectors[(0, 0, 1)][:2], sectors[(1, 1, 1)][:2]],
        "011_plus_101": [sectors[(0, 1, 1)][:2], sectors[(1, 0, 1)][:2]],
        "001_plus_011": [sectors[(0, 0, 1)][:2], sectors[(0, 1, 1)][:2]],
    }
    image_observed = {
        label: raw_word_minor_clean(blocks, exact["literal_word_image_minors_Fp"][label])
        for label, blocks in image_blocks.items()
    }
    joint_observed = {
        label: raw_word_minor_clean(blocks, exact["literal_word_linkage_minors_Fp"][label])
        for label, blocks in joint_blocks.items()
    }
    image_ranks = {label: value["rank_Fp"] for label, value in image_observed.items()}
    joint_ranks = {label: value["rank_Fp"] for label, value in joint_observed.items()}
    check(
        "all eight actual selected image minors",
        image_ranks == {"000": 861, "001": 606, "010": 351, "011": 799,
                        "100": 351, "101": 799, "110": 351, "111": 606},
    )
    check(
        "all five actual selected linkage minors",
        joint_ranks == {"010_plus_100": 702, "100_plus_110": 351,
                        "001_plus_111": 606, "011_plus_101": 799,
                        "001_plus_011": 1389},
    )

    container = exact["containing_algebra_Q"]
    check(
        "complete exact-Q upper container",
        container["dimension_Q"] == 2952
        and container["derived_dimension_Q"] == 2951
        and container["factor_dimensions"] == [861, 351, 351, 575, 783, 15, 15]
        and container["factor_ranks"] == [21, 13, 13, 23, 27, 3, 3]
        and container["exact_Q_image_and_joint_ranks"]["sector_images"] == image_ranks
        and container["exact_Q_image_and_joint_ranks"]["joint_images"] == joint_ranks,
    )

    # Independently canonicalize the full and derived raw Pauli-word trees.
    pauli = PauliRecipeClosure(GOOD_PRIME)
    full = pauli.canonical(derived=False)
    derived = pauli.canonical(derived=True)
    stored_modular = data["modular_word_certificate"]
    for label, observed, expected_rank in (
        ("full_closure", full, 2952),
        ("derived_ideal", derived, 2951),
    ):
        stored = stored_modular[label]
        check(
            f"canonical raw Pauli tree {label}",
            observed["rank"] == stored["rank_Fp"] == expected_rank
            and observed["tree"] == stored["tree"]
            and observed["grading_block_ranks"] == stored["grading_block_ranks"]
            and observed["maximum_depth"] == stored["maximum_depth"]
            and observed["checksum"] == stored["tree_checksum_mod_1000000007"]
            and observed["support_orbits"] == stored_modular["support_orbits"],
        )

    # Rebuild the one global scalar projection of physical B.  No block is
    # primitive-normalized separately: one exact rational z is assembled first,
    # then compared with the artifact's single globally normalized integer.
    scalar_pieces = []
    traceless_checks = []
    for key in CHARACTERS:
        label = "".join(map(str, key))
        a, b, _gram, _inclusion = sectors[key]
        if label in components:
            projectors = multiplicity_free_projectors((a, b))
            scalar = sp.zeros(a.rows)
            for (_, projector), actions in zip(projectors, components[label]):
                scalar_value = actions[1].trace() / actions[1].rows
                scalar += scalar_value * projector
                traceless_checks.append(actions[1].trace() - actions[1].rows * scalar_value == 0)
        else:
            scalar = sp.zeros(a.rows)
            traceless_checks.append(active[label][1].trace() == 0)
        scalar_pieces.append(scalar)
    z = sp.diag(*scalar_pieces)
    z_integer = primitive_integer_matrix(z)
    stored_z_integer = decode_integer_certificate(exact["central_witness_global_integer"])
    scale = sp.Rational(exact["central_projection_equals_integer_times"])
    full_a = sp.diag(*(sectors[key][0] for key in CHARACTERS))
    full_b = sp.diag(*(sectors[key][1] for key in CHARACTERS))
    check("global central witness exact equality", z_integer == stored_z_integer and z == scale * stored_z_integer)
    check("global central witness nonzero", z != sp.zeros(256))
    check("global central witness commutes A", z * full_a == full_a * z)
    check("global central witness commutes B", z * full_b == full_b * z)
    check("B-z traceless on every simple constituent", all(traceless_checks))
    check(
        "B-z exact derived membership by same-basis equality",
        derived["rank"] == container["derived_dimension_Q"] == 2951
        and all(traceless_checks),
        "derived image is contained in the rebuilt block-traceless S and has dim S",
    )
    check(
        "z belongs to g and factor kernel is Qz",
        full["rank"] == container["dimension_Q"] == 2952
        and data["factor_map"]["joint_rank_Q"] == 2951
        and data["factor_map"]["kernel_dimension_Q"] == 1
        and data["factor_map"]["kernel_equals_center"] is True,
        "z=B-(B-z), with B-z in g'=S",
    )

    check(
        "absolute types without split-form overclaim",
        data["Levi_type"]["type"] == "D21 + B13 + B13 + A23 + A27 + A3 + A3"
        and data["Levi_type"]["rank"] == 103
        and "splitness" in data["Killing"]["reconstruction"]["rational_form_scope"]
        and data["Killing"]["rank_Q"] == 2951
        and data["Killing"]["nullity_Q"] == 1
        and "structural_2951_minor" not in data["Killing"]["reconstruction"],
    )
    check("all stored producer checks", all(item["passed"] for item in artifact["checks"]))
    print(f"WALL_SECONDS={time.monotonic() - started:.3f}")
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
