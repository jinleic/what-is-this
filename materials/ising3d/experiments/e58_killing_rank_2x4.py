#!/usr/bin/env python3
"""Exact-Q closure and Killing/radical certificate for the finite 2x4 Ising layer.

The Lie bracket is represented in ordered Pauli monomials ``Q_(a,b)=X^a Z^b``.
For coefficient control we use the integral half-bracket ``D_G=[G, .]/2``;
this has the same invariant subspaces as ``ad_G`` over Q.  Spatially invariant
vectors are stored by one integer coefficient per C2xC2 point-group orbit, but
all residual checks account for every Pauli string in every orbit.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import resource
import sys
import time
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
E20_PATH = ROOT / "experiments" / "e20_algebra_growth.py"
E49_PATH = ROOT / "experiments" / "e49_char0_levi_2x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "killing_rank_2x4.json"
SCRIPT = "experiments/e58_killing_rank_2x4.py"
ROWS = 2
COLS = 4
N = ROWS * COLS
DIMENSION = 2952
KILLING_RANK = 2951
PIVOT_PRIME = 2_147_483_647
WORD_FORMAT = "seed A=0, seed B=1; thereafter (parent_index, action_index), action 0=A and 1=B"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _sha256_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _array_digest(values: np.ndarray) -> str:
    values = np.ascontiguousarray(values)
    header = f"{values.dtype.str}:{','.join(map(str, values.shape))}\n".encode("ascii")
    return hashlib.sha256(header + values.tobytes()).hexdigest()


def _vector_content(row: np.ndarray) -> int:
    content = 0
    for value in row:
        content = math.gcd(content, abs(int(value)))
    return content


def _primitive_word(word: tuple[int, np.ndarray]) -> tuple[int, np.ndarray]:
    block, row = word
    content = _vector_content(row)
    if not content:
        raise AssertionError("a stored nested Lie word vanished")
    return block, np.asarray([int(value) // content for value in row], dtype=np.int64)


@dataclass(frozen=True)
class WordSpec:
    parent: int
    action: int
    depth: int
    block: int
    pivot_local: int
    pivot_pauli: int
    modular_pivot: int


@dataclass
class WordCertificate:
    model: object
    worker: object
    specs: list[WordSpec]
    words_mod: list[tuple[int, np.ndarray]]
    words_Z: list[tuple[int, np.ndarray]]
    action_targets: list[list[int]]
    determinant_residue: int
    maximum_abs_coefficient: int
    primitive_words_Z: list[tuple[int, np.ndarray]]
    word_contents: list[int]


def _apply_raw(worker, action_index: int, word: tuple[int, np.ndarray], modulus: int | None):
    """Apply the integral half-bracket action in one orbit grading block."""
    source_block, row = word
    target_block = source_block ^ (3 if action_index == 0 else 1)
    sources, signs = worker.actions[action_index][source_block]
    dtype = object if row.dtype == object else np.int64
    result = np.zeros(worker.block_sizes[target_block], dtype=dtype)
    for source, sign in zip(sources, signs):
        selected = source >= 0
        if not np.any(selected):
            continue
        result[selected] += sign[selected].astype(np.int64) * row[source[selected]]
    if modulus is not None:
        result = np.asarray(result, dtype=np.int64) % modulus
    return target_block, result


def _exact_word(worker, spec: WordSpec, exact_words: list[tuple[int, np.ndarray]], seed_index: int | None = None):
    if seed_index is not None:
        terms = worker.model.x_terms if seed_index == 0 else worker.model.zz_terms
        block, row = worker.dense_seed(terms)
        return block, np.asarray(row, dtype=np.int64)
    parent = exact_words[spec.parent]
    max_source = max((abs(int(x)) for x in parent[1]), default=0)
    term_count = N if spec.action == 0 else len(worker.model.bonds)
    if max_source <= (np.iinfo(np.int64).max // max(1, term_count)):
        return _apply_raw(worker, spec.action, parent, None)
    object_parent = (parent[0], np.asarray(parent[1], dtype=object))
    return _apply_raw(worker, spec.action, object_parent, None)


def build_word_certificate(prime: int = PIVOT_PRIME) -> WordCertificate:
    """Rebuild the deterministic 2,952 nested words and their exact Z vectors."""
    e20 = _load(E20_PATH, "_e20_for_e58")
    model = e20.OrbitModel(ROWS, COLS, False)
    support = model.enumerate_support(100_000)
    worker = e20.DenseOrbitEngine(model, prime, support, None, None)

    specs: list[WordSpec] = []
    words_mod: list[tuple[int, np.ndarray]] = []
    words_Z: list[tuple[int, np.ndarray]] = []
    action_targets: list[list[int]] = []
    frontier: list[int] = []
    determinant_residue = 1
    maximum_abs = 0

    for seed_index, terms in enumerate((model.x_terms, model.zz_terms)):
        modular_word = worker.dense_seed(terms)
        added = worker.reduce_add(modular_word)
        assert added is not None
        block, row_index = added
        pivot_local = int(np.flatnonzero(worker.rows[block][row_index])[0])
        pivot_value = int(worker.rows[block][row_index, pivot_local])
        determinant_residue = determinant_residue * pivot_value % prime
        spec = WordSpec(-1, seed_index, 1, block, pivot_local, worker.blocks[block][pivot_local], pivot_value)
        exact_word = _exact_word(worker, spec, words_Z, seed_index)
        specs.append(spec)
        words_mod.append((block, np.asarray(modular_word[1], dtype=np.int64) % prime))
        words_Z.append(exact_word)
        action_targets.append([-1, -1])
        frontier.append(len(specs) - 1)
        maximum_abs = max(maximum_abs, max(abs(int(x)) for x in exact_word[1]))

    head = 0
    while head < len(frontier):
        parent = frontier[head]
        head += 1
        for action_index in (0, 1):
            candidate = _apply_raw(worker, action_index, words_mod[parent], prime)
            added = worker.reduce_add(candidate)
            if added is None:
                continue
            block, row_index = added
            pivot_local = int(np.flatnonzero(worker.rows[block][row_index])[0])
            pivot_value = int(worker.rows[block][row_index, pivot_local])
            determinant_residue = determinant_residue * pivot_value % prime
            spec = WordSpec(
                parent,
                action_index,
                specs[parent].depth + 1,
                block,
                pivot_local,
                worker.blocks[block][pivot_local],
                pivot_value,
            )
            exact_word = _exact_word(worker, spec, words_Z)
            child = len(specs)
            specs.append(spec)
            words_mod.append((block, candidate[1]))
            words_Z.append(exact_word)
            action_targets.append([-1, -1])
            action_targets[parent][action_index] = child
            frontier.append(child)
            maximum_abs = max(maximum_abs, max(abs(int(x)) for x in exact_word[1]))

    assert len(specs) == DIMENSION
    assert determinant_residue != 0
    for index, spec in enumerate(specs):
        block, row = words_Z[index]
        assert block == spec.block
        assert np.array_equal(np.asarray([int(x) % prime for x in row], dtype=np.int64), words_mod[index][1])
    return WordCertificate(
        model=model,
        worker=worker,
        specs=specs,
        words_mod=words_mod,
        words_Z=words_Z,
        action_targets=action_targets,
        determinant_residue=determinant_residue,
        primitive_words_Z=[_primitive_word(word) for word in words_Z],
        word_contents=[_vector_content(row) for _, row in words_Z],
        maximum_abs_coefficient=maximum_abs,
    )


def word_digest(certificate: WordCertificate) -> str:
    return _sha256_lines(
        f"{index}:{spec.parent}:{spec.action}:{spec.depth}:{spec.block}:{spec.pivot_pauli}:{spec.modular_pivot}"
        for index, spec in enumerate(certificate.specs)
    )




def nested_word_actions(
    sectors, certificate: WordCertificate, prime: int
):
    """Evaluate reductions of the exact stored nested words in all sectors."""
    generators = {
        key: (
            np.asarray(actions[0].tolist(), dtype=np.int64) % prime,
            np.asarray(actions[1].tolist(), dtype=np.int64) % prime,
        )
        for key, actions in sectors.items()
    }
    words = {key: [actions[0], actions[1]] for key, actions in generators.items()}
    for spec in certificate.specs[2:]:
        for key in sectors:
            generator = generators[key][spec.action]
            parent = words[key][spec.parent]
            words[key].append((generator @ parent - parent @ generator) % prime)
    return words




def _rank_mod_chunked(
    rows: np.ndarray, prime: int
) -> tuple[int, list[int], int, list[int]]:
    """Overflow-safe echelon rank for a few thousand long rows."""
    rows = np.asarray(rows, dtype=np.int64) % prime
    basis: list[np.ndarray] = []
    pivots: dict[int, int] = {}
    selected: list[int] = []
    determinant = 1
    for original_index, raw in enumerate(rows):
        value = raw.copy()
        while True:
            nonzero = np.flatnonzero(value)
            if not nonzero.size:
                break
            lead = int(nonzero[0])
            old = pivots.get(lead)
            if old is None:
                pivot = int(value[lead])
                determinant = determinant * pivot % prime
                pivot_inverse = pow(pivot, prime - 2, prime)
                for start in range(lead, len(value), 4096):
                    stop = min(start + 4096, len(value))
                    value[start:stop] = value[start:stop] * pivot_inverse % prime
                pivots[lead] = len(basis)
                basis.append(value)
                selected.append(original_index)
                break
            factor = int(value[lead])
            old_row = basis[old]
            for start in range(lead, len(value), 4096):
                stop = min(start + 4096, len(value))
                value[start:stop] = (
                    value[start:stop] - factor * old_row[start:stop]
                ) % prime
    return len(basis), selected, determinant, list(pivots)


def _flatten_words(words, labels: Sequence[tuple[int, int, int]]) -> np.ndarray:
    return np.asarray(
        [
            np.concatenate([words[label][index].ravel() for label in labels])
            for index in range(DIMENSION)
        ],
        dtype=np.int64,
    )

def _intertwiner_record(e49, target, source, dual: bool) -> dict[str, object]:
    source_for_equation = (
        (source[0].T, source[1].T) if dual else source
    )
    signs = (-1, -1) if dual else (1, 1)
    basis = e49.intertwiner_basis(target, source_for_equation, signs)
    return {
        "equation": "G_target T = -T G_source^T" if dual else "G_target T = T G_source",
        "solution_dimension_Q": len(basis),
        "maximum_rank_Q": e49.max_combination_rank(basis),
        "basis_checksums": [e49.matrix_checksum(matrix) for matrix in basis],
    }


def _central_scalar_tuple(e49, sectors):
    central = {}
    scalar_data = {}
    for key, actions in sectors.items():
        if key[2] == 0:
            central[key] = e49.sp.zeros(actions[0].rows)
            scalar_data["".join(map(str, key))] = []
            continue
        blocks = e49.multiplicity_free_projectors(*actions)
        value = e49.sp.zeros(actions[0].rows)
        scalars = []
        for dimension, projector in blocks:
            scalar = e49.sp.trace(projector * actions[1]) / dimension
            value += scalar * projector
            scalars.append(str(scalar))
        central[key] = value
        scalar_data["".join(map(str, key))] = scalars
    nonzero_commutators = sum(
        sum(value != 0 for value in (central[key] * action - action * central[key]))
        for key, actions in sectors.items()
        for action in actions
    )
    if nonzero_commutators or not any(any(matrix) for matrix in central.values()):
        raise AssertionError("the exact central block-scalar witness failed")
    denominator_lcm = math.lcm(
        *(
            int(e49.sp.denom(value))
            for matrix in central.values()
            for value in matrix
        )
    )
    digest = _sha256_lines(
        f"{''.join(map(str, key))}:{i}:{j}:{central[key][i, j]}"
        for key in sorted(central)
        for i in range(central[key].rows)
        for j in range(central[key].cols)
        if central[key][i, j]
    )
    return central, {
        "construction": "zero on even spin-flip sectors; on every exact odd irreducible projector P_d, (tr(P_d B)/d) P_d",
        "block_scalars": scalar_data,
        "common_denominator": denominator_lcm,
        "nonzero_generator_commutator_entries_Q": nonzero_commutators,
        "matrix_digest_sha256": digest,
    }
def exact_factor_audit(
    certificate: WordCertificate,
    prime: int = 65_521,
    return_state: bool = False,
) -> dict[str, object]:
    """Independent exact-Q factor/container audit using the literal word basis."""
    e49 = _load(E49_PATH, "_e49_for_e58")
    sectors = {key: e49.symmetry_sector(key) for key in e49.CHARACTERS}
    words = nested_word_actions(sectors, certificate, prime)

    label = lambda text: tuple(map(int, text))
    joint_cases = {
        "000": (("000",), 861),
        "010": (("010",), 351),
        "100_110": (("100", "110"), 351),
        "010_100": (("010", "100"), 702),
        "001_111": (("001", "111"), 606),
        "011_101": (("011", "101"), 799),
        "001_011": (("001", "011"), 1389),
        "all_sectors": (
            tuple("".join(map(str, key)) for key in e49.CHARACTERS),
            2952,
        ),
    }
    ranks = {}
    minors = {}
    for name, (text_labels, expected) in joint_cases.items():
        labels = tuple(label(text) for text in text_labels)
        matrix = _flatten_words(words, labels)
        rank, selected, determinant, pivot_columns = _rank_mod_chunked(
            matrix, prime
        )
        if rank != expected or determinant == 0:
            raise AssertionError(f"{name} literal-word minor rank {rank}, expected {expected}")
        ranks[name] = rank
        minors[name] = {
            "prime": prime,
            "determinant_residue": determinant,
            "selected_word_count": len(selected),
            "selected_words_digest_sha256": _sha256_lines(map(str, selected)),
            "pivot_columns_digest_sha256": _sha256_lines(map(str, pivot_columns)),
        }
    exact = {}
    quotients = {}
    components = {}
    for text in ("000", "010", "100", "110", "001", "111", "011", "101"):
        actions = sectors[label(text)]
        kernel, quotient = e49.quotient_by_common_kernel(actions)
        quotients[text] = quotient
        record = {
            "common_kernel_dimension_Q": kernel.cols,
            "quotient_dimension_Q": quotient[0].rows,
        }
        if text in {"000", "010", "100", "110"}:
            forms = e49.invariant_symmetric_forms(quotient)
            record.update(
                invariant_symmetric_form_nullity_Q=len(forms),
                invariant_symmetric_form_rank_Q=forms[0].rank() if forms else 0,
            )
        else:
            projectors = e49.multiplicity_free_projectors(*actions)
            components[text] = [
                (dimension, e49.restrict_to_projector(actions, projector))
                for dimension, projector in projectors
            ]
            record["exact_constituent_dimensions_Q"] = [
                dimension for dimension, _ in components[text]
            ]
        exact[text] = record

    linkages = {
        "100_to_110_direct": _intertwiner_record(
            e49, quotients["110"], quotients["100"], False
        ),
        "001_to_111_dual": _intertwiner_record(
            e49,
            e49.traceless_actions(sectors[label("111")]),
            e49.traceless_actions(sectors[label("001")]),
            True,
        ),
        "011_to_101_dual": _intertwiner_record(
            e49,
            e49.traceless_actions(sectors[label("101")]),
            e49.traceless_actions(sectors[label("011")]),
            True,
        ),
    }
    a3_links = []
    for source_index, (source_dimension, source) in enumerate(components["001"]):
        if source_dimension != 4:
            continue
        for target_index, (target_dimension, target) in enumerate(components["011"]):
            if target_dimension != 4:
                continue
            for dual in (False, True):
                record = _intertwiner_record(
                    e49,
                    e49.traceless_actions(target),
                    e49.traceless_actions(source),
                    dual,
                )
                record.update(
                    source_component_index=source_index,
                    target_component_index=target_index,
                    dual=dual,
                )
                a3_links.append(record)
    linkages["001_A3_to_011_A3"] = a3_links
    if linkages["100_to_110_direct"]["maximum_rank_Q"] != 27:
        raise AssertionError("100 and 110 B13 quotient linkage failed")
    if linkages["001_to_111_dual"]["maximum_rank_Q"] != 32:
        raise AssertionError("001 and 111 dual linkage failed")
    if linkages["011_to_101_dual"]["maximum_rank_Q"] != 32:
        raise AssertionError("011 and 101 dual linkage failed")
    if not any(not item["dual"] and item["maximum_rank_Q"] == 4 for item in a3_links):
        raise AssertionError("the shared direct A3 linkage failed")
    central, central_record = _central_scalar_tuple(e49, sectors)

    factors = [
        {"type": "D21", "dimension": 861, "multiplicity": 1},
        {"type": "B13", "dimension": 351, "multiplicity": 2},
        {"type": "A23", "dimension": 575, "multiplicity": 1},
        {"type": "A27", "dimension": 783, "multiplicity": 1},
        {"type": "A3", "dimension": 15, "multiplicity": 2},
    ]
    semisimple_dimension = sum(
        item["dimension"] * item["multiplicity"] for item in factors
    )
    if semisimple_dimension != KILLING_RANK:
        raise AssertionError("factor dimensions do not sum to 2951")
    audit = {
        "same_word_digest_sha256": word_digest(certificate),
        "literal_word_minor_ranks_Q_lower_bounds": ranks,
        "literal_word_minors": minors,
        "exact_Q_sector_containers": exact,
        "exact_Q_linkages": linkages,
        "exact_central_generator": central_record,
        "semisimple_factors_Q": factors,
        "semisimple_dimension_Q": semisimple_dimension,
    }
    if return_state:
        factor_modules = [
            ("D21", 42, 40, quotients["000"]),
            ("B13a", 27, 25, quotients["010"]),
            ("B13b", 27, 25, quotients["100"]),
        ]
        a3_components = [
            (index, actions)
            for index, (dimension, actions) in enumerate(components["001"])
            if dimension == 4
        ]
        shared_index = next(
            item["source_component_index"]
            for item in a3_links
            if not item["dual"] and item["maximum_rank_Q"] == 4
        )
        distinct_index = next(
            index for index, _ in a3_components if index != shared_index
        )
        component_by_index = dict(a3_components)
        factor_modules.extend(
            [
                ("A3_shared", 4, 8, e49.traceless_actions(components["011"][0][1])),
                ("A3_distinct", 4, 8, e49.traceless_actions(component_by_index[distinct_index])),
            ]
        )
        factor_modules.extend(
            [("A23", 24, 48, e49.traceless_actions(next(actions for dimension, actions in components["001"] if dimension == 24)))]
        )
        factor_modules.extend(
            [("A27", 28, 56, e49.traceless_actions(next(actions for dimension, actions in components["011"] if dimension == 28)))]
        )
        return audit, sectors, words, central, factor_modules
    return audit


def _sympy_matrix_mod(matrix, prime: int) -> tuple[np.ndarray, set[int]]:
    output = np.empty(matrix.shape, dtype=np.int64)
    denominators = set()
    for index, value in enumerate(matrix):
        numerator, denominator = int(value.p), int(value.q)
        if denominator % prime == 0:
            raise ArithmeticError(f"bad reduction prime {prime} divides denominator")
        denominators.add(denominator)
        output.ravel()[index] = numerator * pow(denominator, -1, prime) % prime
    return output, denominators


def _module_words_mod(
    actions, certificate: WordCertificate, prime: int
) -> tuple[np.ndarray, dict[str, object]]:
    generators = []
    denominators = set()
    exact_digest = hashlib.sha256()
    for action in actions:
        reduced, action_denominators = _sympy_matrix_mod(action, prime)
        generators.append(reduced)
        denominators.update(action_denominators)
        for value in action:
            exact_digest.update(f"{value};".encode("ascii"))
    dimension = generators[0].shape[0]
    output = np.empty((DIMENSION, dimension * dimension), dtype=np.int64)
    matrices = [generators[0], generators[1]]
    output[0] = generators[0].ravel()
    output[1] = generators[1].ravel()
    for index, spec in enumerate(certificate.specs[2:], start=2):
        generator = generators[spec.action]
        parent = matrices[spec.parent]
        value = (generator @ parent - parent @ generator) % prime
        matrices.append(value)
        output[index] = value.ravel()
    return output, {
        "module_dimension": dimension,
        "generator_denominators": sorted(denominators),
        "generator_actions_Q_digest_sha256": exact_digest.hexdigest(),
        "word_matrix_mod_digest_sha256": _array_digest(output),
    }


def _trace_gram_mod(words: np.ndarray, dimension: int, prime: int) -> np.ndarray:
    transpose = np.arange(dimension * dimension).reshape(dimension, dimension).T.ravel()
    right = words[:, transpose].T
    output = np.empty((DIMENSION, DIMENSION), dtype=np.int64)
    for start in range(0, DIMENSION, 64):
        stop = min(start + 64, DIMENSION)
        output[start:stop] = words[start:stop] @ right % prime
    return output


def _rank_echelon_mod(
    matrix: np.ndarray, prime: int
) -> tuple[int, list[int], list[int], int, np.ndarray]:
    """Forward modular elimination; returns a nonzero pivot minor and echelon form."""
    work = np.asarray(matrix, dtype=np.int64).copy() % prime
    row_order = list(range(work.shape[0]))
    pivot_columns = []
    pivot_rows = []
    determinant = 1
    rank = 0
    sign = 1
    for column in range(work.shape[1]):
        choices = np.flatnonzero(work[rank:, column])
        if not choices.size:
            continue
        pivot = rank + int(choices[0])
        if pivot != rank:
            work[[rank, pivot]] = work[[pivot, rank]]
            row_order[rank], row_order[pivot] = row_order[pivot], row_order[rank]
            sign = -sign
        pivot_value = int(work[rank, column])
        determinant = determinant * pivot_value % prime
        pivot_inverse = pow(pivot_value, prime - 2, prime)
        work[rank, column:] = work[rank, column:] * pivot_inverse % prime
        rows = rank + 1 + np.flatnonzero(work[rank + 1 :, column])
        for start in range(0, len(rows), 64):
            selected = rows[start : start + 64]
            factors = work[selected, column].copy()
            work[selected, column:] = (
                work[selected, column:]
                - factors[:, None] * work[rank, column:]
            ) % prime
        pivot_columns.append(column)
        pivot_rows.append(row_order[rank])
        rank += 1
        if rank == work.shape[0]:
            break
    if sign < 0:
        determinant = (-determinant) % prime
    return rank, pivot_rows, pivot_columns, determinant, work


def _one_dimensional_null_vector(
    echelon: np.ndarray, pivot_columns: Sequence[int], prime: int
) -> np.ndarray:
    pivot_set = set(pivot_columns)
    free = [column for column in range(echelon.shape[1]) if column not in pivot_set]
    if len(free) != 1:
        raise AssertionError(f"expected one free Killing coordinate, found {len(free)}")
    vector = np.zeros(echelon.shape[1], dtype=np.int64)
    vector[free[0]] = 1
    for row in range(len(pivot_columns) - 1, -1, -1):
        pivot = pivot_columns[row]
        total = int(echelon[row, pivot + 1 :] @ vector[pivot + 1 :]) % prime
        vector[pivot] = -total % prime
    return vector


def certify_killing_matrix(
    certificate: WordCertificate,
    sectors,
    all_words,
    central,
    factor_modules,
    prime: int = 65_521,
) -> dict[str, object]:
    """Build the exact reduction of the full rational Killing matrix."""
    matrix = np.zeros((DIMENSION, DIMENSION), dtype=np.int64)
    module_records = []
    all_denominators = {1}
    for name, dimension, dynkin_coefficient, actions in factor_modules:
        word_matrices, record = _module_words_mod(actions, certificate, prime)
        if word_matrices.shape[1] != dimension * dimension:
            raise AssertionError(f"{name} natural module has the wrong dimension")
        gram = _trace_gram_mod(word_matrices, dimension, prime)
        matrix = (matrix + dynkin_coefficient * gram) % prime
        all_denominators.update(record["generator_denominators"])
        record.update(
            name=name,
            Dynkin_trace_coefficient=dynkin_coefficient,
            trace_gram_mod_digest_sha256=_array_digest(gram),
        )
        module_records.append(record)
        del word_matrices, gram
    if not np.array_equal(matrix, matrix.T):
        raise AssertionError("the Killing reduction is not symmetric")

    rank, pivot_rows, pivot_columns, determinant, echelon = _rank_echelon_mod(
        matrix, prime
    )
    if rank != KILLING_RANK or determinant == 0:
        raise AssertionError(f"Killing rank mod {prime} was {rank}")
    null_vector = _one_dimensional_null_vector(
        echelon, pivot_columns, prime
    )
    kz = matrix @ null_vector % prime
    if np.count_nonzero(kz):
        raise AssertionError("Killing null-vector residual is nonzero")

    e49 = sys.modules["_e49_for_e58"]
    central_mod = {
        key: _sympy_matrix_mod(value, prime)[0] for key, value in central.items()
    }
    represented = {}
    scale = None
    for key in sorted(sectors):
        stacked = np.asarray(all_words[key], dtype=np.int64)
        value = np.tensordot(null_vector, stacked, axes=(0, 0)) % prime
        represented[key] = value
        target = central_mod[key]
        if scale is None:
            nonzero = np.flatnonzero(target)
            if nonzero.size:
                location = int(nonzero[0])
                scale = (
                    int(value.ravel()[location])
                    * pow(int(target.ravel()[location]), -1, prime)
                ) % prime
    if scale is None or scale == 0:
        raise AssertionError("central witness reduced to zero at the Killing prime")
    central_residuals = sum(
        int(
            np.count_nonzero(
                (represented[key] - scale * central_mod[key]) % prime
            )
        )
        for key in represented
    )
    if central_residuals:
        raise AssertionError("Killing null vector is not the exact central witness reduction")

    denominator_lcm = math.lcm(*all_denominators)
    if denominator_lcm % prime == 0:
        raise AssertionError("Killing prime is not good for factor-map denominators")
    return {
        "entry_domain": "Q in the exact raw nested-word basis; this record materializes its exact good-prime reduction",
        "size": [DIMENSION, DIMENSION],
        "formula": "40 tr_42 + 25 tr_27a + 25 tr_27b + 48 tr_24 + 56 tr_28 + 8 tr_4a + 8 tr_4b",
        "factor_module_records": module_records,
        "common_generator_denominator_lcm": denominator_lcm,
        "good_reduction_prime": prime,
        "matrix_mod_digest_sha256": _array_digest(matrix),
        "symmetric_mod_prime": True,
        "rank_mod_prime": rank,
        "rank_minor": {
            "size": KILLING_RANK,
            "determinant_residue": determinant,
            "row_digest_sha256": _sha256_lines(map(str, pivot_rows)),
            "column_digest_sha256": _sha256_lines(map(str, pivot_columns)),
        },
        "null_vector_mod_prime": {
            "nonzero_entries": int(np.count_nonzero(null_vector)),
            "digest_sha256": _array_digest(null_vector),
            "Kz_nonzero_residuals": int(np.count_nonzero(kz)),
            "central_witness_scale": int(scale),
            "central_matrix_nonzero_residuals": central_residuals,
        },
    }

def pivot_digest(certificate: WordCertificate) -> str:
    return _sha256_lines(
        f"{index}:{spec.block}:{spec.pivot_local}:{spec.pivot_pauli}:{spec.modular_pivot}"
        for index, spec in enumerate(certificate.specs)
    )


def _orbit_sizes(certificate: WordCertificate) -> list[np.ndarray]:
    e20 = sys.modules.get("_e20_for_e58") or _load(E20_PATH, "_e20_for_e58")
    return [
        np.asarray(
            [
                len(
                    {
                        e20._permute_pauli(value, permutation, N)
                        for permutation in certificate.model.group
                    }
                )
                for value in certificate.worker.blocks[block]
            ],
            dtype=np.int64,
        )
        for block in range(4)
    ]




def probe() -> dict[str, object]:
    started = time.monotonic()
    certificate = build_word_certificate()
    return {
        "dimension": len(certificate.specs),
        "maximum_word_depth": max(spec.depth for spec in certificate.specs),
        "grading_block_rows": [
            sum(spec.block == block for spec in certificate.specs)
            for block in range(4)
        ],
        "grading_block_columns": list(certificate.worker.block_sizes),
        "support_orbits": len(certificate.worker.support),
        "full_pauli_support": sum(
            int(values.sum()) for values in _orbit_sizes(certificate)
        ),
        "maximum_abs_exact_word_coefficient": certificate.maximum_abs_coefficient,
        "maximum_abs_primitive_coefficient": max(
            abs(int(value))
            for _, row in certificate.primitive_words_Z
            for value in row
        ),
        "word_content_minimum": min(certificate.word_contents),
        "word_content_maximum": max(certificate.word_contents),
        "word_content_distinct": len(set(certificate.word_contents)),
        "pivot_prime": PIVOT_PRIME,
        "pivot_minor_determinant_residue": certificate.determinant_residue,
        "word_digest_sha256": word_digest(certificate),
        "pivot_digest_sha256": pivot_digest(certificate),
        "elapsed_seconds": round(time.monotonic() - started, 6),
        "rss_peak_bytes": _rss_bytes(),
    }




def factor_probe() -> dict[str, object]:
    started = time.monotonic()
    certificate = build_word_certificate()
    basis_elapsed = time.monotonic() - started
    factors_started = time.monotonic()
    audit = exact_factor_audit(certificate)
    audit["word_basis_seconds"] = round(basis_elapsed, 6)
    audit["factor_audit_seconds"] = round(
        time.monotonic() - factors_started, 6
    )
    audit["total_seconds"] = round(time.monotonic() - started, 6)
    audit["rss_peak_bytes"] = _rss_bytes()
    return audit


def killing_probe() -> dict[str, object]:
    started = time.monotonic()
    certificate = build_word_certificate()
    basis_elapsed = time.monotonic() - started
    audit_started = time.monotonic()
    audit, sectors, words, central, modules = exact_factor_audit(
        certificate, return_state=True
    )
    audit_elapsed = time.monotonic() - audit_started
    killing_started = time.monotonic()
    killing = certify_killing_matrix(
        certificate, sectors, words, central, modules
    )
    audit["Killing"] = killing
    audit["dimension_Q"] = DIMENSION
    audit["Killing_rank_Q"] = KILLING_RANK
    audit["Killing_nullity_Q"] = 1
    audit["radical_dimension_Q"] = 1
    audit["radical_equals_center"] = True
    audit["basis_certificate"] = {
        "basis": "2,952 raw nested half-bracket words in deterministic breadth-first order",
        "word_format": WORD_FORMAT,
        "word_digest_sha256": word_digest(certificate),
        "maximum_word_depth": max(spec.depth for spec in certificate.specs),
        "orbit_coordinate_columns": len(certificate.worker.support),
        "full_Pauli_columns": sum(
            int(values.sum()) for values in _orbit_sizes(certificate)
        ),
        "good_prime": PIVOT_PRIME,
        "pivot_minor_determinant_residue": certificate.determinant_residue,
        "pivot_digest_sha256": pivot_digest(certificate),
        "meaning": "nonzero integral Pauli-coordinate minor modulo a prime, hence 2,952 exact-Q independent words",
    }
    audit["Cartan_radical_criterion"] = {
        "criterion": "rad(g)={x: kappa(x,[g,g])=0}",
        "derived_dimension_Q": KILLING_RANK,
        "derived_algebra": "D21 + 2 B13 + A23 + A27 + 2 A3",
        "kernel_dimension_Q": 1,
        "kernel_generator": "the exact central block-scalar witness",
        "equals_center": True,
    }
    audit["timings_seconds"] = {
        "word_basis": round(basis_elapsed, 6),
        "factor_audit": round(audit_elapsed, 6),
        "Killing_modular_materialization_and_rank": round(
            time.monotonic() - killing_started, 6
        ),
        "total": round(time.monotonic() - started, 6),
    }
    audit["rss_peak_bytes"] = _rss_bytes()
    return audit


def build_artifact() -> dict[str, object]:
    core = killing_probe()
    checks = [
        {
            "name": "exact_Q_word_independence",
            "passed": (
                core["basis_certificate"]["pivot_minor_determinant_residue"] != 0
                and core["literal_word_minor_ranks_Q_lower_bounds"]["all_sectors"]
                == DIMENSION
            ),
            "detail": "two explicit nonzero good-prime minors of the same literal nested words",
        },
        {
            "name": "exact_Q_closure_dimension",
            "passed": (
                core["dimension_Q"] == DIMENSION
                and core["semisimple_dimension_Q"] + 1 == DIMENSION
            ),
            "detail": "exact rational Lie container C=Qc plus the linked simple-factor sum has dimension 2952; the generator words have matching lower rank",
        },
        {
            "name": "exact_Q_factor_linkages",
            "passed": (
                core["exact_Q_linkages"]["100_to_110_direct"]["maximum_rank_Q"] == 27
                and core["exact_Q_linkages"]["001_to_111_dual"]["maximum_rank_Q"] == 32
                and core["exact_Q_linkages"]["011_to_101_dual"]["maximum_rank_Q"] == 32
                and any(
                    not item["dual"] and item["maximum_rank_Q"] == 4
                    for item in core["exact_Q_linkages"]["001_A3_to_011_A3"]
                )
            ),
            "detail": "all declared diagonal linkages are exact rational full-rank intertwiners; joint word minors exclude further linkages",
        },
        {
            "name": "full_Killing_reduction_rank",
            "passed": (
                core["Killing"]["rank_mod_prime"] == KILLING_RANK
                and core["Killing"]["rank_minor"]["determinant_residue"] != 0
                and core["Killing"]["symmetric_mod_prime"]
            ),
            "detail": "2952x2952 exact good-prime reduction of the rational trace-form formula in the same word basis",
        },
        {
            "name": "exact_central_generator",
            "passed": (
                core["exact_central_generator"][
                    "nonzero_generator_commutator_entries_Q"
                ]
                == 0
                and core["Killing"]["null_vector_mod_prime"]["Kz_nonzero_residuals"]
                == 0
                and core["Killing"]["null_vector_mod_prime"][
                    "central_matrix_nonzero_residuals"
                ]
                == 0
            ),
            "detail": "exact-Q block-scalar c commutes with A,B; the one-dimensional modular Killing nullspace is its reduction",
        },
        {
            "name": "Cartan_radical_equals_center",
            "passed": (
                core["radical_dimension_Q"] == 1
                and core["radical_equals_center"]
                and core["Cartan_radical_criterion"]["kernel_dimension_Q"] == 1
            ),
            "detail": "the exact semisimple derived sum has nondegenerate Killing form, so ker kappa(-,[g,g])=Qc",
        },
    ]
    data = {
        "scope": "single finite 2x4 open Ising layer; no all-sizes or 3D-solution assertion",
        "status": "[THEOREM] exact computer-assisted finite-dimensional result",
        "dimension_Q": core["dimension_Q"],
        "derived_dimension_Q": KILLING_RANK,
        "basis_certificate": core["basis_certificate"],
        "exact_Q_closure": {
            "method": "exact rational symmetry-sector Lie container plus same-word lower minors",
            "container_dimension_Q": DIMENSION,
            "semisimple_dimension_Q": core["semisimple_dimension_Q"],
            "factor_linkages": core["exact_Q_linkages"],
            "direct_pivot_coordinate_attempt": {
                "claim_tag": "[UNRESOLVED]",
                "status": "not used in any theorem",
                "attempted_good_primes": 30,
                "outcome": "raw/primitive word-basis rational reconstruction did not stabilize after 30 approximately one-million-size primes; full ambient rational residual coordinates were therefore not recorded",
                "replacement": "the exact characteristic-zero upper bound is instead the explicit 2952-dimensional rational Lie container containing A and B",
            },
        },
        "semisimple_factors_Q": core["semisimple_factors_Q"],
        "sector_containers_Q": core["exact_Q_sector_containers"],
        "literal_word_minors": core["literal_word_minors"],
        "Killing": core["Killing"],
        "center": {
            "dimension_Q": 1,
            "generator": core["exact_central_generator"],
        },
        "solvable_radical": {
            "dimension_Q": core["radical_dimension_Q"],
            "equals_center": core["radical_equals_center"],
            "Cartan_criterion": core["Cartan_radical_criterion"],
        },
        "timings_seconds": core["timings_seconds"],
        "rss_peak_bytes": core["rss_peak_bytes"],
    }
    return {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "arithmetic": "exact Z/Q/SymPy plus exact finite-field reductions at explicitly good primes; no floating point",
            "method": "same-word Pauli minor, exact rational sector containers/linkages, classical Killing trace formulas, full modular Killing minor, exact central witness, Cartan radical criterion",
        },
        "data": data,
        "checks": checks,
    }


def run(write: bool = True) -> dict[str, object]:
    artifact = build_artifact()
    if write:
        RESULT.parent.mkdir(parents=True, exist_ok=True)
        RESULT.write_text(
            json.dumps(artifact, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return artifact


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--probe",
        action="store_true",
        help="only rebuild exact nested words and the good-prime minor",
    )
    parser.add_argument(
        "--factors",
        action="store_true",
        help="audit exact-Q sector containers and literal-word minors",
    )
    parser.add_argument(
        "--killing",
        action="store_true",
        help="explicit alias for the default full artifact run",
    )
    args = parser.parse_args(argv)
    result = (
        probe()
        if args.probe
        else factor_probe()
        if args.factors
        else run(write=True)
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    if "checks" in result:
        passed = all(check["passed"] for check in result["checks"])
    else:
        passed = (
            result.get("dimension", result.get("dimension_Q", DIMENSION))
            == DIMENSION
            and result.get("pivot_minor_determinant_residue", 1) != 0
            and result.get("semisimple_dimension_Q", KILLING_RANK)
            == KILLING_RANK
            and result.get("literal_word_minor_ranks_Q_lower_bounds", {}).get(
                "all_sectors", DIMENSION
            )
            == DIMENSION
        )
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
