#!/usr/bin/env python3
"""Exact characteristic-zero quotient images of the open 2x4 Ising layer.

The decisive lower certificates are modular minors of literal right-normed Lie
words evaluated in the same rational quotient bases used for the exact upper
bounds.  A modular rank is used only as a lower bound over Q.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e59_levi_images_2x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "levi_images_2x4.json"
E49_PATH = ROOT / "experiments" / "e49_char0_levi_2x4.py"
OA_RESULT = ROOT / "results" / "algebra_structure" / "oa_quotient.json"
GOOD_PRIME = 65_521
ORTHOGONAL_TARGETS = {"000": (42, 861, "D21"), "010": (27, 351, "B13"),
                      "100": (27, 351, "B13"), "110": (27, 351, "B13")}
REDUCTIVE_TARGETS = {
    "001": ([4, 4, 24], 606, ["A3", "A3", "A23"]),
    "111": ([4, 4, 24], 606, ["A3", "A3", "A23"]),
    "011": ([4, 28], 799, ["A3", "A27"]),
    "101": ([4, 28], 799, ["A3", "A27"]),
}


def load_e49():
    spec = importlib.util.spec_from_file_location("_e49_for_e59", E49_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


E49 = load_e49()


def label_key(label: str) -> tuple[int, int, int]:
    return tuple(map(int, label))


def matrix_denominator_lcm(matrix: sp.Matrix) -> int:
    return int(sp.ilcm(*(sp.denom(value) for value in matrix))) if matrix else 1


def rational_matrix_mod(matrix: sp.Matrix, prime: int) -> np.ndarray:
    output = np.empty(matrix.shape, dtype=np.int64)
    for index, value in enumerate(matrix):
        value = sp.Rational(value)
        numerator, denominator = int(value.p), int(value.q)
        if denominator % prime == 0:
            raise AssertionError("rational matrix has a denominator divisible by the good prime")
        output.flat[index] = numerator * pow(denominator, -1, prime) % prime
    return output


def rational_matrix_checksum(matrix: sp.Matrix, store_entries: bool = False) -> dict:
    entries = [sp.Rational(value) for value in matrix]
    encoded = ";".join(f"{int(value.p)}/{int(value.q)}" for value in entries).encode("ascii")
    answer = {
        "rows": matrix.rows,
        "columns": matrix.cols,
        "nonzero_entries": sum(value != 0 for value in entries),
        "denominator_lcm": matrix_denominator_lcm(matrix),
        "sha256_row_major_rationals": hashlib.sha256(encoded).hexdigest(),
    }
    if store_entries:
        answer["entries_row_major"] = [
            int(value) if value.q == 1 else f"{int(value.p)}/{int(value.q)}"
            for value in entries
        ]
    return answer




def rank_mod(matrix: np.ndarray, prime: int) -> int:
    work = np.asarray(matrix, dtype=np.int64).copy() % prime
    row = 0
    for column in range(work.shape[1]):
        choices = np.flatnonzero(work[row:, column])
        if not choices.size:
            continue
        pivot = row + int(choices[0])
        work[[row, pivot]] = work[[pivot, row]]
        work[row] = work[row] * pow(int(work[row, column]), -1, prime) % prime
        indices = np.flatnonzero(work[:, column])
        indices = indices[indices != row]
        for start in range(0, len(indices), 128):
            block = indices[start:start + 128]
            work[block] = (work[block] - work[block, column, None] * work[row]) % prime
        row += 1
        if row == work.shape[0]:
            break
    return row


def block_diagonal(actions: list[tuple[sp.Matrix, sp.Matrix]]) -> tuple[sp.Matrix, sp.Matrix]:
    return tuple(sp.diag(*(pair[which] for pair in actions)) for which in range(2))


def audited_common_kernel_quotient(actions: tuple[sp.Matrix, sp.Matrix]):
    dimension = actions[0].rows
    stacked = sp.Matrix.vstack(*actions)
    kernel_columns = stacked.nullspace()
    kernel = sp.Matrix.hstack(*kernel_columns) if kernel_columns else sp.zeros(dimension, 0)
    assert stacked * kernel == sp.zeros(2 * dimension, kernel.cols)

    columns = list(kernel_columns)
    current_rank = kernel.cols
    identity = sp.eye(dimension)
    for column in range(dimension):
        candidate = sp.Matrix.hstack(*columns, identity[:, column])
        candidate_rank = candidate.rank()
        if candidate_rank > current_rank:
            columns.append(identity[:, column])
            current_rank = candidate_rank
        if current_rank == dimension:
            break
    change = sp.Matrix.hstack(*columns)
    assert change.rank() == dimension
    inverse = change.inv()
    assert change * inverse == identity == inverse * change

    conjugated = tuple(inverse * generator * change for generator in actions)
    kernel_dimension = kernel.cols
    for generator, transformed in zip(actions, conjugated):
        assert generator * kernel == sp.zeros(dimension, kernel_dimension)
        assert change * transformed == generator * change
        assert transformed[:, :kernel_dimension] == sp.zeros(dimension, kernel_dimension)
    quotient = tuple(
        transformed[kernel_dimension:, kernel_dimension:] for transformed in conjugated
    )

    prior_kernel, prior_quotient = E49.quotient_by_common_kernel(actions)
    assert prior_kernel.columnspace() == kernel.columnspace()
    assert all(left == right for left, right in zip(quotient, prior_quotient))
    return kernel, change, inverse, conjugated, quotient


def lie_word_minor_certificate(
    rational_generators: tuple[sp.Matrix, sp.Matrix],
    target_dimension: int,
    prime: int = GOOD_PRIME,
) -> dict:
    """Find a nonzero modular minor of literal rational Lie-word matrices.

    Accepted raw rows are reductions of A, B, or [A,w], [B,w] for an earlier
    accepted raw word w.  Echelon rows select the words but never replace the
    stored raw rows, so every selected row has a prime-independent rational
    recipe.  The pivot columns, in discovery order, form a square minor whose
    determinant is the product of the unnormalised insertion pivots.
    """
    dimensions = {generator.rows for generator in rational_generators}
    assert len(dimensions) == 1
    dimension = dimensions.pop()
    assert all(generator.cols == dimension for generator in rational_generators)
    denominator_lcms = [matrix_denominator_lcm(generator) for generator in rational_generators]
    assert all(value % prime for value in denominator_lcms)
    assert dimension * (prime - 1) ** 2 < np.iinfo(np.int64).max
    generators = tuple(rational_matrix_mod(generator, prime) for generator in rational_generators)

    echelon_rows: list[np.ndarray] = []
    pivot_to_row: dict[int, int] = {}
    raw_rows: list[np.ndarray] = []
    recipes: list[list[int | str]] = []
    depths: list[int] = []
    pivot_coordinates: list[int] = []
    insertion_pivots: list[int] = []
    digest = hashlib.sha256()

    def add(raw: np.ndarray, recipe: list[int | str], depth: int) -> bool:
        raw = np.asarray(raw, dtype=np.int64).reshape(-1) % prime
        work = raw.copy()
        while True:
            nonzero = np.flatnonzero(work)
            if not nonzero.size:
                return False
            lead = int(nonzero[0])
            previous = pivot_to_row.get(lead)
            if previous is None:
                pivot = int(work[lead])
                normalised = work * pow(pivot, -1, prime) % prime
                pivot_to_row[lead] = len(echelon_rows)
                echelon_rows.append(normalised)
                raw_rows.append(raw)
                recipes.append(recipe)
                depths.append(depth)
                pivot_coordinates.append(lead)
                insertion_pivots.append(pivot)
                digest.update(np.asarray(raw, dtype="<u4").tobytes())
                return True
            work = (work - int(work[lead]) * echelon_rows[previous]) % prime

    assert add(generators[0], ["seed", 0], 1)
    assert add(generators[1], ["seed", 1], 1)
    head = 0
    while head < len(raw_rows) and len(raw_rows) < target_dimension:
        parent = raw_rows[head].reshape(dimension, dimension)
        parent_depth = depths[head]
        for action_index, generator in enumerate(generators):
            bracket = (generator @ parent - parent @ generator) % prime
            add(bracket, ["ad", action_index, head], parent_depth + 1)
            if len(raw_rows) == target_dimension:
                break
        head += 1
    if len(raw_rows) != target_dimension:
        raise AssertionError(
            f"modular word span stopped at {len(raw_rows)}, below target {target_dimension}"
        )

    minor_determinant = math.prod(insertion_pivots) % prime
    assert minor_determinant != 0
    recipe_payload = json.dumps(recipes, separators=(",", ":")).encode("ascii")
    return {
        "claim_tag": "[COMPUTATION]",
        "prime": prime,
        "field": f"F_{prime}",
        "matrix_dimension": dimension,
        "coordinate_count": dimension * dimension,
        "coordinate_order": "row-major entries in the displayed exact-Q quotient/component basis",
        "rank_lower_bound_Q": len(raw_rows),
        "word_count": len(raw_rows),
        "maximum_right_normed_word_depth": max(depths),
        "generator_denominator_lcms": denominator_lcms,
        "all_generator_denominators_nonzero_mod_prime": True,
        "same_basis_reduction": True,
        "word_recipes": recipes,
        "word_recipe_sha256": hashlib.sha256(recipe_payload).hexdigest(),
        "raw_word_rows_mod_prime_sha256": digest.hexdigest(),
        "selected_minor_pivot_coordinates": pivot_coordinates,
        "selected_minor_insertion_pivots_mod_prime": insertion_pivots,
        "selected_minor_determinant_mod_prime": minor_determinant,
        "lower_bound_reason": (
            "the selected raw rows are literal rational Lie words; their square minor on the "
            "listed columns reduces to the nonzero displayed determinant, so its determinant "
            "is nonzero over Q"
        ),
    }


def exact_orthogonal_image(label: str) -> tuple[dict, list[dict]]:
    quotient_dimension, target_dimension, absolute_type = ORTHOGONAL_TARGETS[label]
    actions = E49.symmetry_sector(label_key(label))
    kernel, change, inverse, conjugated, quotient = audited_common_kernel_quotient(actions)
    assert quotient[0].rows == quotient_dimension
    quotient_common_kernel_dimension = len(sp.Matrix.vstack(*quotient).nullspace())
    assert quotient_common_kernel_dimension == 0

    forms = E49.invariant_symmetric_forms(quotient)
    assert len(forms) == 1
    form = forms[0]
    zero = sp.zeros(quotient_dimension)
    form_identities = [generator.T * form + form * generator for generator in quotient]
    assert all(identity == zero for identity in form_identities)
    determinant = sp.factor(form.det())
    assert determinant != 0

    change_mod = rational_matrix_mod(change, GOOD_PRIME)
    inverse_mod = rational_matrix_mod(inverse, GOOD_PRIME)
    identity_mod = np.eye(change.rows, dtype=np.int64)
    same_basis_inverse = (
        np.array_equal(change_mod @ inverse_mod % GOOD_PRIME, identity_mod)
        and np.array_equal(inverse_mod @ change_mod % GOOD_PRIME, identity_mod)
    )
    assert same_basis_inverse
    assert rank_mod(rational_matrix_mod(form, GOOD_PRIME), GOOD_PRIME) == quotient_dimension

    modular = lie_word_minor_certificate(quotient, target_dimension)
    ambient_dimension = quotient_dimension * (quotient_dimension - 1) // 2
    assert ambient_dimension == target_dimension == modular["rank_lower_bound_Q"]
    kernel_expected = 2 if label == "000" else 1
    assert kernel.cols == kernel_expected

    record = {
        "claim_tag": "[THEOREM]",
        "sector": label,
        "source_module_dimension_Q": actions[0].rows,
        "common_zero_kernel_dimension_Q": kernel.cols,
        "quotient_dimension_Q": quotient_dimension,
        "kernel_audit": {
            "stacked_generator_kernel_rank_Q": kernel.cols,
            "annihilated_by_A_exactly": actions[0] * kernel == sp.zeros(actions[0].rows, kernel.cols),
            "annihilated_by_B_exactly": actions[1] * kernel == sp.zeros(actions[0].rows, kernel.cols),
            "invariant_under_both_generators": True,
            "kernel_basis": rational_matrix_checksum(kernel),
        },
        "quotient_coordinate_audit": {
            "change_of_basis": rational_matrix_checksum(change),
            "change_inverse": rational_matrix_checksum(inverse),
            "conjugated_generators": [rational_matrix_checksum(matrix) for matrix in conjugated],
            "quotient_generators": [rational_matrix_checksum(matrix) for matrix in quotient],
            "agrees_exactly_with_e49_quotient": True,
            "change_and_inverse_reduce_to_inverses_mod_prime": same_basis_inverse,
            "same_exact_basis_reduced_mod_prime": True,
            "quotient_common_zero_kernel_dimension_Q": quotient_common_kernel_dimension,
        },
        "symmetric_form_Q": {
            "solution_space_dimension_Q": len(forms),
            "rank_Q": form.rank(),
            "determinant": str(determinant),
            "determinant_nonzero_mod_prime": int(determinant) % GOOD_PRIME != 0,
            "generator_containment_identities_zero": [identity == zero for identity in form_identities],
            "witness": rational_matrix_checksum(form, store_entries=True),
        },
        "image_Q": {
            "exact_dimension": target_dimension,
            "upper_bound": target_dimension,
            "upper_bound_reason": (
                f"exact generator identities put the image in so(F_{label}); for a nonsingular "
                f"{quotient_dimension}-dimensional symmetric form this ambient algebra has "
                f"dimension {quotient_dimension}*({quotient_dimension}-1)/2={target_dimension}"
            ),
            "lower_bound": target_dimension,
            "lower_certificate": modular,
            "equality": f"image = so(F_{label}, Q)",
            "absolute_Dynkin_type": absolute_type,
            "complexification": f"so({quotient_dimension}, C)",
            "simple_after_complexification": True,
            "split_Q_form_claimed": False,
            "Q_form_scope": "the full orthogonal algebra of the displayed rational form; only the absolute type is asserted",
        },
    }
    checks = [
        {"name": f"{label}_kernel_and_quotient_exact_Q", "passed": True,
         "detail": f"kernel {kernel.cols}, quotient {quotient_dimension}; exact e49 basis identity"},
        {"name": f"{label}_unique_nondegenerate_symmetric_form_Q", "passed": True,
         "detail": f"form nullity 1, rank {quotient_dimension}, both containment residuals zero"},
        {"name": f"{label}_{absolute_type}_same_basis_minor", "passed": True,
         "detail": f"good-prime nonzero {target_dimension}-minor plus exact orthogonal upper bound"},
    ]
    return record, checks


def component_actions_from_projectors(actions, projectors):
    inclusions = []
    blocks = []
    for _, projector in projectors:
        columns = projector.columnspace()
        inclusion = sp.Matrix.hstack(*columns)
        inclusions.append(inclusion)
        blocks.append(E49.restrict_to_projector(actions, projector))
    change = sp.Matrix.hstack(*inclusions)
    inverse = change.inv()
    conjugated = tuple(inverse * generator * change for generator in actions)
    expected = block_diagonal(blocks)
    assert all(left == right for left, right in zip(conjugated, expected))
    assert change * inverse == sp.eye(change.rows) == inverse * change
    return inclusions, blocks, change, inverse, conjugated


def invariant_bilinear_forms(actions: tuple[sp.Matrix, sp.Matrix]) -> list[sp.Matrix]:
    assert all(generator == sp.diag(*(generator[i, i] for i in range(generator.rows)))
               for generator in actions[1:2])
    target = tuple(-generator.T for generator in actions)
    forms = E49.intertwiner_basis(target, actions)
    zero = sp.zeros(actions[0].rows)
    assert all(
        all(generator.T * form + form * generator == zero for generator in actions)
        for form in forms
    )
    return forms


def exact_reductive_image(label: str) -> tuple[dict, list[dict]]:
    expected_dimensions, target_dimension, levi_types = REDUCTIVE_TARGETS[label]
    actions = E49.symmetry_sector(label_key(label))
    projectors = E49.multiplicity_free_projectors(*actions)
    dimensions = [dimension for dimension, _ in projectors]
    assert dimensions == expected_dimensions
    identity = sp.eye(actions[0].rows)
    assert sum((projector for _, projector in projectors), sp.zeros(actions[0].rows)) == identity
    for _, projector in projectors:
        assert projector * projector == projector
        assert all(projector * generator == generator * projector for generator in actions)
    for i, (_, left) in enumerate(projectors):
        for j, (_, right) in enumerate(projectors):
            assert left * right == (left if i == j else sp.zeros(actions[0].rows))

    inclusions, blocks, change, inverse, conjugated = component_actions_from_projectors(actions, projectors)
    assert all(block[1] == sp.diag(*(block[1][i, i] for i in range(block[1].rows))) for block in blocks)

    hom_dimensions = []
    for target in blocks:
        row = []
        for source in blocks:
            hom = E49.intertwiner_basis(target, source)
            zero = sp.zeros(target[0].rows, source[0].rows)
            assert all(
                target_generator * witness - witness * source_generator == zero
                for witness in hom
                for target_generator, source_generator in zip(target, source)
            )
            row.append(len(hom))
        hom_dimensions.append(row)
    assert hom_dimensions == [
        [1 if i == j else 0 for j in range(len(blocks))] for i in range(len(blocks))
    ]

    component_records = []
    for dimension, (projector_dimension, projector), inclusion, block in zip(
        dimensions, projectors, inclusions, blocks
    ):
        assert dimension == projector_dimension
        traceless = E49.traceless_actions(block)
        bilinear_forms = invariant_bilinear_forms(traceless)
        symmetric = sum(form.T == form for form in bilinear_forms)
        alternating = sum(form.T == -form for form in bilinear_forms)
        other = len(bilinear_forms) - symmetric - alternating
        assert not bilinear_forms
        component_records.append({
            "dimension_Q": dimension,
            "projector": rational_matrix_checksum(projector),
            "inclusion_basis": rational_matrix_checksum(inclusion),
            "generator_blocks": [rational_matrix_checksum(generator) for generator in block],
            "generator_traces_Q": [str(generator.trace()) for generator in block],
            "traceless_generator_invariant_bilinear_forms_Q": {
                "total_dimension": len(bilinear_forms),
                "symmetric_dimension": symmetric,
                "alternating_dimension": alternating,
                "other_dimension": other,
            },
        })

    block_generators = block_diagonal(blocks)
    trace_vectors = [sp.Matrix([block[which].trace() for block in blocks]) for which in range(2)]
    trace_rank = sp.Matrix.hstack(*trace_vectors).rank()
    assert trace_rank == 1 and trace_vectors[0] == sp.zeros(len(blocks), 1)
    upper_dimension = 1 + sum(dimension * dimension - 1 for dimension in dimensions)
    assert upper_dimension == target_dimension

    change_mod = rational_matrix_mod(change, GOOD_PRIME)
    inverse_mod = rational_matrix_mod(inverse, GOOD_PRIME)
    identity_mod = np.eye(change.rows, dtype=np.int64)
    assert np.array_equal(change_mod @ inverse_mod % GOOD_PRIME, identity_mod)
    assert np.array_equal(inverse_mod @ change_mod % GOOD_PRIME, identity_mod)
    modular = lie_word_minor_certificate(block_generators, target_dimension)
    assert modular["rank_lower_bound_Q"] == upper_dimension

    scalar_blocks = []
    for block in blocks:
        scalar_blocks.append(block[1].trace() * sp.eye(block[1].rows) / block[1].rows)
    central_trace_line = sp.diag(*scalar_blocks)
    assert all(central_trace_line * generator == generator * central_trace_line
               for generator in block_generators)
    assert all((block[0] - block[0].trace() * sp.eye(block[0].rows) / block[0].rows).trace() == 0
               for block in blocks)

    record = {
        "claim_tag": "[THEOREM]",
        "sector": label,
        "module_dimension_Q": actions[0].rows,
        "constituent_dimensions_Q": dimensions,
        "projector_audit": {
            "sum_to_identity": True,
            "pairwise_orthogonal_idempotents": True,
            "commute_with_both_generators": True,
            "rational_component_change_basis": rational_matrix_checksum(change),
            "rational_component_change_inverse": rational_matrix_checksum(inverse),
            "exact_block_diagonal_generator_identity": True,
            "same_basis_reduces_invertibly_mod_prime": True,
        },
        "pairwise_Hom_dimensions_Q": hom_dimensions,
        "components": component_records,
        "block_trace_space": {
            "generator_trace_vectors_Q": [[str(value) for value in vector] for vector in trace_vectors],
            "dimension_Q": trace_rank,
            "central_trace_line": rational_matrix_checksum(central_trace_line),
        },
        "image_Q": {
            "exact_dimension": target_dimension,
            "upper_bound": upper_dimension,
            "upper_bound_reason": (
                "exact rational projectors make the action block diagonal; every commutator has "
                "zero trace on every block, while the two generator trace vectors span one line"
            ),
            "lower_bound": target_dimension,
            "lower_certificate": modular,
            "equality": "Q*z plus direct sum of the full traceless block algebras",
            "centre_dimension_Q": 1,
            "Levi_types_Q": levi_types,
            "Levi_algebras_Q": [f"sl({dimension}, Q)" for dimension in dimensions],
            "reductive": True,
        },
    }
    checks = [
        {"name": f"{label}_exact_constituent_decomposition_Q", "passed": True,
         "detail": f"dimensions {dimensions}; rational orthogonal idempotents and exact block identity"},
        {"name": f"{label}_pairwise_Hom_and_form_audit_Q", "passed": True,
         "detail": f"Hom matrix is identity; all traceless constituent invariant-form spaces vanish"},
        {"name": f"{label}_exact_reductive_image_Q", "passed": True,
         "detail": f"same-basis nonzero {target_dimension}-minor saturates trace-space upper bound"},
    ]
    return record, checks


def build_artifact() -> dict:
    started = time.monotonic()
    orthogonal_images = {}
    reductive_images = {}
    assert sp.isprime(GOOD_PRIME)
    checks = [{
        "name": "good_prime_is_prime",
        "passed": True,
        "detail": f"{GOOD_PRIME} is prime; every represented denominator is checked separately",
    }]
    timings = {}

    for label in ORTHOGONAL_TARGETS:
        stage = time.monotonic()
        orthogonal_images[label], new_checks = exact_orthogonal_image(label)
        timings[f"orthogonal_{label}"] = round(time.monotonic() - stage, 6)
        checks.extend(new_checks)
    for label in REDUCTIVE_TARGETS:
        stage = time.monotonic()
        reductive_images[label], new_checks = exact_reductive_image(label)
        timings[f"reductive_{label}"] = round(time.monotonic() - stage, 6)
        checks.extend(new_checks)

    oa_artifact = json.loads(OA_RESULT.read_text(encoding="utf-8"))
    levi_corollary = oa_artifact["data"]["levi_corollary"]
    assert levi_corollary["claim_tag"] == "[THEOREM]"
    assert "direct sum of sl2 copies" in levi_corollary["statement"]
    assert orthogonal_images["000"]["image_Q"]["absolute_Dynkin_type"] == "D21"
    # Content digest over the mathematical payload only (sorted keys, no meta/checks):
    # the byte-level file digest is NOT stable, because the e47 producer refreshes
    # timestamps in `meta` on every rerun while `data` is unchanged.
    oa_content_sha256 = hashlib.sha256(
        json.dumps(oa_artifact["data"], sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    oa_no_go = {
        "claim_tag": "[THEOREM]",
        "layer": "single finite open 2x4 layer",
        "base_field": "C after exact-Q construction and scalar extension",
        "further_quotient": "so(42,C), simple type D21",
        "classification_input": levi_corollary["statement"],
        "classification_artifact": "results/algebra_structure/oa_quotient.json",
        "classification_artifact_content_sha256": oa_content_sha256,
        "conclusion": (
            "the complexified open-2x4 layer Lie algebra is not a finite-dimensional "
            "quotient of the Onsager algebra; no generating pair in it can satisfy the "
            "normalized Dolan-Grady presentation"
        ),
        "proof_chain": [
            "the 000 sector gives an exact characteristic-zero surjection onto a Q-form of D21",
            "scalar extension gives a surjection onto simple so(42,C)",
            "a quotient of an OA quotient would itself be an OA quotient",
            "the audited Date-Roan corollary permits only sl2^n as a semisimple OA quotient",
            "simple D21 is not sl2, giving a contradiction",
        ],
    }
    checks.append({
        "name": "open_2x4_Onsager_quotient_no_go",
        "passed": True,
        "detail": "exact D21 further quotient contradicts the audited OA semisimple-quotient corollary",
    })
    timings["total"] = round(time.monotonic() - started, 6)

    artifact = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "arithmetic": "exact Q via SymPy and exact F_p at one explicitly audited good prime; no floating point",
            "source_e49": "experiments/e49_char0_levi_2x4.py",
            "method": (
                "exact rational invariant quotients/projectors and forms; orthogonal or block-trace "
                "upper bounds; nonzero modular minors of literal rational Lie-word matrices in the same bases"
            ),
        },
        "data": {
            "scope": "single finite open 2x4 layer only; no all-sizes or 3D exact-solution claim",
            "good_prime": GOOD_PRIME,
            "orthogonal_quotient_images_Q": orthogonal_images,
            "other_sector_images_Q": reductive_images,
            "complete_2x4_Levi_sum": {
                "claim_tag": "[UNRESOLVED]",
                "value": None,
                "reason": "sector images overlap; exact quotient images are not summed into the full layer Levi decomposition",
            },
            "Onsager_quotient_no_go": oa_no_go,
            "timings_seconds": timings,
        },
        "checks": checks,
    }
    assert all(check["passed"] for check in checks)
    return artifact


def run(write: bool = True) -> dict:
    artifact = build_artifact()
    if write:
        RESULT.parent.mkdir(parents=True, exist_ok=True)
        RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return artifact


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-write", action="store_true", help="build and verify without replacing the JSON artifact")
    arguments = parser.parse_args(argv)
    artifact = run(write=not arguments.no_write)
    for check in artifact["checks"]:
        print(f"PASS: {check['name']} ({check['detail']})")
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
