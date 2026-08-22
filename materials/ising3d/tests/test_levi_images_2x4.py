#!/usr/bin/env python3
"""Standalone verifier for the exact open-2x4 quotient-image theorems."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "e59_levi_images_2x4.py"
E49_PATH = ROOT / "experiments" / "e49_char0_levi_2x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "levi_images_2x4.json"
OA_RESULT = ROOT / "results" / "algebra_structure" / "oa_quotient.json"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


e49 = load("e49_test_levi_images", E49_PATH)
e59 = load("e59_test_levi_images", EXPERIMENT)


def check(name, passed, detail=""):
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        raise AssertionError(name)


def rational_mod_independent(matrix, prime):
    output = np.empty(matrix.shape, dtype=np.int64)
    for index, entry in enumerate(matrix):
        entry = sp.Rational(entry)
        denominator = int(entry.q)
        if denominator % prime == 0:
            raise AssertionError("bad denominator")
        output.flat[index] = int(entry.p) * pow(denominator, -1, prime) % prime
    return output


def checksum_independent(matrix):
    values = [sp.Rational(value) for value in matrix]
    payload = ";".join(f"{int(value.p)}/{int(value.q)}" for value in values).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def rebuild_000_quotient():
    actions = e49.symmetry_sector((0, 0, 0))
    dimension = actions[0].rows
    kernel_columns = sp.Matrix.vstack(*actions).nullspace()
    kernel = sp.Matrix.hstack(*kernel_columns)
    assert len(kernel_columns) == 2
    assert all(generator * kernel == sp.zeros(dimension, 2) for generator in actions)

    columns = list(kernel_columns)
    rank = 2
    identity = sp.eye(dimension)
    for j in range(dimension):
        trial = sp.Matrix.hstack(*columns, identity[:, j])
        new_rank = trial.rank()
        if new_rank > rank:
            columns.append(identity[:, j])
            rank = new_rank
        if rank == dimension:
            break
    change = sp.Matrix.hstack(*columns)
    inverse = change.inv()
    transformed = tuple(inverse * generator * change for generator in actions)
    assert all(matrix[:, :2] == sp.zeros(dimension, 2) for matrix in transformed)
    quotient = tuple(matrix[2:, 2:] for matrix in transformed)
    return actions, kernel, change, inverse, quotient


def independent_word_minor(generators_q, target, prime):
    dimension = generators_q[0].rows
    generators = tuple(rational_mod_independent(generator, prime) for generator in generators_q)
    reduced = []
    pivot_rows = {}
    raw = []
    recipes = []
    depths = []
    pivot_columns = []
    insertion_pivots = []
    raw_digest = hashlib.sha256()

    def insert(candidate, recipe, depth):
        candidate = np.asarray(candidate, dtype=np.int64).reshape(-1) % prime
        remainder = candidate.copy()
        while True:
            indices = np.flatnonzero(remainder)
            if not indices.size:
                return False
            column = int(indices[0])
            old = pivot_rows.get(column)
            if old is None:
                pivot = int(remainder[column])
                pivot_rows[column] = len(reduced)
                reduced.append(remainder * pow(pivot, prime - 2, prime) % prime)
                raw.append(candidate)
                recipes.append(recipe)
                depths.append(depth)
                pivot_columns.append(column)
                insertion_pivots.append(pivot)
                raw_digest.update(np.asarray(candidate, dtype="<u4").tobytes())
                return True
            remainder = (remainder - int(remainder[column]) * reduced[old]) % prime

    assert insert(generators[0], ["seed", 0], 1)
    assert insert(generators[1], ["seed", 1], 1)
    cursor = 0
    while cursor < len(raw) and len(raw) < target:
        parent = raw[cursor].reshape(dimension, dimension)
        for which, generator in enumerate(generators):
            candidate = (generator @ parent - parent @ generator) % prime
            insert(candidate, ["ad", which, cursor], depths[cursor] + 1)
            if len(raw) == target:
                break
        cursor += 1
    determinant = math.prod(insertion_pivots) % prime
    return {
        "rank": len(raw),
        "recipes": recipes,
        "pivot_columns": pivot_columns,
        "insertion_pivots": insertion_pivots,
        "determinant": determinant,
        "maximum_depth": max(depths),
        "raw_digest": raw_digest.hexdigest(),
    }


def main():
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(artifact) == {"provenance", "data", "checks"}
        and artifact["provenance"]["script"] == "experiments/e59_levi_images_2x4.py"
        and artifact["provenance"]["interpreter"] == ".venv/bin/python"
        and all(set(item) == {"name", "passed", "detail"} for item in artifact["checks"]),
    )
    data = artifact["data"]
    check(
        "scope discipline",
        data["complete_2x4_Levi_sum"]["claim_tag"] == "[UNRESOLVED]"
        and data["complete_2x4_Levi_sum"]["value"] is None
        and "single finite open 2x4" in data["scope"],
    )

    actions, kernel, change, inverse, quotient = rebuild_000_quotient()
    record = data["orthogonal_quotient_images_Q"]["000"]
    coordinate_audit = record["quotient_coordinate_audit"]
    check(
        "independent exact-Q 000 quotient",
        actions[0].rows == 44
        and kernel.cols == 2
        and quotient[0].rows == 42
        and not sp.Matrix.vstack(*quotient).nullspace()
        and coordinate_audit["quotient_common_zero_kernel_dimension_Q"] == 0
        and checksum_independent(change) == coordinate_audit["change_of_basis"]["sha256_row_major_rationals"]
        and checksum_independent(inverse) == coordinate_audit["change_inverse"]["sha256_row_major_rationals"]
        and all(
            checksum_independent(generator) == witness["sha256_row_major_rationals"]
            for generator, witness in zip(quotient, coordinate_audit["quotient_generators"])
        ),
    )

    form_witness = record["symmetric_form_Q"]["witness"]
    form = sp.Matrix(
        form_witness["rows"], form_witness["columns"],
        [sp.Rational(value) for value in form_witness["entries_row_major"]],
    )
    exact_forms = e49.invariant_symmetric_forms(quotient)
    check(
        "independent exact orthogonal containment",
        len(exact_forms) == 1
        and form.rank() == 42
        and form.det() != 0
        and all(generator.T * form + form * generator == sp.zeros(42) for generator in quotient)
        and checksum_independent(form) == form_witness["sha256_row_major_rationals"],
    )

    prime = data["good_prime"]
    change_mod = rational_mod_independent(change, prime)
    inverse_mod = rational_mod_independent(inverse, prime)
    identity_mod = np.eye(44, dtype=np.int64)
    check(
        "same quotient basis has good reduction",
        np.array_equal(change_mod @ inverse_mod % prime, identity_mod)
        and np.array_equal(inverse_mod @ change_mod % prime, identity_mod)
        and int(form.det()) % prime != 0
        and all(e59.matrix_denominator_lcm(generator) % prime for generator in quotient),
    )

    rebuilt = independent_word_minor(quotient, 861, prime)
    stored_minor = record["image_Q"]["lower_certificate"]
    check(
        "independent D21 word minor",
        rebuilt["rank"] == 861
        and rebuilt["determinant"] != 0
        and rebuilt["determinant"] == stored_minor["selected_minor_determinant_mod_prime"]
        and rebuilt["pivot_columns"] == stored_minor["selected_minor_pivot_coordinates"]
        and rebuilt["insertion_pivots"] == stored_minor["selected_minor_insertion_pivots_mod_prime"]
        and rebuilt["recipes"] == stored_minor["word_recipes"]
        and rebuilt["maximum_depth"] == stored_minor["maximum_right_normed_word_depth"]
        and rebuilt["raw_digest"] == stored_minor["raw_word_rows_mod_prime_sha256"],
        f"det={rebuilt['determinant']} mod {prime}",
    )
    check(
        "D21 exact dimension squeeze",
        record["image_Q"]["upper_bound"] == 42 * 41 // 2 == 861
        and record["image_Q"]["lower_bound"] == 861
        and record["image_Q"]["exact_dimension"] == 861
        and record["image_Q"]["absolute_Dynkin_type"] == "D21"
        and record["image_Q"]["simple_after_complexification"],
    )

    expected_orthogonal = {"000": (861, "D21"), "010": (351, "B13"),
                           "100": (351, "B13"), "110": (351, "B13")}
    check(
        "all promoted orthogonal images",
        all(
            data["orthogonal_quotient_images_Q"][label]["image_Q"]["exact_dimension"] == dimension
            and data["orthogonal_quotient_images_Q"][label]["image_Q"]["absolute_Dynkin_type"] == kind
            and data["orthogonal_quotient_images_Q"][label]["image_Q"]["lower_certificate"]["selected_minor_determinant_mod_prime"] != 0
            for label, (dimension, kind) in expected_orthogonal.items()
        ),
    )
    expected_other = {
        "001": (606, ["A3", "A3", "A23"]),
        "111": (606, ["A3", "A3", "A23"]),
        "011": (799, ["A3", "A27"]),
        "101": (799, ["A3", "A27"]),
    }
    check(
        "exact 606 and 799 image identifications",
        all(
            data["other_sector_images_Q"][label]["image_Q"]["exact_dimension"] == dimension
            and data["other_sector_images_Q"][label]["image_Q"]["Levi_types_Q"] == types
            and data["other_sector_images_Q"][label]["image_Q"]["centre_dimension_Q"] == 1
            and data["other_sector_images_Q"][label]["pairwise_Hom_dimensions_Q"]
            == [[1 if i == j else 0 for j in range(len(types))] for i in range(len(types))]
            for label, (dimension, types) in expected_other.items()
        ),
    )

    oa_artifact = json.loads(OA_RESULT.read_text(encoding="utf-8"))
    oa_corollary = oa_artifact["data"]["levi_corollary"]
    no_go = data["Onsager_quotient_no_go"]
    # Content digest over the mathematical payload only: the e47 producer refreshes
    # `meta` timestamps on every rerun, so byte-level digests of the file are unstable.
    oa_digest = hashlib.sha256(
        json.dumps(oa_artifact["data"], sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    # A further quotient of an OA quotient is again an OA quotient.  The rebuilt
    # simple D21 quotient is semisimple and is not a direct sum of A1 factors,
    # so the audited Date--Roan corollary gives the stated contradiction.
    check(
        "Onsager-quotient implication",
        oa_corollary["claim_tag"] == "[THEOREM]"
        and "direct sum of sl2 copies" in oa_corollary["statement"]
        and no_go["claim_tag"] == "[THEOREM]"
        and no_go["classification_artifact_content_sha256"] == oa_digest
        and no_go["further_quotient"] == "so(42,C), simple type D21"
        and "not a finite-dimensional quotient" in no_go["conclusion"],
    )

    check("all generated checks", all(item["passed"] for item in artifact["checks"]))
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
