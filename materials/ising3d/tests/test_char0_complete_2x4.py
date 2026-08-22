#!/usr/bin/env python3
"""Standalone exact verifier for the complete 2x4 characteristic-zero theorem."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "e76_char0_complete_2x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "char0_complete_2x4.json"
SPEC = importlib.util.spec_from_file_location("e76_tested", SCRIPT)
e76 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = e76
SPEC.loader.exec_module(e76)


def check(name, passed, detail=""):
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        raise AssertionError(name)


def verify_sparse(certificate):
    matrix = e76.decode_sparse_matrix(certificate)
    observed = e76.sparse_matrix_certificate(matrix)
    return observed == certificate, matrix


def verify_raw_minor(certificate):
    prime = certificate["prime"]
    pivots = certificate["unnormalised_insertion_pivots"]
    determinant = 1
    for pivot in pivots:
        determinant = determinant * pivot % prime
    return (
        certificate["rank_Fp"] == len(certificate["raw_parent_recipes"])
        == len(certificate["pivot_columns"])
        == len(pivots)
        and len(set(certificate["pivot_columns"])) == certificate["rank_Fp"]
        and determinant == certificate["selected_minor_determinant_mod_p"]
        and determinant != 0
        and certificate["all_denominators_nonzero_mod_p"]
    )


def main():
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(artifact) == {"provenance", "data", "checks"}
        and artifact["provenance"]["script"]
        == "experiments/e76_char0_complete_2x4.py"
        and artifact["provenance"]["interpreter"] == ".venv/bin/python"
        and all(set(item) == {"name", "passed", "detail"} for item in artifact["checks"]),
    )
    data = artifact["data"]
    exact = data["exact_Q_certificate"]
    modular = data["modular_word_certificate"]

    # Independently replay the full and derived literal Pauli-word closures at
    # the stored good prime.  This rebuilds both nonzero modular minors from
    # their actual nested-word construction rather than trusting JSON ranks.
    replayed = e76.modular_word_certificate()
    for part, expected in (("full_closure", 2952), ("derived_ideal", 2951)):
        check(
            f"replayed {part} word minor",
            replayed[part]["rank_Fp"] == expected
            and replayed[part]["tree"] == modular[part]["tree"]
            and replayed[part]["grading_block_ranks"]
            == modular[part]["grading_block_ranks"],
        )

    # Recompute every stored raw-word selected-minor determinant, including all
    # five joint linkage certificates that control repeated factor identities.
    image_minors = exact["literal_word_image_minors_Fp"]
    linkage_minors = exact["literal_word_linkage_minors_Fp"]
    expected_images = {
        "000": 861, "001": 606, "010": 351, "011": 799,
        "100": 351, "101": 799, "110": 351, "111": 606,
    }
    expected_links = {
        "010_plus_100": 702,
        "100_plus_110": 351,
        "001_plus_111": 606,
        "011_plus_101": 799,
        "001_plus_011": 1389,
    }
    check(
        "all literal sector image minors",
        {key: value["rank_Fp"] for key, value in image_minors.items()}
        == expected_images
        and all(verify_raw_minor(value) for value in image_minors.values()),
    )
    check(
        "all literal linkage minors",
        {key: value["rank_Fp"] for key, value in linkage_minors.items()}
        == expected_links
        and all(verify_raw_minor(value) for value in linkage_minors.values()),
    )

    # Independently rebuild exact-Q sectors and verify the stored physical Gram
    # splittings.  Invariance makes common kernels direct zero summands, closing
    # the extension/cocycle gap in classical image upper bounds.
    sectors = {key: e76.sector_with_gram(key) for key in e76.CHARACTERS}
    for key, (a, b, gram) in sectors.items():
        label = "".join(map(str, key))
        record = exact["sector_certificates_Q"][label]
        kernel_ok, kernel = verify_sparse(record["kernel_basis"])
        complement_ok, complement = verify_sparse(record["gram_orthogonal_complement"])
        change = sp.Matrix.hstack(kernel, complement)
        inverse = change.inv()
        conjugated = (inverse * a * change, inverse * b * change)
        k = kernel.cols
        split = all(
            generator[:k, :] == sp.zeros(k, generator.cols)
            and generator[:, :k] == sp.zeros(generator.rows, k)
            for generator in conjugated
        )
        check(
            f"{label} exact Gram kernel split",
            kernel_ok and complement_ok and change.rank() == a.rows and split,
        )

    # Re-verify at least three exact invariant-form witnesses from stored sparse
    # entries; all four are checked here, including nondegeneracy and equations.
    for label, expected_rank in {"000": 42, "010": 27, "100": 27, "110": 27}.items():
        record = exact["sector_certificates_Q"][label]
        form_ok, form = verify_sparse(record["invariant_symmetric_form"])
        a, b, gram = sectors[tuple(map(int, label))]
        _, _, active = e76.gram_split_common_kernel((a, b), gram)
        check(
            f"{label} exact invariant form witness",
            form_ok
            and form.rank() == expected_rank
            and form == form.T
            and all(generator.T * form + form * generator == sp.zeros(expected_rank) for generator in active),
        )

    # Rebuild exact intertwiner equations and compare dimensions/ranks.  These
    # checks use exact rational arithmetic, not stored checksums alone.
    rebuilt = e76.exact_q_certificate()
    stored_intertwiners = exact["intertwiners_Q"]
    rebuilt_intertwiners = rebuilt["intertwiners_Q"]
    for family in stored_intertwiners:
        stored_items = stored_intertwiners[family]
        rebuilt_items = rebuilt_intertwiners[family]
        if isinstance(stored_items, dict):
            stored_items = [stored_items]
            rebuilt_items = [rebuilt_items]
        check(
            f"{family} exact intertwiners",
            [(item["space_dimension_Q"], item["rank_Q"], item["dual"]) for item in stored_items]
            == [(item["space_dimension_Q"], item["rank_Q"], item["dual"]) for item in rebuilt_items]
            and all(verify_sparse(item["witness"])[0] for item in stored_items),
        )

    central_ok, central_integer = verify_sparse(
        exact["central_witness_global_integer"]
    )
    central_scale = sp.Rational(
        exact["central_projection_equals_integer_times"]
    )
    central = central_scale * central_integer
    full_a = sp.diag(*(sectors[key][0] for key in e76.CHARACTERS))
    full_b = sp.diag(*(sectors[key][1] for key in e76.CHARACTERS))
    check("central sparse global validity", central_ok)
    check(
        "central independent rebuild equality",
        exact["central_witness_global_integer"]
        == rebuilt["central_witness_global_integer"]
        and exact["central_projection_equals_integer_times"]
        == rebuilt["central_projection_equals_integer_times"],
    )
    check("central witness nonzero", central != sp.zeros(256))
    check("central commutes with A", central * full_a == full_a * central)
    check("central commutes with B", central * full_b == full_b * central)
    check(
        "factor map kernel dimension one",
        data["factor_map"]["kernel_dimension_Q"] == 1
        and data["factor_map"]["kernel_equals_center"] is True,
    )

    check(
        "dimension and Levi squeeze",
        data["dimension_Q"]["value"] == 2952
        and data["derived_algebra"]["dimension_Q"] == 2951
        and sum(data["Levi_type"]["factor_dimensions"]) == 2951
        and sum(data["Levi_type"]["factor_ranks"]) == 103
        and data["Levi_type"]["type"]
        == "D21 + B13 + B13 + A23 + A27 + A3 + A3",
    )
    check(
        "radical and structural Killing conclusions",
        data["solvable_radical"]["equals_center"] is True
        and data["Killing"]["rank_Q"] == 2951
        and data["Killing"]["nullity_Q"] == 1
        and data["Killing"]["direct_adjoint_matrix"]["materialized"] is False
        and data["Killing"]["direct_adjoint_matrix"]["claim_tag"] == "[UNRESOLVED]",
    )
    check("all artifact checks", all(item["passed"] for item in artifact["checks"]))
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
