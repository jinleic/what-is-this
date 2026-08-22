#!/usr/bin/env python3
"""Standalone exact verifier for the partial 2x4 characteristic-zero result."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "e49_char0_levi_2x4.py"
RESULT = ROOT / "results" / "algebra_structure" / "char0_levi_2x4.json"
SPEC = importlib.util.spec_from_file_location("e49", SCRIPT)
e49 = importlib.util.module_from_spec(SPEC)
sys.modules["e49"] = e49
assert SPEC.loader is not None
SPEC.loader.exec_module(e49)


def check(name, passed, detail=""):
    print(("PASS" if passed else "FAIL") + f": {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        raise AssertionError(name)


def matrix_from_witness(witness):
    return sp.Matrix(witness["rows"], witness["columns"], witness["entries_row_major"])


def main():
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    check(
        "artifact envelope",
        set(artifact) == {"provenance", "data", "checks"}
        and artifact["provenance"]["script"] == "experiments/e49_char0_levi_2x4.py"
        and artifact["provenance"]["interpreter"] == ".venv/bin/python"
        and all(set(item) == {"name", "passed", "detail"} for item in artifact["checks"]),
    )
    data = artifact["data"]
    check(
        "scope and full theorem",
        data["scope"].startswith("single finite 2x4")
        and data["Killing"]["rank_Q"] == 2951
        and data["Killing"]["nullity_Q"] == 1
        and data["solvable_radical"]["equals_center"] is True
        and data["Levi_type"]["complete_direct_sum"]
        == "D21 + B13 + B13 + A23 + A27 + A3 + A3"
        and data["status"].startswith("[THEOREM]"),
    )

    sectors = {key: e49.symmetry_sector(key) for key in e49.CHARACTERS}
    expected_dimensions = [44, 32, 28, 32, 28, 32, 28, 32]
    check(
        "independent exact-Q sector modules",
        [sectors[key][0].rows for key in e49.CHARACTERS] == expected_dimensions,
    )

    # Independently re-solve four quotient common kernels and verify at least
    # three stored invariant-form witnesses against the exact generator action.
    for label in ("000", "010", "100", "110"):
        key = tuple(map(int, label))
        actions = sectors[key]
        kernel, quotient = e49.quotient_by_common_kernel(actions)
        record = data["sector_certificates_Q"][label]
        witness = matrix_from_witness(record["invariant_symmetric_form_witness"])
        expected_kernel = 2 if label == "000" else 1
        expected_quotient = 42 if label == "000" else 27
        check(f"{label} exact common kernel", kernel.cols == expected_kernel)
        check(
            f"{label} stored invariant form",
            witness.rows == expected_quotient
            and witness.rank() == expected_quotient
            and all(
                generator.T * witness + witness * generator == sp.zeros(expected_quotient)
                for generator in quotient
            )
            and e49.matrix_checksum(witness, store_entries=True)
            == record["invariant_symmetric_form_witness"],
        )

    # Independently recompute the rational spectral projectors behind the four
    # unidentified sector images.
    expected_components = {
        "001": [4, 4, 24],
        "011": [4, 28],
        "101": [4, 28],
        "111": [4, 4, 24],
    }
    for label, expected in expected_components.items():
        a, b = sectors[tuple(map(int, label))]
        projectors = e49.multiplicity_free_projectors(a, b)
        check(
            f"{label} constituent projectors",
            [dimension for dimension, _ in projectors] == expected
            and sum((projector for _, projector in projectors), sp.zeros(a.rows)) == sp.eye(a.rows)
            and all(projector * projector == projector for _, projector in projectors),
        )

    # Independently reconstruct exact-Q linkage witnesses and recheck stored
    # full-rank intertwiner matrices.
    linkage = e49.exact_linkage_certificates(sectors)
    stored_linkage = data["linkage_certificates_Q"]
    check(
        "exact factor linkages",
        linkage["B13_linkage"]["100_and_110_same_factor"]["full_rank"] == 28
        and linkage["B13_linkage"]["010_and_100_independent"]["prior_modular_joint_rank"] == 702
        and [item["full_rank"] for item in linkage["odd_dual_pairings"]]
        == [4, 4, 24, 4, 28]
        and linkage["shared_A3_between_001_and_011"]["full_rank"] == 4
        and linkage["dimension_upper_bounds_Q"]["all_eight_sectors"] == 2952
        and linkage["dimension_upper_bounds_Q"]["derived_all_eight_sectors"] == 2951
        and linkage == stored_linkage,
    )
    check(
        "Levi dimensions and rank",
        sum(data["Levi_type"]["factor_dimensions"]) == 2951
        and data["Levi_type"]["rank"] == 103,
    )

    check("all generated artifact checks", all(item["passed"] for item in artifact["checks"]))
    print("PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
