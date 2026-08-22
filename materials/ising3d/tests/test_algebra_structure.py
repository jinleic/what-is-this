"""Standalone regression tests for e29 algebra classification."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "e29_algebra_structure.py"
# Quick-mode producer artifacts (owned by this test's subprocess run):
RESULT = ROOT / "results" / "algebra_structure" / "structure_quick.json"
CONTROL = ROOT / "results" / "algebra_structure" / "onsager_controls_quick.json"
# Canonical full artifacts (READ-ONLY here; only the non-quick experiment writes them):
CANONICAL_RESULT = ROOT / "results" / "algebra_structure" / "structure.json"
CANONICAL_CONTROL = ROOT / "results" / "algebra_structure" / "onsager_controls.json"


def fail(message: str) -> int:
    print(f"FAIL: {message}")
    return 1


def main() -> int:
    # Exercise the producer, including the exact-Q 4-cycle path and one-prime
    # sparse 2x4 classification.  Quick mode writes only *_quick.json; the canonical
    # structure.json retains the full second-prime certificate and must not be touched.
    canonical_before = {
        path: path.read_bytes() for path in (CANONICAL_RESULT, CANONICAL_CONTROL)
    }
    run = subprocess.run(
        [sys.executable, str(SCRIPT), "--quick"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=3600,  # calibrated for a contended machine: e29 --quick took 26.6 min wall
        # (8.4 min CPU) alongside external SAT/qec jobs on 2026-08-12; 900 s was idle-machine only
        check=False,
    )
    if run.returncode != 0 or run.stdout.strip().splitlines()[-1:] != ["PASS"]:
        return fail(f"experiment failed: stdout={run.stdout!r} stderr={run.stderr!r}")

    for path, before in canonical_before.items():
        if path.read_bytes() != before:
            return fail(f"quick run mutated the canonical artifact {path.name}")

    # Canonical full artifact must retain the two-prime 2x4 certificate.
    canonical = json.loads(CANONICAL_RESULT.read_text(encoding="utf-8"))
    full24 = canonical["data"]["grid_2x4"]
    if full24["classification_prime_2"] is None or full24["prime_2_invariants"] is None:
        return fail("canonical structure.json lost the second-prime 2x4 certificate")
    if (full24["classification_prime_1"]["sector_image_dimensions"]
            != full24["classification_prime_2"]["sector_image_dimensions"]):
        return fail("canonical 2x4 two-prime sector dimensions disagree")

    structure = json.loads(RESULT.read_text(encoding="utf-8"))
    control = json.loads(CONTROL.read_text(encoding="utf-8"))
    for name, artifact in (("structure", structure), ("controls", control)):
        if set(artifact) != {"provenance", "data", "checks"}:
            return fail(f"{name} envelope keys are wrong")
        if not artifact["checks"] or not all(check["passed"] for check in artifact["checks"]):
            return fail(f"{name} contains a failed check")
        if artifact["provenance"]["script"] != "experiments/e29_algebra_structure.py":
            return fail(f"{name} provenance script is wrong")

    four = structure["data"]["four_cycle_exact_Q"]
    expected_four = {
        "dimension": 11,
        "killing_rank": 9,
        "killing_radical_dimension": 2,
        "centre_dimension": 2,
        "solvable_radical_dimension": 2,
        "levi_dimension": 9,
        "levi_type": "A1^3",
        "derived_series_dimensions": [11, 9, 9],
        "lower_central_dimensions": [11, 9, 9],
    }
    for key, value in expected_four.items():
        if four[key] != value:
            return fail(f"four-cycle {key}: got {four[key]!r}, expected {value!r}")
    if four["root_system"]["simple_root_cartan_matrix"] != [[2]]:
        return fail("four-cycle root certificate is not A1")
    if four["root_system"]["simple_factors"][0]["multiplicity"] != 3:
        return fail("four-cycle does not exhibit three A1 factors")
    if len(four["killing_form"]) != 11 or any(len(row) != 11 for row in four["killing_form"]):
        return fail("four-cycle Killing matrix has wrong shape")

    grid23 = structure["data"]["grid_2x3"]
    inv23 = [grid23["prime_1_invariants"], grid23["prime_2_invariants"]]
    if [(x["dimension"], x["centre_dimension"], x["derived_dimension"]) for x in inv23] != [(263, 1, 262), (263, 1, 262)]:
        return fail("2x3 computed modular invariants changed")
    if (grid23["solvable_radical"] != "NOT_COMPUTED_BY_THIS_EXPERIMENT"
            or grid23["levi_decomposition"] != "NOT_COMPUTED_BY_THIS_EXPERIMENT"):
        return fail("2x3 e29-scope fields changed")
    # Cross-artifact SSOT: e29 points at e45's artifact, which owns the exact
    # characteristic-zero resolution.  Validate the pointer and the resolution.
    pointer23 = grid23["characteristic_zero_resolution"]
    if (not pointer23["resolved_elsewhere"]
            or pointer23["artifact"] != "results/algebra_structure/char0_levi.json"):
        return fail("2x3 characteristic-zero resolution pointer is wrong")
    char0 = json.loads((ROOT / pointer23["artifact"]).read_text(encoding="utf-8"))
    if not char0["checks"] or not all(c["passed"] for c in char0["checks"]):
        return fail("char0_levi.json contains a failed check")
    c0 = char0["data"]
    if not (c0["dimension_Q"] == 263
            and c0["Killing"]["rank_Q"] == 262
            and c0["solvable_radical"]["dimension_Q"] == 1
            and c0["solvable_radical"]["equals_center"]
            and [f["dimension"] for f in c0["semisimple_quotient"]["factors"]] == [105, 21, 21, 80, 35]):
        return fail("char0_levi.json headline invariants changed")
    inv23_dims = (inv23[0]["dimension"], inv23[0]["centre_dimension"], inv23[0]["derived_dimension"])
    if inv23_dims != (c0["dimension_Q"], c0["solvable_radical"]["dimension_Q"], c0["derived_dimension_Q"]):
        return fail("modular invariants disagree with the characteristic-zero certification")
    if grid23["classification_prime_1"]["levi_factors"] is not None:
        return fail("2x3 quotient images were promoted to Levi factors")
    certified23 = grid23["classification_prime_1"]["certified_quotient_images"]
    if [(x["field_type"], x["dimension"]) for x in certified23] != [("C7", 105), ("C3", 21), ("gl6", 36)]:
        return fail("2x3 F_p quotient-image certificates changed")
    form23 = grid23["classification_prime_1"]["linkage_witnesses"]["invariant_form_witnesses"]
    if form23["000"] != {"kernel_dimension": 0, "symmetric": [0, None], "alternating": [1, 14]}:
        return fail("2x3 C7 alternating-form certificate changed")
    if grid23["classification_prime_1"]["uncertified_quotient_images"][0]["dimension"] != 81:
        return fail("2x3 uncertified 81-dimensional image was not retained")

    grid24 = structure["data"]["grid_2x4"]
    inv24 = [grid24["prime_1_invariants"]]
    if grid24["prime_2_invariants"] is not None:
        inv24.append(grid24["prime_2_invariants"])
    if any((x["dimension"], x["centre_dimension"], x["derived_dimension"]) != (2952, 1, 2951) for x in inv24):
        return fail("2x4 computed modular invariants changed")
    if (grid24["solvable_radical"] != "NOT_COMPUTED_BY_THIS_EXPERIMENT"
            or grid24["levi_decomposition"] != "NOT_COMPUTED_BY_THIS_EXPERIMENT"):
        return fail("2x4 e29-scope fields changed")
    if grid24["characteristic_zero_resolution"]["resolved_elsewhere"] is not False:
        return fail("2x4 characteristic-zero resolution pointer must remain unresolved until a certifying artifact exists")
    if grid24["classification_prime_1"]["levi_factors"] is not None:
        return fail("2x4 quotient images were incorrectly promoted to Levi factors")
    witnesses = grid24["classification_prime_1"]["linkage_witnesses"]
    if witnesses["D21_form_witness"] != {"symmetric": [1, 42], "alternating": [0, None]}:
        return fail("2x4 D21 F_p image certificate changed")
    if witnesses["B13_form_witness"] != {"symmetric": [1, 27], "alternating": [0, None]}:
        return fail("2x4 B13 F_p image certificate changed")
    if "205.8 GB" not in grid24["resource_wall"]:
        return fail("2x4 resource wall was not recorded")
    for grid in (grid23, grid24):
        if "solvable radical" not in grid["field_provenance"]["not_computed"]:
            return fail("not-computed provenance guard missing")

    p1_cases = control["data"]["prime_1_cases"]
    p2_cases = control["data"]["prime_2_cases"]
    if [case["dimension"] for case in p1_cases] != [case["dimension"] for case in p2_cases]:
        return fail("Onsager controls disagree between primes")
    for case in p1_cases:
        n = int(case["graph"].split("_")[-1])
        if case["graph"].startswith("open"):
            expected = (n * n, 1, n * n - 1)
        else:
            expected = (3 * n - 1, 2, 3 * n - 3)
        actual = (case["dimension"], case["centre_dimension"], case["derived_dimension"])
        if actual != expected:
            return fail(f"Onsager control {case['graph']}: got {actual}, expected {expected}")
    if "OA / ker(rho_open,n)" not in control["data"]["open_chain_exact_quotient"]:
        return fail("open-chain exact quotient was not named")
    if "OA / ker(rho_per,n)" not in control["data"]["ring_exact_quotient"]:
        return fail("ring exact quotient was not named")

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
