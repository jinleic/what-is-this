"""Zero-level CCZ provenance and common-unit Gate-B accounting sensitivity.

pre_statement.md Revision 4 (2026-09-04). Deterministic arithmetic, no
sampling, no campaign. The frozen campaign
`campaigns/20260829T223540Z_dcf693ab/` remains untouched.

The frozen Gate-B resource ratio divided 22 physical qubits * 24 circuit
layers by Litinski physical qubits * syndrome-extraction cycles. That
mixed-unit value is retained here only as a historical reproduction row.

The current comparison uses one unit on both sides:
physical qubits * surface-code syndrome-extraction cycles per accepted output.
The zero-level paper reports that depth 24 corresponds to three syndrome
rounds. Its per-accepted-state cost is therefore

    (22 + three output-patch footprints) * 3 / acceptance.

For the two paper-simulated output variants, a full rotated distance-d patch is
priced as 2d^2-1 physical qubits for d in {3, 7}. This static full-footprint
pricing is a derived sensitivity, not a paper resource estimate.

The QLOPS Table-6 `cycles` values already include Litinski postselection, so
the accepted-output comparator is

    7 * unit_qubits * reported_cycles_including_postselection.

All resource products, the seven-T-to-CCZ conversion, patch pricing, ratios,
and verdicts are [DERIVED]. Reported inputs remain explicitly identified.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

TARGET = Path(__file__).resolve().parents[1]
SRC = TARGET / "src"
sys.path.insert(0, str(SRC))

import formulas as F  # noqa: E402
import paper_data as P  # noqa: E402

ARTIFACT = SRC / "evidence" / "zero_level_sensitivity.json"
PATCH_DISTANCES = (3, 7)          # the two output variants the paper simulates
FALSIFIER = P.GATE_B["comparability_falsifier"]


def rotated_patch_qubits(d: int) -> int:
    """Rotated surface code, distance d: d^2 data + (d^2 - 1) syndrome ancillas."""
    return 2 * d * d - 1


def frozen_gate_b_legacy_row() -> dict[str, object]:
    """Reproduce the frozen mixed-unit Gate-B magnitude without endorsing it."""
    z = P.ZERO_LEVEL_CCZ
    baselines = []
    for row in P.T6:
        m = re.match(r"15to1_(\d+),(\d+),(\d+)", row["proto"])
        dX, dZ, dm = map(int, m.groups())
        unit = F.litinski_15to1_unit_qubits(dX, dZ, dm)
        cycles = F.litinski_15to1_plain_cycles(dm)
        baselines.append(7 * unit * cycles)
    zero_cost = float(z["legacy_gate_b_spacetime_mixed_units"])
    ratios = [zero_cost / baseline for baseline in baselines]
    return {
        "basis": (
            "LEGACY MIXED UNITS: zero-level physical-qubits*circuit-layers "
            "divided by Litinski physical-qubits*syndrome-cycles"),
        "zero_level_spacetime_per_attempt_mixed_units": zero_cost,
        "acceptance_assumption": 1.0,
        "ratio_zero_vs_litinski_ccz": [min(ratios), max(ratios)],
        "advantage_factor": [1.0 / max(ratios), 1.0 / min(ratios)],
        "usable_for_common_basis_magnitude": False,
        "tag": "[DERIVED; LEGACY-MIXED-UNITS]",
    }


def litinski_accepted_ccz_spacetime() -> list[dict[str, object]]:
    """Accepted 7-T CCZ-equivalent cost from reported QLOPS Table-6 rows."""
    rows = []
    for row in P.T6:
        per_t = row["unit"] * row["cycles"]
        rows.append({
            "target": f"{row['code']}/{row['proto']}",
            "unit_qubits": row["unit"],
            "cycles_including_postselection": row["cycles"],
            "spacetime_per_accepted_T": per_t,
            "spacetime_per_accepted_CCZ_equiv": 7 * per_t,
            "tag": "[DERIVED from REPORTED Table-6 unit and cycles]",
        })
    return rows




def zero_level_variants() -> dict[str, dict[str, object]]:
    """Common-time-unit zero-level costs per attempt for the sourced variants."""
    z = P.ZERO_LEVEL_CCZ
    rounds = z["syndrome_rounds"]
    out: dict[str, dict[str, object]] = {
        "distillation_only_output_patches_omitted": {
            "basis": "22 distillation qubits; three output patches omitted",
            "patch_distance": None,
            "qubits": z["phys_qubits"],
            "syndrome_rounds": rounds,
            "spacetime_per_attempt": z["phys_qubits"] * rounds,
            "scope": "diagnostic only; omits required output-patch footprint",
            "tag": "[DERIVED from REPORTED qubits and syndrome rounds]",
        }
    }
    for d in PATCH_DISTANCES:
        patch = rotated_patch_qubits(d)
        qubits = z["phys_qubits"] + z["log_qubits"] * patch
        out[f"with_output_patches_d{d}"] = {
            "basis": (
                f"22 distillation qubits + {z['log_qubits']} full rotated "
                f"d={d} patches ({patch} qubits each), held for three rounds"),
            "patch_distance": d,
            "patch_qubits_each": patch,
            "qubits": qubits,
            "syndrome_rounds": rounds,
            "spacetime_per_attempt": qubits * rounds,
            "scope": (
                "static full-footprint sensitivity; not a time-resolved "
                "paper resource estimate"),
            "tag": "[DERIVED from REPORTED qubits and syndrome rounds]",
        }
    return out


def sensitivity() -> dict[str, object]:
    z = P.ZERO_LEVEL_CCZ
    baseline = litinski_accepted_ccz_spacetime()
    variants = zero_level_variants()
    acceptances = {
        "p=1e-3_optimistic": {
            "probability": z["acceptance_p1e3"][1],
            "physical_error": 1e-3,
            "tag": "[REPORTED approximate range endpoint]",
            "scope": (
                "range endpoint applied to both paper variants as a "
                "sensitivity arm; the paper does not map 0.40 to one variant"),
        },
        "p=1e-3_pessimistic": {
            "probability": z["acceptance_p1e3"][0],
            "physical_error": 1e-3,
            "tag": "[REPORTED approximate range endpoint]",
            "scope": (
                "range endpoint applied to both paper variants as a "
                "sensitivity arm; the paper does not map 0.30 to one variant"),
        },
        "p=1e-4": {
            "probability": z["acceptance_p1e4"],
            "physical_error": 1e-4,
            "tag": "[REPORTED approximate value]",
            "scope": "approximate paper value applied to both variants",
        },
    }
    combos: dict[str, dict[str, object]] = {}
    for variant_name, variant in variants.items():
        for acceptance_name, acceptance in acceptances.items():
            probability = float(acceptance["probability"])
            zero_cost = float(variant["spacetime_per_attempt"]) / probability
            endpoint_rows = [
                (row["spacetime_per_accepted_CCZ_equiv"] / zero_cost,
                 row["target"])
                for row in baseline
            ]
            lower = min(endpoint_rows)
            upper = max(endpoint_rows)
            ratios = [
                zero_cost / row["spacetime_per_accepted_CCZ_equiv"]
                for row in baseline
            ]
            combos[f"{variant_name}|accept_{acceptance_name}"] = {
                "zero_level_spacetime_per_accepted_CCZ": zero_cost,
                "acceptance": probability,
                "ratio_zero_vs_litinski_ccz": [min(ratios), max(ratios)],
                "advantage_factor": [lower[0], upper[0]],
                "lower_endpoint_target": lower[1],
                "upper_endpoint_target": upper[1],
                "still_exceeds_falsifier": bool(lower[0] > FALSIFIER),
                "tag": "[DERIVED]",
            }

    worst_name, worst_combo = min(
        combos.items(), key=lambda item: item[1]["advantage_factor"][0])
    worst = {
        "combination": worst_name,
        "advantage_factor": worst_combo["advantage_factor"],
        "lower_endpoint_target": worst_combo["lower_endpoint_target"],
        "upper_endpoint_target": worst_combo["upper_endpoint_target"],
        "tag": "[DERIVED]",
    }
    return {
        "revision": "Revision 4 (qlops)",
        "kind": "deterministic re-derivation; no sampling or campaign",
        "paper": {
            "arxiv": z["arxiv"],
            "full_text_read_utc": "2026-09-03",
            "constants_located_in_body": {
                "c=300": (
                    "Sec. IV least-squares fit + Fig. 10 caption "
                    "(also abstract, Sec. I, Sec. V)"),
                "phys_qubits=22": (
                    "abstract, Sec. I, Sec. III, Sec. V "
                    "(distillation circuit only)"),
                "depth=24_and_syndrome_rounds=3": (
                    "Sec. III and Sec. V"),
                "log_qubits=3": (
                    "abstract, Sec. I, Sec. III, Sec. V (output patches)"),
            },
            "fit_p_range": list(z["fit_p_range"]),
            "fit_trials_per_p": z["fit_trials_per_p"],
            "fit_coefficient_uncertainty": z["fit_coefficient_uncertainty"],
            "fit_residuals": z["fit_residuals"],
            "pointwise_error_bars": z["pointwise_error_bars"],
            "decoder_named": z["decoder"],
            "acceptance_p1e3": list(z["acceptance_p1e3"]),
            "acceptance_p1e4": z["acceptance_p1e4"],
            "machine_readable_circuit": False,
            "paper_own_spacetime_metric": (
                "logical qubits x surface-code syndrome rounds "
                "(Sec. IV, Fig. 13)"),
        },
        "conflicts_recorded": [
            "The frozen 22*24 resource row mixes circuit layers with "
            "Litinski syndrome-extraction cycles and has no common-unit "
            "magnitude interpretation.",
            "22 excludes the three output surface-code patches, which the "
            "Litinski comparator's footprint includes.",
            "300 p^2 is acceptance-conditional: acceptance is approximately "
            "30-40% at p=1e-3 and 90% at p=1e-4.",
            "QLOPS Table-6 cycles already include Litinski postselection; "
            "the accepted-output comparison must use the reported cycles.",
        ],
        "comparison_basis": {
            "metric": (
                "physical qubits x surface-code syndrome-extraction cycles "
                "per accepted CCZ"),
            "zero_level_time": {
                "syndrome_rounds": z["syndrome_rounds"],
                "tag": "[REPORTED]",
            },
            "litinski_time": (
                "QLOPS Table-6 reported cycles, including postselection"),
            "ccz_conversion_T_states": 7,
            "ccz_conversion_tag": "[DERIVED assumption]",
            "scope": (
                "static full-patch footprint using the reported equivalence "
                "depth 24 = 3 syndrome rounds; not an end-to-end paper estimate"),
        },
        "falsifier": FALSIFIER,
        "legacy_frozen_gate_b_row": frozen_gate_b_legacy_row(),
        "litinski_accepted_baseline_rows": baseline,
        "zero_level_variants": variants,
        "acceptance_scenarios": acceptances,
        "corrected_combinations": combos,
        "worst_case": worst,
        "verdict_survives_corrections": bool(
            worst["advantage_factor"][0] > FALSIFIER),
    }


def main() -> None:
    payload = sensitivity()
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {ARTIFACT}")
    print("sha256", hashlib.sha256(ARTIFACT.read_bytes()).hexdigest())
    legacy = payload["legacy_frozen_gate_b_row"]["advantage_factor"]
    print(f"legacy mixed-unit frozen row: {legacy[0]:.0f}x .. {legacy[1]:.0f}x")
    for name, combination in payload["corrected_combinations"].items():
        advantage = combination["advantage_factor"]
        print(
            f"  {name:67s} {advantage[0]:8.1f}x .. "
            f"{advantage[1]:10.1f}x")
    worst = payload["worst_case"]
    print(
        f"common-basis worst case {worst['combination']}: "
        f"{worst['advantage_factor'][0]:.1f}x .. "
        f"{worst['advantage_factor'][1]:.1f}x; verdict survives: "
        f"{payload['verdict_survives_corrections']}")


if __name__ == "__main__":
    main()
