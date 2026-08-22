"""Standalone regression checks for experiments/e32_series_extend2.py."""

from __future__ import annotations

import json
import hashlib
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SERIES = ROOT / "results" / "series"


def load(name: str) -> dict:
    payload = json.loads((SERIES / name).read_text(encoding="utf-8"))
    assert set(payload) == {"provenance", "data", "checks"}, name
    assert payload["provenance"]["script"] == "experiments/e32_series_extend2.py", name
    assert payload["provenance"]["timestamp"], name
    assert payload["provenance"]["precision"], name
    assert payload["checks"] and all(check["passed"] for check in payload["checks"]), name
    return payload


def coefficients(payload: dict) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in payload["data"]["coefficients"])


def formal_log(values: list[int]) -> tuple[Fraction, ...]:
    result = [Fraction(0) for _ in values]
    assert values[0] == 1
    for degree in range(1, len(values)):
        convolution = sum(
            (
                index * result[index] * values[degree - index]
                for index in range(1, degree)
            ),
            Fraction(0),
        )
        result[degree] = Fraction(values[degree]) - convolution / degree
    return tuple(result)


def main() -> None:
    ht = load("extended2_sc_ht_free_energy.json")
    lt = load("extended2_sc_lt_free_energy.json")
    analysis = load("dfinite_holdout.json")
    frozen_ht = json.loads((SERIES / "extended_sc_ht_free_energy.json").read_text(encoding="utf-8"))
    frozen_lt = json.loads((SERIES / "extended_sc_lt_free_energy.json").read_text(encoding="utf-8"))

    new_ht = coefficients(ht)
    new_lt = coefficients(lt)
    old_ht = coefficients(frozen_ht)
    old_lt = coefficients(frozen_lt)
    assert new_ht[: len(old_ht)] == old_ht
    assert new_lt[: len(old_lt)] == old_lt

    assert ht["data"]["achieved_order"] == 22
    assert lt["data"]["achieved_order"] == 32
    assert len(new_ht) == 23 and len(new_lt) == 33
    assert new_lt[30] and new_lt[32] and new_ht[22]
    assert ht["data"]["coefficient_provenance"] == "derived_locally"
    assert lt["data"]["coefficient_provenance"] == "derived_locally"

    profiles = {
        ("HT", 22): ht["data"]["resource_profile_v22"],
        ("HT", 24): ht["data"]["resource_wall_v24"],
        ("LT", 32): lt["data"]["resource_profile"],
    }
    assert profiles[("HT", 22)]["memory_wall_box"] == {
        "shape": [4, 5, 5],
        "sites": 100,
        "cross_section_sites": 20,
        "full_polynomial_degree": 235,
        "estimated_peak_bytes": 5939134464,
        "within_crt_cross_section_guard": True,
    }
    assert profiles[("HT", 24)]["memory_wall_box"]["shape"] == [5, 5, 5]
    assert profiles[("HT", 24)]["memory_wall_box"]["cross_section_sites"] == 25
    assert profiles[("HT", 24)]["memory_wall_box"]["estimated_peak_bytes"] == 242397216768
    assert not profiles[("HT", 24)]["reachable_by_current_guard"]
    assert profiles[("LT", 32)]["memory_wall_box"]["shape"] == [3, 3, 3]
    assert profiles[("LT", 32)]["memory_wall_box"]["estimated_peak_bytes"] == 1339392
    witness = lt["data"]["external_witness"]
    assert witness["manifest_key"] == "guttmann_enting1993"
    source = ROOT / witness["local_path"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == witness["sha256"]
    manifest = (ROOT / "sources" / "manifest.yaml").read_text(encoding="utf-8")
    assert "key: guttmann_enting1993" in manifest
    assert witness["sha256"] in manifest


    # Independently transform the exact integers printed in Guttmann--Enting Table 2.
    lambda_values = [0] * 17
    for degree, value in {
        0: 1, 3: 1, 5: 3, 6: -3, 7: 15, 8: -30, 9: 101,
        10: -261, 11: 807, 12: -2308, 13: 7065, 14: -21171,
        15: 65337, 16: -200934,
    }.items():
        lambda_values[degree] = value
    log_lambda = formal_log(lambda_values)
    assert new_lt[30] == log_lambda[15] == Fraction(346966, 5)
    assert new_lt[32] == log_lambda[16] == Fraction(-427509, 2)

    phi_values = [0] * 23
    for degree, value in {
        0: 1, 4: 3, 6: 22, 8: 192, 10: 2046, 12: 24853,
        14: 329334, 16: 4649601, 18: 68884356, 20: 1059830112,
        22: 16809862992,
    }.items():
        phi_values[degree] = value
    log_phi = formal_log(phi_values)
    assert new_ht[22] == log_phi[22] + Fraction(3, 22) == Fraction(360734896503, 22)

    dfinite = analysis["data"]["d_finiteness"]
    assert dfinite["known_coefficient_count"] == 11
    assert dfinite["largest_fully_nonvacuous_parameter_budget"] == 11
    assert dfinite["smallest_vacuous_parameter_budget"] == 12
    assert dfinite["consistent_relation_count"] == 0
    assert dfinite["conclusion"] == (
        "no D-finite recurrence of budget <= 11 is consistent with the known coefficients"
    )
    assert dfinite["tested_pairs"]
    for item in dfinite["tested_pairs"]:
        assert item["parameter_budget"] <= 11
        assert item["rank"] == item["parameter_budget"]
        assert item["nullity"] == 0

    strict = analysis["data"]["strict_holdout"]
    ode_rows = strict["d_finite_interpolants"]
    da_rows = strict["differential_approximants"]
    assert len(ode_rows) >= 5 and len(da_rows) >= 4
    assert all(row["held_out_coefficients"] >= 1 for row in ode_rows + da_rows)
    assert all(row["selected_pair"] is not None for row in ode_rows)
    assert sum(row["selected_degrees"] is not None for row in da_rows) >= 3
    assert all(row["selected_absolute_errors"] for row in ode_rows)
    assert sum(bool(row["selected_absolute_errors"]) for row in da_rows) >= 3

    critical = analysis["data"]["critical_estimate"]
    assert not critical["benchmark_used_in_fit_or_selection"]
    assert Fraction(critical["K_c_uncertainty"]) > 0
    assert not analysis["data"]["final_comparison_only"]["used_in_fit_or_model_selection"]

    note = (ROOT / "notes" / "series_extension2.md").read_text(encoding="utf-8")
    assert "no D-finite recurrence of budget <= 11 is consistent with the known coefficients" in note
    assert "cross-section 25" in note and "242,397,216,768" in note
    assert "[COMPUTATION]" in note and "[LEMMA]" in note and "[EXTERNAL]" in note
    assert note.rstrip().splitlines()[-1].startswith("[EXTERNAL] Final comparison only:")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_series_extend2: {type(error).__name__}: {error}")
        raise SystemExit(1)
    print("PASS test_series_extend2")
