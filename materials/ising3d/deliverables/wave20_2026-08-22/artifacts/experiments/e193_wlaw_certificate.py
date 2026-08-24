#!/usr/bin/env python3
"""Assemble the exact W-law structural route-closure certificate.

The candidate and its obstruction are constructed before the seven historical
ladder rows are loaded.  Those rows are used only as held-out controls.  This
script never launches a ladder closure and never computes W_10.
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import e191_wlaw_structure as structure
import e192_wlaw_module as module

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e193_wlaw_certificate.py"
OUT_PATH = ROOT / "results" / "algebra_growth" / "wlaw.json"
REGRESSION_PATH = ROOT / "results" / "algebra_growth" / "ladder_w8_regression.json"
W9_PATH = ROOT / "results" / "ladder" / "w9_saturation.json"


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def mirror_sector_values(L: int, low: dict[int, int]) -> dict[int, int]:
    full: dict[int, int] = {}
    for sector, value in low.items():
        full[sector] = value
        full[2 * L - sector] = value
    return dict(sorted(full.items()))


def control_record(
    L: int,
    cyclic_dim: int,
    K_dim: int,
    W_dim: int,
    low_ranks: dict[int, int],
    deficits: dict[int, int],
    source: str,
) -> dict:
    ranks = mirror_sector_values(L, low_ranks)
    full_deficits = {sector: deficits.get(sector, 0) for sector in range(0, 2 * L + 1, 2)}
    K_sectors = {sector: ranks[sector] + full_deficits[sector] for sector in ranks}
    predicted = {
        2 * m: structure.structure_record(L, m)["candidate_dimension"]
        for m in range(L + 1)
    }
    sector_matches = {sector: ranks[sector] == predicted[sector] for sector in ranks}
    integrity = (
        sum(ranks.values()) == cyclic_dim
        and sum(full_deficits.values()) == W_dim
        and sum(K_sectors.values()) == K_dim
        and cyclic_dim + W_dim == K_dim
    )
    return {
        "L": L,
        "source": source,
        "role": "[COMPUTATION] held-out exact control; no entry is an input to the structural candidate",
        "cyclic_dimension": cyclic_dim,
        "annihilator_dimension": W_dim,
        "K_dimension": K_dim,
        "certified_rank_by_sector": {str(k): v for k, v in ranks.items()},
        "certified_deficit_by_sector": {str(k): v for k, v in full_deficits.items()},
        "certified_K_dimension_by_sector": {str(k): v for k, v in K_sectors.items()},
        "parameter_free_candidate_rank_by_sector": {str(k): v for k, v in predicted.items()},
        "sector_matches": {str(k): v for k, v in sector_matches.items()},
        "all_sector_matches": all(sector_matches.values()),
        "integrity_passed": integrity,
        "tag": "[COMPUTATION] prior exact-Q totals and sectors compared after construction",
    }


def load_controls() -> tuple[list[dict], dict]:
    regression = json.loads(REGRESSION_PATH.read_text())
    w9 = json.loads(W9_PATH.read_text())
    controls = []
    for row in regression["data"]["regression_closures"]:
        controls.append(control_record(
            L=int(row["L"]),
            cyclic_dim=int(row["cyclic_dim"]),
            K_dim=int(row["K_dim"]),
            W_dim=int(row["W_dim"]),
            low_ranks={int(k): int(v) for k, v in row["low_ranks"].items()},
            deficits={int(k): int(v) for k, v in row["W_sectors"].items()},
            source="results/algebra_growth/ladder_w8_regression.json",
        ))
    w9_data = w9["data"]
    controls.append(control_record(
        L=9,
        cyclic_dim=int(w9_data["W9"]["cyclic_Q_dim"]),
        K_dim=int(w9_data["K9_dimension_record"]["K_dim_formula"]),
        W_dim=int(w9_data["W9"]["value"]),
        low_ranks={int(k): int(v) for k, v in w9_data["closures"]["q1"]["low_ranks"].items()},
        deficits={int(k): int(v) for k, v in w9_data["W9"]["W_sectors"].items()},
        source="results/ladder/w9_saturation.json",
    ))
    source_status = {
        "regression_checks_passed": all(item["passed"] for item in regression["data"]["checks"]),
        "w9_checks_passed": all(item["passed"] for item in w9["checks"]),
        "w9_two_prime_low_ranks_equal": (
            w9_data["closures"]["q1"]["low_ranks"] == w9_data["closures"]["p1"]["low_ranks"]
        ),
    }
    return controls, source_status


def main() -> int:
    started = time.process_time()

    # Parameter-free construction.  Historical ranks are deliberately loaded
    # only after these records have been built.
    structure_records = structure.build_structure_records()
    module_records = module.build_module_records()
    controls, source_status = load_controls()

    by_Lm = {(row["L"], row["m"]): row for row in structure_records}
    control_by_L = {row["L"]: row for row in controls}
    l3m1 = by_Lm[(3, 1)]
    earliest = {
        "L": 3,
        "m": 1,
        "particle_sector": 2,
        "candidate_dimension": l3m1["candidate_dimension"],
        "natural_reflection_equivariant_rank_ceiling": l3m1["balanced_polarization_rank"],
        "forced_kernel_dimension": l3m1["balanced_polarization_kernel"],
        "certified_U_sector_rank": control_by_L[3]["certified_rank_by_sector"]["2"],
        "conclusion": "[THEOREM] the natural reflection-equivariant balanced map is already too small at the earliest controlled row; this is not a counterexample to the rank law",
    }

    checks: list[dict] = []

    def add_check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    formula_rows_ok = all(row["formula_matches_explicit"] for row in structure_records)
    add_check(
        "C1_structure_formulas_exact_L2_9",
        formula_rows_ok,
        f"explicit subset/pair enumeration agreed with every formula on {len(structure_records)} (L,m) rows",
    )
    delta_rows = [
        (row["L"], row["m"])
        for row in structure_records
        if row["top_wedge"]["symmetric_invariant_functional_dim"] == 1
    ]
    expected_delta_rows = [
        (L, L // 2) for L in structure.SMALL_L if L % 4 == 0
    ]
    add_check(
        "C2_top_wedge_delta_location",
        delta_rows == expected_delta_rows,
        f"nonzero symmetric top-wedge trace occurs exactly at {delta_rows} in the exact small range",
    )
    interior = [row for row in structure_records if 0 < row["m"] < row["L"]]
    interior_kernels_ok = all(row["balanced_polarization_kernel"] > 0 for row in interior)
    interior_ceiling_strict = all(
        row["balanced_polarization_rank"] < row["candidate_dimension"] for row in interior
    )
    add_check(
        "C3_reflection_kernel_and_strict_ceiling",
        interior_kernels_ok and interior_ceiling_strict,
        f"all {len(interior)} interior rows have a nonzero reflection-odd kernel and equivariant "
        "ceiling strictly below the candidate dimension",
    )
    earliest_ok = (
        earliest["candidate_dimension"] == earliest["certified_U_sector_rank"] == 6
        and earliest["natural_reflection_equivariant_rank_ceiling"] == 4
        and earliest["forced_kernel_dimension"] == 2
    )
    add_check(
        "C4_earliest_exact_countercertificate",
        earliest_ok,
        "L=3,m=1: candidate/certified rank 6, but every natural reflection-equivariant realization has rank at most 4",
    )
    for piece in ("D", "F", "D_plus"):
        rows = [row for row in module_records if row["piece"] == piece]
        add_check(
            f"C5_module_leakage_{piece}",
            all(row["all_checks_passed"] and row["unbalanced_leakage_support"] > 0 for row in rows),
            f"exact positive-coefficient leakage found for {piece} at every L=2..9",
        )
    exceptional_row = next(row for row in module_records if row["L"] == 4 and row["piece"] == "D_plus")
    add_check(
        "C6_exceptional_witness_is_traceless",
        exceptional_row["exceptional_top_wedge_trace"] == 0,
        "the L=4 D+ source lies in the top-wedge traceless hyperplane and still leaks",
    )
    control_integrity = all(row["integrity_passed"] for row in controls)
    control_matches = all(row["all_sector_matches"] for row in controls)
    add_check(
        "C7_old_controls_integrity",
        len(controls) == 7 and control_integrity and all(source_status.values()),
        "seven prior rows L=3..9 have consistent totals/sectors and all cited source checks pass",
    )
    add_check(
        "C8_old_controls_parameter_free_match",
        control_matches,
        f"all {sum(len(row['sector_matches']) for row in controls)} certified sector ranks match the independently constructed candidate dimensions",
    )

    rss = peak_rss_bytes()
    add_check(
        "C9_resource_ceiling",
        rss < 2 * 1024**3,
        f"producer peak RSS {rss} bytes is below 2 GiB",
    )

    data = {
        "headline": {
            "status": "[THEOREM] rigorous route closure; [UNRESOLVED] the symmetric-square rank law itself",
            "claim": "[THEOREM] the delta has a unique top-wedge trace model, but every natural reflection-equivariant balanced realization has a forced kernel and its image is not invariant under D, F, or D+",
            "law_status": "[UNRESOLVED] finite controls L=3..9 do not prove an all-L rank formula",
        },
        "general_formulas": {
            "candidate": "[LEMMA] dim Sym^2(exterior^m Q^L) - delta = C(L,m)(C(L,m)+1)/2 - [4|L and 2m=L]",
            "delta": "[THEOREM] the unique SL_L-invariant top-wedge bilinear form is symmetric exactly when L=2m with m even",
            "reflection_fixed_subsets": "[LEMMA] f(L,m)=[x^m](1+x)^(L mod 2)(1+x^2)^floor(L/2)",
            "reflection_kernel": "[THEOREM] ker(balanced polarization) is exactly the reflection-odd part, of dimension (C(L,m)^2-f(L,m)^2)/4; it is nonzero for 0<m<L",
            "reflection_rank": "[THEOREM] rank(balanced polarization) = (C(L,m)^2+2C(L,m)+f(L,m)^2)/4 exactly (= the balanced configuration orbit count); the same number bounds every diagonal sign gauge equivariant map",
        },
        "structure_records_exact_L2_9": structure_records,
        "module_records_exact_L2_9": module_records,
        "earliest_countercertificate": earliest,
        "old_exact_controls_L3_9": controls,
        "control_protocol": {
            "tag": "[COMPUTATION] holdout-only comparison",
            "candidate_constructed_before_artifact_load": True,
            "fit_parameters": 0,
            "fit_totals": 0,
            "fit_sector_values": 0,
            "control_totals": 7,
            "control_sector_values": sum(len(row["sector_matches"]) for row in controls),
        },
        "scope": [
            "[THEOREM] the balanced polarization/Reynolds route is closed, including its exceptional traceless restriction",
            "[THEOREM] explicit all-L witnesses in proofs/wlaw.md show non-invariance under each of D, F, D+ for L>=2",
            "[COMPUTATION] exact enumerations L=2..9 check every dimension and leakage lemma without fitting",
            "[UNRESOLVED] no isomorphism between the traceless symmetric-square candidate and U_L^(2m) is proved",
            "[UNRESOLVED] no rank, W_L value, or sector value is claimed for L>=10; W_10 was not computed",
        ],
    }
    envelope = {
        "meta": {
            "script": SCRIPT,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python": sys.version.split()[0],
            "artifact_schema": "meta/data/checks-v1",
            "exact_arithmetic": "Python int; explicit finite set/orbit and sparse-vector enumeration",
            "cpu_seconds": round(time.process_time() - started, 6),
            "peak_rss_bytes": rss,
            "source_artifacts": {
                str(REGRESSION_PATH.relative_to(ROOT)): file_sha256(REGRESSION_PATH),
                str(W9_PATH.relative_to(ROOT)): file_sha256(W9_PATH),
            },
        },
        "data": data,
        "checks": checks,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n")
    passed = all(item["passed"] for item in checks)
    print(
        f"[e193] wrote {OUT_PATH.relative_to(ROOT)}: "
        f"{sum(item['passed'] for item in checks)}/{len(checks)} checks, "
        f"cpu={envelope['meta']['cpu_seconds']}s rss={rss}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
