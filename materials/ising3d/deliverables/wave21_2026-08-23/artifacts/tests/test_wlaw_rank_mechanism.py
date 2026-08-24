#!/usr/bin/env python3
"""Clean-room verifier for the wave-21 W-law rank-mechanism certificate.

No producer module is imported.  Configuration orbits are rebuilt in a
(top-mask, bottom-mask) encoding, the annihilator restriction is ranked with a
different pivot convention, and the rung basis change is recomputed from exact
local matrices.  Coverage is contractual: omitting any (L,m) row for L=3..8,
either exceptional row, or any of the 16 local horizontal sources is a failure.

Run only when explicitly requested:
    .venv/bin/python tests/test_wlaw_rank_mechanism.py
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "wlaw_rank_mechanism.json"
SOURCE = ROOT / "results" / "ladder" / "l8_saturation.json"
L_VALUES = tuple(range(3, 9))
LABELS = ("0", "b", "a", "d")
WEIGHT = {"0": 0, "b": 1, "a": 1, "d": 2}
EXPECTED_CHECKS = {
    "C1_source_exact_annihilators_L3_8",
    "C2_good_primes",
    "C3_complete_sector_coverage_L3_8",
    "C4_direct_K_orbit_census",
    "C5_omitted_pairing_full_column_rank_mod_999983",
    "C5_omitted_pairing_full_column_rank_mod_1000003",
    "C6_two_slice_restriction_injective_exact_Q",
    "C7_exceptional_sector_4_2",
    "C8_exceptional_sector_8_4",
    "C9_bonding_action_complete_and_integral",
    "C10_bonding_action_block_tridiagonal",
    "C11_bonding_truncation_first_possible_row",
    "C12_bonding_witness_exact_symmetry_and_sector",
    "C13_bonding_witness_in_U_by_exact_duality",
    "C14_bonding_truncation_exact_kernel",
    "C15_resource_ceiling",
}

CHECK_COUNT = 0
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    global CHECK_COUNT
    CHECK_COUNT += 1
    status = "ok" if passed else "FAIL"
    print(f"[test_wlaw_rank_mechanism] {status:4} {name} {detail}")
    if not passed:
        FAILURES.append(name)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def reverse_mask(mask: int, L: int) -> int:
    image = 0
    for rung in range(L):
        image |= ((mask >> rung) & 1) << (L - 1 - rung)
    return image


def interleaved(top: int, bottom: int, L: int) -> int:
    value = 0
    for rung in range(L):
        value |= ((top >> rung) & 1) << (2 * rung)
        value |= ((bottom >> rung) & 1) << (2 * rung + 1)
    return value


def decode_interleaved(value: int, L: int) -> tuple[int, int]:
    top = 0
    bottom = 0
    for rung in range(L):
        top |= ((value >> (2 * rung)) & 1) << rung
        bottom |= ((value >> (2 * rung + 1)) & 1) << rung
    return top, bottom


def block_orbit(top: int, bottom: int, L: int) -> set[tuple[int, int]]:
    reflected = (reverse_mask(top, L), reverse_mask(bottom, L))
    return {(top, bottom), (bottom, top), reflected, reflected[::-1]}


def block_orbit_records(L: int) -> list[dict[str, int]]:
    seen: set[tuple[int, int]] = set()
    records: list[dict[str, int]] = []
    for top in range(1 << L):
        for bottom in range(1 << L):
            if (top.bit_count() + bottom.bit_count()) % 2 or (top, bottom) in seen:
                continue
            orbit = block_orbit(top, bottom, L)
            seen.update(orbit)
            representatives = [interleaved(t, b, L) for t, b in orbit]
            records.append(
                {
                    "representative": min(representatives),
                    "orbit_size": len(orbit),
                    "particle_sector": top.bit_count() + bottom.bit_count(),
                    "absolute_leg_imbalance": abs(top.bit_count() - bottom.bit_count()),
                }
            )
    records.sort(key=lambda row: row["representative"])
    return records


def source_records(source: dict[str, object]) -> dict[int, dict[str, object]]:
    data = source["data"]
    records = {int(row["L"]): row for row in data["regression_records"]}
    records[7] = data["L7_record"]
    records[8] = data["L8_record"]
    return records


def independently_exact_source(record: dict[str, object]) -> bool:
    validation = record["validation"]
    return (
        len(record["basis"]) == int(record["W_dim"])
        and sum(int(value) for value in record["W_sectors"].values())
        == int(record["W_dim"])
        and int(record["cyclic_Q_dim"]) + int(record["W_dim"])
        == int(record["K_dim"])
        and validation["A_invariant_over_Q"] is True
        and validation["basis_independent_over_Q"] is True
        and validation["four_B_invariant_over_Q"] is True
        and validation["sector_homogeneous"] is True
        and validation["vacuum_orthogonal"] is True
        and validation["outside_image_indices"] == []
    )


def modular_rank_smallest_pivot(rows: list[dict[int, int]], prime: int) -> int:
    """Sparse exact F_p rank, deliberately opposite to the producer's pivot order."""
    pivots: dict[int, dict[int, int]] = {}
    for original in rows:
        vector = {
            int(column): int(coefficient) % prime
            for column, coefficient in original.items()
            if int(coefficient) % prime
        }
        while vector:
            pivot = min(vector)
            coefficient = vector[pivot]
            existing = pivots.get(pivot)
            if existing is None:
                inverse = pow(coefficient, prime - 2, prime)
                normalized: dict[int, int] = {}
                for column, entry in vector.items():
                    value = entry * inverse % prime
                    if value:
                        normalized[column] = value
                pivots[pivot] = normalized
                break
            for column, entry in existing.items():
                value = (vector.get(column, 0) - coefficient * entry) % prime
                if value:
                    vector[column] = value
                else:
                    vector.pop(column, None)
    return len(pivots)


def restricted_W_rows(
    record: dict[str, object],
    high_orbits: list[dict[str, int]],
    sector: int,
) -> list[dict[int, int]]:
    index = {row["representative"]: position for position, row in enumerate(high_orbits)}
    sizes = {row["representative"]: row["orbit_size"] for row in high_orbits}
    output: list[dict[int, int]] = []
    for vector in record["basis"]:
        if int(vector[0]["representative"]).bit_count() != sector:
            continue
        row: dict[int, int] = {}
        for term in vector:
            representative = int(term["representative"])
            if representative in index:
                row[index[representative]] = (
                    int(term["coefficient"]) * sizes[representative]
                )
        output.append(row)
    return output


def expected_sector_record(
    L: int,
    m: int,
    orbits: list[dict[str, int]],
    source_record: dict[str, object],
    primes: tuple[int, ...],
) -> dict[str, object]:
    sector = 2 * m
    sector_orbits = [row for row in orbits if row["particle_sector"] == sector]
    high = [row for row in sector_orbits if row["absolute_leg_imbalance"] >= 4]
    low = [row for row in sector_orbits if row["absolute_leg_imbalance"] <= 2]
    rows = restricted_W_rows(source_record, high, sector)
    W_by_sector = {int(key): int(value) for key, value in source_record["W_sectors"].items()}
    W_dimension = W_by_sector.get(sector, 0)
    K_dimension = len(sector_orbits)
    U_dimension = K_dimension - W_dimension
    ranks = {prime: modular_rank_smallest_pivot(rows, prime) for prime in primes}
    n = comb(L, m)
    comparison = n * (n + 1) // 2 - int(L % 4 == 0 and 2 * m == L)
    return {
        "K": K_dimension,
        "low": len(low),
        "high": len(high),
        "W": W_dimension,
        "W_rows": len(rows),
        "U": U_dimension,
        "comparison": comparison,
        "ranks": ranks,
    }


NEW_VECTORS: dict[str, tuple[Fraction, ...]] = {
    "0": (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    "b": (Fraction(0), Fraction(1), Fraction(1), Fraction(0)),
    "a": (Fraction(0), Fraction(1), Fraction(-1), Fraction(0)),
    "d": (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
}
NEW_DUALS: dict[str, tuple[Fraction, ...]] = {
    "0": (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    "b": (Fraction(0), Fraction(1, 2), Fraction(1, 2), Fraction(0)),
    "a": (Fraction(0), Fraction(1, 2), Fraction(-1, 2), Fraction(0)),
    "d": (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
}


def independent_horizontal_action(left: str, right: str) -> dict[tuple[str, str], int]:
    old_source = [
        NEW_VECTORS[left][x] * NEW_VECTORS[right][y]
        for x in range(4)
        for y in range(4)
    ]
    old_target = [Fraction(0) for _ in range(16)]
    for x in range(4):
        for y in range(4):
            coefficient = old_source[4 * x + y]
            if not coefficient:
                continue
            old_target[4 * (x ^ 1) + (y ^ 1)] += coefficient
            old_target[4 * (x ^ 2) + (y ^ 2)] += coefficient
    answer: dict[tuple[str, str], int] = {}
    for new_left in LABELS:
        for new_right in LABELS:
            coefficient = Fraction(0)
            for x in range(4):
                for y in range(4):
                    coefficient += (
                        NEW_DUALS[new_left][x]
                        * NEW_DUALS[new_right][y]
                        * old_target[4 * x + y]
                    )
            coefficient *= 2
            if coefficient:
                if coefficient.denominator != 1:
                    raise AssertionError((left, right, coefficient))
                answer[(new_left, new_right)] = int(coefficient)
    return answer


def independent_rung_action(label: str) -> dict[str, int]:
    old_target = [Fraction(0) for _ in range(4)]
    for old, coefficient in enumerate(NEW_VECTORS[label]):
        old_target[old ^ 3] += coefficient
    answer: dict[str, int] = {}
    for new in LABELS:
        coefficient = 2 * sum(
            NEW_DUALS[new][old] * old_target[old] for old in range(4)
        )
        if coefficient:
            if coefficient.denominator != 1:
                raise AssertionError((label, coefficient))
            answer[new] = int(coefficient)
    return answer


def normalized_stored_horizontal(record: dict[str, object]) -> dict[tuple[str, str], int]:
    return {
        tuple(output["target"]): int(output["coefficient"])
        for output in record["outputs"]
    }


def normalized_stored_rung(record: dict[str, object]) -> dict[str, int]:
    return {
        str(output["target"]): int(output["coefficient"])
        for output in record["outputs"]
    }


def transformed_census(L: int, sector: int) -> tuple[int, int]:
    seen: set[tuple[str, ...]] = set()
    omitted = 0
    for word in product(LABELS, repeat=L):
        if sum(WEIGHT[label] for label in word) != sector:
            continue
        if sum(label == "a" for label in word) % 2:
            continue
        canonical = min(word, tuple(reversed(word)))
        if canonical in seen:
            continue
        seen.add(canonical)
        if sum(label == "a" for label in word) >= 4:
            omitted += 1
    return len(seen), omitted


def witness_block() -> dict[tuple[int, int], int]:
    vector = {(0, 0): 1}
    for rung in range(4):
        next_vector: dict[tuple[int, int], int] = {}
        for (top, bottom), coefficient in vector.items():
            next_vector[(top | (1 << rung), bottom)] = coefficient
            next_vector[(top, bottom | (1 << rung))] = -coefficient
        vector = next_vector
    return vector


def pairing_from_block_witness(
    witness: dict[tuple[int, int], int], basis_vector: list[dict[str, int]]
) -> int:
    answer = 0
    for term in basis_vector:
        representative = int(term["representative"])
        top, bottom = decode_interleaved(representative, 4)
        orbit = block_orbit(top, bottom, 4)
        answer += int(term["coefficient"]) * sum(
            witness.get(config, 0) for config in orbit
        )
    return answer


def stage_schema_and_coverage(envelope: dict[str, object]) -> None:
    check("V1_top_level_schema", set(envelope) == {"meta", "data", "checks"})
    checks = envelope["checks"]
    names = [str(row["name"]) for row in checks]
    check(
        "V1_stored_checks_unique_complete_passed",
        set(names) == EXPECTED_CHECKS
        and len(names) == len(EXPECTED_CHECKS)
        and len(names) == len(set(names))
        and all(row["passed"] is True for row in checks),
    )
    contract = envelope["data"]["coverage_contract"]
    required = {(L, m) for L in L_VALUES for m in range(L + 1)}
    stored_rows = {
        (int(row["L"]), int(row["m"]))
        for row in envelope["data"]["two_slice_observability_exact_L3_8"]
    }
    check(
        "V1_all_finite_rows_covered",
        set(map(tuple, contract["required_sector_pairs"])) == required
        and stored_rows == required
        and contract["required_L"] == list(L_VALUES),
        f"{len(stored_rows)}/{len(required)} sectors",
    )
    check(
        "V1_universal_local_coverage_and_scope",
        contract["required_horizontal_local_sources"] == 16
        and contract["required_rung_local_sources"] == 4
        and contract["required_exceptional_pairs"] == [[4, 2], [8, 4]]
        and contract["all_L_computation_claimed"] is False,
    )


def stage_source_and_observability(
    envelope: dict[str, object], source: dict[str, object]
) -> dict[int, dict[str, object]]:
    source_hash = envelope["meta"]["source_artifacts"][
        "results/ladder/l8_saturation.json"
    ]
    check("V2_source_sha256", source_hash == sha256(SOURCE))
    records = source_records(source)
    check(
        "V2_source_exact_Q_contract",
        all(row["passed"] is True for row in source["checks"])
        and set(records) == set(L_VALUES)
        and all(independently_exact_source(records[L]) for L in L_VALUES),
    )
    primes = tuple(int(value) for value in envelope["meta"]["primes"])
    stored = {
        (int(row["L"]), int(row["m"])): row
        for row in envelope["data"]["two_slice_observability_exact_L3_8"]
    }
    all_rows_ok = True
    full_rank_ok = True
    first_bad = ""
    direct_K: dict[int, int] = {}
    for L in L_VALUES:
        orbits = block_orbit_records(L)
        direct_K[L] = len(orbits)
        for m in range(L + 1):
            expected = expected_sector_record(L, m, orbits, records[L], primes)
            row = stored[(L, m)]
            modular = row["weighted_annihilator_restriction"]
            row_ok = (
                row["particle_sector"] == 2 * m
                and row["K_sector_dimension"] == expected["K"]
                and row["retained_abs_leg_imbalance_0_or_2_dimension"] == expected["low"]
                and row["omitted_abs_leg_imbalance_at_least_4_dimension"] == expected["high"]
                and row["annihilator_W_sector_dimension_exact_Q"] == expected["W"]
                and row["annihilator_basis_vectors_in_sector"] == expected["W_rows"]
                and row["cyclic_U_sector_dimension_exact_Q"] == expected["U"]
                and row["prior_symmetric_square_value_for_comparison_only"] == expected["comparison"]
                and row["prior_value_matches_exact_U_rank"] is True
                and row["retained_codomain_excess_over_U"] == expected["low"] - expected["U"]
            )
            ranks_ok = all(
                expected["ranks"][prime] == expected["high"]
                and modular[str(prime)]["rank_Fp_lower_bound_for_rank_Q"] == expected["high"]
                and modular[str(prime)]["rational_rank_ceiling_omitted_coordinate_dimension"]
                == expected["high"]
                and modular[str(prime)]["lower_equals_ceiling"] is True
                for prime in primes
            )
            exact_ok = (
                row["annihilator_restriction_rank_exact_Q"] == expected["high"]
                and row["restriction_kernel_on_U_dimension_exact_Q"] == 0
                and row["retained_restriction_rank_exact_Q"] == expected["U"]
                and row["all_checks_passed"] is True
            )
            if not (row_ok and ranks_ok and exact_ok) and not first_bad:
                first_bad = f"L={L},m={m}"
            all_rows_ok = all_rows_ok and row_ok and exact_ok
            full_rank_ok = full_rank_ok and ranks_ok
    check(
        "V2_direct_K_orbit_census",
        all(direct_K[L] == int(records[L]["K_dim"]) for L in L_VALUES),
        str(direct_K),
    )
    check("V2_all_sector_fields_recomputed", all_rows_ok, first_bad)
    check(
        "V2_modular_lower_meets_rational_ceiling_everywhere",
        full_rank_ok,
        "no modular rank is used without the omitted-column ceiling",
    )
    return records


def stage_exceptional(envelope: dict[str, object]) -> None:
    exceptional = envelope["data"]["exceptional_rows"]
    row4 = exceptional["L4_m2"]
    row8 = exceptional["L8_m4"]
    check(
        "V3_exceptional_4_2",
        (row4["L"], row4["m"]) == (4, 2)
        and row4["K_sector_dimension"] == 22
        and row4["omitted_abs_leg_imbalance_at_least_4_dimension"] == 1
        and row4["annihilator_restriction_rank_exact_Q"] == 1
        and row4["retained_restriction_rank_exact_Q"] == 20,
    )
    check(
        "V3_exceptional_8_4",
        (row8["L"], row8["m"]) == (8, 4)
        and row8["K_sector_dimension"] == 3270
        and row8["omitted_abs_leg_imbalance_at_least_4_dimension"] == 433
        and row8["annihilator_restriction_rank_exact_Q"] == 433
        and row8["retained_restriction_rank_exact_Q"] == 2484,
    )


def stage_bonding(
    envelope: dict[str, object], records: dict[int, dict[str, object]]
) -> None:
    block = envelope["data"]["bonding_antibonding"]
    horizontal = {tuple(row["source"]): row for row in block["horizontal_actions"]}
    rung = {str(row["source"]): row for row in block["rung_actions"]}
    expected_horizontal_sources = set(product(LABELS, repeat=2))
    check(
        "V4_complete_local_source_table",
        set(horizontal) == expected_horizontal_sources and set(rung) == set(LABELS),
    )
    horizontal_ok = all(
        normalized_stored_horizontal(horizontal[source])
        == independent_horizontal_action(*source)
        for source in expected_horizontal_sources
    )
    rung_ok = all(
        normalized_stored_rung(rung[label]) == independent_rung_action(label)
        for label in LABELS
    )
    check("V4_exact_local_basis_change", horizontal_ok and rung_ok)
    deltas_ok = all(
        output["delta_antibonding_count"] in {-2, 0, 2}
        and output["delta_particle_number"] in {-2, 0, 2}
        for row in list(horizontal.values()) + list(rung.values())
        for output in row["outputs"]
    )
    check(
        "V4_all_local_cases_block_tridiagonal",
        deltas_ok and horizontal[("b", "a")]["outputs"] == []
        and horizontal[("a", "b")]["outputs"] == [],
    )

    failure = block["first_truncation_failure"]
    early_specs = ((3, 0), (3, 2), (3, 4), (3, 6), (4, 0), (4, 2), (4, 4))
    early = {
        f"L{L}_k{sector}": {
            "tau_plus_rho_orbit_dimension": transformed_census(L, sector)[0],
            "antibonding_at_least_4_dimension": transformed_census(L, sector)[1],
        }
        for L, sector in early_specs
    }
    check(
        "V4_lexicographically_first_omitted_line",
        failure["earlier_transformed_census"] == early
        and all(
            row["antibonding_at_least_4_dimension"] == 0
            for key, row in early.items()
            if key != "L4_k4"
        )
        and early["L4_k4"]["antibonding_at_least_4_dimension"] == 1,
    )

    witness = witness_block()
    stored_witness = {
        decode_interleaved(int(term["configuration"]), 4): int(term["coefficient"])
        for term in failure["physical_witness"]
    }
    invariant = all(
        witness.get(image, 0) == coefficient
        for config, coefficient in witness.items()
        for image in block_orbit(*config, 4)
    )
    W4_sector = [
        vector
        for vector in records[4]["basis"]
        if int(vector[0]["representative"]).bit_count() == 4
    ]
    pairings = [pairing_from_block_witness(witness, vector) for vector in W4_sector]
    check(
        "V4_antibonding_witness_rebuilt",
        stored_witness == witness
        and len(witness) == 16
        and invariant
        and all(top.bit_count() + bottom.bit_count() == 4 for top, bottom in witness),
    )
    check(
        "V4_witness_exactly_in_U",
        len(W4_sector) == 2
        and pairings == [0, 0]
        and failure["pairings_with_complete_exact_W4_sector_basis"] == [0, 0]
        and failure["witness_in_U_exact_Q"] is True,
    )
    check(
        "V4_exact_kernel_and_route_closure",
        failure["omitted_K_dimension"] == 1
        and failure["kernel_dimension_exact_Q"] == 1
        and failure["restriction_rank_exact_Q"] == 19
        and failure["omitted_basis_word"] == ["a", "a", "a", "a"],
    )


def stage_scope_and_resources(envelope: dict[str, object]) -> None:
    missing = envelope["data"]["single_missing_all_L_lemma"]
    scope = envelope["data"]["scope"]
    check(
        "V5_no_all_L_promotion",
        missing["tag"] == "[UNRESOLVED]"
        and any("[UNRESOLVED]" in line and "L>=9" in line for line in scope)
        and envelope["data"]["coverage_contract"]["all_L_computation_claimed"] is False,
    )
    rss = int(envelope["meta"]["peak_rss_bytes"])
    cap = int(envelope["meta"]["rss_cap_bytes"])
    check("V5_resource_ceiling", rss < cap == 2 * 1024**3, f"{rss} < {cap}")
    current_rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform != "darwin":
        current_rss *= 1024
    check("V5_verifier_resource_ceiling", current_rss < cap, str(current_rss))


def main() -> int:
    envelope = json.loads(ARTIFACT.read_text())
    source = json.loads(SOURCE.read_text())
    stage_schema_and_coverage(envelope)
    records = stage_source_and_observability(envelope, source)
    stage_exceptional(envelope)
    stage_bonding(envelope, records)
    stage_scope_and_resources(envelope)
    print(
        f"[test_wlaw_rank_mechanism] {CHECK_COUNT} checks, "
        f"{len(FAILURES)} failures"
    )
    if FAILURES:
        print("[test_wlaw_rank_mechanism] FAILED: " + ", ".join(FAILURES))
        return 1
    print("[test_wlaw_rank_mechanism] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
