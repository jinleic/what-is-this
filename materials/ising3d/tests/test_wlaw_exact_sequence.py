#!/usr/bin/env python3
"""Clean-room verifier for the e237 integral exact-sequence certificate.

This file does not import the producer.  It enumerates ladder orbits in
(top-mask, bottom-mask) coordinates, reconstructs the weighted annihilator
restriction, replays the explicit leaf certificates, peels again in the
opposite order, evaluates the residual determinants with Fraction Gaussian
elimination, and uses the smallest-column modular pivot convention.
"""
from __future__ import annotations

import hashlib
import heapq
import json
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "wlaw_exact_sequence.json"
L8_SOURCE = ROOT / "results" / "ladder" / "l8_saturation.json"
W9_METADATA = ROOT / "results" / "ladder" / "w9_saturation.json"
W9_BASIS = ROOT / "results" / "ladder" / "w9_basis.json"
L_VALUES = tuple(range(3, 10))
PRIMES = (999_983, 1_000_003)
RSS_CAP_BYTES = 2 * 1024**3
EXPECTED_PRODUCER_CHECKS = {
    "C1_exact_annihilator_sources_L3_9",
    "C2_control_primes",
    "C3_complete_sector_coverage_L3_9",
    "C4_direct_K_orbit_census",
    "C5_leaf_certificates_replay",
    "C6_exact_integer_core_determinants",
    "C7_integral_full_column_minors_L3_9",
    "C8_modular_control_full_rank_999983",
    "C8_modular_control_full_rank_1000003",
    "C9_two_slice_restriction_injective_exact_Q_L3_9",
    "C10_first_exact_row_beyond_L8_is_positive",
    "C11_raw_leaf_method_frontier",
    "C12_resource_ceiling",
}
EXPECTED_CORE_DETERMINANTS = {
    (8, 4): -(2**32) * 3,
    (9, 4): -(2**288) * 3**5,
    (9, 5): -(2**288) * 3**5,
}

CHECK_COUNT = 0
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    global CHECK_COUNT
    CHECK_COUNT += 1
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] {name}{': ' + detail if detail else ''}")
    if not passed:
        FAILURES.append(name)


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 18), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_basis_digest(basis: list[list[dict[str, int]]]) -> str:
    encoded = json.dumps(basis, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(encoded.encode()).hexdigest()


def reverse_mask(mask: int, L: int) -> int:
    answer = 0
    for rung in range(L):
        answer |= ((mask >> rung) & 1) << (L - 1 - rung)
    return answer


def interleave(top: int, bottom: int, L: int) -> int:
    answer = 0
    for rung in range(L):
        answer |= ((top >> rung) & 1) << (2 * rung)
        answer |= ((bottom >> rung) & 1) << (2 * rung + 1)
    return answer


def pair_orbit_records(L: int) -> list[dict[str, int]]:
    seen: set[tuple[int, int]] = set()
    records: list[dict[str, int]] = []
    for top in range(1 << L):
        for bottom in range(1 << L):
            pair = (top, bottom)
            if (top.bit_count() + bottom.bit_count()) % 2 or pair in seen:
                continue
            reflected = (reverse_mask(top, L), reverse_mask(bottom, L))
            orbit = {
                pair,
                (bottom, top),
                reflected,
                (reflected[1], reflected[0]),
            }
            seen.update(orbit)
            representative = min(
                interleave(left, right, L) for left, right in orbit
            )
            records.append(
                {
                    "representative": representative,
                    "orbit_size": len(orbit),
                    "particle_sector": top.bit_count() + bottom.bit_count(),
                    "absolute_leg_imbalance": abs(
                        top.bit_count() - bottom.bit_count()
                    ),
                }
            )
    records.sort(key=lambda row: row["representative"])
    return records


def source_records(
    l8: dict[str, object],
    w9_metadata: dict[str, object],
    w9_basis: list[list[dict[str, int]]],
) -> dict[int, dict[str, object]]:
    data = l8["data"]
    records = {
        int(record["L"]): record for record in data["regression_records"]
    }
    records[7] = data["L7_record"]
    records[8] = data["L8_record"]
    nine = w9_metadata["data"]["W9"]
    records[9] = {
        "L": 9,
        "K_dim": w9_metadata["data"]["K9_dimension_record"][
            "K_dim_direct_orbits"
        ],
        "W_dim": nine["basis_size"],
        "W_sectors": nine["W_sectors"],
        "cyclic_Q_dim": nine["cyclic_Q_dim"],
        "basis": w9_basis,
        "validation": nine["validation"],
    }
    return records


def exact_source(record: dict[str, object]) -> bool:
    validation = record["validation"]
    return (
        len(record["basis"]) == int(record["W_dim"])
        and sum(int(value) for value in record["W_sectors"].values())
        == int(record["W_dim"])
        and int(record["cyclic_Q_dim"]) + int(record["W_dim"])
        == int(record["K_dim"])
        and all(
            bool(validation[key])
            for key in (
                "A_invariant_over_Q",
                "basis_independent_over_Q",
                "four_B_invariant_over_Q",
                "sector_homogeneous",
                "vacuum_orthogonal",
            )
        )
        and not validation["outside_image_indices"]
    )


def restriction_rows(
    L: int,
    m: int,
    orbits: list[dict[str, int]],
    record: dict[str, object],
) -> tuple[list[dict[str, int]], dict[int, dict[int, int]]]:
    sector = 2 * m
    high = [
        row
        for row in orbits
        if row["particle_sector"] == sector
        and row["absolute_leg_imbalance"] >= 4
    ]
    sizes = {row["representative"]: row["orbit_size"] for row in high}
    rows: dict[int, dict[int, int]] = {}
    for basis_index, vector in enumerate(record["basis"]):
        if not vector or int(vector[0]["representative"]).bit_count() != sector:
            continue
        entries = {
            representative: int(term["coefficient"]) * sizes[representative]
            for term in vector
            if (representative := int(term["representative"])) in sizes
        }
        rows[basis_index] = entries
    return high, rows


def replay_stored_leaf(
    stored: dict[str, object],
    rows: dict[int, dict[int, int]],
    high: list[dict[str, int]],
) -> set[int]:
    remaining = {int(row["representative"]) for row in high}
    for step in stored["steps"]:
        basis_index = int(step["source_basis_index"])
        pivot = int(step["pivot_representative"])
        active = set(rows[basis_index]) & remaining
        if active != {pivot}:
            raise AssertionError(
                f"stored leaf row {basis_index} has active support {active}"
            )
        if rows[basis_index][pivot] != int(step["diagonal_entry"]):
            raise AssertionError("stored leaf diagonal changed")
        remaining.remove(pivot)
    return remaining


def reverse_leaf_core(
    rows: dict[int, dict[int, int]], high: list[dict[str, int]]
) -> set[int]:
    remaining = {int(row["representative"]) for row in high}
    by_column: dict[int, list[int]] = {column: [] for column in remaining}
    for basis_index, row in rows.items():
        for column in row:
            by_column[column].append(basis_index)
    queue = [
        -basis_index
        for basis_index, row in rows.items()
        if len(set(row) & remaining) == 1
    ]
    heapq.heapify(queue)
    used: set[int] = set()
    while queue:
        basis_index = -heapq.heappop(queue)
        if basis_index in used:
            continue
        active = set(rows[basis_index]) & remaining
        if len(active) != 1:
            continue
        pivot = next(iter(active))
        used.add(basis_index)
        remaining.remove(pivot)
        for neighbour in by_column[pivot]:
            if neighbour not in used and len(
                set(rows[neighbour]) & remaining
            ) == 1:
                heapq.heappush(queue, -neighbour)
    return remaining


def fraction_determinant(
    rows: list[dict[int, int]], columns: list[int]
) -> int:
    dimension = len(columns)
    if len(rows) != dimension:
        raise ValueError("square determinant expected")
    if not columns:
        return 1
    position = {column: index for index, column in enumerate(columns)}
    matrix = [
        [Fraction(0) for _ in range(dimension)] for _ in range(dimension)
    ]
    for row_index, row in enumerate(rows):
        for column, coefficient in row.items():
            target = position.get(column)
            if target is not None:
                matrix[row_index][target] = Fraction(coefficient)

    determinant = Fraction(1)
    for pivot_column in range(dimension):
        pivot_row = next(
            (
                row
                for row in range(pivot_column, dimension)
                if matrix[row][pivot_column]
            ),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_column:
            matrix[pivot_column], matrix[pivot_row] = (
                matrix[pivot_row],
                matrix[pivot_column],
            )
            determinant = -determinant
        pivot = matrix[pivot_column][pivot_column]
        determinant *= pivot
        for column in range(pivot_column, dimension):
            matrix[pivot_column][column] /= pivot
        for row_index in range(pivot_column + 1, dimension):
            multiplier = matrix[row_index][pivot_column]
            if not multiplier:
                continue
            for column in range(pivot_column, dimension):
                matrix[row_index][column] -= (
                    multiplier * matrix[pivot_column][column]
                )
    if determinant.denominator != 1:
        raise AssertionError("integer matrix determinant was not integral")
    return determinant.numerator


def modular_rank_smallest_pivot(
    rows: list[dict[int, int]], prime: int
) -> int:
    pivots: dict[int, dict[int, int]] = {}
    for source in rows:
        vector = {
            column: coefficient % prime
            for column, coefficient in source.items()
            if coefficient % prime
        }
        while vector:
            pivot = min(vector)
            coefficient = vector[pivot]
            old = pivots.get(pivot)
            if old is None:
                inverse = pow(coefficient, prime - 2, prime)
                pivots[pivot] = {
                    column: entry * inverse % prime
                    for column, entry in vector.items()
                    if entry * inverse % prime
                }
                break
            for column, entry in old.items():
                updated = (vector.get(column, 0) - coefficient * entry) % prime
                if updated:
                    vector[column] = updated
                else:
                    vector.pop(column, None)
    return len(pivots)


def main() -> int:
    started = time.process_time()
    envelope = json.loads(ARTIFACT.read_text())
    l8 = json.loads(L8_SOURCE.read_text())
    w9_metadata = json.loads(W9_METADATA.read_text())
    w9_basis_envelope = json.loads(W9_BASIS.read_text())
    w9_basis = w9_basis_envelope["basis"]
    records = source_records(l8, w9_metadata, w9_basis)

    check(
        "V1 artifact schema and producer checks",
        set(envelope) == {"meta", "data", "checks"}
        and {row["name"] for row in envelope["checks"]}
        == EXPECTED_PRODUCER_CHECKS
        and all(row["passed"] for row in envelope["checks"]),
    )
    expected_hashes = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (L8_SOURCE, W9_METADATA, W9_BASIS)
    }
    check(
        "V2 pinned source hashes and exact annihilators",
        envelope["meta"]["source_artifacts"] == expected_hashes
        and set(records) == set(L_VALUES)
        and all(exact_source(records[L]) for L in L_VALUES)
        and canonical_basis_digest(w9_basis)
        == w9_metadata["data"]["W9"]["basis_sha256"]
        == envelope["meta"]["canonical_w9_basis_sha256"],
    )

    stored_rows = {
        (int(row["L"]), int(row["m"])): row
        for row in envelope["data"]["two_slice_observability_exact_L3_9"]
    }
    complete_pairs = {(L, m) for L in L_VALUES for m in range(L + 1)}
    check(
        "V3 complete finite coverage",
        set(stored_rows) == complete_pairs and len(stored_rows) == 49,
    )

    replay_ok = True
    reverse_core_ok = True
    determinant_ok = True
    modular_ok = True
    dimensions_ok = True
    orbit_totals: dict[int, int] = {}
    observed_cores: dict[tuple[int, int], tuple[int, int]] = {}
    for L in L_VALUES:
        orbits = pair_orbit_records(L)
        orbit_totals[L] = len(orbits)
        for m in range(L + 1):
            stored = stored_rows[(L, m)]
            high, rows = restriction_rows(L, m, orbits, records[L])
            high_representatives = {
                int(row["representative"]) for row in high
            }
            try:
                replay_core = replay_stored_leaf(
                    stored["leaf_certificate"], rows, high
                )
            except (AssertionError, KeyError):
                replay_ok = False
                replay_core = set()
            stored_core = {
                int(value)
                for value in stored["residual_core"][
                    "remaining_coordinate_representatives"
                ]
            }
            replay_ok &= replay_core == stored_core
            reverse_core_ok &= reverse_leaf_core(rows, high) == stored_core

            selected = [
                rows[int(index)]
                for index in stored["residual_core"][
                    "selected_minor_source_basis_indices"
                ]
            ]
            columns = sorted(stored_core)
            try:
                determinant = fraction_determinant(selected, columns)
            except (AssertionError, KeyError, ValueError):
                determinant = 0
            determinant_ok &= determinant == int(
                stored["residual_core"]["determinant_Z"]
            )
            if columns:
                observed_cores[(L, m)] = (len(columns), determinant)

            raw_rows = list(rows.values())
            for prime in PRIMES:
                rank = modular_rank_smallest_pivot(raw_rows, prime)
                control = stored[
                    "weighted_annihilator_restriction_modular_controls"
                ][str(prime)]
                modular_ok &= (
                    rank
                    == int(control["rank_Fp_lower_bound_for_rank_Q"])
                    == len(high)
                    == int(
                        control[
                            "rational_rank_ceiling_omitted_coordinate_dimension"
                        ]
                    )
                )

            sector = 2 * m
            sector_K = sum(
                row["particle_sector"] == sector for row in orbits
            )
            W_dimension = int(records[L]["W_sectors"].get(str(sector), 0))
            dimensions_ok &= (
                sector_K == int(stored["K_sector_dimension"])
                and len(high)
                == int(
                    stored[
                        "omitted_abs_leg_imbalance_at_least_4_dimension"
                    ]
                )
                and len(rows)
                == W_dimension
                == int(stored["annihilator_W_sector_dimension_exact_Q"])
                and sector_K - W_dimension
                == int(stored["cyclic_U_sector_dimension_exact_Q"])
                == int(stored["retained_restriction_rank_exact_Q"])
                and int(stored["restriction_kernel_on_U_dimension_exact_Q"])
                == 0
                and int(stored["leaf_certificate"]["dimension"])
                + len(stored_core)
                == len(high)
                and high_representatives
                == {
                    int(row["representative"])
                    for row in high
                }
            )

    check(
        "V4 independent top-bottom orbit census",
        orbit_totals
        == {3: 14, 4: 44, 5: 152, 6: 560, 7: 2144, 8: 8384, 9: 33152},
        str(orbit_totals),
    )
    check(
        "V5 explicit leaf certificates replay",
        replay_ok,
        "each selected source row has one active pivot with the stored integer diagonal",
    )
    check(
        "V6 opposite-order leaf peeling has the same cores",
        reverse_core_ok,
        "largest source row first, independently reconstructed",
    )
    check(
        "V7 Fraction core determinants",
        determinant_ok
        and observed_cores
        == {
            key: (
                17 if key == (8, 4) else 133,
                determinant,
            )
            for key, determinant in EXPECTED_CORE_DETERMINANTS.items()
        },
        str(observed_cores),
    )
    check(
        "V8 opposite-pivot modular controls",
        modular_ok,
        "smallest-column ranks meet the omitted-coordinate ceiling at both primes",
    )
    check(
        "V9 exact two-sided sector ranks",
        dimensions_ok,
        "integral lower minors equal H and source exact-Q dimensions give the retained ranks",
    )

    L9 = [stored_rows[(9, m)] for m in range(10)]
    check(
        "V10 L9 first-beyond-L8 positive row",
        [
            int(row["omitted_abs_leg_imbalance_at_least_4_dimension"])
            for row in L9
        ]
        == [0, 0, 66, 614, 1689, 1689, 614, 66, 0, 0]
        and [int(row["retained_restriction_rank_exact_Q"]) for row in L9]
        == [1, 45, 666, 3570, 8001, 8001, 3570, 666, 45, 1],
    )
    scope = envelope["data"]["coverage_contract"]
    check(
        "V11 honest finite scope and resources",
        scope["required_L"] == list(L_VALUES)
        and scope["new_exact_L"] == 9
        and scope["all_L_claimed"] is False
        and scope["closure_run"] is False
        and scope["sequence_fit_used"] is False
        and scope["benchmark_coupling_used"] is False
        and int(envelope["meta"]["peak_rss_bytes"]) < RSS_CAP_BYTES
        and peak_rss_bytes() < RSS_CAP_BYTES,
        f"cpu={time.process_time() - started:.6f}s rss={peak_rss_bytes()}",
    )

    if FAILURES:
        print(f"[test_wlaw_exact_sequence] {len(FAILURES)} failures: {FAILURES}")
        return 1
    print(
        f"[test_wlaw_exact_sequence] {CHECK_COUNT}/{CHECK_COUNT} checks passed; "
        f"cpu={time.process_time() - started:.6f}s rss={peak_rss_bytes()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
