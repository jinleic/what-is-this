#!/usr/bin/env python3
"""Integral block-minor certificate for W-law two-slice observability.

The exact annihilator bases already stored for L=3..9 define the weighted
restriction map W_L^(2m) -> H_Lm^*, where H_Lm contains the orbit coordinates
with absolute leg imbalance at least four.  This producer extracts a
field-independent triangular leaf minor and then computes the determinant of
the small residual core exactly over Z.  Modular ranks are controls only; every
rational rank claim is certified by a nonzero integral minor and the omitted
coordinate count.

No ladder closure is run, no dimension formula is fitted, and no benchmark
coupling enters the calculation.
"""
from __future__ import annotations

import hashlib
import heapq
import json
import resource
import sys
import time
from itertools import combinations
from math import comb
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e237_wlaw_exact_sequence.py"
OUT_PATH = ROOT / "results" / "algebra_growth" / "wlaw_exact_sequence.json"
L8_PATH = ROOT / "results" / "ladder" / "l8_saturation.json"
W9_META_PATH = ROOT / "results" / "ladder" / "w9_saturation.json"
W9_BASIS_PATH = ROOT / "results" / "ladder" / "w9_basis.json"
L_VALUES = tuple(range(3, 10))
PRIMES = (999_983, 1_000_003)
RSS_CAP_BYTES = 2 * 1024**3
EXPECTED_CORE_DETERMINANTS = {
    (8, 4): -(2**32) * 3,
    (9, 4): -(2**288) * 3**5,
    (9, 5): -(2**288) * 3**5,
}


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def integer_list_sha256(values: Iterable[int]) -> str:
    text = ",".join(str(int(value)) for value in values)
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def canonical_basis_sha256(basis: list[list[dict[str, int]]]) -> str:
    payload = json.dumps(basis, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return value == divisor
        divisor += 1 if divisor == 2 else 2
    return True


def tau_config(config: int, L: int) -> int:
    image = 0
    for rung in range(L):
        image |= ((config >> (2 * rung)) & 1) << (2 * rung + 1)
        image |= ((config >> (2 * rung + 1)) & 1) << (2 * rung)
    return image


def rho_config(config: int, L: int) -> int:
    image = 0
    for rung in range(L):
        image |= ((config >> (2 * rung)) & 3) << (2 * (L - 1 - rung))
    return image


def absolute_leg_imbalance(config: int, L: int) -> int:
    top_count = sum((config >> (2 * rung)) & 1 for rung in range(L))
    return abs(2 * top_count - config.bit_count())


def configuration_orbits(L: int) -> list[dict[str, int]]:
    seen: set[int] = set()
    records: list[dict[str, int]] = []
    for config in range(1 << (2 * L)):
        if config.bit_count() % 2 or config in seen:
            continue
        reflected = rho_config(config, L)
        orbit = {
            config,
            tau_config(config, L),
            reflected,
            tau_config(reflected, L),
        }
        seen.update(orbit)
        representative = min(orbit)
        records.append(
            {
                "representative": representative,
                "orbit_size": len(orbit),
                "particle_sector": representative.bit_count(),
                "absolute_leg_imbalance": absolute_leg_imbalance(
                    representative, L
                ),
            }
        )
    return records


def l3_l8_records(source: dict[str, object]) -> dict[int, dict[str, object]]:
    data = source["data"]
    records = {
        int(record["L"]): record for record in data["regression_records"]
    }
    records[7] = data["L7_record"]
    records[8] = data["L8_record"]
    return records


def l9_record(
    metadata: dict[str, object], basis: list[list[dict[str, int]]]
) -> dict[str, object]:
    data = metadata["data"]
    source = data["W9"]
    return {
        "L": 9,
        "K_dim": data["K9_dimension_record"]["K_dim_formula"],
        "W_dim": source["value"],
        "W_sectors": source["W_sectors"],
        "cyclic_Q_dim": source["cyclic_Q_dim"],
        "basis": basis,
        "validation": source["validation"],
        "claim_tag": source["claim_tag"],
    }


def source_exact(record: dict[str, object]) -> bool:
    validation = record["validation"]
    return (
        len(record["basis"]) == int(record["W_dim"])
        and sum(int(value) for value in record["W_sectors"].values())
        == int(record["W_dim"])
        and int(record["cyclic_Q_dim"]) + int(record["W_dim"])
        == int(record["K_dim"])
        and bool(validation["A_invariant_over_Q"])
        and bool(validation["basis_independent_over_Q"])
        and bool(validation["four_B_invariant_over_Q"])
        and bool(validation["sector_homogeneous"])
        and bool(validation["vacuum_orthogonal"])
        and not validation["outside_image_indices"]
    )


def weighted_restriction_rows(
    L: int,
    m: int,
    orbit_records: list[dict[str, int]],
    annihilator_record: dict[str, object],
) -> tuple[list[dict[str, int]], list[dict[str, object]]]:
    sector = 2 * m
    high_orbits = [
        row
        for row in orbit_records
        if row["particle_sector"] == sector
        and row["absolute_leg_imbalance"] >= 4
    ]
    high_index = {
        row["representative"]: index for index, row in enumerate(high_orbits)
    }
    high_sizes = {
        row["representative"]: row["orbit_size"] for row in high_orbits
    }
    rows: list[dict[str, object]] = []
    for basis_index, vector in enumerate(annihilator_record["basis"]):
        if not vector or int(vector[0]["representative"]).bit_count() != sector:
            continue
        entries: dict[int, int] = {}
        for term in vector:
            representative = int(term["representative"])
            column = high_index.get(representative)
            if column is not None:
                entries[column] = (
                    int(term["coefficient"]) * high_sizes[representative]
                )
        rows.append({"basis_index": basis_index, "entries": entries})
    return high_orbits, rows


def leaf_certificate(
    rows: list[dict[str, object]], high_orbits: list[dict[str, int]]
) -> tuple[list[dict[str, int]], set[int], set[int]]:
    supports = [set(row["entries"]) for row in rows]
    by_column: list[list[int]] = [[] for _ in high_orbits]
    for row_index, support in enumerate(supports):
        for column in support:
            by_column[column].append(row_index)

    remaining = set(range(len(high_orbits)))
    used_rows: set[int] = set()
    queue = [index for index, support in enumerate(supports) if len(support) == 1]
    heapq.heapify(queue)
    selected: list[dict[str, int]] = []
    while queue:
        row_index = heapq.heappop(queue)
        if row_index in used_rows:
            continue
        active = supports[row_index] & remaining
        if len(active) != 1:
            continue
        pivot = next(iter(active))
        entries = rows[row_index]["entries"]
        diagonal = int(entries[pivot])
        if diagonal == 0:
            raise AssertionError("leaf diagonal unexpectedly vanished")
        selected.append(
            {
                "source_basis_index": int(rows[row_index]["basis_index"]),
                "pivot_representative": int(
                    high_orbits[pivot]["representative"]
                ),
                "diagonal_entry": diagonal,
            }
        )
        used_rows.add(row_index)
        remaining.remove(pivot)
        for neighbour in by_column[pivot]:
            if neighbour not in used_rows and len(
                supports[neighbour] & remaining
            ) == 1:
                heapq.heappush(queue, neighbour)
    return selected, remaining, used_rows


def replay_leaf_certificate(
    certificate: list[dict[str, int]],
    rows: list[dict[str, object]],
    high_orbits: list[dict[str, int]],
) -> set[int]:
    by_basis_index = {
        int(row["basis_index"]): row["entries"] for row in rows
    }
    by_representative = {
        int(row["representative"]): index
        for index, row in enumerate(high_orbits)
    }
    remaining = set(range(len(high_orbits)))
    for step in certificate:
        entries = by_basis_index[int(step["source_basis_index"])]
        pivot = by_representative[int(step["pivot_representative"])]
        active = set(entries) & remaining
        if active != {pivot} or int(entries[pivot]) != int(
            step["diagonal_entry"]
        ):
            raise AssertionError("invalid stored leaf step")
        remaining.remove(pivot)
    return remaining


def bareiss_determinant(
    rows: list[dict[int, int]], columns: list[int]
) -> int:
    dimension = len(columns)
    if len(rows) != dimension:
        raise ValueError("determinant requires a square row selection")
    if dimension == 0:
        return 1
    position = {column: index for index, column in enumerate(columns)}
    matrix = [[0] * dimension for _ in range(dimension)]
    for row_index, row in enumerate(rows):
        for column, coefficient in row.items():
            target = position.get(column)
            if target is not None:
                matrix[row_index][target] = int(coefficient)

    sign = 1
    previous_pivot = 1
    for pivot_column in range(dimension - 1):
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
            sign = -sign
        pivot = matrix[pivot_column][pivot_column]
        pivot_data = matrix[pivot_column]
        for row_index in range(pivot_column + 1, dimension):
            row = matrix[row_index]
            entry = row[pivot_column]
            for column in range(pivot_column + 1, dimension):
                numerator = row[column] * pivot - entry * pivot_data[column]
                quotient, remainder = divmod(numerator, previous_pivot)
                if remainder:
                    raise AssertionError("Bareiss division was not exact")
                row[column] = quotient
            row[pivot_column] = 0
        previous_pivot = pivot
    return sign * matrix[-1][-1]


def first_nonzero_core_minor(
    rows: list[dict[str, object]],
    remaining: set[int],
    used_rows: set[int],
) -> dict[str, object]:
    columns = sorted(remaining)
    if not columns:
        return {
            "active_source_basis_indices": [],
            "selected_minor_source_basis_indices": [],
            "determinant_Z": 1,
            "active_row_degree_histogram": {},
        }
    active_local_indices = [
        index
        for index, row in enumerate(rows)
        if index not in used_rows and set(row["entries"]) & remaining
    ]
    degree_histogram: dict[int, int] = {}
    for index in active_local_indices:
        degree = len(set(rows[index]["entries"]) & remaining)
        degree_histogram[degree] = degree_histogram.get(degree, 0) + 1
    if min(degree_histogram) < 2:
        raise AssertionError("leaf peeling stopped with a singleton row")
    candidate_count = comb(len(active_local_indices), len(columns))
    if candidate_count > 64:
        raise AssertionError(
            f"core minor search would require {candidate_count} candidates"
        )
    selected: tuple[int, ...] | None = None
    determinant = 0
    for choice in combinations(active_local_indices, len(columns)):
        determinant = bareiss_determinant(
            [rows[index]["entries"] for index in choice], columns
        )
        if determinant:
            selected = choice
            break
    if selected is None:
        raise AssertionError("residual core has no full-size integral minor")
    return {
        "active_source_basis_indices": [
            int(rows[index]["basis_index"]) for index in active_local_indices
        ],
        "selected_minor_source_basis_indices": [
            int(rows[index]["basis_index"]) for index in selected
        ],
        "determinant_Z": determinant,
        "active_row_degree_histogram": {
            str(degree): count
            for degree, count in sorted(degree_histogram.items())
        },
    }


def factor_2_3(value: int) -> dict[str, int]:
    sign = -1 if value < 0 else 1
    residual = abs(value)
    exponent_2 = 0
    exponent_3 = 0
    while residual and residual % 2 == 0:
        residual //= 2
        exponent_2 += 1
    while residual and residual % 3 == 0:
        residual //= 3
        exponent_3 += 1
    return {
        "sign": sign,
        "exponent_2": exponent_2,
        "exponent_3": exponent_3,
        "residual_after_2_and_3": residual,
    }


def sparse_modular_rank(
    rows: Iterable[dict[int, int]], prime: int
) -> tuple[int, tuple[int, ...]]:
    pivots: dict[int, dict[int, int]] = {}
    for source in rows:
        vector = {
            int(column): int(coefficient) % prime
            for column, coefficient in source.items()
            if int(coefficient) % prime
        }
        while vector:
            pivot = max(vector)
            coefficient = vector[pivot]
            old = pivots.get(pivot)
            if old is None:
                inverse = pow(coefficient, prime - 2, prime)
                vector = {
                    column: (entry * inverse) % prime
                    for column, entry in vector.items()
                    if (entry * inverse) % prime
                }
                pivots[pivot] = vector
                break
            for column, entry in old.items():
                updated = (vector.get(column, 0) - coefficient * entry) % prime
                if updated:
                    vector[column] = updated
                else:
                    vector.pop(column, None)
    return len(pivots), tuple(sorted(pivots))


def sector_record(
    L: int,
    m: int,
    orbit_records: list[dict[str, int]],
    annihilator_record: dict[str, object],
) -> dict[str, object]:
    started = time.process_time()
    sector = 2 * m
    sector_orbits = [
        row for row in orbit_records if row["particle_sector"] == sector
    ]
    high_orbits, rows = weighted_restriction_rows(
        L, m, orbit_records, annihilator_record
    )
    leaf, remaining, used_rows = leaf_certificate(rows, high_orbits)
    if replay_leaf_certificate(leaf, rows, high_orbits) != remaining:
        raise AssertionError("leaf certificate replay changed the residual core")
    core = first_nonzero_core_minor(rows, remaining, used_rows)
    core["remaining_coordinate_representatives"] = [
        int(high_orbits[index]["representative"])
        for index in sorted(remaining)
    ]
    core["dimension"] = len(remaining)
    core["determinant_factorization"] = factor_2_3(
        int(core["determinant_Z"])
    )

    modular: dict[str, dict[str, object]] = {}
    raw_rows = [row["entries"] for row in rows]
    for prime in PRIMES:
        rank, pivots = sparse_modular_rank(raw_rows, prime)
        representatives = [
            int(high_orbits[index]["representative"]) for index in pivots
        ]
        modular[str(prime)] = {
            "tag": "[COMPUTATION] control lower bound",
            "rank_Fp_lower_bound_for_rank_Q": rank,
            "rational_rank_ceiling_omitted_coordinate_dimension": len(
                high_orbits
            ),
            "lower_equals_ceiling": rank == len(high_orbits),
            "pivot_representatives_sha256": integer_list_sha256(
                representatives
            ),
        }

    W_by_sector = {
        int(key): int(value)
        for key, value in annihilator_record["W_sectors"].items()
    }
    W_dimension = W_by_sector.get(sector, 0)
    K_dimension = len(sector_orbits)
    H_dimension = len(high_orbits)
    E_dimension = K_dimension - H_dimension
    U_dimension = K_dimension - W_dimension
    integral_minor_rank = len(leaf) + int(core["dimension"])
    exact = (
        integral_minor_rank == H_dimension
        and int(core["determinant_Z"]) != 0
        and all(int(step["diagonal_entry"]) != 0 for step in leaf)
    )
    return {
        "tag": "[COMPUTATION] exact finite-row integral rank certificate",
        "L": L,
        "m": m,
        "particle_sector": sector,
        "K_sector_dimension": K_dimension,
        "retained_abs_leg_imbalance_0_or_2_dimension": E_dimension,
        "omitted_abs_leg_imbalance_at_least_4_dimension": H_dimension,
        "annihilator_W_sector_dimension_exact_Q": W_dimension,
        "annihilator_basis_vectors_in_sector": len(rows),
        "cyclic_U_sector_dimension_exact_Q": U_dimension,
        "weighted_restriction_nnz": sum(len(row["entries"]) for row in rows),
        "weighted_restriction_max_row_nnz": max(
            (len(row["entries"]) for row in rows), default=0
        ),
        "leaf_certificate": {
            "dimension": len(leaf),
            "steps": leaf,
            "step_digest": integer_list_sha256(
                value
                for step in leaf
                for value in (
                    step["source_basis_index"],
                    step["pivot_representative"],
                    step["diagonal_entry"],
                )
            ),
        },
        "residual_core": core,
        "integral_block_minor_rank_lower_bound_for_rank_Q": integral_minor_rank,
        "rational_rank_ceiling_omitted_coordinate_dimension": H_dimension,
        "integral_lower_equals_rational_ceiling": exact,
        "weighted_annihilator_restriction_modular_controls": modular,
        "annihilator_restriction_rank_exact_Q": H_dimension if exact else None,
        "restriction_kernel_on_U_dimension_exact_Q": 0 if exact else None,
        "retained_restriction_rank_exact_Q": U_dimension if exact else None,
        "row_cpu_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes_after_row": peak_rss_bytes(),
        "all_checks_passed": (
            exact
            and len(rows) == W_dimension
            and E_dimension + H_dimension == K_dimension
            and all(
                control["lower_equals_ceiling"]
                for control in modular.values()
            )
        ),
    }


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def main() -> int:
    started = time.process_time()
    l8_source = json.loads(L8_PATH.read_text())
    w9_metadata = json.loads(W9_META_PATH.read_text())
    w9_basis_envelope = json.loads(W9_BASIS_PATH.read_text())
    w9_basis = w9_basis_envelope["basis"]
    records = l3_l8_records(l8_source)
    records[9] = l9_record(w9_metadata, w9_basis)
    checks: list[dict[str, object]] = []

    canonical_w9_hash = canonical_basis_sha256(w9_basis)
    expected_w9_hash = str(w9_metadata["data"]["W9"]["basis_sha256"])
    add_check(
        checks,
        "C1_exact_annihilator_sources_L3_9",
        set(records) == set(L_VALUES)
        and all(source_exact(records[L]) for L in L_VALUES)
        and canonical_w9_hash == expected_w9_hash
        and len(w9_basis) == int(w9_metadata["data"]["W9"]["basis_size"]),
        "all cited W bases are exact-Q invariant annihilators; the separate L9 basis matches its canonical digest",
    )
    add_check(
        checks,
        "C2_control_primes",
        all(is_prime(prime) and prime > 4 for prime in PRIMES),
        f"odd prime controls {PRIMES}; their ranks are not the rational proof",
    )

    records_by_L: dict[int, list[dict[str, object]]] = {}
    orbit_counts: dict[int, int] = {}
    for L in L_VALUES:
        orbits = configuration_orbits(L)
        orbit_counts[L] = len(orbits)
        records_by_L[L] = [
            sector_record(L, m, orbits, records[L]) for m in range(L + 1)
        ]
    sector_records = [
        row for L in L_VALUES for row in records_by_L[L]
    ]

    add_check(
        checks,
        "C3_complete_sector_coverage_L3_9",
        {(row["L"], row["m"]) for row in sector_records}
        == {(L, m) for L in L_VALUES for m in range(L + 1)},
        f"all {len(sector_records)} sectors L=3..9 are present",
    )
    add_check(
        checks,
        "C4_direct_K_orbit_census",
        all(orbit_counts[L] == int(records[L]["K_dim"]) for L in L_VALUES),
        f"direct even-parity orbit dimensions {orbit_counts}",
    )
    add_check(
        checks,
        "C5_leaf_certificates_replay",
        all(
            row["leaf_certificate"]["dimension"]
            + row["residual_core"]["dimension"]
            == row["omitted_abs_leg_imbalance_at_least_4_dimension"]
            for row in sector_records
        ),
        "each stored raw-basis leaf sequence is triangular and leaves exactly the declared core columns",
    )

    nonempty_cores = {
        (int(row["L"]), int(row["m"])): row
        for row in sector_records
        if int(row["residual_core"]["dimension"]) > 0
    }
    add_check(
        checks,
        "C6_exact_integer_core_determinants",
        set(nonempty_cores) == set(EXPECTED_CORE_DETERMINANTS)
        and all(
            int(nonempty_cores[key]["residual_core"]["determinant_Z"])
            == determinant
            for key, determinant in EXPECTED_CORE_DETERMINANTS.items()
        ),
        "the only residual cores are (8,4) of size 17 with determinant -3*2^32 and (9,4),(9,5) of size 133 with determinant -3^5*2^288",
    )
    add_check(
        checks,
        "C7_integral_full_column_minors_L3_9",
        all(
            row["integral_lower_equals_rational_ceiling"]
            and row["annihilator_restriction_rank_exact_Q"]
            == row["omitted_abs_leg_imbalance_at_least_4_dimension"]
            for row in sector_records
        ),
        "nonzero integral block minors meet the omitted-coordinate rational ceiling in every sector",
    )
    for prime in PRIMES:
        add_check(
            checks,
            f"C8_modular_control_full_rank_{prime}",
            all(
                row["weighted_annihilator_restriction_modular_controls"][
                    str(prime)
                ]["lower_equals_ceiling"]
                for row in sector_records
            ),
            "the independent finite-field lower bounds agree with the already exact integral ranks",
        )
    add_check(
        checks,
        "C9_two_slice_restriction_injective_exact_Q_L3_9",
        all(
            row["all_checks_passed"]
            and row["restriction_kernel_on_U_dimension_exact_Q"] == 0
            and row["retained_restriction_rank_exact_Q"]
            == row["cyclic_U_sector_dimension_exact_Q"]
            for row in sector_records
        ),
        "annihilator duality and the integral full-column minors give zero two-slice kernel in all 49 sectors",
    )

    L9_rows = records_by_L[9]
    expected_L9_H = [0, 0, 66, 614, 1689, 1689, 614, 66, 0, 0]
    expected_L9_U = [1, 45, 666, 3570, 8001, 8001, 3570, 666, 45, 1]
    add_check(
        checks,
        "C10_first_exact_row_beyond_L8_is_positive",
        [
            int(row["omitted_abs_leg_imbalance_at_least_4_dimension"])
            for row in L9_rows
        ]
        == expected_L9_H
        and [int(row["retained_restriction_rank_exact_Q"]) for row in L9_rows]
        == expected_L9_U,
        "L9 has no two-slice counterexample: omitted ranks are [0,0,66,614,1689,1689,614,66,0,0] and the retained exact ranks are the cited exact U9 sectors",
    )
    add_check(
        checks,
        "C11_raw_leaf_method_frontier",
        int(nonempty_cores[(8, 4)]["residual_core"]["dimension"]) == 17
        and int(nonempty_cores[(9, 4)]["residual_core"]["dimension"])
        == 133
        and int(nonempty_cores[(9, 5)]["residual_core"]["dimension"])
        == 133
        and all(
            min(
                int(degree)
                for degree in row["residual_core"][
                    "active_row_degree_histogram"
                ]
            )
            >= 2
            for row in nonempty_cores.values()
        ),
        "uncombined source-row leaf peeling first leaves a 17-column core at (8,4) and leaves 133-column cores at (9,4),(9,5); every surviving active row has degree at least two",
    )

    rss = peak_rss_bytes()
    add_check(
        checks,
        "C12_resource_ceiling",
        rss < RSS_CAP_BYTES,
        f"producer peak RSS {rss} bytes is below {RSS_CAP_BYTES}",
    )
    names = [str(check["name"]) for check in checks]
    if len(names) != len(set(names)):
        raise AssertionError("check names are not unique")
    if not all(bool(check["passed"]) for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        raise AssertionError(f"failed checks: {failed}")

    data = {
        "headline": {
            "tag": "[LEMMA][COMPUTATION][UNRESOLVED]",
            "result": (
                "[COMPUTATION] The two-slice restriction is injective exactly over Q "
                "in every sector L=3..9. The new L9 row is proved by integral block "
                "minors: raw leaf pivots close m=2,3,6,7, while the m=4,5 residual "
                "133x133 minors have determinant -2^288*3^5. [UNRESOLVED] This is "
                "not an all-L proof; the exact source frontier ends at L9."
            ),
        },
        "integral_block_minor_lemma": {
            "tag": "[LEMMA]",
            "statement": (
                "For an integer matrix with columns H, repeatedly choose a source row "
                "having one remaining nonzero column. In selection order those rows "
                "form a triangular leaf block with nonzero diagonal and vanish on the "
                "final core columns. If a square minor of the residual core has nonzero "
                "integer determinant, adjoining it gives a nonzero |H|-minor. Hence "
                "rank_Q equals the rational ceiling |H|."
            ),
            "duality_consequence": (
                "For W=U^perp in K=E direct_sum H, full column rank of the weighted "
                "restriction W->H* gives dim ker(R_E|U)=0."
            ),
        },
        "two_slice_observability_exact_L3_9": sector_records,
        "new_L9_sector_records": L9_rows,
        "raw_source_leaf_method_frontier": {
            "tag": "[COMPUTATION][UNRESOLVED]",
            "scope": (
                "The method class uses the stored exact annihilator source vectors "
                "without rational row combinations and permits only successive raw "
                "singleton pivots. Failure to peel a core is a limitation of this "
                "certificate family, not a counterexample to another triangular family."
            ),
            "nonempty_cores": [
                {
                    "L": L,
                    "m": m,
                    "dimension": row["residual_core"]["dimension"],
                    "active_row_degree_histogram": row["residual_core"][
                        "active_row_degree_histogram"
                    ],
                    "determinant_Z": row["residual_core"]["determinant_Z"],
                    "determinant_factorization": row["residual_core"][
                        "determinant_factorization"
                    ],
                }
                for (L, m), row in sorted(nonempty_cores.items())
            ],
            "exact_limitation": (
                "Raw singleton peeling is complete in every controlled sector except "
                "(8,4), (9,4), and (9,5). Thus that route cannot be promoted to an "
                "all-L proof without a separate block-core argument or row combinations."
            ),
        },
        "all_L_status": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "The first exact annihilator row beyond e230, L=9, is positive rather "
                "than a counterexample. No exact annihilator source is available here "
                "for L>=10, and this producer launches no closure. The all-L two-slice "
                "exact-sequence lemma and the first possible counterexample at L>=10 "
                "remain unresolved."
            ),
        },
        "coverage_contract": {
            "tag": "[COMPUTATION]",
            "required_L": list(L_VALUES),
            "required_sector_pairs": [
                [L, m] for L in L_VALUES for m in range(L + 1)
            ],
            "new_exact_L": 9,
            "all_L_claimed": False,
            "closure_run": False,
            "sequence_fit_used": False,
            "benchmark_coupling_used": False,
        },
        "scope": [
            "[LEMMA] the integral leaf-plus-core block-minor criterion is general finite-dimensional linear algebra",
            "[COMPUTATION] exact finite observability rows only, L=3..9, reconstructed from pre-existing exact-Q annihilators",
            "[COMPUTATION] modular ranks are controls and each meets the omitted-coordinate rational ceiling already reached by an integral minor",
            "[UNRESOLVED] no all-size injectivity statement, no L>=10 rank, and no counterexample beyond the exact L9 frontier",
        ],
    }
    envelope = {
        "meta": {
            "script": SCRIPT,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python": sys.version.split()[0],
            "artifact_schema": "meta/data/checks-v1",
            "exact_arithmetic": (
                "Python integers; Bareiss determinants over Z; two finite-field "
                "ranks as controls. Every rational rank equals an integral-minor "
                "lower bound and the omitted-coordinate ceiling."
            ),
            "source_artifacts": {
                str(path.relative_to(ROOT)): file_sha256(path)
                for path in (L8_PATH, W9_META_PATH, W9_BASIS_PATH)
            },
            "canonical_w9_basis_sha256": canonical_w9_hash,
            "primes": list(PRIMES),
            "cpu_seconds": round(time.process_time() - started, 6),
            "peak_rss_bytes": rss,
            "rss_cap_bytes": RSS_CAP_BYTES,
        },
        "data": data,
        "checks": checks,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUT_PATH.with_suffix(".tmp")
    temporary.write_text(json.dumps(envelope, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUT_PATH)
    print(
        f"[e237] wrote {OUT_PATH.relative_to(ROOT)}: "
        f"{len(checks)}/{len(checks)} checks, {len(sector_records)} exact sectors, "
        f"cpu={envelope['meta']['cpu_seconds']}s rss={rss}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
