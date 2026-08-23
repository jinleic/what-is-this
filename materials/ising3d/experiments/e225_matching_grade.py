"""Test candidate invariants and Jordan--Wigner grade profiles exactly."""

from __future__ import annotations

import itertools
import math
import time
from collections import Counter
from typing import Iterable, Sequence

from e224_nonhamiltonian_census import (
    CPU_BUDGET_SECONDS,
    RSS_CAP_BYTES,
    budget_tick,
    edges_from_code,
    exact_pauli_closure,
    max_rss_bytes,
    run_census,
)


def path_majorana_masks(n: int, order: Sequence[int]) -> tuple[int, ...]:
    if len(order) != n or set(order) != set(range(n)):
        raise ValueError("order must be a vertex permutation")
    answer: list[int] = []
    prefix_x = 0
    for vertex in order:
        answer.append(prefix_x | (1 << (n + vertex)))
        answer.append(prefix_x | (1 << vertex) | (1 << (n + vertex)))
        prefix_x ^= 1 << vertex
    return tuple(answer)


def inverse_majorana_basis(
    n: int, order: Sequence[int]
) -> dict[int, tuple[int, int]]:
    rows: dict[int, tuple[int, int]] = {}
    for index, pauli in enumerate(path_majorana_masks(n, order)):
        coefficients = 1 << index
        while pauli:
            pivot = pauli.bit_length() - 1
            if pivot not in rows:
                rows[pivot] = (pauli, coefficients)
                break
            row, row_coefficients = rows[pivot]
            pauli ^= row
            coefficients ^= row_coefficients
        else:
            raise AssertionError("Jordan--Wigner Majoranas are dependent")
    return rows


def majorana_grade(pauli: int, inverse: dict[int, tuple[int, int]]) -> int:
    coefficients = 0
    while pauli:
        pivot = pauli.bit_length() - 1
        row, row_coefficients = inverse[pivot]
        pauli ^= row
        coefficients ^= row_coefficients
    return coefficients.bit_count()


def histogram_key(
    values: Iterable[int], grade_lookup: dict[int, int]
) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(grade_lookup[value] for value in values).items()))


def graph_summary(row: dict[str, object]) -> dict[str, object]:
    keys = (
        "n",
        "canonical_code",
        "canonical_code_hex",
        "edges",
        "bipartition_sizes",
        "colour_imbalance",
        "smaller_colour_parity",
        "hamiltonian_path_exists",
        "maximum_matching_size",
        "matching_deficiency",
        "minimum_path_cover",
        "dimension",
        "closure_sha256",
    )
    return {key: row[key] for key in keys}


def exact_grade_profiles(
    graph_rows: Sequence[dict[str, object]], started: float
) -> list[dict[str, object]]:
    nonhamiltonian = [
        row for row in graph_rows if row["hamiltonian_path_exists"] is False
    ]
    closures: dict[tuple[int, int], frozenset[int]] = {}
    for row in nonhamiltonian:
        n = int(row["n"])
        code = int(row["canonical_code"])
        closures[(n, code)] = exact_pauli_closure(n, edges_from_code(n, code))

    profile_counts: dict[
        tuple[int, int], Counter[tuple[tuple[int, int], ...]]
    ] = {key: Counter() for key in closures}
    profile_examples: dict[
        tuple[int, int], dict[tuple[tuple[int, int], ...], tuple[int, ...]]
    ] = {key: {} for key in closures}

    for n in sorted({key[0] for key in closures}):
        keys = tuple(key for key in closures if key[0] == n)
        union = frozenset().union(*(closures[key] for key in keys))
        for order in itertools.permutations(range(n)):
            inverse = inverse_majorana_basis(n, order)
            grade_lookup = {
                value: majorana_grade(value, inverse) for value in union
            }
            for key in keys:
                histogram = histogram_key(closures[key], grade_lookup)
                profile_counts[key][histogram] += 1
                profile_examples[key].setdefault(histogram, order)
        budget_tick(started, f"all-order grade profiles n={n}")

    rows_by_key = {
        (int(row["n"]), int(row["canonical_code"])): row
        for row in nonhamiltonian
    }
    answer: list[dict[str, object]] = []
    for key in sorted(closures):
        n, code = key
        profiles = []
        for histogram, count in sorted(profile_counts[key].items()):
            profiles.append(
                {
                    "histogram": {str(grade): amount for grade, amount in histogram},
                    "grades": [grade for grade, _ in histogram],
                    "order_count": count,
                    "example_order": list(profile_examples[key][histogram]),
                }
            )
        source = rows_by_key[key]
        answer.append(
            {
                "tag": "[COMPUTATION]",
                "n": n,
                "canonical_code": code,
                "canonical_code_hex": source["canonical_code_hex"],
                "edges": source["edges"],
                "bipartition_sizes": source["bipartition_sizes"],
                "dimension": source["dimension"],
                "closure_sha256": source["closure_sha256"],
                "orders_audited": math.factorial(n),
                "distinct_histogram_count": len(profiles),
                "profiles": profiles,
            }
        )
    return answer


def matching_path_counterexample(
    rows: Sequence[dict[str, object]],
) -> dict[str, object]:
    nonhamiltonian = [
        row for row in rows if row["hamiltonian_path_exists"] is False
    ]
    for left, right in itertools.combinations(nonhamiltonian, 2):
        if (
            left["n"] == right["n"]
            and left["matching_deficiency"] == right["matching_deficiency"]
            and left["minimum_path_cover"] == right["minimum_path_cover"]
            and left["dimension"] != right["dimension"]
        ):
            return {
                "tag": "[COMPUTATION]",
                "claim": "matching deficiency together with path-cover number does not determine the dimension",
                "left": graph_summary(left),
                "right": graph_summary(right),
            }
    raise AssertionError("no matching/path-cover counterexample found")


def colour_magnitude_counterexample(
    rows: Sequence[dict[str, object]],
) -> dict[str, object]:
    nonhamiltonian = [
        row for row in rows if row["hamiltonian_path_exists"] is False
    ]
    for left, right in itertools.combinations(nonhamiltonian, 2):
        if (
            int(left["n"]) % 2 == 0
            and left["n"] == right["n"]
            and left["smaller_colour_parity"] == right["smaller_colour_parity"]
            and left["colour_imbalance"] != right["colour_imbalance"]
            and left["dimension"] == right["dimension"]
        ):
            return {
                "tag": "[COMPUTATION]",
                "claim": "colour-imbalance magnitude is finer than necessary; only smaller-class parity enters at even n",
                "left": graph_summary(left),
                "right": graph_summary(right),
            }
    raise AssertionError("no colour-magnitude counterexample found")


def order_dependence_counterexample(
    profiles: Sequence[dict[str, object]],
) -> dict[str, object]:
    for row in profiles:
        grade_sets = {tuple(profile["grades"]) for profile in row["profiles"]}
        if len(grade_sets) > 1:
            first, second = next(
                (left, right)
                for left, right in itertools.combinations(row["profiles"], 2)
                if left["grades"] != right["grades"]
            )
            return {
                "tag": "[COMPUTATION]",
                "claim": "reachable Clifford grades depend on the Jordan--Wigner vertex order",
                "graph": {
                    "n": row["n"],
                    "canonical_code": row["canonical_code"],
                    "canonical_code_hex": row["canonical_code_hex"],
                    "edges": row["edges"],
                    "dimension": row["dimension"],
                },
                "first_profile": first,
                "second_profile": second,
            }
    raise AssertionError("no order-dependent grade set found")


def grade_set_dimension_counterexample(
    profiles: Sequence[dict[str, object]],
) -> dict[str, object]:
    buckets: dict[tuple[int, tuple[int, ...]], list[tuple[dict[str, object], dict[str, object]]]] = {}
    for row in profiles:
        for profile in row["profiles"]:
            key = (int(row["n"]), tuple(profile["grades"]))
            buckets.setdefault(key, []).append((row, profile))
    for (n, grades), entries in sorted(buckets.items()):
        for left, right in itertools.combinations(entries, 2):
            left_row, left_profile = left
            right_row, right_profile = right
            if left_row["dimension"] != right_row["dimension"]:
                return {
                    "tag": "[COMPUTATION]",
                    "claim": "the set of occupied Clifford grades does not determine the dimension",
                    "n": n,
                    "common_grades": list(grades),
                    "left": {
                        "canonical_code": left_row["canonical_code"],
                        "canonical_code_hex": left_row["canonical_code_hex"],
                        "edges": left_row["edges"],
                        "dimension": left_row["dimension"],
                        "profile": left_profile,
                    },
                    "right": {
                        "canonical_code": right_row["canonical_code"],
                        "canonical_code_hex": right_row["canonical_code_hex"],
                        "edges": right_row["edges"],
                        "dimension": right_row["dimension"],
                        "profile": right_profile,
                    },
                }
    raise AssertionError("no same-grade-set/different-dimension pair found")


def run_matching_grade(census: dict[str, object]) -> dict[str, object]:
    started = time.process_time()
    rows = census["graphs"]
    if not isinstance(rows, list):
        raise TypeError("census graph rows must be a list")
    profiles = exact_grade_profiles(rows, started)
    matching_path = matching_path_counterexample(rows)
    colour_magnitude = colour_magnitude_counterexample(rows)
    order_dependence = order_dependence_counterexample(profiles)
    grade_set_dimension = grade_set_dimension_counterexample(profiles)

    nonhamiltonian = [
        row for row in rows if row["hamiltonian_path_exists"] is False
    ]
    checks = [
        {
            "name": "every non-Hamiltonian order was audited",
            "passed": len(profiles) == len(nonhamiltonian)
            and all(
                sum(profile["order_count"] for profile in row["profiles"])
                == math.factorial(int(row["n"]))
                for row in profiles
            ),
            "detail": f"{len(profiles)} graphs, all vertex orders",
        },
        {
            "name": "matching/path-cover candidate is decisively refuted",
            "passed": matching_path["left"]["n"] == matching_path["right"]["n"]
            and matching_path["left"]["matching_deficiency"]
            == matching_path["right"]["matching_deficiency"]
            and matching_path["left"]["minimum_path_cover"]
            == matching_path["right"]["minimum_path_cover"]
            and matching_path["left"]["dimension"]
            != matching_path["right"]["dimension"],
            "detail": "same n, matching deficiency, and path cover; different exact dimensions",
        },
        {
            "name": "grade set is presentation data, not an intrinsic classifier",
            "passed": order_dependence["first_profile"]["grades"]
            != order_dependence["second_profile"]["grades"],
            "detail": "one graph has two occupied-grade sets under different orders",
        },
        {
            "name": "grade set alone does not determine dimension",
            "passed": grade_set_dimension["left"]["dimension"]
            != grade_set_dimension["right"]["dimension"]
            and grade_set_dimension["left"]["profile"]["grades"]
            == grade_set_dimension["right"]["profile"]["grades"],
            "detail": "same n and occupied grades; different exact dimensions",
        },
        {
            "name": "colour magnitude is narrowed to parity",
            "passed": colour_magnitude["left"]["dimension"]
            == colour_magnitude["right"]["dimension"]
            and colour_magnitude["left"]["colour_imbalance"]
            != colour_magnitude["right"]["colour_imbalance"]
            and colour_magnitude["left"]["smaller_colour_parity"]
            == colour_magnitude["right"]["smaller_colour_parity"],
            "detail": "same even n and smaller-class parity, different imbalance, same dimension",
        },
        {
            "name": "resource cap",
            "passed": max_rss_bytes() < RSS_CAP_BYTES,
            "detail": f"peak_rss_bytes={max_rss_bytes()}",
        },
    ]
    budget_tick(started, "candidate invariant audit")
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError(checks)
    return {
        "tag": "[COMPUTATION]",
        "grade_profiles": profiles,
        "counterexamples": {
            "matching_deficiency_and_path_cover": matching_path,
            "colour_imbalance_magnitude": colour_magnitude,
            "order_dependent_grades": order_dependence,
            "same_grade_set_different_dimension": grade_set_dimension,
        },
        "checks": checks,
        "process_time_seconds": time.process_time() - started,
        "peak_rss_bytes": max_rss_bytes(),
        "budgets": {
            "process_time_seconds": CPU_BUDGET_SECONDS,
            "rss_cap_bytes": RSS_CAP_BYTES,
        },
    }


def main() -> int:
    census = run_census()
    result = run_matching_grade(census)
    print(
        "PASS "
        f"(graphs={len(result['grade_profiles'])}, "
        f"process_time={result['process_time_seconds']:.6f}s, "
        f"peak_rss={result['peak_rss_bytes']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
