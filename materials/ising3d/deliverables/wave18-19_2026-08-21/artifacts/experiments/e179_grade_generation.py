"""Exact support calculus for the Clifford-grade generation claim.

The originally proposed lemma is deliberately tested at its endpoint.  The
calculation isolates the exceptional seed ``k=N-2`` and proves the corrected
``k<=N-4`` statement by explicit overlap witnesses.
"""

from __future__ import annotations

import json
import time
from itertools import combinations
from typing import Iterable

from e178_clifford_setup import budget_tick, max_rss_bytes


def bracket_grade_outputs(N: int, a: int, b: int) -> tuple[int, ...]:
    """Grades occurring for brackets of coordinate a- and b-monomials."""
    if not (0 <= a <= N and 0 <= b <= N):
        raise ValueError((N, a, b))
    lower = max(0, a + b - N)
    upper = min(a, b)
    return tuple(
        sorted(
            {
                a + b - 2 * overlap
                for overlap in range(lower, upper + 1)
                if (a * b - overlap) % 2
            }
        )
    )


def grade_closure(N: int, seeds: Iterable[int]) -> tuple[int, ...]:
    grades = set(int(grade) for grade in seeds)
    changed = True
    while changed:
        changed = False
        snapshot = tuple(sorted(grades))
        for a in snapshot:
            for b in snapshot:
                for output in bracket_grade_outputs(N, a, b):
                    if output not in grades:
                        grades.add(output)
                        changed = True
    return tuple(sorted(grades))


def target_noncentral_grades(N: int) -> tuple[int, ...]:
    if N % 2:
        raise ValueError("N must be even")
    return tuple(
        grade
        for grade in range(2, N + 1, 4)
        if grade != N
    )


def monomial_product_sign(left: Iterable[int], right: Iterable[int]) -> int:
    """Sign in gamma_left gamma_right = sign gamma_(left symmetric-difference right)."""
    left_tuple = tuple(sorted(left))
    right_tuple = tuple(sorted(right))
    inversions = sum(1 for a in left_tuple for b in right_tuple if a > b)
    return -1 if inversions % 2 else 1


def clifford_bracket(
    left: Iterable[int], right: Iterable[int]
) -> tuple[int, tuple[int, ...]]:
    left_set, right_set = set(left), set(right)
    coefficient = monomial_product_sign(left_set, right_set) - monomial_product_sign(
        right_set, left_set
    )
    return coefficient, tuple(sorted(left_set ^ right_set))


def overlap_witness(
    N: int, a: int, b: int, overlap: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if not (0 <= overlap <= min(a, b)) or a + b - overlap > N:
        raise ValueError((N, a, b, overlap))
    left = tuple(range(a))
    right = tuple(range(overlap)) + tuple(range(a, a + b - overlap))
    assert len(set(left) & set(right)) == overlap
    return left, right


def orbit_swap_witness(
    N: int, seed: Iterable[int], target: Iterable[int]
) -> dict[str, object]:
    """Give grade-2 brackets that walk between two coordinate k-monomials."""
    current = set(seed)
    wanted = set(target)
    if len(current) != len(wanted) or not current <= set(range(N)) or not wanted <= set(range(N)):
        raise ValueError((N, current, wanted))
    steps: list[dict[str, object]] = []
    while current != wanted:
        removed = min(current - wanted)
        added = min(wanted - current)
        bilinear = {removed, added}
        coefficient, output = clifford_bracket(bilinear, current)
        next_set = (current - {removed}) | {added}
        if coefficient == 0 or set(output) != next_set:
            raise AssertionError("invalid Johnson-graph orbit step")
        steps.append(
            {
                "bilinear_support": sorted(bilinear),
                "input_support": sorted(current),
                "bracket_coefficient": coefficient,
                "output_support": list(output),
            }
        )
        current = next_set
    return {
        "N": N,
        "grade": len(wanted),
        "seed_support": sorted(seed),
        "target_support": sorted(target),
        "steps": steps,
        "passed": current == wanted,
    }


def corrected_generation_witness(N: int, k: int) -> dict[str, object]:
    """Witness grade k -> grade 6, then grade 6 -> every higher target grade."""
    if k < 6 or k % 4 != 2 or k > N - 4:
        raise ValueError((N, k))
    first_left, first_right = overlap_witness(N, k, k, k - 3)
    first_coefficient, first_output = clifford_bracket(first_left, first_right)
    if first_coefficient == 0 or len(first_output) != 6:
        raise AssertionError("the k-to-6 witness failed")
    raises: list[dict[str, object]] = []
    for grade in range(6, N, 4):
        target = grade + 4
        if target not in target_noncentral_grades(N):
            continue
        left, right = overlap_witness(N, 6, grade, 1)
        coefficient, output = clifford_bracket(left, right)
        if coefficient == 0 or len(output) != target:
            raise AssertionError("the grade-6 raising witness failed")
        raises.append(
            {
                "from_grades": [6, grade],
                "overlap": 1,
                "union_size": len(set(left) | set(right)),
                "output_grade": len(output),
                "bracket_coefficient": coefficient,
            }
        )
    return {
        "N": N,
        "seed_grade": k,
        "k_to_6": {
            "overlap": k - 3,
            "union_size": len(set(first_left) | set(first_right)),
            "output_grade": len(first_output),
            "bracket_coefficient": first_coefficient,
        },
        "raising_steps": raises,
        "grade_closure": list(grade_closure(N, (2, k))),
        "target_noncentral_grades": list(target_noncentral_grades(N)),
    }


def run_generation_audit() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    witnesses: list[dict[str, object]] = []

    # Exhaust every admissible corrected seed over a range far beyond the graph cases.
    for N in range(8, 82, 2):
        target = target_noncentral_grades(N)
        for k in range(6, N - 3, 4):
            closure = grade_closure(N, (2, k))
            passed = closure == target
            checks.append(
                {
                    "name": f"corrected generation N={N} k={k}",
                    "passed": passed,
                    "detail": f"closure={closure}, target={target}",
                }
            )
            if N in {10, 12, 14, 16, 24} and k in {6, 10, 14, 18}:
                witnesses.append(corrected_generation_witness(N, k))

    # The endpoint grade is the load-bearing failure in the proposed lemma.
    endpoint_failures: list[dict[str, object]] = []
    for N in range(12, 65, 4):
        k = N - 2
        closure = grade_closure(N, (2, k))
        expected_exception = (2, k)
        target = target_noncentral_grades(N)
        row = {
            "N": N,
            "seed_grade": k,
            "closure": list(closure),
            "full_target": list(target),
            "proper_for_N_at_least_12": closure != target,
        }
        endpoint_failures.append(row)
        checks.append(
            {
                "name": f"endpoint seed exception N={N}",
                "passed": closure == expected_exception and closure != target,
                "detail": f"grade closure {closure}, not {target}",
            }
        )

    # N=8 is the accidental equality: {2,N-2} already is the whole class.
    small_cycle = grade_closure(8, (2, 6))
    checks.append(
        {
            "name": "C4 accidental full grade class",
            "passed": small_cycle == target_noncentral_grades(8) == (2, 6),
            "detail": f"closure={small_cycle}",
        }
    )

    # The literal k<=N proposal also included a central top-grade endpoint when N=2 mod 4.
    top_failure = grade_closure(14, (2, 14))
    checks.append(
        {
            "name": "top-grade seed is central",
            "passed": top_failure == (2, 14),
            "detail": f"N=14, closure={top_failure}",
        }
    )

    # Coordinate orbit certificate, including the split middle exterior power.
    middle_orbit = orbit_swap_witness(12, range(6), range(6, 12))
    checks.append(
        {
            "name": "middle-grade Johnson orbit crosses to Hodge complement",
            "passed": bool(middle_orbit["passed"]) and len(middle_orbit["steps"]) == 6,
            "detail": "six exact grade-2 swaps take I to its disjoint complement",
        }
    )
    general_orbits = [
        orbit_swap_witness(18, (1, 2, 3, 4, 5, 6), (0, 3, 7, 10, 14, 17)),
        orbit_swap_witness(24, range(10), range(14, 24)),
    ]
    checks.append(
        {
            "name": "non-middle coordinate orbit samples",
            "passed": all(bool(row["passed"]) for row in general_orbits),
            "detail": f"step counts={[len(row['steps']) for row in general_orbits]}",
        }
    )

    # Exhaust the support formula itself for N<=12.
    support_formula_passed = True
    support_cases = 0
    for N in range(2, 11):
        universe = range(N)
        for a in range(N + 1):
            predicted_by_b: dict[int, set[int]] = {}
            for b in range(N + 1):
                predicted_by_b[b] = set(bracket_grade_outputs(N, a, b))
            for left in combinations(universe, a):
                for b in range(N + 1):
                    actual: set[int] = set()
                    for right in combinations(universe, b):
                        coefficient, output = clifford_bracket(left, right)
                        if coefficient:
                            actual.add(len(output))
                    support_cases += 1
                    if actual != predicted_by_b[b]:
                        support_formula_passed = False
                        break
                if not support_formula_passed:
                    break
            if not support_formula_passed:
                break
        if not support_formula_passed:
            break
        budget_tick(started, f"support formula N={N}")
    checks.append(
        {
            "name": "commutator overlap formula exhaustive",
            "passed": support_formula_passed,
            "detail": f"{support_cases} fixed-left/grade cases through N=10",
        }
    )

    if not all(bool(row["passed"]) for row in checks):
        raise AssertionError("grade-generation audit failed")
    return {
        "tag": "[COMPUTATION]",
        "original_lemma_verdict": "REFUTED at k=N-2 (and at a central top grade if read literally)",
        "corrected_lemma": (
            "grade 2 plus one full coordinate orbit of grade k generates every noncentral "
            "grade 2 mod 4 iff 6<=k<=N-4; the graph-relevant exception k=N-2 closes on {2,N-2}"
        ),
        "witnesses": witnesses,
        "endpoint_failures": endpoint_failures,
        "middle_orbit": middle_orbit,
        "general_orbits": general_orbits,
        "checks": checks,
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }


def main() -> int:
    result = run_generation_audit()
    print(
        json.dumps(
            {
                "checks": len(result["checks"]),
                "endpoint_failures": len(result["endpoint_failures"]),
                "process_time_seconds": result["process_time_seconds"],
            }
        )
    )
    print("PASS e179 grade generation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
