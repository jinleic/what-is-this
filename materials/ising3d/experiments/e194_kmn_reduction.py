#!/usr/bin/env python3
"""[THEOREM] Exact collective-spin reduction and structural upper bounds for K_{m,n}."""

from __future__ import annotations

import math
import platform
import resource
import time
from dataclasses import dataclass
from math import comb
from typing import Iterable

import numpy as np

CPU_BUDGET_SECONDS = 300.0
RSS_CAP_BYTES = 1_900_000_000
GOOD_PRIMES = (1_000_003, 1_000_033)


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str) -> None:
    used = time.process_time() - started
    rss = max_rss_bytes()
    if used > CPU_BUDGET_SECONDS:
        raise TimeoutError(
            f"process-time budget exceeded at {stage}: {used:.3f}s > "
            f"{CPU_BUDGET_SECONDS:.3f}s"
        )
    if rss >= RSS_CAP_BYTES:
        raise MemoryError(
            f"RSS cap reached at {stage}: {rss} >= {RSS_CAP_BYTES}"
        )


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def spin_labels(part_size: int, *, positive_only: bool = False) -> tuple[int, ...]:
    """Return twice-spin labels r=m mod 2,m mod 2+2,...,m."""
    if part_size < 1:
        raise ValueError("part size must be positive")
    labels = tuple(range(part_size % 2, part_size + 1, 2))
    return tuple(label for label in labels if label or not positive_only)


def spin_matrix_square_sum(part_size: int, *, positive_only: bool = False) -> int:
    return sum((label + 1) ** 2 for label in spin_labels(part_size, positive_only=positive_only))


def schur_commutant_dimension(left_size: int, right_size: int) -> int:
    return comb(left_size + 3, 3) * comb(right_size + 3, 3)


def representation_blocks(left_size: int, right_size: int) -> tuple[tuple[int, int], ...]:
    """Multiplicity-free spin blocks; balanced transpose-equivalent blocks occur once."""
    left = spin_labels(left_size)
    right = spin_labels(right_size)
    if left_size != right_size:
        return tuple((r, s) for r in left for s in right if r or s)
    return tuple(
        (r, s)
        for index, r in enumerate(left)
        for s in right[index:]
        if r or s
    )


@dataclass(frozen=True)
class SpinMatrices:
    label: int
    x: np.ndarray
    h: np.ndarray
    reversal: np.ndarray
    sign: np.ndarray
    metric: np.ndarray


def spin_matrices(label: int) -> SpinMatrices:
    """Integral Sym^r basis for sum X, sum Z and its exact contravariant metric."""
    if label < 0:
        raise ValueError("twice-spin label must be nonnegative")
    dimension = label + 1
    x = np.zeros((dimension, dimension), dtype=np.int64)
    for column in range(dimension):
        if column:
            x[column - 1, column] = column
        if column < label:
            x[column + 1, column] = label - column
    h = np.diag([label - 2 * column for column in range(dimension)]).astype(np.int64)
    reversal = np.fliplr(np.eye(dimension, dtype=np.int64))
    sign = np.diag([1 if column % 2 == 0 else -1 for column in range(dimension)]).astype(
        np.int64
    )
    binomials = [comb(label, column) for column in range(dimension)]
    scale = math.lcm(*binomials)
    metric = np.diag([scale // value for value in binomials]).astype(np.int64)
    return SpinMatrices(label, x, h, reversal, sign, metric)


@dataclass(frozen=True)
class SectorMatrices:
    left_label: int
    right_label: int
    a: np.ndarray
    b: np.ndarray
    parity: np.ndarray
    charge: np.ndarray
    form: np.ndarray
    swap: np.ndarray | None


def tensor_swap(dimension: int) -> np.ndarray:
    swap = np.zeros((dimension * dimension, dimension * dimension), dtype=np.int64)
    for left in range(dimension):
        for right in range(dimension):
            swap[right * dimension + left, left * dimension + right] = 1
    return swap


def sector_matrices(left_label: int, right_label: int) -> SectorMatrices:
    left = spin_matrices(left_label)
    right = spin_matrices(right_label)
    identity_left = np.eye(left_label + 1, dtype=np.int64)
    identity_right = np.eye(right_label + 1, dtype=np.int64)
    a = np.kron(left.x, identity_right) + np.kron(identity_left, right.x)
    b = np.kron(left.h, right.h)
    parity = np.kron(left.reversal, right.reversal)
    charge = np.kron(left.sign @ left.reversal, right.sign)
    metric = np.kron(left.metric, right.metric)
    form = metric @ charge
    swap = tensor_swap(left_label + 1) if left_label == right_label else None
    return SectorMatrices(left_label, right_label, a, b, parity, charge, form, swap)


def sector_identity_audit(left_label: int, right_label: int) -> dict[str, object]:
    sector = sector_matrices(left_label, right_label)
    left = spin_matrices(left_label)
    right = spin_matrices(right_label)
    sign = -1 if (left_label + right_label) % 2 else 1
    checks = {
        "left_spin_metric": np.array_equal(left.x.T @ left.metric, left.metric @ left.x)
        and np.array_equal(left.h.T @ left.metric, left.metric @ left.h),
        "right_spin_metric": np.array_equal(right.x.T @ right.metric, right.metric @ right.x)
        and np.array_equal(right.h.T @ right.metric, right.metric @ right.h),
        "parity_commutes_with_generators": np.array_equal(
            sector.parity @ sector.a, sector.a @ sector.parity
        )
        and np.array_equal(sector.parity @ sector.b, sector.b @ sector.parity),
        "charge_anticommutes_with_generators": np.array_equal(
            sector.charge @ sector.a, -sector.a @ sector.charge
        )
        and np.array_equal(sector.charge @ sector.b, -sector.b @ sector.charge),
        "generators_preserve_form": np.count_nonzero(
            sector.a.T @ sector.form + sector.form @ sector.a
        )
        == 0
        and np.count_nonzero(sector.b.T @ sector.form + sector.form @ sector.b) == 0,
        "parity_charge_relation": np.array_equal(
            sector.parity @ sector.charge,
            sign * sector.charge @ sector.parity,
        ),
        "form_transpose_type": np.array_equal(
            sector.form.T,
            (1 if left_label % 2 == 0 else -1) * sector.form,
        ),
    }
    if sector.swap is not None:
        checks["swap_relations"] = (
            np.array_equal(sector.swap @ sector.a, sector.a @ sector.swap)
            and np.array_equal(sector.swap @ sector.b, sector.b @ sector.swap)
            and np.array_equal(sector.swap @ sector.parity, sector.parity @ sector.swap)
            and np.array_equal(
                sector.swap.T @ sector.form @ sector.swap,
                sector.form @ sector.parity,
            )
        )
    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "labels": [left_label, right_label],
        "matrix_dimension": int(sector.a.shape[0]),
        "checks": checks,
        "passed": all(bool(value) for value in checks.values()),
    }


def so_dimension(vector_dimension: int) -> int:
    return vector_dimension * (vector_dimension - 1) // 2


def sp_dimension_on(vector_dimension: int) -> int:
    if vector_dimension % 2:
        raise ValueError("a nondegenerate alternating form requires even dimension")
    return vector_dimension * (vector_dimension + 1) // 2


def sector_container(left_label: int, right_label: int) -> dict[str, object]:
    """Dimension of the P-, K-, and (when applicable) swap-fixed sector container."""
    if left_label <= 0 or right_label <= 0:
        raise ValueError("sector container is only defined for two positive spins")
    matrix_dimension = (left_label + 1) * (right_label + 1)
    if (left_label + right_label) % 2:
        half = matrix_dimension // 2
        upper = half * half
        return {
            "labels": [left_label, right_label],
            "kind": f"gl({half}) paired across parity",
            "matrix_dimension": matrix_dimension,
            "component_dimensions": [half, half],
            "upper_dimension": upper,
            "center_dimension": 1,
            "semisimple_dimension": upper - 1,
        }
    if left_label != right_label:
        if left_label % 2 == 0:
            plus = (matrix_dimension + 1) // 2
            minus = (matrix_dimension - 1) // 2
            upper = so_dimension(plus) + so_dimension(minus)
            kind = f"so({plus}) + so({minus})"
            components = [plus, minus]
        else:
            half = matrix_dimension // 2
            upper = 2 * sp_dimension_on(half)
            kind = f"sp(on {half}) + sp(on {half})"
            components = [half, half]
        return {
            "labels": [left_label, right_label],
            "kind": kind,
            "matrix_dimension": matrix_dimension,
            "component_dimensions": components,
            "upper_dimension": upper,
            "center_dimension": 0,
            "semisimple_dimension": upper,
        }

    single_dimension = left_label + 1
    if left_label % 2 == 0:
        plus_swap_plus = ((single_dimension + 1) // 2) ** 2
        plus_swap_minus = ((single_dimension - 1) // 2) ** 2
        paired = (single_dimension * single_dimension - 1) // 4
        upper = (
            so_dimension(plus_swap_plus)
            + so_dimension(plus_swap_minus)
            + paired * paired
        )
        kind = (
            f"so({plus_swap_plus}) + so({plus_swap_minus}) + gl({paired})"
        )
    else:
        plus_swap_plus = single_dimension * (single_dimension + 2) // 4
        plus_swap_minus = single_dimension * (single_dimension - 2) // 4
        paired = single_dimension * single_dimension // 4
        upper = (
            sp_dimension_on(plus_swap_plus)
            + sp_dimension_on(plus_swap_minus)
            + paired * paired
        )
        kind = (
            f"sp(on {plus_swap_plus}) + sp(on {plus_swap_minus}) + gl({paired})"
        )
    return {
        "labels": [left_label, right_label],
        "kind": kind,
        "matrix_dimension": matrix_dimension,
        "component_dimensions": [plus_swap_plus, plus_swap_minus, paired, paired],
        "upper_dimension": upper,
        "center_dimension": 1,
        "semisimple_dimension": upper - 1,
    }


def opposite_parity_family_upper(left_size: int, right_size: int) -> dict[str, object]:
    if (left_size - right_size) % 2 == 0:
        raise ValueError("part sizes must have opposite parity")
    sectors = [
        sector_container(left_label, right_label)
        for left_label in spin_labels(left_size, positive_only=True)
        for right_label in spin_labels(right_size, positive_only=True)
    ]
    semisimple = sum(int(row["semisimple_dimension"]) for row in sectors)
    positive_left_count = len(spin_labels(left_size, positive_only=True))
    positive_right_count = len(spin_labels(right_size, positive_only=True))
    positive_left_square_sum = spin_matrix_square_sum(left_size, positive_only=True)
    positive_right_square_sum = spin_matrix_square_sum(right_size, positive_only=True)
    closed_form = (
        1
        + positive_left_square_sum * positive_right_square_sum // 4
        - positive_left_count * positive_right_count
    )
    upper = semisimple + 1
    if upper != closed_form:
        raise AssertionError((upper, closed_form))
    return {
        "tag": "[THEOREM]",
        "part_sizes": [left_size, right_size],
        "positive_spin_counts": [positive_left_count, positive_right_count],
        "positive_spin_square_sums": [
            positive_left_square_sum,
            positive_right_square_sum,
        ],
        "sector_containers": sectors,
        "semisimple_upper": semisimple,
        "global_center_upper": 1,
        "upper_dimension": upper,
        "closed_form": (
            "1 + S_m*S_n/4 - a_m*a_n, where "
            "S_t=C(t+3,3)-1_{t even} and a_t=ceil(t/2)"
        ),
        "closed_form_value": closed_form,
    }


def balanced_family_upper(part_size: int) -> dict[str, object]:
    positive = spin_labels(part_size, positive_only=True)
    sectors = [
        sector_container(left_label, right_label)
        for index, left_label in enumerate(positive)
        for right_label in positive[index:]
    ]
    semisimple = sum(int(row["semisimple_dimension"]) for row in sectors)
    center_upper = 2 if part_size % 2 == 0 else 1
    upper = semisimple + center_upper
    return {
        "tag": "[THEOREM]",
        "part_sizes": [part_size, part_size],
        "sector_containers": sectors,
        "semisimple_upper": semisimple,
        "global_center_upper": center_upper,
        "upper_dimension": upper,
        "center_explanation": (
            "odd parts: only the B center signature; even parts: that signature plus "
            "the common A-only zero-spin signature"
        ),
    }


def structural_upper(left_size: int, right_size: int) -> dict[str, object]:
    if left_size == right_size:
        return balanced_family_upper(left_size)
    if (left_size - right_size) % 2:
        return opposite_parity_family_upper(left_size, right_size)
    raise ValueError("implemented family bounds require balanced or opposite-parity parts")


def local_term_branching_dimension(left_size: int, right_size: int) -> int:
    if min(left_size, right_size) < 2 or max(left_size, right_size) < 3:
        raise ValueError("the cited local-term branching formula requires maximum degree >= 3")
    vertices = left_size + right_size
    if vertices % 2:
        return (1 << (2 * vertices - 2)) - 1
    return (1 << (2 * vertices - 2)) - ((-1) ** (vertices // 2)) * (
        1 << (vertices - 1)
    )


def hamiltonian_complete_bipartite(left_size: int, right_size: int) -> bool:
    return abs(left_size - right_size) <= 1


def reduction_audit(part_sizes: Iterable[int] = range(2, 6)) -> dict[str, object]:
    started = time.process_time()
    sizes = tuple(part_sizes)
    square_sum_rows = []
    for part_size in sizes:
        direct = spin_matrix_square_sum(part_size)
        closed = comb(part_size + 3, 3)
        square_sum_rows.append(
            {
                "part_size": part_size,
                "spin_labels": list(spin_labels(part_size)),
                "direct_square_sum": direct,
                "binomial_value": closed,
                "passed": direct == closed,
            }
        )
    sector_pairs = sorted(
        {
            block
            for left_size, right_size in ((2, 3), (3, 3), (3, 4), (4, 4), (4, 5))
            for block in representation_blocks(left_size, right_size)
        }
    )
    identity_rows = [sector_identity_audit(*block) for block in sector_pairs]
    budget_tick(started, "collective-spin reduction audit")
    return {
        "tag": "[COMPUTATION]",
        "square_sum_rows": square_sum_rows,
        "sector_identity_rows": identity_rows,
        "all_square_sums_pass": all(row["passed"] for row in square_sum_rows),
        "all_sector_identities_pass": all(row["passed"] for row in identity_rows),
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }


def main() -> int:
    audit = reduction_audit()
    passed = bool(audit["all_square_sums_pass"] and audit["all_sector_identities_pass"])
    if not passed:
        raise AssertionError(audit)
    print(
        f"PASS e194: {len(audit['sector_identity_rows'])} sectors, "
        f"cpu={audit['process_time_seconds']}s, rss={audit['peak_rss_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
