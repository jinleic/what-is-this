#!/usr/bin/env python3
"""Exact finite anchors for the rank-one extra-row restriction obstruction.

The all-width obstruction is a symbolic Pauli/linear-algebra argument recorded in
proofs/rectangle_twogen.md.  This producer only supplies bounded integer and
rational anchors; it never promotes their finite range to an induction.
"""

from __future__ import annotations

import json
from fractions import Fraction


def rational_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    matrix = [[Fraction(value) for value in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot_value for value in matrix[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def walsh_columns(width: int) -> list[list[int]]:
    """Evaluation matrix of I,Z_0,...,Z_(width-1) on the Z basis."""
    return [
        [1] + [1 if ((state >> site) & 1) == 0 else -1 for site in range(width)]
        for state in range(1 << width)
    ]


def transverse_field_image(state: int, width: int) -> dict[int, int]:
    image: dict[int, int] = {}
    for site in range(width):
        flipped = state ^ (1 << site)
        image[flipped] = image.get(flipped, 0) + 1
    return image


def _record(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append(
        {
            "name": name,
            "passed": bool(passed),
            "detail": f"[COMPUTATION] {detail}",
        }
    )


def audit_width(width: int) -> tuple[dict[str, object], list[dict[str, object]]]:
    if width < 1:
        raise ValueError("width must be positive")

    schmidt_rank = rational_rank(walsh_columns(width))
    images = [transverse_field_image(state, width) for state in range(1 << width)]
    support_sizes = [len(image) for image in images]
    squared_norms = [sum(coefficient * coefficient for coefficient in image.values()) for image in images]
    original_overlaps = [image.get(state, 0) for state, image in enumerate(images)]

    seam_scalar = width
    removed_horizontal_scalar = 2 * (width - 1)
    anomaly_scalar = seam_scalar + removed_horizontal_scalar

    checks: list[dict[str, object]] = []
    _record(
        checks,
        f"width_{width}_core_schmidt_rank",
        schmidt_rank == width + 1,
        f"rank_Q(I,Z_0,...,Z_{width - 1})={schmidt_rank}",
    )
    _record(
        checks,
        f"width_{width}_field_support",
        min(support_sizes) == max(support_sizes) == width,
        "sum_j X_j sends every joint-Z eigenstate to width distinct basis states",
    )
    _record(
        checks,
        f"width_{width}_field_norm",
        min(squared_norms) == max(squared_norms) == width,
        "the exact squared norm of the transverse-field image is width",
    )
    _record(
        checks,
        f"width_{width}_field_overlap",
        all(overlap == 0 for overlap in original_overlaps),
        "the transverse-field image is orthogonal to its joint-Z eigenstate",
    )
    _record(
        checks,
        f"width_{width}_scalar_accounting",
        anomaly_scalar == 3 * width - 2,
        "width seam self-terms plus two endpoints for each removed-row horizontal edge give 3w-2",
    )

    data = {
        "tag": "[COMPUTATION]",
        "width": width,
        "joint_z_eigenline_count": 1 << width,
        "core_operator_evaluation_rank_over_Q": schmidt_rank,
        "transverse_field_support_size_range": [min(support_sizes), max(support_sizes)],
        "transverse_field_squared_norm_range": [min(squared_norms), max(squared_norms)],
        "transverse_field_original_overlap_values": sorted(set(original_overlaps)),
        "depth_three_scalar_accounting": {
            "seam_self_terms": seam_scalar,
            "removed_horizontal_endpoint_terms": removed_horizontal_scalar,
            "total": anomaly_scalar,
        },
    }
    return data, checks


def run_induction_audit(
    widths: tuple[int, ...] = (1, 2, 3, 4, 5, 6),
) -> tuple[dict[str, object], list[dict[str, object]]]:
    anchors: list[dict[str, object]] = []
    checks: list[dict[str, object]] = []
    for width in widths:
        anchor, width_checks = audit_width(width)
        anchors.append(anchor)
        checks.extend(width_checks)

    data = {
        "rank_one_tensor_restriction": {
            "tag": "[THEOREM]",
            "statement": (
                "For every positive width, no nonzero extra-row vector Omega makes "
                "H_core tensor span{Omega} invariant under both rectangle generators."
            ),
            "proof_kernel": (
                "B-invariance and independence of I,Z_s0,... force every Z_j Omega into "
                "span{Omega}; the common joint-Z eigenlines are computational basis lines, "
                "but sum_j X_j maps each such line to a nonzero orthogonal vector."
            ),
        },
        "plus_compression_counterexample": {
            "tag": "[THEOREM]",
            "statement": (
                "The |+>-row compression matches A up to width times identity, matches B, "
                "and matches {A,B}, but fails on {B,{A,B}} by the noncentral boundary field "
                "plus (3 width-2) times identity."
            ),
        },
        "finite_anchors": anchors,
        "finite_scope": {
            "tag": "[COMPUTATION]",
            "statement": (
                "The enumeration stops at width 6; it checks the ingredients but is not the "
                "source of either all-width theorem."
            ),
        },
        "remaining_routes": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "Higher-rank ancilla codes, non-tensor representation embeddings, and direct "
                "rectangle Lie-word inductions are not excluded."
            ),
        },
    }
    return data, checks


def main() -> int:
    data, checks = run_induction_audit()
    print(json.dumps(data, indent=2, sort_keys=True))
    if all(bool(item["passed"]) for item in checks):
        print("PASS e198 rectangle induction obstruction")
        return 0
    for item in checks:
        if not item["passed"]:
            print(f"FAIL {item['name']}: {item['detail']}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
