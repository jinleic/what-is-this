#!/usr/bin/env python3
"""Clean-room verifier for the rectangle two-generator route closure.

This file intentionally imports none of e197, e198, or e199.  Its Pauli
implementation uses symplectic bit masks rather than producer word tuples.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "rectangle_twogen.json"

Operator = dict[int, int]


def require(condition: bool, name: str, detail: str) -> None:
    if not condition:
        raise AssertionError(f"{name}: {detail}")
    print(f"[PASS] {name}: {detail}", flush=True)


def accumulate(operator: Operator, code: int, coefficient: int) -> None:
    value = operator.get(code, 0) + coefficient
    if value:
        operator[code] = value
    else:
        operator.pop(code, None)


def add(*operators: Operator) -> Operator:
    result: Operator = {}
    for operator in operators:
        for code, coefficient in operator.items():
            accumulate(result, code, coefficient)
    return result


def scale(operator: Operator, scalar: int) -> Operator:
    return {code: scalar * coefficient for code, coefficient in operator.items() if scalar}


def unpack(code: int, sites: int) -> tuple[int, int]:
    mask = (1 << sites) - 1
    return code & mask, code >> sites


def pack(x_mask: int, z_mask: int, sites: int) -> int:
    return x_mask | (z_mask << sites)


def multiplication_phase(left: int, right: int, sites: int) -> int:
    """Exponent q with P_left P_right=i^q P_(left xor right)."""
    left_x, left_z = unpack(left, sites)
    right_x, right_z = unpack(right, sites)
    out_x = left_x ^ right_x
    out_z = left_z ^ right_z
    exponent = (
        (left_x & left_z).bit_count()
        + (right_x & right_z).bit_count()
        - (out_x & out_z).bit_count()
        + 2 * (left_z & right_x).bit_count()
    )
    return exponent % 4


def bracket(left: Operator, right: Operator, sites: int) -> Operator:
    """Exact [left,right]/(2i) in the canonical Hermitian Pauli basis."""
    result: Operator = {}
    for left_code, left_coefficient in left.items():
        for right_code, right_coefficient in right.items():
            phase = multiplication_phase(left_code, right_code, sites)
            if phase == 1:
                sign = 1
            elif phase == 3:
                sign = -1
            else:
                continue
            accumulate(
                result,
                left_code ^ right_code,
                sign * left_coefficient * right_coefficient,
            )
    return result


def edges(rows: int, width: int) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for row in range(rows):
        for column in range(width):
            site = row * width + column
            if column + 1 < width:
                result.append((site, site + 1))
            if row + 1 < rows:
                result.append((site, site + width))
    return result


def generators(rows: int, width: int) -> tuple[Operator, Operator]:
    sites = rows * width
    field: Operator = {}
    bonds: Operator = {}
    for site in range(sites):
        accumulate(field, pack(1 << site, 0, sites), 1)
    for left, right in edges(rows, width):
        accumulate(bonds, pack(0, (1 << left) | (1 << right), sites), 1)
    return field, bonds


def identity(sites: int, coefficient: int = 1) -> Operator:
    return {pack(0, 0, sites): coefficient} if coefficient else {}


def last_row_field(rows: int, width: int) -> Operator:
    sites = rows * width
    result: Operator = {}
    offset = (rows - 1) * width
    for column in range(width):
        accumulate(result, pack(1 << (offset + column), 0, sites), 1)
    return result


def compress_plus(operator: Operator, full_sites: int, core_sites: int) -> Operator:
    core_mask = (1 << core_sites) - 1
    extra_mask = ((1 << full_sites) - 1) ^ core_mask
    result: Operator = {}
    for code, coefficient in operator.items():
        x_mask, z_mask = unpack(code, full_sites)
        if z_mask & extra_mask:
            continue
        compressed = pack(x_mask & core_mask, z_mask & core_mask, core_sites)
        accumulate(result, compressed, coefficient)
    return result


def reflect_rows(operator: Operator, rows: int, width: int) -> Operator:
    sites = rows * width
    result: Operator = {}
    for code, coefficient in operator.items():
        x_mask, z_mask = unpack(code, sites)
        reflected_x = 0
        reflected_z = 0
        for row in range(rows):
            for column in range(width):
                source = row * width + column
                target = (rows - 1 - row) * width + column
                reflected_x |= ((x_mask >> source) & 1) << target
                reflected_z |= ((z_mask >> source) & 1) << target
        accumulate(result, pack(reflected_x, reflected_z, sites), coefficient)
    return result


def pauli_string(code: int, sites: int) -> str:
    x_mask, z_mask = unpack(code, sites)
    letters: list[str] = []
    for site in range(sites):
        key = ((x_mask >> site) & 1, (z_mask >> site) & 1)
        letters.append({(0, 0): "I", (1, 0): "X", (1, 1): "Y", (0, 1): "Z"}[key])
    return "".join(letters)


def encoded(operator: Operator, sites: int) -> list[dict[str, object]]:
    return [
        {"pauli": pauli_string(code, sites), "coefficient": coefficient}
        for code, coefficient in sorted(
            operator.items(), key=lambda item: pauli_string(item[0], sites)
        )
    ]


def digest(operator: Operator, sites: int) -> str:
    payload = json.dumps(encoded(operator, sites), separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def fraction_rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        divisor = work[rank][column]
        for entry in range(column, columns):
            work[rank][entry] /= divisor
        for row in range(rank + 1, rows):
            multiplier = work[row][column]
            if not multiplier:
                continue
            for entry in range(column, columns):
                work[row][entry] -= multiplier * work[rank][entry]
        rank += 1
        if rank == rows:
            break
    return rank


def verify_anchor(anchor: dict[str, object]) -> None:
    core_rows, width = (int(value) for value in anchor["core_shape"])
    full_rows = core_rows + 1
    core_sites = core_rows * width
    full_sites = full_rows * width
    core_a, core_b = generators(core_rows, width)
    full_a, full_b = generators(full_rows, width)

    compressed_a = compress_plus(full_a, full_sites, core_sites)
    compressed_b = compress_plus(full_b, full_sites, core_sites)
    core_first = bracket(core_a, core_b, core_sites)
    full_first = bracket(full_a, full_b, full_sites)
    compressed_first = compress_plus(full_first, full_sites, core_sites)
    core_depth_three = bracket(core_b, core_first, core_sites)
    full_depth_three = bracket(full_b, full_first, full_sites)
    compressed_depth_three = compress_plus(full_depth_three, full_sites, core_sites)
    defect = add(compressed_depth_three, scale(core_depth_three, -1))
    expected = add(last_row_field(core_rows, width), identity(core_sites, 3 * width - 2))
    reflected = reflect_rows(compressed_depth_three, core_rows, width)
    residual = add(compressed_depth_three, scale(reflected, -1))

    stem = f"clean_r{core_rows}_w{width}"
    require(
        compressed_a == add(core_a, identity(core_sites, width)),
        f"{stem}_A",
        "generator compression including the scalar shift",
    )
    require(compressed_b == core_b, f"{stem}_B", "bond generator compression")
    require(compressed_first == core_first, f"{stem}_first", "first-bracket compression")
    require(defect == expected, f"{stem}_defect", "all coefficients in the depth-three anomaly")
    require(residual, f"{stem}_reflection", "nonzero reflection-odd residual")
    require(encoded(defect, core_sites) == anchor["defect"], f"{stem}_stored_vector", "stored exact defect")
    require(digest(defect, core_sites) == anchor["defect_sha256"], f"{stem}_defect_hash", "defect digest")
    require(
        digest(residual, core_sites) == anchor["reflection_residual_sha256"],
        f"{stem}_reflection_hash",
        "reflection-residual digest",
    )


def verify_induction_anchor(anchor: dict[str, object]) -> None:
    width = int(anchor["width"])
    walsh = [
        [1] + [1 if ((state >> site) & 1) == 0 else -1 for site in range(width)]
        for state in range(1 << width)
    ]
    rank = fraction_rank(walsh)
    supports: list[int] = []
    norms: list[int] = []
    overlaps: set[int] = set()
    for state in range(1 << width):
        image: dict[int, int] = {}
        for site in range(width):
            target = state ^ (1 << site)
            image[target] = image.get(target, 0) + 1
        supports.append(len(image))
        norms.append(sum(value * value for value in image.values()))
        overlaps.add(image.get(state, 0))

    stem = f"clean_width_{width}"
    require(rank == width + 1, f"{stem}_rank", "exact rank of I and the seam Z operators")
    require(min(supports) == max(supports) == width, f"{stem}_support", "transverse image support")
    require(min(norms) == max(norms) == width, f"{stem}_norm", "transverse image squared norm")
    require(overlaps == {0}, f"{stem}_overlap", "orthogonality to each joint-Z eigenstate")
    require(
        anchor["core_operator_evaluation_rank_over_Q"] == rank,
        f"{stem}_stored_rank",
        "stored exact rational rank",
    )
    require(
        anchor["depth_three_scalar_accounting"]["total"] == 3 * width - 2,
        f"{stem}_scalar",
        "stored scalar anomaly accounting",
    )


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    require(set(artifact) == {"meta", "data", "checks"}, "artifact_shape", "top-level meta/data/checks")
    require(artifact["checks"], "stored_checks_nonempty", "producer stored computed checks")
    require(
        all(bool(item["passed"]) for item in artifact["checks"]),
        "stored_checks_pass",
        f"all {len(artifact['checks'])} producer checks passed",
    )

    for relative, expected in artifact["meta"]["source_sha256"].items():
        observed = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        require(observed == expected, f"source_hash_{relative}", "source digest matches artifact")

    word_data = artifact["data"]["word_compression"]
    require(word_data["finite_anchors"], "word_anchors_nonempty", "bounded word anchors are present")
    for anchor in word_data["finite_anchors"]:
        verify_anchor(anchor)

    induction_data = artifact["data"]["rank_one_restriction"]
    require(induction_data["finite_anchors"], "restriction_anchors_nonempty", "bounded restriction anchors are present")
    for anchor in induction_data["finite_anchors"]:
        verify_induction_anchor(anchor)

    conclusion = artifact["data"]["conclusion"]
    require(conclusion["tag"] == "[THEOREM]", "conclusion_tag", "route closure is theorem-tagged")
    require(
        artifact["data"]["dimension_status"]["tag"] == "[UNRESOLVED]",
        "dimension_scope_tag",
        "no all-rectangle dimension theorem is promoted",
    )
    require(
        artifact["meta"]["modular_arithmetic"] == "not used",
        "arithmetic_scope",
        "the certificate uses exact integer/rational arithmetic only",
    )
    print("PASS test_rectangle_twogen", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FAIL test_rectangle_twogen: {exc}", file=sys.stderr, flush=True)
        raise
