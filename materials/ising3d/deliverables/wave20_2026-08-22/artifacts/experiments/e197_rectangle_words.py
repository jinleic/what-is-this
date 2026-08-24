#!/usr/bin/env python3
"""Exact Pauli-word audit of boundary-row compression for open rectangles.

The Hermitian bracket is {S,T} = [S,T]/(2i).  It keeps integral Hermitian
Pauli sums integral and differs from the anti-Hermitian Lie bracket only by a
nonzero scalar.  No producer numbered e191 or later is imported.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable

Word = tuple[str, ...]
Operator = dict[Word, int]

# phase is the exponent of i in left * right = i**phase * product.
_LOCAL_PRODUCT: dict[tuple[str, str], tuple[int, str]] = {
    ("I", "I"): (0, "I"),
    ("I", "X"): (0, "X"),
    ("I", "Y"): (0, "Y"),
    ("I", "Z"): (0, "Z"),
    ("X", "I"): (0, "X"),
    ("Y", "I"): (0, "Y"),
    ("Z", "I"): (0, "Z"),
    ("X", "X"): (0, "I"),
    ("Y", "Y"): (0, "I"),
    ("Z", "Z"): (0, "I"),
    ("X", "Y"): (1, "Z"),
    ("Y", "X"): (3, "Z"),
    ("Y", "Z"): (1, "X"),
    ("Z", "Y"): (3, "X"),
    ("Z", "X"): (1, "Y"),
    ("X", "Z"): (3, "Y"),
}


def _accumulate(target: Operator, word: Word, coefficient: int) -> None:
    if not coefficient:
        return
    value = target.get(word, 0) + coefficient
    if value:
        target[word] = value
    else:
        target.pop(word, None)


def add_operators(*operators: Operator) -> Operator:
    result: Operator = {}
    for operator in operators:
        for word, coefficient in operator.items():
            _accumulate(result, word, coefficient)
    return result


def scale_operator(operator: Operator, scalar: int) -> Operator:
    return {word: scalar * coefficient for word, coefficient in operator.items() if scalar}


def single_word(sites: int, entries: dict[int, str]) -> Word:
    letters = ["I"] * sites
    for site, letter in entries.items():
        letters[site] = letter
    return tuple(letters)


def _product_phase(left: Word, right: Word) -> tuple[int, Word]:
    phase = 0
    product: list[str] = []
    for left_letter, right_letter in zip(left, right, strict=True):
        local_phase, local_product = _LOCAL_PRODUCT[(left_letter, right_letter)]
        phase = (phase + local_phase) % 4
        product.append(local_product)
    return phase, tuple(product)


def hermitian_bracket(left: Operator, right: Operator) -> Operator:
    """Return [left,right]/(2i) in the Hermitian Pauli basis over Z."""
    result: Operator = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            phase, product = _product_phase(left_word, right_word)
            if phase == 1:
                sign = 1
            elif phase == 3:
                sign = -1
            else:
                # Products with real phase commute, so their bracket vanishes.
                continue
            _accumulate(result, product, sign * left_coefficient * right_coefficient)
    return result


def grid_edges(rows: int, width: int) -> list[tuple[int, int]]:
    edges: list[tuple[int, int]] = []
    for row in range(rows):
        for column in range(width):
            site = row * width + column
            if column + 1 < width:
                edges.append((site, site + 1))
            if row + 1 < rows:
                edges.append((site, site + width))
    return edges


def grid_generators(rows: int, width: int) -> tuple[Operator, Operator]:
    sites = rows * width
    field: Operator = {}
    bonds: Operator = {}
    for site in range(sites):
        _accumulate(field, single_word(sites, {site: "X"}), 1)
    for left, right in grid_edges(rows, width):
        _accumulate(bonds, single_word(sites, {left: "Z", right: "Z"}), 1)
    return field, bonds


def identity_operator(sites: int, coefficient: int = 1) -> Operator:
    return {single_word(sites, {}): coefficient} if coefficient else {}


def boundary_field(core_rows: int, width: int) -> Operator:
    sites = core_rows * width
    result: Operator = {}
    first = (core_rows - 1) * width
    for column in range(width):
        _accumulate(result, single_word(sites, {first + column: "X"}), 1)
    return result


def compress_plus_row(operator: Operator, core_sites: int) -> Operator:
    """Compute <+|^(extra) operator |+>^(extra), removing the final row."""
    result: Operator = {}
    for word, coefficient in operator.items():
        if all(letter in ("I", "X") for letter in word[core_sites:]):
            _accumulate(result, word[:core_sites], coefficient)
    return result


def reflect_core_rows(operator: Operator, rows: int, width: int) -> Operator:
    result: Operator = {}
    for word, coefficient in operator.items():
        reflected = ["I"] * len(word)
        for row in range(rows):
            for column in range(width):
                source = row * width + column
                target = (rows - 1 - row) * width + column
                reflected[target] = word[source]
        _accumulate(result, tuple(reflected), coefficient)
    return result


def encoded_operator(operator: Operator) -> list[dict[str, object]]:
    return [
        {"pauli": "".join(word), "coefficient": coefficient}
        for word, coefficient in sorted(operator.items())
    ]


def operator_sha256(operator: Operator) -> str:
    payload = json.dumps(encoded_operator(operator), separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()


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


def audit_case(core_rows: int, width: int) -> tuple[dict[str, object], list[dict[str, object]]]:
    if core_rows < 2 or width < 1:
        raise ValueError("the reflection obstruction requires core_rows >= 2 and width >= 1")

    core_sites = core_rows * width
    full_rows = core_rows + 1
    core_a, core_b = grid_generators(core_rows, width)
    full_a, full_b = grid_generators(full_rows, width)

    compressed_a = compress_plus_row(full_a, core_sites)
    compressed_b = compress_plus_row(full_b, core_sites)
    core_first = hermitian_bracket(core_a, core_b)
    compressed_first = compress_plus_row(hermitian_bracket(full_a, full_b), core_sites)
    core_depth_three = hermitian_bracket(core_b, core_first)
    compressed_depth_three = compress_plus_row(
        hermitian_bracket(full_b, hermitian_bracket(full_a, full_b)), core_sites
    )
    defect = add_operators(compressed_depth_three, scale_operator(core_depth_three, -1))
    expected_defect = add_operators(
        boundary_field(core_rows, width), identity_operator(core_sites, 3 * width - 2)
    )
    reflected_compressed = reflect_core_rows(compressed_depth_three, core_rows, width)

    expected_a = add_operators(core_a, identity_operator(core_sites, width))
    checks: list[dict[str, object]] = []
    stem = f"r{core_rows}_w{width}"
    _record(checks, f"{stem}_A_compression", compressed_a == expected_a, "Phi(A)=A_core+w I")
    _record(checks, f"{stem}_B_compression", compressed_b == core_b, "Phi(B)=B_core")
    _record(
        checks,
        f"{stem}_first_bracket",
        compressed_first == core_first,
        "Phi({A,B})={A_core,B_core}",
    )
    _record(
        checks,
        f"{stem}_depth_three_defect",
        defect == expected_defect,
        "Phi({B,{A,B}})-{B_core,{A_core,B_core}}=S_boundary+(3w-2)I",
    )
    _record(
        checks,
        f"{stem}_reflection_obstruction",
        compressed_depth_three != reflected_compressed,
        "the compressed depth-three word is not fixed by core row reflection",
    )

    case = {
        "tag": "[COMPUTATION]",
        "core_shape": [core_rows, width],
        "full_shape": [full_rows, width],
        "core_sites": core_sites,
        "full_sites": full_rows * width,
        "core_edge_count": len(grid_edges(core_rows, width)),
        "full_edge_count": len(grid_edges(full_rows, width)),
        "depth_three_core_support": len(core_depth_three),
        "depth_three_compressed_support": len(compressed_depth_three),
        "defect": encoded_operator(defect),
        "defect_sha256": operator_sha256(defect),
        "reflection_residual_sha256": operator_sha256(
            add_operators(compressed_depth_three, scale_operator(reflected_compressed, -1))
        ),
    }
    return case, checks


def run_word_audit(
    cases: Iterable[tuple[int, int]] = ((2, 1), (2, 2), (2, 3), (2, 4), (3, 2), (3, 3)),
) -> tuple[dict[str, object], list[dict[str, object]]]:
    anchors: list[dict[str, object]] = []
    checks: list[dict[str, object]] = []
    for core_rows, width in cases:
        anchor, case_checks = audit_case(core_rows, width)
        anchors.append(anchor)
        checks.extend(case_checks)
    data = {
        "tag": "[COMPUTATION]",
        "bracket_convention": "{S,T}=[S,T]/(2i) in the Hermitian Pauli basis",
        "compression_state": "the removed final row is |+>^width",
        "finite_anchors": anchors,
        "anchor_scope": "core rows 2,3 and widths at most 4; the all-size statement is proved symbolically, not inferred from these anchors",
    }
    return data, checks


def main() -> int:
    data, checks = run_word_audit()
    print(json.dumps(data, indent=2, sort_keys=True))
    if all(bool(item["passed"]) for item in checks):
        print("PASS e197 rectangle word compression")
        return 0
    for item in checks:
        if not item["passed"]:
            print(f"FAIL {item['name']}: {item['detail']}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
