"""Standalone exact checks for the finite series-structure certificates."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from ising.interlayer import anisotropic_box_even_subgraph, anisotropic_flm_c2_series
from ising.lattices import square

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "series" / "structure_certificates.json"
HT = ROOT / "results" / "series" / "extended2_sc_ht_free_energy.json"
LT = ROOT / "results" / "series" / "extended2_sc_lt_free_energy.json"
ORDER = 8


def _check(name: str, passed: bool, detail: str) -> None:
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def _load_even(path: Path) -> tuple[Fraction, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    full = tuple(Fraction(value) for value in payload["data"]["coefficients"])
    assert all(full[degree] == 0 for degree in range(1, len(full), 2))
    return tuple(full[::2])


def _convolve(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    out = [Fraction(0) for _ in range(order + 1)]
    for i, first in enumerate(left[: order + 1]):
        if first:
            for j, second in enumerate(right[: order + 1 - i]):
                out[i + j] += first * second
    return tuple(out)


def _powers(
    series: tuple[Fraction, ...], maximum: int, order: int
) -> tuple[tuple[Fraction, ...], ...]:
    values = [(Fraction(1),) + (Fraction(0),) * order]
    for _ in range(maximum):
        values.append(_convolve(values[-1], series, order))
    return tuple(values)


def _rank(rows: Sequence[Sequence[Fraction]], columns: int) -> int:
    """Independent exact Fraction Gaussian elimination."""

    matrix = [[Fraction(value) for value in row] for row in rows]
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(pivot_row + 1, len(matrix)):
            if matrix[row][column]:
                multiple = matrix[row][column]
                matrix[row] = [
                    value - multiple * pivot_value
                    for value, pivot_value in zip(
                        matrix[row], matrix[pivot_row], strict=True
                    )
                ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def _algebraic_rows(
    series: tuple[Fraction, ...], degree_x: int, degree_f: int
) -> tuple[list[list[Fraction]], list[tuple[int, int]]]:
    count = len(series)
    power = _powers(series, degree_f, count - 1)
    monomials = [
        (a, b)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    return (
        [
            [power[b][order - a] if order >= a else Fraction(0) for a, b in monomials]
            for order in range(count)
        ],
        monomials,
    )


def _da_rows(
    series: tuple[Fraction, ...],
    degree_x: int,
    degree_f: int,
    degree_derivative: int,
) -> tuple[list[list[Fraction]], list[tuple[int, int, int]]]:
    count = len(series) - 1
    derivative = tuple((degree + 1) * series[degree + 1] for degree in range(count))
    f_power = _powers(tuple(series[:count]), degree_f, count - 1)
    derivative_power = _powers(derivative, degree_derivative, count - 1)
    monomials = [
        (a, b, c)
        for c in range(degree_derivative + 1)
        for b in range(degree_f + 1)
        for a in range(degree_x + 1)
    ]
    columns = []
    for a, b, c in monomials:
        product = _convolve(f_power[b], derivative_power[c], count - 1)
        columns.append(
            [product[order - a] if order >= a else Fraction(0) for order in range(count)]
        )
    return (
        [
            [columns[column][order] for column in range(len(columns))]
            for order in range(count)
        ],
        monomials,
    )


def _dot(row: Sequence[Fraction], vector: Sequence[int]) -> Fraction:
    return sum(
        (value * coefficient for value, coefficient in zip(row, vector, strict=True)),
        Fraction(0),
    )


def _parity_polynomials(width: int, height: int):
    lattice = square(width, height, periodic=False)
    table = [[0] * (ORDER + 1) for _ in range(1 << lattice.n_sites)]
    table[0][0] = 1
    for left, right in lattice.bonds:
        toggle = (1 << left) | (1 << right)
        previous = [row[:] for row in table]
        for mask, polynomial in enumerate(previous):
            target = table[mask ^ toggle]
            for degree in range(ORDER):
                target[degree + 1] += polynomial[degree]
    return lattice, table


def _divide(numerator: Sequence[int], denominator: Sequence[int]) -> tuple[Fraction, ...]:
    quotient = [Fraction(0) for _ in range(ORDER + 1)]
    for degree in range(ORDER + 1):
        lower = sum(
            (
                denominator[shift] * quotient[degree - shift]
                for shift in range(1, degree + 1)
            ),
            Fraction(0),
        )
        quotient[degree] = (numerator[degree] - lower) / denominator[0]
    return tuple(quotient)


def _box_overlap(width: int, height: int) -> tuple[Fraction, ...]:
    lattice, parity = _parity_polynomials(width, height)
    total = [Fraction(0) for _ in range(ORDER + 1)]
    for source in range(lattice.n_sites):
        for target in range(lattice.n_sites):
            if source == target:
                correlation = (Fraction(1),) + (Fraction(0),) * ORDER
            else:
                correlation = _divide(
                    parity[(1 << source) | (1 << target)], parity[0]
                )
            square_series = _convolve(correlation, correlation, ORDER)
            total = [left + right for left, right in zip(total, square_series, strict=True)]
    return tuple(total)


def _overlap_flm() -> tuple[Fraction, ...]:
    boxes = sorted(
        (
            (width, height)
            for width in range(1, ORDER // 2 + 2)
            for height in range(1, ORDER // 2 + 2)
            if width + height - 2 <= ORDER // 2
        ),
        key=lambda shape: (sum(shape), shape),
    )
    weights: dict[tuple[int, int], tuple[Fraction, ...]] = {}
    bulk = [Fraction(0) for _ in range(ORDER + 1)]
    for width, height in boxes:
        weight = list(_box_overlap(width, height))
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                weight = [
                    value - placements * subvalue
                    for value, subvalue in zip(weight, subweight, strict=True)
                ]
        weights[(width, height)] = tuple(weight)
        bulk = [value + delta for value, delta in zip(bulk, weight, strict=True)]
    return tuple(bulk)


def _recheck_survivors(
    result: dict, ht: tuple[Fraction, ...], lt: tuple[Fraction, ...]
) -> tuple[int, int]:
    checked_rows = 0
    checked_vectors = 0
    sections = (
        ("algebraicity", "HT", ht, False),
        ("algebraicity", "LT", lt, False),
        ("differential_algebraicity_order_1", "HT", ht, True),
        ("differential_algebraicity_order_1", "LT", lt, True),
    )
    for section, phase, series, differential in sections:
        for record in result["data"][section][phase]["frontier"]:
            if record["verdict"] != "CANDIDATE_SURVIVES(!)":
                continue
            if differential:
                rows, monomials = _da_rows(
                    series,
                    record["degree_x"],
                    record["degree_f"],
                    record["degree_f_prime"],
                )
            else:
                rows, monomials = _algebraic_rows(
                    series, record["degree_x"], record["degree_f"]
                )
            assert [list(monomial) for monomial in monomials] == record["monomials"]
            for vector in record["full_kernel_basis"]:
                assert any(vector)
                assert all(_dot(row, vector) == 0 for row in rows)
                checked_vectors += 1
            assert record["all_survivors_are_truncation_unobservable_monomials"]
            checked_rows += 1
    return checked_rows, checked_vectors


def main() -> None:
    ht = _load_even(HT)
    lt = _load_even(LT)
    result = json.loads(RESULT.read_text(encoding="utf-8"))

    ht_rows, _ = _algebraic_rows(ht, 2, 2)
    _check(
        "HT representative algebraic certificate",
        _rank(ht_rows[:11], 9) == 9,
        "the (dx,df)=(2,2) prefix matrix has exact rank 9 for 9 unknowns",
    )
    lt_rows, _ = _algebraic_rows(lt, 2, 3)
    _check(
        "LT representative algebraic certificate",
        _rank(lt_rows[:14], 12) == 12,
        "the (dx,df)=(2,3) prefix matrix has exact rank 12 for 12 unknowns",
    )

    ht_da_rows, _ = _da_rows(ht, 1, 1, 1)
    _check(
        "HT representative differential-algebraic certificate",
        _rank(ht_da_rows[:10], 8) == 8,
        "the (dx,df,df')=(1,1,1) prefix matrix has exact rank 8",
    )
    lt_da_rows, _ = _da_rows(lt, 1, 1, 1)
    _check(
        "LT representative differential-algebraic certificate",
        _rank(lt_da_rows[:10], 8) == 8,
        "the (dx,df,df')=(1,1,1) prefix matrix has exact rank 8",
    )

    survivor_rows, survivor_vectors = _recheck_survivors(result, ht, lt)
    _check(
        "all surviving candidates independently rechecked",
        survivor_rows == result["data"]["survivor_parent_audit"]["row_count"]
        and survivor_vectors > 0,
        f"{survivor_vectors} basis vectors in {survivor_rows} survivor rows vanish exactly on every available equation",
    )

    overlap = _overlap_flm()
    expected_overlap = (
        Fraction(1),
        Fraction(0),
        Fraction(4),
        Fraction(0),
        Fraction(36),
        Fraction(0),
        Fraction(236),
        Fraction(0),
        Fraction(1556),
    )
    _check(
        "independent 2D overlap extension",
        overlap == expected_overlap,
        "open-box parity enumeration and independent Fraction inversion give [v^8] sum_r G(r)^2=1556",
    )

    anisotropic = anisotropic_flm_c2_series(ORDER)
    residual = list(value / 2 for value in overlap)
    residual[0] -= Fraction(1, 2)
    _check(
        "interlayer c2 v8 crosscheck",
        tuple(residual) == anisotropic.residual and residual[8] == 778,
        "independent 2D and anisotropic 3D routes agree through v^8",
    )
    odd_zero = all(
        all(row[1] == 0 and row[3] == 0 for row in anisotropic_box_even_subgraph(shape, ORDER, 3))
        for shape in anisotropic.boxes
    )
    _check(
        "interlayer c3 parity check",
        odd_zero,
        "all 15 required anisotropic box polynomials have zero w and w^3 columns through v^8",
    )

    _check(
        "result schema and checks",
        set(result) == {"meta", "data", "checks"}
        and all(check["passed"] for check in result["checks"]),
        "machine-readable artifact has the required top-level schema and only computed passing checks",
    )
    print("PASS test_series_structure")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_series_structure: {error}")
        raise
