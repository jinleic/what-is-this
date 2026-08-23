#!/usr/bin/env python3
"""Clean-room verifier for the c8 closed formula and c10 structure certificate.

The verifier imports no experiment module.  It rebuilds Bell partitions with a
restricted-growth-string generator (different from the producers' insertion
recursion), reconstructs every c8 connected polynomial, recounts all c10
profiles, and redoes the exact formal-series checks.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import time
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "interlayer" / "c10_structure.json"
LEGACY_C8 = ROOT / "results" / "interlayer" / "c8_series.json"
CPU_LIMIT = 900.0
RSS_LIMIT = 1_900_000_000
STARTED = time.process_time()
FAILURES: list[str] = []

Atom = tuple[str, tuple[int, ...]]
Monomial = tuple[Atom, ...]
Polynomial = dict[Monomial, Fraction]
KINDS = {2: "G", 4: "U", 6: "W", 8: "W8"}


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def guard(stage: str) -> None:
    elapsed = time.process_time() - STARTED
    rss = max_rss_bytes()
    if elapsed >= CPU_LIMIT:
        raise RuntimeError(f"process-time cap at {stage}: {elapsed:.3f}s")
    if rss >= RSS_LIMIT:
        raise RuntimeError(f"RSS cap at {stage}: {rss}")


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"{'PASS' if passed else 'FAIL'}: {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def integer_partitions_rgs(n: int) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Set partitions of range(1,n+1) from restricted growth strings."""

    if n == 0:
        return ((),)
    rows: list[tuple[tuple[int, ...], ...]] = []

    def walk(labels: tuple[int, ...], maximum: int) -> None:
        if len(labels) == n:
            blocks = [[] for _ in range(maximum + 1)]
            for slot, label in enumerate(labels, start=1):
                blocks[label].append(slot)
            rows.append(tuple(tuple(block) for block in blocks))
            return
        for label in range(maximum + 2):
            walk(labels + (label,), max(maximum, label))

    walk((0,), 0)
    return tuple(rows)


@lru_cache(maxsize=None)
def partitions_of_slots(
    slots: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    rows = []
    for partition in integer_partitions_rgs(len(slots)):
        rows.append(
            tuple(tuple(slots[index - 1] for index in block) for block in partition)
        )
    return tuple(rows)


def atom(kind: str, slots: Iterable[int]) -> Atom:
    return kind, tuple(sorted(slots))


def monomial(*atoms: Atom) -> Monomial:
    return tuple(sorted(atoms))


@lru_cache(maxsize=None)
def moment_to_connected(slots: tuple[int, ...]) -> Polynomial:
    if not slots:
        return {(): Fraction(1)}
    out: Polynomial = {}
    for partition in partitions_of_slots(slots):
        if any(len(block) % 2 for block in partition):
            continue
        key = monomial(*(atom(KINDS[len(block)], block) for block in partition))
        out[key] = out.get(key, Fraction(0)) + 1
    return out


def polynomial_product(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for first, coefficient in left.items():
        for second, other in right.items():
            key = monomial(*first, *second)
            out[key] = out.get(key, Fraction(0)) + coefficient * other
    return {key: value for key, value in out.items() if value}


def compositions(total: int) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    return tuple(
        (first,) + tail
        for first in range(total, 0, -1)
        for tail in compositions(total - first)
    )


def gap_map(composition: Sequence[int]) -> dict[int, int]:
    result: dict[int, int] = {}
    slot = 1
    for gap, part in enumerate(composition):
        for _ in range(2 * int(part)):
            result[slot] = gap
            slot += 1
    return result


def profile_denominator(composition: Sequence[int]) -> int:
    result = 1
    for part in composition:
        result *= factorial(2 * int(part))
    return result


def block_is_even(gaps: dict[int, int], block: Sequence[int]) -> bool:
    parity = 0
    for slot in block:
        parity ^= 1 << gaps[slot]
    return parity == 0


def induced_layer_sets(
    gaps: dict[int, int], block: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    by_gap: dict[int, list[int]] = defaultdict(list)
    for slot in block:
        by_gap[gaps[slot]].append(slot)
    rows = []
    for layer in range(max(gaps.values()) + 2):
        slots = tuple(sorted(by_gap.get(layer - 1, []) + by_gap.get(layer, [])))
        if slots:
            rows.append(slots)
    return tuple(rows)


def clean_c8_profile(
    composition: tuple[int, ...], partitions8: Sequence[tuple[tuple[int, ...], ...]]
) -> tuple[Polynomial, int, dict[tuple[int, ...], int]]:
    gaps = gap_map(composition)
    placement = Fraction(1, profile_denominator(composition))
    output: Polynomial = {}
    count = 0
    class_weights: Counter[tuple[int, ...]] = Counter()
    for partition in partitions8:
        if not all(block_is_even(gaps, block) for block in partition):
            continue
        count += 1
        mu = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        class_weights[
            tuple(sorted((len(block) for block in partition), reverse=True))
        ] += mu
        term: Polynomial = {(): Fraction(1)}
        for block in partition:
            for slots in induced_layer_sets(gaps, block):
                term = polynomial_product(term, moment_to_connected(slots))
        for key, coefficient in term.items():
            output[key] = output.get(key, Fraction(0)) + placement * mu * coefficient
    return (
        {key: value for key, value in output.items() if value},
        count,
        dict(class_weights),
    )


def polynomial_digest(polynomials: dict[str, Polynomial]) -> str:
    digest = hashlib.sha256()
    for name in sorted(polynomials):
        digest.update(name.encode("utf-8"))
        for key, coefficient in sorted(polynomials[name].items()):
            digest.update(repr((key, coefficient)).encode("utf-8"))
    return digest.hexdigest()


def parse(values: Sequence[str]) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in values)


def series_multiply(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    output = [Fraction(0) for _ in range(order + 1)]
    for left_degree in range(min(len(left), order + 1)):
        for right_degree in range(min(len(right), order + 1 - left_degree)):
            output[left_degree + right_degree] += left[left_degree] * right[right_degree]
    return tuple(output)


def series_power(base: Sequence[Fraction], exponent: int, order: int) -> tuple[Fraction, ...]:
    output = (Fraction(1),) + (Fraction(0),) * order
    for _ in range(exponent):
        output = series_multiply(output, base, order)
    return output


def series_divide(
    numerator: Sequence[Fraction], denominator: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    output = [Fraction(0) for _ in range(order + 1)]
    for degree in range(order + 1):
        value = numerator[degree] if degree < len(numerator) else Fraction(0)
        for previous in range(degree):
            denominator_degree = degree - previous
            if denominator_degree < len(denominator):
                value -= output[previous] * denominator[denominator_degree]
        output[degree] = value / denominator[0]
    return tuple(output)


def series_log(base: Sequence[Fraction], order: int) -> tuple[Fraction, ...]:
    derivative = tuple(
        Fraction((degree + 1) * base[degree + 1])
        for degree in range(order)
    ) + (Fraction(0),)
    quotient = series_divide(derivative, base, order)
    return (Fraction(0),) + tuple(
        quotient[degree - 1] / degree for degree in range(1, order + 1)
    )


def add_rows(rows: Iterable[Sequence[Fraction]]) -> tuple[Fraction, ...]:
    materialized = tuple(rows)
    return tuple(
        sum((row[index] for row in materialized), Fraction(0))
        for index in range(len(materialized[0]))
    )


def verify_c8(
    artifact: dict[str, object], legacy: dict[str, object]
) -> None:
    partitions8 = integer_partitions_rgs(8)
    check("clean Bell(8)", len(partitions8) == 4140, str(len(partitions8)))
    stored = artifact["data"]["c8_closed"]["connected_formula"]
    stored_rows = {tuple(row["composition"]): row for row in stored["profiles"]}
    legacy_rows = {
        tuple(row["pattern"]): row
        for row in legacy["data"]["bell_eighth_order_inventory"]["patterns"]
    }
    polynomials: dict[str, Polynomial] = {}
    top = atom("W8", range(1, 9))
    square = monomial(top, top)
    square_coefficients = {}
    for composition in compositions(4):
        polynomial, count, classes = clean_c8_profile(composition, partitions8)
        polynomials[str(composition)] = polynomial
        digest = polynomial_digest({str(composition): polynomial})
        w8_count = sum(any(kind == "W8" for kind, _slots in key) for key in polynomial)
        square_coefficients[composition] = polynomial.get(square, Fraction(0))
        row = stored_rows[composition]
        legacy_row = legacy_rows[composition]
        check(
            f"clean c8 profile {composition}",
            count == row["admissible_bell_partitions"] == legacy_row["admissible_partitions"]
            and len(polynomial) == row["connected_monomials"] == legacy_row["connected_monomials"]
            and digest == row["connected_polynomial_sha256"] == legacy_row["connected_polynomial_sha256"]
            and w8_count == row["w8_monomials"]
            and {
                "+".join(map(str, key)): str(value)
                for key, value in sorted(classes.items(), reverse=True)
            }
            == row["moebius_by_block_class"],
            f"partitions={count}, monomials={len(polynomial)}, W8={w8_count}",
        )
        guard(f"c8 {composition}")
    combined = polynomial_digest(polynomials)
    check(
        "clean complete c8 digest",
        combined
        == stored["complete_sha256"]
        == legacy["data"]["bell_eighth_order_inventory"]["collected_digest_sha256"],
        combined,
    )
    expected_square = {
        composition: (Fraction(1, factorial(8)) if composition == (4,) else Fraction(0))
        for composition in compositions(4)
    }
    check(
        "clean W8 square obstruction",
        square_coefficients == expected_square,
        "unique coefficient 1/8!",
    )

    legacy_data = legacy["data"]
    lower = legacy_data["reproduced_lower_orders"]
    c2 = parse(lower["c2_total_q_or_w"])
    c4 = parse(lower["c4_direct_q"])
    c6 = parse(lower["c6_direct_q"])
    c8 = parse(legacy_data["c8_direct_coupling_q"]["full_infinite_stack"])
    total = parse(legacy_data["c8_wave_variable_w_tanh_Kz"]["total"])
    residual = parse(
        legacy_data["c8_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"]
    )
    atanh = tuple(
        Fraction(1, degree) if degree % 2 else Fraction(0) for degree in range(9)
    )
    weights = {power: series_power(atanh, power, 8)[8] for power in (2, 4, 6, 8)}
    rebuilt = tuple(
        c8[index] + weights[6] * c6[index] + weights[4] * c4[index] + weights[2] * c2[index]
        for index in range(13)
    )
    rebuilt_residual = list(rebuilt)
    rebuilt_residual[0] -= Fraction(1, 8)
    rebuilt_residual = tuple(rebuilt_residual)
    q_heights = tuple(
        parse(row)
        for row in legacy_data["c8_direct_coupling_q"]["by_exact_vertical_extent"].values()
    )
    w_heights = tuple(
        parse(row)
        for row in legacy_data["c8_wave_variable_w_tanh_Kz"][
            "residual_by_exact_vertical_extent"
        ].values()
    )
    check(
        "clean c8 v0-v12 formal audit",
        add_rows(q_heights) == c8
        and rebuilt == total
        and rebuilt_residual == residual
        and add_rows(w_heights) == residual
        and parse(legacy_data["second_route_w_log"]["residual_w8"]) == residual
        and parse(legacy_data["independent_spin_dos_route"]["residual_w8"]) == residual[:7],
        "all 13 stored degrees; independent overlap through v6",
    )
    stored_audit = artifact["data"]["c8_closed"]["formal_v12_audit"]
    check(
        "artifact c8 coefficient rows",
        parse(stored_audit["direct_q"]) == c8
        and parse(stored_audit["total_w"]) == total
        and parse(stored_audit["residual_w"]) == residual
        and stored_audit["degrees_checked"] == list(range(13)),
        "artifact agrees with clean exact reconstruction",
    )


def c10_profile_inventory(
    composition: tuple[int, ...], partitions10: Sequence[tuple[tuple[int, ...], ...]]
) -> tuple[int, dict[str, str]]:
    gaps = gap_map(composition)
    count = 0
    classes: Counter[tuple[int, ...]] = Counter()
    for partition in partitions10:
        if not all(block_is_even(gaps, block) for block in partition):
            continue
        count += 1
        mu = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        classes[tuple(sorted((len(block) for block in partition), reverse=True))] += mu
    return count, {
        "+".join(map(str, key)): str(value)
        for key, value in sorted(classes.items(), reverse=True)
    }


def layer_degrees(composition: Sequence[int]) -> tuple[int, ...]:
    parts = (0,) + tuple(2 * part for part in composition) + (0,)
    return tuple(parts[index] + parts[index + 1] for index in range(len(parts) - 1))


def exact_c10_low_columns() -> dict[str, str]:
    order = 10
    exp_plus = tuple(Fraction(1, factorial(degree)) for degree in range(order + 1))
    exp_minus = tuple(
        Fraction((-1) ** degree, factorial(degree)) for degree in range(order + 1)
    )
    cosh = tuple((left + right) / 2 for left, right in zip(exp_plus, exp_minus, strict=True))
    sinh = tuple((left - right) / 2 for left, right in zip(exp_plus, exp_minus, strict=True))
    tanh = series_divide(sinh, cosh, order)
    log_cosh = series_log(cosh, order)
    tanh2 = series_multiply(tanh, tanh, order)
    denominator = list(-value for value in tanh2)
    denominator[0] += 1
    ladder_q = tuple(2 * value for value in series_divide(tanh2, denominator, order))
    atanh = tuple(
        Fraction(1, degree) if degree % 2 else Fraction(0)
        for degree in range(order + 1)
    )
    weights = {
        even: series_power(atanh, even, order)[10]
        for even in range(2, 11, 2)
    }
    return {
        "direct_q_v0": str(log_cosh[10]),
        "direct_q_v2": str(ladder_q[10]),
        "total_w_v0": str(sum(weights[n] * log_cosh[n] for n in weights)),
        "residual_w_v2": str(sum(weights[n] * ladder_q[n] for n in weights)),
        **{f"weight_{n}": str(value) for n, value in weights.items()},
    }


def verify_c10(artifact: dict[str, object]) -> None:
    partitions10 = integer_partitions_rgs(10)
    check("clean Bell(10)", len(partitions10) == 115975, str(len(partitions10)))
    structure = artifact["data"]["c10_structure"]
    rows = {
        tuple(row["composition"]): row
        for row in structure["connected_formula"]["profiles"]
    }
    reflected = True
    w10_support = set()
    for composition in compositions(5):
        count, classes = c10_profile_inventory(composition, partitions10)
        row = rows[composition]
        degrees = layer_degrees(composition)
        if len(composition) == 1:
            w10_count = sum(
                all(len(block) % 2 == 0 for block in partition)
                for partition in integer_partitions_rgs(10)
            )
        elif len(composition) == 2:
            w10_count = len(moment_to_connected(tuple(range(1, 2 * composition[0] + 1)))) * len(
                moment_to_connected(tuple(range(1, 2 * composition[1] + 1)))
            )
        else:
            w10_count = 0
        if w10_count:
            w10_support.add(composition)
        check(
            f"clean c10 profile {composition}",
            count == row["admissible_bell_partitions"]
            and classes == row["moebius_by_block_class"]
            and factorial(10) // profile_denominator(composition) == row["ordered_gap_assignments"]
            and str(Fraction(1, profile_denominator(composition))) == row["placement"]
            and list(degrees) == row["layer_slot_degrees"]
            and max(degrees) == row["maximum_connected_rank"]
            and w10_count == row["w10_sector_monomials"],
            f"partitions={count}, W10={w10_count}",
        )
        reverse = rows[tuple(reversed(composition))]
        reflected &= count == reverse["admissible_bell_partitions"] and classes == reverse[
            "moebius_by_block_class"
        ]
        guard(f"c10 {composition}")
    check("clean c10 reflection", reflected)
    check(
        "clean W10 support and square",
        w10_support == {(5,), (4, 1), (3, 2), (2, 3), (1, 4)}
        and structure["w10_sector"]["w10_square_coefficient"] == str(Fraction(1, factorial(10))),
        str(sorted(w10_support)),
    )
    low = exact_c10_low_columns()
    stored_low = structure["exact_low_order_c10"]
    check(
        "clean c10 v0 and v2 columns",
        stored_low["direct_q_v0"] == low["direct_q_v0"]
        and stored_low["direct_q_v2"] == low["direct_q_v2"]
        and stored_low["total_w_v0"] == low["total_w_v0"]
        and stored_low["residual_w_v2"] == low["residual_w_v2"]
        and all(
            stored_low["q_to_w_weights"][str(n)] == low[f"weight_{n}"]
            for n in range(2, 11, 2)
        ),
        str(low),
    )
    check(
        "artifact preserves finite c10 scope",
        structure["finite_series_scope"]["tag"] == "[UNRESOLVED]"
        and "No new c10 finite-lattice boxes" in structure["finite_series_scope"]["statement"],
    )


def main() -> int:
    if not ARTIFACT.is_file() or not LEGACY_C8.is_file():
        print("FAIL test_interlayer_c10: required artifact missing")
        return 1
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    legacy = json.loads(LEGACY_C8.read_text(encoding="utf-8"))
    check("artifact top-level schema", set(artifact) == {"meta", "data", "checks"})
    check(
        "all producer checks passed",
        bool(artifact["checks"]) and all(row.get("passed") is True for row in artifact["checks"]),
        f"checks={len(artifact['checks'])}",
    )
    verify_c8(artifact, legacy)
    guard("after c8 verification")
    verify_c10(artifact)
    source_hashes = artifact["meta"]["source_sha256"]
    check(
        "artifact source hashes",
        all(
            (ROOT / relative).is_file()
            and sha256_file(ROOT / relative) == expected
            for relative, expected in source_hashes.items()
        ),
        f"files={len(source_hashes)}",
    )
    check(
        "resource limits",
        time.process_time() - STARTED < CPU_LIMIT and max_rss_bytes() < RSS_LIMIT,
        f"cpu={time.process_time()-STARTED:.3f}s, rss={max_rss_bytes()}",
    )
    if FAILURES:
        print("FAIL test_interlayer_c10: " + ", ".join(FAILURES))
        return 1
    print("PASS test_interlayer_c10")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
