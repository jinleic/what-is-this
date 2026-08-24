#!/usr/bin/env python3
"""Exact all-v connected-correlation formula and obstruction at interlayer order eight.

This module is a new, pure-stdlib derivation.  It does not import the earlier
c8 producer.  The previous artifact is used only after reconstruction, as a
regression target for exact polynomial digests and finite v-series data.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
LEGACY_C8 = ROOT / "results" / "interlayer" / "c8_series.json"

Atom = tuple[str, tuple[int, ...]]
Monomial = tuple[Atom, ...]
Polynomial = dict[Monomial, Fraction]

ATOM_KIND = {2: "G", 4: "U", 6: "W", 8: "W8"}


def record(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def atom(kind: str, slots: Iterable[int]) -> Atom:
    return kind, tuple(sorted(slots))


def monomial(*atoms: Atom) -> Monomial:
    return tuple(sorted(atoms))


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for first, coefficient in left.items():
        for second, other in right.items():
            key = monomial(*first, *second)
            out[key] = out.get(key, Fraction(0)) + coefficient * other
    return {key: value for key, value in out.items() if value}


@lru_cache(maxsize=None)
def set_partitions(
    items: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Return each unlabeled set partition exactly once."""

    if not items:
        return ((),)
    first, rest = items[0], items[1:]
    rows: list[tuple[tuple[int, ...], ...]] = []
    for partition in set_partitions(rest):
        rows.append(((first,),) + partition)
        for index in range(len(partition)):
            rows.append(
                partition[:index]
                + ((first,) + partition[index],)
                + partition[index + 1 :]
            )
    return tuple(rows)


@lru_cache(maxsize=None)
def even_partitions(
    slots: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Return partitions whose blocks all have positive even size."""

    if not slots:
        return ((),)
    first, rest = slots[0], slots[1:]
    rows: list[tuple[tuple[int, ...], ...]] = []
    for size in range(2, len(slots) + 1, 2):
        for partners in combinations(rest, size - 1):
            chosen = set(partners)
            remaining = tuple(slot for slot in rest if slot not in chosen)
            for tail in even_partitions(remaining):
                rows.append(((first,) + partners,) + tail)
    return tuple(rows)


@lru_cache(maxsize=None)
def connected_expansion(slots: tuple[int, ...]) -> Polynomial:
    """Expand a zero-field moment into even connected spin cumulants."""

    if not slots:
        return {(): Fraction(1)}
    if len(slots) % 2:
        return {}
    out: Polynomial = {}
    for partition in even_partitions(slots):
        key = monomial(
            *(atom(ATOM_KIND[len(block)], block) for block in partition)
        )
        out[key] = out.get(key, Fraction(0)) + 1
    return out


@lru_cache(maxsize=None)
def compositions(total: int) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    return tuple(
        (first,) + tail
        for first in range(total, 0, -1)
        for tail in compositions(total - first)
    )


def gap_map(composition: Sequence[int]) -> dict[int, int]:
    return {
        slot: gap
        for slot, gap in enumerate(
            (
                gap
                for gap, part in enumerate(composition)
                for _ in range(2 * int(part))
            ),
            start=1,
        )
    }


def placement(composition: Sequence[int]) -> Fraction:
    denominator = 1
    for part in composition:
        denominator *= factorial(2 * int(part))
    return Fraction(1, denominator)


def ordered_gap_assignments(composition: Sequence[int]) -> int:
    return factorial(2 * sum(composition)) * placement(composition).numerator // placement(
        composition
    ).denominator


def block_gap_even(gaps: dict[int, int], block: Sequence[int]) -> bool:
    counts = Counter(gaps[slot] for slot in block)
    return all(count % 2 == 0 for count in counts.values())


def block_layer_sets(
    gaps: dict[int, int], block: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    per_gap: dict[int, list[int]] = defaultdict(list)
    for slot in block:
        per_gap[gaps[slot]].append(slot)
    rows = []
    for layer in range(max(gaps.values()) + 2):
        slots = tuple(sorted(per_gap.get(layer - 1, []) + per_gap.get(layer, [])))
        if slots:
            rows.append(slots)
    return tuple(rows)


def profile_polynomial(
    composition: tuple[int, ...]
) -> tuple[Polynomial, int, dict[tuple[int, ...], int]]:
    """Build the full connected-correlation polynomial for one c8 profile."""

    gaps = gap_map(composition)
    scale = placement(composition)
    out: Polynomial = {}
    admissible = 0
    moebius: Counter[tuple[int, ...]] = Counter()
    for partition in set_partitions(tuple(range(1, 9))):
        if not all(block_gap_even(gaps, block) for block in partition):
            continue
        admissible += 1
        mu = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        moebius[tuple(sorted((len(block) for block in partition), reverse=True))] += mu
        term: Polynomial = {(): Fraction(1)}
        for block in partition:
            for slots in block_layer_sets(gaps, block):
                term = multiply(term, connected_expansion(slots))
        for key, coefficient in term.items():
            out[key] = out.get(key, Fraction(0)) + scale * mu * coefficient
    return {key: value for key, value in out.items() if value}, admissible, dict(moebius)


def polynomial_digest(polynomials: dict[str, Polynomial]) -> str:
    digest = hashlib.sha256()
    for name in sorted(polynomials):
        digest.update(name.encode("utf-8"))
        for key, coefficient in sorted(polynomials[name].items()):
            digest.update(repr((key, coefficient)).encode("utf-8"))
    return digest.hexdigest()


def top_atom_sector(polynomial: Polynomial, kind: str) -> Polynomial:
    return {
        key: coefficient
        for key, coefficient in polynomial.items()
        if any(atom_kind == kind for atom_kind, _slots in key)
    }


def predicted_w8_sector(composition: tuple[int, ...]) -> Polynomial:
    """Closed W8 sector forced by the unique one-block Bell partition."""

    full_slots = tuple(range(1, 9))
    full_atom = atom("W8", full_slots)
    scale = placement(composition)
    out: Polynomial = {}
    if len(composition) == 1:
        for lower, coefficient in connected_expansion(full_slots).items():
            key = monomial(full_atom, *lower)
            multiplicity = 1 if lower == (full_atom,) else 2
            out[key] = scale * multiplicity * coefficient
    elif len(composition) == 2:
        split = 2 * composition[0]
        left = connected_expansion(tuple(range(1, split + 1)))
        right = connected_expansion(tuple(range(split + 1, 9)))
        for left_key, left_coefficient in left.items():
            for right_key, right_coefficient in right.items():
                key = monomial(full_atom, *left_key, *right_key)
                out[key] = out.get(key, Fraction(0)) + (
                    scale * left_coefficient * right_coefficient
                )
    return out


def fractions(values: Sequence[str]) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in values)


def sum_rows(rows: Iterable[Sequence[Fraction]]) -> tuple[Fraction, ...]:
    rows = tuple(rows)
    if not rows:
        return ()
    return tuple(sum((row[index] for row in rows), Fraction(0)) for index in range(len(rows[0])))


def multiply_series(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    out = [Fraction(0) for _ in range(order + 1)]
    for i, first in enumerate(left):
        if not first:
            continue
        for j, second in enumerate(right[: order + 1 - i]):
            if second:
                out[i + j] += first * second
    return tuple(out)


def power_series(base: Sequence[Fraction], exponent: int, order: int) -> tuple[Fraction, ...]:
    out = (Fraction(1),) + (Fraction(0),) * order
    factor = tuple(base)
    for _ in range(exponent):
        out = multiply_series(out, factor, order)
    return out


def formal_series_audit(
    legacy: dict[str, object], checks: list[dict[str, object]]
) -> dict[str, object]:
    """Recheck every stored c8 coefficient through v^12 by exact identities."""

    data = legacy["data"]
    lower = data["reproduced_lower_orders"]
    c2 = fractions(lower["c2_total_q_or_w"])
    c4 = fractions(lower["c4_direct_q"])
    c6 = fractions(lower["c6_direct_q"])
    c8 = fractions(data["c8_direct_coupling_q"]["full_infinite_stack"])
    total = fractions(data["c8_wave_variable_w_tanh_Kz"]["total"])
    residual = fractions(
        data["c8_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"]
    )

    atanh = [Fraction(0) for _ in range(9)]
    for degree in range(1, 9, 2):
        atanh[degree] = Fraction(1, degree)
    conversion = {
        power: power_series(atanh, power, 8)[8] for power in (2, 4, 6, 8)
    }
    rebuilt_total = tuple(
        c8[index]
        + conversion[6] * c6[index]
        + conversion[4] * c4[index]
        + conversion[2] * c2[index]
        for index in range(13)
    )
    rebuilt_residual = list(rebuilt_total)
    rebuilt_residual[0] -= Fraction(1, 8)
    rebuilt_residual = tuple(rebuilt_residual)

    q_heights = {
        key: fractions(row)
        for key, row in data["c8_direct_coupling_q"]["by_exact_vertical_extent"].items()
    }
    w_heights = {
        key: fractions(row)
        for key, row in data["c8_wave_variable_w_tanh_Kz"][
            "residual_by_exact_vertical_extent"
        ].items()
    }
    second_route = fractions(data["second_route_w_log"]["residual_w8"])
    spin_route = fractions(data["independent_spin_dos_route"]["residual_w8"])

    record(
        checks,
        "c8 exact atanh conversion weights",
        conversion
        == {2: Fraction(44, 105), 4: Fraction(22, 15), 6: Fraction(2), 8: Fraction(1)},
        ", ".join(f"q^{power}:{value}" for power, value in conversion.items()),
    )
    record(
        checks,
        "c8 all thirteen q coefficients equal exact-height sum",
        sum_rows(q_heights.values()) == c8,
        "v^0 through v^12, including every zero odd degree",
    )
    record(
        checks,
        "c8 all thirteen w coefficients rebuilt from q",
        rebuilt_total == total and rebuilt_residual == residual,
        "exact q=atanh(w) composition and 1/8 prefactor subtraction",
    )
    record(
        checks,
        "c8 all thirteen residual coefficients equal exact-height sum",
        sum_rows(w_heights.values()) == residual,
        "heights two through five",
    )
    record(
        checks,
        "c8 independent stored routes cover every available coefficient",
        second_route == residual and spin_route == residual[:7],
        "w-log through v^12 and spin-DOS through v^6",
    )
    record(
        checks,
        "legacy c8 certificate checks passed",
        all(row.get("passed") is True for row in legacy["checks"]),
        f"{len(legacy['checks'])} producer checks",
    )
    return {
        "tag": "[COMPUTATION]",
        "scope": (
            "[COMPUTATION] Every stored coefficient at degrees v^0,...,v^12 is "
            "rechecked by exact height summation and q-to-w composition; the independent "
            "w-log row covers v^12 and the spin-DOS row covers v^6."
        ),
        "atanh_w8_weights": {str(key): str(value) for key, value in conversion.items()},
        "direct_q": [str(value) for value in c8],
        "total_w": [str(value) for value in rebuilt_total],
        "residual_w": [str(value) for value in rebuilt_residual],
        "degrees_checked": list(range(13)),
    }


def build_c8_closed() -> tuple[dict[str, object], list[dict[str, object]]]:
    if not LEGACY_C8.is_file():
        raise FileNotFoundError(LEGACY_C8)
    legacy = json.loads(LEGACY_C8.read_text(encoding="utf-8"))
    legacy_inventory = legacy["data"]["bell_eighth_order_inventory"]
    legacy_rows = {
        tuple(row["pattern"]): row for row in legacy_inventory["patterns"]
    }

    checks: list[dict[str, object]] = []
    polynomials: dict[str, Polynomial] = {}
    rows: list[dict[str, object]] = []
    square = monomial(atom("W8", range(1, 9)), atom("W8", range(1, 9)))
    square_coefficients: dict[tuple[int, ...], Fraction] = {}

    for composition in compositions(4):
        polynomial, admissible, moebius = profile_polynomial(composition)
        polynomials[str(composition)] = polynomial
        digest = polynomial_digest({str(composition): polynomial})
        sector = top_atom_sector(polynomial, "W8")
        predicted = predicted_w8_sector(composition)
        square_coefficients[composition] = polynomial.get(square, Fraction(0))
        legacy_row = legacy_rows[composition]
        record(
            checks,
            f"c8 profile {composition} exact legacy reconstruction",
            admissible == legacy_row["admissible_partitions"]
            and len(polynomial) == legacy_row["connected_monomials"]
            and digest == legacy_row["connected_polynomial_sha256"],
            f"partitions={admissible}, monomials={len(polynomial)}, sha256={digest}",
        )
        record(
            checks,
            f"c8 profile {composition} W8 sector theorem",
            sector == predicted,
            f"W8 monomials={len(sector)}",
        )
        histogram = Counter(tuple(sorted(kind for kind, _slots in key)) for key in polynomial)
        rows.append(
            {
                "composition": list(composition),
                "gap_multiplicities": [2 * part for part in composition],
                "layers": len(composition) + 1,
                "placement": str(placement(composition)),
                "ordered_gap_assignments": ordered_gap_assignments(composition),
                "admissible_bell_partitions": admissible,
                "connected_monomials": len(polynomial),
                "w8_monomials": len(sector),
                "connected_polynomial_sha256": digest,
                "moebius_by_block_class": {
                    "+".join(map(str, key)): str(value)
                    for key, value in sorted(moebius.items(), reverse=True)
                },
                "atom_profile_histogram": {
                    ",".join(key): value for key, value in sorted(histogram.items())
                },
            }
        )

    combined_digest = polynomial_digest(polynomials)
    record(
        checks,
        "c8 complete connected polynomial digest",
        combined_digest == legacy_inventory["collected_digest_sha256"],
        combined_digest,
    )
    expected_square = {
        composition: (Fraction(1, factorial(8)) if composition == (4,) else Fraction(0))
        for composition in compositions(4)
    }
    record(
        checks,
        "c8 irreducible W8 square coefficient",
        square_coefficients == expected_square,
        "coefficient 1/8! occurs only in the one-gap profile",
    )
    record(
        checks,
        "c8 composition count and Bell number",
        len(rows) == 8 and len(set_partitions(tuple(range(1, 9)))) == 4140,
        "2^(4-1)=8 and Bell(8)=4140",
    )

    formal = formal_series_audit(legacy, checks)
    if not all(bool(row["passed"]) for row in checks):
        failed = [str(row["name"]) for row in checks if not row["passed"]]
        raise AssertionError("c8 closed derivation failed: " + ", ".join(failed))

    data = {
        "tag": "[THEOREM][LEMMA][COMPUTATION]",
        "connected_formula": {
            "tag": "[THEOREM]",
            "statement": (
                "[THEOREM] For c_(8,q), sum over the eight compositions a of 4 with "
                "weight 1/prod_j(2a_j)!. For each composition sum the Bell Moebius "
                "coefficient over partitions whose every block has even count in every "
                "gap; replace every induced layer moment M(S) by the even-partition "
                "expansion E(S)=sum_{rho in EP(S)} prod_{D in rho} C(D). Anchoring x1=0 "
                "and summing x2,...,x8 gives an exact coefficientwise all-v identity."
            ),
            "notation": {
                "G": "[LEMMA] C(D) for |D|=2",
                "U": "[LEMMA] C(D) for |D|=4",
                "W": "[LEMMA] C(D) for |D|=6",
                "W8": "[LEMMA] C(D) for |D|=8",
                "coincident_sites": (
                    "[LEMMA] Slots remain labelled; coincident physical sites are handled by "
                    "the ordinary joint-cumulant identity for the repeated spin variables."
                ),
            },
            "profiles": rows,
            "total_admissible_bell_partitions": sum(
                row["admissible_bell_partitions"] for row in rows
            ),
            "total_connected_monomials_before_orbit_collapse": sum(
                row["connected_monomials"] for row in rows
            ),
            "complete_sha256": combined_digest,
        },
        "w8_sector": {
            "tag": "[THEOREM]",
            "statement": (
                "[THEOREM] Write E_8=W8+R_8. The complete W8 sector is "
                "(W8^2+2 W8 R_8)/8! for composition (4), plus "
                "E(I_a) W8 E(J_a)/((2a)!(8-2a)!) for compositions (a,4-a), "
                "a=1,2,3; profiles with at least three gaps contain no W8."
            ),
            "profile_w8_monomial_counts": {
                str(tuple(row["composition"])): row["w8_monomials"] for row in rows
            },
            "w8_square_coefficient": "1/40320",
            "obstruction": (
                "[THEOREM] In the universal polynomial ring of connected layer cumulants, "
                "the monomial W8(1,...,8)^2 has coefficient 1/8! and cannot belong to the "
                "subring generated by G,U,W. Thus partition-cumulant algebra cannot close "
                "c8 with the previously solved objects alone."
            ),
            "scope": (
                "[UNRESOLVED] This algebraic-independence obstruction does not exclude a "
                "special square-lattice identity that evaluates W8 by other non-polynomial "
                "or model-specific data."
            ),
        },
        "formal_v12_audit": formal,
    }
    return data, checks


def main() -> int:
    data, checks = build_c8_closed()
    print(
        json.dumps(
            {
                "patterns": len(data["connected_formula"]["profiles"]),
                "connected_monomials": data["connected_formula"][
                    "total_connected_monomials_before_orbit_collapse"
                ],
                "checks": len(checks),
            },
            sort_keys=True,
        )
    )
    print("PASS e221 c8 all-v connected formula and W8 obstruction")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
