"""Exact closed 2D-correlation reduction for the sixth interlayer cumulant.

The producer derives the finite G/U/W identity from all Bell(6)=203 set
partitions, then independently re-expands its underlying layer-moment form by
an exact 2D boundary-mask transfer and planar finite-lattice inversion.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import hashlib
import json
import sys
import time
from typing import Iterable, Sequence


SCRIPT = "experiments/e99_c6_closed.py"
ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "results" / "interlayer" / "c6_series.json"
RESULT = ROOT / "results" / "interlayer" / "c6_closed.json"
ORDER = 12
MAX_LAYERS = 4

IntegerSeries = tuple[int, ...]
RationalSeries = tuple[Fraction, ...]
Series = tuple[int | Fraction, ...]
Atom = tuple[str, tuple[int, ...]]
Monomial = tuple[Atom, ...]
Polynomial = dict[Monomial, Fraction]


# ---------------------------------------------------------------------------
# Exact formal series and open-layer boundary transfer
# ---------------------------------------------------------------------------
def _int_zero(order: int = ORDER) -> list[int]:
    return [0] * (order + 1)


def _zero(order: int = ORDER) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def _add_scaled(
    target: list[int | Fraction], source: Sequence[int | Fraction], scale: int | Fraction = 1
) -> None:
    for degree, value in enumerate(source):
        if value:
            target[degree] += scale * value


def _mul(
    left: Sequence[int | Fraction], right: Sequence[int | Fraction], order: int = ORDER
) -> Series:
    out: list[int | Fraction] = [0] * (order + 1)
    for first, left_value in enumerate(left[: order + 1]):
        if not left_value:
            continue
        for second, right_value in enumerate(right[: order + 1 - first]):
            if right_value:
                out[first + second] += left_value * right_value
    return tuple(out)


def _pow(base: Sequence[int | Fraction], exponent: int, order: int = ORDER) -> Series:
    result: Series = (1,) + (0,) * order
    for _ in range(exponent):
        result = _mul(result, base, order)
    return result


def _divide_unit(
    numerator: Sequence[int], denominator: Sequence[int], order: int = ORDER
) -> IntegerSeries:
    if denominator[0] != 1:
        raise ValueError("formal division requires unit constant denominator")
    quotient = _int_zero(order)
    for degree in range(order + 1):
        value = numerator[degree] if degree < len(numerator) else 0
        for lower in range(degree):
            value -= quotient[lower] * denominator[degree - lower]
        quotient[degree] = value
    return tuple(quotient)


def _square_bonds(width: int, height: int) -> tuple[tuple[int, int], ...]:
    bonds: list[tuple[int, int]] = []
    for y in range(height):
        for x in range(width):
            site = y * width + x
            if x + 1 < width:
                bonds.append((site, site + 1))
            if y + 1 < height:
                bonds.append((site, site + width))
    return tuple(bonds)


@lru_cache(maxsize=None)
def _boundary_polynomials(width: int, height: int, order: int = ORDER) -> tuple[IntegerSeries, ...]:
    """P_S(v) for every odd-boundary mask S of one open 2D rectangle."""
    sites = width * height
    values = [[0] * (order + 1) for _ in range(1 << sites)]
    values[0][0] = 1
    for first, second in _square_bonds(width, height):
        toggle = (1 << first) | (1 << second)
        for degree in range(order - 1, -1, -1):
            for mask in range(1 << sites):
                coefficient = values[mask][degree]
                if coefficient:
                    values[mask ^ toggle][degree + 1] += coefficient
    return tuple(tuple(row) for row in values)


@lru_cache(maxsize=None)
def _layer_moments(width: int, height: int, order: int = ORDER) -> tuple[IntegerSeries, ...]:
    """m_S=P_S/P_empty, exact open-layer zero-field spin moments."""
    parity = _boundary_polynomials(width, height, order)
    return tuple(_divide_unit(value, parity[0], order) for value in parity)


# ---------------------------------------------------------------------------
# Complete partition algebra and the G/U/W connected reduction
# ---------------------------------------------------------------------------
def _set_partitions(items: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Return every unlabeled set partition once, in insertion canonical order."""
    @lru_cache(maxsize=None)
    def build(remaining: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
        if not remaining:
            return ((),)
        first = remaining[0]
        rows: list[tuple[tuple[int, ...], ...]] = []
        for partition in build(remaining[1:]):
            rows.append(((first,),) + partition)
            for index, block in enumerate(partition):
                rows.append(partition[:index] + ((first,) + block,) + partition[index + 1 :])
        return tuple(rows)

    return build(items)


PARTITIONS_6 = _set_partitions((1, 2, 3, 4, 5, 6))


def _mu(partition: tuple[tuple[int, ...], ...]) -> int:
    count = len(partition)
    return (-1) ** (count - 1) * factorial(count - 1)


def _atom(kind: str, slots: Iterable[int]) -> Atom:
    return kind, tuple(sorted(slots))


def _monomial(*atoms: Atom) -> Monomial:
    return tuple(sorted(atoms))


def _poly_add(target: Polynomial, source: Polynomial, scale: Fraction = Fraction(1)) -> Polynomial:
    out: defaultdict[Monomial, Fraction] = defaultdict(Fraction)
    for mono, coefficient in target.items():
        out[mono] += coefficient
    for mono, coefficient in source.items():
        out[mono] += scale * coefficient
    return {mono: coefficient for mono, coefficient in out.items() if coefficient}


def _poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    out: defaultdict[Monomial, Fraction] = defaultdict(Fraction)
    for first, first_coefficient in left.items():
        for second, second_coefficient in right.items():
            out[_monomial(*first, *second)] += first_coefficient * second_coefficient
    return {mono: coefficient for mono, coefficient in out.items() if coefficient}


def _pair_partitions(slots: tuple[int, ...]) -> tuple[tuple[tuple[int, int], ...], ...]:
    if not slots:
        return ((),)
    first = slots[0]
    rows: list[tuple[tuple[int, int], ...]] = []
    for index in range(1, len(slots)):
        second = slots[index]
        rest = slots[1:index] + slots[index + 1 :]
        for pairing in _pair_partitions(rest):
            rows.append(((first, second),) + pairing)
    return tuple(rows)


@lru_cache(maxsize=None)
def _moment_connected_expansion(slots: tuple[int, ...]) -> Polynomial:
    """M_I expanded as connected 2D G/U/W objects for an even slot set I.

    M_I denotes the layer moment of the product of the spins at the listed
    slots.  Repeated physical positions are allowed later: this is the
    ordinary moment--cumulant identity for labeled variables, so it remains
    exact after slots coincide.
    """
    length = len(slots)
    if length == 0:
        return {(): Fraction(1)}
    if length % 2:
        return {}
    if length == 2:
        return {_monomial(_atom("G", slots)): Fraction(1)}
    if length == 4:
        out: Polynomial = {_monomial(_atom("U", slots)): Fraction(1)}
        for pairing in _pair_partitions(slots):
            mono = _monomial(*(_atom("G", pair) for pair in pairing))
            out[mono] = out.get(mono, Fraction(0)) + 1
        return out
    if length == 6:
        out = {_monomial(_atom("W", slots)): Fraction(1)}
        for pair in combinations(slots, 2):
            rest = tuple(slot for slot in slots if slot not in pair)
            mono = _monomial(_atom("G", pair), _atom("U", rest))
            out[mono] = out.get(mono, Fraction(0)) + 1
        for pairing in _pair_partitions(slots):
            mono = _monomial(*(_atom("G", pair) for pair in pairing))
            out[mono] = out.get(mono, Fraction(0)) + 1
        return out
    raise ValueError(f"sixth-order reduction received {length} slots")


_PROFILE_GAPS: dict[str, dict[int, int]] = {
    "same_gap": {slot: 0 for slot in range(1, 7)},
    "two_adjacent": {1: 0, 2: 0, 3: 0, 4: 0, 5: 1, 6: 1},
    "three_adjacent": {1: 0, 2: 0, 3: 1, 4: 1, 5: 2, 6: 2},
}
_PROFILE_WEIGHTS = {
    "same_gap": Fraction(1, factorial(6)),
    # This includes (4,2) and its reflected (2,4) profile.
    "two_adjacent": Fraction(1, factorial(4) * factorial(2)) * 2,
    "three_adjacent": Fraction(1, factorial(2) ** 3),
}


def _block_layer_slot_sets(profile: str, block: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Layer-moment slot sets carried by one vertical cumulant block.

    The product of spins in each layer is reduced mod 2 automatically by the
    spin product itself.  Thus the moment of a slot set is the moment at the
    XOR of its physical site masks; no separate coincidence case is needed.
    """
    gaps = _PROFILE_GAPS[profile]
    block_set = set(block)
    if profile == "same_gap":
        return tuple(sorted(block_set)), tuple(sorted(block_set))
    if profile == "two_adjacent":
        first = tuple(sorted(slot for slot in block_set if gaps[slot] == 0))
        second = tuple(sorted(slot for slot in block_set if gaps[slot] == 1))
        return first, second, tuple(sorted(first + second))
    first = tuple(sorted(slot for slot in block_set if gaps[slot] == 0))
    middle = tuple(sorted(slot for slot in block_set if gaps[slot] == 1))
    last = tuple(sorted(slot for slot in block_set if gaps[slot] == 2))
    return first, last, tuple(sorted(first + middle)), tuple(sorted(middle + last))


def _block_connected_polynomial(profile: str, block: tuple[int, ...]) -> Polynomial:
    out: Polynomial = {(): Fraction(1)}
    for layer_slots in _block_layer_slot_sets(profile, block):
        out = _poly_mul(out, _moment_connected_expansion(layer_slots))
        if not out:
            break
    return out


def _profile_connected_polynomial(profile: str) -> tuple[Polynomial, int]:
    out: Polynomial = {}
    admissible = 0
    for partition in PARTITIONS_6:
        term: Polynomial = {(): Fraction(1)}
        for block in partition:
            term = _poly_mul(term, _block_connected_polynomial(profile, block))
            if not term:
                break
        if not term:
            continue
        admissible += 1
        out = _poly_add(out, term, _PROFILE_WEIGHTS[profile] * _mu(partition))
    return out, admissible


_S5_PERMUTATIONS = tuple(
    {1: 1, **dict(zip((2, 3, 4, 5, 6), permutation))}
    for permutation in permutations((2, 3, 4, 5, 6))
)


def _map_monomial(monomial: Monomial, mapping: dict[int, int]) -> Monomial:
    return _monomial(*(_atom(kind, (mapping[slot] for slot in slots)) for kind, slots in monomial))


def _canonical_monomial(monomial: Monomial) -> Monomial:
    return min(_map_monomial(monomial, mapping) for mapping in _S5_PERMUTATIONS)


def _connected_to_anchor(monomial: Monomial) -> bool:
    reached = {1}
    changed = True
    while changed:
        changed = False
        for _kind, slots in monomial:
            if reached.intersection(slots) and not set(slots).issubset(reached):
                reached.update(slots)
                changed = True
    return reached == set(range(1, 7))


def _even_slot_incidence(monomial: Monomial) -> bool:
    incidence = Counter(slot for _kind, slots in monomial for slot in slots)
    return all(incidence[slot] % 2 == 0 for slot in range(1, 7))


def _format_atom(atom: Atom) -> str:
    kind, slots = atom
    return f"{kind}({','.join(map(str, slots))})"


def _format_monomial(monomial: Monomial) -> str:
    return " \\cdot ".join(_format_atom(atom) for atom in monomial)


def _derive_connected_identity() -> tuple[list[dict[str, object]], dict[str, object]]:
    """Expand all 203 partitions and collect the exact 13-term G/U/W form."""
    raw: Polynomial = {}
    admissible: dict[str, int] = {}
    for profile in _PROFILE_GAPS:
        local, count = _profile_connected_polynomial(profile)
        raw = _poly_add(raw, local)
        admissible[profile] = count

    canonical: defaultdict[Monomial, Fraction] = defaultdict(Fraction)
    for monomial, coefficient in raw.items():
        canonical[_canonical_monomial(monomial)] += coefficient
    survivors = {mono: coefficient for mono, coefficient in canonical.items() if coefficient}

    if admissible != {"same_gap": 31, "two_adjacent": 11, "three_adjacent": 5}:
        raise AssertionError(f"unexpected admissible partition counts: {admissible}")
    if len(PARTITIONS_6) != 203 or len(raw) != 376 or len(survivors) != 13:
        raise AssertionError("sixth-order connected algebra did not reach its exact normal form")
    if not all(_connected_to_anchor(mono) and _even_slot_incidence(mono) for mono in survivors):
        raise AssertionError("a disconnected or odd-incidence term survived the cumulant reduction")

    rows = [
        {
            "coefficient": str(coefficient),
            "atoms": [{"kind": kind, "slots": list(slots)} for kind, slots in monomial],
            "latex": _format_monomial(monomial),
        }
        for monomial, coefficient in sorted(survivors.items(), key=lambda item: (len(item[0]), item[0]))
    ]
    expected_first = Fraction(1, 720)
    w_square = _monomial(_atom("W", (1, 2, 3, 4, 5, 6)), _atom("W", (1, 2, 3, 4, 5, 6)))
    if survivors.get(w_square) != expected_first:
        raise AssertionError("W^2 coefficient is not 1/720")
    certificate = {
        "bell_number_B6": len(PARTITIONS_6),
        "raw_monomial_count_before_S5_collection": len(raw),
        "canonical_structure_count": len(rows),
        "admissible_partition_counts": admissible,
        "all_survivors_anchor_connected": True,
        "all_survivors_even_incidence": True,
        "normal_form_sha256": hashlib.sha256(
            json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
    }
    return rows, certificate


# ---------------------------------------------------------------------------
# Independent exact re-expansion of the raw partition form by 2D FLM
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _profile_coefficients(order_n: int) -> tuple[tuple[tuple[int, ...], int], ...]:
    weights: Counter[tuple[int, ...]] = Counter()
    for partition in _set_partitions(tuple(range(order_n))):
        profile = tuple(sorted((len(block) for block in partition), reverse=True))
        weights[profile] += _mu(partition)
    return tuple(sorted(weights.items(), reverse=True))


@lru_cache(maxsize=None)
def _q_transform_weights(vertical_bonds: int) -> tuple[tuple[Fraction, ...], ...]:
    """n![q^n] cosh(q)^B tanh(q)^d for 0<=n,d<=6, exactly."""
    q_order = 6
    cosh_q: RationalSeries = (
        Fraction(1), Fraction(0), Fraction(1, 2), Fraction(0),
        Fraction(1, 24), Fraction(0), Fraction(1, 720),
    )
    tanh_q: RationalSeries = (
        Fraction(0), Fraction(1), Fraction(0), Fraction(-1, 3),
        Fraction(0), Fraction(2, 15), Fraction(0),
    )
    prefactor = _pow(cosh_q, vertical_bonds, q_order)
    rows: list[tuple[Fraction, ...]] = []
    for order_n in range(q_order + 1):
        row: list[Fraction] = []
        for vertical_degree in range(q_order + 1):
            factor = _mul(prefactor, _pow(tanh_q, vertical_degree, q_order), q_order)
            row.append(Fraction(factor[order_n]) * factorial(order_n))
        rows.append(tuple(row))
    return tuple(rows)


def _active_masks(parity: Sequence[IntegerSeries], max_weight: int = 6) -> tuple[tuple[tuple[int, int], ...], bool]:
    active: list[tuple[int, int]] = []
    odd_zero = True
    for mask, series in enumerate(parity):
        weight = mask.bit_count()
        if weight <= max_weight and any(series):
            active.append((mask, weight))
            if weight % 2:
                odd_zero = False
    return tuple(active), odd_zero


def _slab_subset_columns(
    parity: Sequence[IntegerSeries], layers: int, order: int = ORDER, max_weight: int = 6
) -> tuple[IntegerSeries, ...]:
    """Exact vertical-subset columns for an open h-layer slab.

    A gap-mask sequence S_1,...,S_(h-1) carries layer factors
    P[S_1] P[S_1 xor S_2] ... P[S_(h-1)], the raw 2D correlation
    factorization used by the closed form.
    """
    if not 1 <= layers <= MAX_LAYERS:
        raise ValueError("layers must lie between one and four")
    if layers == 1:
        columns = [tuple(_int_zero(order)) for _ in range(max_weight + 1)]
        columns[0] = tuple(parity[0])
        return tuple(columns)

    active, odd_zero = _active_masks(parity, max_weight)
    if not odd_zero:
        raise AssertionError("odd planar boundary mask unexpectedly survived")
    state: dict[tuple[int, int], IntegerSeries] = {
        (mask, weight): tuple(parity[mask]) for mask, weight in active
    }
    for _gap in range(1, layers - 1):
        updated: dict[tuple[int, int], list[int]] = {}
        for (previous, previous_weight), accumulated in state.items():
            for current, current_weight in active:
                total_weight = previous_weight + current_weight
                if total_weight > max_weight:
                    continue
                middle = parity[previous ^ current]
                if not any(middle):
                    continue
                contribution = _mul(accumulated, middle, order)
                if not any(contribution):
                    continue
                key = current, total_weight
                if key not in updated:
                    updated[key] = _int_zero(order)
                _add_scaled(updated[key], contribution)
        state = {key: tuple(value) for key, value in updated.items()}

    columns = [_int_zero(order) for _ in range(max_weight + 1)]
    for (last, vertical_degree), accumulated in state.items():
        _add_scaled(columns[vertical_degree], _mul(accumulated, parity[last], order))
    return tuple(tuple(column) for column in columns)


def _raw_moments(
    columns: Sequence[IntegerSeries], partition: IntegerSeries, layers: int, area: int, order: int = ORDER
) -> dict[int, RationalSeries]:
    denominator = _pow(partition, layers, order)
    normalized = tuple(
        tuple(Fraction(value) for value in _divide_unit(column, denominator, order))
        for column in columns
    )
    weights = _q_transform_weights(area * (layers - 1))
    moments: dict[int, RationalSeries] = {}
    for order_n in range(1, 7):
        values = _zero(order)
        for vertical_degree in range(7):
            if weights[order_n][vertical_degree]:
                _add_scaled(values, normalized[vertical_degree], weights[order_n][vertical_degree])
        moments[order_n] = tuple(values)
    return moments


def _cumulant(moments: dict[int, RationalSeries], order_n: int, order: int = ORDER) -> RationalSeries:
    values = _zero(order)
    unit: RationalSeries = (Fraction(1),) + (Fraction(0),) * order
    for profile, coefficient in _profile_coefficients(order_n):
        term: Series = unit
        for size in profile:
            term = _mul(term, moments[size], order)
        _add_scaled(values, term, coefficient)
    return tuple(values)


@lru_cache(maxsize=None)
def _rectangle_components(width: int, height: int, order: int = ORDER) -> tuple[tuple[str, RationalSeries], ...]:
    first, second = sorted((width, height))
    parity = _boundary_polynomials(first, second, order)
    area = first * second
    finite: dict[int, dict[int, RationalSeries]] = {
        1: {order_n: tuple(_zero(order)) for order_n in range(1, 7)}
    }
    for layers in range(2, MAX_LAYERS + 1):
        raw = _raw_moments(
            _slab_subset_columns(parity, layers, order), parity[0], layers, area, order
        )
        finite[layers] = {
            order_n: tuple(value / factorial(order_n) for value in _cumulant(raw, order_n, order))
            for order_n in range(1, 7)
        }

    vertical_weights: dict[tuple[int, int], RationalSeries] = {}
    for order_n in range(1, 7):
        for layers in range(1, MAX_LAYERS + 1):
            values = list(finite[layers][order_n])
            for lower in range(1, layers):
                _add_scaled(values, vertical_weights[(order_n, lower)], -(layers - lower + 1))
            vertical_weights[(order_n, layers)] = tuple(values)
    return tuple(
        sorted(
            (f"c{order_n}_h{layers}", vertical_weights[(order_n, layers)])
            for order_n in range(1, 7)
            for layers in range(1, MAX_LAYERS + 1)
        )
    )


def _plane_shapes(order: int = ORDER) -> tuple[tuple[int, int], ...]:
    span = order // 2
    return tuple(
        sorted(
            (
                (width, height)
                for width in range(1, span + 2)
                for height in range(1, span + 2)
                if (width - 1) + (height - 1) <= span
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )


def _cumulant_flm(order: int = ORDER) -> tuple[
    dict[str, RationalSeries],
    dict[tuple[int, int], dict[str, RationalSeries]],
    tuple[tuple[int, int], ...],
]:
    shapes = _plane_shapes(order)
    names = tuple(f"c{n}_h{layers}" for n in range(1, 7) for layers in range(1, MAX_LAYERS + 1))
    weights: dict[tuple[int, int], dict[str, RationalSeries]] = {}
    bulk = {name: _zero(order) for name in names}
    for width, height in shapes:
        finite = dict(_rectangle_components(width, height, order))
        local = {name: list(finite[name]) for name in names}
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                for name in names:
                    _add_scaled(local[name], subweight[name], -placements)
        exact = {name: tuple(values) for name, values in local.items()}
        weights[(width, height)] = exact
        for name in names:
            _add_scaled(bulk[name], exact[name])
    return {name: tuple(values) for name, values in bulk.items()}, weights, shapes


# ---------------------------------------------------------------------------
# Artifact helpers and verification
# ---------------------------------------------------------------------------
def _fractions(values: Sequence[str]) -> RationalSeries:
    return tuple(Fraction(value) for value in values)


def _strings(values: Sequence[int | Fraction]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _sum_series(*terms: Sequence[Fraction]) -> RationalSeries:
    return tuple(sum((term[degree] for term in terms), Fraction(0)) for degree in range(ORDER + 1))


def _record(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def main() -> None:
    started = time.monotonic()
    checks: list[dict[str, object]] = []
    identity, identity_certificate = _derive_connected_identity()
    _record(
        checks,
        "complete G/U/W sixth-cumulant reduction",
        identity_certificate["canonical_structure_count"] == 13
        and identity_certificate["all_survivors_anchor_connected"]
        and identity_certificate["all_survivors_even_incidence"],
        "[THEOREM] all 203 partitions reduce exactly to 13 anchored connected G/U/W lattice sums",
    )

    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))["data"]
    bulk, weights, shapes = _cumulant_flm(ORDER)
    c2_q = bulk["c2_h2"]
    c4_q = _sum_series(bulk["c4_h2"], bulk["c4_h3"])
    c6_h2 = bulk["c6_h2"]
    c6_h3 = bulk["c6_h3"]
    c6_h4 = bulk["c6_h4"]
    c6_q = _sum_series(c6_h2, c6_h3, c6_h4)

    expected_q = _fractions(primary["c6_direct_coupling_q"]["full_infinite_stack"])
    expected_h2 = _fractions(primary["c6_direct_coupling_q"]["same_gap_height2"])
    expected_h3 = _fractions(primary["c6_direct_coupling_q"]["two_adjacent_gaps_height3"])
    expected_h4 = _fractions(primary["c6_direct_coupling_q"]["three_adjacent_gaps_height4"])
    _record(
        checks,
        "exact c6(q) prefix through v^12",
        c6_q == expected_q,
        "[COMPUTATION] independent 2D boundary-mask cumulant FLM matches every stored c6(q) coefficient",
    )
    _record(
        checks,
        "exact vertical-profile decomposition",
        c6_h2 == expected_h2 and c6_h3 == expected_h3 and c6_h4 == expected_h4,
        "[COMPUTATION] profiles (6), (4+2/2+4), and (2+2+2) match height 2/3/4 certificate rows",
    )

    expected_c2 = _fractions(primary["reproduced_lower_orders"]["c2_total_q_or_w"])
    expected_c4 = _fractions(primary["reproduced_lower_orders"]["c4_direct_q"])
    _record(
        checks,
        "lower cumulants required by q-to-w conversion",
        c2_q == expected_c2 and c4_q == expected_c4,
        "[COMPUTATION] the same raw correlation engine reproduces c2(q) and c4(q) through v^12",
    )

    c6_w_total = tuple(
        c6_q[degree] + Fraction(4, 3) * c4_q[degree] + Fraction(23, 45) * c2_q[degree]
        for degree in range(ORDER + 1)
    )
    c6_w_residual = list(c6_w_total)
    c6_w_residual[0] -= Fraction(1, 6)
    expected_w_total = _fractions(primary["c6_wave_variable_w_tanh_Kz"]["total"])
    expected_w_residual = _fractions(primary["c6_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"])
    _record(
        checks,
        "exact q-to-w conversion through v^12",
        c6_w_total == expected_w_total and tuple(c6_w_residual) == expected_w_residual,
        "[COMPUTATION] c6,w=c6,q+(4/3)c4,q+(23/45)c2,q, with the 1/6 log-cosh constant removed",
    )

    _record(
        checks,
        "finite-support and parity controls",
        len(shapes) == 28
        and all(c6_q[degree] == 0 for degree in range(1, ORDER + 1, 2))
        and c6_q[0] == Fraction(1, 45),
        "[LEMMA] every surviving connected normal-form term has even slot incidence; span six suffices through v^12",
    )
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError("one or more c6 closed-form checks failed")

    weights_rows = [
        {
            "rectangle": list(shape),
            "c6_h2": _strings(series["c6_h2"]),
            "c6_h3": _strings(series["c6_h3"]),
            "c6_h4": _strings(series["c6_h4"]),
        }
        for shape, series in sorted(weights.items())
    ]
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "actual_executable": sys.executable,
            "method": (
                "complete Bell(6) partition reduction to connected 2D G/U/W moments, then independent "
                "open-rectangle boundary-mask layer transfer, q-cumulants, vertical Mobius inversion, and planar FLM"
            ),
        },
        "data": {
            "claim": (
                "[THEOREM] c6(q) is the displayed finite linear combination of 13 anchored products "
                "of 2D G, connected U4, and connected W6 correlations; [COMPUTATION] its exact v^12 "
                "prefix is certified below."
            ),
            "variables": {
                "v": "tanh(K), the in-plane high-temperature variable",
                "q": "K_z, direct interlayer coupling for the cumulant expansion",
                "w": "tanh(K_z)",
                "slot_convention": (
                    "x1=0 and x2,...,x6 are summed over Z^2. M_I is the layer moment of the product "
                    "of spins at slots I; repeated site positions cancel mod 2 inside that product."
                ),
            },
            "connected_correlation_identity": {
                "tag": "[THEOREM]",
                "definition": (
                    "For every even slot set I, M_I=sum over even set partitions rho of I of product_C K_C, "
                    "where K_ij=G_ij, K_ijkl=U_ijkl, and K_ijklmn=W_ijklmn. "
                    "The 13 rows below give c6,q=sum_{x2,...,x6} coefficient*product(atoms)."
                ),
                "formula": identity,
                "certificate": identity_certificate,
            },
            "vertical_partition_formula": {
                "tag": "[THEOREM]",
                "statement": (
                    "c6,q=sum_{x2,...,x6}[kappa6(Y0,...,Y5)/6! + "
                    "kappa(Y0,Y1,Y2,Y3,Z4,Z5)/(4!2!)*2 + "
                    "kappa(Y0,Y1,Z2,Z3,T4,T5)/(2!2!2!)], with Y,Z,T on consecutive gaps. "
                    "Each block factor is M_B^2 (one gap), M_B0 M_B1 M_(B0 union B1) (two gaps), or "
                    "M_B0 M_B2 M_(B0 union B1) M_(B1 union B2) (three gaps)."
                ),
                "profile_mobius_weight_sums": {
                    "6": {"6": 1, "4+2": -15, "2+2+2": 30},
                    "4+2": {"6": 1, "4+2": -7, "2+2+2": 6},
                    "2+4": {"6": 1, "4+2": -7, "2+2+2": 6},
                    "2+2+2": {"6": 1, "4+2": -3, "2+2+2": 2},
                },
                "placement_weights": {"6": "1/720", "4+2_and_2+4_combined": "1/24", "2+2+2": "1/8"},
            },
            "series_verification": {
                "tag": "[COMPUTATION]",
                "v_order": ORDER,
                "c2_direct_q": _strings(c2_q),
                "c4_direct_q": _strings(c4_q),
                "c6_direct_q": _strings(c6_q),
                "c6_direct_q_by_vertical_profile": {
                    "same_gap_height2": _strings(c6_h2),
                    "two_adjacent_gaps_height3": _strings(c6_h3),
                    "three_adjacent_gaps_height4": _strings(c6_h4),
                },
                "c6_total_w": _strings(c6_w_total),
                "c6_residual_w": _strings(c6_w_residual),
                "rectangle_count": len(shapes),
                "rectangles": [list(shape) for shape in shapes],
                "rectangle_mobius_weights": weights_rows,
            },
            "scope": (
                "[UNRESOLVED] The theorem is a coefficientwise formal high-temperature identity. It does not "
                "evaluate the infinite G/U/W lattice sums at criticality, prove a convergence radius, or solve 3D Ising. "
                "The 2D W6 objects have standard fermionic/Pfaffian evaluation in principle but no scalar closed "
                "evaluation of the required infinite lattice sums is supplied here."
            ),
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    elapsed = time.monotonic() - started
    print(f"PASS e99_c6_closed ({elapsed:.1f}s)")


if __name__ == "__main__":
    main()
