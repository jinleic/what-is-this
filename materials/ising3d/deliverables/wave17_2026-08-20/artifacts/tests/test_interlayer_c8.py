"""Clean-room exact checks for the eighth interlayer coefficient.

This file intentionally does not import ``experiments.e126_interlayer_c8``.  Its
first route builds raw open-box even-subgraph polynomials with an independent
boundary-mask transfer and performs the [w^8] logarithm plus 3D rectangular
finite-lattice inversion.  Its second route starts from the repository's
independent spin density of states for every box through v^6.
"""

from __future__ import annotations

import hashlib
import json
import resource
import time
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
from typing import Iterable, Sequence

from ising.interlayer import anisotropic_box_even_subgraph

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "interlayer" / "c8_series.json"
V_ORDER = 12
W_ORDER = 8
TEST_WALL_BUDGET_SECONDS = 900   # CPU-seconds budget (process_time), not wall time

EXPECTED_C8_RESIDUAL = tuple(
    Fraction(value)
    for value in (
        "0",
        "0",
        "2",
        "0",
        "243",
        "0",
        "13594",
        "0",
        "1046015/2",
        "0",
        "14894834",
        "0",
        "341377607",
    )
)
EXPECTED_PATTERNS = (
    (4,),
    (3, 1),
    (2, 2),
    (2, 1, 1),
    (1, 3),
    (1, 2, 1),
    (1, 1, 2),
    (1, 1, 1, 1),
)
EXPECTED_ADMISSIBLE = (379, 107, 83, 35, 107, 35, 35, 15)
EXPECTED_RAW_MONOMIALS = (379, 107, 83, 22, 107, 31, 22, 8)
EXPECTED_CONNECTED_MONOMIALS = (56932, 8718, 4458, 270, 8718, 621, 270, 27)

Series = tuple[int | Fraction, ...]
Atom = tuple[str, tuple[int, ...]]
Monomial = tuple[Atom, ...]
Polynomial = dict[Monomial, Fraction]
ATOM_KIND = {2: "G", 4: "U", 6: "W", 8: "W8"}


def zeros(order: int) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def add_scaled(target: list[int | Fraction], source: Sequence[int | Fraction], scale: int | Fraction = 1) -> None:
    for degree, value in enumerate(source[: len(target)]):
        if value:
            target[degree] += scale * value


def multiply(left: Sequence[int | Fraction], right: Sequence[int | Fraction], order: int) -> list[int | Fraction]:
    out: list[int | Fraction] = [0] * (order + 1)
    for i, first in enumerate(left[: order + 1]):
        if not first:
            continue
        for j, second in enumerate(right[: order + 1 - i]):
            if second:
                out[i + j] += first * second
    return out


def divide(numerator: Sequence[int | Fraction], denominator: Sequence[int | Fraction], order: int) -> tuple[Fraction, ...]:
    if not denominator[0]:
        raise AssertionError("zero constant term in formal quotient")
    quotient = zeros(order)
    for degree in range(order + 1):
        value = Fraction(numerator[degree]) if degree < len(numerator) else Fraction(0)
        for lower in range(degree):
            if quotient[lower]:
                value -= quotient[lower] * Fraction(denominator[degree - lower])
        quotient[degree] = value / Fraction(denominator[0])
    return tuple(quotient)


def square_bonds(width: int, height: int) -> tuple[tuple[int, int], ...]:
    edges: list[tuple[int, int]] = []
    for y in range(height):
        for x in range(width):
            site = y * width + x
            if x + 1 < width:
                edges.append((site, site + 1))
            if y + 1 < height:
                edges.append((site, site + width))
    assert len(edges) == (width - 1) * height + width * (height - 1)
    return tuple(edges)


@lru_cache(maxsize=None)
def boundary_masks(width: int, height: int, order: int) -> tuple[tuple[int, ...], ...]:
    """Raw edge-subset DP for P_S(v), indexed by its odd-boundary mask."""

    sites = width * height
    values = [[0] * (order + 1) for _ in range(1 << sites)]
    values[0][0] = 1
    for first, second in square_bonds(width, height):
        toggle = (1 << first) | (1 << second)
        for degree in range(order - 1, -1, -1):
            for mask in range(1 << sites):
                if values[mask][degree]:
                    values[mask ^ toggle][degree + 1] += values[mask][degree]
    assert all(not any(values[mask]) for mask in range(1 << sites) if mask.bit_count() % 2)
    return tuple(tuple(row) for row in values)


def rectangle_maps(width: int, height: int) -> tuple[tuple[int, ...], ...]:
    """All graph automorphisms used only to accelerate a raw mask sum."""

    def site(x: int, y: int) -> int:
        return y * width + x

    coordinate_maps = [
        lambda x, y: (x, y),
        lambda x, y: (width - 1 - x, y),
        lambda x, y: (x, height - 1 - y),
        lambda x, y: (width - 1 - x, height - 1 - y),
    ]
    if width == height:
        coordinate_maps.extend(
            (
                lambda x, y: (y, x),
                lambda x, y: (y, width - 1 - x),
                lambda x, y: (height - 1 - y, x),
                lambda x, y: (height - 1 - y, width - 1 - x),
            )
        )
    edges = set(square_bonds(width, height))
    rows: list[tuple[int, ...]] = []
    for transform in coordinate_maps:
        mapping = tuple(site(*transform(x, y)) for y in range(height) for x in range(width))
        transformed = {
            (min(mapping[first], mapping[second]), max(mapping[first], mapping[second]))
            for first, second in edges
        }
        assert transformed == edges
        if mapping not in rows:
            rows.append(mapping)
    return tuple(rows)


def permute_mask(mask: int, mapping: Sequence[int]) -> int:
    result = 0
    for source, target in enumerate(mapping):
        if (mask >> source) & 1:
            result |= 1 << target
    return result


@lru_cache(maxsize=None)
def orbits(width: int, height: int) -> tuple[tuple[int, ...], tuple[int, ...], dict[int, tuple[int, ...]]]:
    maps = rectangle_maps(width, height)
    count = 1 << (width * height)
    representatives = [0] * count
    sizes = [0] * count
    members: dict[int, list[int]] = defaultdict(list)
    for mask in range(count):
        orbit = {permute_mask(mask, mapping) for mapping in maps}
        representative = min(orbit)
        representatives[mask] = representative
        sizes[mask] = len(orbit)
        members[representative].append(mask)
    return tuple(representatives), tuple(sizes), {key: tuple(value) for key, value in members.items()}


def active_masks(parity: Sequence[Series], max_weight: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (mask, mask.bit_count())
        for mask, series in enumerate(parity)
        if mask.bit_count() <= max_weight and any(series)
    )


def plain_slab_columns(parity: Sequence[Series], layers: int, order: int, max_weight: int) -> tuple[Series, ...]:
    active = active_masks(parity, max_weight)
    columns = [[0] * (order + 1) for _ in range(max_weight + 1)]
    if layers == 1:
        columns[0] = list(parity[0])
        return tuple(tuple(row) for row in columns)
    state: dict[tuple[int, int], list[int | Fraction]] = {
        (mask, weight): list(parity[mask]) for mask, weight in active
    }
    for _ in range(1, layers - 1):
        updated: dict[tuple[int, int], list[int | Fraction]] = {}
        for (previous, used), accumulated in state.items():
            for current, weight in active:
                if used + weight > max_weight:
                    continue
                middle = parity[previous ^ current]
                if not any(middle):
                    continue
                term = multiply(accumulated, middle, order)
                if not any(term):
                    continue
                key = (current, used + weight)
                if key not in updated:
                    updated[key] = term
                else:
                    add_scaled(updated[key], term)
        state = updated
    for (last, degree), accumulated in state.items():
        add_scaled(columns[degree], multiply(accumulated, parity[last], order))
    return tuple(tuple(row) for row in columns)


def orbit_slab_columns(width: int, height: int, layers: int, order: int, max_weight: int) -> tuple[Series, ...]:
    """Exact direct-w transfer, summed over every vertical-edge mask orbit."""

    parity = boundary_masks(width, height, order)
    active = active_masks(parity, max_weight)
    columns = [[0] * (order + 1) for _ in range(max_weight + 1)]
    if layers == 1:
        columns[0] = list(parity[0])
        return tuple(tuple(row) for row in columns)
    representative, orbit_size, members = orbits(width, height)
    targets: dict[int, set[int]] = defaultdict(set)
    state: dict[tuple[int, int], list[int | Fraction]] = {}
    for mask, weight in active:
        rep = representative[mask]
        targets[weight].add(rep)
        state[(rep, weight)] = list(parity[rep])
        assert parity[mask] == parity[rep]
    for _ in range(1, layers - 1):
        updated: dict[tuple[int, int], list[int | Fraction]] = {}
        for (source, used), accumulated in state.items():
            for mask in members[source]:
                for target_weight, target_reps in targets.items():
                    if used + target_weight > max_weight:
                        continue
                    for target in target_reps:
                        middle = parity[mask ^ target]
                        if not any(middle):
                            continue
                        term = multiply(accumulated, middle, order)
                        if not any(term):
                            continue
                        key = (target, used + target_weight)
                        if key not in updated:
                            updated[key] = term
                        else:
                            add_scaled(updated[key], term)
        state = updated
    for (last, degree), accumulated in state.items():
        add_scaled(columns[degree], multiply(accumulated, parity[last], order), orbit_size[last])
    return tuple(tuple(row) for row in columns)


def logarithmic_columns(columns: Sequence[Series], order: int) -> dict[int, tuple[Fraction, ...]]:
    """Compute [w^d] log(P/P_0) from raw open-box even-subgraph columns."""

    ratios = {degree: divide(columns[degree], columns[0], order) for degree in range(1, W_ORDER + 1)}
    logs: dict[int, tuple[Fraction, ...]] = {}
    for degree in range(1, W_ORDER + 1):
        value = [Fraction(degree) * coefficient for coefficient in ratios[degree]]
        for lower in range(1, degree):
            add_scaled(value, multiply(logs[lower], ratios[degree - lower], order), -lower)
        logs[degree] = tuple(coefficient / degree for coefficient in value)
    return logs


def plane_shapes(order: int) -> tuple[tuple[int, int], ...]:
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


def direct_wlog_flm(order: int, max_layers: int) -> tuple[tuple[Fraction, ...], dict[int, tuple[Fraction, ...]], dict[tuple[int, int, int], tuple[Fraction, ...]]]:
    """Independent direct-w open-slab FLM for the c8 residual."""

    boxes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in plane_shapes(order)
                for layers in range(1, max_layers + 1)
            ),
            key=lambda box: (sum(box), box),
        )
    )
    weights: dict[tuple[int, int, int], tuple[Fraction, ...]] = {}
    bulk = zeros(order)
    by_height = {layers: zeros(order) for layers in range(1, max_layers + 1)}
    for width, height, layers in boxes:
        columns = orbit_slab_columns(width, height, layers, order, W_ORDER)
        logs = logarithmic_columns(columns, order)
        assert all(not any(logs[degree]) for degree in (1, 3, 5, 7))
        local = list(logs[8])
        for (sub_width, sub_height, sub_layers), subweight in weights.items():
            if sub_width <= width and sub_height <= height and sub_layers <= layers:
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sub_layers + 1)
                )
                add_scaled(local, subweight, -placements)
        exact = tuple(Fraction(value) for value in local)
        weights[(width, height, layers)] = exact
        add_scaled(bulk, exact)
        add_scaled(by_height[layers], exact)
    return tuple(bulk), {layers: tuple(values) for layers, values in by_height.items()}, weights


def spin_dos_wlog_flm(order: int) -> tuple[tuple[Fraction, ...], dict[int, tuple[Fraction, ...]], dict[tuple[int, int, int], tuple[Fraction, ...]]]:
    """Independent spin-DOS/3D-FLM check for every coefficient through v^6."""

    boxes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in plane_shapes(order)
                for layers in range(1, 6)
            ),
            key=lambda box: (sum(box), box),
        )
    )
    weights: dict[tuple[int, int, int], tuple[Fraction, ...]] = {}
    bulk = zeros(order)
    by_height = {layers: zeros(order) for layers in range(1, 6)}
    for width, height, layers in boxes:
        first, second = sorted((width, height))
        polynomial = anisotropic_box_even_subgraph((first, second, layers), order, W_ORDER)
        columns = tuple(
            tuple(polynomial[v_degree][w_degree] for v_degree in range(order + 1))
            for w_degree in range(W_ORDER + 1)
        )
        logs = logarithmic_columns(columns, order)
        assert all(not any(logs[degree]) for degree in (1, 3, 5, 7))
        local = list(logs[8])
        for (sub_width, sub_height, sub_layers), subweight in weights.items():
            if sub_width <= width and sub_height <= height and sub_layers <= layers:
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sub_layers + 1)
                )
                add_scaled(local, subweight, -placements)
        exact = tuple(Fraction(value) for value in local)
        weights[(width, height, layers)] = exact
        add_scaled(bulk, exact)
        add_scaled(by_height[layers], exact)
    return tuple(bulk), {layers: tuple(values) for layers, values in by_height.items()}, weights

def spin_dos_height6_control_weights() -> dict[tuple[int, int, int], dict[int, tuple[Fraction, ...]]]:
    """Rebuild every weight in the producer's explicitly scoped height-six family."""

    shapes = ((1, 1), (2, 1), (3, 1), (4, 1), (2, 2))
    boxes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in shapes
                for layers in range(1, 7)
            ),
            key=lambda box: (sum(box), box),
        )
    )
    degrees = (2, 4, 6, 8)
    weights: dict[tuple[int, int, int], dict[int, tuple[Fraction, ...]]] = {}
    for width, height, layers in boxes:
        first, second = sorted((width, height))
        polynomial = anisotropic_box_even_subgraph((first, second, layers), 6, W_ORDER)
        columns = tuple(
            tuple(polynomial[v_degree][w_degree] for v_degree in range(7))
            for w_degree in range(W_ORDER + 1)
        )
        logs = logarithmic_columns(columns, 6)
        assert all(not any(logs[degree]) for degree in (1, 3, 5, 7))
        local = {degree: list(logs[degree]) for degree in degrees}
        for (sub_width, sub_height, sub_layers), subweight in weights.items():
            if sub_width <= width and sub_height <= height and sub_layers <= layers:
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sub_layers + 1)
                )
                for degree in degrees:
                    add_scaled(local[degree], subweight[degree], -placements)
        weights[(width, height, layers)] = {
            degree: tuple(Fraction(value) for value in values)
            for degree, values in local.items()
        }
    return weights


@lru_cache(maxsize=None)
def set_partitions(slots: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    if not slots:
        return ((),)
    first, tail = slots[0], slots[1:]
    result: list[tuple[tuple[int, ...], ...]] = []
    for partition in set_partitions(tail):
        result.append(((first,),) + partition)
        for index, block in enumerate(partition):
            result.append(partition[:index] + ((first,) + block,) + partition[index + 1 :])
    return tuple(result)


@lru_cache(maxsize=None)
def even_partitions(slots: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    if not slots:
        return ((),)
    first, tail = slots[0], slots[1:]
    result: list[tuple[tuple[int, ...], ...]] = []
    for size in range(2, len(slots) + 1, 2):
        for partners in combinations(tail, size - 1):
            block = (first,) + partners
            remaining = tuple(slot for slot in tail if slot not in partners)
            for rest in even_partitions(remaining):
                result.append((block,) + rest)
    return tuple(result)


@lru_cache(maxsize=None)
def moment_to_connected(slots: tuple[int, ...]) -> tuple[tuple[Monomial, Fraction], ...]:
    if not slots:
        return (((), Fraction(1)),)
    result: Polynomial = {}
    for partition in even_partitions(slots):
        monomial = tuple(
            sorted((ATOM_KIND[len(block)], tuple(sorted(block))) for block in partition)
        )
        result[monomial] = result.get(monomial, Fraction(0)) + 1
    return tuple(result.items())


def gap_map(pattern: Sequence[int]) -> dict[int, int]:
    result: dict[int, int] = {}
    slot = 1
    for gap, count_half in enumerate(pattern):
        for _ in range(2 * count_half):
            result[slot] = gap
            slot += 1
    return result


def block_layer_slots(gaps: dict[int, int], block: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    per_gap: dict[int, list[int]] = defaultdict(list)
    for slot in block:
        per_gap[gaps[slot]].append(slot)
    return tuple(
        tuple(sorted(per_gap.get(layer - 1, []) + per_gap.get(layer, [])))
        for layer in range(max(gaps.values()) + 2)
        if per_gap.get(layer - 1, []) or per_gap.get(layer, [])
    )


def indexed_block_layer_slots(
    gaps: dict[int, int], block: Sequence[int]
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    per_gap: dict[int, list[int]] = defaultdict(list)
    for slot in block:
        per_gap[gaps[slot]].append(slot)
    return tuple(
        (layer, tuple(sorted(per_gap.get(layer - 1, []) + per_gap.get(layer, []))))
        for layer in range(max(gaps.values()) + 2)
        if per_gap.get(layer - 1, []) or per_gap.get(layer, [])
    )


def gap_even(gaps: dict[int, int], block: Sequence[int]) -> bool:
    return all(count % 2 == 0 for count in Counter(gaps[slot] for slot in block).values())


def normalized_gap_multisets(count: int) -> tuple[tuple[int, ...], ...]:
    rows: list[tuple[int, ...]] = []

    def walk(prefix: tuple[int, ...]) -> None:
        if len(prefix) == count:
            rows.append(prefix)
            return
        if not prefix:
            walk((0,))
            return
        for value in range(prefix[-1], count):
            walk(prefix + (value,))

    walk(())
    return tuple(rows)


def literal_profile_moments(
    gap_values: Sequence[int],
) -> dict[tuple[tuple[int, tuple[int, ...]], ...], int]:
    """Literal layer-labelled Bell polynomial; no profile predicate is used."""

    gaps = {slot + 1: value for slot, value in enumerate(gap_values)}
    raw: dict[tuple[tuple[int, tuple[int, ...]], ...], int] = {}
    for partition in set_partitions(tuple(range(1, len(gap_values) + 1))):
        if not all(gap_even(gaps, block) for block in partition):
            continue
        coefficient = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        monomial = tuple(
            sorted(
                factor
                for block in partition
                for factor in indexed_block_layer_slots(gaps, block)
            )
        )
        raw[monomial] = raw.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in raw.items() if coefficient}


def placement_from_ordered_assignments(pattern: Sequence[int]) -> tuple[Fraction, int]:
    labels = tuple(gap for gap, half_count in enumerate(pattern) for _ in range(2 * half_count))
    ordered_assignments = len(set(permutations(labels)))
    return Fraction(ordered_assignments, factorial(len(labels))), ordered_assignments


def polynomial_digest(polynomials: dict[str, Polynomial]) -> str:
    digest = hashlib.sha256()
    for name in sorted(polynomials):
        digest.update(name.encode("utf-8"))
        for monomial, coefficient in sorted(polynomials[name].items()):
            digest.update(repr((monomial, coefficient)).encode("utf-8"))
    return digest.hexdigest()


def profile_inventory(
    pattern: Sequence[int], placement: Fraction
) -> tuple[
    int,
    int,
    int,
    set[str],
    Polynomial,
    Counter[tuple[str, ...]],
    Counter[tuple[int, ...]],
]:
    """A separately written, placement-scaled Bell/XOR expansion for one c8 row."""

    gaps = gap_map(pattern)
    raw: dict[tuple[tuple[int, tuple[int, ...]], ...], Fraction] = {}
    connected: Polynomial = {}
    admissible = 0
    moebius_by_class: Counter[tuple[int, ...]] = Counter()
    for partition in set_partitions(tuple(range(1, 9))):
        if not all(gap_even(gaps, block) for block in partition):
            continue
        admissible += 1
        primitive = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        moebius_by_class[
            tuple(sorted((len(block) for block in partition), reverse=True))
        ] += primitive
        moebius = placement * primitive
        raw_monomial: list[tuple[int, tuple[int, ...]]] = []
        terms: Polynomial = {(): moebius}
        for block in partition:
            for slots in block_layer_slots(gaps, block):
                raw_monomial.append((0, slots))
                expanded = moment_to_connected(slots)
                next_terms: Polynomial = {}
                for monomial, coefficient in terms.items():
                    for additional, other in expanded:
                        key = tuple(sorted(monomial + additional))
                        next_terms[key] = next_terms.get(key, Fraction(0)) + coefficient * other
                terms = {key: value for key, value in next_terms.items() if value}
        raw_key = tuple(sorted(raw_monomial))
        raw[raw_key] = raw.get(raw_key, Fraction(0)) + moebius
        for monomial, coefficient in terms.items():
            connected[monomial] = connected.get(monomial, Fraction(0)) + coefficient
    raw = {key: value for key, value in raw.items() if value}
    connected = {key: value for key, value in connected.items() if value}
    kinds = {kind for monomial in connected for kind, _slots in monomial}
    histogram = Counter(
        tuple(sorted(kind for kind, _slots in monomial)) for monomial in connected
    )
    return (
        admissible,
        len(raw),
        len(connected),
        kinds,
        connected,
        histogram,
        moebius_by_class,
    )

def assert_connected_structure(polynomial: Polynomial) -> None:
    for monomial in polynomial:
        incidences = Counter(slot for _kind, slots in monomial for slot in slots)
        assert all(incidences[slot] % 2 == 0 for slot in range(1, 9))
        reached = {1}
        changed = True
        while changed:
            changed = False
            for _kind, slots in monomial:
                if reached.intersection(slots) and not set(slots).issubset(reached):
                    reached.update(slots)
                    changed = True
        assert reached == set(range(1, 9))


def fractions(values: Iterable[str]) -> tuple[Fraction, ...]:
    return tuple(Fraction(value) for value in values)


def main() -> None:
    # CPU-time budget: on this shared machine wall clock tracks unrelated
    # background load, which must not change the verdict.
    started = time.process_time()

    # The symmetry optimization is only an optimization: this catches a
    # dropped orientation, orbit-size, or vertical-transition multiplicity.
    for width, height in ((1, 1), (2, 1), (3, 1), (2, 2), (3, 2)):
        parity = boundary_masks(width, height, V_ORDER)
        for layers in range(2, 5):
            assert plain_slab_columns(parity, layers, V_ORDER, W_ORDER) == orbit_slab_columns(
                width, height, layers, V_ORDER, W_ORDER
            )

    direct_bulk, direct_by_height, direct_weights = direct_wlog_flm(V_ORDER, 6)
    assert direct_bulk == EXPECTED_C8_RESIDUAL
    assert direct_bulk[2] == 2
    assert all(direct_by_height[1][degree] == 0 for degree in range(V_ORDER + 1))
    assert all(direct_by_height[6][degree] == 0 for degree in range(V_ORDER + 1))
    assert len(direct_weights) == 28 * 6

    # This route begins with raw 3D spin configurations, not boundary masks.
    spin_bulk, spin_by_height, spin_weights = spin_dos_wlog_flm(6)
    assert spin_bulk == EXPECTED_C8_RESIDUAL[:7]
    assert spin_by_height[2] == tuple(Fraction(0) for _ in range(7))
    assert spin_by_height[3] == tuple(Fraction(value) for value in (0, 0, 0, 0, -5, 0, -212))
    assert spin_by_height[4] == tuple(Fraction(value) for value in (0, 0, 0, 0, -40, 0, -1696))
    assert spin_by_height[5] == tuple(Fraction(value) for value in (0, 0, 2, 0, 288, 0, 15502))
    assert len(spin_weights) == 50

    control_weights = spin_dos_height6_control_weights()
    control_boxes = ((1, 1, 6), (2, 1, 6), (3, 1, 6), (4, 1, 6), (2, 2, 6))
    assert set(control_boxes).issubset(control_weights)
    assert all(
        not any(control_weights[box][degree])
        for box in control_boxes
        for degree in (2, 4, 6, 8)
    )

    # The literal, layer-labelled Bell/moment polynomial is evaluated for
    # every normalized eight-gap multiset; selection is not pre-filtered.
    all_gap_multisets = normalized_gap_multisets(8)
    assert len(all_gap_multisets) == 3432
    literal_nonzero = {
        values for values in all_gap_multisets if literal_profile_moments(values)
    }
    expected_gap_multisets = {
        tuple(sorted(gap_map(pattern).values())) for pattern in EXPECTED_PATTERNS
    }
    assert literal_nonzero == expected_gap_multisets
    partitions = set_partitions(tuple(range(1, 9)))
    assert len(partitions) == 4140
    total_connected = 0
    polynomials: dict[str, Polynomial] = {}
    row_reports = []
    for index, pattern in enumerate(EXPECTED_PATTERNS):
        placement, ordered_assignments = placement_from_ordered_assignments(pattern)
        (
            admissible,
            raw_count,
            connected_count,
            kinds,
            polynomial,
            histogram,
            moebius_by_class,
        ) = profile_inventory(pattern, placement)
        assert admissible == EXPECTED_ADMISSIBLE[index]
        assert raw_count == EXPECTED_RAW_MONOMIALS[index]
        assert connected_count == EXPECTED_CONNECTED_MONOMIALS[index]
        assert kinds <= {"G", "U", "W", "W8"}
        assert_connected_structure(polynomial)
        polynomials[str(pattern)] = polynomial
        row_reports.append(
            (
                pattern,
                placement,
                ordered_assignments,
                {",".join(kinds): count for kinds, count in sorted(histogram.items())},
                {
                    "+".join(map(str, block_class)): str(value)
                    for block_class, value in sorted(moebius_by_class.items())
                },
                polynomial_digest({str(pattern): polynomial}),
            )
        )
        total_connected += connected_count
    assert total_connected == 80014

    # Only after all decisive quantities have been independently rebuilt is
    # the producer artifact checked as a report of those quantities.
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    data = payload["data"]
    assert payload["provenance"]["script"] == "experiments/e126_interlayer_c8.py"
    assert all(check["passed"] for check in payload["checks"])
    assert fractions(data["c8_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"]) == direct_bulk
    assert fractions(data["second_route_w_log"]["residual_w8"]) == direct_bulk
    assert fractions(data["independent_spin_dos_route"]["residual_w8"]) == spin_bulk
    assert data["c8_wave_variable_w_tanh_Kz"]["ladder_identity_v2"] == "2"
    height6_control = data["independent_spin_dos_route"]["height6_control"]
    assert height6_control["scope"] == (
        "five-shape control family only; not a claim about every open box shape"
    )
    assert height6_control["shapes"] == [list(box) for box in control_boxes]
    rebuilt_height6 = {
        f"{width}x{height}x{layers}": {
            str(degree): not any(control_weights[(width, height, layers)][degree])
            for degree in (2, 4, 6, 8)
        }
        for width, height, layers in control_boxes
    }
    assert height6_control["per_box_w_columns_zero"] == rebuilt_height6
    assert height6_control["all_per_box_w_columns_zero"] is True
    inventory = data["bell_eighth_order_inventory"]
    assert inventory["bell_number_B8"] == 4140
    assert inventory["pattern_count"] == len(EXPECTED_PATTERNS)
    assert inventory["total_admissible_partitions"] == sum(EXPECTED_ADMISSIBLE)
    assert inventory["total_connected_monomials"] == total_connected
    assert inventory["atom_kinds_present"] == ["G", "U", "W", "W8"]
    assert inventory["placement_factors_independently_derived"] is True
    assert inventory["normalized_gap_multisets_tested"] == len(all_gap_multisets)
    assert inventory["literal_nonzero_gap_multisets"] == [
        list(values) for values in sorted(literal_nonzero)
    ]
    assert inventory["selection_rule_verified"] is True
    assert inventory["inadmissible_patterns_vanish"] is True
    assert inventory["collected_digest_sha256"] == polynomial_digest(polynomials)
    assert [tuple(row["pattern"]) for row in inventory["patterns"]] == list(EXPECTED_PATTERNS)
    assert [row["admissible_partitions"] for row in inventory["patterns"]] == list(EXPECTED_ADMISSIBLE)
    assert [row["raw_layer_moment_monomials"] for row in inventory["patterns"]] == list(EXPECTED_RAW_MONOMIALS)
    assert [row["connected_monomials"] for row in inventory["patterns"]] == list(EXPECTED_CONNECTED_MONOMIALS)
    for row, report in zip(inventory["patterns"], row_reports, strict=True):
        pattern, placement, ordered_assignments, histogram, moebius_by_class, digest = report
        assert tuple(row["pattern"]) == pattern
        assert row["placement_weight"] == str(placement)
        assert row["ordered_gap_assignments"] == ordered_assignments
        assert row["atom_kind_histogram"] == histogram
        assert row["moebius_by_block_class"] == moebius_by_class
        assert row["connected_polynomial_sha256"] == digest
        assert row["all_even_incidence"] and row["all_anchor_connected"]
    resources = payload["provenance"]["resources"]
    assert resources["wall_budget_seconds"] == 1800
    assert resources["observed_wall_is_measured"] is True
    assert resources["peak_rss_budget_is_preflight"] is True

    elapsed = time.process_time() - started
    assert elapsed <= TEST_WALL_BUDGET_SECONDS, (elapsed, TEST_WALL_BUDGET_SECONDS)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(
        "PASS test_interlayer_c8 "
        f"({elapsed:.1f}s CPU; ru_maxrss={rss}; cpu_budget={TEST_WALL_BUDGET_SECONDS}s)"
    )


if __name__ == "__main__":
    main()
