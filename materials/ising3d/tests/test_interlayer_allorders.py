"""Independent exact controls for the all-order interlayer theorem.

Run from the repository root with
    .venv/bin/python tests/test_interlayer_allorders.py

Every control here is exact.  The set-partition combinatorics is recomputed from
scratch, and the c2/c4/c6 series are rebuilt from full spin enumeration of open
3D boxes (``ising.interlayer.anisotropic_box_even_subgraph`` transforms the
integer joint density of states), a route that shares no code with the
edge-subset generator under test.
"""

from __future__ import annotations

import importlib.util
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import factorial
from pathlib import Path
from typing import Sequence

from ising.interlayer import anisotropic_box_even_subgraph

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "e75_interlayer_allorders.py"
ARTIFACT = ROOT / "results" / "interlayer" / "allorders.json"
TEST_V_ORDER = 4
TEST_W_ORDER = 6

RationalSeries = tuple[Fraction, ...]


def _load_experiment():
    spec = importlib.util.spec_from_file_location("e75_interlayer_allorders", EXPERIMENT)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load the all-orders experiment module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------
# independent exact series helpers
# --------------------------------------------------------------------------
def _zero(order: int) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def _add_scaled(target: list[Fraction], source: Sequence[Fraction], scale: int) -> None:
    for degree, coefficient in enumerate(source):
        if coefficient:
            target[degree] += scale * coefficient


def _mul(left: Sequence[Fraction], right: Sequence[Fraction], order: int) -> RationalSeries:
    out = _zero(order)
    for left_degree, left_value in enumerate(left):
        if not left_value or left_degree > order:
            continue
        for right_degree, right_value in enumerate(right):
            if left_degree + right_degree > order:
                break
            if right_value:
                out[left_degree + right_degree] += left_value * right_value
    return tuple(out)


def _divide(
    numerator: Sequence[int | Fraction], denominator: Sequence[int | Fraction], order: int
) -> RationalSeries:
    den = [Fraction(value) for value in denominator]
    while len(den) <= order:
        den.append(Fraction(0))
    quotient = _zero(order)
    for degree in range(order + 1):
        value = Fraction(numerator[degree]) if degree < len(numerator) else Fraction(0)
        for lower in range(degree):
            value -= quotient[lower] * den[degree - lower]
        quotient[degree] = value / den[0]
    return tuple(quotient)


def _log_columns(columns: dict[int, RationalSeries], order: int, w_order: int):
    """Exact ``[w^d] log P`` from the ``[w^d] P`` columns, by Newton recursion."""

    ratios = {
        degree: _divide(columns[degree], columns[0], order)
        for degree in range(1, w_order + 1)
    }
    out: dict[int, RationalSeries] = {}
    for degree in range(1, w_order + 1):
        accumulated = _zero(order)
        _add_scaled(accumulated, ratios[degree], degree)
        for lower in range(1, degree):
            _add_scaled(accumulated, _mul(out[lower], ratios[degree - lower], order), -lower)
        out[degree] = tuple(value / degree for value in accumulated)
    return out


def _box_log_columns(width: int, height: int, layers: int):
    """Spin-enumeration route: exact log columns of one open box."""

    table = anisotropic_box_even_subgraph((width, height, layers), TEST_V_ORDER, TEST_W_ORDER)
    columns = {
        w_degree: tuple(
            Fraction(table[v_degree][w_degree]) for v_degree in range(TEST_V_ORDER + 1)
        )
        for w_degree in range(TEST_W_ORDER + 1)
    }
    odd_zero = all(
        not any(columns[degree]) for degree in range(1, TEST_W_ORDER + 1, 2)
    )
    return _log_columns(columns, TEST_V_ORDER, TEST_W_ORDER), odd_zero


def _brute_force_bulk():
    """Independent 3D finite-lattice inversion of the spin-enumerated boxes."""

    span_budget = TEST_V_ORDER // 2
    layer_budget = TEST_W_ORDER // 2 + 1
    shapes = tuple(
        sorted(
            (
                (width, height, layers)
                for width in range(1, span_budget + 2)
                for height in range(1, span_budget + 2)
                if (width - 1) + (height - 1) <= span_budget
                for layers in range(1, layer_budget + 1)
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    degrees = tuple(range(2, TEST_W_ORDER + 1, 2))
    weights: dict[tuple[int, int, int], dict[int, RationalSeries]] = {}
    bulk = {degree: _zero(TEST_V_ORDER) for degree in degrees}
    by_layers = {
        layers: {degree: _zero(TEST_V_ORDER) for degree in degrees}
        for layers in range(1, layer_budget + 1)
    }
    odd_zero = True
    for width, height, layers in shapes:
        columns, box_odd_zero = _box_log_columns(width, height, layers)
        odd_zero = odd_zero and box_odd_zero
        values = {degree: list(columns[degree]) for degree in degrees}
        for (sub_width, sub_height, sub_layers), subweight in weights.items():
            if sub_width <= width and sub_height <= height and sub_layers <= layers:
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sub_layers + 1)
                )
                for degree in degrees:
                    _add_scaled(values[degree], subweight[degree], -placements)
        exact = {degree: tuple(values[degree]) for degree in degrees}
        weights[(width, height, layers)] = exact
        for degree in degrees:
            _add_scaled(bulk[degree], exact[degree], 1)
            _add_scaled(by_layers[layers][degree], exact[degree], 1)
    return (
        {degree: tuple(bulk[degree]) for degree in degrees},
        {layers: {degree: tuple(values[degree]) for degree in values} for layers, values in by_layers.items()},
        odd_zero,
    )


# --------------------------------------------------------------------------
# independent set-partition combinatorics
# --------------------------------------------------------------------------
def _partitions(items: tuple[int, ...]):
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for partition in _partitions(rest):
        yield ((first,),) + partition
        for index in range(len(partition)):
            yield partition[:index] + ((first,) + partition[index],) + partition[index + 1 :]


def _mobius(blocks: int) -> int:
    return (-1) ** (blocks - 1) * factorial(blocks - 1)


def _odd_positions(values) -> tuple[int, ...]:
    counts = Counter(values)
    return tuple(sorted(value for value, count in counts.items() if count % 2))


def _direct_cumulant(gaps, positions, moment) -> Fraction:
    total = Fraction(0)
    for partition in _partitions(tuple(range(len(gaps)))):
        term = Fraction(_mobius(len(partition)))
        for block in partition:
            per_layer = defaultdict(list)
            for bond in block:
                per_layer[gaps[bond]].append(positions[bond])
                per_layer[gaps[bond] + 1].append(positions[bond])
            for layer, layer_positions in per_layer.items():
                term *= moment(layer, _odd_positions(layer_positions))
        total += term
    return total


def _independent_patterns(order: int) -> tuple[tuple[int, ...], ...]:
    """Exhaustive normalized gap multisets with the two selection properties."""

    found: set[tuple[int, ...]] = set()
    stack: list[tuple[int, ...]] = [()]
    while stack:
        current = stack.pop()
        if len(current) == order:
            counts = Counter(current)
            interval = sorted(counts) == list(range(max(counts) + 1))
            even = all(multiplicity % 2 == 0 for multiplicity in counts.values())
            if interval and even:
                found.add(current)
            continue
        start = current[-1] if current else 0
        for gap in range(start, order // 2):
            stack.append(current + (gap,))
    return tuple(sorted(found))


def _slab_partition_odd_coefficients(
    layers: int, sites_per_layer: int, edges, weight: Fraction, max_degree: int
) -> dict[int, Fraction]:
    """Exact q-Taylor coefficients of a finite open slab partition function."""

    spin_count = layers * sites_per_layer
    coefficients = {degree: Fraction(0) for degree in range(max_degree + 1)}
    for state in range(1 << spin_count):
        spins = tuple(1 if (state >> index) & 1 else -1 for index in range(spin_count))
        horizontal = 0
        for layer in range(layers):
            offset = layer * sites_per_layer
            horizontal += sum(
                spins[offset + left] * spins[offset + right] for left, right in edges
            )
        vertical = sum(
            spins[layer * sites_per_layer + site] * spins[(layer + 1) * sites_per_layer + site]
            for layer in range(layers - 1)
            for site in range(sites_per_layer)
        )
        boltzmann = weight**horizontal
        term = Fraction(1)
        for degree in range(max_degree + 1):
            if degree:
                term *= Fraction(vertical, degree)
            coefficients[degree] += boltzmann * term
    return coefficients


def main() -> None:
    module = _load_experiment()
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e75_interlayer_allorders.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert payload["provenance"]["arithmetic"] == "Python int and Fraction only"
    assert payload["checks"] and all(row["passed"] for row in payload["checks"])
    data = payload["data"]

    # ---- 1. independent set-partition counts at m=1,2,3 ----
    bell = {size: sum(1 for _ in _partitions(tuple(range(size)))) for size in (2, 4, 6)}
    assert bell == {2: 2, 4: 15, 6: 203}
    assert data["structural_counts"]["bell_numbers"] == {"2": 2, "4": 15, "6": 203}

    # ---- 2. independent selection rule and composition count ----
    for m in (1, 2, 3, 4):
        independent = _independent_patterns(2 * m)
        assert independent == module.admissible_gap_patterns(2 * m)
        assert len(independent) == 2 ** (m - 1)
        assert max(max(pattern) + 2 for pattern in independent) == m + 1
        for pattern in independent:
            counts = Counter(pattern)
            assert sorted(counts) == list(range(max(counts) + 1))
            assert all(value % 2 == 0 and value > 0 for value in counts.values())
    assert data["theorem_layer_gap_selection"]["patterns_m1_m2_m3"]["3"] == [
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 1],
        [0, 0, 1, 1, 1, 1],
        [0, 0, 1, 1, 2, 2],
    ]

    # ---- 3. generator versus literal cumulant definition, m=1,2,3 ----
    def moment(layer: int, positions: tuple[int, ...]) -> Fraction:
        if len(positions) % 2:
            return Fraction(0)
        if not positions:
            return Fraction(1)
        code = sum((index + 5) * (position + 3) for index, position in enumerate(positions))
        return Fraction((layer + 3) * (code + 2), code + 11)

    for m in (1, 2, 3):
        positions = tuple(range(2 * m))
        for pattern in module.admissible_gap_patterns(2 * m):
            generated = module.evaluate_polynomial(
                module.cumulant_polynomial(pattern, positions), moment
            )
            assert generated == _direct_cumulant(pattern, positions, moment)

    # An inadmissible pattern must be identically zero as a formal polynomial,
    # both when parity fails and when the gap support is disconnected.
    for pattern in (
        (0, 0, 0, 1),
        (0, 1, 1, 2),
        (0, 0, 2, 2),
        (0, 2, 2, 4),
        (0, 0, 0, 0, 2, 2),
        (0, 0, 1, 1, 3, 3),
        (0, 0, 0, 1, 1, 1),
    ):
        assert module.cumulant_polynomial(pattern, tuple(range(len(pattern)))) == {}
        assert _direct_cumulant(pattern, tuple(range(len(pattern))), moment) == 0

    # ---- 4. certified c2/c4 formulas ----
    assert module.cumulant_polynomial((0, 0), (0, 1)) == {((0, (0, 1)), (1, (0, 1))): 1}
    assert module.symbolic_key(
        module.cumulant_polynomial((0, 0, 0, 0), (0, 1, 2, 3))
    ) == module.symbolic_key(module.same_gap_c4_expected())
    assert module.symbolic_key(
        module.cumulant_polynomial((0, 0, 1, 1), (0, 1, 2, 3))
    ) == module.symbolic_key(module.adjacent_c4_expected())

    # ---- 5. spin-enumeration rebuild of c2, c4 and c6 in the w variable ----
    bulk, by_layers, odd_zero = _brute_force_bulk()
    assert odd_zero
    stored = data["exact_series"]
    constants = {2: Fraction(1, 2), 4: Fraction(1, 4), 6: Fraction(1, 6)}
    keys = {2: "c2_direct_q", 4: "c4_total_w", 6: "c6_total_w"}
    for degree in (2, 4, 6):
        rebuilt = tuple(
            bulk[degree][index] + (constants[degree] if index == 0 else 0)
            for index in range(TEST_V_ORDER + 1)
        )
        recorded = tuple(Fraction(value) for value in stored[keys[degree]][: TEST_V_ORDER + 1])
        assert rebuilt == recorded, (degree, rebuilt, recorded)

    # ---- 6. vertical extent bound at w^2, w^4, w^6 ----
    for m in (1, 2, 3):
        for layers, columns in by_layers.items():
            if layers > m + 1:
                assert not any(columns[2 * m])
        assert any(by_layers[m + 1][2 * m])

    # ---- 7. maximal-extent c6 pattern equals the four-layer weight ----
    maximal = tuple(
        Fraction(value)
        for value in stored["c6_pattern_components_direct_q"]["0,0,1,1,2,2"][: TEST_V_ORDER + 1]
    )
    assert by_layers[4][6] == maximal

    # ---- 8. reflection pairing of the two three-layer c6 patterns ----
    components = stored["c6_pattern_components_direct_q"]
    assert components["0,0,0,0,1,1"] == components["0,0,1,1,1,1"]

    # ---- 9. ladder identity at second in-plane order ----
    for m in (1, 2, 3):
        assert bulk[2 * m][2] == 2
    assert data["structural_counts"]["ladder_identity_v2_coefficients"] == {
        str(m): "2" for m in range(1, 7)
    }

    # ---- 10. decoupled normalization ----
    assert Fraction(stored["c6_direct_q"][0]) == Fraction(1, 45)
    assert Fraction(stored["c4_direct_q"][0]) == Fraction(-1, 12)
    assert Fraction(stored["c2_direct_q"][0]) == Fraction(1, 2)

    # ---- 11. finite slab evenness, and its failure without layer parity ----
    slab = _slab_partition_odd_coefficients(
        layers=3,
        sites_per_layer=4,
        edges=((0, 1), (1, 2), (2, 3), (3, 0)),
        weight=Fraction(5, 4),
        max_degree=7,
    )
    assert all(slab[degree] == 0 for degree in (1, 3, 5, 7))
    assert slab[0] > 0 and slab[2] > 0

    triangle_ring = {
        degree: Fraction(0) for degree in range(4)
    }
    for state in range(1 << 3):
        spins = tuple(1 if (state >> index) & 1 else -1 for index in range(3))
        energy = spins[0] * spins[1] + spins[1] * spins[2] + spins[2] * spins[0]
        term = Fraction(1)
        for degree in range(4):
            if degree:
                term *= Fraction(energy, degree)
            triangle_ring[degree] += term
    assert triangle_ring[1] == 0 and triangle_ring[3] == 8
    assert data["structural_counts"]["odd_periodic_ring_q3"] == "8"

    # ---- 12. honest limitations ----
    limitations = data["limitations"]
    assert limitations["tag"] == "[UNRESOLVED]"
    assert "No sign or positivity theorem" in limitations["sign"]
    assert "No radius of convergence" in limitations["radius"]
    print("PASS test_interlayer_allorders")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL {exc}")
        raise
