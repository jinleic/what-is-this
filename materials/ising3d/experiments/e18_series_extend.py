"""Extend the simple-cubic HT/LT series with the modular CRT transfer engine.

Run from the repository root:
    .venv/bin/python experiments/e18_series_extend.py
"""

from __future__ import annotations

import hashlib
import json
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path

import ising.series as flm
from ising.series import FLMSeries, high_temperature_free_energy, low_temperature_free_energy
from ising.transfer_matrix import (
    box_broken_bond_poly as int64_box_broken_bond_poly,
)
from ising.transfer_matrix import (
    torus_broken_bond_poly as int64_torus_broken_bond_poly,
)
from ising.transfer_matrix.crt import (
    box_broken_bond_poly as crt_box_broken_bond_poly,
)
from ising.transfer_matrix.crt import (
    estimate_box_peak_bytes,
    primes_for_sites,
    torus_broken_bond_poly as crt_torus_broken_bond_poly,
)

SCRIPT = "experiments/e18_series_extend.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "series"
HT_ORDER = 20
LT_ORDER = 28
DIRECT_ORDER = 10


def _fraction_strings(values):
    return [str(value) for value in values]


def _nonzero_coefficients(values):
    return {str(index): str(value) for index, value in enumerate(values) if value}


def _canonical_boxes(boxes):
    return sorted(
        {tuple(sorted(shape)) for shape in boxes}, key=lambda shape: (sum(shape), shape)
    )


def _series_payload(series: FLMSeries) -> dict:
    canonical = _canonical_boxes(series.boxes)
    return {
        "dimension": series.dimension,
        "variable": series.variable,
        "achieved_order": series.order,
        "coefficients": _fraction_strings(series.coefficients),
        "nonzero_coefficients": _nonzero_coefficients(series.coefficients),
        "interaction_coefficients": _fraction_strings(
            series.interaction_coefficients
        ),
        "nonzero_interaction_coefficients": _nonzero_coefficients(
            series.interaction_coefficients
        ),
        "ordered_box_count": len(series.boxes),
        "canonical_box_count": len(canonical),
        "canonical_box_list": [list(shape) for shape in canonical],
        "order_bound": {
            "statement": series.order_bound,
            "budget": series.bound_budget,
            "bound_slack": series.bound_slack,
        },
    }


def _record_check(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{name}: {detail}")
    print(f"PASS {name}: {detail}", flush=True)


def _clear_flm_transfer_caches() -> None:
    flm._free_broken_polynomial.cache_clear()
    flm._plus_broken_polynomial.cache_clear()
    flm._ht_box_log.cache_clear()
    flm._lt_box_log.cache_clear()


@contextmanager
def _hybrid_flm_engine():
    """Temporarily route only N>62 boxes through CRT, preserving existing FLM."""

    original = flm.box_broken_bond_poly
    records: dict[tuple[tuple[int, ...], bool], dict] = {}

    def hybrid(shape, periodic=None, plus_boundary=False):
        canonical = tuple(int(side) for side in shape)
        sites = prod(canonical)
        started = time.perf_counter()
        if sites <= 62:
            coefficients = int64_box_broken_bond_poly(
                canonical,
                periodic=periodic,
                plus_boundary=plus_boundary,
            )
            engine = "int64"
            moduli = None
        else:
            coefficients = crt_box_broken_bond_poly(
                canonical,
                periodic=periodic,
                plus_boundary=plus_boundary,
            )
            engine = "31-bit-prime CRT"
            moduli = primes_for_sites(sites)
        elapsed = time.perf_counter() - started
        record = {
            "shape": list(canonical),
            "sites": sites,
            "cross_section_sites": prod(canonical[:-1]),
            "engine": engine,
            "plus_boundary": bool(plus_boundary),
            "degree": len(coefficients) - 1,
            "nonzero_coefficient_count": sum(value != 0 for value in coefficients),
            "coefficient_sum": str(sum(coefficients)),
            "elapsed_seconds": f"{elapsed:.6f}",
        }
        if moduli is not None:
            record.update(
                {
                    "primes": list(moduli.primes),
                    "prime_product": str(moduli.product),
                    "prime_product_bits": moduli.product_bits,
                    "coefficient_sha256": hashlib.sha256(
                        json.dumps(coefficients, separators=(",", ":")).encode("ascii")
                    ).hexdigest(),
                    "estimated_peak_bytes": estimate_box_peak_bytes(
                        canonical,
                        periodic=periodic,
                        plus_boundary=plus_boundary,
                    ),
                }
            )
        records[(canonical, bool(plus_boundary))] = record
        return coefficients

    _clear_flm_transfer_caches()
    flm.box_broken_bond_poly = hybrid
    try:
        yield records
    finally:
        flm.box_broken_bond_poly = original
        _clear_flm_transfer_caches()


def _load_coefficients(filename: str) -> tuple[Fraction, ...]:
    payload = json.loads((RESULT_DIR / filename).read_text(encoding="utf-8"))
    return tuple(Fraction(value) for value in payload["data"]["coefficients"])


def _crt_int64_crosschecks() -> list[dict]:
    cases = [
        ("box", (2, 2), {}),
        ("box", (3, 3), {}),
        ("box", (2, 2, 2), {}),
        ("box", (2, 3, 2), {}),
        ("box", (3, 3, 2), {}),
        ("box", (2, 2, 3), {"plus_boundary": True}),
        ("box", (2, 3, 1), {"plus_boundary": True}),
        ("box", (3, 2, 2), {"periodic": (True, False)}),
        ("box", (2, 31), {}),
        ("torus", (4, 4), {}),
        ("torus", (2, 2, 2), {}),
        ("torus", (2, 3, 2), {"block": 2}),
    ]
    records = []
    for kind, shape, kwargs in cases:
        started = time.perf_counter()
        if kind == "box":
            expected = int64_box_broken_bond_poly(shape, **kwargs)
            actual = crt_box_broken_bond_poly(shape, **kwargs)
        else:
            expected = int64_torus_broken_bond_poly(shape, **kwargs)
            actual = crt_torus_broken_bond_poly(shape, **kwargs)
        if actual != expected:
            raise AssertionError(f"CRT/int64 mismatch for {kind} {shape} {kwargs}")
        records.append(
            {
                "kind": kind,
                "shape": list(shape),
                "options": kwargs,
                "sites": prod(shape),
                "degree": len(actual) - 1,
                "coefficient_sum": str(sum(actual)),
                "primes": list(primes_for_sites(prod(shape)).primes),
                "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
            }
        )
    return records


def _box_edges(shape: tuple[int, ...]) -> tuple[int, tuple[tuple[int, int], ...]]:
    strides = []
    stride = 1
    for side in reversed(shape):
        strides.append(stride)
        stride *= side
    strides.reverse()

    def index(coordinate):
        return sum(value * step for value, step in zip(coordinate, strides))

    edges = []
    for coordinate in product(*(range(side) for side in shape)):
        left = index(coordinate)
        for axis, side in enumerate(shape):
            if coordinate[axis] + 1 >= side:
                continue
            neighbour = list(coordinate)
            neighbour[axis] += 1
            edges.append((left, index(neighbour)))
    return prod(shape), tuple(edges)


def _fundamental_cycle_basis(
    vertex_count: int, edges: tuple[tuple[int, int], ...]
) -> tuple[int, ...]:
    adjacency: list[list[tuple[int, int]]] = [[] for _ in range(vertex_count)]
    for edge_index, (left, right) in enumerate(edges):
        adjacency[left].append((right, edge_index))
        adjacency[right].append((left, edge_index))

    parent = [-1] * vertex_count
    parent_edge = [-1] * vertex_count
    depth = [0] * vertex_count
    tree_edges: set[int] = set()
    stack = [0]
    parent[0] = 0
    while stack:
        vertex = stack.pop()
        for neighbour, edge_index in adjacency[vertex]:
            if parent[neighbour] != -1:
                continue
            parent[neighbour] = vertex
            parent_edge[neighbour] = edge_index
            depth[neighbour] = depth[vertex] + 1
            tree_edges.add(edge_index)
            stack.append(neighbour)
    if any(value == -1 for value in parent):
        raise AssertionError("free box graph must be connected")

    basis = []
    for chord, (initial_left, initial_right) in enumerate(edges):
        if chord in tree_edges:
            continue
        left, right = initial_left, initial_right
        cycle = 1 << chord
        while depth[left] > depth[right]:
            cycle ^= 1 << parent_edge[left]
            left = parent[left]
        while depth[right] > depth[left]:
            cycle ^= 1 << parent_edge[right]
            right = parent[right]
        while left != right:
            cycle ^= 1 << parent_edge[left]
            cycle ^= 1 << parent_edge[right]
            left = parent[left]
            right = parent[right]
        basis.append(cycle)
    expected_rank = len(edges) - vertex_count + 1
    if len(basis) != expected_rank:
        raise AssertionError("fundamental-cycle rank mismatch")
    return tuple(basis)


def _vertex_coordinates(vertex: int, shape: tuple[int, ...]) -> tuple[int, ...]:
    coordinates = []
    for side in reversed(shape):
        coordinates.append(vertex % side)
        vertex //= side
    return tuple(reversed(coordinates))


def _connected_vertices(
    edge_mask: int, edges: tuple[tuple[int, int], ...]
) -> frozenset[int] | None:
    adjacency: dict[int, list[int]] = {}
    remaining = edge_mask
    while remaining:
        lowest_bit = remaining & -remaining
        edge_index = lowest_bit.bit_length() - 1
        remaining ^= lowest_bit
        left, right = edges[edge_index]
        adjacency.setdefault(left, []).append(right)
        adjacency.setdefault(right, []).append(left)
    if not adjacency:
        return None
    start = next(iter(adjacency))
    visited = {start}
    stack = [start]
    while stack:
        vertex = stack.pop()
        for neighbour in adjacency[vertex]:
            if neighbour not in visited:
                visited.add(neighbour)
                stack.append(neighbour)
    if len(visited) != len(adjacency):
        return None
    return frozenset(visited)


@lru_cache(maxsize=None)
def _connected_even_polymers(
    shape: tuple[int, ...], order: int
) -> tuple[tuple[tuple[int, frozenset[tuple[int, ...]]], ...], int]:
    """List connected even edge sets whose exact bounding box is ``shape``."""

    vertices, edges = _box_edges(shape)
    basis = _fundamental_cycle_basis(vertices, edges)
    polymers: list[tuple[int, frozenset[tuple[int, ...]]]] = []
    edge_mask = 0
    previous_gray = 0
    for counter in range(1, 1 << len(basis)):
        gray = counter ^ (counter >> 1)
        changed = gray ^ previous_gray
        edge_mask ^= basis[changed.bit_length() - 1]
        previous_gray = gray
        edge_count = edge_mask.bit_count()
        if edge_count > order:
            continue
        used_vertices = _connected_vertices(edge_mask, edges)
        if used_vertices is None:
            continue
        coordinates = tuple(
            _vertex_coordinates(vertex, shape) for vertex in used_vertices
        )
        minima = tuple(
            min(coordinate[axis] for coordinate in coordinates)
            for axis in range(len(shape))
        )
        maxima = tuple(
            max(coordinate[axis] for coordinate in coordinates)
            for axis in range(len(shape))
        )
        if minima != (0,) * len(shape):
            continue
        if tuple(value + 1 for value in maxima) != shape:
            continue
        polymers.append((edge_count, frozenset(coordinates)))
    return tuple(polymers), len(basis)


def _direct_connected_cluster_series(order: int = DIRECT_ORDER) -> dict:
    """Direct polymer/overlap expansion, independent of finite-lattice inversion."""

    if order >= 12:
        raise ValueError("the direct two-polymer formula is valid only below order 12")
    span_budget = order // 2
    shapes = [
        tuple(sides)
        for sides in product(range(1, span_budget + 2), repeat=3)
        if sum(side - 1 for side in sides) <= span_budget
    ]
    polymers_by_order: dict[
        int, list[frozenset[tuple[int, ...]]]
    ] = {degree: [] for degree in range(order + 1)}
    rank_by_shape: dict[tuple[int, ...], int] = {}
    for shape in shapes:
        polymers, cycle_rank = _connected_even_polymers(shape, order)
        rank_by_shape[shape] = cycle_rank
        for edge_count, vertices in polymers:
            polymers_by_order[edge_count].append(vertices)

    interaction = [Fraction(0) for _ in range(order + 1)]
    overlap_counts = [0] * (order + 1)
    for degree in range(1, order + 1):
        interaction[degree] = Fraction(len(polymers_by_order[degree]))
        ordered_overlaps = 0
        for left_order in range(1, degree):
            right_order = degree - left_order
            for left_vertices in polymers_by_order[left_order]:
                for right_vertices in polymers_by_order[right_order]:
                    translations = {
                        tuple(
                            left_coordinate[axis] - right_coordinate[axis]
                            for axis in range(3)
                        )
                        for left_coordinate in left_vertices
                        for right_coordinate in right_vertices
                    }
                    ordered_overlaps += len(translations)
        overlap_counts[degree] = ordered_overlaps
        interaction[degree] -= Fraction(ordered_overlaps, 2)

    coefficients = list(interaction)
    for degree in range(2, order + 1, 2):
        coefficients[degree] += Fraction(3, degree)
    maximum_rank = max(rank_by_shape.values(), default=0)
    return {
        "method": (
            "Enumerate every connected even edge set (polymer) in its exact normalized "
            "bounding box using an independent GF(2) cycle-space Gray code.  Count "
            "translations of every ordered overlapping polymer pair directly.  Below "
            "order 12 at most two polymers occur, so the pressure coefficient is exactly "
            "single-polymer embeddings minus one half the ordered overlap count.  This "
            "uses neither spin transfer, broken-bond conversion, finite-lattice "
            "Moebius inversion, nor CRT residues."
        ),
        "cluster_formula": (
            "b_n = number of normalized connected even n-edge polymers "
            "- (1/2)*number of ordered overlapping polymer pairs of total size n"
        ),
        "achieved_order": order,
        "coefficients": _fraction_strings(coefficients),
        "interaction_coefficients": _fraction_strings(interaction),
        "ordered_box_count": len(shapes),
        "maximum_cycle_rank": maximum_rank,
        "largest_enumeration_size": 1 << maximum_rank,
        "connected_even_subgraph_counts": {
            str(degree): len(polymers)
            for degree, polymers in polymers_by_order.items()
            if polymers
        },
        "ordered_overlapping_pair_counts": {
            str(degree): count
            for degree, count in enumerate(overlap_counts)
            if count
        },
    }


def _cross_section_twenty_probe() -> dict:
    shape = (4, 5, 2)
    started = time.perf_counter()
    crt_coefficients = crt_box_broken_bond_poly(shape)
    crt_seconds = time.perf_counter() - started
    started = time.perf_counter()
    int64_coefficients = int64_box_broken_bond_poly(shape)
    int64_seconds = time.perf_counter() - started
    if crt_coefficients != int64_coefficients:
        raise AssertionError("4x5 cross-section CRT/int64 mismatch")
    return {
        "shape": list(shape),
        "sites": prod(shape),
        "cross_section": [4, 5],
        "cross_section_sites": 20,
        "degree": len(crt_coefficients) - 1,
        "coefficient_sum": str(sum(crt_coefficients)),
        "primes": list(primes_for_sites(prod(shape)).primes),
        "estimated_peak_bytes": estimate_box_peak_bytes(shape),
        "crt_elapsed_seconds": f"{crt_seconds:.6f}",
        "int64_elapsed_seconds": f"{int64_seconds:.6f}",
        "matched_int64": True,
    }


def main() -> None:
    timestamp = datetime.now(timezone.utc).isoformat()
    provenance = {
        "script": SCRIPT,
        "timestamp": timestamp,
        "precision": (
            "exact Python int/Fraction; numpy.int64 residues modulo deterministic "
            "31-bit primes; CRT product strictly exceeds 2^N"
        ),
    }

    crt_crosschecks = _crt_int64_crosschecks()
    print(f"PASS CRT/int64 on {len(crt_crosschecks)} lattices", flush=True)

    with _hybrid_flm_engine() as transfer_records:
        sc_ht = high_temperature_free_energy(3, HT_ORDER)
        sc_ht_enlarged = high_temperature_free_energy(
            3, HT_ORDER, bound_slack=1
        )
        sc_lt = low_temperature_free_energy(3, LT_ORDER)
        sc_lt_enlarged = low_temperature_free_energy(
            3, LT_ORDER, bound_slack=1
        )

    direct = _direct_connected_cluster_series()
    cross_section_probe = _cross_section_twenty_probe()

    old_ht = _load_coefficients("sc_ht_free_energy.json")
    old_lt = _load_coefficients("sc_lt_free_energy.json")
    ht_extra = set(sc_ht_enlarged.boxes) - set(sc_ht.boxes)
    lt_extra = set(sc_lt_enlarged.boxes) - set(sc_lt.boxes)
    crt_boxes = sorted(
        (
            record
            for record in transfer_records.values()
            if record["engine"] == "31-bit-prime CRT"
        ),
        key=lambda record: (sum(record["shape"]), record["shape"]),
    )
    free_4_cube = transfer_records.get(((4, 4, 4), False))
    largest_sites = max((record["sites"] for record in crt_boxes), default=0)
    largest_cross_section = max(
        (record["cross_section_sites"] for record in crt_boxes), default=0
    )

    ht_checks: list[dict] = []
    lt_checks: list[dict] = []
    _record_check(
        ht_checks,
        "crt_matches_int64_at_least_8_lattices",
        len(crt_crosschecks) >= 8,
        f"modular propagation exactly matches int64 on {len(crt_crosschecks)} lattices",
    )
    _record_check(
        ht_checks,
        "crt_4x4x4_free_box",
        free_4_cube is not None
        and free_4_cube["coefficient_sum"] == str(1 << 64),
        (
            "computed the full 4x4x4 free-box polynomial with three 31-bit primes; "
            f"degree={free_4_cube['degree'] if free_4_cube else 'missing'} and sum=2^64"
        ),
    )
    _record_check(
        ht_checks,
        "crt_cross_section_4x5",
        cross_section_probe["matched_int64"]
        and cross_section_probe["cross_section_sites"] == 20,
        "propagated a 4x5=20-site cross-section for the 4x5x2 box and matched int64 exactly",
    )
    _record_check(
        ht_checks,
        "sc_ht_existing_prefix_exact",
        sc_ht.coefficients[: len(old_ht)] == old_ht,
        f"every existing coefficient through v^{len(old_ht) - 1} is reproduced exactly",
    )
    _record_check(
        ht_checks,
        "sc_ht_enlarged_box_self_consistency",
        bool(ht_extra)
        and sc_ht.coefficients == sc_ht_enlarged.coefficients
        and all(not any(sc_ht_enlarged.box_weights[shape]) for shape in ht_extra),
        (
            f"adding {len(ht_extra)} ordered boxes at span budget {sc_ht.bound_budget + 1} "
            f"changes no coefficient through v^{HT_ORDER}; every added W(A) truncates to zero"
        ),
    )
    direct_coefficients = tuple(Fraction(value) for value in direct["coefficients"])
    _record_check(
        ht_checks,
        "direct_even_subgraph_crosscheck",
        direct_coefficients == sc_ht.coefficients[: DIRECT_ORDER + 1],
        (
            f"direct GF(2) even-subgraph enumeration and connected embedding weights agree "
            f"through v^{DIRECT_ORDER} (maximum cycle-space enumeration "
            f"2^{direct['maximum_cycle_rank']})"
        ),
    )
    _record_check(
        lt_checks,
        "sc_lt_existing_prefix_exact",
        sc_lt.coefficients[: len(old_lt)] == old_lt,
        f"every existing coefficient through x^{len(old_lt) - 1} is reproduced exactly",
    )
    _record_check(
        lt_checks,
        "sc_lt_enlarged_box_self_consistency",
        bool(lt_extra)
        and sc_lt.coefficients == sc_lt_enlarged.coefficients
        and all(not any(sc_lt_enlarged.box_weights[shape]) for shape in lt_extra),
        (
            f"adding {len(lt_extra)} ordered boxes at side-sum budget {sc_lt.bound_budget + 1} "
            f"changes no coefficient through x^{LT_ORDER}; every added W(A) truncates to zero"
        ),
    )
    _record_check(
        lt_checks,
        "crt_matches_int64_at_least_8_lattices",
        len(crt_crosschecks) >= 8,
        f"shared CRT correctness gate passed on {len(crt_crosschecks)} lattices",
    )

    ht_data = _series_payload(sc_ht)
    ht_data.update(
        {
            "lattice": "simple cubic",
            "normalization": "phi = log(2) + sum_n coefficients[n] * v^n",
            "interaction_normalization": (
                "phi = log(2) + 3*log(cosh K) + "
                "sum_n interaction_coefficients[n] * v^n"
            ),
            "order_bound_proof": (
                "A connected even graph spanning side lengths a_i has an Euler circuit "
                "that traverses every coordinate span out and back, so its edge count is "
                "at least 2*sum_i(a_i-1)."
            ),
            "crt_required_box_count": len(crt_boxes),
            "crt_box_records": crt_boxes,
            "largest_total_sites_reached": largest_sites,
            "largest_required_cross_section_sites": largest_cross_section,
            "cross_section_20_probe": cross_section_probe,
            "crt_int64_crosschecks": crt_crosschecks,
            "independent_direct_crosscheck": direct,
            "self_consistency": {
                "extra_ordered_boxes": len(ht_extra),
                "extra_budget": sc_ht.bound_budget + 1,
                "all_extra_weights_zero": True,
            },
        }
    )
    lt_data = _series_payload(sc_lt)
    lt_data.update(
        {
            "lattice": "simple cubic",
            "normalization": "phi = 3*K + sum_n coefficients[n] * x^n",
            "order_bound_proof": (
                "A connected droplet cluster spanning a,b,c has surface at least "
                "4*(a+b+c)-6; boxes above the retained side-sum budget cannot contribute."
            ),
            "self_consistency": {
                "extra_ordered_boxes": len(lt_extra),
                "extra_budget": sc_lt.bound_budget + 1,
                "all_extra_weights_zero": True,
            },
        }
    )

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        RESULT_DIR / "extended_sc_ht_free_energy.json": {
            "provenance": provenance,
            "data": ht_data,
            "checks": ht_checks,
        },
        RESULT_DIR / "extended_sc_lt_free_energy.json": {
            "provenance": provenance,
            "data": lt_data,
            "checks": lt_checks,
        },
    }
    for path, payload in outputs.items():
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"WROTE {path.relative_to(ROOT)}", flush=True)

    print(
        f"PASS e18_series_extend: HT v^{HT_ORDER}; LT x^{LT_ORDER}; "
        f"largest required box N={largest_sites}, required cross-section={largest_cross_section}; "
        "4x5 probe cross-section=20",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e18_series_extend: {error}", flush=True)
        raise
