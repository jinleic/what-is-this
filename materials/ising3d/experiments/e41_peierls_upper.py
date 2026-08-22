"""Exact small-area 3D Peierls census and a rigorous rational tail certificate.

The finite combinatorics use only Python integers.  The sole transcendental
quantity, exp(-2K), is enclosed by an alternating Taylor series with exact
``Fraction`` arithmetic and then rounded outward to a decimal rational interval.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import defaultdict, deque
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path


SCRIPT = "experiments/e41_peierls_upper.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "peierls_upper.json"
FRIEDLI_VELENIK = ROOT / "sources" / "fulltext" / "friedli_velenik_ising_ch3.pdf"
FRIEDLI_VELENIK_SHA256 = "4809d134deb0a97bb37bcc306345d3ace659f7543bbcd8cfedde53c5dbe47045"
BONATI = ROOT / "sources" / "fulltext" / "bonati2014_peierls.pdf"
BONATI_SHA256 = "2bd9daa39aad312c9a1141b404a9de19020c85e97d0b779644ddfe71567acd36"

A_MAX = 24
N_MAX = 8
DIRECT_A_MAX = 10
TAYLOR_ORDER = 60
EXP_DECIMAL_PLACES = 50
K0 = Fraction(848449, 500000)  # 1.696898
K_PREVIOUS_GRID = Fraction(1696897, 1000000)  # 1.696897
INCUMBENT_UPPER = Fraction(25273100985867, 100000000000000)
LAMBDA = Fraction(11**11, 10**10)
TAIL_COEFFICIENT = Fraction(78, 625)
DIRECTIONS = (
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
)
EXPECTED_LEVEL_COUNTS = {1: 1, 2: 3, 3: 15, 4: 86, 5: 534, 6: 3481, 7: 23502, 8: 162913}
EXPECTED_HEAD = {6: 1, 8: 0, 10: 6, 12: 0, 14: 45, 16: 12, 18: 332, 20: 240, 22: 2538, 24: 3040}

Cell = tuple[int, int, int]
Polycube = tuple[Cell, ...]
Plaquette = tuple[int, int, int, int]
DualEdge = tuple[int, int, int, int]


def canonical_polycube(cells: set[Cell]) -> Polycube:
    """Quotient translations only; rotations and reflections remain distinct."""
    minima = tuple(min(cell[axis] for cell in cells) for axis in range(3))
    return tuple(
        sorted(
            (cell[0] - minima[0], cell[1] - minima[1], cell[2] - minima[2])
            for cell in cells
        )
    )


def boundary_area(shape: Polycube) -> int:
    cells = set(shape)
    return sum(
        (cell[0] + step[0], cell[1] + step[1], cell[2] + step[2]) not in cells
        for cell in shape
        for step in DIRECTIONS
    )


def complement_connected(shape: Polycube) -> bool:
    """Exact finite flood fill in a one-cell-expanded bounding box."""
    cells = set(shape)
    minima = [min(cell[axis] for cell in shape) - 1 for axis in range(3)]
    maxima = [max(cell[axis] for cell in shape) + 1 for axis in range(3)]
    universe = set(
        itertools.product(*(range(minima[axis], maxima[axis] + 1) for axis in range(3)))
    )
    start = tuple(minima)
    exterior = {start}
    queue: deque[Cell] = deque([start])
    while queue:
        cell = queue.popleft()
        for step in DIRECTIONS:
            neighbor = (
                cell[0] + step[0],
                cell[1] + step[1],
                cell[2] + step[2],
            )
            if neighbor in universe and neighbor not in cells and neighbor not in exterior:
                exterior.add(neighbor)
                queue.append(neighbor)
    return cells | exterior == universe


def shape_digest(shapes: list[Polycube]) -> str:
    encoded = "\n".join(repr(shape) for shape in sorted(shapes)).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def enumerate_polycube_head() -> dict[str, object]:
    """Enumerate fixed polycubes through eight cells by canonical BFS growth."""
    current: set[Polycube] = {((0, 0, 0),)}
    level_counts: dict[int, int] = {}
    rooted_counts: defaultdict[int, int] = defaultdict(int)
    unrooted_counts: defaultdict[int, int] = defaultdict(int)
    by_cells_area: defaultdict[int, defaultdict[int, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    relevant_shapes: defaultdict[int, list[Polycube]] = defaultdict(list)
    relevant_with_cavities = 0

    for size in range(1, N_MAX + 1):
        level_counts[size] = len(current)
        for shape in current:
            area = boundary_area(shape)
            by_cells_area[size][area] += 1
            if area <= A_MAX:
                if complement_connected(shape):
                    unrooted_counts[area] += 1
                    rooted_counts[area] += size
                    relevant_shapes[area].append(shape)
                else:
                    relevant_with_cavities += 1
        if size == N_MAX:
            break
        following: set[Polycube] = set()
        for shape in current:
            cells = set(shape)
            for cell in shape:
                for step in DIRECTIONS:
                    neighbor = (
                        cell[0] + step[0],
                        cell[1] + step[1],
                        cell[2] + step[2],
                    )
                    if neighbor not in cells:
                        following.add(canonical_polycube(cells | {neighbor}))
        current = following

    return {
        "level_counts": level_counts,
        "rooted_counts": {area: rooted_counts.get(area, 0) for area in range(6, A_MAX + 1, 2)},
        "unrooted_counts": {
            area: unrooted_counts.get(area, 0) for area in range(6, A_MAX + 1, 2)
        },
        "by_cells_area": {
            size: dict(sorted(counts.items())) for size, counts in sorted(by_cells_area.items())
        },
        "shape_sha256_by_area": {
            area: shape_digest(shapes) for area, shapes in sorted(relevant_shapes.items())
        },
        "relevant_shape_count": sum(len(shapes) for shapes in relevant_shapes.values()),
        "relevant_with_cavities": relevant_with_cavities,
    }


def plaquette_edges(plaquette: Plaquette) -> tuple[DualEdge, ...]:
    normal, *center = plaquette
    tangent = [axis for axis in range(3) if axis != normal]
    first, second = tangent
    edges: list[DualEdge] = []
    for sign in (-1, 1):
        edge_center = center.copy()
        edge_center[second] += sign
        edges.append((first, *edge_center))
        edge_center = center.copy()
        edge_center[first] += sign
        edges.append((second, *edge_center))
    return tuple(edges)


def incident_plaquettes(edge: DualEdge) -> tuple[Plaquette, ...]:
    direction, *center = edge
    transverse = [axis for axis in range(3) if axis != direction]
    plaquettes: list[Plaquette] = []
    for normal, tangent in ((transverse[0], transverse[1]), (transverse[1], transverse[0])):
        for sign in (-1, 1):
            plaquette_center = center.copy()
            plaquette_center[tangent] += sign
            plaquettes.append((normal, *plaquette_center))
    return tuple(plaquettes)


def positive_face(cell: Cell, axis: int) -> Plaquette:
    center = [2 * coordinate for coordinate in cell]
    center[axis] += 1
    return (axis, *center)


def cell_boundary(cells: set[Cell]) -> frozenset[Plaquette]:
    result: set[Plaquette] = set()
    for cell in cells:
        for axis in range(3):
            lower = tuple(cell[index] - (index == axis) for index in range(3))
            upper = tuple(cell[index] + (index == axis) for index in range(3))
            if lower not in cells:
                result.add(positive_face(lower, axis))
            if upper not in cells:
                result.add(positive_face(cell, axis))
    return frozenset(result)


def plaquette_interior(surface: frozenset[Plaquette]) -> set[Cell]:
    adjacent_cells: list[Cell] = []
    for plaquette in surface:
        normal, *center = plaquette
        lower = [coordinate // 2 for coordinate in center]
        lower[normal] = (center[normal] - 1) // 2
        upper = lower.copy()
        upper[normal] += 1
        adjacent_cells.extend((tuple(lower), tuple(upper)))
    minima = [min(cell[axis] for cell in adjacent_cells) - 1 for axis in range(3)]
    maxima = [max(cell[axis] for cell in adjacent_cells) + 1 for axis in range(3)]
    universe = set(
        itertools.product(*(range(minima[axis], maxima[axis] + 1) for axis in range(3)))
    )
    start = tuple(minima)
    exterior = {start}
    queue: deque[Cell] = deque([start])
    while queue:
        cell = queue.popleft()
        for axis in range(3):
            for sign in (-1, 1):
                neighbor_list = list(cell)
                neighbor_list[axis] += sign
                neighbor = tuple(neighbor_list)
                if neighbor not in universe or neighbor in exterior:
                    continue
                lower = neighbor if sign < 0 else cell
                if positive_face(lower, axis) not in surface:
                    exterior.add(neighbor)
                    queue.append(neighbor)
    return universe - exterior


def direct_plaquette_counts(max_area: int = DIRECT_A_MAX) -> tuple[dict[int, int], dict[str, int]]:
    """Independent closure search: branch on the first odd-incidence dual edge."""
    roots = [(0, 2 * offset + 1, 0, 0) for offset in range(max_area)]
    stack: list[tuple[frozenset[Plaquette], frozenset[DualEdge]]] = [
        (frozenset({root}), frozenset(plaquette_edges(root))) for root in roots
    ]
    seen: set[frozenset[Plaquette]] = set()
    completed: set[frozenset[Plaquette]] = set()
    popped = 0

    while stack:
        selected, odd_edges = stack.pop()
        popped += 1
        if selected in seen:
            continue
        seen.add(selected)
        if not odd_edges:
            completed.add(selected)
            # A closed proper subset of a connected target can be extended
            # across one of its shared-edge adjacencies.  Adding just that
            # plaquette starts a new open frontier with its four odd edges.
            # The one-plaquette state itself is already present from the ray
            # roots and does not need this extension.
            if 1 < len(selected) < max_area:
                adjacent: set[Plaquette] = set()
                for plaquette in selected:
                    for edge in plaquette_edges(plaquette):
                        adjacent.update(incident_plaquettes(edge))
                for plaquette in adjacent - set(selected):
                    stack.append(
                        (
                            selected | {plaquette},
                            frozenset(plaquette_edges(plaquette)),
                        )
                    )
            continue
        if len(selected) >= max_area:
            continue
        edge = min(odd_edges)
        for plaquette in incident_plaquettes(edge):
            if plaquette in selected:
                continue
            next_odd = set(odd_edges)
            for boundary_edge in plaquette_edges(plaquette):
                if boundary_edge in next_odd:
                    next_odd.remove(boundary_edge)
                else:
                    next_odd.add(boundary_edge)
            stack.append((selected | {plaquette}, frozenset(next_odd)))

    enclosing: defaultdict[int, set[frozenset[Plaquette]]] = defaultdict(set)
    origin = (0, 0, 0)
    for surface in completed:
        interior = plaquette_interior(surface)
        if origin in interior and cell_boundary(interior) == surface:
            enclosing[len(surface)].add(surface)
    counts = {area: len(enclosing.get(area, set())) for area in range(6, max_area + 1, 2)}
    stats = {
        "stack_states_popped": popped,
        "distinct_partial_sets": len(seen),
        "closed_sets_before_enclosure_filter": len(completed),
        "enclosing_closed_sets": sum(counts.values()),
    }
    return counts, stats


def exp_negative_taylor_bounds(x: Fraction, order: int = TAYLOR_ORDER) -> tuple[Fraction, Fraction]:
    """Exact alternating-series enclosure of exp(-x), with even upper partial sum."""
    if x < 0 or order % 2 or Fraction(order + 2) <= x:
        raise ValueError("need x>=0, even order, and a decreasing Taylor tail")
    term = Fraction(1)
    partial = term
    upper: Fraction | None = None
    for index in range(1, order + 2):
        term *= -x / index
        partial += term
        if index == order:
            upper = partial
    assert upper is not None and 0 < partial < upper
    return partial, upper


def outward_decimal_interval(
    lower: Fraction, upper: Fraction, places: int = EXP_DECIMAL_PLACES
) -> tuple[Fraction, Fraction]:
    scale = 10**places
    lower_integer = lower.numerator * scale // lower.denominator
    upper_integer = (upper.numerator * scale + upper.denominator - 1) // upper.denominator
    return Fraction(lower_integer, scale), Fraction(upper_integer, scale)


def criterion_terms(y: Fraction, head: dict[int, int]) -> tuple[Fraction, Fraction, Fraction]:
    exact_head = sum(Fraction(count) * y**area for area, count in head.items())
    ratio = LAMBDA * y
    if ratio >= 1:
        raise ValueError("tail majorant does not converge")
    tail = TAIL_COEFFICIENT * ratio**26 / (1 - ratio**2)
    return exact_head, tail, exact_head + tail


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def scaled_decimal(integer: int, places: int) -> str:
    sign = "-" if integer < 0 else ""
    digits = str(abs(integer)).rjust(places + 1, "0")
    return f"{sign}{digits[:-places]}.{digits[-places:]}"


def decimal_enclosure(value: Fraction, places: int = 40) -> list[str]:
    scale = 10**places
    lower = value.numerator * scale // value.denominator
    upper = (value.numerator * scale + value.denominator - 1) // value.denominator
    return [scaled_decimal(lower, places), scaled_decimal(upper, places)]


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")
        if not passed:
            raise AssertionError(name)

    census = enumerate_polycube_head()
    rooted_head = census["rooted_counts"]
    assert isinstance(rooted_head, dict)
    direct_counts, direct_stats = direct_plaquette_counts()

    raw_lower, raw_upper = exp_negative_taylor_bounds(2 * K0)
    y_lower, y_upper = outward_decimal_interval(raw_lower, raw_upper)
    head_upper, tail_upper, criterion_upper = criterion_terms(y_upper, rooted_head)
    margin = Fraction(1, 2) - criterion_upper

    previous_raw_lower, previous_raw_upper = exp_negative_taylor_bounds(2 * K_PREVIOUS_GRID)
    previous_y_lower, previous_y_upper = outward_decimal_interval(
        previous_raw_lower, previous_raw_upper
    )
    previous_head_lower, previous_tail_lower, previous_criterion_lower = criterion_terms(
        previous_y_lower, rooted_head
    )

    fv_hash = hashlib.sha256(FRIEDLI_VELENIK.read_bytes()).hexdigest()
    bonati_hash = hashlib.sha256(BONATI.read_bytes()).hexdigest()

    check(
        "fixed_polycube_level_counts",
        census["level_counts"] == EXPECTED_LEVEL_COUNTS,
        f"computed translation-quotiented fixed counts {census['level_counts']}",
    )
    check(
        "exact_head_counts",
        rooted_head == EXPECTED_HEAD,
        f"computed rooted contour counts {rooted_head}",
    )
    check("unit_cube", rooted_head[6] == 1, "N(6)=1")
    check(
        "direct_plaquette_crosscheck",
        all(direct_counts.get(area, 0) == rooted_head.get(area, 0) for area in (6, 8, 10)),
        f"independent direct closure search gives {direct_counts}; stats={direct_stats}",
    )
    check(
        "isoperimetric_completeness",
        A_MAX**3 == 216 * N_MAX**2 and A_MAX**3 < 216 * (N_MAX + 1) ** 2,
        "A^3>=216*n^2 implies every area<=24 interior has n<=8 cells",
    )
    check(
        "no_relevant_cavities",
        census["relevant_with_cavities"] == 0,
        f"exact flood fill found {census['relevant_with_cavities']} cavity-bearing candidates at area<=24",
    )
    check(
        "tail_growth_constant",
        11**11 < 29 * 10**10 and LAMBDA > 0,
        f"lambda=11^11/10^10={fraction_text(LAMBDA)}<29",
    )
    check(
        "exact_exp_enclosure",
        y_lower <= raw_lower < raw_upper <= y_upper and y_upper - y_lower <= Fraction(1, 10**49),
        f"50-place outward rational interval {fraction_text(y_lower)} .. {fraction_text(y_upper)}",
    )
    check(
        "tail_converges_at_K0",
        LAMBDA * y_upper < 1,
        f"lambda*y_upper decimal enclosure={decimal_enclosure(LAMBDA * y_upper, 30)}",
    )
    check(
        "certified_peierls_inequality",
        criterion_upper < Fraction(1, 2),
        f"upper enclosure={decimal_enclosure(criterion_upper)}; positive margin={decimal_enclosure(margin)}",
    )
    check(
        "previous_millionth_fails",
        previous_criterion_lower > Fraction(1, 2)
        and K0 - K_PREVIOUS_GRID == Fraction(1, 1000000),
        f"at 1.696897 lower enclosure={decimal_enclosure(previous_criterion_lower)}",
    )
    check(
        "does_not_improve_incumbent",
        K0 > INCUMBENT_UPPER,
        f"{fraction_text(K0)} is weaker than incumbent decimal 0.25273100985867",
    )
    check(
        "friedli_velenik_source_sha256",
        fv_hash == FRIEDLI_VELENIK_SHA256,
        f"sha256={fv_hash}",
    )
    check(
        "bonati_source_sha256",
        bonati_hash == BONATI_SHA256,
        f"sha256={bonati_hash}",
    )

    result = {
        "meta": {
            "script": SCRIPT,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "classification": "[THEOREM] finite-volume Peierls inequality; [COMPUTATION] exact finite census and rational certificate",
            "arithmetic": {
                "combinatorics": "exact Python integers",
                "transcendental": f"alternating Taylor enclosure through orders {TAYLOR_ORDER} and {TAYLOR_ORDER + 1}, exact Fraction arithmetic",
                "stored_exp_interval": f"rounded outward to denominator 10^{EXP_DECIMAL_PLACES}",
                "floating_point_used_in_certificate": False,
            },
            "sources": [
                {
                    "key": "friedli_velenik2017_ising_ch3",
                    "path": str(FRIEDLI_VELENIK.relative_to(ROOT)),
                    "sha256": fv_hash,
                    "location": "Section 3.7.2, especially Lemma 3.37 and Exercise 3.20",
                },
                {
                    "key": "bonati2014_peierls",
                    "path": str(BONATI.relative_to(ROOT)),
                    "sha256": bonati_hash,
                    "location": "Section 3, equations (17)-(23), comparison only",
                },
            ],
        },
        "data": {
            "model": "nearest-neighbor ferromagnetic Ising model on Z^3, K=beta*J",
            "contour_definition": "edge-connected outer dual-plaquette boundary of the filled minus component containing the origin",
            "exact_head": {
                "A_max": A_MAX,
                "isoperimetric_cell_cutoff": N_MAX,
                "N_of_A": {str(area): count for area, count in sorted(rooted_head.items())},
                "translation_quotiented_shape_counts": {
                    str(area): count
                    for area, count in sorted(census["unrooted_counts"].items())
                },
                "fixed_polycubes_by_cell_count": {
                    str(size): count for size, count in census["level_counts"].items()
                },
                "shape_sha256_by_area": {
                    str(area): digest
                    for area, digest in census["shape_sha256_by_area"].items()
                },
                "shape_counts_by_cells_and_area": {
                    str(size): {str(area): count for area, count in counts.items()}
                    for size, counts in census["by_cells_area"].items()
                },
                "direct_plaquette_search_A_le_10": {
                    "counts": {str(area): count for area, count in direct_counts.items()},
                    "stats": direct_stats,
                },
            },
            "tail": {
                "plaquette_edge_adjacency_max_degree": 12,
                "root_choices_bound": "at most A positive-x ray roots for an area-A contour",
                "universal_cover_animal_bound": "B_A = 12/(A-1)*binom(11*A,A-2) for A>=2",
                "binomial_majorant_derivation": "binom(11A,A-2) <= (11/10)^(11A)*10^(A-2) = lambda^A/100",
                "lambda_exact": fraction_text(LAMBDA),
                "lambda_decimal_enclosure": decimal_enclosure(LAMBDA, 20),
                "N_bound_A_ge_26": "N(A) <= (78/625)*lambda^A",
                "even_area_tail_at_y_exp_minus_2K": "(78/625)*(lambda*y)^26/(1-(lambda*y)^2)",
                "convergence_barrier": "requires lambda*exp(-2K)<1, i.e. K>(1/2)log(lambda); this barrier alone exceeds the incumbent upper endpoint",
            },
            "certificate": {
                "criterion": "sum_A N(A)*exp(-2*K*A) < 1/2",
                "K0_exact": fraction_text(K0),
                "K0_decimal": "1.696898",
                "grid_optimality_scope": "least multiple of 10^-6 certified by this exact-head plus stated-tail majorant",
                "exp_minus_2K0_interval_exact": [fraction_text(y_lower), fraction_text(y_upper)],
                "exp_minus_2K0_interval_decimal": [
                    decimal_enclosure(y_lower, 50)[0],
                    decimal_enclosure(y_upper, 50)[1],
                ],
                "head_upper_exact": fraction_text(head_upper),
                "head_upper_decimal": decimal_enclosure(head_upper),
                "tail_upper_exact": fraction_text(tail_upper),
                "tail_upper_decimal": decimal_enclosure(tail_upper),
                "sum_upper_exact": fraction_text(criterion_upper),
                "sum_upper_decimal": decimal_enclosure(criterion_upper),
                "margin_below_half_exact": fraction_text(margin),
                "margin_below_half_decimal": decimal_enclosure(margin),
                "previous_grid_point": {
                    "K_exact": fraction_text(K_PREVIOUS_GRID),
                    "criterion_lower_exact": fraction_text(previous_criterion_lower),
                    "criterion_lower_decimal": decimal_enclosure(previous_criterion_lower),
                    "head_lower_exact": fraction_text(previous_head_lower),
                    "tail_lower_exact": fraction_text(previous_tail_lower),
                },
                "theorem": "K_c <= 848449/500000",
                "improves_incumbent": False,
                "incumbent_upper_decimal": "0.25273100985867",
                "outcome": "certified negative result: exact areas through 24 do not overcome the rigorous tail growth constant",
            },
            "scope_and_limitations": [
                "The N(A) table is exact only for even A<=24; the isoperimetric cutoff makes that finite census complete.",
                "For A>=26 the table is replaced by a proved overcount of all connected plaquette animals; no equality is claimed.",
                "The resulting K0 is a valid upper bound but is much weaker than both the incumbent infrared bound and the previously recorded Bonati Peierls bound.",
                "The benchmark K_c approximation was not read or used by this experiment.",
            ],
        },
        "checks": checks,
    }
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULT_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL: {error}")
        raise
