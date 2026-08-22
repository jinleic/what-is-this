"""Standalone exact checks for the quantitative 3D Peierls certificate."""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import defaultdict, deque
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "bounds" / "peierls_upper.json"
FRIEDLI_VELENIK = ROOT / "sources" / "fulltext" / "friedli_velenik_ising_ch3.pdf"
FRIEDLI_VELENIK_SHA256 = "4809d134deb0a97bb37bcc306345d3ace659f7543bbcd8cfedde53c5dbe47045"
A_MAX = 24
N_MAX = 8
K0 = Fraction(848449, 500000)
K_PREVIOUS = Fraction(1696897, 1000000)
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
EXPECTED_HEAD = {6: 1, 8: 0, 10: 6, 12: 0, 14: 45, 16: 12, 18: 332, 20: 240, 22: 2538, 24: 3040}


def canonical(cells: set[tuple[int, int, int]]) -> tuple[tuple[int, int, int], ...]:
    minima = [min(cell[axis] for cell in cells) for axis in range(3)]
    return tuple(
        sorted(
            tuple(cell[axis] - minima[axis] for axis in range(3)) for cell in cells
        )
    )


def area(shape: tuple[tuple[int, int, int], ...]) -> int:
    cells = set(shape)
    return sum(
        tuple(cell[axis] + step[axis] for axis in range(3)) not in cells
        for cell in cells
        for step in DIRECTIONS
    )


def has_connected_complement(shape: tuple[tuple[int, int, int], ...]) -> bool:
    cells = set(shape)
    low = [min(cell[axis] for cell in cells) - 1 for axis in range(3)]
    high = [max(cell[axis] for cell in cells) + 1 for axis in range(3)]
    box = set(itertools.product(*(range(low[axis], high[axis] + 1) for axis in range(3))))
    outside = {tuple(low)}
    queue = deque(outside)
    while queue:
        cell = queue.popleft()
        for step in DIRECTIONS:
            neighbor = tuple(cell[axis] + step[axis] for axis in range(3))
            if neighbor in box and neighbor not in cells and neighbor not in outside:
                outside.add(neighbor)
                queue.append(neighbor)
    return outside | cells == box


def recompute_polycube_head() -> tuple[dict[int, int], dict[int, int]]:
    level = {((0, 0, 0),)}
    level_sizes: dict[int, int] = {}
    rooted: defaultdict[int, int] = defaultdict(int)
    for size in range(1, N_MAX + 1):
        level_sizes[size] = len(level)
        for shape in level:
            boundary = area(shape)
            if boundary <= A_MAX and has_connected_complement(shape):
                rooted[boundary] += size
        if size == N_MAX:
            break
        next_level = set()
        for shape in level:
            cells = set(shape)
            for cell in shape:
                for step in DIRECTIONS:
                    neighbor = tuple(cell[axis] + step[axis] for axis in range(3))
                    if neighbor not in cells:
                        next_level.add(canonical(cells | {neighbor}))
        level = next_level
    return (
        {boundary: rooted.get(boundary, 0) for boundary in range(6, A_MAX + 1, 2)},
        level_sizes,
    )


def face_edges(face: tuple[int, int, int, int]) -> tuple[tuple[int, int, int, int], ...]:
    normal, *center = face
    tangent = [axis for axis in range(3) if axis != normal]
    result = []
    for sign in (-1, 1):
        edge = center.copy()
        edge[tangent[1]] += sign
        result.append((tangent[0], *edge))
        edge = center.copy()
        edge[tangent[0]] += sign
        result.append((tangent[1], *edge))
    return tuple(result)


def faces_at_edge(edge: tuple[int, int, int, int]) -> tuple[tuple[int, int, int, int], ...]:
    direction, *center = edge
    transverse = [axis for axis in range(3) if axis != direction]
    result = []
    for normal, tangent in ((transverse[0], transverse[1]), (transverse[1], transverse[0])):
        for sign in (-1, 1):
            face = center.copy()
            face[tangent] += sign
            result.append((normal, *face))
    return tuple(result)


def plus_face(cell: tuple[int, int, int], axis: int) -> tuple[int, int, int, int]:
    center = [2 * coordinate for coordinate in cell]
    center[axis] += 1
    return (axis, *center)


def boundary_faces(cells: set[tuple[int, int, int]]) -> frozenset[tuple[int, int, int, int]]:
    result = set()
    for cell in cells:
        for axis in range(3):
            lower = tuple(cell[index] - (index == axis) for index in range(3))
            upper = tuple(cell[index] + (index == axis) for index in range(3))
            if lower not in cells:
                result.add(plus_face(lower, axis))
            if upper not in cells:
                result.add(plus_face(cell, axis))
    return frozenset(result)


def enclosed_cells(surface: frozenset[tuple[int, int, int, int]]) -> set[tuple[int, int, int]]:
    adjacent = []
    for face in surface:
        normal, *center = face
        lower = [coordinate // 2 for coordinate in center]
        lower[normal] = (center[normal] - 1) // 2
        upper = lower.copy()
        upper[normal] += 1
        adjacent.extend((tuple(lower), tuple(upper)))
    low = [min(cell[axis] for cell in adjacent) - 1 for axis in range(3)]
    high = [max(cell[axis] for cell in adjacent) + 1 for axis in range(3)]
    box = set(itertools.product(*(range(low[axis], high[axis] + 1) for axis in range(3))))
    outside = {tuple(low)}
    queue = deque(outside)
    while queue:
        cell = queue.popleft()
        for axis in range(3):
            for sign in (-1, 1):
                neighbor_list = list(cell)
                neighbor_list[axis] += sign
                neighbor = tuple(neighbor_list)
                if neighbor not in box or neighbor in outside:
                    continue
                lower = neighbor if sign < 0 else cell
                if plus_face(lower, axis) not in surface:
                    outside.add(neighbor)
                    queue.append(neighbor)
    return box - outside


def direct_counts(max_area: int = 10) -> dict[int, int]:
    roots = [(0, 2 * offset + 1, 0, 0) for offset in range(max_area)]
    stack = [(frozenset({root}), frozenset(face_edges(root))) for root in roots]
    seen = set()
    closed = set()
    while stack:
        chosen, odd = stack.pop()
        if chosen in seen:
            continue
        seen.add(chosen)
        if not odd:
            closed.add(chosen)
            if 1 < len(chosen) < max_area:
                neighbors = set()
                for face in chosen:
                    for edge in face_edges(face):
                        neighbors.update(faces_at_edge(edge))
                for face in neighbors - set(chosen):
                    stack.append((chosen | {face}, frozenset(face_edges(face))))
            continue
        if len(chosen) >= max_area:
            continue
        edge = min(odd)
        for face in faces_at_edge(edge):
            if face in chosen:
                continue
            updated = set(odd)
            for item in face_edges(face):
                if item in updated:
                    updated.remove(item)
                else:
                    updated.add(item)
            stack.append((chosen | {face}, frozenset(updated)))
    valid: defaultdict[int, set[frozenset[tuple[int, int, int, int]]]] = defaultdict(set)
    for surface in closed:
        inside = enclosed_cells(surface)
        if (0, 0, 0) in inside and boundary_faces(inside) == surface:
            valid[len(surface)].add(surface)
    return {boundary: len(valid.get(boundary, set())) for boundary in (6, 8, 10)}


def exp_bounds(x: Fraction, order: int = 62) -> tuple[Fraction, Fraction]:
    """Independent exact Taylor enclosure, using two more terms than experiment."""
    assert order % 2 == 0 and Fraction(order + 2) > x
    term = Fraction(1)
    partial = term
    upper = None
    for index in range(1, order + 2):
        term *= -x / index
        partial += term
        if index == order:
            upper = partial
    assert upper is not None and 0 < partial < upper
    return partial, upper


def round_out(lower: Fraction, upper: Fraction, places: int = 55) -> tuple[Fraction, Fraction]:
    scale = 10**places
    return (
        Fraction(lower.numerator * scale // lower.denominator, scale),
        Fraction((upper.numerator * scale + upper.denominator - 1) // upper.denominator, scale),
    )


def majorant(y: Fraction, head: dict[int, int]) -> Fraction:
    exact_head = sum(Fraction(count) * y**boundary for boundary, count in head.items())
    ratio = LAMBDA * y
    assert ratio < 1
    return exact_head + TAIL_COEFFICIENT * ratio**26 / (1 - ratio**2)


def main() -> None:
    failures = 0

    def check(name: str, condition: bool, detail: str) -> None:
        nonlocal failures
        print(f"{'PASS' if condition else 'FAIL'} {name}: {detail}")
        failures += not condition

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    recomputed_head, levels = recompute_polycube_head()
    check("polycube_head", recomputed_head == EXPECTED_HEAD, str(recomputed_head))
    check(
        "polycube_levels",
        levels == {1: 1, 2: 3, 3: 15, 4: 86, 5: 534, 6: 3481, 7: 23502, 8: 162913},
        str(levels),
    )
    independently_direct = direct_counts()
    check(
        "direct_plaquette_search",
        independently_direct == {6: 1, 8: 0, 10: 6},
        str(independently_direct),
    )
    check(
        "two_enumerators_agree",
        all(independently_direct[boundary] == recomputed_head[boundary] for boundary in (6, 8, 10)),
        "cell-BFS and direct even-edge closure search agree through area 10",
    )
    check(
        "isoperimetric_cutoff",
        A_MAX**3 == 216 * N_MAX**2 < 216 * (N_MAX + 1) ** 2,
        "24^3=216*8^2<216*9^2",
    )

    representative = (0, 1, 0, 0)
    adjacency = set()
    for edge in face_edges(representative):
        adjacency.update(faces_at_edge(edge))
    adjacency.remove(representative)
    check("plaquette_degree", len(adjacency) == 12, f"computed degree {len(adjacency)}")

    animal_coefficients_ok = True
    for number in range(2, 101):
        universal_cover_count = Fraction(12, number - 1) * comb(11 * number, number - 2)
        exact_exponential_bound = Fraction(12, 100 * (number - 1)) * LAMBDA**number
        animal_coefficients_ok &= universal_cover_count <= exact_exponential_bound
    check(
        "tail_coefficient_sanity",
        animal_coefficients_ok and TAIL_COEFFICIENT == Fraction(12 * 26, 100 * 25),
        "exact coefficient inequalities checked for A=2..100; tail prefactor is 78/625",
    )

    k_lower, k_upper = round_out(*exp_bounds(2 * K0))
    previous_lower, previous_upper = round_out(*exp_bounds(2 * K_PREVIOUS))
    candidate_upper = majorant(k_upper, recomputed_head)
    previous_bound_lower = majorant(previous_lower, recomputed_head)
    check(
        "exact_final_inequality",
        candidate_upper < Fraction(1, 2),
        f"exact positive margin has {(Fraction(1, 2) - candidate_upper).numerator.bit_length()} numerator bits",
    )
    check(
        "millionth_grid_minimum",
        previous_bound_lower > Fraction(1, 2) and K0 - K_PREVIOUS == Fraction(1, 1000000),
        "1.696897 fails while 1.696898 passes the same majorant",
    )
    check(
        "stored_counts",
        payload["data"]["exact_head"]["N_of_A"]
        == {str(boundary): count for boundary, count in recomputed_head.items()},
        "JSON head equals independently recomputed head",
    )
    check(
        "stored_K0",
        payload["data"]["certificate"]["K0_exact"] == "848449/500000",
        payload["data"]["certificate"]["K0_exact"],
    )
    check(
        "stored_checks",
        all(item["passed"] for item in payload["checks"]),
        f"{len(payload['checks'])} generated checks are true",
    )
    source_hash = hashlib.sha256(FRIEDLI_VELENIK.read_bytes()).hexdigest()
    check("source_hash", source_hash == FRIEDLI_VELENIK_SHA256, source_hash)

    if failures:
        raise SystemExit(f"FAIL: {failures} Peierls checks failed")
    print("PASS: quantitative Peierls certificate independently verified")


if __name__ == "__main__":
    main()
