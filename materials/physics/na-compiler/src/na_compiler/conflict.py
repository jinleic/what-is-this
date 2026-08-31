"""Conflict-graph model of the strict AOD move-compatibility relation.

Relation pinned from QMAP IndependentSetRouter.cpp (read 2026-08-29):
moves m1=(xs1,ys1,xt1,yt1), m2=(xs2,ys2,xt2,yt2) are *incompatible* (cannot
share a batch) unless ALL of:
  (xs1 == xs2) <=> (xt1 == xt2)
  (ys1 == ys2) <=> (yt1 == yt2)
  (xs1 <  xs2) <=> (xt1 <  xt2)
  (ys1 <  ys2) <=> (yt1 <  yt2)
Note: equality collapsing is forbidden by 1-2 (two atoms starting in the same
row must end in the same row — 'preservation' — and distinct rows must not be
merged); order inversions are forbidden by 3-4 ('non-crossing').
Ghost-spot constraints are NOT part of QMAP's strict router; optionally added
by ghost_spots=True (pick-up/drop-off row/column disjointness).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .architecture import Architecture, Site


@dataclass(frozen=True)
class Move:
    atom: int
    start: tuple[float, float]
    target: tuple[float, float]
    start_site: Site
    target_site: Site

    @property
    def vec(self) -> tuple[float, float, float, float]:
        x0, y0 = self.start
        x1, y1 = self.target
        return (x0, y0, x1, y1)


def compatible(m1: Move, m2: Move) -> bool:
    x0a, y0a, x1a, y1a = m1.vec
    x0b, y0b, x1b, y1b = m2.vec
    if (x0a == x0b) != (x1a == x1b):
        return False
    if (y0a == y0b) != (y1a == y1b):
        return False
    if (x0a < x0b) != (x1a < x1b):
        return False
    if (y0a < y0b) != (y1a < y1b):
        return False
    return True


def moves_from_trajectory(
    start: dict[str, tuple[float, float]],
    target: dict[str, tuple[float, float]],
    arch: Architecture | None = None,
) -> list[Move]:
    atoms = sorted(start)
    assert set(atoms) == set(target), "atom sets differ across transition"
    out = []
    for i, a in enumerate(atoms):
        if start[a] != target[a]:
            out.append(
                Move(
                    atom=i,
                    start=start[a],
                    target=target[a],
                    start_site=Site("s", 0, -1, -1),  # zone sites resolved by caller
                    target_site=Site("s", 0, -1, -1),
                )
            )
    return out


def conflict_graph(moves: list[Move], ghost_spots: bool = False) -> list[tuple[int, int]]:
    """Edges between incompatible moves; optionally add ghost-spot conflicts:
    two moves conflict if they share a start grid row/col with an
    unloading/loading atom elsewhere on that row/col (disjoint pick-up and
    drop-off row/column sets batch-wide is the paper's constraint; modeled
    pairwise here as in NALAC)."""
    edges = set()
    for i, j in combinations(range(len(moves)), 2):
        if not compatible(moves[i], moves[j]):
            edges.add((i, j))
    if ghost_spots:
        # pairwise ghost-spot: move i and j conflict if start rows/cols or
        # target rows/cols cross between the two moves' endpoints
        for i, j in combinations(range(len(moves)), 2):
            (x0a, y0a, x1a, y1a) = moves[i].vec
            (x0b, y0b, x1b, y1b) = moves[j].vec
            row_or_col_cross = (
                (y0a == y0b) != (y1a == y1b)  # handled by compatibility
            )
            # ghost spot addition: x1a == x0b or x1b == x0a etc. is allowed by
            # QMAP but excluded under ghost-spot discipline
            if (x1a == x0b and y0a == y0b) or (x1b == x0a and y0a == y0b):
                edges.add((i, j))
            if row_or_col_cross:
                edges.add((i, j))
    return sorted(edges)


def chromatic_number_upper_bound(n_moves: int, edges: list[tuple[int, int]]) -> int:
    """Greedy DSATUR upper bound, for solver warm starts / batch counting."""
    if n_moves == 0:
        return 0
    adj = [set() for _ in range(n_moves)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    colors: dict[int, int] = {}
    import heapq

    sat = [set() for _ in range(n_moves)]
    heap = [(-len(adj[v]), -v, v) for v in range(n_moves)]
    heapq.heapify(heap)
    while heap:
        _, _, v = heapq.heappop(heap)
        if v in colors:
            continue
        c = 0
        while c in sat[v]:
            c += 1
        colors[v] = c
        for w in adj[v]:
            if w not in colors:
                sat[w].add(c)
    return max(colors.values()) + 1 if colors else 0
