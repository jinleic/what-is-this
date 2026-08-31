"""Gate B: exhaustive-BFS hardest permutation instances on small grids.

Semantics (pinned in pre_statement.md section 4):
- fixed S x S grid of sites, n atoms, one atom per site;
- a batch = any set of singleton moves (atom -> currently free site) that is
  pairwise compatible under the strict relation AND has pairwise-distinct
  targets; all targets must be free in the CURRENT state (conservative
  drop-into-free-trap semantics, matching pre_statement section 1);
- BFS over states from the identity placement gives exact min-batches for
  every reachable state; the exhaustion log (node/visit counts + the dist
  dict) is the certificate.

Correctness invariants tested in smoke:
- reachable states are always permutations of DISTINCT sites (no two atoms
  share a site), computed on the FIXED grid (never growing);
- grid-side is passed explicitly and never re-derived from state contents.
"""

from __future__ import annotations

from itertools import combinations

from .conflict import Move, compatible


def _rc(site: int, grid_side: int) -> tuple[int, int]:
    return divmod(site, grid_side)


def legal_batches(
    state: tuple[int, ...], grid_side: int
) -> list[tuple[tuple[int, int], ...]]:
    """All legal batches for one state on the fixed grid. Alpha-sorted moves
    for determinism. n_sites = grid_side^2."""
    nsites = grid_side * grid_side
    occupied = set(state)
    # candidate singleton moves (atom, target), target free and != own site
    cands: list[tuple[int, int]] = []
    for atom, site in enumerate(state):
        for tgt in range(nsites):
            if tgt not in occupied:
                cands.append((atom, tgt))

    def rc_move(a: int, tgt: int) -> Move:
        x0, y0 = _rc(state[a], grid_side)
        xt, yt = _rc(tgt, grid_side)
        return Move(a, (x0, y0), (xt, yt), None, None)

    res: list[tuple[tuple[int, int], ...]] = []

    def extend(cur: list[tuple[int, int]], start_idx: int) -> None:
        res.append(tuple(cur))
        for idx in range(start_idx, len(cands)):
            a, tgt = cands[idx]
            if any(a == a2 for a2, _ in cur):
                continue
            if any(tgt == t2 for _, t2 in cur):
                continue
            m1 = rc_move(a, tgt)
            ok = True
            for a2, t2 in cur:
                m2 = rc_move(a2, t2)
                if not compatible(m1, m2):
                    ok = False
                    break
            if ok:
                cur.append(cands[idx])
                extend(cur, idx + 1)
                cur.pop()

    # NOTE: generating ALL compatible subsets is exponential in |cands|;
    # for gate B sizes (n<=8, <=16 free targets) this is the dominant cost
    # and is itself part of the scaling table.
    extend([], 0)
    return res


def apply_batch(
    state: tuple[int, ...], batch: tuple[tuple[int, int], ...]
) -> tuple[int, ...]:
    new = list(state)
    for atom, tgt in batch:
        new[atom] = tgt
    return tuple(new)


def bfs_all_distances(
    start: tuple[int, ...], grid_side: int, max_depth: int | None = None,
    max_states: int | None = None,
):
    """Exhaustive BFS from start. Returns (dist, frontier_sizes, visited,
    aborted_flag). dist maps state -> min batches. The returned structure is
    the certificate for the reached subspace."""
    dist = {start: 0}
    frontier = [start]
    frontier_sizes = [1]
    depth = 0
    aborted = False
    batch_cache: dict[tuple[int, ...], list] = {}
    while frontier:
        if max_depth is not None and depth >= max_depth:
            aborted = True
            break
        if max_states is not None and len(dist) > max_states:
            aborted = True
            break
        nxt = []
        seen_next = set()
        for state in frontier:
            batches = batch_cache.get(state)
            if batches is None:
                batches = legal_batches(state, grid_side)
                batch_cache[state] = batches
            for batch in batches:
                s2 = apply_batch(state, batch)
                if s2 not in dist and s2 not in seen_next:
                    seen_next.add(s2)
                    dist[s2] = depth + 1
                    nxt.append(s2)
        frontier = nxt
        depth += 1
        frontier_sizes.append(len(nxt))
    return dist, frontier_sizes, len(dist), aborted


def enumerate_reachable(
    n_atoms: int, grid_side: int, max_depth: int | None = None,
    max_states: int | None = None,
):
    start = tuple(range(n_atoms))
    return bfs_all_distances(start, grid_side, max_depth, max_states)


def hardest_states(dist: dict, k: int = 5) -> list[tuple[tuple[int, ...], int]]:
    items = sorted(dist.items(), key=lambda kv: (-kv[1], kv[0]))
    return items[:k]
