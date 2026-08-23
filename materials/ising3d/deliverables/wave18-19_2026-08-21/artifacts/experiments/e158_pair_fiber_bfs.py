#!/usr/bin/env python3
"""Exact affine Pauli-string BFS for every weight-two fiber at L=2,...,7.

The BFS uses only sound commutator productions.  A complete source pair fiber
is transferred across one graph bond by the anticommuting half-fiber; toggling
with the X generator at the newly occupied endpoint restores the one direction
lost by the hyperplane filter.  The initial complete fiber is supplied by the
bounded square move proved and checked in e157.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import time
from collections import deque

from e157_saturation_local_move import (
    checkerboard_mask,
    grid,
    parity,
    square_triple_witnesses,
)


def gf2_rank(vectors, stop_at: int | None = None) -> int:
    basis: dict[int, int] = {}
    for value in vectors:
        row = value
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                if stop_at is not None and len(basis) == stop_at:
                    return stop_at
                break
            row ^= basis[pivot]
    return len(basis)


def pair_law(x: int, z: int, L: int) -> bool:
    n = 2 * L
    J = (1 << n) - 1
    d = checkerboard_mask(L)
    return parity((J ^ z) & x) == (1 ^ parity(z & d))


def pair_digest(fibers: dict[int, set[int]], n: int) -> str:
    width = max(1, (2 * n + 7) // 8)
    digest = hashlib.sha256()
    for z in sorted(fibers):
        for x in sorted(fibers[z]):
            digest.update((x | (z << n)).to_bytes(width, "little"))
    return digest.hexdigest()


def base_fiber_L2() -> tuple[int, set[int], dict[int, dict[str, int]]]:
    """Fill one L=2 edge fiber from its edge seed and one local triple."""
    L = 2
    n = 4
    p, q, r, s = 0, 1, 2, 3
    witnesses = square_triple_witnesses(p, q, r, s, n)
    z = (1 << p) ^ (1 << r)
    edge_points = {
        0,
        1 << p,
        1 << r,
        (1 << p) ^ (1 << r),
    }
    # The triple omitting r is p+q+s and intersects z only at p.
    toggle = witnesses[r]["output"] & ((1 << n) - 1)
    if parity(toggle & z) != 1:
        raise AssertionError("the base square triple must anticommute with the edge fiber")
    points = edge_points | {x ^ toggle for x in edge_points}
    if len(points) != 1 << (n - 1):
        raise AssertionError("the L=2 base fiber did not acquire its third direction")
    if not all(pair_law(x, z, L) for x in points):
        raise AssertionError("the L=2 base construction escaped the law fiber")
    if gf2_rank([x ^ min(points) for x in points]) != n - 1:
        raise AssertionError("the L=2 base tangent has deficient rank")
    return z, points, witnesses


def extend_seed_fiber(
    old_L: int, old_z: int, old_points: set[int]
) -> tuple[int, set[int], dict[int, dict[str, int]]]:
    """Raise a full old boundary-edge fiber by two directions at L+1."""
    if old_L < 2:
        raise ValueError("the induction starts from L=2")
    new_L = old_L + 1
    n = 2 * new_L
    t = 2 * (old_L - 2)
    p, q = 2 * (old_L - 1), 2 * (old_L - 1) + 1
    r, s = 2 * old_L, 2 * old_L + 1
    expected_old_z = (1 << t) ^ (1 << p)
    if old_z != expected_old_z:
        raise AssertionError("the carried fiber is not the old top boundary edge")
    if len(old_points) != 1 << (2 * old_L - 1):
        raise AssertionError("the carried old fiber is not complete")

    witnesses = square_triple_witnesses(p, q, r, s, n)
    # q1=p+q+r (omit s), q2=p+q+s (omit r).  Both meet {t,p} once.
    q1 = witnesses[s]["output"] & ((1 << n) - 1)
    q2 = witnesses[r]["output"] & ((1 << n) - 1)
    if parity(q1 & old_z) != 1 or parity(q2 & old_z) != 1:
        raise AssertionError("both extension triples must toggle the carried fiber")

    points = {
        x ^ extra
        for x in old_points
        for extra in (0, q1, q2, q1 ^ q2)
    }
    expected_size = 1 << (n - 1)
    if len(points) != expected_size:
        raise AssertionError("the two square directions are not independent modulo old support")
    if not all(pair_law(x, old_z, new_L) for x in points):
        raise AssertionError("the extension construction escaped the new law fiber")
    anchor = min(points)
    if gf2_rank([x ^ anchor for x in points]) != n - 1:
        raise AssertionError("the extended seed fiber has deficient tangent rank")
    return old_z, points, witnesses


def spread_pair_fibers(
    L: int, seed_z: int, seed_points: set[int], deadline: float
) -> tuple[dict[int, set[int]], list[dict[str, int]]]:
    """BFS on the two-token graph, retaining every exact Pauli x-label."""
    n = 2 * L
    _, adjacency = grid(L)
    expected_size = 1 << (n - 1)
    if len(seed_points) != expected_size:
        raise AssertionError("pair BFS needs one complete seed fiber")
    if not all(pair_law(x, seed_z, L) for x in seed_points):
        raise AssertionError("seed points violate the target pair law")

    fibers: dict[int, set[int]] = {seed_z: set(seed_points)}
    queue: deque[int] = deque([seed_z])
    transitions: list[dict[str, int]] = []
    while queue:
        if time.process_time() > deadline:
            raise TimeoutError("NON-DECISIVE: e158 process-time budget expired")
        source_z = queue.popleft()
        source_vertices = [i for i in range(n) if (source_z >> i) & 1]
        if len(source_vertices) != 2:
            raise AssertionError("pair BFS encountered a non-pair support")
        for position, removed in enumerate(source_vertices):
            fixed = source_vertices[1 - position]
            for added in adjacency[removed]:
                if added == fixed:
                    continue
                target_z = source_z ^ (1 << removed) ^ (1 << added)
                if target_z in fibers:
                    continue
                bond = (1 << removed) ^ (1 << added)
                transferred = {
                    x for x in fibers[source_z] if parity(x & bond) == 1
                }
                if len(transferred) != expected_size // 2:
                    raise AssertionError("bond filter did not cut a complete source fiber in half")
                # The bond commutator keeps x and changes z.  X_added then supplies
                # the missing tangent direction in the target fiber.
                target_points = transferred | {
                    x ^ (1 << added) for x in transferred
                }
                if len(target_points) != expected_size:
                    raise AssertionError("endpoint-X recovery did not fill the target fiber")
                if not all(pair_law(x, target_z, L) for x in target_points):
                    raise AssertionError("a token-graph transfer escaped the law fiber")
                fibers[target_z] = target_points
                queue.append(target_z)
                transitions.append(
                    {
                        "source_z": source_z,
                        "target_z": target_z,
                        "bond": bond,
                        "removed": removed,
                        "added": added,
                        "transferred_half": len(transferred),
                    }
                )

    expected_pairs = math.comb(n, 2)
    if len(fibers) != expected_pairs:
        raise AssertionError(
            f"two-token BFS covered {len(fibers)}/{expected_pairs} pair supports"
        )
    if len(transitions) != expected_pairs - 1:
        raise AssertionError("the pair-fiber BFS certificate is not a spanning tree")
    return fibers, transitions


def run_pair_bfs(L_max: int, cpu_budget: float) -> list[dict]:
    if L_max < 2:
        raise ValueError("L_max must be at least 2")
    if cpu_budget <= 0:
        raise ValueError("cpu_budget must be positive")
    deadline = time.process_time() + cpu_budget
    seed_z, seed_points, _ = base_fiber_L2()
    rows: list[dict] = []
    for L in range(2, L_max + 1):
        if L > 2:
            seed_z, seed_points, _ = extend_seed_fiber(L - 1, seed_z, seed_points)
        fibers, transitions = spread_pair_fibers(L, seed_z, seed_points, deadline)
        n = 2 * L
        expected_size = 1 << (n - 1)
        ranks = set()
        for points in fibers.values():
            anchor = min(points)
            ranks.add(
                gf2_rank(
                    (x ^ anchor for x in points),
                    stop_at=n - 1,
                )
            )
        if ranks != {n - 1}:
            raise AssertionError(f"L={L} has incomplete pair-fiber tangent ranks: {ranks}")
        row = {
            "L": L,
            "n": n,
            "pair_supports": len(fibers),
            "pair_fiber_size": expected_size,
            "pair_strings": len(fibers) * expected_size,
            "tangent_rank": n - 1,
            "transport_edges": len(transitions),
            "pair_labels_sha256": pair_digest(fibers, n),
        }
        rows.append(row)

        # Carry the new top boundary edge into the next induction step.
        if L < L_max:
            top_left = 2 * (L - 2)
            top_right = 2 * (L - 1)
            seed_z = (1 << top_left) ^ (1 << top_right)
            seed_points = set(fibers[seed_z])
    return rows


def main() -> bool:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-L", type=int, default=7)
    parser.add_argument("--cpu-budget", type=float, default=120.0)
    args = parser.parse_args()
    rows = run_pair_bfs(args.max_L, args.cpu_budget)
    for row in rows:
        print(
            f"L={row['L']} pairs={row['pair_supports']} "
            f"fiber={row['pair_fiber_size']} rank={row['tangent_rank']} "
            f"strings={row['pair_strings']}"
        )
    print("PASS e158 exact pair-fiber Pauli BFS")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
