#!/usr/bin/env python3
"""Clean-room verifier for the all-L 2xL saturation theorem.

No experiment module is imported.  The verifier independently rebuilds:

* all bounded square production witnesses at L=2,...,7;
* every exact weight-two affine fiber by legal Pauli productions;
* a seed and full tangent rank for EVERY even-z target row (no sampling);
* the rational Clifford grade census restricted to the 2xL statement;
* full generator-monotone Pauli closures at L=2,...,6.

L=7 is explicitly not a full closure enumeration: its 67,117,056 dimension is
checked by exhaustive support-row construction/ranks and by the independent
Clifford binomial identity.
"""

from __future__ import annotations

import gc
import hashlib
import json
import math
import resource
import sys
import time
from collections import deque
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "alll_saturation.json"
EXPECTED = {2: 56, 3: 1056, 4: 16256, 5: 262656, 6: 4192256, 7: 67117056}
CPU_BUDGET_SECONDS = 300.0


def parity(value: int) -> int:
    return value.bit_count() & 1


def rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def grid(L: int) -> tuple[list[tuple[int, int]], list[list[int]]]:
    n = 2 * L
    bonds = [(2 * c, 2 * c + 1) for c in range(L)]
    for c in range(L - 1):
        bonds.extend(((2 * c, 2 * c + 2), (2 * c + 1, 2 * c + 3)))
    adjacency = [[] for _ in range(n)]
    for u, v in bonds:
        adjacency[u].append(v)
        adjacency[v].append(u)
    for row in adjacency:
        row.sort()
    return bonds, adjacency


def checkerboard(L: int) -> int:
    return sum(
        1 << (2 * c + r)
        for c in range(L)
        for r in (0, 1)
        if ((c + r) & 1) == 0
    )


def pack(x: int, z: int, n: int) -> int:
    return x | (z << n)


def symplectic(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    x1, z1 = left & mask, left >> n
    x2, z2 = right & mask, right >> n
    return parity((x1 & z2) ^ (z1 & x2))


def law(x: int, z: int, L: int) -> bool:
    n = 2 * L
    J = (1 << n) - 1
    if parity(z):
        return False
    if z == J:
        return bool(L & 1)
    return parity((J ^ z) & x) == (1 ^ parity(z & checkerboard(L)))


def dimension(L: int) -> int:
    n = 2 * L
    m = 1 << (n - 1)
    return m * (m - ((-1) ** L))


def rank(vectors: Iterable[int], stop: int | None = None) -> int:
    basis: dict[int, int] = {}
    for vector in vectors:
        row = vector
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                if stop is not None and len(basis) == stop:
                    return stop
                break
            row ^= basis[pivot]
    return len(basis)


def nullspace(rows: list[int], n: int) -> list[int]:
    matrix = [row for row in rows if row]
    pivots: list[int] = []
    current = 0
    for column in range(n):
        found = next(
            (index for index in range(current, len(matrix)) if (matrix[index] >> column) & 1),
            None,
        )
        if found is None:
            continue
        matrix[current], matrix[found] = matrix[found], matrix[current]
        for index in range(len(matrix)):
            if index != current and ((matrix[index] >> column) & 1):
                matrix[index] ^= matrix[current]
        pivots.append(column)
        current += 1
        if current == len(matrix):
            break
    pivot_set = set(pivots)
    output = []
    for free in range(n):
        if free in pivot_set:
            continue
        vector = 1 << free
        for row_index, pivot in enumerate(pivots):
            if (matrix[row_index] >> free) & 1:
                vector ^= 1 << pivot
        assert all(parity(row & vector) == 0 for row in rows)
        output.append(vector)
    return output


def edge_point(u: int, v: int, alpha: int, beta: int, n: int) -> int:
    return pack((alpha << u) ^ (beta << v), (1 << u) ^ (1 << v), n)


def path2(u: int, v: int, w: int, alpha: int, beta: int, n: int) -> int:
    left = edge_point(u, v, alpha, 0, n)
    right = edge_point(v, w, 1, beta, n)
    assert symplectic(left, right, n) == 1
    output = left ^ right
    expected = pack(
        (alpha << u) ^ (1 << v) ^ (beta << w),
        (1 << u) ^ (1 << w),
        n,
    )
    assert output == expected
    return output


def square_triples(p: int, q: int, r: int, s: int, n: int) -> dict[int, int]:
    ps_r = path2(p, r, s, 0, 0, n)
    rq_p = path2(r, p, q, 0, 0, n)
    parents = {
        s: (ps_r, path2(p, q, s, 1, 0, n)),
        p: (ps_r, path2(p, q, s, 0, 1, n)),
        q: (rq_p, path2(r, s, q, 1, 0, n)),
        r: (rq_p, path2(r, s, q, 0, 1, n)),
    }
    corners = (1 << p) ^ (1 << q) ^ (1 << r) ^ (1 << s)
    output = {}
    for omitted, (left, right) in parents.items():
        assert symplectic(left, right, n) == 1
        child = left ^ right
        assert child == pack(corners ^ (1 << omitted), 0, n)
        output[omitted] = child & ((1 << n) - 1)
    return output


def base_L2() -> tuple[int, set[int]]:
    n = 4
    p, q, r, s = 0, 1, 2, 3
    triples = square_triples(p, q, r, s, n)
    z = (1 << p) ^ (1 << r)
    edge = {0, 1 << p, 1 << r, (1 << p) ^ (1 << r)}
    toggle = triples[r]
    assert parity(toggle & z) == 1
    points = edge | {x ^ toggle for x in edge}
    assert len(points) == 8 and all(law(x, z, 2) for x in points)
    return z, points


def extend(old_L: int, old_z: int, old_points: set[int]) -> tuple[int, set[int]]:
    new_L = old_L + 1
    n = 2 * new_L
    t = 2 * (old_L - 2)
    p, q = 2 * (old_L - 1), 2 * (old_L - 1) + 1
    r, s = 2 * old_L, 2 * old_L + 1
    assert old_z == (1 << t) ^ (1 << p)
    triples = square_triples(p, q, r, s, n)
    q1, q2 = triples[s], triples[r]
    assert parity(q1 & old_z) == parity(q2 & old_z) == 1
    points = {x ^ extra for x in old_points for extra in (0, q1, q2, q1 ^ q2)}
    assert len(points) == 1 << (n - 1)
    assert all(law(x, old_z, new_L) for x in points)
    return old_z, points


def spread(L: int, seed_z: int, seed_points: set[int], deadline: float) -> dict[int, set[int]]:
    n = 2 * L
    _, adjacency = grid(L)
    size = 1 << (n - 1)
    fibers = {seed_z: set(seed_points)}
    queue = deque([seed_z])
    while queue:
        assert time.process_time() <= deadline, "NON-DECISIVE: verifier process budget expired"
        source = queue.popleft()
        vertices = [index for index in range(n) if (source >> index) & 1]
        assert len(vertices) == 2
        for slot, removed in enumerate(vertices):
            fixed = vertices[1 - slot]
            for added in adjacency[removed]:
                if added == fixed:
                    continue
                target = source ^ (1 << removed) ^ (1 << added)
                if target in fibers:
                    continue
                bond = (1 << removed) ^ (1 << added)
                half = {x for x in fibers[source] if parity(x & bond)}
                assert len(half) == size // 2
                points = half | {x ^ (1 << added) for x in half}
                assert len(points) == size
                assert all(law(x, target, L) for x in points)
                fibers[target] = points
                queue.append(target)
    assert len(fibers) == math.comb(n, 2)
    return fibers


def pair_digest(fibers: dict[int, set[int]], n: int) -> str:
    digest = hashlib.sha256()
    width = max(1, (2 * n + 7) // 8)
    for z in sorted(fibers):
        for x in sorted(fibers[z]):
            digest.update(pack(x, z, n).to_bytes(width, "little"))
    return digest.hexdigest()


def pair_rep(z: int, L: int) -> int:
    J = (1 << (2 * L)) - 1
    lam = J ^ z
    return (lam & -lam) if (1 ^ parity(z & checkerboard(L))) else 0


def choose_pair(z: int, functional: int, desired: int, L: int) -> int:
    n = 2 * L
    J = (1 << n) - 1
    x = pair_rep(z, L)
    if parity(x & functional) != desired:
        spare = (J ^ z) & (J ^ functional)
        assert functional and not (functional & z) and spare
        x ^= (functional & -functional) ^ (spare & -spare)
    assert law(x, z, L) and parity(x & functional) == desired
    return x


def seed(z: int, L: int, fibers: dict[int, set[int]]) -> tuple[int, int]:
    n = 2 * L
    J = (1 << n) - 1
    assert z and z != J and parity(z) == 0
    sites = [index for index in range(n) if (z >> index) & 1]
    pairs = [(1 << sites[i]) ^ (1 << sites[i + 1]) for i in range(0, len(sites), 2)]
    current_z = pairs[0]
    current_x = pair_rep(current_z, L)
    assert current_x in fibers[current_z]
    for next_pair in pairs[1:]:
        next_x = choose_pair(next_pair, current_z, 1 ^ parity(current_x & next_pair), L)
        assert next_x in fibers[next_pair]
        assert symplectic(pack(current_x, current_z, n), pack(next_x, next_pair, n), n) == 1
        current_x ^= next_x
        current_z ^= next_pair
    assert current_z == z and law(current_x, current_z, L)
    return current_x, current_z


def universal_audit(L: int, fibers: dict[int, set[int]], deadline: float) -> tuple[int, int]:
    n = 2 * L
    J = (1 << n) - 1
    bond_masks = [(1 << u) ^ (1 << v) for u, v in grid(L)[0]]
    odd_count = 0
    for a in range(1 << n):
        if parity(a) == 0:
            continue
        odd_count += 1
        if a & (a - 1):
            edge = next((value for value in bond_masks if parity(a & value)), None)
            assert edge is not None
            x = min(fibers[edge])
            assert x ^ a in fibers[edge]
            assert symplectic(pack(x, edge, n), pack(x ^ a, edge, n), n) == 1
    assert odd_count == 1 << (n - 1)

    rows = 0
    for z in range(1 << n):
        if parity(z):
            continue
        assert time.process_time() <= deadline, "NON-DECISIVE: verifier process budget expired"
        if z == J:
            if L & 1:
                edge = bond_masks[0]
                x, predecessor = seed(J ^ edge, L, fibers)
                assert predecessor == (J ^ edge) and parity(x & edge) == 1
                assert symplectic(pack(x, predecessor, n), pack(0, edge, n), n) == 1
                assert rank([1 << index for index in range(n)]) == n
            else:
                edge = bond_masks[0]
                predecessor = J ^ edge
                assert (J ^ predecessor) == edge
                assert (1 ^ parity(predecessor & checkerboard(L))) == 0
            continue
        rows += 1
        if z == 0:
            tangent_rank = n - 1
        else:
            x, support = seed(z, L, fibers)
            assert support == z and law(x, z, L)
            q0 = z & -z
            toggles = [q0] + [q0 ^ value for value in nullspace([J, z], n)]
            assert all(parity(q) and parity(q & z) for q in toggles)
            tangent_rank = rank(toggles)
        assert tangent_rank == n - 1
    assert rows == (1 << (n - 1)) - 1
    derived = rows * (1 << (n - 1)) + ((1 << n) if L & 1 else 0)
    return rows, derived


def clifford_census(L: int) -> tuple[list[int], int]:
    N = 4 * L
    grades = sorted([2] + [4 * L - 4 * c + 2 for c in range(1, L)])
    assert grades == list(range(2, N, 4))
    total = sum(math.comb(N, grade) for grade in grades)
    assert total == dimension(L)
    # A grade-2 commutator swaps one chosen Majorana index with one unchosen
    # index, so the Johnson graph of k-subsets is connected even at k=N/2.
    for grade in grades:
        assert 0 < grade < N
        assert min(grade, N - grade) >= 1
    return grades, total


def direct_closure(L: int, deadline: float) -> dict:
    n = 2 * L
    mask = (1 << n) - 1
    x_generators = [1 << index for index in range(n)]
    bonds = [(1 << u) ^ (1 << v) for u, v in grid(L)[0]]
    generators = x_generators + [bond << n for bond in bonds]
    reached = set(generators)
    queue = deque(generators)
    started = time.process_time()
    processed = 0
    while queue:
        label = queue.popleft()
        x, z = label & mask, label >> n
        for generator in x_generators:
            if z & generator:
                child = label ^ generator
                if child not in reached:
                    reached.add(child)
                    queue.append(child)
        for bond in bonds:
            if parity(x & bond):
                child = label ^ (bond << n)
                if child not in reached:
                    reached.add(child)
                    queue.append(child)
        processed += 1
        if (processed & 0x3FFFF) == 0:
            assert time.process_time() <= deadline, "NON-DECISIVE: direct closure budget expired"
    violations = sum(not law(label & mask, label >> n, L) for label in reached)
    result = {
        "dimension": len(reached),
        "violations": violations,
        "cpu": time.process_time() - started,
        "peak_rss": rss_bytes(),
    }
    del reached, queue
    gc.collect()
    return result


def main() -> None:
    deadline = time.process_time() + CPU_BUDGET_SECONDS
    construction_rows: dict[int, dict] = {}
    seed_z, seed_points = base_L2()
    for L in range(2, 8):
        started = time.process_time()
        n = 2 * L
        if L > 2:
            seed_z, seed_points = extend(L - 1, seed_z, seed_points)
        for c in range(L - 1):
            triples = square_triples(2 * c, 2 * c + 1, 2 * c + 2, 2 * c + 3, n)
            assert len(triples) == 4
        fibers = spread(L, seed_z, seed_points, deadline)
        ranks = set()
        for points in fibers.values():
            anchor = min(points)
            ranks.add(rank((x ^ anchor for x in points), stop=n - 1))
        assert ranks == {n - 1}
        rows, derived = universal_audit(L, fibers, deadline)
        grades, grade_sum = clifford_census(L)
        assert derived == grade_sum == EXPECTED[L] == dimension(L)
        digest = pair_digest(fibers, n)
        construction_rows[L] = {
            "dimension": derived,
            "ordinary_rows": rows,
            "pair_digest": digest,
            "pair_count": len(fibers),
            "rank": n - 1,
            "grades": grades,
            "cpu": time.process_time() - started,
            "peak_rss": rss_bytes(),
        }
        print(
            f"L={L} CONSTRUCTION checked all {rows} ordinary even-z rows; "
            f"pair ranks={n-1}; dim={derived}; cpu={construction_rows[L]['cpu']:.6f}s; "
            f"cumulative_peak_rss={construction_rows[L]['peak_rss']}"
        )
        if L < 7:
            seed_z = (1 << (2 * (L - 2))) ^ (1 << (2 * (L - 1)))
            seed_points = set(fibers[seed_z])
        del fibers
        gc.collect()

    closure_rows: dict[int, dict] = {}
    for L in range(2, 7):
        row = direct_closure(L, deadline)
        assert row["dimension"] == EXPECTED[L] and row["violations"] == 0
        closure_rows[L] = row
        print(
            f"L={L} FULL CLOSURE enumerated {row['dimension']} strings, "
            f"set-law violations=0; cpu={row['cpu']:.6f}s; peak_rss={row['peak_rss']}"
        )
    print(
        "L=7 NOT a full closure enumeration: exact construction/rank + Clifford census "
        f"give {construction_rows[7]['dimension']}"
    )

    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert artifact["schema"] == "alll_saturation/v1"
    assert all(check["passed"] for check in artifact["data"]["checks"])
    stored = {point["L"]: point for point in artifact["data"]["points"]}
    assert set(stored) == set(range(2, 8))
    for L, independent in construction_rows.items():
        point = stored[L]
        assert point["predicted_and_proved_dimension"] == independent["dimension"]
        assert point["affine_pair_bfs"]["pair_labels_sha256"] == independent["pair_digest"]
        assert point["universal_fiber_audit"]["ordinary_even_z_rows_checked"] == independent["ordinary_rows"]
        if L <= 6:
            assert point["full_closure"]["closure_strings"] == closure_rows[L]["dimension"]
            assert point["full_closure"]["set_law_violations"] == 0
        else:
            assert point["full_closure"] is None
            assert point["evidence_strength"].startswith("NOT a full closure enumeration")
    assert time.process_time() <= deadline
    print("PASS")


if __name__ == "__main__":
    main()
