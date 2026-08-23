#!/usr/bin/env python3
"""Exact local production schemas for the all-L 2xL saturation proof.

This is the bounded local part of the certificate.  It verifies, using only
GF(2) Pauli labels, that edge fibers telescope along a length-two path and
that the two routes across a ladder square produce every three-corner pure-X
string.  No conclusion from a stored artifact is used.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def parity(value: int) -> int:
    return value.bit_count() & 1


def pack(x: int, z: int, n: int) -> int:
    return x | (z << n)


def unpack(label: int, n: int) -> tuple[int, int]:
    mask = (1 << n) - 1
    return label & mask, label >> n


def symplectic(left: int, right: int, n: int) -> int:
    x1, z1 = unpack(left, n)
    x2, z2 = unpack(right, n)
    return parity((x1 & z2) ^ (z1 & x2))


def grid(L: int) -> tuple[list[tuple[int, int]], list[list[int]]]:
    """Return open-ladder bonds and adjacency in column-major labeling."""
    if L < 2:
        raise ValueError("the saturation theorem starts at L=2")
    n = 2 * L
    bonds: list[tuple[int, int]] = []
    for column in range(L):
        bonds.append((2 * column, 2 * column + 1))
    for column in range(L - 1):
        bonds.append((2 * column, 2 * (column + 1)))
        bonds.append((2 * column + 1, 2 * (column + 1) + 1))
    adjacency = [[] for _ in range(n)]
    for u, v in bonds:
        adjacency[u].append(v)
        adjacency[v].append(u)
    for neighbors in adjacency:
        neighbors.sort()
    return bonds, adjacency


def checkerboard_mask(L: int) -> int:
    return sum(
        1 << (2 * column + row)
        for column in range(L)
        for row in (0, 1)
        if ((column + row) & 1) == 0
    )


def edge_label(u: int, v: int, alpha: int, beta: int, n: int) -> int:
    """One of the four labels in the completely generated edge fiber."""
    if alpha not in (0, 1) or beta not in (0, 1):
        raise ValueError("edge endpoint bits must lie in GF(2)")
    x = (alpha << u) ^ (beta << v)
    z = (1 << u) ^ (1 << v)
    return pack(x, z, n)


def path2_witness(
    u: int, v: int, w: int, alpha: int, beta: int, n: int
) -> dict[str, int]:
    """Produce the forced-interior path label on u-v-w from two edge labels."""
    if len({u, v, w}) != 3:
        raise ValueError("a length-two witness needs three distinct vertices")
    left = edge_label(u, v, alpha, 0, n)
    right = edge_label(v, w, 1, beta, n)
    if symplectic(left, right, n) != 1:
        raise AssertionError("length-two parents must anticommute")
    output = left ^ right
    expected = pack(
        (alpha << u) ^ (1 << v) ^ (beta << w),
        (1 << u) ^ (1 << w),
        n,
    )
    if output != expected:
        raise AssertionError("length-two telescope has the wrong label")
    return {"left": left, "right": right, "output": output}


def square_triple_witnesses(
    p: int, q: int, r: int, s: int, n: int
) -> dict[int, dict[str, int]]:
    """Produce all pure-X triples of the square p-r-s-q-p.

    The returned key is the omitted corner.  Every output is obtained from two
    path labels with the same diagonal Z support and odd symplectic pairing.
    """
    if len({p, q, r, s}) != 4:
        raise ValueError("a square needs four distinct vertices")

    # Diagonal {p,s}: interiors r and q.
    ps_r = path2_witness(p, r, s, 0, 0, n)["output"]
    ps_q_with_p = path2_witness(p, q, s, 1, 0, n)["output"]
    ps_q_with_s = path2_witness(p, q, s, 0, 1, n)["output"]

    # Diagonal {r,q}: interiors p and s.
    rq_p = path2_witness(r, p, q, 0, 0, n)["output"]
    rq_s_with_r = path2_witness(r, s, q, 1, 0, n)["output"]
    rq_s_with_q = path2_witness(r, s, q, 0, 1, n)["output"]

    parent_pairs = {
        s: (ps_r, ps_q_with_p),
        p: (ps_r, ps_q_with_s),
        q: (rq_p, rq_s_with_r),
        r: (rq_p, rq_s_with_q),
    }
    all_corners = (1 << p) ^ (1 << q) ^ (1 << r) ^ (1 << s)
    witnesses: dict[int, dict[str, int]] = {}
    for omitted, (left, right) in parent_pairs.items():
        if symplectic(left, right, n) != 1:
            raise AssertionError("the two diagonal-path labels must anticommute")
        output = left ^ right
        expected = pack(all_corners ^ (1 << omitted), 0, n)
        if output != expected:
            raise AssertionError("square move did not produce the requested triple")
        witnesses[omitted] = {"left": left, "right": right, "output": output}
    return witnesses


def audit_local_schemas(L: int, deadline: float) -> dict[str, int | str]:
    n = 2 * L
    bonds, _ = grid(L)
    bond_set = {tuple(sorted(edge)) for edge in bonds}
    witness_count = 0
    digest_accumulator = 0
    for column in range(L - 1):
        if time.process_time() > deadline:
            raise TimeoutError("NON-DECISIVE: e157 process-time budget expired")
        p, q = 2 * column, 2 * column + 1
        r, s = 2 * (column + 1), 2 * (column + 1) + 1
        required = {(p, q), (p, r), (q, s), (r, s)}
        if not all(tuple(sorted(edge)) in bond_set for edge in required):
            raise AssertionError("square boundary is not contained in the ladder")
        witnesses = square_triple_witnesses(p, q, r, s, n)
        if set(witnesses) != {p, q, r, s}:
            raise AssertionError("not every omitted-corner triple was produced")
        witness_count += len(witnesses)
        for row in witnesses.values():
            digest_accumulator ^= row["output"]

    d = checkerboard_mask(L)
    for u, v in bonds:
        if parity(d & ((1 << u) ^ (1 << v))) != 1:
            raise AssertionError("every ladder bond must cross the checkerboard")

    return {
        "L": L,
        "n": n,
        "squares": L - 1,
        "triple_witnesses": witness_count,
        "xor_control": str(digest_accumulator),
    }


def main() -> bool:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu-budget", type=float, default=30.0)
    args = parser.parse_args()
    if args.cpu_budget <= 0:
        raise ValueError("--cpu-budget must be positive")
    deadline = time.process_time() + args.cpu_budget
    rows = [audit_local_schemas(L, deadline) for L in range(2, 8)]
    for row in rows:
        print(
            f"L={row['L']} squares={row['squares']} "
            f"triple_witnesses={row['triple_witnesses']}"
        )
    expected = sum(4 * (L - 1) for L in range(2, 8))
    if sum(int(row["triple_witnesses"]) for row in rows) != expected:
        raise AssertionError("local witness census is incomplete")
    print("PASS e157 local square move")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
