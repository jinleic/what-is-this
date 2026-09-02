#!/usr/bin/env python3
"""Exact off-diagonal census for the remaining public mm3 decompositions.

This campaign instrument extends the corrected ternary-unimodular sandwich census
from paper55/sun56 to mws59/stapleton60.  It has three deliberately separate
entry points:

  --init                 create a producer-named campaign containing ONLY the
                         pre-statement (no target computation)
  --run CAMPAIGN         execute anchors, controls, census, and bound decisions
  --freeze CAMPAIGN      freeze .asrun inputs, manifest, report, and checksums

The standard-matrix convention is fixed throughout.  Each factor row is reshaped
row-major as a 3x3 block M[i][j] = row[3*i+j].  The tensor-preserving action is

  U' = P^-1 U Q^-T,   V' = Q^T V R^-T,   W' = P^T W R.

This is the family for F(A,B,C)=sum_ikj A_ik B_kj C_ij = tr((AB)^T C), not
tr(ABC).  The frozen gatec_sweep.sandwich helper uses transposed-L letters; at
monomials its (X,Y,Z) output equals this action at (P,Q,R)=(X^T,Y^T,Z^T).
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import platform
import random
import shutil
import socket
import sys
import time
import uuid
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

ROOT = Path("/Users/jinleic/jinleic-workspace/cs/mm3")
SRC = ROOT / "src"
SCRATCH = ROOT / "scratch"
CAMPAIGNS = ROOT / "campaigns"
THIS = Path(__file__).resolve()

sys.path.insert(0, str(SRC))
sys.path.insert(0, str(SCRATCH))

N = 3
FLAT = 9
RANK = 23
GAP = RANK - FLAT
SEED = 20260831
TARGET_DECOMPS = ("mws59", "stapleton60")
CONTROL_DECOMPS = ("paper55", "sun56", "mws59", "stapleton60")
ALL_DECOMPS = ("paper55", "perminov58", "sun56", "mws59", "stapleton60")

# Frozen campaign 113838Z values.  These are startup anchors, never campaign
# outputs: every d/floor tuple must reproduce before the new tables are built.
EXPECTED_CLASSES = {
    "paper55": [
        ((12, 13, 13), (False, False, False), 55),
        ((13, 13, 12), (False, False, False), 55),
        ((13, 12, 13), (False, False, False), 55),
    ],
    "sun56": [
        ((12, 11, 16), (False, False, True), 55),
        ((11, 16, 12), (False, True, False), 55),
        ((16, 12, 11), (True, False, False), 55),
    ],
    "mws59": [
        ((13, 12, 14), (False, False, False), 56),
        ((12, 14, 13), (False, False, False), 56),
        ((14, 13, 12), (False, False, False), 56),
    ],
    "stapleton60": [
        ((14, 14, 13), (False, False, False), 58),
        ((14, 13, 14), (False, False, False), 58),
        ((13, 14, 14), (False, False, False), 58),
    ],
}

PHASE_A_CAP_S = 2 * 3600
TABLE_CAP_S = 8 * 3600
DECISION_CAP_S = 24 * 3600
CAMPAIGN_CAP_S = 48 * 3600


class CampaignError(RuntimeError):
    pass


class RecordFlag(CampaignError):
    """Raised before an apparent <=54 row is written to any artifact."""


class Logger:
    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("a", encoding="utf-8", buffering=1)

    def __call__(self, message: str) -> None:
        stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        line = f"[{stamp}] {message}"
        print(line, flush=True)
        self.handle.write(line + "\n")

    def close(self) -> None:
        self.handle.close()


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def mat_tuple(rows: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(x) for x in row) for row in rows)


def flat_mat(rows: Sequence[Sequence[int]]) -> tuple[int, ...]:
    return tuple(int(x) for row in rows for x in row)


def block(row: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    if len(row) != FLAT:
        raise ValueError(f"expected 9 entries, got {len(row)}")
    return tuple(tuple(int(row[3 * i + j]) for j in range(3)) for i in range(3))


def flat_block(rows: Sequence[Sequence[int]]) -> tuple[int, ...]:
    return tuple(int(rows[i][j]) for i in range(3) for j in range(3))


def transpose(a: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(a[j][i]) for j in range(3)) for i in range(3))


def matmul(a: Sequence[Sequence[int]], b: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(sum(int(a[i][k]) * int(b[k][j]) for k in range(3)) for j in range(3))
        for i in range(3)
    )


def det3(a: Sequence[Sequence[int]]) -> int:
    return (
        int(a[0][0]) * (int(a[1][1]) * int(a[2][2]) - int(a[1][2]) * int(a[2][1]))
        - int(a[0][1]) * (int(a[1][0]) * int(a[2][2]) - int(a[1][2]) * int(a[2][0]))
        + int(a[0][2]) * (int(a[1][0]) * int(a[2][1]) - int(a[1][1]) * int(a[2][0]))
    )


def inverse_unimodular(a: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    d = det3(a)
    if d not in (-1, 1):
        raise ValueError(f"matrix is not unimodular, determinant={d}")
    cof = (
        (
            int(a[1][1]) * int(a[2][2]) - int(a[1][2]) * int(a[2][1]),
            -(int(a[1][0]) * int(a[2][2]) - int(a[1][2]) * int(a[2][0])),
            int(a[1][0]) * int(a[2][1]) - int(a[1][1]) * int(a[2][0]),
        ),
        (
            -(int(a[0][1]) * int(a[2][2]) - int(a[0][2]) * int(a[2][1])),
            int(a[0][0]) * int(a[2][2]) - int(a[0][2]) * int(a[2][0]),
            -(int(a[0][0]) * int(a[2][1]) - int(a[0][1]) * int(a[2][0])),
        ),
        (
            int(a[0][1]) * int(a[1][2]) - int(a[0][2]) * int(a[1][1]),
            -(int(a[0][0]) * int(a[1][2]) - int(a[0][2]) * int(a[1][0])),
            int(a[0][0]) * int(a[1][1]) - int(a[0][1]) * int(a[1][0]),
        ),
    )
    inv = tuple(tuple(transpose(cof)[i][j] // d for j in range(3)) for i in range(3))
    ident = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    if matmul(a, inv) != ident or matmul(inv, a) != ident:
        raise AssertionError("unimodular inverse audit failed")
    return inv


def is_monomial(a: Sequence[Sequence[int]]) -> bool:
    return (
        all(sum(int(x) != 0 for x in row) == 1 for row in a)
        and all(sum(int(a[i][j]) != 0 for i in range(3)) == 1 for j in range(3))
        and all(int(x) in (-1, 0, 1) for row in a for x in row)
    )


def signed_monomials() -> tuple[tuple[tuple[int, ...], ...], ...]:
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            a = [[0] * 3 for _ in range(3)]
            for i, j in enumerate(perm):
                a[i][j] = signs[i]
            out.append(mat_tuple(a))
    out = sorted(set(out), key=flat_mat)
    if len(out) != 48 or any(not is_monomial(a) for a in out):
        raise AssertionError("signed-monomial construction failed")
    return tuple(out)


def enumerate_ternary_unimodular() -> tuple[list[tuple[tuple[int, ...], ...]], dict[str, int]]:
    unimodular = []
    invertible = 0
    inverse_ternary = 0
    inverse_ternary_nonmonomial = 0
    for values in itertools.product((-1, 0, 1), repeat=9):
        a = block(values)
        d = det3(a)
        if d:
            invertible += 1
        if d in (-1, 1):
            unimodular.append(a)
            inv = inverse_unimodular(a)
            if all(x in (-1, 0, 1) for row in inv for x in row):
                inverse_ternary += 1
                if not is_monomial(a):
                    inverse_ternary_nonmonomial += 1
    monomial = sum(is_monomial(a) for a in unimodular)
    stats = {
        "ternary_matrices": 3**9,
        "invertible": invertible,
        "unimodular": len(unimodular),
        "monomial": monomial,
        "nonmonomial_unimodular": len(unimodular) - monomial,
        "inverse_ternary": inverse_ternary,
        "inverse_ternary_nonmonomial": inverse_ternary_nonmonomial,
    }
    expected = {
        "ternary_matrices": 19683,
        "invertible": 11808,
        "unimodular": 6960,
        "monomial": 48,
        "nonmonomial_unimodular": 6912,
        "inverse_ternary": 4656,
        "inverse_ternary_nonmonomial": 4608,
    }
    if stats != expected:
        raise AssertionError(f"ternary census anchor mismatch: {stats} != {expected}")
    return unimodular, stats


def build_data_orbits(
    unimodular: Sequence[tuple[tuple[int, ...], ...]],
    monomials: Sequence[tuple[tuple[int, ...], ...]],
) -> tuple[list[tuple[tuple[int, ...], ...]], dict[tuple[int, ...], int], dict[str, object]]:
    universe = {flat_mat(a): a for a in unimodular}
    unassigned = set(universe)
    raw_orbits = []
    while unassigned:
        seed = universe[min(unassigned)]
        orbit = {flat_mat(matmul(seed, m)) for m in monomials}
        if len(orbit) != 48 or not orbit <= universe.keys():
            raise AssertionError("right-monomial action did not give a 48-element T-orbit")
        raw_orbits.append(tuple(sorted(orbit)))
        unassigned -= orbit
    raw_orbits.sort(key=lambda orb: orb[0])
    reps = [universe[orb[0]] for orb in raw_orbits]
    member_to_rep: dict[tuple[int, ...], int] = {}
    factorization_uniqueness = True
    for idx, orb in enumerate(raw_orbits):
        rep = reps[idx]
        for member in orb:
            member_to_rep[member] = idx
            matches = sum(flat_mat(matmul(rep, m)) == member for m in monomials)
            factorization_uniqueness &= matches == 1
    ident = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    mono_idx = member_to_rep[ident]
    if len(reps) != 145 or len(member_to_rep) != 6960 or not factorization_uniqueness:
        raise AssertionError("T=A*M factorization failed")
    if not is_monomial(reps[mono_idx]):
        raise AssertionError("identity orbit representative is not monomial")
    summary = {
        "data_nodes": len(reps),
        "right_orbit_size": 48,
        "partition_size": len(member_to_rep),
        "unique_factorization": factorization_uniqueness,
        "monomial_data_node": mono_idx,
        "representative_rule": "lexicographically least flattened matrix in each right-monomial orbit",
    }
    return reps, member_to_rep, summary


def map_blocks(
    blocks: Sequence[Sequence[int]],
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
) -> list[tuple[int, ...]]:
    return [flat_block(matmul(matmul(left, block(row)), right)) for row in blocks]


def honest_sandwich(
    U: Sequence[Sequence[int]],
    V: Sequence[Sequence[int]],
    W: Sequence[Sequence[int]],
    P: Sequence[Sequence[int]],
    Q: Sequence[Sequence[int]],
    R: Sequence[Sequence[int]],
) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], list[tuple[int, ...]]]:
    """Standard semantics: U'=P^-1 U Q^-T, V'=Q^T V R^-T, W'=P^T W R."""
    Pi = inverse_unimodular(P)
    Qi = inverse_unimodular(Q)
    Ri = inverse_unimodular(R)
    U2 = map_blocks(U, Pi, transpose(Qi))
    V2 = map_blocks(V, transpose(Q), transpose(Ri))
    W2 = map_blocks(W, transpose(P), R)
    return U2, V2, W2


def wrong_w_variant(
    U: Sequence[Sequence[int]],
    V: Sequence[Sequence[int]],
    W: Sequence[Sequence[int]],
    G: Sequence[Sequence[int]],
) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], list[tuple[int, ...]]]:
    """Counterfactual plant: honest U,V but wrong standard W'=G W G^T."""
    U2, V2, _ = honest_sandwich(U, V, W, G, G, G)
    W_bad = map_blocks(W, G, transpose(G))
    return U2, V2, W_bad


def transpose_blocks(blocks: Sequence[Sequence[int]]) -> list[tuple[int, ...]]:
    return [flat_block(transpose(block(row))) for row in blocks]


def sigma_orbit(
    U: Sequence[Sequence[int]], V: Sequence[Sequence[int]], W: Sequence[Sequence[int]]
) -> list[tuple[list[tuple[int, ...]], list[tuple[int, ...]], list[tuple[int, ...]]]]:
    U0, V0, W0 = [list(map(tuple, x)) for x in (U, V, W)]
    return [
        (U0, V0, W0),
        (V0, transpose_blocks(W0), transpose_blocks(U0)),
        (transpose_blocks(W0), U0, transpose_blocks(V0)),
    ]


def all_ternary(*sides: Sequence[Sequence[int]]) -> bool:
    return all(int(x) in (-1, 0, 1) for side in sides for row in side for x in row)


def brent_failures_int(
    U: Sequence[Sequence[int]], V: Sequence[Sequence[int]], W: Sequence[Sequence[int]]
) -> int:
    failures = 0
    for i, j, k, ip, jp, kp in itertools.product(range(3), repeat=6):
        got = sum(
            int(U[r][3 * i + k]) * int(V[r][3 * kp + j]) * int(W[r][3 * ip + jp])
            for r in range(RANK)
        )
        expected = int(i == ip and j == jp and k == kp)
        failures += got != expected
    return failures


def brent_failures_fmpz(
    U: Sequence[Sequence[int]], V: Sequence[Sequence[int]], W: Sequence[Sequence[int]]
) -> int:
    from flint import fmpz

    failures = 0
    for i, j, k, ip, jp, kp in itertools.product(range(3), repeat=6):
        got = fmpz(0)
        for r in range(RANK):
            got += fmpz(U[r][3 * i + k]) * fmpz(V[r][3 * kp + j]) * fmpz(W[r][3 * ip + jp])
        expected = int(i == ip and j == jp and k == kp)
        failures += int(got) != expected
    return failures


def canon_vec(v: Sequence[int]) -> tuple[int, ...]:
    t = tuple(int(x) for x in v)
    neg = tuple(-x for x in t)
    return t if t <= neg else neg


def input_dirs() -> tuple[tuple[int, ...], ...]:
    return tuple(canon_vec(tuple(int(i == j) for j in range(FLAT))) for i in range(FLAT))


def prep_floor(targets: Sequence[Sequence[int]]) -> tuple[
    list[tuple[int, ...]], list[tuple[int, ...]], list[list[tuple[int, int]]]
]:
    inputs = set(input_dirs())
    classes = sorted({canon_vec(t) for t in targets if any(int(x) for x in t)} - inputs)
    base = sorted(set(classes) | inputs)
    reps: list[list[tuple[int, int]]] = []
    for c in classes:
        pairs = []
        for ia, a in enumerate(base):
            for ib in range(ia, len(base)):
                b = base[ib]
                found = False
                for sa, sb in itertools.product((-1, 1), repeat=2):
                    s = tuple(sa * a[i] + sb * b[i] for i in range(FLAT))
                    if canon_vec(s) == c:
                        found = True
                        break
                if found:
                    pairs.append((ia, ib))
        reps.append(pairs)
    return classes, base, reps


def subset_floor_decide(
    classes: Sequence[tuple[int, ...]],
    base: Sequence[tuple[int, ...]],
    reps: Sequence[Sequence[tuple[int, int]]],
) -> tuple[bool, int, list[int]]:
    c_index = {c: i for i, c in enumerate(classes)}
    b_index = {c: i for i, c in enumerate(base)}
    inputs = set(input_dirs())
    init_av = 0
    for c in inputs:
        init_av |= 1 << b_index[c]
    class_base_bits = [1 << b_index[c] for c in classes]
    full = (1 << len(classes)) - 1
    states = 0
    choice: dict[int, int] = {}

    @lru_cache(maxsize=None)
    def can(mask: int) -> bool:
        nonlocal states
        states += 1
        if mask == full:
            return True
        av = init_av
        for ci, bit in enumerate(class_base_bits):
            if mask & (1 << ci):
                av |= bit
        for ci, pairs in enumerate(reps):
            bit = 1 << ci
            if mask & bit:
                continue
            if any((av & (1 << ia)) and (av & (1 << ib)) for ia, ib in pairs):
                if can(mask | bit):
                    choice[mask] = ci
                    return True
        return False

    ok = can(0)
    order = []
    if ok:
        mask = 0
        while mask != full:
            ci = choice[mask]
            order.append(ci)
            mask |= 1 << ci
        if len(order) != len(classes) or set(order) != set(range(len(classes))):
            raise AssertionError("floor witness order is not a class bijection")
        available = set(inputs)
        for ci in order:
            c = classes[ci]
            if c in available:
                raise AssertionError("floor witness repeats an available target")
            if not any(base[ia] in available and base[ib] in available for ia, ib in reps[ci]):
                raise AssertionError("floor witness transition is unavailable")
            available.add(c)
        if not set(classes) <= available:
            raise AssertionError("floor witness does not cover all classes")
    return ok, states, order


def side_decision(targets: Sequence[Sequence[int]], include_witness: bool = False) -> dict[str, object]:
    classes, base, reps = prep_floor(targets)
    ok, states, order = subset_floor_decide(classes, base, reps)
    out: dict[str, object] = {
        "d": len(classes),
        "floor_exists": ok,
        "lb": len(classes) + int(not ok),
        "states": states,
    }
    if include_witness and ok:
        out["order"] = [list(classes[i]) for i in order]
    return out


def certified_total(
    U: Sequence[Sequence[int]], V: Sequence[Sequence[int]], W: Sequence[Sequence[int]]
) -> dict[str, object]:
    sides = [side_decision(side) for side in (U, V, W)]
    return {
        "sides": sides,
        "total_lb": sum(int(x["lb"]) for x in sides) + GAP,
    }


def rank_q(rows: Sequence[Sequence[int]]) -> int:
    a = [[Fraction(int(x)) for x in row] for row in rows]
    rank = 0
    col = 0
    while rank < len(a) and col < FLAT:
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            col += 1
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        p = a[rank][col]
        a[rank] = [x / p for x in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][col]:
                q = a[i][col]
                a[i] = [a[i][j] - q * a[rank][j] for j in range(FLAT)]
        rank += 1
        col += 1
    return rank


def activity_audit(U: Sequence[Sequence[int]], V: Sequence[Sequence[int]], W: Sequence[Sequence[int]]) -> dict[str, object]:
    sides = []
    for side in (U, V, W):
        sides.append({
            "nonzero_rows": sum(any(int(x) for x in row) for row in side),
            "rank_Q": rank_q(side),
            "nonzero_columns": sum(any(int(side[r][j]) for r in range(RANK)) for j in range(FLAT)),
        })
    ok = all(x == {"nonzero_rows": 23, "rank_Q": 9, "nonzero_columns": 9} for x in sides)
    if not ok:
        raise AssertionError(f"transposition activity precondition failed: {sides}")
    return {"sides": sides, "ok": True}


def load_decompositions() -> tuple[dict[str, tuple[list[tuple[int, ...]], ...]], dict[str, object]]:
    from gatec_decomps import LOADERS, META

    data = {}
    for name in ALL_DECOMPS:
        U, V, W = LOADERS[name]()
        data[name] = (list(map(tuple, U)), list(map(tuple, V)), list(map(tuple, W)))
    return data, META


def anchor_all(
    data: dict[str, tuple[list[tuple[int, ...]], ...]],
    meta: dict[str, object],
    monomial_rep: Sequence[Sequence[int]],
) -> dict[str, object]:
    import gate_b_floor as frozen_floor
    import gatec_sweep as frozen_sweep
    from verify_anchors import brent_failures as frozen_brent

    published = {}
    for name in ALL_DECOMPS:
        U, V, W = data[name]
        got_int = brent_failures_int(U, V, W)
        got_fmpz = brent_failures_fmpz(U, V, W)
        got_frozen = frozen_brent(U, V, W)
        ternary = all_ternary(U, V, W)
        anchor = int(meta[name]["anchor"])
        split_meta = meta[name]["split"]
        count_provenance: dict[str, object]
        if split_meta is None:
            if name != "perminov58":
                raise AssertionError(f"unexpected missing split for {name}")
            import gate_a

            U9 = [[int(U[r][i]) for r in range(RANK)] for i in range(FLAT)]
            V9 = [[int(V[r][i]) for r in range(RANK)] for i in range(FLAT)]
            W9 = [[int(W[r][i]) for r in range(RANK)] for i in range(FLAT)]
            mismatches, complexity = gate_a.check_perminov(
                SCRATCH / "cr58_cn122_ZT_reduced.json", U9, V9, W9
            )
            if mismatches or int(complexity["reduced"]) != anchor:
                raise AssertionError(
                    f"Perminov JSON anchor failed: mismatches={len(mismatches)} "
                    f"complexity={complexity}"
                )
            published_split = None
            count_provenance = {
                "kind": "Perminov pinned JSON complexity field",
                "complexity": complexity,
                "decoded_factor_mismatches": len(mismatches),
            }
        else:
            published_split = tuple(int(x) for x in split_meta)
            if sum(published_split) != anchor:
                raise AssertionError(f"published split anchor failed for {name}")
            count_provenance = {"kind": "recounted left/right/output split"}
        if (got_int, got_fmpz, got_frozen) != (0, 0, 0) or not ternary:
            raise AssertionError(f"published Brent/ternary anchor failed for {name}")
        published[name] = {
            "brent_int": 729 - got_int,
            "brent_fmpz": 729 - got_fmpz,
            "brent_frozen": 729 - got_frozen,
            "ternary": ternary,
            "published_split": list(published_split) if published_split is not None else None,
            "published_total": anchor,
            "count_provenance": count_provenance,
        }

    certified = {}
    for name in CONTROL_DECOMPS:
        U, V, W = data[name]
        rows = []
        for sigma, triple in enumerate(sigma_orbit(U, V, W)):
            dec = certified_total(*triple)
            expected_d, expected_floor, expected_total = EXPECTED_CLASSES[name][sigma]
            got_d = tuple(int(s["d"]) for s in dec["sides"])
            got_floor = tuple(bool(s["floor_exists"]) for s in dec["sides"])
            if got_d != expected_d or got_floor != expected_floor or dec["total_lb"] != expected_total:
                raise AssertionError(
                    f"certified anchor mismatch {name} sigma={sigma}: "
                    f"{got_d}/{got_floor}/{dec['total_lb']}"
                )
            frozen = []
            for side in triple:
                classes, reps = frozen_floor.prep(list(side))
                ok, _stats, _order = frozen_floor.subset_dfs(classes, reps)
                frozen.append((len(classes), bool(ok)))
            if tuple(x[0] for x in frozen) != got_d or tuple(x[1] for x in frozen) != got_floor:
                raise AssertionError(f"independent frozen floor mismatch {name} sigma={sigma}")
            rows.append({
                "sigma": sigma,
                "d": list(got_d),
                "floor_exists": list(got_floor),
                "total_lb": int(dec["total_lb"]),
                "states_own": [int(s["states"]) for s in dec["sides"]],
            })

        oriented = honest_sandwich(U, V, W, monomial_rep, monomial_rep, monomial_rep)
        mono_total = int(certified_total(*oriented)["total_lb"])
        expected_total = EXPECTED_CLASSES[name][0][2]
        if mono_total != expected_total:
            raise AssertionError(f"all-monomial data triple anchor failed for {name}: {mono_total}")
        certified[name] = {"sigma_classes": rows, "all_monomial_data_total_lb": mono_total}

    # Mandatory startup regression: the all-monomial data triple must be 55 in
    # both completed 55-lower-bound landscapes.  This is intentionally separate
    # from the new 56/58 anchors: it catches a missing +14 transposition gap.
    if certified["paper55"]["all_monomial_data_total_lb"] != 55:
        raise AssertionError("MANDATORY paper55 all-monomial anchor is not 55")
    if certified["sun56"]["all_monomial_data_total_lb"] != 55:
        raise AssertionError("MANDATORY sun56 all-monomial anchor is not 55")
    if certified["mws59"]["all_monomial_data_total_lb"] != 56:
        raise AssertionError("mws59 all-monomial anchor is not 56")
    if certified["stapleton60"]["all_monomial_data_total_lb"] != 58:
        raise AssertionError("stapleton60 all-monomial anchor is not 58")

    # Letter-convention control: frozen (X,Y,Z) equals honest
    # (P,Q,R)=(X^T,Y^T,Z^T) on signed monomials, bit-for-bit.
    monos = signed_monomials()
    rng = random.Random(SEED)
    convention_checks = 0
    for name in ALL_DECOMPS:
        U, V, W = data[name]
        for _ in range(5):
            X, Y, Z = (monos[rng.randrange(48)] for _ in range(3))
            frozen = frozen_sweep.sandwich(U, V, W, X, Y, Z)
            honest = honest_sandwich(U, V, W, transpose(X), transpose(Y), transpose(Z))
            if tuple(tuple(map(tuple, side)) for side in frozen) != tuple(tuple(map(tuple, side)) for side in honest):
                raise AssertionError(f"frozen/honest monomial convention mismatch for {name}")
            convention_checks += 1

    return {
        "published": published,
        "certified": certified,
        "mandatory_55_anchor_asserted_before_new_counts": True,
        "frozen_vs_standard_letter_controls": {
            "checks": convention_checks,
            "rule": "frozen (X,Y,Z) == standard honest (P,Q,R)=(X^T,Y^T,Z^T) on monomials",
        },
    }


def semantic_controls(data: dict[str, tuple[list[tuple[int, ...]], ...]]) -> dict[str, object]:
    import gatec_sweep as frozen_sweep

    # Fixed witness from frozen session-7 artifact witness_sun56_nonmono_diagonal.json.
    G = ((-1, -1, -1), (0, 1, 1), (0, 0, 1))
    if det3(G) != -1 or is_monomial(G):
        raise AssertionError("fixed nonmonomial control matrix changed")
    U, V, W = data["sun56"]
    honest = honest_sandwich(U, V, W, G, G, G)
    honest_int = brent_failures_int(*honest)
    honest_fmpz = brent_failures_fmpz(*honest)
    honest_ternary = all_ternary(*honest)
    if (honest_int, honest_fmpz, honest_ternary) != (0, 0, True):
        raise AssertionError("fixed final-family witness no longer passes")

    swapped = wrong_w_variant(U, V, W, G)
    swapped_fail = brent_failures_int(*swapped)
    swapped_max = max(abs(int(x)) for side in swapped for row in side for x in row)
    if swapped_fail == 0 and all_ternary(*swapped):
        raise AssertionError("wrong-W counterfactual plant was not detected")

    frozen_bad = frozen_sweep.sandwich(U, V, W, G, G, G)
    frozen_bad_fail = brent_failures_int(*frozen_bad)
    if frozen_bad_fail == 0:
        raise AssertionError("transpose-as-inverse counterfactual plant was not detected")

    # Direct Brent corruption plant: alter a nonzero W coefficient of an active summand.
    W_corrupt = [list(row) for row in W]
    changed = None
    for r, row in enumerate(W_corrupt):
        for j, x in enumerate(row):
            if x and any(U[r]) and any(V[r]):
                row[j] = -x
                changed = (r, j, x, -x)
                break
        if changed:
            break
    corrupt_fail = brent_failures_int(U, V, W_corrupt)
    if not changed or corrupt_fail == 0:
        raise AssertionError("Brent corruption plant was not detected")

    nonunimodular = ((2, 0, 0), (0, 1, 0), (0, 0, 1))
    rejected_nonunimodular = False
    try:
        inverse_unimodular(nonunimodular)
    except ValueError:
        rejected_nonunimodular = True
    if not rejected_nonunimodular:
        raise AssertionError("non-unimodular counterfactual plant entered the domain")

    # Exact-rank/activity checks justify the +14 direction for every orientation:
    # each honest block action is an invertible 9-coordinate map.
    activity = {}
    for name in CONTROL_DECOMPS:
        activity[name] = []
        for sigma, triple in enumerate(sigma_orbit(*data[name])):
            activity[name].append({"sigma": sigma, **activity_audit(*triple)})

    return {
        "fixed_G": [list(row) for row in G],
        "honest_final_family": {
            "brent_int": 729 - honest_int,
            "brent_fmpz": 729 - honest_fmpz,
            "ternary": honest_ternary,
        },
        "wrong_W_plant": {"brent_failures": swapped_fail, "max_abs_entry": swapped_max},
        "transpose_as_inverse_plant": {"brent_failures": frozen_bad_fail},
        "direct_brent_corruption_plant": {"changed": list(changed), "brent_failures": corrupt_fail},
        "nonunimodular_det2_rejected": rejected_nonunimodular,
        "activity": activity,
    }


def numpy_pair_targets(
    blocks: np.ndarray,
    side: str,
    a: int,
    b: int,
    reps: np.ndarray,
    inverses: np.ndarray,
) -> np.ndarray:
    if side == "U":       # pair is (P,Q)
        return inverses[a] @ blocks @ inverses[b].T
    if side == "V":       # pair is (Q,R)
        return reps[a].T @ blocks @ inverses[b].T
    if side == "W":       # pair is (P,R)
        return reps[a].T @ blocks @ reps[b]
    raise ValueError(side)


def python_pair_targets(
    blocks: Sequence[Sequence[int]],
    side: str,
    a: Sequence[Sequence[int]],
    b: Sequence[Sequence[int]],
) -> list[tuple[int, ...]]:
    if side == "U":
        return map_blocks(blocks, inverse_unimodular(a), transpose(inverse_unimodular(b)))
    if side == "V":
        return map_blocks(blocks, transpose(a), transpose(inverse_unimodular(b)))
    if side == "W":
        return map_blocks(blocks, transpose(a), b)
    raise ValueError(side)


def build_pair_tables(
    data: tuple[list[tuple[int, ...]], ...],
    reps_list: Sequence[tuple[tuple[int, ...], ...]],
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    U, V, W = data
    reps = np.asarray(reps_list, dtype=np.int64)
    inverses = np.asarray([inverse_unimodular(a) for a in reps_list], dtype=np.int64)
    blocks = {
        "U": np.asarray([block(row) for row in U], dtype=np.int64),
        "V": np.asarray([block(row) for row in V], dtype=np.int64),
        "W": np.asarray([block(row) for row in W], dtype=np.int64),
    }
    tables = {side: np.zeros((len(reps), len(reps)), dtype=np.bool_) for side in ("U", "V", "W")}
    max_abs = {side: 0 for side in tables}
    for side in ("U", "V", "W"):
        for a in range(len(reps)):
            for b in range(len(reps)):
                targets = numpy_pair_targets(blocks[side], side, a, b, reps, inverses)
                max_abs[side] = max(max_abs[side], int(np.max(np.abs(targets))))
                tables[side][a, b] = bool(np.all((targets >= -1) & (targets <= 1)))
    # int64 bound: inverse entries <=2, so U<=36, V<=18, W<=9.  These
    # maxima are checked, not merely stated; overflow is impossible.
    hard_bounds = {"U": 36, "V": 18, "W": 9}
    if any(max_abs[s] > hard_bounds[s] for s in hard_bounds):
        raise AssertionError(f"exact int64 bound derivation failed: {max_abs}")
    info = {
        "positive_pairs": {side: int(np.count_nonzero(tables[side])) for side in tables},
        "max_abs_observed": max_abs,
        "int64_hard_bounds": hard_bounds,
        "int64_overflow_impossible": True,
    }
    return tables, info


def direct_table_controls(
    data: tuple[list[tuple[int, ...]], ...],
    reps: Sequence[tuple[tuple[int, ...], ...]],
    tables: dict[str, np.ndarray],
    mono_idx: int,
    seed: int,
) -> dict[str, object]:
    U, V, W = data
    mismatches = 0
    direct_positives = 0
    table_positives = 0
    p = mono_idx
    for q in range(len(reps)):
        for r in range(len(reps)):
            direct = all_ternary(*honest_sandwich(U, V, W, reps[p], reps[q], reps[r]))
            table = bool(tables["U"][p, q] and tables["V"][q, r] and tables["W"][p, r])
            mismatches += direct != table
            direct_positives += direct
            table_positives += table
    if mismatches:
        raise AssertionError(f"full 145^2 direct/table subcube mismatch count={mismatches}")

    rng = random.Random(seed)
    random_mismatches = 0
    for _ in range(200):
        p, q, r = (rng.randrange(len(reps)) for _ in range(3))
        direct = all_ternary(*honest_sandwich(U, V, W, reps[p], reps[q], reps[r]))
        table = bool(tables["U"][p, q] and tables["V"][q, r] and tables["W"][p, r])
        random_mismatches += direct != table
    if random_mismatches:
        raise AssertionError(f"random direct/table mismatch count={random_mismatches}")

    # Numpy vs pure-Python transform checks all entries, not just predicates.
    np_reps = np.asarray(reps, dtype=np.int64)
    np_inv = np.asarray([inverse_unimodular(a) for a in reps], dtype=np.int64)
    base = {"U": np.asarray([block(x) for x in U], dtype=np.int64),
            "V": np.asarray([block(x) for x in V], dtype=np.int64),
            "W": np.asarray([block(x) for x in W], dtype=np.int64)}
    entry_checks = 0
    for _ in range(24):
        side = ("U", "V", "W")[_ % 3]
        a, b = rng.randrange(len(reps)), rng.randrange(len(reps))
        got_np = numpy_pair_targets(base[side], side, a, b, np_reps, np_inv)
        got_py = python_pair_targets((U, V, W)[("U", "V", "W").index(side)], side, reps[a], reps[b])
        if got_np.tolist() != [[list(row) for row in block(x)] for x in got_py]:
            raise AssertionError("numpy/pure-Python pair-transform mismatch")
        entry_checks += 1

    # Structural fiber control: right-monomial factors must not change the
    # ternarity predicate.  This probes the exact P=a*p factorization, not only
    # the data representatives used by the table builder.
    monomials = signed_monomials()
    fiber_mismatches = 0
    for _ in range(25):
        p, q, r = (rng.randrange(len(reps)) for _ in range(3))
        mp, mq, mr = (monomials[rng.randrange(48)] for _ in range(3))
        P, Q, R = matmul(reps[p], mp), matmul(reps[q], mq), matmul(reps[r], mr)
        data_predicate = all_ternary(
            *honest_sandwich(U, V, W, reps[p], reps[q], reps[r])
        )
        fiber_predicate = all_ternary(*honest_sandwich(U, V, W, P, Q, R))
        fiber_mismatches += data_predicate != fiber_predicate
    if fiber_mismatches:
        raise AssertionError(f"right-monomial fiber predicate mismatches={fiber_mismatches}")

    return {
        "anchored_subcube": {
            "p_index": mono_idx,
            "triples": len(reps) ** 2,
            "mismatches": mismatches,
            "direct_positives": direct_positives,
            "table_positives": table_positives,
        },
        "random_concordance": {"seed": seed, "trials": 200, "mismatches": random_mismatches},
        "numpy_vs_python_entry_checks": entry_checks,
        "right_monomial_fiber_predicate": {
            "trials": 25,
            "mismatches": fiber_mismatches,
        },
    }


def enumerate_survivors(tables: dict[str, np.ndarray], mono_idx: int) -> tuple[np.ndarray, dict[str, object]]:
    U, V, W = tables["U"], tables["V"], tables["W"]
    triples = []
    for p in range(U.shape[0]):
        for q in np.flatnonzero(U[p]).tolist():
            for r in np.flatnonzero(V[q] & W[p]).tolist():
                triples.append((p, q, r))
    arr = np.asarray(triples, dtype=np.uint16).reshape((-1, 3))
    path1 = len(triples)
    dense = U.astype(np.int64) @ V.astype(np.int64)
    path2 = int(np.sum(dense * W.astype(np.int64), dtype=np.int64))
    if path1 != path2:
        raise AssertionError(f"independent survivor count mismatch {path1} != {path2}")
    mono_counts = Counter(sum(int(x == mono_idx) for x in t) for t in triples)
    return arr, {
        "survivor_data_triples": path1,
        "independent_dense_count": path2,
        "by_monomial_node_count": {str(k): int(mono_counts.get(k, 0)) for k in range(4)},
        "lexicographic_order": True,
    }


def decide_pairs_and_triples(
    name: str,
    data: tuple[list[tuple[int, ...]], ...],
    reps_list: Sequence[tuple[tuple[int, ...], ...]],
    tables: dict[str, np.ndarray],
    survivors: np.ndarray,
    mono_idx: int,
    log: Logger,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    U, V, W = data
    base = {"U": U, "V": V, "W": W}
    n = len(reps_list)
    decisions: dict[str, np.ndarray] = {}
    state_tables: dict[str, np.ndarray] = {}
    d_tables: dict[str, np.ndarray] = {}
    floor_tables: dict[str, np.ndarray] = {}
    lb_tables: dict[str, np.ndarray] = {}

    started = time.perf_counter()
    for side in ("U", "V", "W"):
        d_tab = np.full((n, n), -1, dtype=np.int16)
        floor_tab = np.full((n, n), -1, dtype=np.int8)
        lb_tab = np.full((n, n), -1, dtype=np.int16)
        states_tab = np.full((n, n), -1, dtype=np.int64)
        positives = np.argwhere(tables[side])
        for offset, (a_raw, b_raw) in enumerate(positives):
            a, b = int(a_raw), int(b_raw)
            targets = python_pair_targets(base[side], side, reps_list[a], reps_list[b])
            if not all_ternary(targets):
                raise AssertionError("positive pair table produced nonternary targets")
            dec = side_decision(targets)
            d_tab[a, b] = int(dec["d"])
            floor_tab[a, b] = int(bool(dec["floor_exists"]))
            lb_tab[a, b] = int(dec["lb"])
            states_tab[a, b] = int(dec["states"])
            if offset and offset % 256 == 0:
                log(f"{name} {side}: decided {offset}/{len(positives)} ternary pairs")
            if time.perf_counter() - started > DECISION_CAP_S:
                raise CampaignError(f"decision cap hit in {name} {side}; no partial claim")
        d_tables[side] = d_tab
        floor_tables[side] = floor_tab
        lb_tables[side] = lb_tab
        state_tables[side] = states_tab

    totals = np.empty(len(survivors), dtype=np.int16)
    side_lbs = np.empty((len(survivors), 3), dtype=np.int16)
    for idx, (p_raw, q_raw, r_raw) in enumerate(survivors):
        p, q, r = int(p_raw), int(q_raw), int(r_raw)
        vals = (lb_tables["U"][p, q], lb_tables["V"][q, r], lb_tables["W"][p, r])
        if min(vals) < 0:
            raise AssertionError("survivor references undecided pair")
        side_lbs[idx] = vals
        totals[idx] = int(sum(int(x) for x in vals) + GAP)

    # Mandatory bookkeeping anchor is asserted before any histogram or minimum
    # is computed or written.  The monomial data triple is guaranteed present.
    mono_rows = np.flatnonzero(np.all(survivors == mono_idx, axis=1))
    if len(mono_rows) != 1:
        raise AssertionError(f"{name}: monomial data triple multiplicity {len(mono_rows)}")
    mono_total = int(totals[int(mono_rows[0])])
    expected = 56 if name == "mws59" else 58
    if mono_total != expected:
        raise AssertionError(f"{name}: mandatory all-monomial total {mono_total}, expected {expected}")

    # The user requires escalation before recording any <=54.  Deliberately do
    # not write the offending row, orientation, or histogram before raising.
    if np.any(totals <= 54):
        raise RecordFlag("A <=54 bound-layer row was detected; notify Main before recording details")

    # Independent frozen d-count/floor cross-check on deterministic positive pairs.
    import gate_b_floor as frozen_floor
    import gatec_sweep as frozen_sweep

    rng = random.Random(SEED + (0 if name == "mws59" else 1000))
    positive_keys = []
    for side in ("U", "V", "W"):
        for a_raw, b_raw in np.argwhere(tables[side]).tolist():
            positive_keys.append((side, int(a_raw), int(b_raw)))
    rng.shuffle(positive_keys)
    cross = []
    for side, a, b in positive_keys[:18]:
        targets = python_pair_targets(base[side], side, reps_list[a], reps_list[b])
        own_d, own_floor = int(d_tables[side][a, b]), bool(floor_tables[side][a, b])
        frozen_d = frozen_sweep.d_count(targets)
        classes, reps = frozen_floor.prep(list(targets))
        frozen_ok, _stats, _order = frozen_floor.subset_dfs(classes, reps)
        if (own_d, own_floor) != (frozen_d, bool(frozen_ok)):
            raise AssertionError(f"frozen pair-decision cross-check failed {name}/{side}/{a}/{b}")
        cross.append({"side": side, "a": a, "b": b, "d": own_d, "floor_exists": own_floor})

    # Recompute random full totals independently from scratch, not table lookup.
    # The first 12 also check the 48^3 right-monomial fiber and all three
    # post-sandwich sigma powers, whose equality supplies the exact multipliers.
    fresh_total_checks = []
    monomials = signed_monomials()
    sampled_indices = rng.sample(range(len(survivors)), min(24, len(survivors)))
    for ordinal, idx in enumerate(sampled_indices):
        p, q, r = (int(x) for x in survivors[idx])
        triple = honest_sandwich(U, V, W, reps_list[p], reps_list[q], reps_list[r])
        fresh = certified_total(*triple)
        expected_total = int(totals[idx])
        if int(fresh["total_lb"]) != expected_total:
            raise AssertionError(f"fresh total mismatch {name} survivor={idx}")
        row: dict[str, object] = {"survivor_index": idx, "total_lb": expected_total}
        if ordinal < 12:
            sigma_totals = [int(certified_total(*x)["total_lb"]) for x in sigma_orbit(*triple)]
            if sigma_totals != [expected_total] * 3:
                raise AssertionError(f"sigma-total invariance failed {name} survivor={idx}")
            mp, mq, mr = (monomials[rng.randrange(48)] for _ in range(3))
            P = matmul(reps_list[p], mp)
            Q = matmul(reps_list[q], mq)
            R = matmul(reps_list[r], mr)
            fiber = honest_sandwich(U, V, W, P, Q, R)
            fiber_total = int(certified_total(*fiber)["total_lb"])
            if fiber_total != expected_total or not all_ternary(*fiber):
                raise AssertionError(f"monomial-fiber total invariance failed {name} survivor={idx}")
            row["sigma_totals"] = sigma_totals
            row["fiber_total_lb"] = fiber_total
        fresh_total_checks.append(row)

    # Exact fmpz Brent rechecks on deterministic survivors; validity follows
    # algebraically for all, but this validates the implementation path.
    brent_checks = []
    for idx in rng.sample(range(len(survivors)), min(12, len(survivors))):
        p, q, r = (int(x) for x in survivors[idx])
        triple = honest_sandwich(U, V, W, reps_list[p], reps_list[q], reps_list[r])
        if not all_ternary(*triple):
            raise AssertionError("survivor full-map recheck left ternary alphabet")
        fi, ff = brent_failures_int(*triple), brent_failures_fmpz(*triple)
        if fi or ff:
            raise AssertionError("survivor full-map Brent recheck failed")
        brent_checks.append({"survivor_index": idx, "brent_int": 729, "brent_fmpz": 729})

    hist = Counter(int(x) for x in totals.tolist())
    sorted_totals = sorted(int(x) for x in totals.tolist())
    median = sorted_totals[len(sorted_totals) // 2] if sorted_totals else None
    min_total = min(sorted_totals) if sorted_totals else None
    max_total = max(sorted_totals) if sorted_totals else None
    minimizer_rows = np.flatnonzero(totals == min_total).tolist() if min_total is not None else []
    nonmono_mask = np.array([not all(int(x) == mono_idx for x in row) for row in survivors], dtype=np.bool_)
    nonmono_min = int(np.min(totals[nonmono_mask])) if np.any(nonmono_mask) else None

    summary = {
        "decomposition": name,
        "decision_formula": "U_lb + V_lb + Wfac_lb + 14, side_lb=d+[floor impossible]",
        "all_monomial_data_triple": {
            "index": int(mono_rows[0]),
            "total_lb": mono_total,
            "asserted_before_histogram": True,
        },
        "survivors_decided": len(survivors),
        "histogram": {str(k): int(hist[k]) for k in sorted(hist)},
        "min_total_lb": min_total,
        "max_total_lb": max_total,
        "upper_median_total_lb": median,
        "minimizer_count_data_triples": len(minimizer_rows),
        "minimizer_survivor_indices": minimizer_rows,
        "nonmonomial_involving_min_total_lb": nonmono_min,
        "flags_le_54": int(np.count_nonzero(totals <= 54)),
        "pair_decisions": {
            side: {
                "ternary_pairs": int(np.count_nonzero(tables[side])),
                "floor_exists": int(np.count_nonzero(floor_tables[side] == 1)),
                "floor_impossible_exact_census": int(np.count_nonzero(floor_tables[side] == 0)),
                "d_histogram": {
                    str(k): int(v)
                    for k, v in sorted(Counter(int(x) for x in d_tables[side][d_tables[side] >= 0].tolist()).items())
                },
                "max_dfs_states": int(np.max(state_tables[side])),
            }
            for side in ("U", "V", "W")
        },
        "independent_pair_crosschecks": cross,
        "fresh_total_crosschecks": fresh_total_checks,
        "survivor_brent_crosschecks": brent_checks,
        "seconds": time.perf_counter() - started,
    }
    arrays = {
        **{f"{side}_d": d_tables[side] for side in ("U", "V", "W")},
        **{f"{side}_floor": floor_tables[side] for side in ("U", "V", "W")},
        **{f"{side}_lb": lb_tables[side] for side in ("U", "V", "W")},
        **{f"{side}_states": state_tables[side] for side in ("U", "V", "W")},
        "survivors": survivors,
        "side_lbs": side_lbs,
        "totals": totals,
    }
    return arrays, summary


def render_pre_statement(campaign: Path, code_hash: str) -> str:
    full_per_decomp = 6960**3 * 3
    full_total = 2 * full_per_decomp
    data_total = 2 * 145**3
    return f"""# Pre-statement — new-decomposition off-diagonal census (`mws59`, `stapleton60`)

**Created {utc_stamp()} by the producer and committed BEFORE any target computation.**
Campaign: `{campaign}`.  Producer SHA-256 prefix: `{code_hash}`.

**Rule-5 correction / prior abort:** producer hash `4dc928ac6dbe` at campaign
`2026-08-31T082412Z_a2d87381-157e-45bb-ad3a-519ddf899bc1_4dc928ac6dbe`
aborted in 0.26 s before any anchor or target count: it attempted to iterate the
intentional `None` in `META["perminov58"]["split"]`.  That failed source and log are
frozen, with no target claim.  This code path instead decodes the pinned Perminov JSON,
asserts zero factor mismatches and `complexity.reduced == 58`, and receives this new hash
and separate pre-registration before rerun.

## 1. Fixed question and falsifiable outcomes

The completed off-diagonal programme covers only `paper55` and `sun56`; its undecided
remainder is empty.  This campaign therefore takes the highest-value listed extension:
a **different decomposition**.  It fixes exactly the two remaining on-disk, previously
729/729-verified public rank-23 decompositions: `mws59` (published 59 additions) and
`stapleton60` (published 60 additions).  `perminov58` is not a new decomposition here:
its factor blocks are the already-covered `paper55` blocks.

Exactly one verdict is admissible per fully decided data triple:

* **NOGO:** exact lower bound at least 55, via the exact d-counter, complete subset-DFS
  floor census, and the +14 transposition gap with its rank/activity precondition checked.
* **LIVE:** lower bound at most 54; halt before recording the row and escalate to Main.
  A low lower bound is not a scheme.  It then enters a separately frozen exact scheduling
  adjudication; any SAT witness is replayed over Z, and any reported UNSAT must have DRAT
  converted to LRAT and accepted by both `drat-trim` and `lrat-check`.
* **RECORD:** an explicit <=54 circuit, only after exact gate-list replay and all 729 Brent
  identities over Z; escalate to Main before writing anywhere.

The primary output is the exact survivor count and exact certified-lower-bound histogram
for each decomposition.  If a wall cap intervenes, only the lexicographic prefix is
reported and the exact remainder named; no total claim is made for it.

## 2. Domain fixed before results (rules 14-16)

Let T be the **set** (not a group) of all 6960 ternary 3x3 integer matrices with
determinant +/-1, exactly re-enumerated from all 3^9=19,683 ternary matrices.  For each
D in {{`mws59`, `stapleton60`}}, enumerate every independent triple
(P,Q,R) in T^3 and then every post-sandwich sigma power k in {{0,1,2}}, where
sigma(U,V,W)=(V,W^T,U^T).  Full raw domain: {full_total:,} scheme-instances
({full_per_decomp:,} per decomposition).  The 6960 set's failure to close under
multiplication is irrelevant: each P,Q,R is used independently and no product is assumed
to remain in T.

Each factor row is a row-major 3x3 block.  In **standard matrix semantics** the fixed
family is

    U' = P^-1 U Q^-T,   V' = Q^T V R^-T,   W' = P^T W R.

This is re-derived from F(A,B,C)=tr((AB)^T C), not tr(ABC).  Convention warning: frozen
`gatec_sweep.sandwich` uses transposed-L letters; on monomials frozen (X,Y,Z) equals the
standard family at (P,Q,R)=(X^T,Y^T,Z^T), checked bit-for-bit.

An instance is admitted only when all 23 rows of all three transformed blocks remain in
the alphabet {{-1,0,1}}^9.  Brent validity follows from the proved automorphism and is
also spot-checked 729/729 in two exact paths.  P^-1,Q^-1,R^-1 need not themselves be
ternary.  Non-ternary factor images are counted as excluded from the fixed alphabet, not
silently dropped.

## 3. Exact factorization and completeness

T is partitioned by free right multiplication by the 48 signed monomials M:
T=A*M, |A|=145, every orbit size 48, representatives chosen lexicographically.  If
P=a1 p1, Q=a2 p2, R=a3 p3, the outer p_i induce a common signed permutation of the
nine input coordinates on each side.  Ternarity, d, floor existence, activity, and
circuit cost are invariant.  Thus every data triple (a1,a2,a3) represents exactly
48^3 sandwich triples, and sigma contributes a further factor 3.  The campaign decides
exactly {data_total:,} (decomposition, data-triple) points across the two decompositions
before the proven fiber multiplier.
The factorization is accepted only after: all 6960 elements partition into 145 size-48
orbits with unique factorization; a full anchored 145^2 direct-vs-table subcube for each
decomposition has zero mismatches; 200 deterministic random direct/table triples per
decomposition agree; 24 numpy-vs-pure-Python entry-level checks agree; and 25 random
right-monomial fiber predicates per decomposition agree.  Any mismatch aborts before a
count is reported.

Counting uses exact int64 only after the hard entry bounds |U'|<=36, |V'|<=18,
|W'|<=9 are proved from |G|<=1, |G^-1|<=2; counts are accumulated in Python integers
and independently as a bounded int64 matrix contraction (<145^3<2^63).  No float,
interval arithmetic, HiGHS, SAT, or tolerance enters the main verdict.

## 4. Bound semantics and direction

For one side F, d(F) is the number of distinct non-input sign classes among its 23
needed forms.  Each addition creates at most one new class, so C(F)>=d(F).  At exactly d
gates every gate must be a needed class; complete memoized subset search enumerates every
ordering and every representation t=+/-a+/-b (including a=b).  Failure is therefore an
exact finite census proving C(F)>=d+1, not a solver UNSAT label.  A found floor schedule
is replayed gate by gate.

For the output map Wfac^T:23->9, reverse-mode transposition of an active A-gate circuit
uses A+9-23=A-14 additions for Wfac:9->23, hence C(output)>=C(Wfac)+14.  The permissive
direction is explicitly excluded.  Startup asserts every slot has 23 nonzero rows and
exact Q-rank 9; the honest action is invertible on the nine coordinates, so the
precondition holds on the full domain.  Total lower bound is

    U_lb + V_lb + Wfac_lb + 14.

Sigma only cycles the three side costs and transposes blocks (a coordinate permutation),
so the total is invariant; it multiplies the orientation count by exactly 3.  Twelve
survivor totals per decomposition are replayed through all three sigma powers and through
independent right-monomial fibers.

## 5. Mandatory anchors and counterfactual controls — before new counts

1. Reproduce 55/58/56/59/60 published totals and 729/729 Brent over Z in independent
   Python-int, fmpz, and frozen paths.
2. Assert the **all-monomial data triple is exactly 55** on `paper55` and `sun56`
   before any new histogram.  Also assert the inherited all-monomial certified totals
   56 for `mws59` and 58 for `stapleton60`, across their sigma classes.
3. Reproduce the frozen d/floor tuples for paper55, sun56, mws59, and stapleton60 with
   both this independent DFS and the frozen DFS.
4. Fixed nonmonomial G=[[-1,-1,-1],[0,1,1],[0,0,1]]: the final family must be ternary
   and Brent 729/729 on sun56; the wrong W'=G W G^T and transpose-as-inverse plants
   must fail.  A deliberate nonzero W-coefficient flip must produce Brent failures;
   det-2 must be rejected.  These are the controls that give a zero-event result meaning.

## 6. Caps, order, escalation, and non-goals

Lexicographic order in (decomposition,p,q,r), with mws59 before stapleton60.  Fixed seed
{SEED} only for controls.  Caps: anchors/controls 2 h, tables 8 h, decisions 24 h,
total 48 h; unfinished work is an exact named remainder.  Any <=54 row halts before the
row is written and is escalated to Main; no computation is adjusted toward a published
value.

Not searched: `paper55`/`sun56` beyond startup controls (their full domain is already
closed); decompositions outside the five named public ones (Laderman, other Smirnov,
Schwartz-Vaknin, unpublished, F2-only); non-ternary factor alphabets; GL(3,Q) or GL(3,Z)
beyond T; actions outside the proved family; anti-cyclic pair swaps (compute B*A); and
upper-bound circuit synthesis unless a <=54 adjudication triggers.  No radii exist: all
objects are exact integers.
"""


def init_campaign() -> Path:
    code_hash = sha256(THIS)[:12]
    campaign = CAMPAIGNS / f"{utc_stamp()}_{uuid.uuid4()}_{code_hash}"
    campaign.mkdir(parents=False, exist_ok=False)
    pre = render_pre_statement(campaign, code_hash)
    (campaign / "pre_statement.md").write_text(pre, encoding="utf-8")
    files = list(campaign.iterdir())
    if files != [campaign / "pre_statement.md"]:
        raise AssertionError("pre_statement.md was not the campaign's first and only init file")
    print(campaign)
    return campaign


def run_campaign(campaign: Path) -> dict[str, object]:
    campaign = campaign.resolve()
    if campaign.parent != CAMPAIGNS.resolve() or not (campaign / "pre_statement.md").exists():
        raise CampaignError("invalid campaign path or missing pre_statement.md")
    if (campaign / "campaign_result.json").exists():
        raise CampaignError("campaign already has a result; frozen campaigns are immutable")
    log = Logger(campaign / "run.log")
    t_campaign = time.perf_counter()
    try:
        log("RUN START; exact new-decomposition off-diagonal census")
        log(f"python={sys.version.split()[0]} numpy={np.__version__} code_sha256={sha256(THIS)}")
        if sha256(THIS)[:12] != campaign.name.rsplit("_", 1)[-1]:
            raise CampaignError("source hash no longer matches campaign directory suffix")

        phase_a_start = time.perf_counter()
        # Rule-17b startup order is literal: no population/table count is
        # computed before the known all-monomial 55 anchor is asserted.
        data, meta = load_decompositions()
        identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        anchors = anchor_all(data, meta, identity)
        log("MANDATORY STARTUP ANCHOR PASS: paper55=55, sun56=55, mws59=56, stapleton60=58")
        semantic = semantic_controls(data)
        unimodular, t_stats = enumerate_ternary_unimodular()
        monomials = signed_monomials()
        reps, _member_map, orbit_stats = build_data_orbits(unimodular, monomials)
        mono_idx = int(orbit_stats["monomial_data_node"])
        canonical_anchors = {}
        for name in CONTROL_DECOMPS:
            oriented = honest_sandwich(
                *data[name], reps[mono_idx], reps[mono_idx], reps[mono_idx]
            )
            got = int(certified_total(*oriented)["total_lb"])
            expected = int(EXPECTED_CLASSES[name][0][2])
            if got != expected:
                raise AssertionError(
                    f"canonical all-monomial data anchor failed for {name}: {got} != {expected}"
                )
            canonical_anchors[name] = got
        phase_a = {
            "ternary_matrix_census": t_stats,
            "data_factorization": orbit_stats,
            "anchors": anchors,
            "canonical_data_node_anchors": canonical_anchors,
            "semantic_controls": semantic,
            "seconds": time.perf_counter() - phase_a_start,
        }
        if phase_a["seconds"] > PHASE_A_CAP_S:
            raise CampaignError("phase A cap exceeded")
        atomic_json(campaign / "phaseA_controls.json", phase_a)
        log("PHASE A PASS: all anchors, counterfactual plants, ranks, and factorization checks")

        full_outputs = {}
        npz_pair_payload: dict[str, np.ndarray] = {
            "data_reps": np.asarray(reps, dtype=np.int8),
            "data_inverses": np.asarray([inverse_unimodular(a) for a in reps], dtype=np.int8),
        }
        for d_index, name in enumerate(TARGET_DECOMPS):
            table_start = time.perf_counter()
            log(f"{name}: building exact 145x145 pair tables")
            tables, table_info = build_pair_tables(data[name], reps)
            controls = direct_table_controls(data[name], reps, tables, mono_idx, SEED + d_index)
            if time.perf_counter() - table_start > TABLE_CAP_S:
                raise CampaignError(f"{name}: table cap exceeded")
            survivors, census = enumerate_survivors(tables, mono_idx)
            multiplier_per_sigma = 48**3
            census.update({
                "pair_tables": table_info,
                "controls": controls,
                "sandwich_instances_per_sigma": int(census["survivor_data_triples"]) * multiplier_per_sigma,
                "sandwich_instances_all_sigma": int(census["survivor_data_triples"]) * multiplier_per_sigma * 3,
                "full_raw_domain_per_decomposition": 6960**3 * 3,
                "sigma_multiplier": 3,
                "fiber_multiplier_per_sigma": 48**3,
                "seconds": time.perf_counter() - table_start,
            })
            for side in ("U", "V", "W"):
                npz_pair_payload[f"{name}_{side}_ternary"] = tables[side]
            npz_pair_payload[f"{name}_survivors"] = survivors
            log(f"{name}: census controls pass; survivor data-triples={len(survivors)}")

            arrays, decision = decide_pairs_and_triples(
                name, data[name], reps, tables, survivors, mono_idx, log
            )
            np.savez_compressed(campaign / f"decision_{name}.npz", **arrays)
            atomic_json(campaign / f"decision_{name}_summary.json", decision)
            full_outputs[name] = {"census": census, "decision": decision}
            log(
                f"{name}: decision complete, min={decision['min_total_lb']} "
                f"max={decision['max_total_lb']} flags<=54={decision['flags_le_54']}"
            )
            if time.perf_counter() - t_campaign > CAMPAIGN_CAP_S:
                raise CampaignError("campaign cap exceeded")

        np.savez_compressed(campaign / "pair_tables.npz", **npz_pair_payload)

        named_valid = sum(int(x["census"]["sandwich_instances_all_sigma"]) for x in full_outputs.values())
        named_data = sum(int(x["census"]["survivor_data_triples"]) for x in full_outputs.values())
        global_min = min(int(x["decision"]["min_total_lb"]) for x in full_outputs.values())
        global_flags = sum(int(x["decision"]["flags_le_54"]) for x in full_outputs.values())
        nonmono_values = [
            int(x["decision"]["nonmonomial_involving_min_total_lb"])
            for x in full_outputs.values()
            if x["decision"]["nonmonomial_involving_min_total_lb"] is not None
        ]
        nonmono_min = min(nonmono_values) if nonmono_values else None
        result = {
            "campaign": campaign.name,
            "code_sha256": sha256(THIS),
            "domain": {
                "decompositions": list(TARGET_DECOMPS),
                "T_size": 6960,
                "data_nodes": 145,
                "full_raw_scheme_instances": 2 * 6960**3 * 3,
                "data_triples_enumerated": 2 * 145**3,
                "fiber_multiplier_per_sigma": 48**3,
                "sigma_multiplier": 3,
                "action_standard": "U'=P^-1 U Q^-T; V'=Q^T V R^-T; W'=P^T W R",
                "alphabet": [-1, 0, 1],
            },
            "outputs": full_outputs,
            "aggregate": {
                "survivor_data_triples": named_data,
                "valid_scheme_instances_all_sigma": named_valid,
                "min_total_lb": global_min,
                "nonmonomial_involving_min_total_lb": nonmono_min,
                "flags_le_54": global_flags,
            },
            "evidence": "MACHINE-VERIFIED exact integer census + complete finite subset DFS; no float",
            "seconds": time.perf_counter() - t_campaign,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        if global_flags:
            raise RecordFlag("A <=54 row was detected; notify Main before recording aggregate")
        atomic_json(campaign / "campaign_result.json", result)
        log(f"RUN PASS: valid orientations={named_valid}, aggregate min_lb={global_min}, zero <=54")
        return result
    finally:
        log.close()


def get_git_commit_for(path: Path) -> str:
    import subprocess

    p = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(path)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    commit = p.stdout.strip()
    if not commit:
        raise CampaignError(f"no committed provenance for {path}")
    return commit


def package_version(distribution: str) -> str | None:
    try:
        import importlib.metadata
        return importlib.metadata.version(distribution)
    except Exception:
        return None


def rule7_sentence(result: dict[str, object]) -> str:
    agg = result["aggregate"]
    return (
        "This campaign swept exactly the two decompositions `mws59` and `stapleton60` "
        "under every independent (P,Q,R) in T^3, where T is the exactly enumerated "
        "6960-element set of ternary 3x3 integer determinant-+/-1 matrices, followed by "
        "all three post-sandwich sigma powers, for 2,022,921,216,000 raw scheme-instances "
        "represented without sampling by the proved 145^3 data-node factorization per "
        f"decomposition; it admitted exactly {agg['survivor_data_triples']:,} ternary "
        f"data-triples ({agg['valid_scheme_instances_all_sigma']:,} orientations after the "
        "48^3 and 3-sigma multiplicities), and decided every admitted data-triple with exact "
        "integer d-counting, complete subset-DFS floor census, and the rank-audited +14 "
        "transposition gap under the standard action U'=P^-1 U Q^-T, V'=Q^T V R^-T, "
        "W'=P^T W R; the alphabet was exactly {-1,0,1}, parameters were the two fixed "
        "decompositions, all 6960 choices for each of P,Q,R, and sigma=0,1,2, radii were "
        "none, and no instance in that named domain was sampled or left undecided; it did "
        "not sweep paper55/sun56 beyond startup controls, any decomposition outside the five "
        "named public ones (including Laderman, other Smirnov schemes, Schwartz-Vaknin, "
        "unpublished, or F2-only constructions), non-ternary factor alphabets, GL(3,Q) or "
        "GL(3,Z) sandwiches beyond T, actions outside the proved family, anti-cyclic swaps "
        "that compute B*A, or upper-bound circuit synthesis, so it establishes no universal "
        "no-54 claim outside the named domain."
    )


def render_report(result: dict[str, object], campaign: Path) -> str:
    mws = result["outputs"]["mws59"]
    sta = result["outputs"]["stapleton60"]
    agg = result["aggregate"]
    rule7 = rule7_sentence(result)
    if int(agg["min_total_lb"]) >= 56:
        headline = "The paper55 all-monomial data triple remains the unique 55 lower-bound minimizer across the four distinct on-disk public decompositions."
    else:
        headline = "A new 55 lower-bound data triple exists; this is not an upper-bound circuit."
    return f"""# Report — new-decomposition off-diagonal census

## Headline

**{headline}**  No <=54 flag occurred.  `[MACHINE-VERIFIED]`

Frozen campaign: `{campaign}`.

## Exact verdicts

* `mws59`: `{mws['census']['survivor_data_triples']:,}` survivor data-triples,
  `{mws['census']['sandwich_instances_all_sigma']:,}` valid orientations after the
  48^3 monomial-fiber and 3-sigma multipliers; certified lower-bound range
  **{mws['decision']['min_total_lb']}..{mws['decision']['max_total_lb']}**, upper median
  **{mws['decision']['upper_median_total_lb']}**, `{mws['decision']['flags_le_54']}` at <=54.
  `[MACHINE-VERIFIED]`
* `stapleton60`: `{sta['census']['survivor_data_triples']:,}` survivor data-triples,
  `{sta['census']['sandwich_instances_all_sigma']:,}` valid orientations; certified
  lower-bound range **{sta['decision']['min_total_lb']}..{sta['decision']['max_total_lb']}**,
  upper median **{sta['decision']['upper_median_total_lb']}**,
  `{sta['decision']['flags_le_54']}` at <=54. `[MACHINE-VERIFIED]`
* Aggregate new domain: `{agg['survivor_data_triples']:,}` survivor data-triples,
  `{agg['valid_scheme_instances_all_sigma']:,}` valid orientations; minimum certified
  lower bound **{agg['min_total_lb']}** and every non-monomial-involving survivor at least
  **{agg['nonmonomial_involving_min_total_lb']}**. `[MACHINE-VERIFIED]`

The 59- and 60-addition published witnesses were reproduced first, 729/729 Brent over Z
in Python-int, fmpz, and frozen paths.  The all-monomial certified-lower-bound anchors
were 56 (`mws59`) and 58 (`stapleton60`); the mandatory paper55/sun56 55 anchors were
asserted before histogram computation. `[REPRODUCED] [MACHINE-VERIFIED]`

These numbers are **lower bounds**, not synthesized circuits.  A row at lower bound 55
would not by itself be a 55-addition upper bound; no <=54 lower-bound row appeared, so no
record adjudication or escalation was triggered.

## Instrument adequacy

The corrected family was re-derived from `F(A,B,C)=tr((AB)^T C)` and the exact direction
of both lower bounds was checked.  The 6960 matrix census, 145 size-48 right-monomial
orbits, and unique factorization were recomputed.  For each decomposition a full 145^2
anchored direct-vs-table subcube, 200 random direct/table comparisons, 24 entrywise
numpy-vs-pure-Python comparisons, and 25 right-monomial fiber-predicate comparisons had
zero mismatches; 12 survivor totals were replayed through independent fibers and every
sigma power.  The fixed nonmonomial final-family
plant passed Brent 729/729; the wrong-W, transpose-as-inverse, direct coefficient-corruption,
and det-2 plants were all rejected.  Exact rank 9 and 23 nonzero rows in every slot justify
the +14 transposition gap throughout. `[MACHINE-VERIFIED]`

No SAT solver was invoked; therefore there is no unproved UNSAT label.  Every
floor-impossibility verdict is a complete finite subset census with an explicit
model<->schedule argument.  No float or tolerance entered any result.

## Rule-5 implementation incident

The first producer hash `4dc928ac6dbe` aborted before its first anchor because it treated
Perminov's intentional missing left/right/output split as iterable.  Its campaign, source,
log, and no-claim verdict remain frozen.  This successful source replays the 58 anchor from
the pinned JSON's `complexity.reduced` field and asserts its decoded factors match exactly.

## Rule-7 scope sentence

{rule7}

## Prompt premises disproved

None.  The prompt's premise that the audited machinery could be reused on a different
on-disk decomposition survived all controls.  The non-closure of the 6960 set was handled
explicitly and never assumed away.

## Named next campaign and cost

`laderman23_source_lock`: lock a primary/canonical factor file for Laderman's original
rank-23 decomposition, reproduce a published addition split and all 729 Brent identities
over Z, then pre-register its ternary-unimodular off-diagonal census if and only if the
factor alphabet fits this model.  Estimated cost: 2-6 h source/provenance and convention
resolution, then 15-60 min census/decision at the observed 145-node rate (longer only if
its floor-DFS state counts are materially larger).
"""


def freeze_campaign(campaign: Path) -> None:
    campaign = campaign.resolve()
    result_path = campaign / "campaign_result.json"
    if not result_path.exists():
        raise CampaignError("cannot freeze without campaign_result.json")
    if (campaign / "checksums.sha256").exists():
        raise CampaignError("campaign already frozen")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    inputs = campaign / "inputs"
    inputs.mkdir(exist_ok=False)
    to_copy = [
        THIS,
        SRC / "gatec_decomps.py",
        SRC / "gate_b_floor.py",
        SRC / "gatec_sweep.py",
        SRC / "verify_anchors.py",
        SRC / "gate_a.py",
        SRC / "tensor_data.py",
        SRC / "stapleton60_data.py",
        SCRATCH / "mws59_layout.txt",
        SCRATCH / "cr58_cn122_ZT_reduced.json",
        SCRATCH / "sun56_verify.py",
    ]
    provenance = {}
    for source in to_copy:
        if not source.exists():
            raise CampaignError(f"missing freeze input {source}")
        dest = inputs / (source.name + ".asrun")
        shutil.copyfile(source, dest)
        if source == THIS and sha256(dest) != result["code_sha256"]:
            raise CampaignError("asrun producer does not match run hash")
        provenance[str(dest.relative_to(campaign))] = {
            "source": str(source),
            "sha256": sha256(dest),
        }

    report = render_report(result, campaign)
    (campaign / "report.md").write_text(report, encoding="utf-8")
    manifest = {
        "campaign": campaign.name,
        "status": "FROZEN",
        "freeze_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "pre_statement_first_file": True,
        "pre_statement_sha256": sha256(campaign / "pre_statement.md"),
        "pre_statement_commit": get_git_commit_for(campaign / "pre_statement.md"),
        "source_sha256": result["code_sha256"],
        "domain": result["domain"],
        "aggregate": result["aggregate"],
        "evidence": result["evidence"],
        "rule7_scope_sentence": rule7_sentence(result),
        "toolchain": {
            "python": sys.version,
            "python_flint": package_version("python-flint"),
            "numpy": np.__version__,
            "sympy": package_version("sympy"),
            "kissat": "not invoked",
            "cadical": "not invoked",
            "drat_trim": "not invoked (no SAT UNSAT claim)",
            "lrat_check": "not invoked (no SAT UNSAT claim)",
        },
        "machine": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "hostname": socket.gethostname(),
            "single_process": True,
        },
        "inputs": provenance,
        "timing_seconds_float_computational_evidence_only": result["seconds"],
        "immutable_after_checksums": True,
    }
    atomic_json(campaign / "manifest.json", manifest)

    files = sorted(p for p in campaign.rglob("*") if p.is_file() and p.name != "checksums.sha256")
    lines = [f"{sha256(p)}  {p.relative_to(campaign)}" for p in files]
    (campaign / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(campaign)


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--init", action="store_true")
    group.add_argument("--run", metavar="CAMPAIGN")
    group.add_argument("--freeze", metavar="CAMPAIGN")
    args = parser.parse_args()
    if args.init:
        init_campaign()
    elif args.run:
        run_campaign(Path(args.run))
    else:
        freeze_campaign(Path(args.freeze))


if __name__ == "__main__":
    main()
