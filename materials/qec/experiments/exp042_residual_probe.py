"""EXP-042: exhaustive small-lattice probe of the residual class T(P) < k_P/2.

Conjecture B-prime lives at T(P) < k_P/2.  Across all 134 EXP-039 certificates
exactly one parent sits there (9a7638586033*, ell=12 x m=6, k_P=12, T=4), and
its lattice has ell*m = 72.  This experiment asks whether *small* lattices
(ell*m <= 24) contain the residual class at all, and if so pounds every such
parent with exhaustive perturbation tests.  Everything is exact GF(2)
arithmetic; no SAT solver, no Stim, single-threaded NumPy.

Arms
----
1. **Exhaustive parent enumeration.**  For every ordered lattice shape
   (ell, m) with 2 <= ell, 2 <= m, ell*m <= 24 and every unordered support
   family (w_A, w_B) in {(2,2), (2,3), (3,3)} (the catalogue is uniformly
   (3,3); (2,2)/(2,3) extend the sweep), enumerate all translation-canonical
   supports: two supports are identified when a lattice translate x^a y^b maps
   one to the other.  k_P, d_Z and T are invariant under independent block
   translations (a translate of [A|B] is a qubit relabelling), so one
   canonical representative per class is exhaustive for classification.  Pairs
   {A, B} are unordered (block swap is again a qubit relabelling).  Filters,
   matching EXP-040's nondegeneracy predicate: k_P >= 2 and a single Tanner
   connected component.

2. **Exact d_Z and T for every surviving parent.**  d_Z by ascending-weight
   meet-in-the-middle enumeration of ker[A B] (side tables and matching
   adapted from EXP-041, which was itself cross-checked bit-exactly against
   full nullspace spans), survivors filtered against S_Z = rowspace(H_Z) by
   echelon reduction.  T by inserting the translation orbits of *every*
   minimum-weight survivor into an echelon basis seeded with S_Z (EXP-040's
   definition); every translate is separately verified to lie in ker[A B].

3. **Residual-class perturbation test.**  For every parent with 2T < k_P:
   all delta>0 perturbations of total [C D] support <= 3 (i) must keep
   dim(Delta_bar) <= T (the catalogue upper law) and (ii) may not produce a
   strict distance increase d_Q > d_Z with dim(Delta_bar) != T (Conjecture
   B-prime).  d_Q > d_Z is decided exactly: if the minimum-logical module is
   not absorbed, a surviving pure-Z logical of weight d_Z caps the child;
   otherwise an exhaustive symplectic meet-in-the-middle search over all
   ternary Pauli errors of weight <= d_Z runs under an explicit pattern
   budget, and exceeding the budget is recorded as unresolved (never silently
   treated as a pass or a violation).

Independence protocol: the module never imports the SAT path; ``pysat`` /
``sat_decide`` must be absent from ``sys.modules`` at write time.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import os
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from math import comb
from pathlib import Path
from typing import Any, Iterable

# Shared workstation: everything here is exact GF(2) enumeration and must stay
# single-threaded (hard cap: one core per agent).
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_env] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SCHEMA = "exp042-residual-probe-v1"
OUT = ROOT / "results" / "processed" / "exp042_residual_probe.json"

MAX_L = 24                 # lattice area cap: ell*m <= 24
MIN_LATTICE_DIM = 2        # both lattice sides >= 2 (catalogue minimum is 3)
SUPPORT_FAMILIES = ((2, 2), (2, 3), (3, 3))  # unordered (w_A, w_B), w_A <= w_B
D_WEIGHT_CAP = 12          # d_Z search cap; exhausting it is recorded loudly
SIDE_WEIGHT_CAP = 3_000_000  # single-side mask table cap for one weight
MATCH_GUARD = 400_000      # max matches per (parent, total weight) split set
PERT_SUPPORT_CAP = 3       # total [C D] support cap for residual probing
SYM_PATTERN_BUDGET = 2e7   # per-side ternary pattern budget for d_Q searches
SPAN_CROSS_CHECK_PARENTS = 8   # full-nullspace cross-checks of the MITM set
SPAN_CROSS_CHECK_MAX_KER = 20

# Reuse the verified EXP-041 kernel machinery read-only (its own module level
# only loads E27 and imports numpy helpers; no solver anywhere).
_SPEC_041 = importlib.util.spec_from_file_location(
    "exp041_t_crosscheck", ROOT / "experiments" / "exp041_t_crosscheck.py"
)
E41 = importlib.util.module_from_spec(_SPEC_041)
sys.modules[_SPEC_041.name] = E41
_SPEC_041.loader.exec_module(E41)
E27 = E41.E27

from qec_research.codes.bicycle import (  # noqa: E402
    BBSpec,
    bb_stabilizer,
    commutation_defect,
    poly_matrix,
)
from qec_research.codes.pbb_survival import translation_orbit  # noqa: E402
from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


# ---------------------------------------------------------------------------
# Small GF(2) / bitset helpers (EXP-040 conventions)
# ---------------------------------------------------------------------------
def row_mask(row: np.ndarray) -> int:
    bits = np.flatnonzero(np.asarray(row, dtype=np.uint8).reshape(-1) & 1)
    out = 0
    for bit in bits:
        out |= 1 << int(bit)
    return out


def mask_row(mask: int, ncols: int) -> np.ndarray:
    return np.fromiter(
        ((int(mask) >> j) & 1 for j in range(ncols)),
        dtype=np.uint8,
        count=ncols,
    )


def rows_to_masks(M: np.ndarray) -> list[int]:
    return [row_mask(row) for row in np.asarray(M, dtype=np.uint8)]


def gf2_product(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return (
        np.asarray(left, dtype=np.uint8).astype(np.int64)
        @ np.asarray(right, dtype=np.uint8).astype(np.int64)
        % 2
    ).astype(np.uint8)


def json_terms(terms: Iterable[tuple[int, int]]) -> list[list[int]]:
    return [[int(a), int(b)] for a, b in terms]


def support_terms(
    combo: Iterable[int], ell: int, m: int
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    block = ell * m
    c_terms: list[tuple[int, int]] = []
    d_terms: list[tuple[int, int]] = []
    for coordinate in combo:
        if coordinate < block:
            c_terms.append((coordinate // m, coordinate % m))
        else:
            q = coordinate - block
            d_terms.append((q // m, q % m))
    return tuple(c_terms), tuple(d_terms)


def rank_of_masks(rows: Iterable[int]) -> int:
    basis = E41.Echelon()
    for row in rows:
        basis.insert(int(row))
    return basis.rank


def matrix_rank(M: np.ndarray) -> int:
    return rank_of_masks(rows_to_masks(M))


# ---------------------------------------------------------------------------
# Translation-canonical support enumeration
# ---------------------------------------------------------------------------
def canonical_supports(ell: int, m: int, weight: int) -> list[tuple[int, ...]]:
    """One representative per orbit of weight-w supports under Z_ell x Z_m."""

    dim = ell * m
    # Precompute translation lookup: position p = i*m + j -> p + (a,b).
    translate: list[list[int]] = []
    for a in range(ell):
        for b in range(m):
            table = [0] * dim
            for i in range(ell):
                for j in range(m):
                    p = i * m + j
                    table[p] = ((i + a) % ell) * m + ((j + b) % m)
            translate.append(table)

    seen: set[frozenset[int]] = set()
    reps: list[tuple[int, ...]] = []
    for combo in itertools.combinations(range(dim), weight):
        if frozenset(combo) in seen:
            continue
        orbit = [tuple(sorted(table[p] for p in combo)) for table in translate]
        best = min(orbit)
        for member in orbit:
            seen.add(frozenset(member))
        if best not in reps:
            reps.append(best)
    return sorted(reps)


def canonical_rep(ell: int, m: int, support: Iterable[tuple[int, int]]) -> tuple[int, ...]:
    """Canonicalise an explicit term list (used for the catalogue coverage check)."""

    base = frozenset((a % ell) * m + (b % m) for a, b in support)
    best: tuple[int, ...] | None = None
    for a in range(ell):
        for b in range(m):
            shifted = tuple(
                sorted(((p // m + a) % ell) * m + ((p % m + b) % m) for p in base)
            )
            if best is None or shifted < best:
                best = shifted
    if best is None:
        raise ValueError("empty support")
    return best


# ---------------------------------------------------------------------------
# Exact d_Z (MITM) and T (translation orbits) for one parent
# ---------------------------------------------------------------------------
class SideTable:
    """Per-parent cache of EXP-041 side enumerations, keyed by (side, weight)."""

    def __init__(self, L: int, cols_a: np.ndarray, cols_b: np.ndarray) -> None:
        self.L = L
        self.cols_a = cols_a
        self.cols_b = cols_b
        self.data: dict[tuple[str, int], tuple[np.ndarray, np.ndarray]] = {}

    def get(
        self, side: str, weight: int, cap: int | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        key = (side, weight, cap)
        cached = self.data.get(key)
        if cached is None:
            cols = self.cols_a if side == "A" else self.cols_b
            limit = SIDE_WEIGHT_CAP if cap is None else cap
            if comb(self.L, weight) > limit:
                raise OverflowError(f"side table C({self.L},{weight}) too large")
            cached = E41.side_enumeration(self.L, weight, cols)
            self.data[key] = cached
        return cached

    def join_matches(self, w1: int, w2: int, cap: int | None = None) -> list[int]:
        """Vectorized MITM join; identical set to matches(), numpy sort-join."""
        masks_a, synd_a = self.get("A", w1, cap)
        masks_b, synd_b = self.get("B", w2, cap)
        order = np.argsort(synd_b, kind="stable")
        sorted_synd = synd_b[order]
        lo = np.searchsorted(sorted_synd, synd_a, side="left")
        hi = np.searchsorted(sorted_synd, synd_a, side="right")
        counts = hi - lo
        left_idx = np.flatnonzero(counts)
        if left_idx.size == 0:
            return []
        sorted_masks = masks_b[order]
        out: list[int] = []
        L = self.L
        for i in left_idx:
            lm = int(masks_a[i])
            start = int(lo[i])
            for k in range(int(counts[i])):
                out.append(lm | (int(sorted_masks[start + k]) << L))
        return out

    def matches(self, w1: int, w2: int) -> list[int]:
        """All z=(z1|z2) with wt(z1)=w1, wt(z2)=w2 and A z1 = B z2 (E41 logic)."""
        masks_a, synd_a = self.get("A", w1)
        masks_b, synd_b = self.get("B", w2)
        keys = np.intersect1d(np.unique(synd_a), np.unique(synd_b))
        if keys.size == 0:
            return []
        left_idx = np.nonzero(np.isin(synd_a, keys))[0]
        right_idx = np.nonzero(np.isin(synd_b, keys))[0]
        right_map: dict[int, list[int]] = {}
        for j in right_idx:
            right_map.setdefault(int(synd_b[j]), []).append(int(masks_b[j]))
        out: list[int] = []
        L = self.L
        for i in left_idx:
            lm = int(masks_a[i])
            for rm in right_map[int(synd_a[i])]:
                out.append(lm | (rm << L))
        return out

COSET_CAP_BITS = 16      # max dim ker(A)/ker(B) for coset enumeration
SPLIT_WORK_CAP = 2e7     # max enumerated masks per split option
JOIN_ROW_CAP = 8_000_000  # per-side rows for the vectorized MITM join


class ColumnSolver:
    """Solve ``M z = s`` over GF(2) for a square matrix M.

    Columns are echelonized (pivot = highest row bit) with combination
    tracking: a solve returns a particular solution plus a basis of the
    kernel, both as column bitmasks (bit j = column j).
    """

    def __init__(self, M: np.ndarray) -> None:
        self.L = int(M.shape[0])
        cols = [int(c) for c in E41.column_masks(np.asarray(M, dtype=np.uint8))]
        basis: dict[int, tuple[int, int]] = {}
        kernel: list[int] = []
        for j in range(self.L):
            residual, combo = cols[j], 1 << j
            while residual:
                pivot = residual.bit_length() - 1
                entry = basis.get(pivot)
                if entry is None:
                    break
                residual ^= entry[0]
                combo ^= entry[1]
            if residual:
                basis[pivot] = (residual, combo)
            else:
                kernel.append(combo)
        self.basis = basis
        self.kernel = kernel

    def solve(self, s: int) -> tuple[int, list[int]] | None:
        residual, combo = s, 0
        while residual:
            pivot = residual.bit_length() - 1
            entry = self.basis.get(pivot)
            if entry is None:
                return None
            residual ^= entry[0]
            combo ^= entry[1]
        return combo, self.kernel


def masks_of_weight(L: int, w: int):
    if w == 0:
        yield 0
        return
    for combo in itertools.combinations(range(L), w):
        mask = 0
        for p in combo:
            mask |= 1 << p
        yield mask


def coset_masks(particular: int, kernel: list[int]):
    """particular XOR u for every u in the span of the kernel basis."""

    yield particular
    u = particular
    for i in range(1, 1 << len(kernel)):
        g = i ^ (i >> 1)
        prev = (i - 1) ^ ((i - 1) >> 1)
        u ^= kernel[(g ^ prev).bit_length() - 1]
        yield u


def kernel_by_weight(solver: ColumnSolver) -> dict[int, list[int]] | None:
    """Kernel masks bucketed by Hamming weight, or None if the kernel is big."""

    if len(solver.kernel) > COSET_CAP_BITS:
        return None
    buckets: dict[int, list[int]] = {}
    for u in coset_masks(0, solver.kernel):
        buckets.setdefault(u.bit_count(), []).append(u)
    return buckets


def split_matches_solved(
    L: int,
    w1: int,
    w2: int,
    cols_b: np.ndarray,
    cols_a: np.ndarray,
    solver_a: ColumnSolver,
    solver_b: ColumnSolver,
    table: "SideTable",
) -> list[int]:
    """All z=(z1|z2) with wt(z1)=w1, wt(z2)=w2 and A z1 = B z2.

    Engine choice per split, cheapest first: enumerate one side and solve the
    other over a coset of its kernel (exact; kernel dims are tiny on these
    lattices), or the vectorized MITM join when both side tables fit the row
    cap.  The EXP-041 dict join is the last-resort fallback.
    """

    options: list[tuple[int, str]] = []
    if len(solver_a.kernel) <= COSET_CAP_BITS and comb(L, w2) <= SIDE_WEIGHT_CAP:
        options.append((comb(L, w2) << len(solver_a.kernel), "left"))
    if len(solver_b.kernel) <= COSET_CAP_BITS and comb(L, w1) <= SIDE_WEIGHT_CAP:
        options.append((comb(L, w1) << len(solver_b.kernel), "right"))
    if comb(L, w1) <= JOIN_ROW_CAP and comb(L, w2) <= JOIN_ROW_CAP:
        options.append((comb(L, w1) + comb(L, w2), "join"))
    options.sort()
    for _, which in options:
        out: list[int] = []
        if which == "left":
            for z2 in masks_of_weight(L, w2):
                s = 0
                y = z2
                while y:
                    low = y & -y
                    s ^= int(cols_b[low.bit_length() - 1])
                    y ^= low
                solved = solver_a.solve(s)
                if solved is None:
                    continue
                particular, kernel = solved
                for z1 in coset_masks(particular, kernel):
                    if z1.bit_count() == w1:
                        out.append(z1 | (z2 << L))
        elif which == "join":
            return table.join_matches(w1, w2, cap=JOIN_ROW_CAP)
        else:
            for z1 in masks_of_weight(L, w1):
                s = 0
                x = z1
                while x:
                    low = x & -x
                    s ^= int(cols_a[low.bit_length() - 1])
                    x ^= low
                solved = solver_b.solve(s)
                if solved is None:
                    continue
                particular, kernel = solved
                for z2 in coset_masks(particular, kernel):
                    if z2.bit_count() == w2:
                        out.append(z1 | (z2 << L))
        return out
    return table.matches(w1, w2)


def syndrome_of(z: int, L: int, cols_a: np.ndarray, cols_b: np.ndarray) -> int:
    s = 0
    x = z & ((1 << L) - 1)
    while x:
        low = x & -x
        s ^= int(cols_a[low.bit_length() - 1])
        x ^= low
    y = z >> L
    while y:
        low = y & -y
        s ^= int(cols_b[low.bit_length() - 1])
        y ^= low
    return s


def collect_by_weight(HX: np.ndarray, weights: range, use_solver: bool) -> dict[int, list[int]]:
    """Kernel masks by total weight via one engine (solver or pure MITM)."""

    L = HX.shape[1] // 2
    A = np.ascontiguousarray(HX[:, :L], dtype=np.uint8)
    B = np.ascontiguousarray(HX[:, L:], dtype=np.uint8)
    cols_a = E41.column_masks(A)
    cols_b = E41.column_masks(B)
    table = SideTable(L, cols_a, cols_b)
    solver_a = ColumnSolver(A) if use_solver else None
    solver_b = ColumnSolver(B) if use_solver else None
    ker_a = kernel_by_weight(solver_a) if use_solver else None
    ker_b = kernel_by_weight(solver_b) if use_solver else None
    by_weight: dict[int, list[int]] = {0: [0]}
    for weight in weights:
        collected: list[int] = []
        for w1 in range(weight + 1):
            w2 = weight - w1
            if use_solver and w1 and w2:
                matches = split_matches_solved(
                    L, w1, w2, cols_b, cols_a, solver_a, solver_b, table
                )
            elif use_solver and w1:
                matches = list(ker_a.get(w1, []))
            elif use_solver:
                matches = [z2 << L for z2 in ker_b.get(w2, [])]
            else:
                matches = table.matches(w1, w2)
            collected.extend(matches)
        by_weight[weight] = sorted(set(collected))
    return by_weight


def analyse_parent(
    HX: np.ndarray,
    HZ: np.ndarray,
    ell: int,
    m: int,
    k_parent: int,
) -> dict[str, Any]:
    """Exact d_Z (MITM, ascending weight) and T (orbit span mod S_Z)."""

    started = time.perf_counter()
    L = ell * m
    A = np.ascontiguousarray(HX[:, :L], dtype=np.uint8)
    B = np.ascontiguousarray(HX[:, L:], dtype=np.uint8)
    cols_a = E41.column_masks(A)
    cols_b = E41.column_masks(B)

    sz_basis = E41.Echelon(rows_to_masks(HZ))

    # Sanity: every H_Z row lies in ker[A B] (exact syndrome check).
    for row in rows_to_masks(HZ):
        if syndrome_of(row, L, cols_a, cols_b) != 0:
            raise AssertionError("H_Z row outside ker[A B]")

    table = SideTable(L, cols_a, cols_b)
    solver_a = ColumnSolver(A)
    solver_b = ColumnSolver(B)
    ker_a = kernel_by_weight(solver_a)
    ker_b = kernel_by_weight(solver_b)
    by_weight: dict[int, list[int]] = {0: [0]}
    d_z: int | None = None
    minimum_logicals: list[int] = []
    overflowed = False
    for weight in range(1, D_WEIGHT_CAP + 1):
        collected: list[int] = []
        for w1 in range(weight + 1):
            w2 = weight - w1
            if w1 and w2:
                try:
                    matches = split_matches_solved(
                        L, w1, w2, cols_b, cols_a, solver_a, solver_b, table
                    )
                except OverflowError:
                    overflowed = True
                    break
            elif w1:  # z = (z1 | 0) with A z1 = 0
                if ker_a is not None:
                    matches = list(ker_a.get(w1, []))
                else:
                    try:
                        matches = table.matches(w1, 0)
                    except OverflowError:
                        overflowed = True
                        break
            else:  # z = (0 | z2) with B z2 = 0
                if ker_b is not None:
                    matches = [z2 << L for z2 in ker_b.get(w2, [])]
                else:
                    try:
                        matches = table.matches(0, w2)
                    except OverflowError:
                        overflowed = True
                        break
            if len(matches) > MATCH_GUARD:
                raise RuntimeError(
                    f"match guard exceeded at weight {weight}: {len(matches)}"
                )
            collected.extend(matches)
        if overflowed:
            break
        distinct = sorted(set(collected))
        by_weight[weight] = distinct
        survivors = [z for z in distinct if not sz_basis.contains(z)]
        if survivors:
            d_z = weight
            minimum_logicals = survivors
            break
    if d_z is None:
        return {
            "d_Z": None,
            "T": None,
            "unresolved": "weight cap exhausted" if not overflowed else "side table cap",
            "wall_time_s": round(time.perf_counter() - started, 6),
        }

    # T: orbits of ALL minimum-weight logicals, inserted mod S_Z; every
    # translate is verified to lie in ker[A B] exactly (not on a prefix).
    module_basis = sz_basis.copy()
    module_rows: list[int] = []
    t_value = 0
    for logical in minimum_logicals:
        for translated in E41.orbit_masks(logical, ell, m):
            if syndrome_of(translated, L, cols_a, cols_b) != 0:
                raise AssertionError("translate of minimum logical left ker[A B]")
            module_rows.append(translated)
            if module_basis.insert(translated):
                t_value += 1
    T = int(module_basis.rank - sz_basis.rank)
    if t_value != T or T > k_parent:
        raise AssertionError("orbit insertion inconsistent with quotient rank")
    if T < 1:
        raise AssertionError("T = 0 despite a nontrivial minimum logical")

    return {
        "d_Z": int(d_z),
        "T": T,
        "minimum_logicals": minimum_logicals,
        "module_rows": module_rows,
        "by_weight": {w: v for w, v in by_weight.items()},
        "rank_hz": int(sz_basis.rank),
        "num_minimum_logicals": len(minimum_logicals),
        "num_module_rows": len(module_rows),
        "wall_time_s": round(time.perf_counter() - started, 6),
    }


# ---------------------------------------------------------------------------
# Exact symplectic threshold decision (d_Q > d_Z?) for candidate children
# ---------------------------------------------------------------------------
def xor_lookup(contributions: list[int]) -> np.ndarray:
    lookup = np.zeros(1 << len(contributions), dtype=np.uint64)
    for mask in range(1, lookup.shape[0]):
        low = mask & -mask
        bit = low.bit_length() - 1
        lookup[mask] = lookup[mask ^ low] ^ np.uint64(contributions[bit])
    return lookup


def ternary_patterns(qubits: int, weight: int) -> tuple[np.ndarray, np.ndarray]:
    """All ternary Pauli patterns of exact symplectic weight on ``qubits``."""

    total = comb(qubits, weight) * (3**weight)
    xs = np.empty(total, dtype=np.uint32)
    zs = np.empty(total, dtype=np.uint32)
    cursor = 0
    for support in itertools.combinations(range(qubits), weight):
        for states in itertools.product((1, 2, 3), repeat=weight):
            x_mask = 0
            z_mask = 0
            for position, state in zip(support, states, strict=True):
                if state & 1:
                    x_mask |= 1 << position
                if state & 2:
                    z_mask |= 1 << position
            xs[cursor] = x_mask
            zs[cursor] = z_mask
            cursor += 1
    if cursor != total:
        raise AssertionError("ternary pattern enumeration mismatch")
    return xs, zs


def symplectic_column_masks(H: np.ndarray) -> tuple[list[int], list[int]]:
    n = H.shape[1] // 2
    x_contrib = [row_mask(H[:, n + q]) for q in range(n)]
    z_contrib = [row_mask(H[:, q]) for q in range(n)]
    return x_contrib, z_contrib


def logical_at_or_below(
    H: np.ndarray, cap: int, *, budget: float = SYM_PATTERN_BUDGET
) -> dict[str, Any]:
    """Exact MITM search for a non-stabilizer centralizer vector of wt <= cap.

    Returns ``logical_exists_at_or_below_cap`` or ``unresolved`` when the
    per-side ternary pattern budget is exceeded (recorded, never guessed).
    """

    H = np.asarray(H, dtype=np.uint8)
    n = H.shape[1] // 2
    n_left = n // 2
    n_right = n - n_left
    x_contrib, z_contrib = symplectic_column_masks(H)
    left_x = xor_lookup(x_contrib[:n_left])
    left_z = xor_lookup(z_contrib[:n_left])
    right_x = xor_lookup(x_contrib[n_left:])
    right_z = xor_lookup(z_contrib[n_left:])
    stabilizers = E41.Echelon(rows_to_masks(H))

    examined = 0
    discarded = 0
    for total_weight in range(1, cap + 1):
        for left_weight in range(total_weight + 1):
            right_weight = total_weight - left_weight
            if left_weight > n_left or right_weight > n_right:
                continue
            if (
                comb(n_left, left_weight) * 3**left_weight > budget
                or comb(n_right, right_weight) * 3**right_weight > budget
            ):
                return {
                    "logical_exists_at_or_below_cap": None,
                    "unresolved": "pattern budget exceeded",
                    "unresolved_weight": total_weight,
                    "centralizer_pairs_examined": examined,
                }
            left_x_p, left_z_p = ternary_patterns(n_left, left_weight)
            right_x_p, right_z_p = ternary_patterns(n_right, right_weight)
            left_syndrome = (
                left_x[left_x_p.astype(np.intp)] ^ left_z[left_z_p.astype(np.intp)]
            )
            right_syndrome = (
                right_x[right_x_p.astype(np.intp)] ^ right_z[right_z_p.astype(np.intp)]
            )
            order = np.argsort(right_syndrome, kind="stable")
            sorted_right = right_syndrome[order]
            lo = np.searchsorted(sorted_right, left_syndrome, side="left")
            hi = np.searchsorted(sorted_right, left_syndrome, side="right")
            matched_left = np.flatnonzero(hi > lo)
            for li in matched_left:
                for ordered_ri in range(int(lo[li]), int(hi[li])):
                    ri = int(order[ordered_ri])
                    examined += 1
                    x_mask = int(left_x_p[li]) | (int(right_x_p[ri]) << n_left)
                    z_mask = int(left_z_p[li]) | (int(right_z_p[ri]) << n_left)
                    vector = x_mask | (z_mask << n)
                    if stabilizers.contains(vector):
                        discarded += 1
                        continue
                    return {
                        "logical_exists_at_or_below_cap": True,
                        "minimum_logical_weight_within_cap": total_weight,
                        "centralizer_pairs_examined": examined,
                        "stabilizers_discarded": discarded,
                    }
    return {
        "logical_exists_at_or_below_cap": False,
        "minimum_logical_weight_within_cap": None,
        "centralizer_pairs_examined": examined,
        "stabilizers_discarded": discarded,
    }


# ---------------------------------------------------------------------------
# Residual-class perturbation probe (total [C D] support <= 3)
# ---------------------------------------------------------------------------
def singleton_defect_signatures(
    ell: int, m: int, A: np.ndarray, B: np.ndarray
) -> list[int]:
    from qec_research.codes.bicycle import monomial_matrix

    dim = ell * m
    zero = np.zeros((dim, dim), dtype=np.uint8)
    signatures: list[int] = []
    for coordinate in range(2 * dim):
        monomial = monomial_matrix(
            ell, m, (coordinate % dim) // m, (coordinate % dim) % m
        )
        if coordinate < dim:
            defect = commutation_defect(A, B, monomial, zero)
        else:
            defect = commutation_defect(A, B, zero, monomial)
        signatures.append(row_mask(defect.reshape(-1)))
    return signatures


def pbb_check_matrix(HX: np.ndarray, HZ: np.ndarray, CD: np.ndarray) -> np.ndarray:
    n = HX.shape[1]
    return np.vstack(
        (
            np.hstack((HX, CD)),
            np.hstack((np.zeros((HZ.shape[0], n), dtype=np.uint8), HZ)),
        )
    ).astype(np.uint8)


def probe_perturbations(parent: dict[str, Any]) -> dict[str, Any]:
    """Exhaust all commuting total-support<=3 perturbations of one parent."""

    started = time.perf_counter()
    ell, m = int(parent["ell"]), int(parent["m"])
    dim, n = ell * m, 2 * ell * m
    HX = np.asarray(parent["HX"], dtype=np.uint8)
    HZ = np.asarray(parent["HZ"], dtype=np.uint8)
    A, B = HX[:, :dim], HX[:, dim:]
    left_kernel = nullspace_np(np.ascontiguousarray(HX.T))
    sz_rank = int(parent["rank_hz"])
    s_z_rows = rows_to_masks(HZ)
    t_parent = int(parent["T"])
    d_z = int(parent["d_Z"])
    k_parent = int(parent["k_P"])

    signatures = singleton_defect_signatures(ell, m, A, B)
    zero = np.zeros((dim, dim), dtype=np.uint8)

    raw = 0
    commuting = 0
    delta_positive = 0
    relation: Counter[str] = Counter()
    containment: Counter[str] = Counter()
    distance_method: Counter[str] = Counter()
    d_q_relation: Counter[str] = Counter()
    significant: list[dict[str, Any]] = []
    upper_law_breaks: list[dict[str, Any]] = []
    b_prime_breaks: list[dict[str, Any]] = []

    for weight in range(1, PERT_SUPPORT_CAP + 1):
        for combo in itertools.combinations(range(2 * dim), weight):
            raw += 1
            signature = 0
            for coordinate in combo:
                signature ^= signatures[coordinate]
            if signature:
                continue
            commuting += 1
            c_terms, d_terms = support_terms(combo, ell, m)
            C = poly_matrix(ell, m, list(c_terms)) if c_terms else zero
            D = poly_matrix(ell, m, list(d_terms)) if d_terms else zero
            if commutation_defect(A, B, C, D).any():
                raise AssertionError("signature filter accepted a noncommuting PBB")
            CD = np.hstack((C, D)).astype(np.uint8)
            delta_rows = gf2_product(left_kernel, CD)
            bar_basis = E41.Echelon(s_z_rows)
            for row in delta_rows:
                bar_basis.insert(row_mask(row))
            dim_bar = int(bar_basis.rank - sz_rank)
            if dim_bar == 0:
                continue
            delta_positive += 1
            rel = (
                "less_than_T"
                if dim_bar < t_parent
                else "equal_T"
                if dim_bar == t_parent
                else "greater_than_T"
            )
            relation[rel] += 1

            H = pbb_check_matrix(HX, HZ, CD)
            k_q = int(n - matrix_rank(H))
            if k_q != k_parent - dim_bar:
                raise AssertionError("PBB dimension identity failed")

            record: dict[str, Any] = {
                "C": json_terms(c_terms),
                "D": json_terms(d_terms),
                "total_C_D_support": weight,
                "dim_Delta": int(matrix_rank(delta_rows)),
                "dim_Delta_bar": dim_bar,
                "T": t_parent,
                "delta_bar_relation_to_T": rel,
                "k_Q": k_q,
            }

            if dim_bar > t_parent:
                record["violation"] = "upper_law_dim_Delta_bar_exceeds_T"
                upper_law_breaks.append(dict(record))

            survivor = next(
                (row for row in parent["_module_rows"] if not bar_basis.contains(row)),
                None,
            )
            absorbed = survivor is None
            containment["contained" if absorbed else "not_contained"] += 1
            record["M_bar_contained_in_Delta_bar"] = absorbed
            if not absorbed:
                if int(survivor).bit_count() != d_z:
                    raise AssertionError("module survivor is not minimum weight")
                record.update(
                    {
                        "d_Q_gt_d_Z": False,
                        "distance_decision": "surviving_minimum_pure_Z_logical",
                        "witness_weight": d_z,
                    }
                )
                distance_method["surviving_minimum_pure_Z_logical"] += 1
                d_q_relation["not_greater"] += 1
            else:
                threshold = logical_at_or_below(H, d_z)
                if threshold["logical_exists_at_or_below_cap"] is None:
                    record.update(
                        {
                            "d_Q_gt_d_Z": None,
                            "distance_decision": "unresolved_pattern_budget",
                            "threshold_search": threshold,
                        }
                    )
                    distance_method["unresolved_pattern_budget"] += 1
                    d_q_relation["unresolved"] += 1
                else:
                    d_q_gt = not bool(threshold["logical_exists_at_or_below_cap"])
                    record.update(
                        {
                            "d_Q_gt_d_Z": d_q_gt,
                            "distance_decision": "exhaustive_symplectic_MITM",
                            "threshold_search": threshold,
                        }
                    )
                    distance_method["exhaustive_symplectic_MITM"] += 1
                    d_q_relation["greater" if d_q_gt else "not_greater"] += 1
                    if d_q_gt and dim_bar != t_parent:
                        record["violation"] = "b_prime_distance_increase_off_T"
                        b_prime_breaks.append(dict(record))

            if absorbed or record.get("d_Q_gt_d_Z") or "violation" in record:
                significant.append(record)

    return {
        "raw_total_supports_enumerated": raw,
        "commutation_valid_supports": commuting,
        "delta_positive_supports": delta_positive,
        "delta_bar_relation_to_T": dict(sorted(relation.items())),
        "minimum_module_containment": dict(sorted(containment.items())),
        "distance_decision_method": dict(sorted(distance_method.items())),
        "d_Q_vs_d_Z": dict(sorted(d_q_relation.items())),
        "significant_records": significant,
        "upper_law_breaks": upper_law_breaks,
        "b_prime_breaks": b_prime_breaks,
        "wall_time_s": round(time.perf_counter() - started, 6),
    }


# ---------------------------------------------------------------------------
# Self-tests (independence protocol)
# ---------------------------------------------------------------------------
def self_tests(shapes: list[tuple[int, int]]) -> dict[str, Any]:
    rng = np.random.default_rng(20260818)
    out: dict[str, Any] = {}

    # 1. Own orbit packing vs the shared reference implementation.
    for ell, m in shapes:
        L = ell * m
        for _ in range(8):
            v = (rng.random(2 * L) < 0.2).astype(np.uint8)
            z = row_mask(v)
            ref = sorted(row_mask(row) for row in translation_orbit(v, ell, m, blocks=2))
            mine = sorted(E41.orbit_masks(z, ell, m))
            if ref != mine:
                raise AssertionError(f"orbit mismatch on {ell}x{m}")
    out["orbit_matches_reference"] = True

    # 2. Echelon rank over packed ints vs rank_np.
    for ell, m in shapes:
        L = ell * m
        for _ in range(8):
            M = (rng.random((L, 2 * L)) < 0.3).astype(np.uint8)
            if matrix_rank(M) != int(rank_np(M)):
                raise AssertionError(f"rank mismatch on {ell}x{m}")
    out["echelon_rank_matches_rank_np"] = True

    # 3. Column-XOR syndrome vs dense GF(2) product.
    for ell, m in shapes[:4]:
        L = ell * m
        M = (rng.random((L, L)) < 0.3).astype(np.uint8)
        cols = E41.column_masks(M)
        for _ in range(50):
            z = int(rng.integers(0, 1 << L))
            if syndrome_of(z, L, cols, cols) != row_mask(
                gf2_product(M, mask_row(z, L)[:, None])[:, 0]
            ):
                raise AssertionError("syndrome mismatch")
    out["column_syndrome_matches_dense"] = True

    # 4. Symplectic MITM threshold vs full centralizer enumeration (tiny PBBs).
    checked = 0
    for ell, m in ((2, 2), (2, 3), (3, 2)):
        L = ell * m
        for seed in range(4):
            rng2 = np.random.default_rng(1000 * ell + 100 * m + seed)
            H = (rng2.random((2 * L, 4 * L)) < 0.2).astype(np.uint8)
            swapped = np.hstack((H[:, 2 * L :], H[:, : 2 * L]))
            centralizer = nullspace_np(swapped)
            if centralizer.shape[0] > 14:
                continue
            stabilizers = E41.Echelon(rows_to_masks(H))
            best = None
            basis = rows_to_masks(centralizer)
            z = 0
            if not stabilizers.contains(0):
                best = 0
            for i in range(1, 1 << len(basis)):
                g = i ^ (i >> 1)
                prev = (i - 1) ^ ((i - 1) >> 1)
                z ^= basis[(g ^ prev).bit_length() - 1]
                if not stabilizers.contains(z):
                    w = z.bit_count()
                    if best is None or w < best:
                        best = w
            for cap in range(0, (best or 0) + 2):
                got = logical_at_or_below(H, cap)
                if got["logical_exists_at_or_below_cap"] is None:
                    raise AssertionError("unexpected unresolved on tiny instance")
                want = best is not None and best <= cap
                if bool(got["logical_exists_at_or_below_cap"]) != want:
                    raise AssertionError("MITM threshold disagrees with enumeration")
            checked += 1
    if checked == 0:
        raise AssertionError("no symplectic self-test instance exercised")
    out["symplectic_mitm_matches_full_enumeration"] = {"instances": checked}

    # 5. Solver-based kernel enumeration vs the pure EXP-041 MITM join on
    #    real enumerated parents (bit-exact by-weight set equality).
    engine_checks = 0
    for ell, m, family in ((2, 3, (2, 3)), (3, 3, (3, 3)), (3, 4, (3, 3)), (2, 6, (3, 3)), (4, 3, (3, 3))):
        L = ell * m
        wa, wb = family
        sup_a = canonical_supports(ell, m, wa)
        sup_b = sup_a if wa == wb else canonical_supports(ell, m, wb)
        mat_a = [poly_matrix(ell, m, [(p // m, p % m) for p in s]) for s in sup_a]
        mat_b = mat_a if wa == wb else [
            poly_matrix(ell, m, [(p // m, p % m) for p in s]) for s in sup_b
        ]
        monomials = [(p // m, p % m) for p in range(L)]
        taken = 0
        for i, (sa, A) in enumerate(zip(sup_a, mat_a, strict=True)):
            for j in range(i if wa == wb else 0, len(sup_b)):
                if taken >= 2:
                    break
                HX = np.hstack((A, mat_b[j])).astype(np.uint8)
                HZ = np.hstack((mat_b[j].T, A.T)).astype(np.uint8)
                k = 2 * L - matrix_rank(HX) - matrix_rank(HZ)
                if k < 2:
                    continue
                code = bb_stabilizer(
                    BBSpec(
                        ell=ell,
                        m=m,
                        A=[monomials[p] for p in sa],
                        B=[monomials[p] for p in sup_b[j]],
                    )
                )
                if len(code.connected_components()) != 1:
                    continue
                analysis = analyse_parent(HX, HZ, ell, m, k)
                if analysis["d_Z"] is None:
                    continue
                solved = collect_by_weight(HX, range(1, analysis["d_Z"] + 1), True)
                legacy = collect_by_weight(HX, range(1, analysis["d_Z"] + 1), False)
                if solved != legacy:
                    raise AssertionError(
                        f"solver/MITM engine mismatch on {ell}x{m} A={sa} B={sup_b[j]}"
                    )
                engine_checks += 1
                taken += 1
            if taken >= 2:
                break
    if engine_checks == 0:
        raise AssertionError("no engine-equality instance exercised")
    out["solver_engine_matches_mitm"] = {"instances": engine_checks}
    return out


def span_cross_check(
    HX: np.ndarray, by_weight: dict[int, list[int]], d_z: int
) -> dict[str, Any] | None:
    """Full nullspace-span equality check of the MITM set (EXP-041 helper)."""

    return E41.nullspace_cross_check(np.ascontiguousarray(HX), by_weight, d_z)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def lattice_shapes() -> list[tuple[int, int]]:
    return sorted(
        (ell, m)
        for ell in range(MIN_LATTICE_DIM, MAX_L + 1)
        for m in range(MIN_LATTICE_DIM, MAX_L + 1)
        if ell * m <= MAX_L
    )


class LatticeDeadlineExceeded(Exception):
    """Raised when a lattice scan passes its wall-clock deadline."""


def scan_shape_family(
    ell: int,
    m: int,
    family: tuple[int, int],
    catalogue_pairs: set | None,
    deadline_ts: float | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], set]:
    L = ell * m
    wa, wb = family
    lists: dict[int, list[tuple[int, ...]]] = {}
    matrices: dict[int, list[np.ndarray]] = {}
    for w in {wa, wb}:
        sup = canonical_supports(ell, m, w)
        lists[w] = sup
        matrices[w] = [
            poly_matrix(ell, m, [(p // m, p % m) for p in support]) for support in sup
        ]

    pairs: list[tuple[tuple[int, ...], np.ndarray, tuple[int, ...], np.ndarray]] = []
    if wa == wb:
        for i, (sa, A) in enumerate(zip(lists[wa], matrices[wa], strict=True)):
            for sb, B in list(zip(lists[wb], matrices[wb]))[i:]:
                pairs.append((sa, A, sb, B))
    else:
        for sa, A in zip(lists[wa], matrices[wa], strict=True):
            for sb, B in zip(lists[wb], matrices[wb], strict=True):
                pairs.append((sa, A, sb, B))

    monomials = [(p // m, p % m) for p in range(L)]
    k_ge_2 = 0
    connected = 0
    analyzed: list[dict[str, Any]] = []
    residual: list[dict[str, Any]] = []
    unresolved: list[str] = []
    k_hist: Counter[int] = Counter()
    covered_catalogue: set = set()

    for sa, A, sb, B in pairs:
        HX = np.hstack((A, B)).astype(np.uint8)
        HZ = np.hstack((B.T, A.T)).astype(np.uint8)
        rank_x = matrix_rank(HX)
        rank_z = matrix_rank(HZ)
        k_parent = 2 * L - rank_x - rank_z
        k_hist[k_parent] += 1
        if k_parent < 2:
            continue
        k_ge_2 += 1
        a_terms = [monomials[p] for p in sa]
        b_terms = [monomials[p] for p in sb]
        code = bb_stabilizer(BBSpec(ell=ell, m=m, A=list(a_terms), B=list(b_terms)))
        if len(code.connected_components()) != 1:
            continue
        connected += 1
        if deadline_ts is not None and time.perf_counter() > deadline_ts:
            raise LatticeDeadlineExceeded(
                f"{ell}x{m}: wall-clock deadline exceeded after {connected} "
                "connected candidates"
            )
        if catalogue_pairs is not None:
            key = (min(sa, sb), max(sa, sb))
            if key in catalogue_pairs:
                covered_catalogue.add(key)

        analysis = analyse_parent(HX, HZ, ell, m, k_parent)
        if analysis["d_Z"] is None:
            unresolved.append(f"{sa}|{sb}: {analysis['unresolved']}")
            continue
        record = {
            "ell": ell,
            "m": m,
            "family": f"{wa}x{wb}",
            "A_terms": json_terms(a_terms),
            "B_terms": json_terms(b_terms),
            "HX": HX,
            "HZ": HZ,
            "k_P": k_parent,
            "rank_hx": rank_x,
            "rank_hz": analysis["rank_hz"],
            "d_Z": analysis["d_Z"],
            "T": analysis["T"],
            "num_minimum_logicals": analysis["num_minimum_logicals"],
            "num_module_rows": analysis["num_module_rows"],
            "analysis_wall_time_s": analysis["wall_time_s"],
            "_module_rows": analysis["module_rows"],
            "_by_weight": analysis["by_weight"],
        }
        analyzed.append(record)
        if 2 * analysis["T"] < k_parent:
            residual.append(record)

    accounting = {
        "family": f"{wa}x{wb}",
        "canonical_A_supports": len(lists[wa]),
        "canonical_B_supports": len(lists[wb]),
        "support_pairs_enumerated": len(pairs),
        "k_parent_histogram": {str(k): int(v) for k, v in sorted(k_hist.items())},
        "k_ge_2_pairs": k_ge_2,
        "connected_positive_k_pairs": connected,
        "parents_analyzed": len(analyzed),
        "unresolved_parents": unresolved,
        "buckets": {},
        "residual_count": len(residual),
    }
    return accounting, analyzed, residual, covered_catalogue


def run() -> dict[str, Any]:
    started = time.perf_counter()
    shapes = lattice_shapes()
    checks = self_tests(shapes[:6] + [(3, 6), (6, 3)])

    # Catalogue coverage expectation: every L=18 catalogue row's canonical
    # parent pair must be rediscovered by the exhaustive (3,3) sweep.
    catalogue_rows = E27.load_catalogue()
    catalogue_pairs: dict[tuple[int, int], set] = {}
    for row in catalogue_rows:
        ell, m = int(row["ell"]), int(row["m"])
        if ell * m != 18:
            continue
        key = (
            canonical_rep(ell, m, [tuple(t) for t in row["A_terms"]]),
            canonical_rep(ell, m, [tuple(t) for t in row["B_terms"]]),
        )
        catalogue_pairs.setdefault((ell, m), set()).add((min(key), max(key)))

    lattice_records: list[dict[str, Any]] = []
    all_residual: list[dict[str, Any]] = []
    totals: Counter[str] = Counter()
    span_checks: list[dict[str, Any]] = []
    fingerprints: set[str] = set()
    catalogue_hits: dict[tuple[int, int], set] = {}
    catalogue_misses: list[str] = []

    for ell, m in shapes:
        L = ell * m
        family_records = []
        for family in SUPPORT_FAMILIES:
            t0 = time.perf_counter()
            accounting, analyzed, residual, covered = scan_shape_family(
                ell, m, family, catalogue_pairs.get((ell, m))
            )
            wall = time.perf_counter() - t0
            buckets: Counter[str] = Counter()
            for rec in analyzed:
                k, t = rec["k_P"], rec["T"]
                if 2 * t > k:
                    buckets["T_gt_half_k"] += 1
                elif 2 * t == k:
                    buckets["T_eq_half_k"] += 1
                else:
                    buckets["T_lt_half_k"] += 1
                if t == k:
                    buckets["T_eq_k"] += 1
                fingerprints.add(E27.matrix_fingerprint(rec["HX"], rec["HZ"]))
            accounting["buckets"] = {k2: int(v) for k2, v in sorted(buckets.items())}
            accounting["wall_time_s"] = round(wall, 3)
            family_records.append(accounting)
            totals["parents_analyzed"] += len(analyzed)
            for key in ("T_gt_half_k", "T_eq_half_k", "T_lt_half_k", "T_eq_k"):
                totals[key] += int(buckets.get(key, 0))
            for rec in residual:
                rec = dict(rec)
                rec["id"] = f"{ell}x{m}_{family[0]}x{family[1]}_{len(all_residual):03d}"
                all_residual.append(rec)
            # Full-nullspace cross-checks of the MITM set (bounded sample).
            for rec in analyzed:
                if len(span_checks) >= SPAN_CROSS_CHECK_PARENTS:
                    break
                ker = nullspace_np(np.ascontiguousarray(rec["HX"]))
                if int(ker.shape[0]) > SPAN_CROSS_CHECK_MAX_KER:
                    continue
                check = span_cross_check(rec["HX"], rec["_by_weight"], rec["d_Z"])
                if check is None:
                    continue
                if not check["set_equal"]:
                    raise AssertionError(
                        f"MITM set disagrees with nullspace span on {ell}x{m} "
                        f"A={rec['A_terms']} B={rec['B_terms']}"
                    )
                span_checks.append(
                    {
                        "lattice": f"{ell}x{m}",
                        "A_terms": rec["A_terms"],
                        "B_terms": rec["B_terms"],
                        **{k2: v for k2, v in check.items()},
                    }
                )
                break  # one per (shape, family) at most
            if catalogue_pairs.get((ell, m)):
                got = catalogue_hits.setdefault((ell, m), set())
                got |= covered
            print(
                f"[exp042] {ell}x{m} family {family[0]}x{family[1]}: "
                f"{accounting['support_pairs_enumerated']} pairs, "
                f"{accounting['parents_analyzed']} analyzed, "
                f"residual={accounting['residual_count']} ({wall:.1f}s)",
                flush=True,
            )

        lattice_records.append(
            {
                "lattice": f"{ell}x{m}",
                "ell": ell,
                "m": m,
                "ell_times_m": L,
                "families": family_records,
            }
        )

    for (ell, m), want in catalogue_pairs.items():
        got = catalogue_hits.get((ell, m), set())
        missing = want - got
        for pair in sorted(missing):
            catalogue_misses.append(f"{ell}x{m}: catalogue pair {pair} not recovered")

    # Residual-class perturbation probe.
    residual_outcomes: list[dict[str, Any]] = []
    upper_law_breaks: list[dict[str, Any]] = []
    b_prime_breaks: list[dict[str, Any]] = []
    for rec in all_residual:
        t0 = time.perf_counter()
        summary = probe_perturbations(rec)
        wall = time.perf_counter() - t0
        outcome = {
            "id": rec["id"],
            "ell": rec["ell"],
            "m": rec["m"],
            "family": rec["family"],
            "A_terms": rec["A_terms"],
            "B_terms": rec["B_terms"],
            "k_P": rec["k_P"],
            "d_Z": rec["d_Z"],
            "T": rec["T"],
            "num_minimum_logicals": rec["num_minimum_logicals"],
            "perturbations": summary,
            "probe_wall_time_s": round(wall, 3),
        }
        residual_outcomes.append(outcome)
        for brk in summary["upper_law_breaks"]:
            upper_law_breaks.append({"parent_id": rec["id"], **brk})
        for brk in summary["b_prime_breaks"]:
            b_prime_breaks.append({"parent_id": rec["id"], **brk})
        print(
            f"[exp042] residual {rec['id']}: k_P={rec['k_P']} T={rec['T']} "
            f"d_Z={rec['d_Z']}; {summary['delta_positive_supports']} delta>0 "
            f"perturbations, breaks={len(summary['upper_law_breaks']) + len(summary['b_prime_breaks'])}",
            flush=True,
        )

    if not all_residual:
        verdict = "SMALL_CLASS_EMPTY"
        verdict_detail = (
            f"No parent with T < k_P/2 exists among the {totals['parents_analyzed']} "
            f"exhaustively enumerated connected k>=2 parents on any lattice with "
            f"ell*m <= {MAX_L} (families "
            + ", ".join(f"{a}x{b}" for a, b in SUPPORT_FAMILIES)
            + "); Conjecture B-prime's residual class is vacuous at small lattices."
        )
    elif upper_law_breaks or b_prime_breaks:
        verdict = "B_PRIME_VIOLATED"
        verdict_detail = (
            f"{len(all_residual)} residual parents found; "
            f"{len(upper_law_breaks)} upper-law breaks and "
            f"{len(b_prime_breaks)} B-prime breaks attached as counterexamples."
        )
    else:
        verdict = "SMALL_CLASS_TESTED"
        unresolved_decisions = sum(
            outcome["perturbations"]["distance_decision_method"].get(
                "unresolved_pattern_budget", 0
            )
            for outcome in residual_outcomes
        )
        verdict_detail = (
            f"{len(all_residual)} residual parents found and exhaustively probed "
            f"over all delta>0 perturbations of total [C D] support <= "
            f"{PERT_SUPPORT_CAP}: zero upper-law breaks, zero B-prime breaks"
            + (
                f" ({unresolved_decisions} distance decisions unresolved by budget)."
                if unresolved_decisions
                else "."
            )
        )

    if catalogue_misses:
        raise AssertionError(
            "exhaustive sweep failed to recover catalogue parents: "
            + "; ".join(catalogue_misses[:5])
        )

    payload = {
        "schema": SCHEMA,
        "experiment": "EXP-042",
        "utc": utc_now(),
        "protocol": {
            "lattice_set": (
                f"all ordered (ell, m) with ell, m >= {MIN_LATTICE_DIM} and "
                f"ell*m <= {MAX_L} ({len(shapes)} shapes)"
            ),
            "support_families": [f"{a}x{b}" for a, b in SUPPORT_FAMILIES],
            "parent_equivalence": (
                "supports identified under independent Z_ell x Z_m translations of "
                "the A and B blocks; {A, B} unordered; k_P, d_Z, T are invariants "
                "of the class"
            ),
            "nondegenerate_predicate": "k_P >= 2 and a single Tanner component",
            "d_Z_method": (
                "ascending-weight exact enumeration of ker[A B]: per split the "
                "cheaper side is enumerated and the other side solved over a "
                "coset of its kernel (GF(2) column solver); one-sided splits "
                "read the kernel spans directly; the EXP-041 MITM join is the "
                "fallback when neither solve orientation is feasible.  Both "
                "engines verified bit-exact on real parents (self_tests) and "
                "against full nullspace spans (span_cross_checks)"
            ),
            "T_method": (
                "translation-orbit span of every minimum-weight Z-logical modulo "
                "S_Z; every translate verified to lie in ker[A B]"
            ),
            "perturbation_class": (
                "every literal [C D] support of total weight 1..3 with zero "
                "AC^T + BD^T commutation defect"
            ),
            "d_Q_threshold_method": (
                "surviving minimum pure-Z logical, else exhaustive symplectic "
                "MITM over all ternary Pauli patterns of weight <= d_Z under "
                f"per-side budget {int(SYM_PATTERN_BUDGET)}"
            ),
            "d_weight_cap": D_WEIGHT_CAP,
            "match_guard": MATCH_GUARD,
        },
        "known_exception_context": {
            "fingerprint": "9a7638586033f4c7ae22623a4b53173fa72492f1b1212e94a47e9d8aa2467a4c",
            "lattice": "12x6",
            "ell_times_m": 72,
            "k_parent": 12,
            "T": 4,
            "d_z_parent": 6,
            "note": "outside the ell*m <= 24 scan range (certificate read-only)",
        },
        "self_tests": checks,
        "span_cross_checks": span_checks,
        "catalogue_coverage": {
            "catalogue_l18_parent_pairs": {
                f"{ell}x{m}": len(pairs) for (ell, m), pairs in sorted(catalogue_pairs.items())
            },
            "recovered_pairs": {
                f"{ell}x{m}": len(catalogue_hits.get((ell, m), set()))
                for (ell, m) in sorted(catalogue_pairs)
            },
            "all_recovered": not catalogue_misses,
        },
        "aggregate": {
            "shapes": len(shapes),
            "parents_analyzed": int(totals["parents_analyzed"]),
            "bucket_T_gt_half_k": int(totals["T_gt_half_k"]),
            "bucket_T_eq_half_k": int(totals["T_eq_half_k"]),
            "bucket_T_lt_half_k": int(totals["T_lt_half_k"]),
            "bucket_T_eq_k": int(totals["T_eq_k"]),
            "distinct_parent_fingerprints": len(fingerprints),
            "residual_parents": len(all_residual),
            "residual_delta_positive_perturbations": sum(
                o["perturbations"]["delta_positive_supports"] for o in residual_outcomes
            ),
            "upper_law_breaks": len(upper_law_breaks),
            "b_prime_breaks": len(b_prime_breaks),
        },
        "lattices": lattice_records,
        "residual_parents": residual_outcomes,
        "counterexamples": {
            "upper_law_dim_Delta_bar_exceeds_T": upper_law_breaks,
            "b_prime_distance_increase_off_T": b_prime_breaks,
        },
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        "sat_solver_used": False,
        "pysat_imported": "pysat" in sys.modules,
        "sat_decide_imported": "qec_research.distance.sat_decide" in sys.modules,
        "stim_imported": "stim" in sys.modules,
        "wall_time_s": round(time.perf_counter() - started, 3),
    }
    atomic_write_json(OUT, payload)
    return payload


def main() -> int:
    payload = run()
    print(f"[exp042] wrote {OUT}")
    print(f"[exp042] verdict: {payload['verdict']} -- {payload['verdict_detail']}")
    print(
        f"[exp042] aggregate: {payload['aggregate']} "
        f"({payload['wall_time_s']:.1f}s)"
    )
    if payload["pysat_imported"] or payload["sat_decide_imported"]:
        print("[exp042] ERROR: SAT modules were imported")
        return 1
    return 0 if payload["verdict"] != "B_PRIME_VIOLATED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
