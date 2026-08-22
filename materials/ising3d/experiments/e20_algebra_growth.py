"""Exact depth-graded growth of the two-generator Ising transfer algebra.

For a layer graph ``Gamma`` this script studies

    A = sum_v X_v,              B = sum_{(u,v) in E(Gamma)} Z_u Z_v

and ``D_k``, the dimension of all left-normed brackets of ``A`` and ``B`` of
length at most ``k``.  Arithmetic is exact over a 31-bit prime.  Every
commutator is divided by two; two is invertible modulo the prime, so this
changes no span or rank.

Two exact engines are used.  The dense engine first quotients Pauli support by
all point-group symmetries and stores echelon rows as int32.  Products are
formed in int64; ``PRIME_31**2 < 2**63`` makes this overflow-safe.  The sparse
engine discovers only the orbit sums reached through the requested depth and
is used when the complete support would itself be too large.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp
import numpy as np


PRIME_31 = 2_147_483_647
SECOND_PRIME_31 = 2_147_483_629
SCRIPT = "experiments/e20_algebra_growth.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "algebra_growth"
_MERSENNE_SHIFT = 31


class SupportLimit(RuntimeError):
    """Raised when complete orbit-support enumeration exceeds its explicit cap."""

    def __init__(self, count: int, cap: int) -> None:
        super().__init__(f"orbit support exceeded cap={cap:,} at count={count:,}")
        self.count = count
        self.cap = cap


def _compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Return the permutation ``left after right`` (old-site to new-site form)."""

    return tuple(left[right[i]] for i in range(len(left)))


def _permutation_group(generators: Iterable[tuple[int, ...]], n: int) -> tuple[tuple[int, ...], ...]:
    identity = tuple(range(n))
    group = [identity]
    seen = {identity}
    head = 0
    generators = tuple(generators)
    while head < len(group):
        current = group[head]
        head += 1
        for generator in generators:
            product = _compose(generator, current)
            if product not in seen:
                seen.add(product)
                group.append(product)
    return tuple(group)


def grid_bonds(rows: int, cols: int, periodic_cols: bool = False) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Return the open rectangular layer, optionally periodic in its column direction."""

    if rows < 1 or cols < 2:
        raise ValueError("rows must be positive and cols must be at least two")
    n = rows * cols

    def site(x: int, y: int) -> int:
        return x * cols + y

    bonds: list[tuple[int, int]] = []
    for x in range(rows):
        for y in range(cols):
            if x + 1 < rows:
                bonds.append((site(x, y), site(x + 1, y)))
            if y + 1 < cols:
                bonds.append((site(x, y), site(x, y + 1)))
            elif periodic_cols:
                if cols == 2:
                    # The repository convention keeps both parallel bonds.  The
                    # duplicate only rescales B over an odd field, but retaining
                    # it keeps the lattice definition exact.
                    bonds.append((site(x, y), site(x, 0)))
                elif cols > 2:
                    bonds.append((site(x, y), site(x, 0)))
    return n, tuple(bonds)


def _point_group(rows: int, cols: int, periodic_cols: bool) -> tuple[tuple[int, ...], ...]:
    n = rows * cols

    def site(x: int, y: int) -> int:
        return x * cols + y

    generators = [
        tuple(site(rows - 1 - x, y) for x in range(rows) for y in range(cols)),
        tuple(site(x, cols - 1 - y) for x in range(rows) for y in range(cols)),
    ]
    if rows == cols and not periodic_cols:
        generators.append(tuple(site(y, x) for x in range(rows) for y in range(cols)))
    if periodic_cols:
        generators.append(tuple(site(x, (y + 1) % cols) for x in range(rows) for y in range(cols)))
    return _permutation_group(generators, n)


def _permute_pauli(value: int, permutation: Sequence[int], n: int) -> int:
    mask = (1 << n) - 1
    a = value & mask
    b = (value >> n) & mask
    new_a = 0
    new_b = 0
    for old, new in enumerate(permutation):
        if (a >> old) & 1:
            new_a |= 1 << new
        if (b >> old) & 1:
            new_b |= 1 << new
    return new_a | (new_b << n)


def _commutator_sign_half(generator: int, value: int, n: int) -> int:
    """Coefficient of Q_(g xor v) in [Q_g,Q_v]/2: zero or +/-1."""

    mask = (1 << n) - 1
    first = ((((generator >> n) & mask) & (value & mask)).bit_count()) & 1
    second = ((((value >> n) & mask) & (generator & mask)).bit_count()) & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1

def _pauli_grade(value: int, n: int) -> int:
    """Joint (X-count parity, transpose parity) grading, encoded in two bits."""

    mask = (1 << n) - 1
    a = value & mask
    b = (value >> n) & mask
    return ((a.bit_count() & 1) << 1) | ((a & b).bit_count() & 1)


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _physical_memory() -> int | None:
    try:
        if sys.platform == "darwin":
            return int(os.popen("sysctl -n hw.memsize").read().strip())
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        return int(pages * page_size)
    except (OSError, ValueError):
        return None


def _case_name(rows: int, cols: int, periodic_cols: bool) -> str:
    if rows == 1:
        return f"ring_{cols}" if periodic_cols else f"chain_open_{cols}"
    return f"grid_{rows}x{cols}"


@dataclass
class OrbitModel:
    rows: int
    cols: int
    periodic_cols: bool = False

    def __post_init__(self) -> None:
        self.n, self.bonds = grid_bonds(self.rows, self.cols, self.periodic_cols)
        self.group = _point_group(self.rows, self.cols, self.periodic_cols)
        self.x_terms = tuple(1 << i for i in range(self.n))
        self.zz_terms = tuple((1 << (self.n + i)) | (1 << (self.n + j)) for i, j in self.bonds)
        self.single_terms = self.x_terms + self.zz_terms
        self._canonical_cache: dict[int, int] = {}

    def canonical(self, value: int) -> int:
        cached = self._canonical_cache.get(value)
        if cached is not None:
            return cached
        result = min(_permute_pauli(value, permutation, self.n) for permutation in self.group)
        self._canonical_cache[value] = result
        return result

    def seed(self, terms: Sequence[int]) -> dict[int, int]:
        # Orbit coordinates are coefficients per Pauli string, not sums over
        # orbit members.  Parallel bonds do contribute their true multiplicity.
        multiplicities: dict[int, int] = {}
        for term in terms:
            multiplicities[term] = multiplicities.get(term, 0) + 1
        result: dict[int, int] = {}
        for term, multiplicity in multiplicities.items():
            orbit = self.canonical(term)
            previous = result.setdefault(orbit, multiplicity)
            if previous != multiplicity:
                raise AssertionError("point-group orbit has nonuniform generator coefficients")
        return result

    def enumerate_support(self, cap: int) -> tuple[int, ...]:
        """Exact DLA support, already quotiented into point-group orbits."""

        seen = {self.canonical(term) for term in self.single_terms}
        queue = sorted(seen)
        head = 0
        while head < len(queue):
            value = queue[head]
            head += 1
            for generator in self.single_terms:
                if _commutator_sign_half(generator, value, self.n):
                    target = self.canonical(value ^ generator)
                    if target not in seen:
                        seen.add(target)
                        queue.append(target)
                        if len(seen) > cap:
                            raise SupportLimit(len(seen), cap)
        mask = (1 << self.n) - 1
        if any((((value >> self.n) & mask).bit_count() & 1) for value in seen):
            raise AssertionError("global spin-flip Z2 support restriction was violated")
        # Reverse BFS order makes shallow-support pivots land at the right of
        # dense rows, sharply reducing exact elimination slice lengths.
        return tuple(reversed(queue))


class DenseOrbitEngine:
    """Dense int32 echelon engine split by two exact Z2 gradings."""

    def __init__(
        self,
        model: OrbitModel,
        prime: int,
        support: tuple[int, ...],
        memory_limit_bytes: int | None,
        max_dimension: int | None,
    ) -> None:
        if prime >= 1 << 31 or prime * prime >= 1 << 63:
            raise ValueError("prime must fit int32 and have prime**2 < 2**63")
        self.model = model
        self.prime = prime
        self.support = support
        self.m = len(support)
        self.max_dimension = max_dimension
        self.blocks = tuple(
            tuple(value for value in support if _pauli_grade(value, model.n) == block)
            for block in range(4)
        )
        self.block_sizes = tuple(len(block) for block in self.blocks)
        self.index = tuple(
            {value: local_index for local_index, value in enumerate(block)}
            for block in self.blocks
        )

        full_bytes = sum(4 * size * size for size in self.block_sizes)
        if memory_limit_bytes is None or full_bytes <= memory_limit_bytes:
            scale = 1.0
        else:
            scale = memory_limit_bytes / full_bytes
        capacities = [
            min(size, max(1, int(size * scale)))
            for size in self.block_sizes
        ]
        if sum(4 * capacity * size for capacity, size in zip(capacities, self.block_sizes)) > (
            memory_limit_bytes if memory_limit_bytes is not None else full_bytes
        ):
            raise MemoryError("memory limit cannot hold the two seed rows")
        self.capacities = tuple(capacities)
        self.rows = tuple(
            np.empty((capacity, size), dtype=np.int32)
            for capacity, size in zip(self.capacities, self.block_sizes)
        )
        self.pivot_row = tuple(np.full(size, -1, dtype=np.int32) for size in self.block_sizes)
        self.pivot_inverse = tuple(np.zeros(size, dtype=np.int32) for size in self.block_sizes)
        self.row_counts = [0, 0, 0, 0]
        self.actions = tuple(
            tuple(self._action_table(terms, block, action_index) for block in range(4))
            for action_index, terms in enumerate((model.x_terms, model.zz_terms))
        )

    @property
    def row_count(self) -> int:
        return sum(self.row_counts)

    @property
    def allocated_basis_bytes(self) -> int:
        return sum(int(rows.nbytes) for rows in self.rows)

    @property
    def used_basis_bytes(self) -> int:
        return sum(
            count * size * np.dtype(np.int32).itemsize
            for count, size in zip(self.row_counts, self.block_sizes)
        )

    def _action_table(
        self,
        terms: Sequence[int],
        source_block: int,
        action_index: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        # ad_A toggles both X-count parity and transpose parity; ad_B
        # preserves X-count parity and toggles transpose parity.
        target_block = source_block ^ (3 if action_index == 0 else 1)
        targets = self.blocks[target_block]
        sources = np.full((len(terms), len(targets)), -1, dtype=np.int32)
        signs = np.zeros((len(terms), len(targets)), dtype=np.int8)
        for term_index, generator in enumerate(terms):
            for target_index, target in enumerate(targets):
                source = target ^ generator
                sign = _commutator_sign_half(generator, source, self.model.n)
                if sign:
                    canonical_source = self.model.canonical(source)
                    source_index = self.index[source_block].get(canonical_source)
                    if source_index is None:
                        raise AssertionError("support or grading block is not closed under an action")
                    sources[term_index, target_index] = source_index
                    signs[term_index, target_index] = sign
        return sources, signs

    def dense_seed(self, terms: Sequence[int]) -> tuple[int, np.ndarray]:
        seed = self.model.seed(terms)
        grades = {_pauli_grade(value, self.model.n) for value in seed}
        if len(grades) != 1:
            raise AssertionError("a seed crossed exact grading blocks")
        block = grades.pop()
        result = np.zeros(self.block_sizes[block], dtype=np.int64)
        for orbit, coefficient in seed.items():
            result[self.index[block][orbit]] = coefficient
        return block, result

    def apply(self, action_index: int, row_reference: tuple[int, int]) -> tuple[int, np.ndarray]:
        source_block, row_index = row_reference
        target_block = source_block ^ (3 if action_index == 0 else 1)
        sources, signs = self.actions[action_index][source_block]
        row = self.rows[source_block][row_index]
        result = np.zeros(self.block_sizes[target_block], dtype=np.int64)
        for source, sign in zip(sources, signs):
            selected = source >= 0
            result[selected] += sign[selected].astype(np.int64) * row[source[selected]].astype(np.int64)
        result %= self.prime
        return target_block, result

    def _subtract_multiple(
        self,
        block: int,
        target: np.ndarray,
        start: int,
        factor: int,
        row_index: int,
    ) -> None:
        basis = self.rows[block][row_index, start:].astype(np.int64)
        if self.prime == PRIME_31:
            product = basis * factor
            product = (product & PRIME_31) + (product >> _MERSENNE_SHIFT)
            product = (product & PRIME_31) + (product >> _MERSENNE_SHIFT)
            product[product >= PRIME_31] -= PRIME_31
            view = target[start:]
            view -= product
            view[view < 0] += PRIME_31
        else:
            target[start:] = (target[start:] - factor * basis) % self.prime

    def reduce_add(self, vector: tuple[int, np.ndarray]) -> tuple[int, int] | None:
        block, raw_value = vector
        value = np.asarray(raw_value, dtype=np.int64) % self.prime
        while True:
            nonzero = np.flatnonzero(value)
            if not nonzero.size:
                return None
            lead = int(nonzero[0])
            existing = int(self.pivot_row[block][lead])
            if existing < 0:
                if self.max_dimension is not None and self.row_count >= self.max_dimension:
                    raise MemoryError(f"dense basis exceeded max_dimension={self.max_dimension:,}")
                if self.row_counts[block] >= self.capacities[block]:
                    raise MemoryError(
                        f"grading block {block} exhausted int32 capacity "
                        f"{self.capacities[block]:,} rows x {self.block_sizes[block]:,} columns; "
                        f"{self.allocated_basis_bytes:,} total allocated bytes"
                    )
                row_index = self.row_counts[block]
                self.rows[block][row_index] = value.astype(np.int32)
                self.pivot_row[block][lead] = row_index
                self.pivot_inverse[block][lead] = pow(
                    int(value[lead]), self.prime - 2, self.prime
                )
                self.row_counts[block] += 1
                return block, row_index
            factor = (
                int(value[lead]) * int(self.pivot_inverse[block][lead])
            ) % self.prime
            self._subtract_multiple(block, value, lead, factor, existing)


class SparseOrbitEngine:
    """Sparse exact engine that discovers Pauli orbit sums only as depths require them."""

    # Python dict storage is implementation-dependent.  This is deliberately a
    # conservative accounting constant used only for an explicit safety wall;
    # ru_maxrss is also recorded as the observed wall.
    ESTIMATED_BYTES_PER_ENTRY = 96

    def __init__(
        self,
        model: OrbitModel,
        prime: int,
        memory_limit_bytes: int | None,
        max_dimension: int | None,
    ) -> None:
        self.model = model
        self.prime = prime
        self.memory_limit_bytes = memory_limit_bytes
        self.max_dimension = max_dimension
        self.basis: dict[int, dict[int, int]] = {}
        self.inverse: dict[int, int] = {}
        self.rows: list[dict[int, int]] = []
        self.nnz = 0
        self.transition_cache: dict[tuple[int, int], dict[int, int]] = {}
        self.discovered_orbits: set[int] = set()

    @property
    def row_count(self) -> int:
        return len(self.rows)

    @property
    def used_basis_bytes(self) -> int:
        return self.nnz * self.ESTIMATED_BYTES_PER_ENTRY

    @property
    def allocated_basis_bytes(self) -> int:
        return self.used_basis_bytes

    def _transitions(self, source: int, action_index: int) -> dict[int, int]:
        key = (source, action_index)
        cached = self.transition_cache.get(key)
        if cached is not None:
            return cached
        terms = self.model.x_terms if action_index == 0 else self.model.zz_terms
        targets = {
            self.model.canonical(source ^ generator)
            for generator in terms
            if _commutator_sign_half(generator, source, self.model.n)
        }
        result: dict[int, int] = {}
        for target in targets:
            coefficient = 0
            for generator in terms:
                preimage = target ^ generator
                sign = _commutator_sign_half(generator, preimage, self.model.n)
                if sign and self.model.canonical(preimage) == source:
                    coefficient += sign
            if coefficient:
                result[target] = coefficient
        self.discovered_orbits.add(source)
        self.discovered_orbits.update(result)
        self.transition_cache[key] = result
        return result

    def dense_seed(self, terms: Sequence[int]) -> dict[int, int]:
        seed = self.model.seed(terms)
        self.discovered_orbits.update(seed)
        return seed

    def apply(self, action_index: int, row_index: int) -> dict[int, int]:
        result: dict[int, int] = {}
        for source, value in self.rows[row_index].items():
            for target, coefficient in self._transitions(source, action_index).items():
                new_value = (result.get(target, 0) + value * coefficient) % self.prime
                if new_value:
                    result[target] = new_value
                else:
                    result.pop(target, None)
        return result

    def reduce_add(self, vector: dict[int, int]) -> int | None:
        value = {key: coefficient % self.prime for key, coefficient in vector.items() if coefficient % self.prime}
        while value:
            lead = min(value)
            existing = self.basis.get(lead)
            if existing is None:
                if self.max_dimension is not None and self.row_count >= self.max_dimension:
                    raise MemoryError(f"sparse basis exceeded max_dimension={self.max_dimension:,}")
                projected_nnz = self.nnz + len(value)
                projected_bytes = projected_nnz * self.ESTIMATED_BYTES_PER_ENTRY
                if self.memory_limit_bytes is not None and projected_bytes > self.memory_limit_bytes:
                    raise MemoryError(
                        f"estimated sparse basis {projected_nnz:,} entries x "
                        f"{self.ESTIMATED_BYTES_PER_ENTRY} bytes = {projected_bytes:,} bytes "
                        f"exceeds limit={self.memory_limit_bytes:,}"
                    )
                row_index = self.row_count
                self.rows.append(value)
                self.basis[lead] = value
                self.inverse[lead] = pow(value[lead], self.prime - 2, self.prime)
                self.nnz = projected_nnz
                return row_index
            factor = (value[lead] * self.inverse[lead]) % self.prime
            for key, coefficient in existing.items():
                new_value = (value.get(key, 0) - factor * coefficient) % self.prime
                if new_value:
                    value[key] = new_value
                else:
                    value.pop(key, None)
        return None


def compute_support_statistics(
    rows: int,
    cols: int,
    *,
    periodic_cols: bool = False,
    cap: int = 2_000_000,
) -> dict[str, object]:
    """Enumerate the exact point-group orbit column space without elimination."""

    started = time.monotonic()
    rss_started = _rss_bytes()
    model = OrbitModel(rows, cols, periodic_cols)
    try:
        support = model.enumerate_support(cap)
    except SupportLimit as error:
        return {
            "case": _case_name(rows, cols, periodic_cols),
            "complete": False,
            "count_lower_bound": error.count,
            "cap": cap,
            "stop_reason": str(error),
            "elapsed_seconds": f"{time.monotonic() - started:.6f}",
            "rss_start_bytes": rss_started,
            "rss_peak_bytes": _rss_bytes(),
        }
    block_sizes = [0, 0, 0, 0]
    for value in support:
        block_sizes[_pauli_grade(value, model.n)] += 1
    return {
        "case": _case_name(rows, cols, periodic_cols),
        "complete": True,
        "support_orbits": len(support),
        "point_group_order": len(model.group),
        "grading_block_columns": block_sizes,
        "unblocked_int32_square_bytes": 4 * len(support) ** 2,
        "graded_int32_square_bytes": 4 * sum(size * size for size in block_sizes),
        "elapsed_seconds": f"{time.monotonic() - started:.6f}",
        "rss_start_bytes": rss_started,
        "rss_peak_bytes": _rss_bytes(),
    }


def _growth_estimate(dimensions: Sequence[int]) -> dict[str, object] | None:
    """Descriptive log-linear and log-log fits; never used as a proof."""

    rising = list(dimensions)
    if len(rising) >= 2 and rising[-1] == rising[-2]:
        rising.pop()
    if len(rising) < 4:
        return None
    end = len(rising)
    taper_depth = None
    for index in range(5, len(rising)):
        if 10 * rising[index] < 11 * rising[index - 1]:
            end = index
            taper_depth = index + 1
            break
    start = max(1, end - 8)
    depths = list(range(start + 1, end + 1))
    values = rising[start:end]
    mp.mp.dps = 50

    def regression(xs: Sequence[mp.mpf], ys: Sequence[mp.mpf]) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        count = mp.mpf(len(xs))
        mean_x = sum(xs) / count
        mean_y = sum(ys) / count
        covariance = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
        variance = sum((x - mean_x) ** 2 for x in xs)
        slope = covariance / variance
        intercept = mean_y - slope * mean_x
        residual = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
        total = sum((y - mean_y) ** 2 for y in ys)
        r_squared = mp.mpf(1) if total == 0 else 1 - residual / total
        return slope, intercept, r_squared

    log_values = [mp.log(value) for value in values]
    exp_slope, _, exp_r2 = regression([mp.mpf(depth) for depth in depths], log_values)
    power_slope, _, power_r2 = regression([mp.log(depth) for depth in depths], log_values)
    return {
        "window_depths": [depths[0], depths[-1]],
        "window_rule": "last up to eight points before the first post-depth-5 ratio below 1.1",
        "excluded_finite_size_taper_from_depth": taper_depth,
        "points": len(depths),
        "exponential_base_per_depth": mp.nstr(mp.e ** exp_slope, 14),
        "exponential_log_linear_r_squared": mp.nstr(exp_r2, 14),
        "power_law_exponent": mp.nstr(power_slope, 14),
        "power_log_log_r_squared": mp.nstr(power_r2, 14),
        "interpretation": "descriptive finite-depth fit only; not an asymptotic bound",
        "mpmath_dps": 50,
    }


def compute_profile(
    rows: int,
    cols: int,
    *,
    periodic_cols: bool = False,
    max_depth: int | None = None,
    prime: int = PRIME_31,
    engine: str = "auto",
    dense_support_cap: int = 100_000,
    memory_limit_bytes: int | None = None,
    max_dimension: int | None = None,
    max_seconds: float | None = None,
) -> dict[str, object]:
    """Compute exact ``D_k`` values through saturation or an explicit resource wall.

    ``dimensions[k-1]`` is exactly ``D_k`` for every completed depth.  If a
    resource wall is hit part-way through the next depth, ``lower_bound`` also
    counts the independent rows already certified at that incomplete depth.
    """

    if max_depth is not None and max_depth < 1:
        raise ValueError("max_depth must be positive or None")
    if engine not in {"auto", "dense", "sparse"}:
        raise ValueError("engine must be 'auto', 'dense', or 'sparse'")
    model = OrbitModel(rows, cols, periodic_cols)
    started = time.monotonic()
    rss_started = _rss_bytes()
    support: tuple[int, ...] | None = None
    support_wall: str | None = None
    chosen_engine = engine
    if engine in {"auto", "dense"}:
        try:
            support = model.enumerate_support(dense_support_cap)
            chosen_engine = "dense"
        except SupportLimit as error:
            support_wall = str(error)
            if engine == "dense":
                raise
            chosen_engine = "sparse"

    if chosen_engine == "dense":
        assert support is not None
        worker: DenseOrbitEngine | SparseOrbitEngine = DenseOrbitEngine(
            model, prime, support, memory_limit_bytes, max_dimension
        )
        support_orbits: int | None = len(support)
        support_complete = True
    else:
        worker = SparseOrbitEngine(model, prime, memory_limit_bytes, max_dimension)
        support_orbits = None
        support_complete = False

    frontier: list[object] = []
    for seed_terms in (model.x_terms, model.zz_terms):
        added = worker.reduce_add(worker.dense_seed(seed_terms))
        if added is not None:
            frontier.append(added)
    if worker.row_count != 2:
        raise AssertionError("A and B must be linearly independent")

    dimensions = [2]
    completed_depth = 1
    saturated = False
    stop_reason = ""
    incomplete_depth: int | None = None
    try:
        while True:
            if max_depth is not None and completed_depth >= max_depth:
                stop_reason = "max_depth"
                break
            if max_seconds is not None and time.monotonic() - started >= max_seconds:
                stop_reason = f"TimeoutError: elapsed time exceeded max_seconds={max_seconds} at a depth boundary"
                break
            next_frontier: list[object] = []
            next_depth = completed_depth + 1
            for row_index in frontier:
                for action_index in (0, 1):
                    if max_seconds is not None and time.monotonic() - started >= max_seconds:
                        incomplete_depth = next_depth
                        raise TimeoutError(f"elapsed time exceeded max_seconds={max_seconds}")
                    candidate = worker.apply(action_index, row_index)
                    added = worker.reduce_add(candidate)
                    if added is not None:
                        next_frontier.append(added)
            completed_depth = next_depth
            dimensions.append(worker.row_count)
            if not next_frontier:
                saturated = True
                stop_reason = "saturated"
                break
            frontier = next_frontier
    except (MemoryError, TimeoutError) as error:
        stop_reason = f"{type(error).__name__}: {error}"
        if incomplete_depth is None:
            incomplete_depth = completed_depth + 1

    elapsed = time.monotonic() - started
    lower_bound = worker.row_count
    record: dict[str, object] = {
        "case": _case_name(rows, cols, periodic_cols),
        "kind": "ring" if rows == 1 and periodic_cols else ("chain_open" if rows == 1 else "grid_open"),
        "rows": rows,
        "cols": cols,
        "n_sites": model.n,
        "n_bonds": len(model.bonds),
        "periodic_cols": periodic_cols,
        "prime": prime,
        "arithmetic": "exact rank over F_p; commutators divided by invertible scalar 2",
        "engine": chosen_engine,
        "point_group_order": len(model.group),
        "global_spin_flip_restriction": "all Pauli support has even Z parity",
        "grading": "(parity of X mask, transpose parity a.b); ad_A and ad_B map fixed blocks",
        "support_orbits": support_orbits,
        "support_complete": support_complete,
        "support_wall": support_wall,
        "dimensions": dimensions,
        "completed_depth": completed_depth,
        "saturated": saturated,
        "exact_dimension": lower_bound if saturated else None,
        "lower_bound": lower_bound,
        "incomplete_depth": incomplete_depth,
        "stop_reason": stop_reason,
        "basis_rows": worker.row_count,
        "basis_columns": (
            support_orbits
            if support_orbits is not None
            else len(worker.discovered_orbits)
        ),
        "grading_block_columns": (
            list(worker.block_sizes) if isinstance(worker, DenseOrbitEngine) else None
        ),
        "grading_block_rows": (
            list(worker.row_counts)
            if isinstance(worker, DenseOrbitEngine)
            else [
                sum(
                    1
                    for row in worker.rows
                    if row and _pauli_grade(next(iter(row)), model.n) == block
                )
                for block in range(4)
            ]
        ),
        "basis_used_bytes": worker.used_basis_bytes,
        "basis_allocated_bytes": worker.allocated_basis_bytes,
        "rss_start_bytes": rss_started,
        "rss_peak_bytes": _rss_bytes(),
        "elapsed_seconds": f"{elapsed:.6f}",
        "growth_estimate": _growth_estimate(dimensions),
    }
    return record


def _artifact(data: dict[str, object], checks: list[dict[str, object]]) -> dict[str, object]:
    return {
        "provenance": {
            "script": SCRIPT,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
            "prime": PRIME_31,
            "second_prime": SECOND_PRIME_31,
            "physical_memory_bytes": _physical_memory(),
        },
        "data": data,
        "checks": checks,
    }


def _write_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run_experiment(quick: bool = False) -> tuple[dict[str, object], dict[str, object]]:
    """Run the requested family and write the two result artifacts."""

    profiles: list[dict[str, object]] = []
    for n in range(2, 33):
        profiles.append(compute_profile(1, n))
    for n in range(3, 33):
        profiles.append(compute_profile(1, n, periodic_cols=True))

    schedules = [
        (2, 2, None, "dense", None),
        (2, 3, None, "dense", None),
        (2, 4, None, "dense", None),
        (3, 3, None, "dense", None),
        (2, 5, None, "dense", 300),
        (2, 6, 14 if quick else 18, "sparse", None if quick else 900),
        (3, 4, 14 if quick else 16, "sparse", None if quick else 1_200),
    ]
    grid_records: list[dict[str, object]] = []
    for rows, cols, depth, selected_engine, time_limit in schedules:
        record = compute_profile(
            rows,
            cols,
            max_depth=depth,
            engine=selected_engine,
            max_seconds=time_limit,
            dense_support_cap=100_000,
            memory_limit_bytes=8 * 1024**3 if selected_engine == "dense" else 4 * 1024**3,
        )
        profiles.append(record)
        grid_records.append(record)

    support_statistics = [
        compute_support_statistics(
            rows,
            cols,
            cap=2_000_000,
        )
        for rows, cols in ((2, 2), (2, 3), (2, 4), (2, 5), (3, 3), (2, 6), (3, 4))
    ]

    known = {"grid_2x2": 11, "grid_2x3": 263, "grid_2x4": 2952}
    second_prime_expected = {**known, "grid_3x3": 8034}
    second_prime_records = [
        compute_profile(rows, cols, prime=SECOND_PRIME_31)
        for rows, cols in ((2, 2), (2, 3), (2, 4), (3, 3))
    ]

    two_by_five_attempt = dict(
        next(row for row in grid_records if row["case"] == "grid_2x5")
    )
    two_by_five_attempt["attempted"] = True

    ladder_certificate = []
    for length in range(2, 7):
        record = next(row for row in profiles if row["case"] == f"grid_2x{length}")
        certificate_depth = 2 * length
        bound = (
            record["dimensions"][certificate_depth - 1]
            if record["completed_depth"] >= certificate_depth
            else None
        )
        ladder_certificate.append(
            {
                "length": length,
                "depth": certificate_depth,
                "D_depth": bound,
                "claimed_finite_bound": 2**length,
                "passed": bound is not None and bound >= 2**length,
            }
        )
    mp.mp.dps = 50
    ladder_size_estimate = {
        "exact_lengths": [2, 3, 4],
        "exact_dimensions": [11, 263, 2952],
        "successive_ratios": ["263/11", "2952/263"],
        "endpoint_base_per_added_rung": mp.nstr(mp.sqrt(mp.mpf(2952) / 11), 14),
        "interpretation": "three-size descriptive estimate only; not an asymptotic theorem",
        "mpmath_dps": 50,
    }

    checks = [
        {
            "name": f"known_dimension_{case}",
            "passed": next(row for row in profiles if row["case"] == case)["exact_dimension"] == expected,
            "detail": f"p={PRIME_31}; expected {expected}",
        }
        for case, expected in known.items()
    ]
    checks.extend(
        {
            "name": f"second_prime_{case}",
            "passed": next(row for row in second_prime_records if row["case"] == case)[
                "exact_dimension"
            ]
            == expected,
            "detail": f"p={SECOND_PRIME_31}; expected {expected}",
        }
        for case, expected in second_prime_expected.items()
    )
    checks.extend(
        {
            "name": f"one_dimensional_formula_{row['case']}",
            "passed": row["saturated"]
            and row["exact_dimension"]
            == (row["n_sites"] ** 2 if row["kind"] == "chain_open" else 3 * row["n_sites"] - 1),
            "detail": "open n^2; ring 3n-1",
        }
        for row in profiles
        if row["kind"] in {"chain_open", "ring"}
    )
    checks.extend(
        [
            {
                "name": "all_requested_profiles_present",
                "passed": len(profiles) == 68,
                "detail": "31 open chains + 30 rings + 7 grids",
            },
            {
                "name": "grid_3x3_completed",
                "passed": next(row for row in grid_records if row["case"] == "grid_3x3")[
                    "saturated"
                ],
                "detail": "dense four-block elimination reached a fixed point",
            },
            {
                "name": "complete_support_counts",
                "passed": quick or all(row["complete"] for row in support_statistics),
                "detail": "non-quick run enumerates every requested point-group orbit support",
            },
            {
                "name": "finite_ladder_exponential_certificate",
                "passed": quick or all(row["passed"] for row in ladder_certificate),
                "detail": "exact modular ranks certify D_(2L) >= 2^L for L=2,...,6",
            },
        ]
    )

    profiles_artifact = _artifact(
        {
            "definition": "D_k = dim span of all left-normed brackets in A,B of length <= k",
            "profiles": profiles,
            "finite_ladder_certificate": ladder_certificate,
            "ladder_size_growth_estimate": ladder_size_estimate,
            "rank_scope": (
                "Every D_k is an exact rank over F_p and therefore a rigorous lower bound over Q. "
                "Agreement at two primes is a cross-check, not a deterministic proof of rational rank."
            ),
        },
        checks,
    )

    exact_keys = (
        "case",
        "n_sites",
        "support_orbits",
        "point_group_order",
        "grading_block_columns",
        "grading_block_rows",
        "engine",
        "exact_dimension",
        "lower_bound",
        "completed_depth",
        "incomplete_depth",
        "stop_reason",
        "basis_used_bytes",
        "basis_allocated_bytes",
        "rss_peak_bytes",
        "elapsed_seconds",
    )
    exact_artifact = _artifact(
        {
            "dimensions": [{key: row[key] for key in exact_keys} for row in grid_records],
            "second_prime_dimensions": [
                {key: row[key] for key in exact_keys} for row in second_prime_records
            ],
            "support_statistics": support_statistics,
            "two_by_five_completion_attempt": two_by_five_attempt,
            "memory_method": (
                "global-spin-flip even-Z support; point-group orbit columns; four exact grading "
                "blocks; int32 rows; int64 products with p^2 < 2^63"
            ),
        },
        [
            check
            for check in checks
            if check["name"].startswith(("known_dimension_", "second_prime_"))
            or check["name"] in {"grid_3x3_completed", "complete_support_counts"}
        ],
    )
    _write_json(RESULT_DIR / "profiles.json", profiles_artifact)
    _write_json(RESULT_DIR / "exact_dimensions.json", exact_artifact)
    return profiles_artifact, exact_artifact


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="use shorter explicit depths for the largest layers")
    arguments = parser.parse_args(argv)
    try:
        profiles, exact = run_experiment(quick=arguments.quick)
        passed = all(check["passed"] for check in profiles["checks"] + exact["checks"])
        for row in exact["data"]["dimensions"]:
            status = row["exact_dimension"] if row["exact_dimension"] is not None else f">= {row['lower_bound']}"
            print(
                f"{row['case']:12s} dim {status!s:>8s} depth {row['completed_depth']:>3d} "
                f"support {str(row['support_orbits']):>8s} {row['elapsed_seconds']}s"
            )
        print("PASS" if passed else "FAIL")
        return 0 if passed else 1
    except Exception as error:  # standalone experiment must report an unambiguous status
        print(f"FAIL: {type(error).__name__}: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
