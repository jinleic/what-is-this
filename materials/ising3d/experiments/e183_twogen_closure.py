#!/usr/bin/env python3
"""Exact two-sum Pauli closures, symmetry containers, and small-size census.

The large anchor closures use exact arithmetic modulo the good Mersenne prime;
each accepted row is an actual nested integer commutator word, hence certifies a
rational lower bound.  The small complete-bipartite cases are additionally
closed over ``fractions.Fraction``, giving exact characteristic-zero dimensions.
The integrated artifact is written by e184.
"""

from __future__ import annotations

import hashlib
import itertools
import platform
import resource
import time
from collections import deque
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

GOOD_PRIME = 2_147_483_647
CPU_BUDGET_SECONDS = 1_800.0
RSS_CAP_BYTES = 1_900_000_000
ROOT = Path(__file__).resolve().parents[1]


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str, budget: float = CPU_BUDGET_SECONDS) -> None:
    used = time.process_time() - started
    if used > budget:
        raise RuntimeError(
            f"NON-DECISIVE process-time expiry at {stage}: {used:.3f}s > {budget:.3f}s"
        )
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise RuntimeError(
            f"NON-DECISIVE RSS expiry at {stage}: {rss} >= {RSS_CAP_BYTES} bytes"
        )


def normalize_edges(n: int, edges: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    answer: set[tuple[int, int]] = set()
    for left, right in edges:
        if not 0 <= left < n or not 0 <= right < n or left == right:
            raise ValueError((n, left, right))
        answer.add(tuple(sorted((left, right))))
    return tuple(sorted(answer))


def x_terms(n: int) -> tuple[int, ...]:
    return tuple(1 << site for site in range(n))


def zz_terms(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    return tuple(((1 << left) | (1 << right)) << n for left, right in edges)


def commutator_sign_half(generator: int, value: int, n: int) -> int:
    """Coefficient of Q_(g xor v) in [Q_g,Q_v]/2 for Q_v=X^a Z^b."""
    mask = (1 << n) - 1
    first = ((((generator >> n) & mask) & (value & mask)).bit_count()) & 1
    second = ((((value >> n) & mask) & (generator & mask)).bit_count()) & 1
    if first == second:
        return 0
    return 1 if first == 0 else -1


def pauli_grade_block(value: int, n: int) -> int:
    """Four-block commutator grading used only to reduce exact row storage."""
    mask = (1 << n) - 1
    x_mask, z_mask = value & mask, value >> n
    return ((x_mask.bit_count() & 1) << 1) | ((x_mask & z_mask).bit_count() & 1)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Apply ``right`` and then ``left``."""
    return tuple(left[right[index]] for index in range(len(left)))


def permute_pauli(value: int, permutation: tuple[int, ...], n: int) -> int:
    mask = (1 << n) - 1
    x_mask, z_mask = value & mask, value >> n

    def move(bits: int) -> int:
        output = 0
        for source, target in enumerate(permutation):
            if (bits >> source) & 1:
                output |= 1 << target
        return output

    return move(x_mask) | (move(z_mask) << n)


def grid_edges(rows: int, columns: int) -> tuple[tuple[int, int], ...]:
    if rows < 1 or columns < 1:
        raise ValueError((rows, columns))
    edges: list[tuple[int, int]] = []
    index = lambda row, column: row * columns + column
    for row in range(rows):
        for column in range(columns):
            if row + 1 < rows:
                edges.append((index(row, column), index(row + 1, column)))
            if column + 1 < columns:
                edges.append((index(row, column), index(row, column + 1)))
    return normalize_edges(rows * columns, edges)


def grid_automorphisms(rows: int, columns: int) -> tuple[tuple[int, ...], ...]:
    n = rows * columns
    transforms = [
        lambda r, c: (r, c),
        lambda r, c: (rows - 1 - r, c),
        lambda r, c: (r, columns - 1 - c),
        lambda r, c: (rows - 1 - r, columns - 1 - c),
    ]
    if rows == columns:
        size = rows
        transforms.extend(
            [
                lambda r, c: (c, r),
                lambda r, c: (size - 1 - c, r),
                lambda r, c: (c, size - 1 - r),
                lambda r, c: (size - 1 - c, size - 1 - r),
            ]
        )
    permutations = set()
    for transform in transforms:
        mapping = []
        for row in range(rows):
            for column in range(columns):
                new_row, new_column = transform(row, column)
                mapping.append(new_row * columns + new_column)
        permutations.add(tuple(mapping))
    identity = tuple(range(n))
    if identity not in permutations:
        raise AssertionError("grid identity missing")
    return tuple(sorted(permutations))


def complete_bipartite_edges(
    n: int, left: Sequence[int], right: Sequence[int]
) -> tuple[tuple[int, int], ...]:
    if set(left) & set(right) or set(left) | set(right) != set(range(n)):
        raise ValueError((n, left, right))
    return normalize_edges(n, ((u, v) for u in left for v in right))


def complete_bipartite_automorphisms(
    n: int, left: Sequence[int], right: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    left = tuple(left)
    right = tuple(right)
    group: set[tuple[int, ...]] = set()
    for left_image in itertools.permutations(left):
        for right_image in itertools.permutations(right):
            mapping = list(range(n))
            for source, target in zip(left, left_image, strict=True):
                mapping[source] = target
            for source, target in zip(right, right_image, strict=True):
                mapping[source] = target
            group.add(tuple(mapping))
    if len(left) == len(right):
        swap = list(range(n))
        for u, v in zip(left, right, strict=True):
            swap[u] = v
            swap[v] = u
        swap_tuple = tuple(swap)
        group |= {compose(swap_tuple, permutation) for permutation in tuple(group)}
    return tuple(sorted(group))


def raw_support(
    n: int, edges: Sequence[tuple[int, int]], *, cap: int = 1 << 22
) -> tuple[int, ...]:
    terms = x_terms(n) + zz_terms(n, edges)
    seen = set(terms)
    queue = deque(terms)
    while queue:
        value = queue.popleft()
        for generator in terms:
            if commutator_sign_half(generator, value, n):
                target = value ^ generator
                if target not in seen:
                    seen.add(target)
                    queue.append(target)
                    if len(seen) > cap:
                        raise MemoryError(f"raw support exceeded {cap}")
    return tuple(sorted(seen, reverse=True))


class OrbitClosure:
    """Closure of invariant Pauli orbit sums under the two generator adjoints."""

    def __init__(
        self,
        n: int,
        edges: Sequence[tuple[int, int]],
        automorphisms: Sequence[tuple[int, ...]],
        *,
        prime: int = GOOD_PRIME,
        started: float | None = None,
    ) -> None:
        self.n = n
        self.edges = normalize_edges(n, edges)
        self.group = tuple(automorphisms)
        if not self.group or any(sorted(permutation) != list(range(n)) for permutation in self.group):
            raise ValueError("invalid automorphism group")
        edge_set = set(self.edges)
        for permutation in self.group:
            moved = {
                tuple(sorted((permutation[left], permutation[right])))
                for left, right in self.edges
            }
            if moved != edge_set:
                raise ValueError("a supplied permutation is not a graph automorphism")
        self.prime = prime
        self.started = time.process_time() if started is None else started
        self.field_terms = x_terms(n)
        self.bond_terms = zz_terms(n, self.edges)
        self.single_terms = self.field_terms + self.bond_terms
        self._canonical_cache: dict[int, int] = {}
        self.support = self._enumerate_orbit_support()
        self.blocks = tuple(
            tuple(value for value in self.support if pauli_grade_block(value, n) == block)
            for block in range(4)
        )
        self.block_sizes = tuple(len(block) for block in self.blocks)
        self.index = tuple(
            {value: position for position, value in enumerate(block)} for block in self.blocks
        )
        self.actions = tuple(
            tuple(
                self._action_table(action_index, terms, source_block)
                for source_block in range(4)
            )
            for action_index, terms in enumerate((self.field_terms, self.bond_terms))
        )
        budget_tick(self.started, f"orbit closure initialization n={n}")

    def canonical(self, value: int) -> int:
        cached = self._canonical_cache.get(value)
        if cached is None:
            cached = min(
                permute_pauli(value, permutation, self.n) for permutation in self.group
            )
            self._canonical_cache[value] = cached
        return cached

    def _enumerate_orbit_support(self) -> tuple[int, ...]:
        seen = {self.canonical(term) for term in self.single_terms}
        queue = deque(sorted(seen))
        while queue:
            value = queue.popleft()
            for generator in self.single_terms:
                if commutator_sign_half(generator, value, self.n):
                    target = self.canonical(value ^ generator)
                    if target not in seen:
                        seen.add(target)
                        queue.append(target)
            if len(seen) % 4096 == 0:
                budget_tick(self.started, f"orbit support n={self.n}, size={len(seen)}")
        return tuple(sorted(seen, reverse=True))

    def _action_table(
        self, action_index: int, terms: tuple[int, ...], source_block: int
    ) -> tuple[np.ndarray, np.ndarray]:
        target_block = source_block ^ (3 if action_index == 0 else 1)
        sources = np.full((len(terms), self.block_sizes[target_block]), -1, dtype=np.int32)
        signs = np.zeros((len(terms), self.block_sizes[target_block]), dtype=np.int8)
        for term_index, generator in enumerate(terms):
            for target_index, target in enumerate(self.blocks[target_block]):
                source = target ^ generator
                sign = commutator_sign_half(generator, source, self.n)
                if sign:
                    source_index = self.index[source_block].get(self.canonical(source))
                    if source_index is None:
                        raise AssertionError("orbit support is not closed")
                    sources[term_index, target_index] = source_index
                    signs[term_index, target_index] = sign
        return sources, signs

    def seed(self, terms: tuple[int, ...]) -> tuple[int, np.ndarray]:
        coefficients: dict[int, int] = {}
        for term in terms:
            orbit = self.canonical(term)
            previous = coefficients.setdefault(orbit, 1)
            if previous != 1:
                raise AssertionError("nonuniform invariant seed")
        blocks = {pauli_grade_block(orbit, self.n) for orbit in coefficients}
        if len(blocks) != 1:
            raise AssertionError("seed crosses grading blocks")
        block = blocks.pop()
        row = np.zeros(self.block_sizes[block], dtype=np.int64)
        for orbit in coefficients:
            row[self.index[block][orbit]] = 1
        return block, row

    def apply_modular(self, action_index: int, word: tuple[int, np.ndarray]) -> tuple[int, np.ndarray]:
        source_block, row = word
        target_block = source_block ^ (3 if action_index == 0 else 1)
        sources, signs = self.actions[action_index][source_block]
        output = np.zeros(self.block_sizes[target_block], dtype=np.int64)
        for source, sign in zip(sources, signs):
            selected = source >= 0
            output[selected] += (
                sign[selected].astype(np.int64) * row[source[selected]].astype(np.int64)
            )
        return target_block, output % self.prime

    def modular_closure(self, *, rank_stop: int | None = None) -> dict[str, object]:
        rows = tuple(np.empty((size, size), dtype=np.int32) for size in self.block_sizes)
        pivot_row = tuple(np.full(size, -1, dtype=np.int32) for size in self.block_sizes)
        pivot_inverse = tuple(np.zeros(size, dtype=np.int32) for size in self.block_sizes)
        row_counts = [0, 0, 0, 0]
        words: list[tuple[int, np.ndarray]] = []
        recipes: list[tuple[int, int, int]] = []

        def reduce_add(vector: tuple[int, np.ndarray]) -> bool:
            block, raw = vector
            value = np.asarray(raw, dtype=np.int64) % self.prime
            while True:
                nonzero = np.flatnonzero(value)
                if not nonzero.size:
                    return False
                lead = int(nonzero[0])
                existing = int(pivot_row[block][lead])
                if existing < 0:
                    row_index = row_counts[block]
                    insertion = int(value[lead])
                    rows[block][row_index] = value.astype(np.int32)
                    pivot_row[block][lead] = row_index
                    pivot_inverse[block][lead] = pow(insertion, self.prime - 2, self.prime)
                    row_counts[block] += 1
                    return True
                factor = int(value[lead]) * int(pivot_inverse[block][lead]) % self.prime
                basis = rows[block][existing, lead:].astype(np.int64)
                product = basis * factor
                product = (product & self.prime) + (product >> 31)
                product = (product & self.prime) + (product >> 31)
                product[product >= self.prime] -= self.prime
                tail = value[lead:]
                tail -= product
                tail[tail < 0] += self.prime

        seeds = (
            (self.seed(self.field_terms), (-1, 0, 1)),
            (self.seed(self.bond_terms), (-1, 1, 1)),
        )
        frontier: list[int] = []
        for word, recipe in seeds:
            if not reduce_add(word):
                raise AssertionError("dependent raw generators")
            words.append((word[0], np.asarray(word[1], dtype=np.int64) % self.prime))
            recipes.append(recipe)
            frontier.append(len(words) - 1)
        stopped = False
        while frontier and not stopped:
            next_frontier: list[int] = []
            for parent in frontier:
                for action_index in (0, 1):
                    candidate = self.apply_modular(action_index, words[parent])
                    if reduce_add(candidate):
                        words.append(candidate)
                        recipes.append((parent, action_index, recipes[parent][2] + 1))
                        next_frontier.append(len(words) - 1)
                        if rank_stop is not None and len(words) >= rank_stop:
                            stopped = True
                            break
                        if len(words) % 256 == 0:
                            budget_tick(
                                self.started,
                                f"modular orbit closure n={self.n}, rank={len(words)}",
                            )
                if stopped:
                    break
            frontier = next_frontier
        recipe_bytes = ";".join(
            f"{parent},{action},{depth}" for parent, action, depth in recipes
        ).encode("ascii")
        return {
            "rank_Fp": len(words),
            "prime": self.prime,
            "saturated": not stopped and not frontier,
            "maximum_depth": max(depth for _, _, depth in recipes),
            "grading_block_ranks": row_counts,
            "orbit_support_size": len(self.support),
            "orbit_block_sizes": list(self.block_sizes),
            "recipe_sha256": hashlib.sha256(recipe_bytes).hexdigest(),
        }

    def _apply_fraction(
        self, action_index: int, word: tuple[int, tuple[Fraction, ...]]
    ) -> tuple[int, tuple[Fraction, ...]]:
        source_block, row = word
        target_block = source_block ^ (3 if action_index == 0 else 1)
        sources, signs = self.actions[action_index][source_block]
        output = [Fraction(0) for _ in range(self.block_sizes[target_block])]
        for source, sign in zip(sources, signs):
            for target_index in np.flatnonzero(source >= 0):
                output[int(target_index)] += (
                    int(sign[target_index]) * row[int(source[target_index])]
                )
        return target_block, tuple(output)

    def rational_closure(self, *, maximum_support: int = 256) -> dict[str, object]:
        if len(self.support) > maximum_support:
            raise ValueError(
                f"rational closure restricted to <= {maximum_support} orbit coordinates"
            )
        pivots: tuple[dict[int, tuple[Fraction, ...]], ...] = tuple({} for _ in range(4))
        row_counts = [0, 0, 0, 0]

        def reduce_add(
            vector: tuple[int, tuple[Fraction, ...]]
        ) -> tuple[int, tuple[Fraction, ...]] | None:
            block, raw = vector
            value = list(raw)
            while True:
                lead = next((index for index, entry in enumerate(value) if entry), None)
                if lead is None:
                    return None
                existing = pivots[block].get(lead)
                if existing is None:
                    scale = value[lead]
                    normalized = tuple(entry / scale for entry in value)
                    pivots[block][lead] = normalized
                    row_counts[block] += 1
                    return block, normalized
                factor = value[lead]
                value = [entry - factor * basis for entry, basis in zip(value, existing, strict=True)]

        seeds: list[tuple[int, tuple[Fraction, ...]]] = []
        for terms in (self.field_terms, self.bond_terms):
            block, row = self.seed(terms)
            seeds.append((block, tuple(Fraction(int(value)) for value in row)))
        frontier: list[tuple[int, tuple[Fraction, ...]]] = []
        for seed in seeds:
            inserted = reduce_add(seed)
            if inserted is None:
                raise AssertionError("dependent rational seeds")
            frontier.append(inserted)
        rank = len(frontier)
        depth = 1
        while frontier:
            next_frontier: list[tuple[int, tuple[Fraction, ...]]] = []
            for word in frontier:
                for action_index in (0, 1):
                    inserted = reduce_add(self._apply_fraction(action_index, word))
                    if inserted is not None:
                        next_frontier.append(inserted)
                        rank += 1
            frontier = next_frontier
            depth += 1
            budget_tick(self.started, f"rational orbit closure n={self.n}, rank={rank}")
        return {
            "dimension_Q": rank,
            "saturated": True,
            "grading_block_ranks": row_counts,
            "orbit_support_size": len(self.support),
            "maximum_frontier_depth": depth - 1,
        }


class CappedRawClosure:
    """Identity-orbit closure stopped as soon as a requested lower bound is met."""

    def __init__(
        self,
        n: int,
        edges: Sequence[tuple[int, int]],
        columns: Sequence[int],
        *,
        prime: int = GOOD_PRIME,
    ) -> None:
        self.n = n
        self.edges = normalize_edges(n, edges)
        self.columns = tuple(columns)
        self.index = {value: position for position, value in enumerate(self.columns)}
        self.prime = prime
        self.groups = (x_terms(n), zz_terms(n, self.edges))
        colarr = np.asarray(self.columns, dtype=np.int64)
        mask = (1 << n) - 1
        x_col, z_col = colarr & mask, colarr >> n
        self.actions: list[tuple[tuple[np.ndarray, np.ndarray], ...]] = []
        for terms in self.groups:
            tables = []
            for generator in terms:
                x_gen, z_gen = generator & mask, generator >> n
                first = np.bitwise_count(np.int64(z_gen) & x_col) & 1
                second = np.bitwise_count(z_col & np.int64(x_gen)) & 1
                anti = (first ^ second).astype(bool)
                targets = np.full(len(self.columns), -1, dtype=np.int32)
                moved = np.bitwise_xor(colarr, np.int64(generator))
                lookup = np.array(
                    [self.index.get(int(value), -1) for value in moved], dtype=np.int32
                )
                targets[anti] = lookup[anti]
                if np.any(targets[anti] < 0):
                    raise AssertionError("universal support is not closed")
                signs = np.zeros(len(self.columns), dtype=np.int8)
                signs[anti] = np.where(first[anti] == 0, 1, -1)
                tables.append((targets, signs))
            self.actions.append(tuple(tables))

    def seed(self, terms: tuple[int, ...]) -> np.ndarray:
        row = np.zeros(len(self.columns), dtype=np.int64)
        for term in terms:
            row[self.index[term]] = 1
        return row

    def apply(self, action_index: int, row: np.ndarray) -> np.ndarray:
        output = np.zeros(len(self.columns), dtype=np.int64)
        for targets, signs in self.actions[action_index]:
            selected = targets >= 0
            output[targets[selected]] += signs[selected].astype(np.int64) * row[selected]
        return output % self.prime

    def close_until(self, rank_stop: int) -> dict[str, object]:
        width = len(self.columns)
        rows = np.empty((rank_stop, width), dtype=np.int32)
        pivot_row = np.full(width, -1, dtype=np.int32)
        pivot_inverse = np.zeros(width, dtype=np.int32)
        rank = 0

        def reduce_add(raw: np.ndarray) -> np.ndarray | None:
            nonlocal rank
            value = np.asarray(raw, dtype=np.int64) % self.prime
            while True:
                nonzero = np.flatnonzero(value)
                if not nonzero.size:
                    return None
                lead = int(nonzero[0])
                existing = int(pivot_row[lead])
                if existing < 0:
                    insertion = int(value[lead])
                    normalized = value * pow(insertion, self.prime - 2, self.prime)
                    normalized = (normalized & self.prime) + (normalized >> 31)
                    normalized = (normalized & self.prime) + (normalized >> 31)
                    normalized[normalized >= self.prime] -= self.prime
                    rows[rank] = normalized.astype(np.int32)
                    pivot_row[lead] = rank
                    pivot_inverse[lead] = 1
                    rank += 1
                    return normalized
                factor = int(value[lead]) * int(pivot_inverse[lead]) % self.prime
                product = rows[existing, lead:].astype(np.int64) * factor
                product = (product & self.prime) + (product >> 31)
                product = (product & self.prime) + (product >> 31)
                product[product >= self.prime] -= self.prime
                tail = value[lead:]
                tail -= product
                tail[tail < 0] += self.prime

        frontier: list[np.ndarray] = []
        for terms in self.groups:
            inserted = reduce_add(self.seed(terms))
            if inserted is None:
                raise AssertionError("dependent raw seeds")
            frontier.append(inserted)
        stopped = rank >= rank_stop
        while frontier and not stopped:
            next_frontier: list[np.ndarray] = []
            for word in frontier:
                for action_index in (0, 1):
                    inserted = reduce_add(self.apply(action_index, word))
                    if inserted is not None:
                        next_frontier.append(inserted)
                        if rank >= rank_stop:
                            stopped = True
                            break
                if stopped:
                    break
            frontier = next_frontier
        return {
            "rank_Fp": rank,
            "saturated": not stopped and not frontier,
            "stopped_at_requested_rank": stopped,
            "support_size": width,
        }


def complete_bipartite_case(left_size: int, right_size: int) -> dict[str, object]:
    n = left_size + right_size
    even_positions = tuple(range(0, n, 2))
    odd_positions = tuple(range(1, n, 2))
    if (len(even_positions), len(odd_positions)) == (left_size, right_size):
        left, right = even_positions, odd_positions
    elif (len(odd_positions), len(even_positions)) == (left_size, right_size):
        left, right = odd_positions, even_positions
    else:
        raise ValueError(
            "K_{a,b} carries an alternating Hamiltonian path only when |a-b| <= 1"
        )
    edges = complete_bipartite_edges(n, left, right)
    group = complete_bipartite_automorphisms(n, left, right)
    closure = OrbitClosure(n, edges, group)
    rational = closure.rational_closure()
    return {
        "tag": "[COMPUTATION]",
        "graph": f"K_{{{left_size},{right_size}}}",
        "n_vertices": n,
        "parts": [list(left), list(right)],
        "edges": [list(edge) for edge in edges],
        "hamiltonian_path": list(range(n)),
        "maximum_degree": max(left_size, right_size),
        "automorphism_group_order": len(group),
        "symmetry_fixed_local_term_upper_bound": len(closure.support),
        "rational_closure": rational,
        "quadratic_ceiling": n * (2 * n - 1),
        "clears_quadratic_ceiling": rational["dimension_Q"] > n * (2 * n - 1),
    }


def small_hamiltonian_census() -> dict[str, object]:
    started = time.process_time()
    rows: list[dict[str, object]] = []
    for n in (4, 5, 6, 7):
        path = {(site, site + 1) for site in range(n - 1)}
        chords = tuple(
            (left, right)
            for left in range(n)
            for right in range(left + 1, n)
            if (right - left) % 2 == 1 and (left, right) not in path
        )
        branching_graphs: list[tuple[tuple[int, int], ...]] = []
        for mask in range(1 << len(chords)):
            edges = normalize_edges(
                n,
                path
                | {
                    chords[index]
                    for index in range(len(chords))
                    if (mask >> index) & 1
                },
            )
            degrees = [0] * n
            for left, right in edges:
                degrees[left] += 1
                degrees[right] += 1
            if max(degrees) >= 3:
                branching_graphs.append(edges)
        if not branching_graphs:
            rows.append(
                {
                    "n_vertices": n,
                    "allowed_chords": [list(edge) for edge in chords],
                    "branching_graph_count": 0,
                    "quadratic_ceiling": n * (2 * n - 1),
                    "cases_clearing_ceiling": 0,
                    "exceptions": [],
                }
            )
            continue
        largest_left = tuple(range(0, n, 2))
        largest_right = tuple(range(1, n, 2))
        universal_edges = complete_bipartite_edges(n, largest_left, largest_right)
        columns = raw_support(n, universal_edges)
        ceiling = n * (2 * n - 1)
        clearing = 0
        exceptions: list[dict[str, object]] = []
        for edges in branching_graphs:
            result = CappedRawClosure(n, edges, columns).close_until(ceiling + 1)
            if result["rank_Fp"] > ceiling:
                clearing += 1
            else:
                exceptions.append(
                    {
                        "edges": [list(edge) for edge in edges],
                        "rank_Fp": result["rank_Fp"],
                        "saturated_mod_prime": result["saturated"],
                    }
                )
        rows.append(
            {
                "n_vertices": n,
                "allowed_chords": [list(edge) for edge in chords],
                "branching_graph_count": len(branching_graphs),
                "quadratic_ceiling": ceiling,
                "cases_clearing_ceiling": clearing,
                "exceptions": exceptions,
                "universal_local_term_support_size": len(columns),
            }
        )
        budget_tick(started, f"small Hamiltonian census n={n}", budget=600.0)
    return {
        "tag": "[COMPUTATION]",
        "scope": (
            "Every bipartite graph carrying the fixed labeled Hamiltonian path at n=4..7; "
            "odd-distance chord subsets exhaust all such relabelled graphs."
        ),
        "rows": rows,
        "process_time_seconds": round(time.process_time() - started, 6),
    }


ANCHOR_SPECS = (
    ("open_2x2", 2, 2, 11, 28),
    ("open_2x3", 2, 3, 263, 66),
    ("open_2x4", 2, 4, 2952, 120),
    ("open_3x3", 3, 3, 8034, 153),
)


def run_anchor_closures() -> tuple[dict[str, object], ...]:
    started = time.process_time()
    rows: list[dict[str, object]] = []
    for label, grid_rows, grid_columns, expected, ceiling in ANCHOR_SPECS:
        n = grid_rows * grid_columns
        edges = grid_edges(grid_rows, grid_columns)
        group = grid_automorphisms(grid_rows, grid_columns)
        closure = OrbitClosure(n, edges, group, started=started)
        modular = closure.modular_closure()
        if not modular["saturated"] or modular["rank_Fp"] != expected:
            raise AssertionError((label, modular, expected))
        rows.append(
            {
                "tag": "[COMPUTATION][EXTERNAL exact-Q upper certificate]",
                "label": label,
                "shape": [grid_rows, grid_columns],
                "n_vertices": n,
                "edges": [list(edge) for edge in edges],
                "automorphism_group_order": len(group),
                "raw_generator_modular_closure": modular,
                "certified_dimension_Q": expected,
                "quadratic_ceiling": ceiling,
                "clears_quadratic_ceiling": expected > ceiling,
                "exact_Q_scope": (
                    "This producer re-derives the good-prime raw-word closure; the opposite "
                    "characteristic-zero inequality is supplied by the cited existing exact-Q "
                    "sector certificate."
                ),
            }
        )
        budget_tick(started, f"anchor {label}")
    return tuple(rows)


def run_closure_audit(*, include_anchors: bool = True) -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    counterexamples = (
        complete_bipartite_case(2, 3),
        complete_bipartite_case(3, 3),
        complete_bipartite_case(3, 4),
    )
    expected = {"K_{2,3}": 44, "K_{3,3}": 63, "K_{3,4}": 167}
    checks: list[dict[str, object]] = []
    for row in counterexamples:
        dimension = row["rational_closure"]["dimension_Q"]
        checks.append(
            {
                "name": f"exact rational closure {row['graph']}",
                "passed": dimension == expected[row["graph"]],
                "detail": (
                    f"dimension={dimension}, symmetry upper="
                    f"{row['symmetry_fixed_local_term_upper_bound']}"
                ),
            }
        )
    census = small_hamiltonian_census()
    census_by_n = {row["n_vertices"]: row for row in census["rows"]}
    checks.extend(
        [
            {
                "name": "n=5 exhaustive exception is K_2,3",
                "passed": (
                    census_by_n[5]["branching_graph_count"] == 3
                    and census_by_n[5]["cases_clearing_ceiling"] == 2
                    and len(census_by_n[5]["exceptions"]) == 1
                    and census_by_n[5]["exceptions"][0]["rank_Fp"] == 44
                ),
                "detail": "3 graphs; 2 clear 45; K_2,3 saturates at 44",
            },
            {
                "name": "n=6 exhaustive exception is K_3,3",
                "passed": (
                    census_by_n[6]["branching_graph_count"] == 14
                    and census_by_n[6]["cases_clearing_ceiling"] == 13
                    and len(census_by_n[6]["exceptions"]) == 1
                    and census_by_n[6]["exceptions"][0]["rank_Fp"] == 63
                ),
                "detail": "14 graphs; 13 clear 66; K_3,3 saturates at 63",
            },
            {
                "name": "all n=7 branching Hamiltonian bipartite graphs clear",
                "passed": (
                    census_by_n[7]["branching_graph_count"] == 63
                    and census_by_n[7]["cases_clearing_ceiling"] == 63
                    and not census_by_n[7]["exceptions"]
                ),
                "detail": "finite n=7 census only; not an all-size induction",
            },
        ]
    )
    anchors = run_anchor_closures() if include_anchors else ()
    if anchors:
        checks.append(
            {
                "name": "four raw grid anchor closures",
                "passed": [
                    row["raw_generator_modular_closure"]["rank_Fp"] for row in anchors
                ]
                == [11, 263, 2952, 8034],
                "detail": "2x2, 2x3, 2x4, 3x3",
            }
        )
    checks.append(
        {
            "name": "resource cap",
            "passed": max_rss_bytes() < RSS_CAP_BYTES,
            "detail": f"peak_rss_bytes={max_rss_bytes()}",
        }
    )
    if not all(bool(row["passed"]) for row in checks):
        raise AssertionError("two-generator closure audit failed")
    budget_tick(started, "completed closure audit")
    data = {
        "tag": "[COMPUTATION][THEOREM — finite counterexamples]",
        "complete_bipartite_cases": list(counterexamples),
        "small_hamiltonian_census": census,
        "grid_anchors": list(anchors),
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }
    return data, checks


def main() -> int:
    # The standalone e183 invocation is intentionally the light theory/counterexample pass.
    # e184 runs the four larger raw grid anchors before writing the integrated artifact.
    data, checks = run_closure_audit(include_anchors=False)
    print(
        f"PASS e183: {len(checks)} checks, "
        f"cpu={data['process_time_seconds']}s, rss={data['peak_rss_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
