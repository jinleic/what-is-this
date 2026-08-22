#!/usr/bin/env python3
"""Clean-room verifier for the wave-19 two-generator obstruction.

This file imports no experiment module and no project Clifford helper.  It
rebuilds the raw packed-Pauli generators, the point-group orbit supports, all
four good-prime anchor closures, the exact small complete-bipartite rational
closures, the n<=7 fixed-path census, and the endpoint ad_(iA) brackets.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import time
from collections import deque
from fractions import Fraction
from math import comb
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "twogen_allsize.json"
GOOD_PRIME = 2_147_483_647
CPU_LIMIT = 1_800.0
RSS_LIMIT = 1_900_000_000
STARTED = time.process_time()
FAILURES: list[str] = []


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def guard(stage: str) -> None:
    elapsed = time.process_time() - STARTED
    if elapsed > CPU_LIMIT:
        raise RuntimeError(f"process-time wall at {stage}: {elapsed:.3f}s")
    if max_rss_bytes() >= RSS_LIMIT:
        raise RuntimeError(f"RSS wall at {stage}: {max_rss_bytes()}")


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"{'PASS' if passed else 'FAIL'}: {name}" + (f" ({detail})" if detail else ""))
    if not passed:
        FAILURES.append(name)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def edges_normalized(n: int, edges: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    result: set[tuple[int, int]] = set()
    for u, v in edges:
        if u == v or not 0 <= u < n or not 0 <= v < n:
            raise ValueError((n, u, v))
        result.add((min(u, v), max(u, v)))
    return tuple(sorted(result))


def grid_bonds(height: int, width: int) -> tuple[tuple[int, int], ...]:
    bonds = []
    for r in range(height):
        for c in range(width):
            site = r * width + c
            if r + 1 < height:
                bonds.append((site, site + width))
            if c + 1 < width:
                bonds.append((site, site + 1))
    return edges_normalized(height * width, bonds)


def grid_group(height: int, width: int) -> tuple[tuple[int, ...], ...]:
    functions = [
        lambda r, c: (r, c),
        lambda r, c: (height - 1 - r, c),
        lambda r, c: (r, width - 1 - c),
        lambda r, c: (height - 1 - r, width - 1 - c),
    ]
    if height == width:
        s = height
        functions += [
            lambda r, c: (c, r),
            lambda r, c: (s - 1 - c, r),
            lambda r, c: (c, s - 1 - r),
            lambda r, c: (s - 1 - c, s - 1 - r),
        ]
    group = set()
    for function in functions:
        group.add(
            tuple(
                function(r, c)[0] * width + function(r, c)[1]
                for r in range(height)
                for c in range(width)
            )
        )
    return tuple(sorted(group))


def bipartite_bonds(left: Sequence[int], right: Sequence[int]) -> tuple[tuple[int, int], ...]:
    n = len(left) + len(right)
    return edges_normalized(n, ((u, v) for u in left for v in right))


def permutation_compose(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a[b[index]] for index in range(len(a)))


def bipartite_group(left: Sequence[int], right: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    left, right = tuple(left), tuple(right)
    n = len(left) + len(right)
    group = set()
    for left_order in itertools.permutations(left):
        for right_order in itertools.permutations(right):
            mapping = list(range(n))
            for source, target in zip(left, left_order, strict=True):
                mapping[source] = target
            for source, target in zip(right, right_order, strict=True):
                mapping[source] = target
            group.add(tuple(mapping))
    if len(left) == len(right):
        swap = list(range(n))
        for u, v in zip(left, right, strict=True):
            swap[u], swap[v] = v, u
        swap = tuple(swap)
        group |= {permutation_compose(swap, value) for value in tuple(group)}
    return tuple(sorted(group))


def fields(n: int) -> tuple[int, ...]:
    return tuple(1 << site for site in range(n))


def bonds_packed(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    return tuple(((1 << u) | (1 << v)) << n for u, v in edges)


def half_bracket_sign(g: int, v: int, n: int) -> int:
    mask = (1 << n) - 1
    left = ((((g >> n) & mask) & (v & mask)).bit_count()) & 1
    right = ((((v >> n) & mask) & (g & mask)).bit_count()) & 1
    return 0 if left == right else (1 if left == 0 else -1)


def move_bits(bits: int, permutation: tuple[int, ...]) -> int:
    output = 0
    for source, target in enumerate(permutation):
        output |= ((bits >> source) & 1) << target
    return output


def move_pauli(value: int, permutation: tuple[int, ...], n: int) -> int:
    mask = (1 << n) - 1
    return move_bits(value & mask, permutation) | (move_bits(value >> n, permutation) << n)


def block_label(value: int, n: int) -> int:
    mask = (1 << n) - 1
    x, z = value & mask, value >> n
    return ((x.bit_count() & 1) << 1) | ((x & z).bit_count() & 1)


class CleanOrbitClosure:
    """Independent orbit-coordinate implementation of the two adjoint actions."""

    def __init__(
        self,
        n: int,
        edges: Sequence[tuple[int, int]],
        group: Sequence[tuple[int, ...]],
    ) -> None:
        self.n = n
        self.edges = edges_normalized(n, edges)
        self.group = tuple(group)
        self.xs = fields(n)
        self.zzs = bonds_packed(n, self.edges)
        self.terms = self.xs + self.zzs
        self.cache: dict[int, int] = {}
        self.representatives = self._support_orbits()
        self.parts = tuple(
            tuple(value for value in self.representatives if block_label(value, n) == label)
            for label in range(4)
        )
        self.positions = tuple(
            {value: index for index, value in enumerate(part)} for part in self.parts
        )
        self.maps = tuple(
            tuple(
                self._map(action, terms, source)
                for source in range(4)
            )
            for action, terms in enumerate((self.xs, self.zzs))
        )
        guard(f"orbit initialization n={n}")

    def canonical(self, value: int) -> int:
        answer = self.cache.get(value)
        if answer is None:
            answer = min(move_pauli(value, permutation, self.n) for permutation in self.group)
            self.cache[value] = answer
        return answer

    def _support_orbits(self) -> tuple[int, ...]:
        seen = {self.canonical(term) for term in self.terms}
        queue = deque(seen)
        while queue:
            value = queue.popleft()
            for generator in self.terms:
                if half_bracket_sign(generator, value, self.n):
                    candidate = self.canonical(value ^ generator)
                    if candidate not in seen:
                        seen.add(candidate)
                        queue.append(candidate)
        return tuple(sorted(seen, reverse=True))

    def _map(
        self, action: int, terms: tuple[int, ...], source_block: int
    ) -> tuple[np.ndarray, np.ndarray]:
        target_block = source_block ^ (3 if action == 0 else 1)
        source_indices = np.full((len(terms), len(self.parts[target_block])), -1, dtype=np.int32)
        coefficients = np.zeros(source_indices.shape, dtype=np.int8)
        for term_index, term in enumerate(terms):
            for target_index, target in enumerate(self.parts[target_block]):
                source = target ^ term
                sign = half_bracket_sign(term, source, self.n)
                if sign:
                    source_indices[term_index, target_index] = self.positions[source_block][
                        self.canonical(source)
                    ]
                    coefficients[term_index, target_index] = sign
        return source_indices, coefficients

    def seed(self, terms: tuple[int, ...]) -> tuple[int, np.ndarray]:
        orbits = {self.canonical(term) for term in terms}
        labels = {block_label(value, self.n) for value in orbits}
        if len(labels) != 1:
            raise AssertionError("seed grading")
        label = labels.pop()
        row = np.zeros(len(self.parts[label]), dtype=np.int64)
        for value in orbits:
            row[self.positions[label][value]] = 1
        return label, row

    def act_mod(self, action: int, word: tuple[int, np.ndarray]) -> tuple[int, np.ndarray]:
        source_block, row = word
        target_block = source_block ^ (3 if action == 0 else 1)
        sources, signs = self.maps[action][source_block]
        result = np.zeros(len(self.parts[target_block]), dtype=np.int64)
        for source, sign in zip(sources, signs):
            valid = source >= 0
            result[valid] += sign[valid].astype(np.int64) * row[source[valid]]
        return target_block, result % GOOD_PRIME

    def rank_mod_prime(self) -> dict[str, object]:
        sizes = tuple(len(part) for part in self.parts)
        echelon = tuple(np.empty((size, size), dtype=np.int32) for size in sizes)
        owner = tuple(np.full(size, -1, dtype=np.int32) for size in sizes)
        inverse = tuple(np.zeros(size, dtype=np.int32) for size in sizes)
        counts = [0, 0, 0, 0]

        def insert(word: tuple[int, np.ndarray]) -> bool:
            label, raw = word
            row = np.asarray(raw, dtype=np.int64) % GOOD_PRIME
            while True:
                occupied = np.flatnonzero(row)
                if not occupied.size:
                    return False
                pivot = int(occupied[0])
                prior = int(owner[label][pivot])
                if prior < 0:
                    location = counts[label]
                    echelon[label][location] = row.astype(np.int32)
                    owner[label][pivot] = location
                    inverse[label][pivot] = pow(int(row[pivot]), GOOD_PRIME - 2, GOOD_PRIME)
                    counts[label] += 1
                    return True
                factor = int(row[pivot]) * int(inverse[label][pivot]) % GOOD_PRIME
                product = echelon[label][prior, pivot:].astype(np.int64) * factor
                product = (product & GOOD_PRIME) + (product >> 31)
                product = (product & GOOD_PRIME) + (product >> 31)
                product[product >= GOOD_PRIME] -= GOOD_PRIME
                row[pivot:] -= product
                row[pivot:][row[pivot:] < 0] += GOOD_PRIME

        pending: list[tuple[int, np.ndarray]] = []
        for terms in (self.xs, self.zzs):
            seed = self.seed(terms)
            if not insert(seed):
                raise AssertionError("dependent generators")
            pending.append(seed)
        head = 0
        while head < len(pending):
            word = pending[head]
            head += 1
            for action in (0, 1):
                candidate = self.act_mod(action, word)
                if insert(candidate):
                    pending.append(candidate)
                    if len(pending) % 256 == 0:
                        guard(f"clean modular closure n={self.n}, rank={len(pending)}")
        return {
            "rank": len(pending),
            "block_ranks": counts,
            "orbit_support": len(self.representatives),
            "block_sizes": list(sizes),
        }

    def act_fraction(
        self, action: int, word: tuple[int, tuple[Fraction, ...]]
    ) -> tuple[int, tuple[Fraction, ...]]:
        source_block, row = word
        target_block = source_block ^ (3 if action == 0 else 1)
        sources, signs = self.maps[action][source_block]
        result = [Fraction(0) for _ in self.parts[target_block]]
        for source, sign in zip(sources, signs):
            for target in np.flatnonzero(source >= 0):
                result[int(target)] += int(sign[target]) * row[int(source[target])]
        return target_block, tuple(result)

    def rank_over_q(self) -> int:
        if len(self.representatives) > 256:
            raise ValueError("small exact-Q verifier only")
        bases: tuple[dict[int, tuple[Fraction, ...]], ...] = tuple({} for _ in range(4))

        def insert(
            word: tuple[int, tuple[Fraction, ...]]
        ) -> tuple[int, tuple[Fraction, ...]] | None:
            label, raw = word
            row = list(raw)
            while True:
                pivot = next((index for index, value in enumerate(row) if value), None)
                if pivot is None:
                    return None
                prior = bases[label].get(pivot)
                if prior is None:
                    scale = row[pivot]
                    normalized = tuple(value / scale for value in row)
                    bases[label][pivot] = normalized
                    return label, normalized
                scale = row[pivot]
                row = [value - scale * old for value, old in zip(row, prior, strict=True)]

        pending: list[tuple[int, tuple[Fraction, ...]]] = []
        for terms in (self.xs, self.zzs):
            label, seed = self.seed(terms)
            accepted = insert((label, tuple(Fraction(int(value)) for value in seed)))
            if accepted is None:
                raise AssertionError("rational seed")
            pending.append(accepted)
        head = 0
        while head < len(pending):
            word = pending[head]
            head += 1
            for action in (0, 1):
                accepted = insert(self.act_fraction(action, word))
                if accepted is not None:
                    pending.append(accepted)
            guard(f"clean rational closure n={self.n}, rank={len(pending)}")
        return len(pending)


def raw_support(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    terms = fields(n) + bonds_packed(n, edges)
    seen = set(terms)
    queue = deque(terms)
    while queue:
        value = queue.popleft()
        for term in terms:
            if half_bracket_sign(term, value, n):
                candidate = value ^ term
                if candidate not in seen:
                    seen.add(candidate)
                    queue.append(candidate)
    return tuple(sorted(seen, reverse=True))


class CleanCappedClosure:
    def __init__(self, n: int, edges: Sequence[tuple[int, int]], support: Sequence[int]) -> None:
        self.n = n
        self.edges = tuple(edges)
        self.support = tuple(support)
        self.position = {value: index for index, value in enumerate(self.support)}
        self.groups = (fields(n), bonds_packed(n, self.edges))
        columns = np.asarray(self.support, dtype=np.int64)
        mask = (1 << n) - 1
        x_columns, z_columns = columns & mask, columns >> n
        self.tables = []
        for terms in self.groups:
            action_tables = []
            for term in terms:
                x_term, z_term = term & mask, term >> n
                first = np.bitwise_count(np.int64(z_term) & x_columns) & 1
                second = np.bitwise_count(z_columns & np.int64(x_term)) & 1
                anti = (first ^ second).astype(bool)
                target = np.full(len(columns), -1, dtype=np.int32)
                moved = np.bitwise_xor(columns, np.int64(term))
                lookup = np.array([self.position.get(int(value), -1) for value in moved], dtype=np.int32)
                target[anti] = lookup[anti]
                sign = np.zeros(len(columns), dtype=np.int8)
                sign[anti] = np.where(first[anti] == 0, 1, -1)
                action_tables.append((target, sign))
            self.tables.append(tuple(action_tables))

    def act(self, action: int, row: np.ndarray) -> np.ndarray:
        result = np.zeros(len(self.support), dtype=np.int64)
        for target, sign in self.tables[action]:
            valid = target >= 0
            result[target[valid]] += sign[valid].astype(np.int64) * row[valid]
        return result % GOOD_PRIME

    def until(self, target_rank: int) -> tuple[int, bool]:
        width = len(self.support)
        basis = np.empty((target_rank, width), dtype=np.int32)
        owners = np.full(width, -1, dtype=np.int32)
        rank = 0

        def insert(raw: np.ndarray) -> np.ndarray | None:
            nonlocal rank
            row = np.asarray(raw, dtype=np.int64) % GOOD_PRIME
            while True:
                occupied = np.flatnonzero(row)
                if not occupied.size:
                    return None
                pivot = int(occupied[0])
                prior = int(owners[pivot])
                if prior < 0:
                    inv = pow(int(row[pivot]), GOOD_PRIME - 2, GOOD_PRIME)
                    row *= inv
                    row = (row & GOOD_PRIME) + (row >> 31)
                    row = (row & GOOD_PRIME) + (row >> 31)
                    row[row >= GOOD_PRIME] -= GOOD_PRIME
                    basis[rank] = row.astype(np.int32)
                    owners[pivot] = rank
                    rank += 1
                    return row
                factor = int(row[pivot])
                product = basis[prior, pivot:].astype(np.int64) * factor
                product = (product & GOOD_PRIME) + (product >> 31)
                product = (product & GOOD_PRIME) + (product >> 31)
                product[product >= GOOD_PRIME] -= GOOD_PRIME
                row[pivot:] -= product
                row[pivot:][row[pivot:] < 0] += GOOD_PRIME

        pending = []
        for terms in self.groups:
            row = np.zeros(width, dtype=np.int64)
            for term in terms:
                row[self.position[term]] = 1
            accepted = insert(row)
            if accepted is None:
                raise AssertionError("capped seed")
            pending.append(accepted)
        head = 0
        while head < len(pending) and rank < target_rank:
            row = pending[head]
            head += 1
            for action in (0, 1):
                accepted = insert(self.act(action, row))
                if accepted is not None:
                    pending.append(accepted)
                    if rank >= target_rank:
                        break
        return rank, head == len(pending) and rank < target_rank


def complete_case(a: int, b: int) -> tuple[CleanOrbitClosure, int]:
    n = a + b
    evens, odds = tuple(range(0, n, 2)), tuple(range(1, n, 2))
    if (len(evens), len(odds)) == (a, b):
        left, right = evens, odds
    else:
        left, right = odds, evens
    closure = CleanOrbitClosure(n, bipartite_bonds(left, right), bipartite_group(left, right))
    return closure, closure.rank_over_q()


def spectrum(n: int, grade: int) -> dict[int, int]:
    result = {}
    for minus in range(max(0, grade - n), min(n, grade) + 1):
        plus = grade - minus
        result[2 * (plus - minus)] = comb(n, minus) * comb(n, plus)
    return result


def add_vector(target: dict[tuple[str, ...], Fraction], word: tuple[str, ...], value: Fraction) -> None:
    new_value = target.get(word, Fraction(0)) + value
    if new_value:
        target[word] = new_value
    else:
        target.pop(word, None)


def d_action(vector: dict[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    result: dict[tuple[str, ...], Fraction] = {}
    for word, coefficient in vector.items():
        for site, letter in enumerate(word):
            if letter == "Z":
                changed = word[:site] + ("Y",) + word[site + 1 :]
                add_vector(result, changed, 2 * coefficient)
            elif letter == "Y":
                changed = word[:site] + ("Z",) + word[site + 1 :]
                add_vector(result, changed, -2 * coefficient)
    return result


def vector_sum(*vectors: dict[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    result: dict[tuple[str, ...], Fraction] = {}
    for vector in vectors:
        for word, coefficient in vector.items():
            add_vector(result, word, coefficient)
    return result


def vector_scale(
    vector: dict[tuple[str, ...], Fraction], coefficient: Fraction
) -> dict[tuple[str, ...], Fraction]:
    return {word: value * coefficient for word, value in vector.items() if value * coefficient}


def bond_sum(n: int, edges: Sequence[tuple[int, int]], left: str, right: str) -> dict[tuple[str, ...], Fraction]:
    result = {}
    for u, v in edges:
        word = ["I"] * n
        word[u], word[v] = left, right
        result[tuple(word)] = Fraction(1)
    return result


def verify_endpoint_brackets(label: str, n: int, edges: Sequence[tuple[int, int]]) -> bool:
    zz = bond_sum(n, edges, "Z", "Z")
    yy = bond_sum(n, edges, "Y", "Y")
    yz = vector_sum(bond_sum(n, edges, "Y", "Z"), bond_sum(n, edges, "Z", "Y"))
    d1 = d_action(zz)
    d2 = d_action(d1)
    d3 = d_action(d2)
    zero = vector_sum(zz, vector_scale(d2, Fraction(1, 16)))
    cosine = vector_scale(d2, Fraction(-1, 16))
    sine = vector_scale(d1, Fraction(1, 4))
    expected_zero = vector_scale(vector_sum(zz, yy), Fraction(1, 2))
    expected_cosine = vector_scale(vector_sum(zz, vector_scale(yy, -1)), Fraction(1, 2))
    expected_sine = vector_scale(yz, Fraction(1, 2))
    passed = (
        vector_sum(d3, vector_scale(d1, 16)) == {}
        and zero == expected_zero
        and cosine == expected_cosine
        and sine == expected_sine
        and d_action(zero) == {}
        and d_action(cosine) == vector_scale(sine, 4)
        and d_action(sine) == vector_scale(cosine, -4)
        and vector_sum(zero, cosine) == zz
    )
    check(f"exact ad_(iA) endpoint brackets {label}", passed, f"n={n}, edges={len(edges)}")
    return passed


def verify_census() -> dict[int, tuple[int, int, list[int]]]:
    result = {}
    for n in (4, 5, 6, 7):
        path = {(site, site + 1) for site in range(n - 1)}
        chords = tuple(
            (u, v)
            for u in range(n)
            for v in range(u + 1, n)
            if (v - u) % 2 and (u, v) not in path
        )
        graphs = []
        for mask in range(1 << len(chords)):
            edges = edges_normalized(
                n,
                path | {chords[index] for index in range(len(chords)) if (mask >> index) & 1},
            )
            degree = [0] * n
            for u, v in edges:
                degree[u] += 1
                degree[v] += 1
            if max(degree) >= 3:
                graphs.append(edges)
        if not graphs:
            result[n] = (0, 0, [])
            continue
        evens, odds = tuple(range(0, n, 2)), tuple(range(1, n, 2))
        support = raw_support(n, bipartite_bonds(evens, odds))
        ceiling = n * (2 * n - 1)
        clearing = 0
        exceptions = []
        for edges in graphs:
            rank, saturated = CleanCappedClosure(n, edges, support).until(ceiling + 1)
            if rank > ceiling:
                clearing += 1
            else:
                if not saturated:
                    raise AssertionError("sub-ceiling case did not saturate")
                exceptions.append(rank)
        result[n] = (len(graphs), clearing, exceptions)
        guard(f"clean census n={n}")
    return result


def verify_external_factor_arithmetic() -> None:
    factor_sums = {
        "2x2": 11,
        "2x3": 1 + 105 + 21 + 21 + 80 + 35,
        "2x4": 1 + 861 + 2 * 351 + 575 + 783 + 2 * 15,
        "3x3": 1 + sum(size * size - 1 for size in (33, 48, 59, 3, 16, 30)),
    }
    check(
        "independent exact-Q factor arithmetic",
        list(factor_sums.values()) == [11, 263, 2952, 8034],
        str(factor_sums),
    )


def main() -> int:
    if not ARTIFACT.is_file():
        raise FileNotFoundError(ARTIFACT)

    # 1. The exact spectrum and the load-bearing edge coincidence.
    for n in (6, 8, 9):
        for grade in range(2, 2 * n, 4):
            row = spectrum(n, grade)
            check(
                f"grade spectrum n={n}, k={grade}",
                sum(row.values()) == comb(2 * n, grade) and {-4, 0, 4} <= set(row),
                f"eigenvalue_i_coefficients={sorted(row)}",
            )
    graph_brackets = (
        ("open 2x3", 6, grid_bonds(2, 3)),
        ("open 2x4", 8, grid_bonds(2, 4)),
        ("open 3x3", 9, grid_bonds(3, 3)),
        ("non-grid K_3,3", 6, bipartite_bonds((0, 2, 4), (1, 3, 5))),
        ("non-grid K_3,4", 7, bipartite_bonds((1, 3, 5), (0, 2, 4, 6))),
        ("2x2 method exception", 4, grid_bonds(2, 2)),
    )
    for label, n, edges in graph_brackets:
        verify_endpoint_brackets(label, n, edges)

    # 2. Exact-Q complete-bipartite counterexamples and non-grid positive control.
    small_dimensions = {}
    small_upper_bounds = {}
    for a, b in ((2, 3), (3, 3), (3, 4)):
        closure, dimension = complete_case(a, b)
        label = f"K_{{{a},{b}}}"
        small_dimensions[label] = dimension
        small_upper_bounds[label] = len(closure.representatives)
    check(
        "exact rational complete-bipartite dimensions",
        small_dimensions == {"K_{2,3}": 44, "K_{3,3}": 63, "K_{3,4}": 167},
        f"dimensions={small_dimensions}, symmetry_upper={small_upper_bounds}",
    )
    check(
        "headline quadratic comparisons",
        44 < 5 * 9 and 63 < 6 * 11 and 167 > 7 * 13,
        "44<45, 63<66, 167>91",
    )

    # 3. Exhaust the entire relabelled fixed-path graph class through n=7.
    census = verify_census()
    check(
        "independent fixed-path census n=4..7",
        census == {
            4: (0, 0, []),
            5: (3, 2, [44]),
            6: (14, 13, [63]),
            7: (63, 63, []),
        },
        str(census),
    )

    # 4. Rebuild all four raw generators and saturate their good-prime closures.
    anchor_specs = (
        ("open_2x2", 2, 2, 11),
        ("open_2x3", 2, 3, 263),
        ("open_2x4", 2, 4, 2952),
        ("open_3x3", 3, 3, 8034),
    )
    anchor_results = {}
    for label, height, width, expected in anchor_specs:
        closure = CleanOrbitClosure(
            height * width,
            grid_bonds(height, width),
            grid_group(height, width),
        )
        result = closure.rank_mod_prime()
        anchor_results[label] = result
        check(
            f"raw modular closure {label}",
            result["rank"] == expected,
            f"rank={result['rank']}, orbit_support={result['orbit_support']}",
        )
    verify_external_factor_arithmetic()

    # 5. Validate, but do not trust, the producer artifact after recomputation.
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    check("artifact top-level schema", set(artifact) == {"meta", "data", "checks"})
    check("all producer checks passed", all(row.get("passed") is True for row in artifact["checks"]))
    headline = artifact["data"]["headline"]
    check(
        "artifact keeps n0 unresolved",
        headline["n0_status"]["exact_n0"] is None
        and headline["n0_status"]["status"] == "[UNRESOLVED]",
    )
    artifact_counterexamples = {
        row["graph"]: (row["dimension_Q"], row["quadratic_ceiling"])
        for row in headline["counterexamples"]
    }
    check(
        "artifact counterexamples match clean-room values",
        artifact_counterexamples == {"K_{2,3}": (44, 45), "K_{3,3}": (63, 66)},
    )
    artifact_anchors = {
        row["label"]: row["raw_generator_modular_closure"]["rank_Fp"]
        for row in artifact["data"]["closure_computations"]["grid_anchors"]
    }
    check(
        "artifact anchor ranks match clean-room values",
        artifact_anchors
        == {"open_2x2": 11, "open_2x3": 263, "open_2x4": 2952, "open_3x3": 8034},
    )
    source_hashes = artifact["meta"]["source_sha256"]
    check(
        "artifact source hashes",
        all((ROOT / relative).is_file() and digest(ROOT / relative) == value for relative, value in source_hashes.items()),
        f"files={len(source_hashes)}",
    )
    check(
        "resource limits",
        time.process_time() - STARTED < CPU_LIMIT and max_rss_bytes() < RSS_LIMIT,
        f"cpu={time.process_time()-STARTED:.3f}s, rss={max_rss_bytes()}",
    )

    if FAILURES:
        print("FAIL test_twogen_allsize: " + ", ".join(FAILURES))
        return 1
    print("PASS test_twogen_allsize")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
