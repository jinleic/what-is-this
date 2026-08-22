"""Exact bracket witnesses for all-L local-term DLA lower bounds.

The open 2 x L ladder is ordered along a perimeter Hamiltonian path. In the
associated Jordan--Wigner Majorana basis, consecutive Majorana bilinears are
literal local terms. This module constructs and evaluates:

* every path bilinear (Clifford grade 2);
* its closing-rung orbit (Clifford grade 4L-2);
* finite central-rung Johnson hypercubes, exhaustively through L=8; and
* a compact Johnson-orbit schema for the larger binomial lower bound.

Every materialised member is evaluated as an actual nested bracket. Arithmetic
is exact integer Pauli arithmetic; phases are represented by nonzero integer
coefficients and packed Pauli vectors use ising.clifford's (x | z) convention.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from dataclasses import dataclass
from math import comb
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from ising.clifford import pauli_x, symplectic_form, zz  # noqa: E402

CPU_BUDGET_SECONDS = 60.0


@dataclass(frozen=True)
class Monomial:
    """An exact scalar multiple ``coefficient * Q_packed``."""

    coefficient: int
    packed: int


class NonDecisiveBudget(RuntimeError):
    """A process-time wall is non-decisive and carries no mathematical verdict."""


def snake_path(L: int) -> tuple[int, ...]:
    """Top row left-to-right, then bottom row right-to-left."""

    if L < 2:
        raise ValueError("the theorem starts at L=2")
    return tuple(range(L)) + tuple(range(2 * L - 1, L - 1, -1))


def site_coordinate(site: int, L: int) -> tuple[str, int]:
    if not 0 <= site < 2 * L:
        raise ValueError(site)
    return ("top", site + 1) if site < L else ("bottom", site - L + 1)


def ladder_bonds(L: int) -> tuple[tuple[int, int], ...]:
    bonds: set[tuple[int, int]] = set()
    for column in range(L - 1):
        bonds.add((column, column + 1))
        bonds.add((L + column, L + column + 1))
    for column in range(L):
        bonds.add((column, L + column))
    return tuple(sorted(tuple(sorted(edge)) for edge in bonds))


def local_term_set(L: int) -> set[int]:
    n = 2 * L
    terms = {pauli_x(n, site) for site in range(n)}
    terms.update(zz(n, u, v) for u, v in ladder_bonds(L))
    return terms


def canonical_majoranas(L: int) -> tuple[int, ...]:
    """Packed canonical Majoranas in snake-path order, with phases suppressed."""

    n = 2 * L
    prefix_x = 0
    values: list[int] = []
    for site in snake_path(L):
        z_here = 1 << (n + site)
        values.append(prefix_x | z_here)
        values.append(prefix_x | (1 << site) | z_here)
        prefix_x ^= 1 << site
    return tuple(values)


def path_edge_terms(L: int) -> tuple[int, ...]:
    """The 4L-1 local terms representing consecutive Majorana bilinears."""

    n = 2 * L
    path = snake_path(L)
    edges: list[int] = []
    for majorana_edge in range(4 * L - 1):
        if majorana_edge % 2 == 0:
            edges.append(pauli_x(n, path[majorana_edge // 2]))
        else:
            j = majorana_edge // 2
            edges.append(zz(n, path[j], path[j + 1]))
    return tuple(edges)


def bracket(left: Monomial, right: Monomial, n: int) -> Monomial:
    """Evaluate ``[left, right]`` exactly in the ordered Pauli convention."""

    mask = (1 << n) - 1
    a_left, b_left = left.packed & mask, (left.packed >> n) & mask
    a_right, b_right = right.packed & mask, (right.packed >> n) & mask
    left_phase = -1 if ((b_left & a_right).bit_count() & 1) else 1
    right_phase = -1 if ((b_right & a_left).bit_count() & 1) else 1
    factor = left_phase - right_phase
    return Monomial(left.coefficient * right.coefficient * factor, left.packed ^ right.packed)


def path_bilinear(
    L: int,
    a: int,
    b: int,
    edges: tuple[int, ...] | None = None,
) -> Monomial:
    """Evaluate H[a,b], with H[a,b]=[G[b-1],H[a,b-1]]."""

    N = 4 * L
    if not 0 <= a < b < N:
        raise ValueError((a, b, N))
    if edges is None:
        edges = path_edge_terms(L)
    current = Monomial(1, edges[a])
    for edge in range(a + 1, b):
        current = bracket(Monomial(1, edges[edge]), current, 2 * L)
        if current.coefficient == 0:
            raise AssertionError(f"zero path bracket at L={L}, H[{a + 1},{b + 1}], edge={edge + 1}")
    return current


def path_bilinear_cache(L: int) -> dict[tuple[int, int], Monomial]:
    edges = path_edge_terms(L)
    return {
        (a, b): path_bilinear(L, a, b, edges)
        for a in range(4 * L - 1)
        for b in range(a + 1, 4 * L)
    }


def closing_rung(L: int) -> int:
    path = snake_path(L)
    return zz(2 * L, path[0], path[-1])


def dual_bilinear(
    L: int,
    a: int,
    b: int,
    cache: dict[tuple[int, int], Monomial],
) -> Monomial:
    """Derive the grade-(4L-2) monomial with Majorana holes ``{a,b}``."""

    N = 4 * L
    if not 0 <= a < b < N:
        raise ValueError((a, b, N))
    current = Monomial(1, closing_rung(L))
    if a > 0:
        current = bracket(cache[(0, a)], current, 2 * L)
        if current.coefficient == 0:
            raise AssertionError(f"zero first dual bracket at L={L}, holes={a + 1,b + 1}")
    if b < N - 1:
        current = bracket(cache[(b, N - 1)], current, 2 * L)
        if current.coefficient == 0:
            raise AssertionError(f"zero second dual bracket at L={L}, holes={a + 1,b + 1}")
    return current


def central_columns(L: int) -> tuple[int, ...]:
    """Zero-based rung columns used by the exhaustively materialised hypercubes."""

    if L < 3:
        return ()
    if L % 2:
        return (L // 2,)
    return (L // 2 - 1, L // 2)


def orbit_column(L: int) -> int:
    """One central rung for the compact binomial family (grade n or n+2)."""

    return L // 2 if L % 2 else L // 2 - 1


def rung_majorana_support(L: int, column: int) -> tuple[int, ...]:
    """Majorana support of the physical rung in ``column`` (all indices zero-based)."""

    if not 0 <= column < L:
        raise ValueError(column)
    n = 2 * L
    left_position = column
    right_position = n - column - 1
    support: list[int] = []
    for path_bond in range(left_position, right_position):
        support.extend((2 * path_bond + 1, 2 * path_bond + 2))
    return tuple(support)


def rung_term(L: int, column: int) -> int:
    return zz(2 * L, column, L + column)


def hypercube_spec(L: int, column: int) -> dict:
    N = 4 * L
    support = rung_majorana_support(L, column)
    support_set = set(support)
    complement = tuple(index for index in range(N) if index not in support_set)
    rank = min(len(support), len(complement))
    pairs = tuple(zip(support[:rank], complement[:rank]))
    return {
        "column_zero_based": column,
        "grade": len(support),
        "rank": rank,
        "support": support,
        "pairs": pairs,
    }


def hypercube_member(
    L: int,
    column: int,
    mask: int,
    cache: dict[tuple[int, int], Monomial],
    spec: dict | None = None,
) -> Monomial:
    """Evaluate one nested-bracket member of a central-rung Johnson hypercube."""

    if spec is None:
        spec = hypercube_spec(L, column)
    rank = int(spec["rank"])
    if not 0 <= mask < (1 << rank):
        raise ValueError((mask, rank))
    current = Monomial(1, rung_term(L, column))
    for bit, (inside, outside) in enumerate(spec["pairs"]):
        if (mask >> bit) & 1:
            a, b = sorted((inside, outside))
            current = bracket(cache[(a, b)], current, 2 * L)
            if current.coefficient == 0:
                raise AssertionError(
                    f"zero hypercube bracket L={L}, column={column + 1}, mask={mask}, bit={bit}"
                )
    return current


def johnson_member(
    L: int,
    column: int,
    target: Iterable[int],
    cache: dict[tuple[int, int], Monomial],
) -> Monomial:
    """Derive an arbitrary same-grade target by deterministic remove/add swaps."""

    seed = set(rung_majorana_support(L, column))
    target_set = set(int(index) for index in target)
    if len(target_set) != len(seed) or not target_set <= set(range(4 * L)):
        raise ValueError("target must be a same-grade Majorana subset")
    removed = sorted(seed - target_set)
    added = sorted(target_set - seed)
    if len(removed) != len(added):
        raise AssertionError("same-grade subsets must have equally many removals and additions")
    current = Monomial(1, rung_term(L, column))
    for inside, outside in zip(removed, added):
        a, b = sorted((inside, outside))
        current = bracket(cache[(a, b)], current, 2 * L)
        if current.coefficient == 0:
            raise AssertionError(
                f"zero Johnson swap L={L}, column={column + 1}, inside={inside + 1}, outside={outside + 1}"
            )
    return current


def compact_orbit_count(L: int) -> int:
    """Conservative claimed count: half of middle grade, full non-middle grade."""

    grade = len(rung_majorana_support(L, orbit_column(L)))
    count = comb(4 * L, grade)
    return count // 2 if grade == 2 * L else count


def gf2_rank(values: Iterable[int]) -> int:
    basis: dict[int, int] = {}
    for value in values:
        row = int(value)
        while row:
            lead = row.bit_length() - 1
            if lead not in basis:
                basis[lead] = row
                break
            row ^= basis[lead]
    return len(basis)


def digest_ints(values: Iterable[int]) -> str:
    payload = "\n".join(format(int(value), "x") for value in values).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def certify_size(L: int, deadline: float | None = None) -> dict:
    """Evaluate every materialised witness for one size and verify the compact schema."""

    def tick() -> None:
        if deadline is not None and time.process_time() > deadline:
            raise NonDecisiveBudget(f"NON-DECISIVE: e160 exceeded its declared process-time budget at L={L}")

    n = 2 * L
    N = 2 * n
    path = snake_path(L)
    majoranas = canonical_majoranas(L)
    edges = path_edge_terms(L)
    local = local_term_set(L)
    if not all(edge in local for edge in edges):
        raise AssertionError(f"Hamiltonian-path edge is not a local term at L={L}")
    if gf2_rank(majoranas) != N:
        raise AssertionError(f"canonical Majoranas do not form a binary basis at L={L}")
    for a in range(N):
        for b in range(N):
            expected = 0 if a == b else 1
            if symplectic_form(majoranas[a], majoranas[b], n) != expected:
                raise AssertionError(f"Majorana CAR failed at L={L}, pair={a + 1,b + 1}")
    for index, edge in enumerate(edges):
        if edge != (majoranas[index] ^ majoranas[index + 1]):
            raise AssertionError(f"Majorana-edge identity failed at L={L}, edge={index + 1}")

    cache = path_bilinear_cache(L)
    path_values: list[int] = []
    for (a, b), value in cache.items():
        tick()
        if value.coefficient == 0 or value.packed != (majoranas[a] ^ majoranas[b]):
            raise AssertionError(f"path witness failed at L={L}, pair={a + 1,b + 1}")
        path_values.append(value.packed)

    close = closing_rung(L)
    if close not in local or close in set(path_values):
        raise AssertionError(f"closing-rung baseline witness failed at L={L}")

    all_majoranas = 0
    for value in majoranas:
        all_majoranas ^= value
    if all_majoranas != (1 << n) - 1:
        raise AssertionError(f"Majorana volume element is not global X at L={L}")
    if close != (all_majoranas ^ majoranas[0] ^ majoranas[-1]):
        raise AssertionError(f"closing rung has wrong Clifford support at L={L}")

    dual_values: list[int] = []
    for a in range(N - 1):
        for b in range(a + 1, N):
            tick()
            value = dual_bilinear(L, a, b, cache)
            expected = all_majoranas ^ majoranas[a] ^ majoranas[b]
            if value.coefficient == 0 or value.packed != expected:
                raise AssertionError(f"dual witness failed at L={L}, holes={a + 1,b + 1}")
            dual_values.append(value.packed)

    path_set, dual_set = set(path_values), set(dual_values)
    expected_orbit = comb(N, 2)
    if len(path_set) != expected_orbit or len(dual_set) != expected_orbit or path_set & dual_set:
        raise AssertionError(f"ring-orbit independence failed at L={L}")

    central_records: list[dict] = []
    materialised_set = path_set | dual_set
    for column in central_columns(L):
        tick()
        spec = hypercube_spec(L, column)
        seed_from_majoranas = 0
        for index in spec["support"]:
            seed_from_majoranas ^= majoranas[index]
        if seed_from_majoranas != rung_term(L, column) or rung_term(L, column) not in local:
            raise AssertionError(f"central rung seed failed at L={L}, column={column + 1}")
        values: list[int] = []
        for mask in range(1 << spec["rank"]):
            tick()
            member = hypercube_member(L, column, mask, cache, spec)
            expected = seed_from_majoranas
            for bit, (inside, outside) in enumerate(spec["pairs"]):
                if (mask >> bit) & 1:
                    expected ^= majoranas[inside] ^ majoranas[outside]
            if member.coefficient == 0 or member.packed != expected:
                raise AssertionError(
                    f"hypercube witness failed at L={L}, column={column + 1}, mask={mask}"
                )
            values.append(member.packed)
        value_set = set(values)
        if len(value_set) != 1 << spec["rank"] or materialised_set & value_set:
            raise AssertionError(f"hypercube independence failed at L={L}, column={column + 1}")
        materialised_set |= value_set
        central_records.append(
            {
                "column": column + 1,
                "grade": spec["grade"],
                "rank": spec["rank"],
                "seed_majorana_indices": [index + 1 for index in spec["support"]],
                "swap_pairs_majorana_indices": [
                    [inside + 1, outside + 1] for inside, outside in spec["pairs"]
                ],
                "count": len(values),
                "members_hex": [format(value, "x") for value in values],
                "sha256": digest_ints(values),
            }
        )

    central_count = sum(record["count"] for record in central_records)
    expected_central = 0 if L == 2 else ((1 << n) if L % 2 else (1 << (n - 1)))
    expected_materialised = 2 * n * (2 * n - 1) + expected_central
    if central_count != expected_central or len(materialised_set) != expected_materialised:
        raise AssertionError(f"materialised-family count failed at L={L}")

    # Compact binomial family: verify its seed and the universal one-index-swap rule.
    column = orbit_column(L)
    support = set(rung_majorana_support(L, column))
    seed_from_majoranas = 0
    for index in support:
        seed_from_majoranas ^= majoranas[index]
    if seed_from_majoranas != rung_term(L, column):
        raise AssertionError(f"compact orbit seed failed at L={L}")
    complement = set(range(N)) - support
    for inside in support:
        for outside in complement:
            tick()
            control = cache[tuple(sorted((inside, outside)))]
            swapped = bracket(control, Monomial(1, seed_from_majoranas), n)
            expected = seed_from_majoranas ^ majoranas[inside] ^ majoranas[outside]
            if swapped.coefficient == 0 or swapped.packed != expected:
                raise AssertionError(
                    f"compact orbit swap schema failed at L={L}, pair={inside + 1,outside + 1}"
                )

    return {
        "L": L,
        "physical_qubits_n": n,
        "majoranas_N": N,
        "local_term_count": len(local),
        "path_sites": [
            {"site": site, "row": site_coordinate(site, L)[0], "column": site_coordinate(site, L)[1]}
            for site in path
        ],
        "physical_gaussian_ceiling": n * (2 * n - 1),
        "path_bilinears": {
            "count": len(path_values),
            "members_hex": [format(value, "x") for value in path_values],
            "sha256": digest_ints(path_values),
        },
        "closing_rung": {
            "endpoints": [["top", 1], ["bottom", 1]],
            "packed_hex": format(close, "x"),
            "outside_path_family": True,
        },
        "polynomial_family_count": len(path_set | {close}),
        "dual_bilinears": {
            "count": len(dual_values),
            "members_hex": [format(value, "x") for value in dual_values],
            "sha256": digest_ints(dual_values),
        },
        "ring_family_count": len(path_set | dual_set),
        "central_hypercubes": central_records,
        "central_hypercube_count": central_count,
        "materialised_family_count": len(materialised_set),
        "materialised_family_sha256": digest_ints(sorted(materialised_set)),
        "compact_central_orbit": {
            "column": column + 1,
            "grade": len(support),
            "middle_grade": len(support) == n,
            "selection_rule": (
                "one lexicographic representative from each complement pair"
                if len(support) == n
                else "all same-grade Majorana subsets"
            ),
            "claimed_count": compact_orbit_count(L),
            "seed_majorana_indices": [index + 1 for index in sorted(support)],
            "seed_packed_hex": format(seed_from_majoranas, "x"),
            "single_swap_schema_checks": len(support) * len(complement),
        },
    }


def main() -> None:
    started = time.process_time()
    deadline = started + CPU_BUDGET_SECONDS
    records = [certify_size(L, deadline) for L in range(2, 9)]
    for record in records:
        print(
            f"L={record['L']}: polynomial={record['polynomial_family_count']} "
            f"ring={record['ring_family_count']} materialised={record['materialised_family_count']} "
            f"compact={record['compact_central_orbit']['claimed_count']}"
        )
    summary = json.dumps(
        [
            (
                row["L"],
                row["materialised_family_count"],
                row["materialised_family_sha256"],
                row["compact_central_orbit"]["claimed_count"],
            )
            for row in records
        ],
        separators=(",", ":"),
    )
    print(f"certificate_sha256={hashlib.sha256(summary.encode('ascii')).hexdigest()}")
    print("PASS")


if __name__ == "__main__":
    main()
