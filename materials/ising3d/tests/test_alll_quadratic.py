"""Independent exact verifier for the all-L local-term quadratic no-go.

No producer module is imported. From raw 2 x L bonds, this test independently
rebuilds and evaluates every materialised bracket word at L=2..6. For the much
larger compact binomial family it checks the universal CAR swap identity, binary
Majorana-basis injectivity, every derived bilinear control, every possible
one-swap seed transition, and canonical paths of every possible swap length.
Those are the complete premises of the uniform Johnson-graph membership proof;
no stored conclusion is trusted.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from math import comb, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from ising.clifford import (  # noqa: E402
    dla_pauli_closure,
    pauli_x,
    symplectic_form,
    tfim_generators,
    zz,
)

ARTIFACT = ROOT / "results" / "algebra_growth" / "alll_quadratic.json"
CPU_BUDGET_SECONDS = 60.0


def ladder_bonds(L: int) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    for c in range(L - 1):
        bonds.extend(((c, c + 1), (L + c, L + c + 1)))
    bonds.extend((c, L + c) for c in range(L))
    return sorted(tuple(sorted(edge)) for edge in bonds)


def snake(L: int) -> list[int]:
    return list(range(L)) + list(range(2 * L - 1, L - 1, -1))


def majoranas(L: int) -> list[int]:
    n = 2 * L
    prefix = 0
    out: list[int] = []
    for site in snake(L):
        z = 1 << (n + site)
        out.extend((prefix | z, prefix | (1 << site) | z))
        prefix ^= 1 << site
    return out


def path_edges(L: int) -> list[int]:
    n = 2 * L
    path = snake(L)
    out: list[int] = []
    for q in range(4 * L - 1):
        if q % 2 == 0:
            out.append(pauli_x(n, path[q // 2]))
        else:
            j = q // 2
            out.append(zz(n, path[j], path[j + 1]))
    return out


def commutator(left: tuple[int, int], right: tuple[int, int], n: int) -> tuple[int, int]:
    """Independent ordered-Pauli evaluation of [left,right]."""

    lc, lv = left
    rc, rv = right
    mask = (1 << n) - 1
    lx, lz = lv & mask, (lv >> n) & mask
    rx, rz = rv & mask, (rv >> n) & mask
    sign_lr = -1 if ((lz & rx).bit_count() & 1) else 1
    sign_rl = -1 if ((rz & lx).bit_count() & 1) else 1
    return lc * rc * (sign_lr - sign_rl), lv ^ rv


def derive_h(a: int, b: int, edges: list[int], n: int) -> tuple[int, int]:
    assert 0 <= a < b < len(edges) + 1
    word = (1, edges[a])
    for q in range(a + 1, b):
        word = commutator((1, edges[q]), word, n)
        assert word[0] != 0, ("zero H word", a, b, q)
    return word


def gf2_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            lead = row.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = row
                break
            row ^= pivots[lead]
    return len(pivots)


def digest(values: list[int]) -> str:
    raw = "\n".join(format(value, "x") for value in values).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


def rung_support(L: int, column: int) -> list[int]:
    n = 2 * L
    support: list[int] = []
    for path_bond in range(column, n - column - 1):
        support.extend((2 * path_bond + 1, 2 * path_bond + 2))
    return support


def central_columns(L: int) -> list[int]:
    if L < 3:
        return []
    return [L // 2] if L % 2 else [L // 2 - 1, L // 2]


def orbit_column(L: int) -> int:
    return L // 2 if L % 2 else L // 2 - 1


def gaussian_dimension(m: int) -> int:
    return m * (2 * m - 1)


def min_modes(d: int) -> int:
    candidate = max(0, (1 + isqrt(1 + 8 * d)) // 4)
    while gaussian_dimension(candidate) < d:
        candidate += 1
    while candidate and gaussian_dimension(candidate - 1) >= d:
        candidate -= 1
    return candidate


def verify_size(L: int, stored: dict, deadline: float) -> set[int]:
    def tick() -> None:
        if time.process_time() > deadline:
            raise RuntimeError(f"NON-DECISIVE: verifier exceeded declared process-time budget at L={L}")

    n, N = 2 * L, 4 * L
    bonds = ladder_bonds(L)
    raw_generators, _, raw_n = tfim_generators(list(range(n)), bonds, n)
    assert raw_n == n
    local = set(raw_generators)
    assert len(local) == 5 * L - 2

    gamma = majoranas(L)
    edges = path_edges(L)
    assert gf2_rank(gamma) == N
    assert all(edge in local for edge in edges)
    for a in range(N):
        for b in range(N):
            assert symplectic_form(gamma[a], gamma[b], n) == (a != b)
    for q, edge in enumerate(edges):
        assert edge == gamma[q] ^ gamma[q + 1]

    h_cache: dict[tuple[int, int], tuple[int, int]] = {}
    path_values: list[int] = []
    for a in range(N - 1):
        for b in range(a + 1, N):
            tick()
            value = derive_h(a, b, edges, n)
            assert value[0] != 0 and value[1] == gamma[a] ^ gamma[b], (L, a, b)
            h_cache[a, b] = value
            path_values.append(value[1])
    assert len(set(path_values)) == comb(N, 2) == n * (2 * n - 1)

    closing = zz(n, snake(L)[0], snake(L)[-1])
    assert closing in local and closing not in set(path_values)
    assert len(set(path_values) | {closing}) == n * (2 * n - 1) + 1

    volume = 0
    for value in gamma:
        volume ^= value
    assert volume == (1 << n) - 1
    assert closing == volume ^ gamma[0] ^ gamma[-1]

    dual_values: list[int] = []
    for a in range(N - 1):
        for b in range(a + 1, N):
            tick()
            current = (1, closing)
            if a:
                current = commutator(h_cache[0, a], current, n)
                assert current[0] != 0, ("dual first", L, a, b)
            if b < N - 1:
                current = commutator(h_cache[b, N - 1], current, n)
                assert current[0] != 0, ("dual second", L, a, b)
            assert current[1] == volume ^ gamma[a] ^ gamma[b]
            dual_values.append(current[1])
    ring = set(path_values) | set(dual_values)
    assert len(ring) == 2 * n * (2 * n - 1)

    assert stored["path_bilinears"]["members_hex"] == [format(value, "x") for value in path_values]
    assert stored["path_bilinears"]["sha256"] == digest(path_values)
    assert stored["dual_bilinears"]["members_hex"] == [format(value, "x") for value in dual_values]
    assert stored["dual_bilinears"]["sha256"] == digest(dual_values)
    assert int(stored["closing_rung"]["packed_hex"], 16) == closing

    materialised = set(ring)
    rebuilt_hypercubes: list[dict] = []
    for column in central_columns(L):
        tick()
        support = rung_support(L, column)
        support_set = set(support)
        outside = [index for index in range(N) if index not in support_set]
        rank = min(len(support), len(outside))
        pairs = list(zip(support[:rank], outside[:rank]))
        seed = zz(n, column, L + column)
        assert seed in local
        seed_from_gamma = 0
        for index in support:
            seed_from_gamma ^= gamma[index]
        assert seed == seed_from_gamma

        values: list[int] = []
        for mask in range(1 << rank):
            tick()
            current = (1, seed)
            expected = seed
            for bit, (inside, out) in enumerate(pairs):
                if (mask >> bit) & 1:
                    control = h_cache[tuple(sorted((inside, out)))]
                    current = commutator(control, current, n)
                    assert current[0] != 0, ("hypercube", L, column, mask, bit)
                    expected ^= gamma[inside] ^ gamma[out]
            assert current[1] == expected
            values.append(current[1])
        assert len(set(values)) == 1 << rank
        assert not (materialised & set(values))
        materialised.update(values)
        rebuilt_hypercubes.append(
            {
                "column": column + 1,
                "grade": len(support),
                "rank": rank,
                "seed_majorana_indices": [index + 1 for index in support],
                "swap_pairs_majorana_indices": [[a + 1, b + 1] for a, b in pairs],
                "count": len(values),
                "members_hex": [format(value, "x") for value in values],
                "sha256": digest(values),
            }
        )
    assert stored["central_hypercubes"] == rebuilt_hypercubes
    assert stored["materialised_family_sha256"] == digest(sorted(materialised))

    # Compact central orbit: re-derive its universal schema, not a stored count.
    column = orbit_column(L)
    support = set(rung_support(L, column))
    outside = set(range(N)) - support
    grade = len(support)
    seed = zz(n, column, L + column)
    expected_seed = 0
    for index in support:
        expected_seed ^= gamma[index]
    assert seed == expected_seed and seed in local

    # Every one-index replacement is a nonzero exact bracket with an already
    # derived H[a,b]. This is the universal Johnson edge rule.
    for inside in sorted(support):
        for out in sorted(outside):
            tick()
            control = h_cache[tuple(sorted((inside, out)))]
            value = commutator(control, (1, seed), n)
            assert value[0] != 0
            assert value[1] == seed ^ gamma[inside] ^ gamma[out]

    # Independently exercise the deterministic word at every possible Johnson
    # distance. Connectivity then covers every same-grade target.
    seed_sorted, outside_sorted = sorted(support), sorted(outside)
    for distance in range(min(len(seed_sorted), len(outside_sorted)) + 1):
        tick()
        removed = seed_sorted[:distance]
        added = outside_sorted[:distance]
        target = (support - set(removed)) | set(added)
        current = (1, seed)
        for inside, out in zip(removed, added):
            current = commutator(h_cache[tuple(sorted((inside, out)))], current, n)
            assert current[0] != 0
        expected = 0
        for index in target:
            expected ^= gamma[index]
        assert current[1] == expected

    compact_count = comb(N, grade) // 2 if grade == n else comb(N, grade)
    compact_stored = stored["compact_central_orbit"]
    assert compact_stored["column"] == column + 1
    assert compact_stored["grade"] == grade
    assert compact_stored["claimed_count"] == compact_count
    assert compact_stored["single_swap_schema_checks"] == len(support) * len(outside)

    ceiling = n * (2 * n - 1)
    materialised_extra = 0 if L == 2 else ((1 << n) if L % 2 else (1 << (n - 1)))
    strongest = 2 * ceiling if L == 2 else 2 * ceiling + compact_count
    assert stored["physical_gaussian_ceiling"] == ceiling
    assert stored["polynomial_family_bound"] == ceiling + 1
    assert stored["ring_family_bound"] == 2 * ceiling
    assert stored["materialised_hypercube_count"] == materialised_extra
    assert stored["materialised_family_bound"] == 2 * ceiling + materialised_extra
    assert stored["strongest_dimension_bound"] == strongest
    modes = min_modes(strongest)
    assert stored["minimum_modes_from_strongest_bound"] == modes
    assert gaussian_dimension(modes - 1) < strongest <= gaussian_dimension(modes)
    assert 5 * (4 * L + 1) * compact_count >= (1 << (4 * L + 1))
    return materialised


def main() -> None:
    started = time.process_time()
    deadline = started + CPU_BUDGET_SECONDS
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert set(artifact) == {"meta", "data"}
    assert artifact["meta"]["status"] == "PASS"
    assert all(item["passed"] for item in artifact["data"]["checks"])
    stored_rows = {row["L"]: row for row in artifact["data"]["sizes"]}
    assert set(stored_rows) == set(range(2, 9))

    rebuilt: dict[int, set[int]] = {}
    for L in range(2, 7):
        rebuilt[L] = verify_size(L, stored_rows[L], deadline)
        print(
            f"L={L}: independently evaluated {len(rebuilt[L])} stored materialised members; "
            f"compact bound={stored_rows[L]['compact_central_orbit_bound']}"
        )

    # Full raw-source closures are cheap at L=2,3 and independently confirm that
    # every constructed Pauli lies in the actual local-term DLA.
    for L, exact_dimension in ((2, 56), (3, 1056)):
        n = 2 * L
        generators, _, _ = tfim_generators(list(range(n)), ladder_bonds(L), n)
        closure, saturated = dla_pauli_closure(generators, n, max_size=1 << 16)
        closure_set = set(closure)
        assert saturated and len(closure_set) == exact_dimension
        assert rebuilt[L] <= closure_set

    # Independently recompute the finite exact-dimension mode refinements.
    expected_exact = {2: (56, 6), 3: (1056, 24), 4: (16256, 91), 5: (262656, 363)}
    refinements = {row["L"]: row for row in artifact["data"]["certified_exact_finite_refinements"]}
    assert set(refinements) == set(expected_exact)
    for L, (dimension, modes) in expected_exact.items():
        row = refinements[L]
        assert row["exact_characteristic_zero_dimension"] == dimension
        assert row["minimum_modes_from_exact_dimension"] == min_modes(dimension) == modes

    print("PASS")


if __name__ == "__main__":
    main()
