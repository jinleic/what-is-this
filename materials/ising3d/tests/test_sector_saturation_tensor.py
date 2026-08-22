#!/usr/bin/env python3
"""Clean-room verifier for the tensor-cut sector-saturation obstruction.

It does not import the producer.  From raw ladder configurations it rebuilds
all exact orbit spaces, the cut-extension map, the 4B action, every stored
invariant-annihilator certificate, the modular cyclic lower ranks, and the
local new-rung identities.  It deliberately flips one certificate sign as a
negative control; that mutation must cease to be 4B-invariant.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "sector_saturation_tensor.json"
P = 2_147_483_647

CHECKS = 0
FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if condition:
        print(f"  ok   {name} {detail}")
    else:
        print(f"  FAIL {name} {detail}")
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Independent ladder, symmetry, and orbit construction (pure Python).
# ---------------------------------------------------------------------------
def tau(config: int, L: int) -> int:
    answer = 0
    for rung in range(L):
        top = (config >> (2 * rung)) & 1
        bottom = (config >> (2 * rung + 1)) & 1
        answer |= (bottom << (2 * rung)) | (top << (2 * rung + 1))
    return answer


def rho(config: int, L: int) -> int:
    answer = 0
    for rung in range(L):
        state = (config >> (2 * rung)) & 3
        answer |= state << (2 * (L - 1 - rung))
    return answer


def group_images(config: int, L: int) -> tuple[int, ...]:
    return tuple(sorted({config, tau(config, L), rho(config, L), tau(rho(config, L), L)}))


def even_orbits(L: int) -> tuple[list[int], list[tuple[int, ...]]]:
    orbit_of = [-1] * (1 << (2 * L))
    members: list[tuple[int, ...]] = []
    for config in range(1 << (2 * L)):
        if orbit_of[config] >= 0 or config.bit_count() & 1:
            continue
        orbit = group_images(config, L)
        index = len(members)
        for member in orbit:
            orbit_of[member] = index
        members.append(orbit)
    return orbit_of, members


def edges(L: int) -> list[tuple[int, int]]:
    answer = [(2 * rung, 2 * rung + 1) for rung in range(L)]
    for rung in range(L - 1):
        answer.extend(((2 * rung, 2 * rung + 2), (2 * rung + 1, 2 * rung + 3)))
    return answer


def burnside_K_dim(L: int) -> int:
    return (2 ** (2 * L - 1) + 3 * 2**L) // 4


def four_B_rows(L: int, orbit_of: list[int], members: list[tuple[int, ...]], modulus: int | None) -> list[dict[int, int]]:
    """The exact integer matrix 4B in the orbit-sum basis."""
    answer: list[dict[int, int]] = []
    for orbit in members:
        counts: dict[int, int] = defaultdict(int)
        for config in orbit:
            for u, v in edges(L):
                counts[orbit_of[config ^ ((1 << u) | (1 << v))]] += 1
        row: dict[int, int] = {}
        for target, count in counts.items():
            target_size = len(members[target])
            if count % target_size:
                raise AssertionError("4B orbit coefficient was not integral")
            coefficient = 4 * (count // target_size)
            if coefficient:
                row[target] = coefficient % modulus if modulus else coefficient
        answer.append(row)
    return answer


def multiply(vector: dict[int, int], rows: list[dict[int, int]], modulus: int | None) -> dict[int, int]:
    answer: dict[int, int] = defaultdict(int)
    for source, coefficient in vector.items():
        for target, entry in rows[source].items():
            answer[target] += coefficient * entry
    if modulus:
        return {index: value % modulus for index, value in answer.items() if value % modulus}
    return {index: value for index, value in answer.items() if value}


class ModEchelon:
    def __init__(self) -> None:
        self.rows: dict[int, dict[int, int]] = {}

    def add(self, vector: dict[int, int]) -> int | None:
        row = {index: value % P for index, value in vector.items() if value % P}
        while row:
            pivot = max(row)
            old = self.rows.get(pivot)
            if old is None:
                inverse = pow(row[pivot], P - 2, P)
                self.rows[pivot] = {index: value * inverse % P for index, value in row.items()}
                return pivot
            scale = row[pivot]
            for index, value in old.items():
                new_value = (row.get(index, 0) - scale * value) % P
                if new_value:
                    row[index] = new_value
                elif index in row:
                    del row[index]
        return None


class QSpan:
    """Sparse exact-Q row span; the certificate bases are deliberately small."""
    def __init__(self) -> None:
        self.rows: dict[int, dict[int, Fraction]] = {}

    def residual(self, vector: dict[int, int | Fraction]) -> dict[int, Fraction]:
        row = {index: value if isinstance(value, Fraction) else Fraction(value)
               for index, value in vector.items() if value}
        while row:
            pivot = max(row)
            old = self.rows.get(pivot)
            if old is None:
                return row
            scale = row[pivot]
            for index, value in old.items():
                new_value = row.get(index, Fraction(0)) - scale * value
                if new_value:
                    row[index] = new_value
                elif index in row:
                    del row[index]
        return {}

    def add(self, vector: dict[int, int | Fraction]) -> bool:
        row = self.residual(vector)
        if not row:
            return False
        pivot = max(row)
        pivot_value = row[pivot]
        self.rows[pivot] = {index: value / pivot_value for index, value in row.items()}
        return True

    @property
    def rank(self) -> int:
        return len(self.rows)


def extension_rows(L: int, old_orbit_of: list[int], old_members: list[tuple[int, ...]],
                   new_orbit_of: list[int], new_members: list[tuple[int, ...]]) -> list[dict[int, int]]:
    """Pi_{L+1}(K_L tensor span{00,11}) in orbit-sum coordinates."""
    answer: list[dict[int, int]] = []
    for old_orbit in old_members:
        for new_rung in (0, 3):
            counts: dict[int, int] = defaultdict(int)
            for old_config in old_orbit:
                full_config = old_config | (new_rung << (2 * L))
                for image in group_images(full_config, L + 1):
                    counts[new_orbit_of[image]] += 1
            row: dict[int, int] = {}
            for target, count in counts.items():
                target_size = len(new_members[target])
                if count % target_size:
                    raise AssertionError("extension orbit coefficient was not integral")
                coefficient = count // target_size
                if coefficient:
                    row[target] = coefficient
            answer.append(row)
    return answer


def end_block(config: int, L: int) -> str:
    left = (config >> 0) & 3
    right = (config >> (2 * (L - 1))) & 3
    left_kind = "M" if left in (1, 2) else "E"
    right_kind = "M" if right in (1, 2) else "E"
    return "".join(sorted((left_kind, right_kind)))


def extension_record(L: int, claimed: dict) -> None:
    old_orbit_of, old_members = even_orbits(L)
    new_orbit_of, new_members = even_orbits(L + 1)
    rows = extension_rows(L, old_orbit_of, old_members, new_orbit_of, new_members)
    span = QSpan()
    for row in rows:
        span.add(row)
    free = [index for index in range(len(new_members)) if index not in span.rows]
    reps = [new_members[index][0] for index in free]
    histogram = Counter(f"{end_block(new_members[index][0], L + 1)}:k{new_members[index][0].bit_count()}" for index in free)
    check(f"K_dim_L{L}", len(old_members) == burnside_K_dim(L), f"{len(old_members)}")
    check(f"K_dim_L{L+1}", len(new_members) == burnside_K_dim(L + 1), f"{len(new_members)}")
    check(f"extension_domain_L{L}", len(rows) == claimed["extension_domain_dim"], str(len(rows)))
    check(f"extension_rank_Q_L{L}", span.rank == claimed["extension_rank_Q"], str(span.rank))
    check(f"correction_dim_L{L}", len(free) == claimed["correction_dim"], str(len(free)))
    check(f"correction_direct_sum_L{L}", span.rank + len(free) == len(new_members), "pivot/free decomposition")
    check(f"correction_reps_L{L}", reps == claimed["correction_representatives"], f"n={len(reps)}")
    check(f"correction_blocks_L{L}", dict(sorted(histogram.items())) == claimed["correction_block_histogram"], f"n={len(histogram)}")


def cyclic_rank_mod_p(L: int) -> int:
    orbit_of, members = even_orbits(L)
    B4 = four_B_rows(L, orbit_of, members, P)
    A = [(2 * L - 2 * orbit[0].bit_count()) % P for orbit in members]
    echelon = ModEchelon()
    pivot = echelon.add({orbit_of[0]: 1})
    assert pivot is not None
    todo = [pivot]
    processed: set[int] = set()
    while todo:
        pivot = todo.pop()
        if pivot in processed:
            continue
        processed.add(pivot)
        vector = echelon.rows[pivot]
        images = (
            {index: A[index] * coefficient % P for index, coefficient in vector.items()},
            multiply(vector, B4, P),
        )
        for image in images:
            new = echelon.add(image)
            if new is not None:
                todo.append(new)
    return len(echelon.rows)


def certificate_record(record: dict) -> tuple[list[dict[int, int]], list[tuple[int, ...]], list[int]]:
    L = record["L"]
    orbit_of, members = even_orbits(L)
    by_rep = {orbit[0]: index for index, orbit in enumerate(members)}
    basis: list[dict[int, int]] = []
    canonical_reps = True
    for stored_vector in record["basis"]:
        vector: dict[int, int] = {}
        for term in stored_vector:
            rep = term["representative"]
            canonical_reps = canonical_reps and rep in by_rep
            if rep not in by_rep:
                continue
            index = by_rep[rep]
            vector[index] = vector.get(index, 0) + term["coefficient"]
        basis.append({index: coefficient for index, coefficient in vector.items() if coefficient})
    span = QSpan()
    for vector in basis:
        span.add(vector)
    B4 = four_B_rows(L, orbit_of, members, None)
    A_sectors = [{members[index][0].bit_count() for index in vector} for vector in basis]
    sector_histogram = Counter(next(iter(sectors)) for sectors in A_sectors if len(sectors) == 1)
    invariant = all(not span.residual(multiply(vector, B4, None)) for vector in basis)
    check(f"certificate_reps_L{L}", canonical_reps, f"{sum(len(vector) for vector in basis)} terms")
    check(f"certificate_basis_dim_L{L}", span.rank == record["unreachable_block_dim"], str(span.rank))
    check(f"certificate_sector_L{L}", all(len(sectors) == 1 for sectors in A_sectors), str(dict(sorted(sector_histogram.items()))))
    check(f"certificate_vacuum_L{L}", all(orbit_of[0] not in vector for vector in basis), "psi orthogonal")
    check(f"certificate_4B_invariant_L{L}", invariant, "exact Q membership")
    q_upper = len(members) - span.rank
    rank_p = cyclic_rank_mod_p(L)
    check(f"cyclic_modular_rank_L{L}", rank_p == record["cyclic_mod_p_rank"], str(rank_p))
    check(f"cyclic_rank_sandwich_L{L}", rank_p == q_upper == record["cyclic_Q_dim"], f"{rank_p} <= dim_Q <= {q_upper}")
    return basis, members, orbit_of


# ---------------------------------------------------------------------------
# Local new-rung action and the rho-breaking edge-local witness.
# ---------------------------------------------------------------------------
def flip_if_mixed(config: int, u: int, v: int) -> int | None:
    if ((config >> u) & 1) == ((config >> v) & 1):
        return None
    return config ^ ((1 << u) | (1 << v))


def create_if_empty(config: int, u: int, v: int) -> int | None:
    if ((config >> u) & 1) or ((config >> v) & 1):
        return None
    return config ^ ((1 << u) | (1 << v))


def local_action_checks() -> None:
    for Lold in (3, 4):
        top_old, bottom_old = 2 * (Lold - 1), 2 * (Lold - 1) + 1
        top_new, bottom_new = 2 * Lold, 2 * Lold + 1
        pair_ok = True
        top_ok = True
        bottom_ok = True
        count = 0
        for old in range(1 << (2 * Lold)):
            count += 1
            full = old
            pair_ok = pair_ok and create_if_empty(full, top_new, bottom_new) == old | (3 << (2 * Lold))
            top = flip_if_mixed(full, top_old, top_new)
            expected_top = None if not ((old >> top_old) & 1) else (old ^ (1 << top_old)) | (1 << top_new)
            top_ok = top_ok and top == expected_top
            bottom = flip_if_mixed(full, bottom_old, bottom_new)
            expected_bottom = None if not ((old >> bottom_old) & 1) else (old ^ (1 << bottom_old)) | (1 << bottom_new)
            bottom_ok = bottom_ok and bottom == expected_bottom
        check(f"pair_create_all_L{Lold}", pair_ok, f"{count} raw configurations")
        check(f"rail_top_hop_all_L{Lold}", top_ok, f"{count} raw configurations")
        check(f"rail_bottom_hop_all_L{Lold}", bottom_ok, f"{count} raw configurations")


def apply_new_F_to_orbit(orbit: tuple[int, ...], Lold: int) -> dict[int, int]:
    Lnew = Lold + 1
    top_old, bottom_old = 2 * (Lold - 1), 2 * (Lold - 1) + 1
    top_new, bottom_new = 2 * Lold, 2 * Lold + 1
    new_edges = ((top_old, top_new), (bottom_old, bottom_new), (top_new, bottom_new))
    answer: dict[int, int] = defaultdict(int)
    for config in orbit:
        for u, v in new_edges:
            image = flip_if_mixed(config, u, v)
            if image is not None:
                answer[image] += 1
    return {config: coefficient for config, coefficient in answer.items() if coefficient}


def cut_witness(Lold: int) -> tuple[int, int]:
    _, members = even_orbits(Lold + 1)
    for orbit in members:
        image = apply_new_F_to_orbit(orbit, Lold)
        if not image:
            continue
        transformed: dict[int, int] = defaultdict(int)
        for config, coefficient in image.items():
            transformed[rho(config, Lold + 1)] += coefficient
        if dict(transformed) != image:
            difference = min(set(transformed) ^ set(image) | {c for c in set(transformed) | set(image) if transformed.get(c, 0) != image.get(c, 0)})
            return orbit[0], difference
    raise AssertionError("expected a rho-breaking cut witness")


def negative_control(L: int, basis: list[dict[int, int]], members: list[tuple[int, ...]], orbit_of: list[int]) -> None:
    mutated = [dict(vector) for vector in basis]
    first = next(iter(mutated[0]))
    mutated[0][first] *= -1
    span = QSpan()
    for vector in mutated:
        span.add(vector)
    B4 = four_B_rows(L, orbit_of, members, None)
    fails = any(span.residual(multiply(vector, B4, None)) for vector in mutated)
    check("negative_control_mutated_certificate_fails", fails, "one sign flip destroys 4B invariance")


def main() -> int:
    if not ARTIFACT.exists():
        raise FileNotFoundError(f"missing producer artifact: {ARTIFACT}")
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact_schema", artifact.get("provenance", {}).get("artifact_schema") == "provenance/data/checks-v1")
    stored_hash = artifact.get("provenance", {}).get("artifact_sha256")
    canonical = json.dumps({"provenance": {k: v for k, v in artifact["provenance"].items() if k != "artifact_sha256"},
                            "data": artifact["data"], "checks": artifact["checks"]},
                           sort_keys=True, separators=(",", ":")).encode()
    check("artifact_hash", stored_hash == hashlib.sha256(canonical).hexdigest())

    factor_rows = artifact["data"]["tensor_factorisation"]
    for row in factor_rows:
        extension_record(row["L"], row)

    certificate_rows = artifact["data"]["unreachable_invariant_blocks"]
    W4_basis = None
    W4_members = None
    W4_orbit_of = None
    for row in certificate_rows:
        basis, members, orbit_of = certificate_record(row)
        if row["L"] == 4:
            W4_basis, W4_members, W4_orbit_of = basis, members, orbit_of
    assert W4_basis is not None and W4_members is not None and W4_orbit_of is not None
    negative_control(4, W4_basis, W4_members, W4_orbit_of)

    local_action_checks()
    witness = artifact["data"]["cut_locality_witness"]
    observed = cut_witness(witness["L_old"])
    check("new_edge_piece_breaks_rho", observed == (witness["orbit_representative"], witness["difference_configuration"]), str(observed))

    if FAILURES:
        print(f"FAILED {len(FAILURES)}/{CHECKS}: {', '.join(FAILURES)}")
        return 1
    print(f"OK {CHECKS} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
