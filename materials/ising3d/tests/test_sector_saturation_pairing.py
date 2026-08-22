#!/usr/bin/env python3
"""Clean-room verifier for experiments/e132_sector_saturation_pairing.py.

Independently rebuilds the <tau,rho> orbit space of the 2xL ladder, the
stabiliser-module columns under {A, 4B} on the vacuum, and the pairing
kernel certificate, and re-checks the stored annihilator bases over Q exactly.
No helper or constant is imported from the producer.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import defaultdict, Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "sector_saturation_pairing.json"
P1 = 2_147_483_647
P2 = 2_147_483_629

EXPECTED_CYCLIC = {3: (14, 14), 4: (42, 44), 5: (142, 152), 6: (494, 560)}

CHECKS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str) -> None:
    CHECKS.append((name, bool(ok), detail))


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
        answer |= ((config >> (2 * rung)) & 3) << (2 * (L - 1 - rung))
    return answer


def group_images(config: int, L: int) -> tuple[int, ...]:
    return tuple(sorted({config, tau(config, L), rho(config, L), tau(rho(config, L), L)}))


def orbits(L: int) -> tuple[list[int], list[tuple[int, ...]]]:
    orbit_of = [-1] * (1 << (2 * L))
    members: list[tuple[int, ...]] = []
    for config in range(1 << (2 * L)):
        if orbit_of[config] >= 0 or config.bit_count() & 1:
            continue
        orbit = group_images(config, L)
        index = len(members)
        for image in orbit:
            orbit_of[image] = index
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
    answer: list[dict[int, int]] = []
    for orbit in members:
        counts: dict[int, int] = defaultdict(int)
        for config in orbit:
            for u, v in edges(L):
                counts[orbit_of[config ^ ((1 << u) | (1 << v))]] += 1
        row: dict[int, int] = {}
        for target, count in counts.items():
            size = len(members[target])
            if count % size:
                raise AssertionError("4B orbit coefficient was not integral")
            coefficient = 4 * (count // size)
            if coefficient:
                row[target] = coefficient % modulus if modulus else coefficient
        answer.append(row)
    return answer


def multiply(vector: dict[int, int], rows: list[dict[int, int]], modulus: int | None) -> dict[int, int]:
    answer: dict[int, int] = defaultdict(int)
    for source, coefficient in vector.items():
        for target, entry in rows[source].items():
            answer[target] += coefficient * entry
    if modulus is not None:
        return {index: value % modulus for index, value in answer.items() if value % modulus}
    return {index: value for index, value in answer.items() if value}


class ModEchelon:
    def __init__(self, P: int) -> None:
        self.P = P
        self.rows: dict[int, dict[int, int]] = {}

    def add(self, vector: dict[int, int]) -> int | None:
        P = self.P
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
                reduced = (row.get(index, 0) - scale * value) % P
                if reduced:
                    row[index] = reduced
                elif index in row:
                    del row[index]
        return None

    @property
    def rank(self) -> int:
        return len(self.rows)


class QSpan:
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
                reduced = row.get(index, Fraction(0)) - scale * value
                if reduced:
                    row[index] = reduced
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


def closure_rank(L: int, P: int) -> int:
    orbit_of, members = orbits(L)
    B4 = four_B_rows(L, orbit_of, members, P)
    A = [(2 * L - 2 * orbit[0].bit_count()) % P for orbit in members]
    echelon = ModEchelon(P)
    if echelon.add({orbit_of[0]: 1}) is None:
        raise AssertionError("vacuum seeding failed")
    todo = [{orbit_of[0]: 1}]
    while todo:
        vector = todo.pop()
        A_image = {index: (A[index] * coefficient) % P
                   for index, coefficient in vector.items() if (A[index] * coefficient) % P}
        B_image = multiply(vector, B4, P)
        for image in (A_image, B_image):
            if image and echelon.add(image) is not None:
                todo.append(image)
    return echelon.rank


def verify_record(record: dict) -> None:
    L = record["L"]
    orbit_of, members = orbits(L)
    check(f"burnside_K_L{L}", len(members) == burnside_K_dim(L), f"{len(members)}")
    rank1 = closure_rank(L, P1)
    rank2 = closure_rank(L, P2)
    expected_cyclic, expected_K = EXPECTED_CYCLIC[L]
    check(f"two_prime_rank_L{L}", rank1 == rank2 == expected_cyclic,
          f"{rank1}/{rank2} == cyclic dim {expected_cyclic}")
    by_rep = {orbit[0]: index for index, orbit in enumerate(members)}
    basis: list[dict[int, int]] = []
    for stored in record["pairing_kernel_basis"]:
        vector: dict[int, int] = {}
        for term in stored:
            index = by_rep.get(term["representative"])
            if index is not None:
                vector[index] = vector.get(index, 0) + term["coefficient"]
        basis.append({i: coefficient for i, coefficient in vector.items() if coefficient})
    span = QSpan()
    for vector in basis:
        span.add(vector)
    check(f"kernel_dim_L{L}", len(basis) == expected_K - expected_cyclic
          and span.rank == len(basis), f"{len(basis)} / K - cyclic = {expected_K - expected_cyclic}")
    B4 = four_B_rows(L, orbit_of, members, None)
    A = [2 * L - 2 * orbit[0].bit_count() for orbit in members]
    b_invariant = all(not span.residual(multiply(vector, B4, None)) for vector in basis)
    a_invariant = all(
        not span.residual({index: A[index] * coefficient for index, coefficient in vector.items() if A[index] * coefficient})
        for vector in basis)
    homogeneous = all(len({members[index][0].bit_count() for index in vector}) == 1 for vector in basis)
    vacuum_orth = all(orbit_of[0] not in vector for vector in basis)
    sectors = Counter(next(iter({members[index][0].bit_count() for index in vector})) for vector in basis if vector)
    check(f"kernel_exact_L{L}", b_invariant and a_invariant and homogeneous and vacuum_orth,
          f"sectors {dict(sorted(sectors.items()))}")


def negative_controls() -> None:
    orbit_of, members = orbits(4)
    record = json.loads(ARTIFACT.read_text())["data"]["pairing_records"][1]
    by_rep = {orbit[0]: index for index, orbit in enumerate(members)}
    basis = []
    for stored in record["pairing_kernel_basis"]:
        vector: dict[int, int] = {}
        for term in stored:
            index = by_rep[term["representative"]]
            vector[index] = vector.get(index, 0) + term["coefficient"]
        basis.append({i: coefficient for i, coefficient in vector.items() if coefficient})
    mutated = [dict(vector) for vector in basis]
    first = next(iter(mutated[0]))
    mutated[0][first] *= -1
    span = QSpan()
    for vector in mutated:
        span.add(vector)
    B4 = four_B_rows(4, orbit_of, members, None)
    broken = any(span.residual(multiply(vector, B4, None)) for vector in mutated)
    check("negative_control_sign_flip_breaks_invariance_4", broken, "one sign flip destroys B-invariance")


def main() -> int:
    started = time.process_time()
    if not ARTIFACT.exists():
        raise FileNotFoundError(f"missing producer artifact: {ARTIFACT}")
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact_schema", artifact.get("provenance", {}).get("artifact_schema") == "provenance/data/checks-v1",
          artifact.get("provenance", {}).get("artifact_schema", "missing"))
    stored_hash = artifact.get("provenance", {}).get("artifact_sha256")
    canonical = json.dumps(
        {"provenance": {k: v for k, v in artifact["provenance"].items() if k != "artifact_sha256"},
         "data": artifact["data"], "checks": artifact["checks"]},
        sort_keys=True, separators=(",", ":")).encode()
    check("artifact_hash", stored_hash == hashlib.sha256(canonical).hexdigest(),
          hashlib.sha256(canonical).hexdigest()[:16])
    for record in artifact["data"]["pairing_records"]:
        verify_record(record)
    negative_controls()
    failed = [name for name, ok, _ in CHECKS if not ok]
    for name, ok, detail in CHECKS:
        print(f"  {name}: {'ok' if ok else 'FAIL'} {detail}")
    if failed:
        print(f"FAILED {len(failed)}/{len(CHECKS)}: {failed}")
        return 1
    print(f"OK {len(CHECKS)} checks (cpu {round(time.process_time() - started, 3)} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
