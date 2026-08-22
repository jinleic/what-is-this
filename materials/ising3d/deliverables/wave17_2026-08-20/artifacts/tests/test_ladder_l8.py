#!/usr/bin/env python3
"""Clean-room verifier for the L = 8 ladder sector-saturation certificate.

Artifact: results/ladder/l8_saturation.json (producer
experiments/e142_sector_saturation_l8.py).

This verifier shares no code with the producer.  It rebuilds the ladder
orbit space, the exact integer `4B` matrix, and the closure/echelon machinery
from scratch, then:

  (V1) recomputes `dim K_L` independently for L = 7, 8 (Burnside closed form,
       direct orbit enumeration of the rank-two reflection group <tau, rho>,
       and the per-sector Burnside sums from directly enumerated fixed sets);
  (V2) recomputes the full modular closures for the fast regression cases
       L = 3..6 over F_p at p = 2147483647 AND p = 2147483629, together with
       the exact-Q pairing-kernel certificates, and checks the wave-14 values
       cyclic 14/42/142/494 and W = 0/2/10/66 with their sector splits;
  (V3) re-validates the STORED L = 8 annihilator basis over Q with exact
       Fractions (4B-invariance, A-invariance, vacuum-orthogonality, sector
       homogeneity, independence) and checks the arithmetic identities
       W_8 = K_8 - cyclic, agreement of every recorded prime rank with the
       recorded exact cyclic dimension, and the sector symmetry k <-> 16-k;
  (V4) checks the artifact's L = 7 record against the frozen wave-14 sibling
       artifact results/algebra_growth/sector_saturation_tensor.json;
  (V5) negative control: a single sign flip in one stored L = 4 kernel vector
       must destroy 4B-invariance.

The full L = 8 closure (rank 6562, ~7 CPU-minutes at the small primes) is NOT
rerun here; its lower-bound half is anchored by the recorded multi-prime
agreement checked in (V3) and by the L <= 6 regressions in (V2), and its
upper-bound half IS independently recomputed in (V3).
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "ladder" / "l8_saturation.json"
SIBLING = ROOT / "results" / "algebra_growth" / "sector_saturation_tensor.json"
P1 = 2_147_483_647
P2 = 2_147_483_629

CHECKS = 0
FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    status = "ok" if condition else "FAIL"
    print(f"[{status}] {name}" + (f" -- {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Independent ladder, symmetry, and orbit construction.
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
        answer |= ((config >> (2 * rung)) & 3) << (2 * (L - 1 - rung))
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
        for image in orbit:
            orbit_of[image] = index
        members.append(orbit)
    return orbit_of, members


def edges(L: int) -> list[tuple[int, int]]:
    answer = [(2 * rung, 2 * rung + 1) for rung in range(L)]
    for rung in range(L - 1):
        answer.extend(((2 * rung, 2 * rung + 2), (2 * rung + 1, 2 * rung + 3)))
    return answer


def burnside(L: int) -> int:
    return (2 ** (2 * L - 1) + 3 * 2**L) // 4


def four_B_rows(L: int, orbit_of: list[int], members: list[tuple[int, ...]], modulus: int | None) -> list[dict[int, int]]:
    answer: list[dict[int, int]] = []
    table = edges(L)
    for orbit in members:
        counts: dict[int, int] = defaultdict(int)
        for config in orbit:
            for u, v in table:
                counts[orbit_of[config ^ ((1 << u) | (1 << v))]] += 1
        row: dict[int, int] = {}
        for target, count in counts.items():
            size = len(members[target])
            if count % size:
                raise AssertionError("orbit integrality failed")
            coefficient = 4 * (count // size)
            if coefficient:
                row[target] = coefficient % modulus if modulus else coefficient
        answer.append(row)
    return answer


def multiply(vector: dict[int, int], rows: list[dict[int, int]], modulus: int | None) -> dict[int, int]:
    answer: dict[int, int] = defaultdict(int)
    for source, coefficient in vector.items():
        for target, matrix_entry in rows[source].items():
            answer[target] += coefficient * matrix_entry
    if modulus is not None:
        return {index: value % modulus for index, value in answer.items() if value % modulus}
    return {index: value for index, value in answer.items() if value}


class ModEchelon:
    """Sparse row echelon over F_P, pivot = maximum coordinate."""

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


def closure_rank(L: int, P: int) -> ModEchelon:
    orbit_of, members = even_orbits(L)
    B4 = four_B_rows(L, orbit_of, members, P)
    A_diagonal = [(2 * L - 2 * orbit[0].bit_count()) % P for orbit in members]
    echelon = ModEchelon(P)
    seed = echelon.add({orbit_of[0]: 1})
    if seed is None:
        raise AssertionError("vacuum failed to seed")
    todo = [seed]
    processed: set[int] = set()
    while todo:
        pivot = todo.pop()
        if pivot in processed:
            continue
        processed.add(pivot)
        vector = echelon.rows[pivot]
        for image in (
            {index: A_diagonal[index] * coefficient % P for index, coefficient in vector.items()},
            multiply(vector, B4, P),
        ):
            new_pivot = echelon.add(image)
            if new_pivot is not None:
                todo.append(new_pivot)
    return echelon


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

    def add(self, vector) -> bool:
        row = self.residual(vector)
        if not row:
            return False
        pivot = max(row)
        value = row[pivot]
        self.rows[pivot] = {index: v / value for index, v in row.items()}
        return True

    @property
    def rank(self) -> int:
        return len(self.rows)


def kernel_basis(L: int, echelon: ModEchelon, members: list[tuple[int, ...]]) -> list[dict[int, int]]:
    sizes = [len(orbit) for orbit in members]
    P = echelon.P
    constraints = {
        pivot: {index: coefficient * sizes[index] % P for index, coefficient in row.items()}
        for pivot, row in echelon.rows.items()
    }
    candidates = []
    free_indices = [index for index in range(len(members)) if index not in constraints]
    for free in free_indices:
        vector: dict[int, int] = {free: 1}
        for pivot in sorted(constraints):
            row = constraints[pivot]
            residual = sum(row.get(index, 0) * coefficient for index, coefficient in vector.items()) % P
            if residual:
                vector[pivot] = -residual * pow(row[pivot], P - 2, P) % P
        lifted = {}
        for index, residue in vector.items():
            integer = residue if residue <= P // 2 else residue - P
            if integer:
                lifted[index] = integer
        candidates.append(lifted)
    return candidates


def validate(L: int, basis: list[dict[int, int]], orbit_of: list[int], members: list[tuple[int, ...]]) -> dict:
    B4 = four_B_rows(L, orbit_of, members, None)
    A_diagonal = [2 * L - 2 * orbit[0].bit_count() for orbit in members]
    span = QSpan()
    for vector in basis:
        span.add(vector)
    sector_lists = [sorted({members[index][0].bit_count() for index in vector}) for vector in basis]
    b_ok = all(not span.residual(multiply(vector, B4, None)) for vector in basis)
    a_ok = all(not span.residual({i: A_diagonal[i] * c for i, c in vector.items() if A_diagonal[i] * c})
               for vector in basis)
    return {
        "span_rank": span.rank,
        "independent": span.rank == len(basis),
        "homogeneous": all(len(values) == 1 for values in sector_lists),
        "vacuum_orthogonal": all(orbit_of[0] not in vector for vector in basis),
        "four_B_invariant": b_ok,
        "A_invariant": a_ok,
        "sectors": Counter(values[0] for values in sector_lists if len(values) == 1),
    }


def deserialise(basis: list[list[dict[str, int]]], orbit_of: list[int]) -> list[dict[int, int]]:
    answer = []
    for vector in basis:
        entry = {orbit_of[term["representative"]]: term["coefficient"] for term in vector}
        answer.append(entry)
    return answer


# ---------------------------------------------------------------------------
# Verification stages.
# ---------------------------------------------------------------------------
def stage_K_dimensions() -> None:
    for L, expect in ((7, 2144), (8, 8384)):
        orbit_of, members = even_orbits(L)
        check(f"K{L}_burnside_formula", burnside(L) == expect, f"2^(2L-1)+3*2^L)/4 = {burnside(L)}")
        check(f"K{L}_direct_orbits", len(members) == expect, f"{len(members)} even <tau,rho>-orbits")
        fixed = [Counter() for _ in range(3)]
        for config in range(1 << (2 * L)):
            if config.bit_count() & 1:
                continue
            weight = config.bit_count()
            if tau(config, L) == config:
                fixed[0][weight] += 1
            if rho(config, L) == config:
                fixed[1][weight] += 1
            if tau(rho(config, L), L) == config:
                fixed[2][weight] += 1
        sector_total = 0
        for k in range(0, 2 * L + 1, 2):
            plain = sum(1 for config in range(1 << (2 * L)) if config.bit_count() == k)
            sector_total += (plain + fixed[0][k] + fixed[1][k] + fixed[2][k]) // 4
        check(f"K{L}_sector_burnside_sum", sector_total == expect, f"per-sector Burnside total {sector_total}")


def stage_negative_control(artifact: dict) -> None:
    record = [r for r in artifact["data"]["regression_records"] if r["L"] == 4][0]
    orbit_of, members = even_orbits(4)
    basis = deserialise(record["basis"], orbit_of)
    check("negative_control_seed", validate(4, basis, orbit_of, members)["four_B_invariant"],
          "stored L=4 basis is 4B-invariant")
    # Flip ONE coefficient: the span is closed under global rescaling, so only a
    # single-entry sign change can destroy invariance.
    victim = min(basis[0])
    basis[0] = dict(basis[0])
    basis[0][victim] = -basis[0][victim]
    check("negative_control_sign_flip", not validate(4, basis, orbit_of, members)["four_B_invariant"],
          "flipping one coefficient of one kernel vector destroys 4B-invariance")


def stage_regressions() -> None:
    expected = {3: (14, 0, {}), 4: (42, 2, {4: 2}), 5: (142, 10, {4: 5, 6: 5}), 6: (494, 66, {4: 15, 6: 36, 8: 15})}
    for L, (cyclic, W, sectors) in expected.items():
        orbit_of, members = even_orbits(L)
        ech1 = closure_rank(L, P1)
        ech2 = closure_rank(L, P2)
        check(f"cyclic_L{L}_two_primes", ech1.rank == ech2.rank == cyclic,
              f"rank p1 {ech1.rank}, p2 {ech2.rank}, expected {cyclic}")
        basis = kernel_basis(L, ech1, members)
        report = validate(L, basis, orbit_of, members)
        check(f"annihilator_L{L}_over_Q",
              report["span_rank"] == W and report["independent"] and report["homogeneous"]
              and report["vacuum_orthogonal"] and report["four_B_invariant"] and report["A_invariant"]
              and dict(sorted(report["sectors"].items())) == sectors,
              f"W dim {report['span_rank']}, sectors {dict(sorted(report['sectors'].items()))}")
        check(f"sandwich_L{L}", cyclic + W == len(members) == burnside(L),
              f"{cyclic} + {W} = {len(members)}")


def stage_L8(artifact: dict) -> None:
    record = artifact["data"]["L8_record"]
    orbit_of, members = even_orbits(8)
    K = len(members)
    cyclic = record["cyclic_Q_dim"]
    W = record["W_dim"]
    check("L8_K_dim", record["K_dim"] == K == 8384, f"K_8 = {K}")
    check("L8_rank_prime_agreement",
          all(r == cyclic for r in record["mod_p_ranks"]) and all(r == cyclic for r in record.get("primes_all_agree", [])),
          f"ranks {record.get('primes_all_agree', record['mod_p_ranks'])} vs cyclic {cyclic}")
    basis = deserialise(record["basis"], orbit_of)
    check("L8_basis_size", len(basis) == W, f"{len(basis)} stored kernel vectors vs W_8 = {W}")
    started = time.process_time()
    report = validate(8, basis, orbit_of, members)
    check("L8_annihilator_over_Q",
          report["span_rank"] == W and report["independent"] and report["homogeneous"]
          and report["vacuum_orthogonal"] and report["four_B_invariant"] and report["A_invariant"],
          f"re-validated over Q in {time.process_time() - started:.1f}s cpu")
    sectors = {int(k): v for k, v in record["W_sectors"].items()}
    check("L8_sector_split_matches_basis", dict(sorted(report["sectors"].items())) == sectors,
          f"{dict(sorted(report['sectors'].items()))}")
    check("L8_sector_symmetry", all(sectors.get(16 - k) == v for k, v in sectors.items()),
          f"split {sectors} symmetric under k -> 16-k")
    check("L8_sandwich_identity", cyclic + W == K, f"{cyclic} + {W} = {K}")
    pattern = artifact["data"]["pattern"]
    check("W_pattern_sequence", pattern["W_sequence"] == [0, 2, 10, 66, 364, W],
          f"W sequence {pattern['W_sequence']}")
    check("L8_growth_verdict", (W > 364) == pattern["W_strictly_increasing"],
          f"W_8 = {W}, strictly increasing = {pattern['W_strictly_increasing']}")
    check("L8_evaluation_bound", cyclic >= 2 ** 8, f"dim g_8 >= {cyclic} >= 256")


def stage_sibling_consistency(artifact: dict) -> None:
    if not SIBLING.exists():
        check("sibling_artifact_present", False, str(SIBLING))
        return
    sibling = json.loads(SIBLING.read_text())
    rows = {row["L"]: row for row in sibling["data"]["unreachable_invariant_blocks"]}
    l7 = artifact["data"]["L7_record"]
    row7 = rows[7]
    check("L7_matches_wave14",
          (l7["cyclic_Q_dim"], l7["K_dim"], l7["W_dim"]) ==
          (row7["cyclic_Q_dim"], row7["K_dim"], row7["unreachable_block_dim"]) ==
          (1780, 2144, 364),
          f"producer L=7 {l7['cyclic_Q_dim']}/{l7['K_dim']} W={l7['W_dim']} vs frozen tensor front")



def main() -> int:
    if not ARTIFACT.exists():
        print(f"missing artifact: {ARTIFACT}", file=sys.stderr)
        return 1
    artifact = json.loads(ARTIFACT.read_text())
    failed_checks = [c for c in artifact["checks"] if not c["passed"]]
    check("artifact_all_checks_passed", not failed_checks, f"{len(artifact['checks'])} checks, {len(failed_checks)} failing")
    stage_K_dimensions()
    stage_regressions()
    stage_L8(artifact)
    stage_sibling_consistency(artifact)
    stage_negative_control(artifact)
    print(f"\n{CHECKS - len(FAILURES)}/{CHECKS} checks passed")
    if FAILURES:
        print("FAILED: " + ", ".join(FAILURES))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
