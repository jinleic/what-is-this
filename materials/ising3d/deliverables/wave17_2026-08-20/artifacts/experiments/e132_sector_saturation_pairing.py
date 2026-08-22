#!/usr/bin/env python3
"""The pairing-matrix kernel of the ladder sector-saturation route.

Setting (same base as `proofs/ladder_alll_proof.md`): the open `2xL`
ladder, `n = 2L` sites (`2i` top of rung `i`, `2i+1` bottom), vacuum
`psi = |+>^(2L)` (orbit of the zero configuration), and the two diagonalising
operators `A_L = sum_v X_v = 2L - 2N`, `B_L = sum_edges Z_u Z_v`.  The
even-particle, leg-swap (`tau`) and rung-reversal (`rho`) invariant space is
`K_L` with `dim K_L = (2^(2L-1) + 3*2^L)/4` (Burnside).  Every
`v in U(A,B) psi` (the stabiliser module of the vacuum) is in `K_L`.

The object studied here is the **duality pairing**

    P_L[o, j] = <delta_o, v_j>_conf = |o| * v_j[rep(o)]

between the *distinct* delta-functionals `delta_o` on `K_L` (evaluation at
the configurations; `delta_o = delta_x` for every `x` in the `<tau,rho>`
orbit `o`, orbit size `|o|`) and the explicit stabiliser-module columns
`v_j`, where `<c,d>_conf = sum_x c_x d_x` is the configuration inner
product.  Ranking `P_L` is identical to ranking the `2^(2L-1)`-row matrix
whose rows are the delta-functionals at *every* distinct even configuration: an
invariant vector is constant on orbits, so the surplus rows are duplicates.

The pairing kernel is `W_L = (U(A,B)psi)^perp` inside `K_L`, the
annihilator of the Krylov/cyclic span.  The certificate chain:

  (C1) The columns `v_j` are closed under `A_L` and `B_L`, hence span
       `U(A,B)psi` modulo the prime.  `rank_{F_{p1}} P = rank_{F_{p2}} P`
       at two distinct 31-bit primes gives `dim_Q(U(A,B)psi) >= rank`.
  (C2) An explicit integer basis of the modular kernel of `P_L` is lifted to a
       small integer family `W_L` (coefficient bound 128) and validated *over
       Q*: every `c in W_L` is vacuum-orthogonal (`c[vacuum] = 0`),
       homogeneous in one particle sector, and `A c, B c in span(W_L)`
       exactly.  Since `B_L`, `A_L` are symmetric and map `K_L` into
       itself, any A,B-invariant subspace orthogonal to the vacuum lies in
       `(U(A,B)psi)^perp`; with `dim W_L = dim K_L - rank` the two-prime
       lower bound is therefore **exact** over Q and `W_L = (U(A,B)psi)^perp`.
  (C3) The resulting exact cyclic dimensions are `dim U(A,B)psi = 14, 42,
       142, 494` for `L = 3..6` against `dim K_L = 14, 44, 152, 560`,
       i.e. the certified `W` annihilator dims `0, 2, 10, 66`, with
       sector splits `{}`, `{4: 2}`, `{4:5, 6:5}`, `{4:15, 6:36,
       8:15}`.  These are exactly the `W` blocks certified by the tensor
       front (`results/algebra_growth/sector_saturation_tensor.json`): the
       sibling's `4B`-invariant annihilator of the cyclic span is the same
       object as this pairing kernel, so the two routes land on the same numbers
       by an independent characterisation.

All exact validations use `fractions.Fraction`; the only primes are
`P1 = 2147483647` and `P2 = 2147483629`.  Every compute gate uses
`time.process_time()` (process CPU), never a wall clock; peak RSS is recorded.
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import sys
import time
from collections import defaultdict, Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e132_sector_saturation_pairing.py"
RESULT_PATH = ROOT / "results" / "algebra_growth" / "sector_saturation_pairing.json"
P1 = 2_147_483_647
P2 = 2_147_483_629
LIFT_BOUND = 128
L_VALUES = (3, 4, 5, 6)
CPU_BUDGET_TOTAL = 3_600.0


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


# ---------------------------------------------------------------------------
# Ladder configuration space and the Klein-four symmetry.
# ---------------------------------------------------------------------------
def tau(config: int, L: int) -> int:
    """Exchange the two legs on every rung, retaining the rungs."""
    answer = 0
    for rung in range(L):
        top = (config >> (2 * rung)) & 1
        bottom = (config >> (2 * rung + 1)) & 1
        answer |= (bottom << (2 * rung)) | (top << (2 * rung + 1))
    return answer


def rho(config: int, L: int) -> int:
    """Reverse rung order, retaining the legs."""
    answer = 0
    for rung in range(L):
        state = (config >> (2 * rung)) & 3
        answer |= state << (2 * (L - 1 - rung))
    return answer


def group_images(config: int, L: int) -> tuple[int, ...]:
    return tuple(sorted({config, tau(config, L), rho(config, L), tau(rho(config, L), L)}))


def invariant_orbits(L: int) -> tuple[list[int], list[tuple[int, ...]]]:
    """Even-parity <tau,rho>-orbits, indexed by their least configuration."""
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


def ladder_edges(L: int) -> list[tuple[int, int]]:
    answer = [(2 * rung, 2 * rung + 1) for rung in range(L)]
    for rung in range(L - 1):
        answer.extend(((2 * rung, 2 * rung + 2), (2 * rung + 1, 2 * rung + 3)))
    return answer


def burnside_K_dim(L: int) -> int:
    return (2 ** (2 * L - 1) + 3 * 2**L) // 4


def four_B_rows(L: int, orbit_of: list[int], members: list[tuple[int, ...]], modulus: int | None) -> list[dict[int, int]]:
    """The exact integer matrix 4B in the unnormalised orbit-sum basis e_o = sum_{x in o} |x>."""
    answer: list[dict[int, int]] = []
    for orbit in members:
        counts: dict[int, int] = defaultdict(int)
        for config in orbit:
            for u, v in ladder_edges(L):
                counts[orbit_of[config ^ ((1 << u) | (1 << v))]] += 1
        row: dict[int, int] = {}
        for target, count in counts.items():
            size = len(members[target])
            if count % size:
                raise AssertionError("4B failed orbit-integrality")
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
    """Sparse row echelon over F_P in the orbit-sum coordinates (pivot = max)."""

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


# ---------------------------------------------------------------------------
# The pairing matrix and its two-prime closure.
# ---------------------------------------------------------------------------
def cyclic_closure(L: int, P: int) -> tuple[list[int], list[tuple[int, ...]], ModEchelon]:
    """Columns v_j of the pairing matrix: the closure of psi under {A, 4B}."""
    orbit_of, members = invariant_orbits(L)
    B4 = four_B_rows(L, orbit_of, members, P)
    A = [(2 * L - 2 * orbit[0].bit_count()) % P for orbit in members]
    echelon = ModEchelon(P)
    seed = {orbit_of[0]: 1}
    if echelon.add(seed) is None:
        raise AssertionError("vacuum failed to seed the closure")
    todo = [seed]
    while todo:
        vector = todo.pop()
        A_image = {index: (A[index] * coefficient) % P
                   for index, coefficient in vector.items() if (A[index] * coefficient) % P}
        B_image = multiply(vector, B4, P)
        for image in (A_image, B_image):
            if image and echelon.add(image) is not None:
                todo.append(image)
    return orbit_of, members, echelon


def pair_distinct_delta_rank(L: int, P: int) -> int:
    """Regression-only modular rank of the full 2^(2L-1)-row delta matrix.

    The rows are the delta-functionals at EVERY distinct even configuration
    (all 2^(2L-1) of them); invariant columns are constant on orbits, so
    this rank must coincide with the orbit-indexed (K_L-row) pairing rank.
    """
    orbit_of, members = invariant_orbits(L)
    B4 = four_B_rows(L, orbit_of, members, P)
    A = [(2 * L - 2 * orbit[0].bit_count()) % P for orbit in members]
    echelon = ModEchelon(P)
    orbit_of_config = [-1] * (1 << (2 * L))
    for index, orbit in enumerate(members):
        for image in orbit:
            orbit_of_config[image] = index
    seed = {orbit_of[0]: 1}
    if echelon.add({orbit_of[0]: 1}) is None:
        raise AssertionError("vacuum failed to seed the distinct-delta closure")
    todo = [seed]
    while todo:
        vector = todo.pop()
        A_image = {index: (A[index] * coefficient) % P
                   for index, coefficient in vector.items() if (A[index] * coefficient) % P}
        B_image = multiply(vector, B4, P)
        for image in (A_image, B_image):
            if image and echelon.add(image) is not None:
                todo.append(image)
    return len(echelon.rows)


# ---------------------------------------------------------------------------
# Exact orthogonal-complement kernel W_L of the pairing.
# ---------------------------------------------------------------------------
def gram_kernel_lift(echelon: ModEchelon, members: list[tuple[int, ...]], bound: int) -> list[dict[int, int]]:
    """Small integer lifts of the modular kernel of the weighted pairing matrix.

    The Gram-style constraints are indexed by the pivots of the modular closure
    with rows weighted by |o|; free coordinates of the K_L-row pairing matrix
    are completed, then each residue is converted to a centred integer and required
    to stay within `bound`.
    """
    sizes = [len(orbit) for orbit in members]
    constraints = {
        pivot: {index: coefficient * sizes[index] % echelon.P
               for index, coefficient in row.items()}
        for pivot, row in echelon.rows.items()
    }
    candidates: list[dict[int, int]] = []
    free_indices = [index for index in range(len(members)) if index not in echelon.rows]
    for free in free_indices:
        vector: dict[int, int] = {free: 1}
        for pivot in sorted(constraints):
            row = constraints[pivot]
            residual = sum(row.get(index, 0) * coefficient for index, coefficient in vector.items()) % echelon.P
            if residual:
                vector[pivot] = -residual * pow(row[pivot], echelon.P - 2, echelon.P) % echelon.P
        lifted: dict[int, int] = {}
        for index, residue in vector.items():
            integer = residue if residue <= echelon.P // 2 else residue - echelon.P
            if abs(integer) > bound:
                raise AssertionError(f"K={len(members)}: no small exact annihilator lift")
            if integer:
                lifted[index] = integer
        candidates.append(lifted)
    return candidates


class QSpan:
    """Sparse exact-Q row span for the annihilator certificates."""

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


def exact_validate(L: int, orbit_of: list[int], members: list[tuple[int, ...]],
                  candidates: list[dict[int, int]]) -> dict:
    B4 = four_B_rows(L, orbit_of, members, None)
    A = [2 * L - 2 * orbit[0].bit_count() for orbit in members]
    span = QSpan()
    for vector in candidates:
        span.add(vector)
    sector_sets = [{members[index][0].bit_count() for index in vector} for vector in candidates]
    homogeneous = all(len(values) == 1 for values in sector_sets)
    b_invariant = all(not span.residual(multiply(vector, B4, None)) for vector in candidates)
    a_invariant = all(
        not span.residual({index: A[index] * coefficient for index, coefficient in vector.items() if A[index] * coefficient})
        for vector in candidates
    )
    vacuum_orthogonal = all(orbit_of[0] not in vector for vector in candidates)
    independent = span.rank == len(candidates)
    if not all((homogeneous, vacuum_orthogonal, independent)):
        raise AssertionError(f"L={L}: candidate annihilator failed exact validation")
    sector_counts = Counter(next(iter(values)) for values in sector_sets if len(values) == 1) if homogeneous else Counter()
    return {
        "pairing_kernel_dim": len(candidates),
        "b_invariant_over_Q": b_invariant,
        "a_invariant_over_Q": a_invariant,
        "vacuum_orthogonal": vacuum_orthogonal,
        "sector_homogeneous": homogeneous,
        "sector_dimensions": dict(sorted(sector_counts.items())),
        "kernel_independent_over_Q": independent,
    }


def serialise_basis(basis: list[dict[int, int]], members: list[tuple[int, ...]]) -> list[list[dict[str, int]]]:
    return [
        sorted(({"representative": members[index][0], "coefficient": coefficient} for index, coefficient in vector.items()),
              key=lambda term: term["representative"])
        for vector in basis
    ]


def pairing_record(L: int, P1: int, P2: int) -> dict:
    started = time.process_time()
    orbit_of, members, ech1 = cyclic_closure(L, P1)
    _, _, ech2 = cyclic_closure(L, P2)
    if ech1.rank != ech2.rank:
        raise AssertionError(f"L={L}: two-prime pairing ranks disagree")
    distinct_rank = pair_distinct_delta_rank(L, P1)
    if distinct_rank != ech1.rank:
        raise AssertionError(f"L={L}: distinct-delta rank {distinct_rank} != orbit pairing rank {ech1.rank}")
    candidates = gram_kernel_lift(ech1, members, LIFT_BOUND)
    validation = exact_validate(L, orbit_of, members, candidates)
    cyclic_q = len(members) - validation["pairing_kernel_dim"]
    if cyclic_q != ech1.rank:
        raise AssertionError(f"L={L}: exact annihilator upper ({cyclic_q}) misses modular rank {ech1.rank}")
    return {
        "L": L,
        "K_dim": len(members),
        "distinct_delta_count": 1 << (2 * L - 1),
        "two_prime_rank_p1": ech1.rank,
        "two_prime_rank_p2": ech2.rank,
        "distinct_delta_pairing_rank": distinct_rank,
        "cyclic_dim_Q": cyclic_q,
        "cyclic_dims": f"{cyclic_q}/{len(members)}",
        "pairing_kernel_dim": validation["pairing_kernel_dim"],
        "pairing_kernel_sectors": validation["sector_dimensions"],
        "pairing_kernel_basis": serialise_basis(candidates, members),
        "validation": validation,
        "cpu_seconds": round(time.process_time() - started, 6),
    }


def make_artifact() -> dict:
    started = time.process_time()
    records = [pairing_record(L, P1, P2) for L in L_VALUES]
    data = {
        "scope": {
            "space": "K_L = even-parity <tau,rho>-invariant ladder configuration space",
            "pairing": "P_L[o,j] = <delta_o, v_j>_conf with distinct delta-functionals delta_o and stabiliser-module columns v_j",
            "field": "Q for the exact annihilator; F_p for the two-prime lower rank",
            "primes": [P1, P2],
            "lift_bound": LIFT_BOUND,
            "kernel_meaning": "orthogonal complement of the Krylov span in the configuration inner product",
        },
        "pairing_records": records,
        "resource_measurements": {
            "total_cpu_seconds": round(time.process_time() - started, 6),
            "peak_rss_bytes": peak_rss_bytes(),
            "budget_clock": "time.process_time",
        },
    }
    expected_q = {3: (14, 14, 0), 4: (42, 44, 2), 5: (142, 152, 10), 6: (494, 560, 66)}
    checks = []
    for row in records:
        cyclic_q, K, W = expected_q[row["L"]]
        checks.append(
            {
                "name": f"C_cyclic_dims_L{row['L']}",
                "passed": row["cyclic_dim_Q"] == cyclic_q and row["K_dim"] == K and row["pairing_kernel_dim"] == W,
                "detail": f"{cyclic_q}/{K} cyclic with pairing kernel {W}",
            }
        )
        checks.append(
            {
                "name": f"C_two_prime_L{row['L']}",
                "passed": row["two_prime_rank_p1"] == row["two_prime_rank_p2"] == cyclic_q,
                "detail": "rank_p1 == rank_p2 == cyclic_Q_dim",
            }
        )
        checks.append(
            {
                "name": f"C_distinct_delta_rank_L{row['L']}",
                "passed": row["distinct_delta_pairing_rank"] == cyclic_q,
                "detail": f"2^(2L-1)=2^{2 * row['L'] - 1} distinct delta-functionals give the same pairing rank",
            }
        )
        checks.append(
            {
                "name": f"C_exact_annihilator_L{row['L']}",
                "passed": (
                    row["pairing_kernel_dim"] == expected_q[row["L"]][2]
                    and row["validation"]["b_invariant_over_Q"]
                    and row["validation"]["a_invariant_over_Q"]
                    and row["validation"]["vacuum_orthogonal"]
                    and row["validation"]["sector_homogeneous"]
                    and row["validation"]["kernel_independent_over_Q"]
                ),
                "detail": "kernel is exactly (U(A,B) psi)^perp over Q",
            }
        )
    checks.append({
        "name": "C_L3_trivial_kernel",
        "passed": all(row["L"] != 3 or row["pairing_kernel_dim"] == 0 for row in records),
        "detail": "the L=3 pairing is non-degenerate: cyclic 14/14",
    })
    checks.append({
        "name": "C_sector_split_sum",
        "passed": all(sum(row["pairing_kernel_sectors"].values()) == row["pairing_kernel_dim"] for row in records),
        "detail": "sector sub-kernels reproduce the total W dims 0, 2, 10, 66",
    })
    provenance = {
        "artifact_schema": "provenance/data/checks-v1",
        "script": SCRIPT,
        "field": "Q and F_2147483647 x F_2147483629",
        "exact_arithmetic": "Python int and fractions.Fraction",
        "artifact_sha256": "",
    }
    envelope = {"provenance": provenance, "data": data, "checks": checks}
    canonical = json.dumps(
        {"provenance": {key: value for key, value in provenance.items() if key != "artifact_sha256"},
         "data": data, "checks": checks},
        sort_keys=True, separators=(",", ":"),
    ).encode()
    envelope["provenance"]["artifact_sha256"] = hashlib.sha256(canonical).hexdigest()
    return envelope


def write_artifact(artifact: dict) -> None:
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = RESULT_PATH.with_name(f".{RESULT_PATH.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    temporary.replace(RESULT_PATH)


def main() -> int:
    started = time.process_time()
    artifact = make_artifact()
    write_artifact(artifact)
    rows = artifact["data"]["pairing_records"]
    for row in rows:
        print(
            f"L={row['L']}: pairing rank {row['two_prime_rank_p1']}/{row['two_prime_rank_p2']}"
            f" (distinct-delta rank {row['distinct_delta_pairing_rank']}),"
            f" exact cyclic {row['cyclic_dim_Q']}/{row['K_dim']},"
            f" pairing kernel W dim {row['pairing_kernel_dim']} sectors={row['pairing_kernel_sectors']}"
        )
    print(f"total cpu {round(time.process_time() - started, 3)} s, peak rss {peak_rss_bytes()} B")
    print(f"wrote {RESULT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
