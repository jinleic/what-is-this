#!/usr/bin/env python3
"""L = 8 extension of the wave-14 ladder sector-saturation certificate.

Setting (identical to `proofs/ladder_alll_proof.md` and the wave-14 fronts
`experiments/e131_sector_saturation_tensor.py` /
`experiments/e132_sector_saturation_pairing.py`): the open `2xL` ladder with
`n = 2L` sites, vacuum `psi = |+>^(2L)` (the zero configuration), and the two
diagonalising operators `A_L = sum_v X_v = 2L - 2N` (diagonal in the
X-eigenbasis) and `B_L = sum_edges Z_u Z_v` (integer in the orbit-sum basis of
the even-particle `<tau, rho>`-invariant space `K_L`).

The wave-14 certificates stop at `L = 7`:

    dim U(A,B) psi = K_L - W_L =  14, 42, 142, 494, 1780   (L = 3..7)
    W_L =                             0,  2,  10,  66,  364

against `dim K_L = 2^(2L-3) + 3*2^(L-2) = 14, 44, 152, 560, 2144`.  This
experiment pushes the same certificate to `L = 8`, where
`dim K_8 = 2^13 + 3*2^6 = 8384` in the `4^8 = 65536`-dimensional
configuration space (`32768` even configurations).

Certificate chain per `L` (the wave-14 sandwich, unchanged):

  (S1) *Modular lower rank.*  The closure of `psi` under `{A, 4B}` computed by
        sparse echelon elimination over `F_p` has rank `r_p`, and
        `dim_Q U(A,B) psi >= r_p` for *every* prime `p` (reduction modulo `p`
        cannot increase the rank of an integral family).  Two engines are
        used: the wave-14 sparse dictionary echelon at the 31-bit primes
        `P1 = 2147483647`, `P2 = 2147483629`, and a dense-vector / sparse-pivot
        numpy engine at the primes `Q1 = 999983`, `Q2 = 999979` (int64-safe:
        row nnz * (Q-1)^2 < 2^63).  The engines agree with each other on
        every `L` where both run.
  (S2) *Exact Q upper bound.*  The pairing kernel `W_L = (U(A,B) psi)^perp` is
        exhibited: for each non-pivot ("free") orbit coordinate the modular
        constraint system is back-substituted (the echelon rows are triangular
        because each row's pivot is its *maximum* coordinate), lifted to
        centred integers, and validated **over Q** with `fractions.Fraction`:
        vacuum-orthogonal, particle-sector homogeneous, linearly independent,
        and `4B`-invariant.  By the annihilator-recursion theorem of
        `proofs/sector_saturation_pairing.md` any `A,B`-invariant subspace of
        `K_L` orthogonal to `psi` lies in `(U(A,B) psi)^perp`, so
        `dim_Q U(A,B) psi <= K_L - dim W_L`.
  (S3) If `r_p + dim W_L = K_L` at a single prime, both inequalities are
        equalities: `dim_Q U(A,B) psi = K_L - dim W_L` **exactly**, and
        `W_L` is exactly the pairing kernel.  Agreement at further primes is
        redundancy against implementation error, not part of the proof.

`dim K_8` is computed independently three ways: the Burnside closed form, the
direct orbit enumeration of the rank-two reflection group `<tau, rho>` on even
configurations, and the per-particle-sector Burnside sum (fixed-set counts by
direct enumeration of the tau-, rho-, and tau rho-fixed even configurations).

If the `L = 8` dictionary-engine closure exceeds its CPU budget the run stops
with an explicit wall record (exact counts of what completed); the fast-engine
certificate (S1 at Q1, Q2 + S2 + S3) stands on its own because the sandwich
needs only one prime.

Every compute gate uses `time.process_time()` (process CPU), never wall time;
peak RSS is recorded.  Run under `nice -n 19` with pinned single-thread
environment.
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e142_sector_saturation_l8.py"
RESULT_PATH = ROOT / "results" / "ladder" / "l8_saturation.json"
P1 = 2_147_483_647
P2 = 2_147_483_629
Q1 = 999_983
Q2 = 999_979
Q3 = 999_953
LIFT_BOUND_LADDER = (128, 4096, 1 << 20)
REGRESSION_L = (3, 4, 5, 6)
EXPECTED = {  # wave-14 certified values (cyclic, K, W, sector split)
    3: (14, 14, 0, {}),
    4: (42, 44, 2, {"4": 2}),
    5: (142, 152, 10, {"4": 5, "6": 5}),
    6: (494, 560, 66, {"4": 15, "6": 36, "8": 15}),
    7: (1780, 2144, 364, {"4": 35, "6": 147, "8": 147, "10": 35}),
}
L7_DICT_CPU_BUDGET = 900.0
L8_FAST_CPU_BUDGET = float(os.environ.get("E142_L8_FAST_BUDGET", "5400"))
L8_DICT_CPU_BUDGET = float(os.environ.get("E142_L8_DICT_BUDGET", "5400"))
STAGE = os.environ.get("E142_STAGE", "all")


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def log(message: str) -> None:
    print(f"[e142 {time.strftime('%H:%M:%S')}] {message}", flush=True)


# ---------------------------------------------------------------------------
# Ladder configuration space and the Klein-four symmetry (wave-14 core).
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
    edges = ladder_edges(L)
    for orbit in members:
        counts: dict[int, int] = defaultdict(int)
        for config in orbit:
            for u, v in edges:
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


def cyclic_module_mod_p(L: int, P: int, cpu_budget: float) -> tuple[ModEchelon, list[int], list[tuple[int, ...]], float]:
    """The wave-14 dictionary-engine closure of psi under {A, 4B} over F_P."""
    orbit_of, members = invariant_orbits(L)
    B4 = four_B_rows(L, orbit_of, members, P)
    A_diagonal = [(2 * L - 2 * orbit[0].bit_count()) % P for orbit in members]
    echelon = ModEchelon(P)
    first = echelon.add({orbit_of[0]: 1})
    if first is None:
        raise AssertionError("vacuum did not seed the modular closure")
    todo = [first]
    processed: set[int] = set()
    started = time.process_time()
    next_report = 256
    while todo:
        pivot = todo.pop()
        if pivot in processed:
            continue
        processed.add(pivot)
        vector = echelon.rows[pivot]
        images = (
            {index: A_diagonal[index] * coefficient % P for index, coefficient in vector.items()},
            multiply(vector, B4, P),
        )
        for image in images:
            new_pivot = echelon.add(image)
            if new_pivot is not None:
                todo.append(new_pivot)
        if len(processed) >= next_report:
            log(f"dict L={L}: processed {len(processed)} rows, rank {echelon.rank}, cpu {time.process_time() - started:.1f}s")
            next_report = len(processed) * 2
        if time.process_time() - started > cpu_budget:
            raise TimeoutError(
                f"L={L}: dictionary closure exceeded {cpu_budget} process CPU seconds "
                f"at rank {echelon.rank} with {len(processed)} rows processed"
            )
    return echelon, orbit_of, members, time.process_time() - started


# ---------------------------------------------------------------------------
# Fast dense-vector / sparse-pivot closure over F_q (q small, int64-safe).
# ---------------------------------------------------------------------------
class FastClosure:
    """Closure of psi under {A, 4B} over F_q with dense numpy vectors.

    Semantics are identical to `cyclic_module_mod_p`: the action is the
    row-vector action `v -> v * 4B` (i.e. `(4B)^T v` in column convention) and
    `v -> v * A` (diagonal), the echelon keeps pivot-normalised sparse rows
    with the pivot the *maximum* coordinate, and the todo queue processes each
    stored row once.  Products stay below 2^63 because
    row_nnz * (q-1)^2 < 2^63 for q <= 999983.
    """

    def __init__(self, L: int, q: int) -> None:
        self.q = q
        self.orbit_of, self.members = invariant_orbits(L)
        n = len(self.members)
        self.n = n
        B4 = four_B_rows(L, self.orbit_of, self.members, None)
        indptr = np.zeros(n + 1, dtype=np.int64)
        indices: list[int] = []
        data: list[int] = []
        for source, row in enumerate(B4):
            for target, coefficient in sorted(row.items()):
                indices.append(target)
                data.append(coefficient % q)
            indptr[source + 1] = len(indices)
        matrix = csr_matrix(
            (np.array(data, dtype=np.int64), np.array(indices, dtype=np.int32), indptr),
            shape=(n, n),
        )
        self.MT = matrix.transpose().tocsr()  # (v * 4B)^T action
        self.A_diagonal = np.array(
            [(2 * L - 2 * orbit[0].bit_count()) % q for orbit in self.members], dtype=np.int64
        )
        self.rows: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        self.processed: set[int] = set()
        self.adds = 0

    def add(self, w: np.ndarray) -> int | None:
        """Reduce the dense vector w against the stored pivot rows.

        Exactly the wave-14 `ModEchelon.add` semantics with the pivot the
        *current maximum* nonzero coordinate re-found after every reduction,
        so every stored row's pivot is its maximum coordinate (the invariant
        the kernel back-substitution relies on).
        """
        q = self.q
        rows = self.rows
        self.adds += 1
        while True:
            nz = np.nonzero(w)[0]
            if nz.size == 0:
                return None
            pivot = int(nz[-1])
            stored = rows.get(pivot)
            if stored is None:
                inverse = pow(int(w[pivot]), -1, q)
                values = w[nz] * inverse % q
                rows[pivot] = (nz.astype(np.int32).copy(), values.copy())
                return pivot
            idx, vals = stored
            w[idx] = (w[idx] - int(w[pivot]) * vals) % q

    def run(self, cpu_budget: float, progress: int = 512) -> tuple[int, float, dict]:
        q = self.q
        w0 = np.zeros(self.n, dtype=np.int64)
        vacuum = self.orbit_of[0]
        w0[vacuum] = 1
        first = self.add(w0)
        if first is None:
            raise AssertionError("vacuum did not seed the fast closure")
        todo = [first]
        started = time.process_time()
        next_report = progress
        while todo:
            pivot = todo.pop()
            if pivot in self.processed:
                continue
            self.processed.add(pivot)
            idx, vals = self.rows[pivot]
            w = np.zeros(self.n, dtype=np.int64)
            w[idx] = vals
            b_image = self.MT.dot(w) % q
            new_pivot = self.add(b_image)
            if new_pivot is not None:
                todo.append(new_pivot)
            a_image = w * self.A_diagonal % q
            new_pivot = self.add(a_image)
            if new_pivot is not None:
                todo.append(new_pivot)
            if len(self.processed) >= next_report:
                log(f"fast q={q}: processed {len(self.processed)} rows, rank {len(self.rows)}, cpu {time.process_time() - started:.1f}s")
                next_report = len(self.processed) * 2
            if time.process_time() - started > cpu_budget:
                raise TimeoutError(
                    f"fast closure q={q} exceeded {cpu_budget} process CPU seconds "
                    f"at rank {len(self.rows)} with {len(self.processed)} rows processed"
                )
        stats = {"adds": self.adds, "rows_processed": len(self.processed)}
        return len(self.rows), time.process_time() - started, stats

    def rows_as_dicts(self) -> dict[int, dict[int, int]]:
        return {
            pivot: {int(i): int(v) for i, v in zip(idx, vals)}
            for pivot, (idx, vals) in self.rows.items()
        }


# ---------------------------------------------------------------------------
# Pairing-kernel lifts and exact Q validation (wave-14 core, prime-generic).
# ---------------------------------------------------------------------------
def kernel_lifts(echelon_rows: dict[int, dict[int, int]], member_count: int, sizes: list[int], modulus: int) -> list[dict[int, int]]:
    """Small integer lifts of the modular pairing-kernel parametrised by the free coordinates.

    Back-substitution runs over pivots in ascending order, which is exactly
    triangular because each stored row's pivot is its maximum coordinate.
    """
    constraints = {
        pivot: {index: coefficient * sizes[index] % modulus for index, coefficient in row.items()}
        for pivot, row in echelon_rows.items()
    }
    pivots = sorted(constraints)
    row_lookup = [(pivot, constraints[pivot]) for pivot in pivots]
    free_indices = [index for index in range(member_count) if index not in constraints]
    candidates: list[dict[int, int]] = []
    half = modulus // 2
    inverses: dict[int, int] = {}
    for free in free_indices:
        vector: dict[int, int] = {free: 1}
        for pivot, row in row_lookup:
            residual = sum(row.get(index, 0) * coefficient for index, coefficient in vector.items()) % modulus
            if residual:
                inverse = inverses.get(pivot)
                if inverse is None:
                    inverse = pow(row[pivot], modulus - 2, modulus)
                    inverses[pivot] = inverse
                vector[pivot] = -residual * inverse % modulus
        lifted: dict[int, int] = {}
        for index, residue in vector.items():
            integer = residue if residue <= half else residue - modulus
            if integer:
                lifted[index] = integer
        candidates.append(lifted)
    return candidates


def kernel_lifts_dense(echelon_rows: dict[int, dict[int, int]], member_count: int, sizes: list[int], modulus: int, bound: int) -> list[dict[int, int]]:
    """Vectorised version of `kernel_lifts` for the L = 8 scale.

    Solves the same ascending triangular system for all free coordinates at
    once with a dense (member_count x n_free) int64 tableau.  Raises if any
    centred coefficient leaves [-bound, bound].
    """
    pivots = sorted(echelon_rows)
    free_indices = [index for index in range(member_count) if index not in echelon_rows]
    n_free = len(free_indices)
    free_column = {index: column for column, index in enumerate(free_indices)}
    X = np.zeros((member_count, n_free), dtype=np.int64)
    for index, column in free_column.items():
        X[index, column] = 1
    row_cache = [
        (
            pivot,
            np.array([index for index in sorted(row) if index != pivot], dtype=np.int64),
            np.array([row[index] * sizes[index] % modulus for index in sorted(row) if index != pivot], dtype=np.int64),
            row[pivot] * sizes[pivot] % modulus,
        )
        for pivot, row in ((p, echelon_rows[p]) for p in pivots)
    ]
    inverses: dict[int, int] = {}
    half = modulus // 2
    for pivot, others, weights, diagonal in row_cache:
        if others.size:
            contribution = weights @ X[others]  # (n_free,) int64; |w|<q, |X|<q, row nnz <= 88
            X[pivot] = (-contribution) % modulus
        else:
            X[pivot] = 0
        inverse = inverses.get(pivot)
        if inverse is None:
            inverse = pow(diagonal, modulus - 2, modulus)
            inverses[pivot] = inverse
        X[pivot] = X[pivot] * inverse % modulus
    centred = np.where(X > half, X - modulus, X)
    if int(np.abs(centred).max(initial=0)) > bound:
        raise AssertionError(
            f"kernel lift exceeded coefficient bound {bound}: max |c| = {int(np.abs(centred).max())}"
        )
    candidates: list[dict[int, int]] = []
    rows, columns = np.nonzero(centred)
    per_column: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for r, c in zip(rows, columns):
        per_column[int(c)].append((int(r), int(centred[r, c])))
    for column in range(n_free):
        candidates.append({index: value for index, value in per_column.get(column, [])})
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


def validate_W_over_Q(L: int, basis: list[dict[int, int]], orbit_of: list[int],
                      members: list[tuple[int, ...]]) -> dict:
    B4 = four_B_rows(L, orbit_of, members, None)
    span = QSpan()
    for vector in basis:
        span.add(vector)
    sector_lists = [sorted({members[index][0].bit_count() for index in vector}) for vector in basis]
    outside = [index for index, vector in enumerate(basis) if span.residual(multiply(vector, B4, None))]
    a_diagonal = [2 * L - 2 * orbit[0].bit_count() for orbit in members]
    a_outside = [
        index for index, vector in enumerate(basis)
        if span.residual({i: a_diagonal[i] * c for i, c in vector.items() if a_diagonal[i] * c})
    ]
    result = {
        "unreachable_block_dim": span.rank,
        "sector_particle_counts": sector_lists,
        "sector_homogeneous": all(len(values) == 1 for values in sector_lists),
        "vacuum_orthogonal": all(orbit_of[0] not in vector for vector in basis),
        "basis_independent_over_Q": span.rank == len(basis),
        "four_B_invariant_over_Q": not outside,
        "A_invariant_over_Q": not a_outside,
        "outside_image_indices": outside,
    }
    if not all((result["sector_homogeneous"], result["vacuum_orthogonal"],
                result["basis_independent_over_Q"], result["four_B_invariant_over_Q"],
                result["A_invariant_over_Q"])):
        raise AssertionError(f"L={L}: candidate W failed exact validation: "
                             f"{ {k: v for k, v in result.items() if k != 'sector_particle_counts'} }")
    return result


def serialise_basis(basis: list[dict[int, int]], members: list[tuple[int, ...]]) -> list[list[dict[str, int]]]:
    return [
        [
            {"representative": members[index][0], "coefficient": coefficient}
            for index, coefficient in sorted(vector.items())
        ]
        for vector in basis
    ]


def lift_with_ladder(echelon_rows, member_count, sizes, modulus, dense: bool):
    """Try increasing lift bounds; validation over Q is what certifies, the
    bound is only an implementation contract (wave-14 used 64/128)."""
    last_error: Exception | None = None
    for bound in LIFT_BOUND_LADDER:
        try:
            if dense:
                return kernel_lifts_dense(echelon_rows, member_count, sizes, modulus, bound), bound
            candidates = kernel_lifts(echelon_rows, member_count, sizes, modulus)
            biggest = max((abs(value) for vector in candidates for value in vector.values()), default=0)
            if biggest > bound:
                raise AssertionError(f"lift coefficient {biggest} exceeds bound {bound}")
            return candidates, bound
        except AssertionError as error:
            last_error = error
            log(f"lift bound {bound} failed ({error}); escalating")
    raise last_error


# ---------------------------------------------------------------------------
# Independent K_L dimension accounting.
# ---------------------------------------------------------------------------
def K_dimension_record(L: int) -> dict:
    orbit_of, members = invariant_orbits(L)
    direct = len(members)
    formula = burnside_K_dim(L)
    sector_direct = Counter(orbit[0].bit_count() for orbit in members)
    fixed_tau = Counter()
    fixed_rho = Counter()
    fixed_tau_rho = Counter()
    for config in range(1 << (2 * L)):
        if config.bit_count() & 1:
            continue
        weight = config.bit_count()
        if tau(config, L) == config:
            fixed_tau[weight] += 1
        if rho(config, L) == config:
            fixed_rho[weight] += 1
        if tau(rho(config, L), L) == config:
            fixed_tau_rho[weight] += 1
    sector_burnside = {
        str(k): (sum(1 for config in range(1 << (2 * L)) if config.bit_count() == k)
                 + fixed_tau[k] + fixed_rho[k] + fixed_tau_rho[k]) // 4
        for k in range(0, 2 * L + 1, 2)
    }
    if direct != formula or dict(sorted(sector_direct.items())) != {int(k): v for k, v in sector_burnside.items()}:
        raise AssertionError(
            f"L={L}: K accounting mismatch: direct {direct}, formula {formula}, "
            f"sector direct {dict(sorted(sector_direct.items()))}, sector Burnside {sector_burnside}"
        )
    return {
        "L": L,
        "K_dim_formula": formula,
        "K_dim_direct_orbits": direct,
        "K_dim_sector_burnside": sector_burnside,
        "K_dim_sector_direct": {str(k): sector_direct[k] for k in sorted(sector_direct)},
        "fixed_tau_by_sector": {str(k): fixed_tau[k] for k in sorted(fixed_tau)},
        "fixed_rho_by_sector": {str(k): fixed_rho[k] for k in sorted(fixed_rho)},
        "fixed_tau_rho_by_sector": {str(k): fixed_tau_rho[k] for k in sorted(fixed_tau_rho)},
    }


# ---------------------------------------------------------------------------
# Per-L certificate records.
# ---------------------------------------------------------------------------
def dict_record(L: int, primes: tuple[int, ...], cpu_budget: float) -> dict:
    started = time.process_time()
    ranks = []
    lift_echelon = None
    orbit_of = None
    members = None
    closure_cpu = 0.0
    for index, prime in enumerate(primes):
        echelon, orbit_of, members, closure_cpu = cyclic_module_mod_p(L, prime, cpu_budget)
        ranks.append(echelon.rank)
        if index == 0:
            lift_echelon = echelon  # lifts must use the SAME prime as the echelon rows
        log(f"dict L={L} p={prime}: rank {echelon.rank} in {closure_cpu:.1f}s cpu")
    sizes = [len(orbit) for orbit in members]
    candidates, bound = lift_with_ladder(lift_echelon.rows, len(members), sizes, primes[0], dense=False)
    validation = validate_W_over_Q(L, candidates, orbit_of, members)
    cyclic_q = len(members) - validation["unreachable_block_dim"]
    if cyclic_q != ranks[0]:
        raise AssertionError(f"L={L}: exact upper {cyclic_q} misses modular rank {ranks[0]}")
    sectors = Counter(values[0] for values in validation["sector_particle_counts"])
    return {
        "L": L,
        "engine": "dict",
        "K_dim": len(members),
        "mod_p_ranks": ranks,
        "cyclic_Q_dim": cyclic_q,
        "W_dim": validation["unreachable_block_dim"],
        "W_sectors": {str(k): sectors[k] for k in sorted(sectors)},
        "lift_bound": bound,
        "basis": serialise_basis(candidates, members),
        "validation": {k: v for k, v in validation.items() if k != "sector_particle_counts"},
        "closure_cpu_seconds": round(closure_cpu, 3),
        "record_cpu_seconds": round(time.process_time() - started, 3),
        "claim_tag": "[COMPUTATION]",
    }


def fast_record(L: int, primes: tuple[int, ...], cpu_budget: float, serialise: bool) -> dict:
    started = time.process_time()
    ranks = []
    lift_closure = None
    for index, prime in enumerate(primes):
        closure = FastClosure(L, prime)
        try:
            rank, cpu, stats = closure.run(cpu_budget)
        except TimeoutError as error:
            log(f"fast L={L} q={prime} walled: {error}")
            raise
        ranks.append(rank)
        if index == 0:
            lift_closure = closure  # lifts must use the SAME prime as the echelon rows
        log(f"fast L={L} q={prime}: rank {rank} in {cpu:.1f}s cpu ({stats})")
    sizes = [len(orbit) for orbit in lift_closure.members]
    candidates, bound = lift_with_ladder(lift_closure.rows_as_dicts(), len(lift_closure.members), sizes, primes[0], dense=True)
    validation = validate_W_over_Q(L, candidates, lift_closure.orbit_of, lift_closure.members)
    cyclic_q = len(lift_closure.members) - validation["unreachable_block_dim"]
    if cyclic_q != ranks[0]:
        raise AssertionError(f"L={L}: exact upper {cyclic_q} misses fast modular rank {ranks[0]}")
    sectors = Counter(values[0] for values in validation["sector_particle_counts"])
    record = {
        "L": L,
        "engine": "dense-numpy",
        "K_dim": len(closure.members),
        "mod_p_ranks": ranks,
        "cyclic_Q_dim": cyclic_q,
        "W_dim": validation["unreachable_block_dim"],
        "W_sectors": {str(k): sectors[k] for k in sorted(sectors)},
        "lift_bound": bound,
        "validation": {k: v for k, v in validation.items() if k != "sector_particle_counts"},
        "record_cpu_seconds": round(time.process_time() - started, 3),
        "claim_tag": "[COMPUTATION]",
    }
    if serialise:
        record["basis"] = serialise_basis(candidates, closure.members)
    return record


def l8_full_record(K8: dict) -> dict:
    """The L = 8 certificate: fast engine at two primes, dictionary engine at P1."""
    record: dict = {
        "L": 8,
        "K_dim": K8["K_dim_formula"],
        "claim_tag": "[COMPUTATION]",
    }
    fast = fast_record(8, (Q1, Q2), L8_FAST_CPU_BUDGET, serialise=True)
    record["fast"] = fast
    wall = None
    try:
        echelon, orbit_of, members, closure_cpu = cyclic_module_mod_p(8, P1, L8_DICT_CPU_BUDGET)
        dict_rank = echelon.rank
        log(f"dict L=8 p={P1}: rank {dict_rank} in {closure_cpu:.1f}s cpu")
        if dict_rank != fast["cyclic_Q_dim"]:
            raise AssertionError(f"L=8: dict rank {dict_rank} != fast rank {fast['cyclic_Q_dim']}")
        record["dict_confirm"] = {
            "prime": P1,
            "rank": dict_rank,
            "closure_cpu_seconds": round(closure_cpu, 3),
            "agrees_with_fast": True,
        }
    except TimeoutError as error:
        wall = {
            "engine": "dict",
            "prime": P1,
            "detail": str(error),
            "budget_cpu_seconds": L8_DICT_CPU_BUDGET,
            "note": "wall record only; the fast-engine sandwich (S1 at multiple primes + S2 exact-Q validation) is self-contained",
        }
        log(f"dict L=8 walled: {error}")
        record["dict_confirm"] = wall
        # Third small prime as replacement redundancy for the walled engine.
        rank3, cpu3, _ = FastClosure(8, Q3).run(L8_FAST_CPU_BUDGET)
        log(f"third-prime confirmation q={Q3}: rank {rank3}")
        if rank3 != fast["cyclic_Q_dim"]:
            raise AssertionError(f"L=8: third prime rank {rank3} != {fast['cyclic_Q_dim']}")
        record["third_prime_confirm"] = {"prime": Q3, "rank": rank3, "cpu_seconds": round(cpu3, 3)}
    # Merge headline numbers from the fast certificate.
    for key in ("engine", "mod_p_ranks", "cyclic_Q_dim", "W_dim", "W_sectors", "lift_bound", "basis", "validation", "record_cpu_seconds"):
        record[key] = fast[key]
    extra_ranks = []
    if "rank" in record.get("dict_confirm", {}):
        extra_ranks.append(record["dict_confirm"]["rank"])
    if "rank" in record.get("third_prime_confirm", {}):
        extra_ranks.append(record["third_prime_confirm"]["rank"])
    record["primes_all_agree"] = fast["mod_p_ranks"] + extra_ranks
    return record


def pattern_analysis(records: list[dict]) -> dict:
    table = []
    for record in sorted(records, key=lambda r: r["L"]):
        K = record["K_dim"]
        W = record["W_dim"]
        table.append({
            "L": record["L"],
            "K_dim": K,
            "cyclic_dim": record["cyclic_Q_dim"],
            "W": W,
            "W_over_K": round(W / K, 6),
            "cyclic_geq_2powL": record["cyclic_Q_dim"] >= 2 ** record["L"],
            "rescue_slack_K_minus_W_minus_2powL": K - W - 2 ** record["L"],
        })
    W_values = [row["W"] for row in table]
    return {
        "table": table,
        "W_sequence": W_values,
        "W_strictly_increasing": all(b > a for a, b in zip(W_values, W_values[1:])),
        "growth_ratios": [round(b / a, 6) for a, b in zip(W_values, W_values[1:]) if a],
        "claim_tag": "[COMPUTATION]",
    }


def make_artifact() -> dict:
    started = time.process_time()
    K_records = {L: K_dimension_record(L) for L in (7, 8)}
    regression = []
    for L in REGRESSION_L:
        regression.append(dict_record(L, (P1, P2), 300.0))
    l7 = dict_record(7, (P1,), L7_DICT_CPU_BUDGET)
    # Cross-validate the fast engine against the certified dictionary engine.
    for L in (5, 6):
        fast_rank, _, _ = FastClosure(L, Q1).run(600.0)
        if fast_rank != EXPECTED[L][0]:
            raise AssertionError(f"fast engine disagrees at L={L}: {fast_rank} != {EXPECTED[L][0]}")
        log(f"fast engine cross-check L={L}: rank {fast_rank} matches wave-14")
    K8 = K_records[8]
    if K8["K_dim_formula"] != 8384:
        raise AssertionError(f"dim K_8 != 8384: {K8}")
    l8 = l8_full_record(K8)

    for record in regression + [l7]:
        expected = EXPECTED[record["L"]]
        if (record["cyclic_Q_dim"], record["K_dim"], record["W_dim"]) != expected[:3]:
            raise AssertionError(f"L={record['L']}: {record['cyclic_Q_dim']}/{record['K_dim']} W={record['W_dim']} != wave-14 {expected}")
    checks = []
    for record in regression + [l7]:
        L = record["L"]
        expected = EXPECTED[L]
        checks.append({
            "name": f"C_regression_L{L}",
            "passed": (record["cyclic_Q_dim"], record["K_dim"], record["W_dim"]) == expected[:3]
                      and {k: v for k, v in record["W_sectors"].items()} == expected[3],
            "detail": f"cyclic {record['cyclic_Q_dim']}/{record['K_dim']}, W {record['W_dim']}, sectors {record['W_sectors']}",
        })
        checks.append({
            "name": f"C_sandwich_L{L}",
            "passed": all(r == record["cyclic_Q_dim"] for r in record["mod_p_ranks"])
                      and record["validation"]["four_B_invariant_over_Q"]
                      and record["validation"]["A_invariant_over_Q"]
                      and record["validation"]["vacuum_orthogonal"]
                      and record["validation"]["sector_homogeneous"]
                      and record["validation"]["basis_independent_over_Q"],
            "detail": "modular lower rank equals exact Q upper bound K - dim W",
        })
    sectors8 = l8["W_sectors"]
    symmetric = all(sectors8.get(str(16 - int(k))) == v for k, v in sectors8.items())
    checks.extend([
        {
            "name": "C_K8_three_ways",
            "passed": K8["K_dim_formula"] == K8["K_dim_direct_orbits"] == 8384
                      and sum(K8["K_dim_sector_burnside"].values()) == 8384
                      and sum(K8["K_dim_sector_direct"].values()) == 8384,
            "detail": "dim K_8 = 2^13 + 3*2^6 = 8384 by Burnside formula, direct orbit count, per-sector Burnside",
        },
        {
            "name": "C_L8_sandwich",
            "passed": all(r == l8["cyclic_Q_dim"] for r in l8["mod_p_ranks"])
                      and l8["K_dim"] - l8["W_dim"] == l8["cyclic_Q_dim"]
                      and l8["validation"]["four_B_invariant_over_Q"]
                      and l8["validation"]["A_invariant_over_Q"]
                      and l8["validation"]["vacuum_orthogonal"]
                      and l8["validation"]["sector_homogeneous"]
                      and l8["validation"]["basis_independent_over_Q"],
            "detail": f"cyclic {l8['cyclic_Q_dim']}/8384 exact over Q; ranks agree at primes {l8['mod_p_ranks']}",
        },
        {
            "name": "C_L8_sector_symmetry",
            "passed": symmetric,
            "detail": f"W_8 sector split {sectors8} is symmetric under k -> 16 - k",
        },
        {
            "name": "C_L8_pattern_verdict_consistent",
            "passed": (l8["W_dim"] > 364) == all(
                b > a for a, b in zip([0, 2, 10, 66, 364, l8["W_dim"]], [2, 10, 66, 364, l8["W_dim"]])
            ),
            "detail": f"W_8 = {l8['W_dim']}; empirical growth "
                      f"{'continues' if l8['W_dim'] > 364 else 'breaks'} the strictly increasing sequence 0, 2, 10, 66, 364",
        },
        {
            "name": "C_L8_evaluation_bound",
            "passed": l8["cyclic_Q_dim"] >= 2 ** 8,
            "detail": f"dim g_8 >= dim U(A,B) psi = {l8['cyclic_Q_dim']} >= 2^8 = 256",
        },
    ])
    if isinstance(l8["dict_confirm"], dict) and "rank" in l8["dict_confirm"]:
        checks.append({
            "name": "C_L8_dict_engine_agreement",
            "passed": l8["dict_confirm"]["rank"] == l8["cyclic_Q_dim"],
            "detail": f"wave-14 dictionary engine at p={P1} reproduces rank {l8['dict_confirm']['rank']}",
        })
    else:
        checks.append({
            "name": "C_L8_dict_engine_wall_recorded",
            "passed": "detail" in l8["dict_confirm"],
            "detail": l8["dict_confirm"].get("detail", ""),
        })
    data = {
        "scope": {
            "space": "K_L = even-parity <tau,rho>-invariant ladder configuration space; X-basis particle grading A_L = 2L-2N, B_L = D+F+Ddag",
            "field": "Q for the W certificates; F_p ranks at p in [P1, P2] = [2147483647, 2147483629] (dictionary engine) and [Q1, Q2] = [999983, 999979] (dense engine)",
            "certificate": "modular lower rank (any one prime) + exact-Q validated pairing kernel W (annihilator recursion) => equality",
            "regression_scope": "L = 3..7 reproduce the wave-14 values 14/14, 42/44, 142/152, 494/560, 1780/2144 with W = 0, 2, 10, 66, 364",
        },
        "K_dimension_records": {str(L): K_records[L] for L in K_records},
        "regression_records": regression,
        "L7_record": l7,
        "L8_record": l8,
        "pattern": pattern_analysis(regression + [l7, {"L": 8, "K_dim": l8["K_dim"], "cyclic_Q_dim": l8["cyclic_Q_dim"], "W_dim": l8["W_dim"]}]),
        "resource_measurements": {
            "total_cpu_seconds": round(time.process_time() - started, 3),
            "peak_rss_bytes": peak_rss_bytes(),
            "budget_clock": "time.process_time",
        },
    }
    provenance = {
        "artifact_schema": "provenance/data/checks-v1",
        "script": SCRIPT,
        "field": "Q and F_2147483647 x F_2147483629 x F_999983 x F_999979",
        "exact_arithmetic": "Python int, numpy int64 (mod q <= 999983), fractions.Fraction",
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


def pilot() -> int:
    """Fast-engine measurement without artifact side effects."""
    for L in (5, 6, 7):
        closure = FastClosure(L, Q1)
        rank, cpu, stats = closure.run(3600.0)
        log(f"pilot fast L={L}: rank {rank} (expected {EXPECTED[L][0]}) in {cpu:.1f}s cpu, {stats}")
        if rank != EXPECTED[L][0]:
            raise AssertionError(f"fast engine disagrees at L={L}")
    closure = FastClosure(8, Q1)
    rank, cpu, stats = closure.run(L8_FAST_CPU_BUDGET)
    log(f"pilot fast L=8 q={Q1}: rank {rank} of K_8=8384 in {cpu:.1f}s cpu, {stats}")
    sizes = [len(orbit) for orbit in closure.members]
    started = time.process_time()
    candidates, bound = lift_with_ladder(closure.rows_as_dicts(), len(closure.members), sizes, Q1, dense=True)
    log(f"pilot L=8 lifts: {len(candidates)} candidates, bound {bound}, "
        f"max |c| = {max((abs(v) for c in candidates for v in c.values()), default=0)}, "
        f"{time.process_time() - started:.1f}s cpu")
    validation = validate_W_over_Q(8, candidates, closure.orbit_of, closure.members)
    log(f"pilot L=8 validation over Q: { {k: v for k, v in validation.items() if k != 'sector_particle_counts'} }")
    return 0


def main() -> int:
    if STAGE == "pilot":
        return pilot()
    started = time.process_time()
    artifact = make_artifact()
    write_artifact(artifact)
    l8 = artifact["data"]["L8_record"]
    print(f"L=8: cyclic {l8['cyclic_Q_dim']}/{l8['K_dim']}, W_8 = {l8['W_dim']}, sectors {l8['W_sectors']}")
    print(f"pattern: {artifact['data']['pattern']['W_sequence']}")
    print(f"total cpu {round(time.process_time() - started, 3)} s, peak rss {peak_rss_bytes()} B")
    print(f"wrote {RESULT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
