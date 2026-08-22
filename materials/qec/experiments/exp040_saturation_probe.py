"""EXP-040: forced-saturation probe for PBB dressing modules.

The certified-reversal arm rebuilds every named catalogue PBB from its stored
polynomials and measures the actual image

    Delta = { lambda [C D] : lambda [A B] = 0 }

rather than inferring it from the child k alone.  It then quotients by the
parent Z-stabilizer space and checks the exact EXP-039 minimum-logical module
reconstructed from the persisted witness-orbit basis.

The small-lattice arm is deliberately solver-free.  It samples connected,
positive-k BB parents on lattices with ell*m <= 12, exhaustively enumerates
commuting perturbations of total C/D support weight at most four, computes the
parent d_Z and T from the full GF(2) nullspace span, and performs an exact
meet-in-the-middle search for every candidate child that absorbs the complete
minimum-logical module.  The latter decides the threshold predicate
``d_Q > d_Z(P)`` without invoking a SAT solver.
"""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import os
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from math import comb
from pathlib import Path
from typing import Any, Iterable

# This experiment shares a workstation with live SAT shards.  All of its work
# is exact GF(2) enumeration, and NumPy/Accelerate must stay single-threaded.
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_env, "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SCHEMA = "exp040-saturation-probe-v1"
OUT = ROOT / "results" / "processed" / "exp040_saturation_probe.json"
SEED = 20260817
MAX_PERTURBATION_SUPPORT = 4
PARENTS_PER_LATTICE = 16
MAX_PARENT_KERNEL_DIMENSION = 20

# The requested trinomial family has no positive-k parent on 5x2.  There the
# binomial fallback is explicitly tagged rather than silently treating it as a
# trinomial sample.
LATTICE_PROTOCOL = (
    {
        "ell": 4,
        "m": 3,
        "primary_weights": (3, 3),
        "fallback_weights": None,
    },
    {
        "ell": 3,
        "m": 4,
        "primary_weights": (3, 3),
        "fallback_weights": None,
    },
    {
        "ell": 5,
        "m": 2,
        "primary_weights": (3, 3),
        "fallback_weights": (2, 2),
    },
    {
        "ell": 6,
        "m": 2,
        "primary_weights": (3, 3),
        "fallback_weights": None,
    },
)

_SPEC_039 = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
)
E39 = importlib.util.module_from_spec(_SPEC_039)
sys.modules[_SPEC_039.name] = E39
_SPEC_039.loader.exec_module(E39)
E27 = E39.E27

from qec_research.codes.bicycle import (  # noqa: E402
    BBSpec,
    bb_stabilizer,
    commutation_defect,
    monomial_matrix,
    poly_matrix,
)
from qec_research.codes.pbb_survival import translation_orbit  # noqa: E402
from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402


# ---------------------------------------------------------------------------
# Small exact GF(2) / bitset helpers
# ---------------------------------------------------------------------------
def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def row_mask(row: np.ndarray) -> int:
    """Pack a short GF(2) row with column j at integer bit j."""

    bits = np.flatnonzero(np.asarray(row, dtype=np.uint8).reshape(-1) & 1)
    out = 0
    for bit in bits:
        out |= 1 << int(bit)
    return out


def mask_row(mask: int, ncols: int) -> np.ndarray:
    return np.fromiter(
        ((int(mask) >> j) & 1 for j in range(ncols)),
        dtype=np.uint8,
        count=ncols,
    )


def rows_to_masks(M: np.ndarray) -> list[int]:
    return [row_mask(row) for row in np.asarray(M, dtype=np.uint8)]


class Echelon:
    """Exact GF(2) row-space basis using Python integer bitsets."""

    def __init__(self, rows: Iterable[int] = ()) -> None:
        self._pivots: dict[int, int] = {}
        for row in rows:
            self.insert(int(row))

    def copy(self) -> "Echelon":
        out = Echelon()
        out._pivots = dict(self._pivots)
        return out

    def reduce(self, value: int) -> int:
        value = int(value)
        while value:
            pivot = value.bit_length() - 1
            basis_row = self._pivots.get(pivot)
            if basis_row is None:
                break
            value ^= basis_row
        return value

    def insert(self, value: int) -> bool:
        reduced = self.reduce(value)
        if not reduced:
            return False
        self._pivots[reduced.bit_length() - 1] = reduced
        return True

    def contains(self, value: int) -> bool:
        return self.reduce(value) == 0

    @property
    def rank(self) -> int:
        return len(self._pivots)


def enumerate_span(basis: np.ndarray) -> np.ndarray:
    """Every vector in the GF(2) span of independent basis rows, packed."""

    masks = np.zeros(1, dtype=np.uint32)
    for row in np.asarray(basis, dtype=np.uint8):
        b = np.uint32(row_mask(row))
        masks = np.concatenate((masks, masks ^ b))
    return masks


def orbit_masks(mask: int, ell: int, m: int) -> list[int]:
    """All diagonal lattice translates of a two-block Z-sector vector."""

    return [
        row_mask(row)
        for row in translation_orbit(mask_row(mask, 2 * ell * m), ell, m)
    ]


def gf2_product(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Small exact GF(2) matrix product without dtype-overflow ambiguity."""

    return (
        np.asarray(left, dtype=np.uint8).astype(np.int64)
        @ np.asarray(right, dtype=np.uint8).astype(np.int64)
        % 2
    ).astype(np.uint8)


def json_terms(terms: Iterable[tuple[int, int]]) -> list[list[int]]:
    return [[int(a), int(b)] for a, b in terms]


def support_terms(
    combo: Iterable[int], ell: int, m: int
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    """Split a total [C D] support into canonical C and D term tuples."""

    block = ell * m
    c_terms: list[tuple[int, int]] = []
    d_terms: list[tuple[int, int]] = []
    for coordinate in combo:
        if coordinate < block:
            c_terms.append((coordinate // m, coordinate % m))
        else:
            q = coordinate - block
            d_terms.append((q // m, q % m))
    return tuple(c_terms), tuple(d_terms)


# ---------------------------------------------------------------------------
# Certified reversal reconstruction
# ---------------------------------------------------------------------------
def pbb_block_from_row(row: dict[str, Any]) -> np.ndarray:
    ell, m = int(row["ell"]), int(row["m"])
    c = poly_matrix(ell, m, E27.terms(row, "C_terms"))
    d = poly_matrix(ell, m, E27.terms(row, "D_terms"))
    return np.hstack((c, d)).astype(np.uint8)


def dressing_image(HX: np.ndarray, CD: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return a left-kernel basis and rows spanning Delta = leftker(HX) CD."""

    left_kernel = nullspace_np(np.asarray(HX, dtype=np.uint8).T)
    if left_kernel.shape[0] == 0:
        return left_kernel, np.zeros((0, CD.shape[1]), dtype=np.uint8)
    return left_kernel, gf2_product(left_kernel, CD)


def reversal_records() -> tuple[list[dict[str, Any]], bool, list[str]]:
    """Rebuild the seven named reversals and test the forced-saturation rank."""

    rows = E27.load_catalogue()
    certs = E39.load_certificates()
    lower_bounds = E39.certified_pbb_lower_bounds()
    targets = set(E39.CERTIFIED_REVERSALS)
    records: list[dict[str, Any]] = []
    problems: list[str] = []

    for index, row in enumerate(rows):
        label = E27.catalogue_label(row, index)
        if label not in targets:
            continue

        pspec, code = E27.pbb_code(row)
        spec, HX, HZ = E27.parent_matrices(row)
        del pspec, spec
        ell, m = int(row["ell"]), int(row["m"])
        dim = ell * m
        n = 2 * dim
        parent_fp = E27.matrix_fingerprint(HX, HZ)
        cert = certs.get(parent_fp)
        if cert is None:
            raise RuntimeError(f"{label}: missing EXP-039 parent certificate")
        if not cert.get("T_is_exact"):
            raise RuntimeError(f"{label}: parent T is not exact")

        CD = pbb_block_from_row(row)
        H = np.asarray(code.H, dtype=np.uint8)
        if H.shape != (2 * dim, 2 * n):
            raise AssertionError(f"{label}: unexpected PBB check shape {H.shape}")
        if not np.array_equal(H[:dim, :n], HX):
            raise AssertionError(f"{label}: rebuilt [A B] disagrees with PBB block")
        if not np.array_equal(H[:dim, n:], CD):
            raise AssertionError(f"{label}: rebuilt [C D] disagrees with PBB block")
        if H[dim:, :n].any() or not np.array_equal(H[dim:, n:], HZ):
            raise AssertionError(f"{label}: rebuilt H_Z disagrees with PBB block")

        left_kernel, delta_rows = dressing_image(HX, CD)
        dim_ker_map = int(left_kernel.shape[0])
        dim_delta = int(rank_np(delta_rows))
        dim_s_z = int(rank_np(HZ))
        span_delta_s_z = np.vstack((HZ, delta_rows))
        dim_delta_bar = int(rank_np(span_delta_s_z) - dim_s_z)
        k_parent = int(n - rank_np(HX) - rank_np(HZ))
        k_child = int(code.n - rank_np(H))

        module_rows = HZ.copy()
        for witness in cert.get("witness_vectors", []):
            witness_array = np.asarray(witness, dtype=np.uint8).reshape(-1)
            if witness_array.shape[0] != n:
                raise AssertionError(f"{label}: malformed certificate witness")
            module_rows = np.vstack((module_rows, translation_orbit(witness_array, ell, m)))
        t_from_persisted_orbits = int(rank_np(module_rows) - dim_s_z)
        module_contained = (
            int(rank_np(np.vstack((span_delta_s_z, module_rows))))
            == int(rank_np(span_delta_s_z))
        )

        t_cert = int(cert["T"])
        equality = dim_delta_bar == t_cert
        orbit_record_consistent = t_from_persisted_orbits == t_cert
        structural_saturation = bool(
            equality and module_contained and orbit_record_consistent
        )
        d_z_parent = int(cert["d_z_parent"])
        d_q_lower_gt = lower_bounds.get(label)
        reversal_hypothesis_certified = bool(
            d_q_lower_gt is not None and int(d_q_lower_gt) >= d_z_parent
        )

        record = {
            "label": label,
            "catalogue_index": int(index),
            "parent_fingerprint": parent_fp,
            "ell": ell,
            "m": m,
            "n": n,
            "C_terms": row["C_terms"],
            "D_terms": row["D_terms"],
            "dim_ker_A_B_map": dim_ker_map,
            "dim_Delta": dim_delta,
            "dim_S_Z": dim_s_z,
            "dim_Delta_bar": dim_delta_bar,
            "dim_Delta_plus_S_Z_minus_dim_S_Z": dim_delta_bar,
            "T": t_cert,
            "T_is_exact": bool(cert["T_is_exact"]),
            "k_P": k_parent,
            "k_Q": k_child,
            "catalogue_k_Q": int(row["k"]),
            "k_drop": k_parent - k_child,
            "d_Z_parent_from_exact_certificate": d_z_parent,
            "certified_d_Q_gt": (
                int(d_q_lower_gt) if d_q_lower_gt is not None else None
            ),
            "reversal_hypothesis_certified": reversal_hypothesis_certified,
            "equality_dim_Delta_bar_eq_T": equality,
            "T_from_persisted_witness_orbits": t_from_persisted_orbits,
            "persisted_orbit_module_matches_T": orbit_record_consistent,
            "M_bar_contained_in_Delta_bar": module_contained,
            "Delta_bar_equals_M_bar": structural_saturation,
        }
        records.append(record)

        if not equality:
            problems.append(
                f"{label}: dim(Delta_bar)={dim_delta_bar} != T={t_cert}"
            )
        if k_child != int(row["k"]):
            problems.append(
                f"{label}: recomputed k_Q={k_child} != catalogue k={row['k']}"
            )
        if k_child != k_parent - dim_delta_bar:
            problems.append(
                f"{label}: k_Q != k_P - dim(Delta_bar)"
            )
        if not orbit_record_consistent:
            problems.append(f"{label}: persisted witness orbits do not reproduce T")
        if not module_contained:
            problems.append(f"{label}: exact minimum-logical module is not absorbed")
        if not reversal_hypothesis_certified:
            problems.append(f"{label}: strict-distance hypothesis is not certified")

    found = {record["label"] for record in records}
    missing = targets - found
    if missing:
        problems.append(f"missing certified reversals: {sorted(missing)}")
    records.sort(key=lambda record: E39.CERTIFIED_REVERSALS.index(record["label"]))
    return records, not problems, problems


# ---------------------------------------------------------------------------
# Small-lattice parent enumeration and exact d_Z, T(P)
# ---------------------------------------------------------------------------
def all_supports(dim: int, weight: int) -> list[tuple[int, ...]]:
    return list(itertools.combinations(range(dim), weight))


def parent_candidates(
    ell: int, m: int, weights: tuple[int, int], family_name: str
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Connected, positive-k BB parents from one unordered support family."""

    dim = ell * m
    monomials = [(a, b) for a in range(ell) for b in range(m)]
    wa, wb = weights
    supports_a = all_supports(dim, wa)
    supports_b = supports_a if wa == wb else all_supports(dim, wb)
    matrices_a = [
        poly_matrix(ell, m, [monomials[position] for position in support])
        for support in supports_a
    ]
    matrices_b = (
        matrices_a
        if wa == wb
        else [
            poly_matrix(ell, m, [monomials[position] for position in support])
            for support in supports_b
        ]
    )

    candidates: list[dict[str, Any]] = []
    k_histogram: Counter[int] = Counter()
    connected_positive_k = 0
    kernel_dimension_cap_skips = 0
    raw_pairs = 0
    for ai, A in enumerate(matrices_a):
        begin = ai if wa == wb else 0
        for bi in range(begin, len(matrices_b)):
            raw_pairs += 1
            B = matrices_b[bi]
            HX = np.hstack((A, B)).astype(np.uint8)
            HZ = np.hstack((B.T, A.T)).astype(np.uint8)
            rank_x, rank_z = int(rank_np(HX)), int(rank_np(HZ))
            k_parent = int(2 * dim - rank_x - rank_z)
            k_histogram[k_parent] += 1
            if k_parent < 2:
                continue

            a_terms = tuple(monomials[position] for position in supports_a[ai])
            b_terms = tuple(monomials[position] for position in supports_b[bi])
            parent_code = bb_stabilizer(
                BBSpec(ell=ell, m=m, A=list(a_terms), B=list(b_terms))
            )
            if len(parent_code.connected_components()) != 1:
                continue
            connected_positive_k += 1

            kernel_dimension = int(nullspace_np(HX).shape[0])
            if kernel_dimension > MAX_PARENT_KERNEL_DIMENSION:
                kernel_dimension_cap_skips += 1
                continue
            candidates.append(
                {
                    "ell": ell,
                    "m": m,
                    "A": a_terms,
                    "B": b_terms,
                    "HX": HX,
                    "HZ": HZ,
                    "k_P": k_parent,
                    "rank_HX": rank_x,
                    "rank_HZ": rank_z,
                    "kernel_dimension": kernel_dimension,
                    "family": family_name,
                }
            )

    accounting = {
        "family": family_name,
        "A_weight": wa,
        "B_weight": wb,
        "raw_unordered_support_pairs": raw_pairs,
        "k_parent_histogram": {
            str(k): int(v) for k, v in sorted(k_histogram.items())
        },
        "connected_positive_k_pairs": connected_positive_k,
        "kernel_dimension_cap_skips": kernel_dimension_cap_skips,
        "enumerable_connected_positive_k_pairs": len(candidates),
        "nondegenerate_predicate": (
            "one Tanner connected component, k_P >= 2, and right-nullspace "
            f"dimension <= {MAX_PARENT_KERNEL_DIMENSION} for exhaustive GF(2) span"
        ),
    }
    return candidates, accounting


def stratified_sample(
    candidates: list[dict[str, Any]], *, ell: int, m: int
) -> tuple[list[dict[str, Any]], int]:
    """Deterministic round-robin over (k_P, family) strata, then fill."""

    seed = int(
        np.random.SeedSequence([SEED, ell, m]).generate_state(1, dtype=np.uint32)[0]
    )
    rng = np.random.default_rng(seed)
    strata: dict[tuple[int, str], list[dict[str, Any]]] = defaultdict(list)
    for candidate in candidates:
        strata[(int(candidate["k_P"]), str(candidate["family"]))].append(candidate)
    buckets: list[list[dict[str, Any]]] = []
    for key in sorted(strata):
        bucket = list(strata[key])
        rng.shuffle(bucket)
        buckets.append(bucket)

    target = min(PARENTS_PER_LATTICE, len(candidates))
    chosen: list[dict[str, Any]] = []
    cursors = [0] * len(buckets)
    while len(chosen) < target:
        progressed = False
        for bucket_index, bucket in enumerate(buckets):
            if len(chosen) >= target:
                break
            cursor = cursors[bucket_index]
            if cursor < len(bucket):
                chosen.append(bucket[cursor])
                cursors[bucket_index] += 1
                progressed = True
        if not progressed:
            break
    return chosen, seed


def analyse_parent_exact(parent: dict[str, Any]) -> dict[str, Any]:
    """Exhaustively derive d_Z and T from the full parent Z-centralizer."""

    started = time.perf_counter()
    HX = np.asarray(parent["HX"], dtype=np.uint8)
    HZ = np.asarray(parent["HZ"], dtype=np.uint8)
    ell, m = int(parent["ell"]), int(parent["m"])
    n = int(HX.shape[1])
    k_parent = int(parent["k_P"])

    ker_basis = nullspace_np(HX)
    if int(ker_basis.shape[0]) > MAX_PARENT_KERNEL_DIMENSION:
        raise RuntimeError("sampled parent exceeded exact-kernel enumeration cap")
    masks = enumerate_span(ker_basis)
    if masks.shape[0] != 1 << int(ker_basis.shape[0]):
        raise AssertionError("nullspace-span enumeration cardinality mismatch")
    weights = np.bitwise_count(masks)

    sz_basis = Echelon(rows_to_masks(HZ))
    if sz_basis.rank != int(rank_np(HZ)):
        raise AssertionError("bitset and NumPy ranks disagree for S_Z")

    d_z: int | None = None
    minimum_logicals: list[int] = []
    for weight in range(1, n + 1):
        at_weight = masks[weights == weight]
        survivors = [int(mask) for mask in at_weight if not sz_basis.contains(int(mask))]
        if survivors:
            d_z = weight
            minimum_logicals = survivors
            break
    if d_z is None:
        raise AssertionError("positive-k CSS parent has no Z logical")

    module_masks: list[int] = []
    module_basis = sz_basis.copy()
    for logical in minimum_logicals:
        for translated in orbit_masks(logical, ell, m):
            module_masks.append(translated)
            module_basis.insert(translated)
    T = int(module_basis.rank - sz_basis.rank)
    if T > k_parent:
        raise AssertionError("minimum-logical module exceeds parent logical dimension")

    left_kernel = nullspace_np(HX.T)
    if int(left_kernel.shape[0]) * 2 != k_parent:
        raise AssertionError("left-kernel dimension is inconsistent with k_P")

    return {
        **parent,
        "left_kernel": left_kernel,
        "s_z_basis_masks": rows_to_masks(HZ),
        "minimum_logical_masks": minimum_logicals,
        "minimum_module_orbit_masks": module_masks,
        "d_Z": int(d_z),
        "T": T,
        "num_kernel_vectors_enumerated": int(masks.shape[0]),
        "num_minimum_Z_logicals": len(minimum_logicals),
        "num_minimum_logical_orbit_rows": len(module_masks),
        "wall_time_s": round(time.perf_counter() - started, 6),
    }


# ---------------------------------------------------------------------------
# Exact small PBB threshold decision: d_Q <= d_Z or d_Q > d_Z
# ---------------------------------------------------------------------------
_HALF_PATTERN_CACHE: dict[
    tuple[int, int], dict[int, tuple[np.ndarray, np.ndarray]]
] = {}


def half_patterns(
    qubits: int, cap: int
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """All ternary Pauli patterns of each exact symplectic weight <= cap."""

    key = (qubits, cap)
    cached = _HALF_PATTERN_CACHE.get(key)
    if cached is not None:
        return cached

    out: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for weight in range(cap + 1):
        total = comb(qubits, weight) * (3**weight)
        xs = np.empty(total, dtype=np.uint16)
        zs = np.empty(total, dtype=np.uint16)
        cursor = 0
        for support in itertools.combinations(range(qubits), weight):
            for states in itertools.product((1, 2, 3), repeat=weight):
                x_mask = 0
                z_mask = 0
                for position, state in zip(support, states, strict=True):
                    if state & 1:  # X or Y
                        x_mask |= 1 << position
                    if state & 2:  # Z or Y
                        z_mask |= 1 << position
                xs[cursor] = x_mask
                zs[cursor] = z_mask
                cursor += 1
        if cursor != total:
            raise AssertionError("ternary half-pattern enumeration mismatch")
        out[weight] = (xs, zs)
    _HALF_PATTERN_CACHE[key] = out
    return out


def xor_lookup(contributions: list[int]) -> np.ndarray:
    """Syndrome contribution of every binary subset of a short qubit half."""

    lookup = np.zeros(1 << len(contributions), dtype=np.uint64)
    for mask in range(1, lookup.shape[0]):
        low = mask & -mask
        bit = low.bit_length() - 1
        lookup[mask] = lookup[mask ^ low] ^ np.uint64(contributions[bit])
    return lookup


def column_masks(M: np.ndarray) -> list[int]:
    return [row_mask(np.asarray(M, dtype=np.uint8)[:, column]) for column in range(M.shape[1])]


def logical_at_or_below(H: np.ndarray, cap: int) -> dict[str, Any]:
    """Exact MITM search for a non-stabilizer centralizer vector of wt <= cap.

    For a check row (a|b), an X on qubit q contributes b_q to its
    symplectic syndrome and a Z contributes a_q.  Splitting qubits into two
    halves turns a zero-syndrome ternary error into an equality match of their
    packed GF(2) syndromes.  Every pattern of every total weight <= cap occurs
    exactly once, so no hit is a proof of d_Q > cap.
    """

    H = np.asarray(H, dtype=np.uint8)
    n = H.shape[1] // 2
    if H.shape[1] != 2 * n:
        raise ValueError("symplectic H must have an even number of columns")
    n_left = n // 2
    n_right = n - n_left
    h_x = H[:, :n]
    h_z = H[:, n:]

    # Syndrome of an X error is the check's Z column; vice versa for Z.
    x_contrib = column_masks(h_z)
    z_contrib = column_masks(h_x)
    left_x_lookup = xor_lookup(x_contrib[:n_left])
    left_z_lookup = xor_lookup(z_contrib[:n_left])
    right_x_lookup = xor_lookup(x_contrib[n_left:])
    right_z_lookup = xor_lookup(z_contrib[n_left:])
    patterns_left = half_patterns(n_left, cap)
    patterns_right = half_patterns(n_right, cap)
    stabilizers = Echelon(rows_to_masks(H))

    centralizer_pairs = 0
    stabilizers_discarded = 0
    for total_weight in range(1, cap + 1):
        for left_weight in range(total_weight + 1):
            right_weight = total_weight - left_weight
            if left_weight > n_left or right_weight > n_right:
                continue
            left_x, left_z = patterns_left[left_weight]
            right_x, right_z = patterns_right[right_weight]
            left_syndrome = (
                left_x_lookup[left_x.astype(np.intp)]
                ^ left_z_lookup[left_z.astype(np.intp)]
            )
            right_syndrome = (
                right_x_lookup[right_x.astype(np.intp)]
                ^ right_z_lookup[right_z.astype(np.intp)]
            )
            order = np.argsort(right_syndrome, kind="stable")
            sorted_right = right_syndrome[order]
            lo = np.searchsorted(sorted_right, left_syndrome, side="left")
            hi = np.searchsorted(sorted_right, left_syndrome, side="right")
            matched_left = np.flatnonzero(hi > lo)
            for li in matched_left:
                for ordered_ri in range(int(lo[li]), int(hi[li])):
                    ri = int(order[ordered_ri])
                    centralizer_pairs += 1
                    x_mask = int(left_x[li]) | (int(right_x[ri]) << n_left)
                    z_mask = int(left_z[li]) | (int(right_z[ri]) << n_left)
                    vector = x_mask | (z_mask << n)
                    if stabilizers.contains(vector):
                        stabilizers_discarded += 1
                        continue
                    kind = (
                        "pure_Z"
                        if x_mask == 0
                        else "pure_X"
                        if z_mask == 0
                        else "mixed"
                    )
                    return {
                        "logical_exists_at_or_below_cap": True,
                        "minimum_logical_weight_within_cap": total_weight,
                        "kind": kind,
                        "centralizer_pairs_examined": centralizer_pairs,
                        "stabilizers_discarded": stabilizers_discarded,
                    }
    return {
        "logical_exists_at_or_below_cap": False,
        "minimum_logical_weight_within_cap": None,
        "kind": None,
        "centralizer_pairs_examined": centralizer_pairs,
        "stabilizers_discarded": stabilizers_discarded,
    }


# ---------------------------------------------------------------------------
# Small-lattice perturbation enumeration
# ---------------------------------------------------------------------------
def defect_signature(defect: np.ndarray) -> int:
    return row_mask(np.asarray(defect, dtype=np.uint8).reshape(-1))


def singleton_defect_signatures(
    ell: int, m: int, A: np.ndarray, B: np.ndarray
) -> list[int]:
    """Linear AC^T+BD^T commutation-defect signature per C/D monomial."""

    dim = ell * m
    zero = np.zeros((dim, dim), dtype=np.uint8)
    signatures: list[int] = []
    for coordinate in range(2 * dim):
        monomial = monomial_matrix(
            ell, m, (coordinate % dim) // m, (coordinate % dim) % m
        )
        if coordinate < dim:
            defect = commutation_defect(A, B, monomial, zero)
        else:
            defect = commutation_defect(A, B, zero, monomial)
        signatures.append(defect_signature(defect))
    return signatures


def pbb_check(HX: np.ndarray, HZ: np.ndarray, CD: np.ndarray) -> np.ndarray:
    n = HX.shape[1]
    return np.vstack(
        (
            np.hstack((HX, CD)),
            np.hstack((np.zeros((HZ.shape[0], n), dtype=np.uint8), HZ)),
        )
    ).astype(np.uint8)


def probe_perturbations(parent: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Exhaust all commuting total-support<=4 perturbations for one parent."""

    started = time.perf_counter()
    ell, m = int(parent["ell"]), int(parent["m"])
    dim, n = ell * m, 2 * ell * m
    HX = np.asarray(parent["HX"], dtype=np.uint8)
    HZ = np.asarray(parent["HZ"], dtype=np.uint8)
    A, B = HX[:, :dim], HX[:, dim:]
    left_kernel = np.asarray(parent["left_kernel"], dtype=np.uint8)
    s_z_basis = Echelon(parent["s_z_basis_masks"])
    t_parent = int(parent["T"])
    d_z = int(parent["d_Z"])

    signatures = singleton_defect_signatures(ell, m, A, B)
    zero = np.zeros((dim, dim), dtype=np.uint8)
    raw_supports = 0
    commutation_valid = 0
    delta_positive = 0
    direct_defect_checks = 0
    delta_relation: Counter[str] = Counter()
    containment_counts: Counter[str] = Counter()
    distance_method: Counter[str] = Counter()
    d_q_relation: Counter[str] = Counter()
    significant_records: list[dict[str, Any]] = []
    counterexamples: list[dict[str, Any]] = []

    for weight in range(1, MAX_PERTURBATION_SUPPORT + 1):
        for combo in itertools.combinations(range(2 * dim), weight):
            raw_supports += 1
            signature = 0
            for coordinate in combo:
                signature ^= signatures[coordinate]
            if signature:
                continue
            commutation_valid += 1
            c_terms, d_terms = support_terms(combo, ell, m)
            C = poly_matrix(ell, m, list(c_terms)) if c_terms else zero
            D = poly_matrix(ell, m, list(d_terms)) if d_terms else zero
            direct_defect_checks += 1
            if commutation_defect(A, B, C, D).any():
                raise AssertionError("linear defect signature accepted a noncommuting PBB")
            CD = np.hstack((C, D)).astype(np.uint8)
            _, delta_rows = dressing_image(HX, CD)
            dim_delta = int(rank_np(delta_rows))
            bar_basis = Echelon(rows_to_masks(HZ))
            for row in delta_rows:
                bar_basis.insert(row_mask(row))
            dim_delta_bar = int(bar_basis.rank - s_z_basis.rank)
            if dim_delta_bar == 0:
                continue
            delta_positive += 1
            if dim_delta_bar < t_parent:
                relation = "less_than_T"
            elif dim_delta_bar == t_parent:
                relation = "equal_T"
            else:
                relation = "greater_than_T"
            delta_relation[relation] += 1

            survivor_mask = next(
                (
                    minimum
                    for minimum in parent["minimum_module_orbit_masks"]
                    if not bar_basis.contains(minimum)
                ),
                None,
            )
            module_contained = survivor_mask is None
            containment_counts["contained" if module_contained else "not_contained"] += 1
            k_q = int(n - rank_np(pbb_check(HX, HZ, CD)))
            if k_q != int(parent["k_P"]) - dim_delta_bar:
                raise AssertionError("PBB dimension identity failed on small lattice")

            record: dict[str, Any] = {
                "C": json_terms(c_terms),
                "D": json_terms(d_terms),
                "total_C_D_support": weight,
                "dim_Delta": dim_delta,
                "dim_Delta_bar": dim_delta_bar,
                "T": t_parent,
                "delta_bar_relation_to_T": relation,
                "k_Q": k_q,
                "M_bar_contained_in_Delta_bar": module_contained,
            }
            if not module_contained:
                # This translated minimum parent Z-logical is an explicit
                # pure-Z child logical of weight d_Z, so no full PBB search is
                # necessary to prove d_Q <= d_Z.
                if int(survivor_mask).bit_count() != d_z:
                    raise AssertionError("minimum-logical module retained a nonminimum row")
                record.update(
                    {
                        "d_Q_gt_d_Z": False,
                        "distance_decision": "surviving_minimum_pure_Z_logical",
                        "witness_weight": d_z,
                    }
                )
                distance_method["surviving_minimum_pure_Z_logical"] += 1
                d_q_relation["not_greater"] += 1
            else:
                H = pbb_check(HX, HZ, CD)
                threshold = logical_at_or_below(H, d_z)
                d_q_gt_d_z = not bool(threshold["logical_exists_at_or_below_cap"])
                record.update(
                    {
                        "d_Q_gt_d_Z": d_q_gt_d_z,
                        "distance_decision": "exhaustive_symplectic_MITM",
                        "threshold_search": threshold,
                    }
                )
                distance_method["exhaustive_symplectic_MITM"] += 1
                d_q_relation["greater" if d_q_gt_d_z else "not_greater"] += 1
                if d_q_gt_d_z and dim_delta_bar != t_parent:
                    counterexamples.append(record)

            # Retain all expensive-gate records plus every counterexample; the
            # summary below accounts for all other literal perturbations.
            if module_contained or record["d_Q_gt_d_Z"]:
                significant_records.append(record)

    summary = {
        "raw_total_supports_enumerated": raw_supports,
        "commutation_valid_supports": commutation_valid,
        "direct_commutation_checks": direct_defect_checks,
        "delta_positive_supports": delta_positive,
        "delta_bar_relation_to_T": dict(sorted(delta_relation.items())),
        "minimum_module_containment": dict(sorted(containment_counts.items())),
        "distance_decision_method": dict(sorted(distance_method.items())),
        "d_Q_vs_d_Z": dict(sorted(d_q_relation.items())),
        "significant_records": significant_records,
        "counterexamples": counterexamples,
        "wall_time_s": round(time.perf_counter() - started, 6),
    }
    return summary, counterexamples


def run_small_lattice_probe() -> dict[str, Any]:
    started = time.perf_counter()
    all_parent_records: list[dict[str, Any]] = []
    global_counterexamples: list[dict[str, Any]] = []
    lattice_records: list[dict[str, Any]] = []

    global_counts: Counter[str] = Counter()
    for protocol in LATTICE_PROTOCOL:
        ell, m = int(protocol["ell"]), int(protocol["m"])
        lattice_name = f"{ell}x{m}"
        primary_name = f"{protocol['primary_weights'][0]}-term x {protocol['primary_weights'][1]}-term"
        primary, primary_accounting = parent_candidates(
            ell, m, protocol["primary_weights"], primary_name
        )
        selected_pool = primary
        fallback_accounting: dict[str, Any] | None = None
        parent_family_used = primary_name
        if not selected_pool and protocol["fallback_weights"] is not None:
            fallback_name = (
                f"fallback {protocol['fallback_weights'][0]}-term x "
                f"{protocol['fallback_weights'][1]}-term"
            )
            selected_pool, fallback_accounting = parent_candidates(
                ell, m, protocol["fallback_weights"], fallback_name
            )
            parent_family_used = fallback_name

        sampled, sample_seed = stratified_sample(selected_pool, ell=ell, m=m)
        parent_summaries: list[dict[str, Any]] = []
        for sample_index, raw_parent in enumerate(sampled):
            parent = analyse_parent_exact(raw_parent)
            perturbations, counterexamples = probe_perturbations(parent)
            parent_id = f"{lattice_name}_{sample_index:03d}"
            parent_summary = {
                "id": parent_id,
                "ell": ell,
                "m": m,
                "family": parent["family"],
                "A_terms": json_terms(parent["A"]),
                "B_terms": json_terms(parent["B"]),
                "k_P": int(parent["k_P"]),
                "d_Z": int(parent["d_Z"]),
                "T": int(parent["T"]),
                "rank_HX": int(parent["rank_HX"]),
                "rank_HZ": int(parent["rank_HZ"]),
                "kernel_dimension": int(parent["kernel_dimension"]),
                "num_kernel_vectors_enumerated": int(parent["num_kernel_vectors_enumerated"]),
                "num_minimum_Z_logicals": int(parent["num_minimum_Z_logicals"]),
                "num_minimum_logical_orbit_rows": int(parent["num_minimum_logical_orbit_rows"]),
                "parent_bruteforce_wall_time_s": parent["wall_time_s"],
                "perturbations": perturbations,
            }
            parent_summaries.append(parent_summary)
            all_parent_records.append(parent_summary)
            for counterexample in counterexamples:
                global_counterexamples.append(
                    {"parent_id": parent_id, **counterexample}
                )
            global_counts["parents_sampled"] += 1
            global_counts["delta_positive_perturbations"] += int(
                perturbations["delta_positive_supports"]
            )
            global_counts["commutation_valid_supports"] += int(
                perturbations["commutation_valid_supports"]
            )
            global_counts["raw_supports"] += int(
                perturbations["raw_total_supports_enumerated"]
            )
            global_counts["d_Q_gt_d_Z"] += int(
                perturbations["d_Q_vs_d_Z"].get("greater", 0)
            )

        lattice_records.append(
            {
                "lattice": lattice_name,
                "ell": ell,
                "m": m,
                "primary_trinomial_family": primary_accounting,
                "fallback_family": fallback_accounting,
                "family_used_for_sample": parent_family_used,
                "sample_seed": sample_seed,
                "parent_sample_target": PARENTS_PER_LATTICE,
                "parents_sampled": len(parent_summaries),
                "parent_records": parent_summaries,
            }
        )

    return {
        "protocol": {
            "lattices": [[p["ell"], p["m"]] for p in LATTICE_PROTOCOL],
            "ell_times_m_max": 12,
            "parent_primary_family": "unordered trinomial support pairs",
            "nondegenerate_predicate": (
                "one Tanner connected component and k_P >= 2; exact enumeration "
                f"also requires kernel dimension <= {MAX_PARENT_KERNEL_DIMENSION}"
            ),
            "fallback_rule": (
                "only 5x2 uses a clearly labelled binomial fallback because its "
                "positive-k trinomial universe is empty"
            ),
            "parents_per_lattice": PARENTS_PER_LATTICE,
            "parent_selection": (
                "deterministic round-robin stratification over (k_P, polynomial "
                "family), seed derived from EXP-040 seed and lattice"
            ),
            "perturbation_class": (
                "every literal [C D] support of total C/D Hamming weight 1..4 "
                "whose exact AC^T + BD^T commutation defect is zero"
            ),
            "d_Z_method": (
                "enumerate the full GF(2) right-nullspace span of [A B], ordered "
                "by Hamming weight, modulo rowspace([B^T A^T])"
            ),
            "T_method": (
                "translation-orbit span of every brute-force minimum-weight "
                "Z-logical, modulo S_Z"
            ),
            "d_Q_threshold_method": (
                "exact meet-in-the-middle enumeration of all ternary Pauli errors "
                "of symplectic weight <= d_Z, filtered by centralizer and stabilizer "
                "rowspace; no SAT solver"
            ),
        },
        "aggregate": {
            **{key: int(value) for key, value in sorted(global_counts.items())},
            "counterexamples": len(global_counterexamples),
            "strict_distance_examples_exist": bool(global_counts["d_Q_gt_d_Z"]),
        },
        "lattices": lattice_records,
        "counterexamples": global_counterexamples,
        "wall_time_s": round(time.perf_counter() - started, 6),
    }


# ---------------------------------------------------------------------------
# Driver / verdict
# ---------------------------------------------------------------------------
def run() -> tuple[dict[str, Any], int]:
    started = time.perf_counter()
    reversals, reversal_clean, reversal_problems = reversal_records()
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "experiment": "EXP-040",
        "utc": utc_now(),
        "certified_reversals": reversals,
    }

    if not reversal_clean:
        paragraph = (
            "REFUTED: a certified reversal violates forced saturation or an "
            "independent reconstruction invariant; small-lattice probing was "
            "stopped as required. "
            + " | ".join(reversal_problems)
        )
        payload.update(
            {
                "verdict": "REFUTED",
                "verdict_paragraph": paragraph,
                "certified_reversal_problems": reversal_problems,
                "small_lattice": {"status": "not_run_after_certified_refutation"},
                "wall_time_s": round(time.perf_counter() - started, 6),
            }
        )
        atomic_write_json(OUT, payload)
        return payload, 1

    small_lattice = run_small_lattice_probe()
    counterexamples = small_lattice["counterexamples"]
    if counterexamples:
        verdict = "REFUTED"
        paragraph = (
            "REFUTED: all seven certified reversals saturate, but the attached "
            "small-lattice perturbation has d_Q > d_Z(P) with "
            "dim(Delta_bar) != T(P)."
        )
        status = 1
    else:
        verdict = "SUPPORTED"
        paragraph = (
            "SUPPORTED: all 7 certified reversals have dim(Delta_bar) = T(P), "
            "and the exact small-lattice sample found no counterexample."
        )
        status = 0

    payload.update(
        {
            "verdict": verdict,
            "verdict_paragraph": paragraph,
            "small_lattice": small_lattice,
            "wall_time_s": round(time.perf_counter() - started, 6),
        }
    )
    atomic_write_json(OUT, payload)
    return payload, status


def main() -> int:
    _, status = run()
    return status


if __name__ == "__main__":
    raise SystemExit(main())
